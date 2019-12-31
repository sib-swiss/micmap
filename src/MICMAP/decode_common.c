/*
 * ------------------------------------------------------------------------------------------------------------------------
 *
 *                                * micmap *
 *      Mapping of short reads in fastq format onto a reference
 *
 *
 *  Copyright (C) SIB  - Swiss Institute of Bioinformatics,  2015-2019 Nicolas Guex, Thierry Schuepbach and Christian Iseli
 *  Copyright (C) UNIL - University of Lausanne, Switzerland      2019 Nicolas Guex and Christian Iseli
 *
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *
 *      Code:       Nicolas Guex, Thierry Schuepbach and Christian Iseli
 *      Contacts:   Nicolas.Guex@unil.ch and Christian.Iseli@unil.ch
 *      Repository: https://github.com/sib-swiss/micmap
 *
 * ------------------------------------------------------------------------------------------------------------------------
 */
#define _GNU_SOURCE
#include "constants.h"
#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/statfs.h>
#include <fcntl.h>
#include <errno.h>
#include <pthread.h>
#include <string.h>
#include <assert.h>
#include <signal.h>
#include <libgen.h>
#include <mm_malloc.h>
#include "Decoder.h"

//--------------------------------------------------------------- 
// DEFINITIONS 
//--------------------------------------------------------------- 
#define FOUR_G (4UL * 1024 * 1024 * 1024)
#define MAX_NANOSEC 999999999
#define CEIL(X) ((X-(int)(X)) > 0 ? (int)(X+1) : (int)(X))
#define MAP_ADDR NULL
//(void*) FOUR_G
#define MAP_FLAGS MAP_SHARED
//(MAP_FIXED|MAP_SHARED)


//--------------------------------------------------------------- 
// STRUCTURE DEFINITIONS 
//---------------------------------------------------------------
typedef struct _thr_args {
	char *bufPtr;
	size_t size;
	unsigned int chunk;
	char fn[512];
} thr_arg_t, *thr_arg_p_t;


//--------------------------------------------------------------- 
// GLOBALS 
//--------------------------------------------------------------- 
// VALIDDATA ValidDataTable = {
// 	.valid_tbl             = MAP_FAILED,
// 	.following_nt_tbl      = MAP_FAILED,
// 	.valid_tbl_size        = kVALID_DATA_SIZE,
// 	.following_nt_tbl_size = kFOLLOWING_NT_DATA_SIZE,
// 	.valid_fileID          = 0,
// 	.following_nt_fileID   = 0
// };

/* contain table of chr,pos of each 16mer tag (note we have 8 successive such tables,
   corresponding to each of the 8 2nt followups). */
// unsigned int *chrpos_tbl[kFollowNtCount] = { [0 ... kFollowNtCount-1] = MAP_FAILED };

/* contain table of strandness of each unique exact match Nmers. */
// unsigned char *strand_tbl = MAP_FAILED;

extern int verbose;

//--------------------------------------------------------------- 
// LOCAL FCTS 
//--------------------------------------------------------------- 
static void * read_buffer(void *arg)
{
	thr_arg_p_t tap = (thr_arg_p_t) arg;
	int fd = open(tap->fn, O_RDONLY);
	size_t total = 0;
	if (lseek(fd, tap->size * tap->chunk, SEEK_SET) == -1) {
		fprintf(stderr, "Could not lseek: %s(%d)\n", strerror(errno), errno);
		return (void *) 1;
	}
	while (total < tap->size) {
		ssize_t rc;
		while ((rc = read(fd, tap->bufPtr, 32768)) == -1) {
			if (errno != EINTR && errno != EAGAIN) {
				fprintf(stderr, "Could not read: %s(%d)\n", strerror(errno), errno);
				return (void *) 1;
			}
		}
		total += rc;
		tap->bufPtr += rc;
	}
	close(fd);
	return 0;
} /* read_4G_buffer */
//---------------------------------------------------------------

static void waitForDecoderPool(DecoderPool_t * const restrict thpool){

	/* Continuous polling */
	double timeout = 1.0;
	time_t start, end;
	double tpassed = 0.0;
	time (&start);
	while (tpassed < timeout && (thpool->PoolJob_q.len || atomic_load(&(thpool->num_threads_working)) > 0) )
	{
		time (&end);
		tpassed = difftime(end,start);
	}

	/* Exponential polling */
	long init_nano =  1; /* MUST be above 0 */
	long new_nano;
	double multiplier = 1.01;
	int  max_secs = 20;

	struct timespec polling_interval;
	polling_interval.tv_sec  = 0;
	polling_interval.tv_nsec = init_nano;

	while (thpool->PoolJob_q.len || atomic_load(&(thpool->num_threads_working)) > 0)
	{
		nanosleep(&polling_interval, NULL);
		if ( polling_interval.tv_sec < max_secs ){
			new_nano = CEIL(polling_interval.tv_nsec * multiplier);
			polling_interval.tv_nsec = new_nano % MAX_NANOSEC;
			if ( new_nano > MAX_NANOSEC ) {
				polling_interval.tv_sec ++;
			}
		}
		else break;
	}

	/* Fall back to max polling */
	while (thpool->PoolJob_q.len || atomic_load(&(thpool->num_threads_working)) > 0 ) sleep(max_secs);
}
//---------------------------------------------------------------

static void terminateDecoderPool(DecoderPool_t * const restrict thpool) 
{
	volatile int threads_total = atomic_load(&(thpool->num_threads_alive));

	/* End each thread 's infinite loop */
	thpool->threads_keepalive = 0;

	/* Give one second to kill idle threads */
	double TIMEOUT = 1.0;
	time_t start, end;
	double tpassed = 0.0;
	time (&start);
	while (tpassed < TIMEOUT && atomic_load(&(thpool->num_threads_alive)) > 0) {
		bsem_post_all(thpool->PoolJob_q.has_items);
		time (&end);
		tpassed = difftime(end,start);
	}

	/* Poll remaining threads */
	while (atomic_load(&(thpool->num_threads_alive)) > 0)  {
		bsem_post_all(thpool->PoolJob_q.has_items);
		sleep(1);
	}	
} 
//---------------------------------------------------------------

//--------------------------------------------------------------- 
// FCTS 
//--------------------------------------------------------------- 

// That one could be inlined in the header
void freeValidTable(VALIDDATA * restrict ValidDataTable)
{
		if (ValidDataTable->following_nt_tbl != MAP_FAILED) {
			munmap(ValidDataTable->following_nt_tbl,kFOLLOWING_NT_DATA_SIZE);
			ValidDataTable->following_nt_tbl = MAP_FAILED;
		}
		if (ValidDataTable->valid_tbl != MAP_FAILED) {
			munmap(ValidDataTable->valid_tbl,kVALID_DATA_SIZE);
			ValidDataTable->following_nt_tbl = MAP_FAILED;
		}		
 		ValidDataTable = NULL;
}
//---------------------------------------------------------------

VALIDDATA * loadValidTable(const char * const restrict prefix, const char NT) 
{
	struct stat64 st;
	struct statfs64 stfs;
	char fn[512];
	
	VALIDDATA * restrict ValidDataTable = (VALIDDATA *) malloc(sizeof(VALIDDATA));
	if (ValidDataTable == NULL) return NULL;
	
	ValidDataTable->valid_tbl        = MAP_FAILED;
	ValidDataTable->following_nt_tbl = MAP_FAILED;

	strncpy(fn, prefix, 512);
	if (statfs64(dirname(fn), &stfs)) {
		perror("stat64");
		goto bail;
	}
	

	sprintf(fn,"%s%c.valid.bin",(char*)prefix,NT);
	int fd = open (fn, O_RDONLY);
	if (fd < 0)	{
		fprintf(stderr, "can't open bin file\n");
		goto bail;
	}
	
	if (fstat64(fd, &st) != 0)	{
		fprintf(stderr, "cannot stat file %s\n", fn);
		goto bail;
	}
	if (st.st_size != kVALID_DATA_SIZE)	{
		fprintf(stderr, "validStat.st_size %lx != %lx\n",st.st_size,kVALID_DATA_SIZE);
		goto bail;
	}

#if 0
	int flags;
	if (stfs.f_type != HUGETLBFS_MAGIC) {
		flags = MAP_PRIVATE|MAP_ANONYMOUS|MAP_HUGETLB|MAP_HUGE_1GB;
		ValidDataTable->valid_tbl = (unsigned char *)mmap (NULL, kVALID_DATA_SIZE, PROT_READ, flags, fd, 0);
		if (ValidDataTable->valid_tbl == MAP_FAILED) {
			fprintf(stderr, "Unable to map Valid Table to Huge TLB 1GB\n");
			flags = MAP_PRIVATE|MAP_ANONYMOUS|MAP_HUGETLB|MAP_HUGE_2MB;
			ValidDataTable->valid_tbl = (unsigned char *)mmap (NULL, kVALID_DATA_SIZE, PROT_READ, flags, fd, 0);
			if (ValidDataTable->valid_tbl == MAP_FAILED) {
				fprintf(stderr, "Unable to map Valid Table to Huge TLB 2MB\n");
				flags = MAP_PRIVATE;
			}
			else goto Done;
		}
		else goto Done;
	}
	else 
		flags = MAP_SHARED;
#else
	int flags = (stfs.f_type == HUGETLBFS_MAGIC) ? MAP_SHARED : MAP_PRIVATE;
#endif
	
	ValidDataTable->valid_tbl = (unsigned char *) mmap (NULL, kVALID_DATA_SIZE, PROT_READ, flags, fd, 0);
	if (ValidDataTable->valid_tbl == MAP_FAILED) {
		fprintf(stderr, "Valid Table mapping failed.\n");
		goto bail;
	}
	close(fd);

	sprintf(fn,"%s%c.%dnt.bin",(char*)prefix,NT,kFollowNtLen);
	fd = open (fn, O_RDONLY);
	if (fd < 0) {
		fprintf(stderr, "can't open bin file\n");
		goto bail;
	}
	
	if (fstat64(fd,&st) != 0) {
		fprintf(stderr, "cannot stat file %s\n", fn);
		goto bail;
	}

	if (st.st_size != kFOLLOWING_NT_DATA_SIZE) {
		fprintf(stderr, "ntStat.st_size %lx != %lx\n",st.st_size,kFOLLOWING_NT_DATA_SIZE);
		goto bail;
	}
#if 0
	if (stfs.f_type != HUGETLBFS_MAGIC) {
		flags = MAP_PRIVATE|MAP_ANONYMOUS|MAP_HUGETLB|MAP_HUGE_1GB;
		ValidDataTable->following_nt_tbl = (unsigned int *)mmap (0, kFOLLOWING_NT_DATA_SIZE, PROT_READ, flags, fd, 0);
		if (ValidDataTable->following_nt_tbl == MAP_FAILED) {
			fprintf(stderr, "Unable to map NT Table to Huge TLB 1GB\n");
			flags = MAP_PRIVATE|MAP_ANONYMOUS|MAP_HUGETLB|MAP_HUGE_2MB;
			ValidDataTable->following_nt_tbl = (unsigned int *)mmap (0, kFOLLOWING_NT_DATA_SIZE, PROT_READ, flags, fd, 0);
			if (ValidDataTable->following_nt_tbl == MAP_FAILED) {
				fprintf(stderr, "Unable to map Valid Table to Huge TLB 2MB\n");
				flags = MAP_PRIVATE;
			}
			else goto Done1;
		}
		else goto Done1;
	}
	else 
		flags = MAP_SHARED;
#else
	flags = (stfs.f_type == HUGETLBFS_MAGIC) ? MAP_SHARED : MAP_PRIVATE;
#endif
	
	ValidDataTable->following_nt_tbl = (unsigned int *)mmap (0, kFOLLOWING_NT_DATA_SIZE, PROT_READ, flags, fd, 0);
	if (ValidDataTable->following_nt_tbl == MAP_FAILED) {
		fprintf(stderr, "NT table mapping failed.\n");
		goto bail;
	}
	close(fd);
	
	return ValidDataTable;
	bail:;
		freeValidTable(ValidDataTable);
		return NULL;
}
//---------------------------------------------------------------

int loadDecodingTable(DecoderPool_t * const restrict thpool, const char * const restrict prefix,
                      const char NT, const _Bool useHugeTLB) 
{
	char temp[512];
	pthread_t t[kFollowNtCount+1];
	thr_arg_t ta[kFollowNtCount+1];
	struct statfs64 sb;
	
	strncpy(temp, prefix, 512);
	const char * const path = dirname(temp);
	if (statfs64(path, &sb)) {
		perror("statfs in loadDecodingTable");
		goto bail;
	}
	
	for (int i = 0; i < kFollowNtCount; i++) thpool->chrpos_tbl[i] = MAP_FAILED;
	thpool->strand_tbl = MAP_FAILED;
	
	const _Bool isHugeTLB = (sb.f_type == HUGETLBFS_MAGIC);
	
	if (isHugeTLB) {
		snprintf(temp, 512, "%s%c.chr.bin",prefix,NT);
		int fd = open(temp, O_RDONLY);
		if (fd < 0) {
			perror("open");
			goto bail;
		}
		
		struct stat64 st;
		if (fstat64(fd, &st)) {
			perror("stat64");
			goto bail;
		}
		
		if (st.st_size != kFollowNtCount*FOUR_G) {
			fprintf(stderr, "Invalid size for file %s\n", temp);
			goto bail;
		}
		
		size_t offset = 0UL;
		for (int i = 0; i < kFollowNtCount; i++) {
			thpool->chrpos_tbl[i] = mmap(NULL, FOUR_G, PROT_READ, MAP_SHARED, fd, offset);
			if (thpool->chrpos_tbl[i]  == MAP_FAILED) {
				fprintf(stderr, "Unable to map chunk %i of %s\n", i, temp);
				perror("mmap");
				goto bail;
			}
			offset += FOUR_G;
		}
		close(fd);
		
		snprintf(temp, 512, "%s%c.strand.bin", prefix, NT);
		fd = open(temp, O_RDONLY);
		if (fd < 0) {
			perror("open");
			goto bail;
		}
		if (fstat64(fd, &st)) {
			perror("stat64");
			goto bail;
		}
		
		if (st.st_size != kSTRAND_DATA_SIZE) {
			fprintf(stderr, "Invalid size for file %s\n", temp);
			goto bail;
		}
		
		thpool->strand_tbl = mmap(NULL, kSTRAND_DATA_SIZE, PROT_READ, MAP_SHARED, fd, 0);
		if (thpool->strand_tbl == MAP_FAILED) {
			fprintf(stderr, "Unable to map %s\n", temp);
			goto bail;
		}
		close(fd);
		
		thpool->SharedMemoryFileDecriptor = fd;
	}
	else {
		if (!useHugeTLB) goto standard;
		{
			for (int i = 0; i < kFollowNtCount; i++) {
				thpool->chrpos_tbl[i] = mmap(NULL, FOUR_G, PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS|MAP_HUGETLB|MAP_HUGE_1GB, -1, 0);
				if (thpool->chrpos_tbl[i] == MAP_FAILED) {
					fprintf(stderr, "Could not mmap on Huge TLB 1GB for loading decoder thread %d: %s(%d), trying 2MB...\n", i, strerror(errno), errno);
					
					thpool->chrpos_tbl[i] = mmap(NULL, FOUR_G, PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS|MAP_HUGETLB|MAP_HUGE_2MB, -1, 0);
					if (thpool->chrpos_tbl[i] == MAP_FAILED) {
						fprintf(stderr, "Could not mmap on Huge TLB 2MB for loading decoder thread %d: %s(%d), trying std...\n", i, strerror(errno), errno);
						goto clean;
					}
				}
			}
			thpool->strand_tbl = mmap(NULL, kSTRAND_DATA_SIZE, PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS|MAP_HUGETLB|MAP_HUGE_1GB, -1, 0);
			if (thpool->strand_tbl == MAP_FAILED) {
				fprintf(stderr, "Could not mmap on Huge 1GB strand for loading decoder thread : %s(%d)\n", strerror(errno), errno);
				thpool->strand_tbl = mmap(NULL, kSTRAND_DATA_SIZE, PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS|MAP_HUGETLB|MAP_HUGE_2MB, -1, 0);
				if (thpool->strand_tbl == MAP_FAILED) {
					fprintf(stderr, "Could not mmap on Huge 2MB strand for loading decoder thread : %s(%d), trying std...\n", strerror(errno), errno);
					goto clean;
				}
			}
			goto load;
		}
		clean:;
		{
			for (int i = 0; i < kFollowNtCount; i++) {
				if (thpool->chrpos_tbl[i] != MAP_FAILED) munmap((void*)thpool->chrpos_tbl[i], FOUR_G);
			}
			if (thpool->strand_tbl != MAP_FAILED) munmap((void*)thpool->strand_tbl, kSTRAND_DATA_SIZE);
		}
		standard:;
		{
			for (int i = 0; i < kFollowNtCount; i++) {
				thpool->chrpos_tbl[i] = mmap(NULL, FOUR_G, PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);
				if (thpool->chrpos_tbl[i] == MAP_FAILED) {
					fprintf(stderr, "Could not mmap for loading decoder thread %d: %s(%d)\n", i, strerror(errno), errno);
					goto bail;
				}
			}
			thpool->strand_tbl = mmap(NULL, kSTRAND_DATA_SIZE, PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, -1, 0);
			if (thpool->strand_tbl == MAP_FAILED) {
				fprintf(stderr, "Could not mmap strand for loading decoder thread : %s(%d)\n", strerror(errno), errno);
				goto bail;
			}
		}
		load:;
		for (int i = 0; i < kFollowNtCount; i++) {
			ta[i].bufPtr = (char*)thpool->chrpos_tbl[i];
			ta[i].chunk = i;
			ta[i].size = FOUR_G;
			sprintf(ta[i].fn,"%s%c.chr.bin",prefix,NT);
			if (pthread_create(t + i, NULL, read_buffer, ta + i) != 0) {
				fprintf(stderr, "Could not create thread %d: %s(%d)\n", i, strerror(errno), errno);
				goto bail;
			}
		}
		ta[kFollowNtCount].bufPtr = (char*) thpool->strand_tbl;
		ta[kFollowNtCount].chunk = 0;
		ta[kFollowNtCount].size = kSTRAND_DATA_SIZE;
		sprintf(ta[kFollowNtCount].fn,"%s%c.strand.bin",prefix,NT);
		if (pthread_create(t + kFollowNtCount, NULL, read_buffer, ta + kFollowNtCount) != 0) {
			fprintf(stderr, "Could not create thread %d: %s(%d)\n", kFollowNtCount, strerror(errno), errno);
			goto bail;
		}
		
		for (int i = 0; i < (kFollowNtCount+1); i++) {
			void *res;
			if (pthread_join(t[i], &res) != 0) {
				fprintf(stderr, "Could not join thread %d: %s(%d)\n", i, strerror(errno), errno);
				return 1;
			}
			free(res);
		}
		thpool->SharedMemoryFileDecriptor = -1;
	}
	
	
	return 0;
	
	bail:
	return 1;
}
//---------------------------------------------------------------

void freeDecodingTable(DecoderPool_t * const restrict thpool) 
{ 
	for (int i = 0; i < kFollowNtCount; i++) { 
		if (thpool->chrpos_tbl[i] != MAP_FAILED) munmap((void*)thpool->chrpos_tbl[i],FOUR_G);
	}
	if (thpool->strand_tbl != MAP_FAILED) munmap((void*)thpool->strand_tbl, kSTRAND_DATA_SIZE);

	for (int i = 0; i < kFollowNtCount; i++) thpool->chrpos_tbl[i] = MAP_FAILED;
	thpool->strand_tbl = MAP_FAILED;
} 
//---------------------------------------------------------------

int DecodingTableToRAM(const char * const restrict path, const char * const restrict prefix,
                       const char NT)
{
	struct statfs64 stfs;
	thr_arg_t ta;
	struct DecodingTableFile {
			size_t Size;
			char FileName[24];
	} DecodingTableFiles[] = { 
			{ (4UL*1024*1024*1024),  "start0.?.2nt.bin" },
			{ (32UL*1024*1024*1024), "start0.?.chr.bin" },
			{ (1UL*1024*1024*1024),  "start0.?.strand.bin" },
			{ (1UL*1024*1024*1024),  "start0.?.valid.bin" },
	};
	int i = -1;
	char fn[512];	
	
	if(statfs64(path, &stfs)) {
		perror("statfs64");
		goto bail;
	}
	
	if (stfs.f_type == HUGETLBFS_MAGIC) {
		i = -1;
		char temp[512];
		strncpy(temp, prefix, 512);
		const char * const restrict basePrefix = dirname(temp);
		while (++i<4) {
			DecodingTableFiles[i].FileName[7] = NT;
			snprintf(fn, 512, "%s/%s", path, DecodingTableFiles[i].FileName);
			int fd = open(fn, O_EXCL | O_CREAT | O_RDWR, S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH);
			if (fd < 0) {
				perror("open");
				goto bail;
			}
			if (ftruncate64(fd, DecodingTableFiles[i].Size)) {
				fprintf(stderr, "Unable to extend file %s to %lu\n", DecodingTableFiles[i].FileName, DecodingTableFiles[i].Size);
				goto bail;
			}
			ta.bufPtr = mmap64(NULL, DecodingTableFiles[i].Size, PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0);
			if (ta.bufPtr == MAP_FAILED) { perror("mmap64"); goto bail; }
			ta.chunk = 0;
			snprintf(ta.fn, 512, "%s/%s", basePrefix, DecodingTableFiles[i].FileName);
			ta.size = DecodingTableFiles[i].Size;
			void * res = read_buffer(&ta);
			munmap(ta.bufPtr,  DecodingTableFiles[i].Size);
			close(fd);
		}
	}
	else {
		fprintf(stderr, "%s does not point to a huge TLB file system\n", path);
		goto bail;
	}
	return 0;
	
	bail:
		while (i>=0) {
			snprintf(fn, 512, "%s/%s", path, DecodingTableFiles[i--].FileName);
			unlink(fn);
		}
		return 1;
}
//---------------------------------------------------------------

int compare_res_element(const DecoderElement_t * const restrict a, const DecoderElement_t * const restrict b)
{
	const register unsigned int reversedA = (a->info & 0x01) ^ ((a->info >> 1) & 0x01);
	const register unsigned int reversedB = (b->info & 0x01) ^ ((b->info >> 1) & 0x01);
	if (reversedA != reversedB)
	{
		// keep forward strand first
		if (reversedA == 0)
			return -1;
		return 1;
	}
	if (a->which16mer < b->which16mer)
		return -1;
	if (a->which16mer > b->which16mer)
		return 1;
	return 0;
}
//---------------------------------------------------------------

DecoderMemoryBulk_t * allocateDecoderMemory(const size_t nBlocks)
{
	/* Allocates the return structure */
	DecoderMemoryBulk_t * const restrict data = (DecoderMemoryBulk_t*) malloc(sizeof(DecoderMemoryBulk_t));
	if (data == NULL) goto err1;
	
	/* Allocates the queue itself */
	jobqueue_t * const restrict queue = (jobqueue_t*) malloc(sizeof(jobqueue_t));
	if (queue == NULL) goto err2;
	
	if (jobqueue_init(queue) == -1) goto err3;

	if (jobqueue_init(&(data->ToBeAligned_q)) == -1) goto err3a;

	
	/* Allocates the queue elements */
	DecoderMemoryElement_t * const restrict elems = (DecoderMemoryElement_t *) malloc(nBlocks*sizeof(DecoderMemoryElement_t));
	if (elems == NULL) goto err3;
	
	/* Allocates the blocks */
	DecoderElement_t (*blocks)[kMaxTagsPerBuffer][kDesiredHitCnt] = (DecoderElement_t (*)[kMaxTagsPerBuffer][kDesiredHitCnt]) malloc(nBlocks*kDesiredHitCnt*kMaxTagsPerBuffer*sizeof(DecoderElement_t));
	if (blocks == NULL) goto err4;
	atomic_uint * counters = (atomic_uint*) malloc(nBlocks*sizeof(atomic_uint));
	if (counters == NULL) goto err5;
	
	assert((kDecoderTagBatch % kDesiredHitCnt) == 0);
	
	/* Populate the queue */
	pthread_mutex_lock(&queue->rwmutex);
	for (size_t i=0; i<nBlocks; i++) {
		elems[i].Anchors = &blocks[i];
		// atomic_init(&counters[i], 0U); Not always recognized by compiler, should be safe to avoid
		elems[i].count = &counters[i];
		jobqueue_push(queue, (job_t*) &elems[i]);
	}
	pthread_mutex_unlock(&queue->rwmutex);
	data->DecoderMemory_q = queue;
	data->Anchors = blocks;
	data->DecoderMemoryElements = elems;
	
	return data;
err5: ;
	free(counters);
err4: ;
	free(elems);
err3a: ;
	jobqueue_destroy(queue);
err3: ;
	free(queue);
err2: ;
	free(data);
err1: ;
	return NULL;
}
//---------------------------------------------------------------

void freeDecoderMemory(DecoderMemoryBulk_t * data)
{
		if (data->Anchors) free(data->Anchors);
		if(data->DecoderMemoryElements) free(data->DecoderMemoryElements);
		if (data->DecoderMemory_q) {
			jobqueue_destroy(data->DecoderMemory_q);
			free(data->DecoderMemory_q);
		}
}
//---------------------------------------------------------------

int createDecoderPool(DecoderPool_t * const restrict thpool, const DecodeTagsPtr Fct,
                      const Affinity_Mask_t * const restrict affinities,
                      DecoderMemoryBulk_t * const restrict InOut, const unsigned int nThreads)
{
	if (thpool == NULL) goto err1;
	
	thpool->threads = (pthread_t*) malloc(nThreads*sizeof(pthread_t));
	if (thpool->threads == NULL) goto err1;
	
	thpool->threads_keepalive = 1;
	
	thpool->num_threads = nThreads;
	atomic_store(&(thpool->num_threads_alive), 0U);
	atomic_store(&(thpool->num_threads_working), 0U);
	thpool->ToBeAligned_q = &(InOut->ToBeAligned_q);
	
	if (jobqueue_init(&(thpool->PoolJob_q)) == -1) goto err2;
	if (jobqueue_init(&(thpool->PoolMemory_q)) == -1) goto err3;
		
	/* Create free memory jobs and place them into donequeue, EXTRA could be NICE ?*/
	char * restrict Jobs = (void*) malloc((2*nThreads)*sizeof(DecoderJobElement_t));
	if (Jobs == NULL) goto err4;
	
	thpool->jobs = (job_t*) Jobs;
	pthread_mutex_lock(&thpool->PoolMemory_q.rwmutex);
	for(int iJob=0; iJob<(2*nThreads); iJob++) {
		jobqueue_push(&thpool->PoolMemory_q, (job_t*) Jobs);
		Jobs += sizeof(DecoderJobElement_t);
	}
	pthread_mutex_unlock(&thpool->PoolMemory_q.rwmutex);
	
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	//pthread_attr_setstacksize(&attr, (size_t) (2*1024*1024));
	if (affinities) {
		for (int iThread = 0; iThread<nThreads; iThread++) {
			if(pthread_attr_setaffinity_np(&attr, sizeof(Affinity_Mask_t), (cpu_set_t*) &affinities[iThread])) {
				perror("pthread_attr_setaffinity_np");
				while (--iThread > 0) {
                                        pthread_kill(thpool->threads[iThread], SIGKILL);
                                }
				goto err5;
			}
			if (pthread_create(&(thpool->threads[iThread]), &attr, (void* (*)(void*)) Fct, thpool) != 0) {
				fprintf(stderr, "Error creating decode worker thread %i\n", iThread);
				while (--iThread > 0) {
					pthread_kill(thpool->threads[iThread], SIGKILL);
				}
				goto err5;
			}
		}
	}
	else {
		for (int iThread = 0; iThread<nThreads; iThread++) {
			if (pthread_create(&(thpool->threads[iThread]), &attr, (void* (*)(void*)) Fct, thpool) != 0) {
				while (--iThread > 0) {
					pthread_kill(thpool->threads[iThread], SIGKILL);
				}
				goto err5;
			}
		}
	}
	return 0;
	
	err5:
		jobqueue_destroy(&thpool->PoolJob_q);
	err4:
		jobqueue_destroy(&thpool->PoolMemory_q);
	err3:
	err2:
		free(thpool->threads);
	err1:
	return 1;
}
//---------------------------------------------------------------

void destroyThreadPool(DecoderPool_t * restrict thpool)
{
	/* No need to destory if it's NULL */
	if (thpool == NULL) return ;

	/* End each thread 's infinite loop */
	thpool->threads_keepalive = 0;

	/* Give one second to kill idle threads */
	const double TIMEOUT = 1.0;
	time_t start, end;
	double tpassed = 0.0;
	time (&start);
	while (tpassed < TIMEOUT && atomic_load(&(thpool->num_threads_alive)) > 0) {
		bsem_post_all(thpool->PoolJob_q.has_items);
		time (&end);
		tpassed = difftime(end,start);
	}

	/* Poll remaining threads */
	while (atomic_load(&(thpool->num_threads_alive)) > 0 ) {
		bsem_post_all(thpool->PoolJob_q.has_items);
		sleep(1);
	}

	/* Job queue cleanup */
	jobqueue_destroy(&thpool->PoolJob_q);
	if (thpool->jobs) jobqueue_destroy(&thpool->PoolMemory_q);
	
	if (thpool->jobs) {
		free(thpool->jobs);
	}
	
	free(thpool->threads);
	free(thpool);
	thpool = NULL;
}
//---------------------------------------------------------------

void* ReaderToDecoderThread(ReaderToDecoderArgs_t * const restrict data)
{
	jobqueue_t * const restrict ReaderJobQueue = &(data->ReaderData->ReaderJobQueue);
	DecoderPool_t * const restrict DecoderPool = data->DecoderPool;
	jobqueue_t * const restrict DecoderMemory = data->DecoderMemory_q;

	// Initialize some other data
	ReaderJob_t * restrict task = NULL;

	//////////////////////////////////////////////////////////////////////////////////////////////
	// Transfer jobs from reader to internal decoders
	while(1) {
		DecoderMemoryElement_t * restrict space;
		bsem_wait(ReaderJobQueue->has_items);
		
		/* Read job from queue and execute it */
		pthread_mutex_lock(&(ReaderJobQueue->rwmutex));
		task = (ReaderJob_t *) jobqueue_pull(ReaderJobQueue);
		pthread_mutex_unlock(&(ReaderJobQueue->rwmutex));
		
		/***************************************************************/
		/*                     START THE TASK                          */
		/***************************************************************/
		if (task) {
			if (verbose & 0x1)
				printf("%s got a new batch of %u tags\n", __FUNCTION__, task->TagCnt);
			//////////////////////////////////////////////////////////////////////////////////////////////
			// Aquire a Batch of memory to accomodate all tags matches
			do {
				bsem_wait(DecoderMemory->has_items);
			
				/* Get memory slot from pool and zero atomic counter */
				pthread_mutex_lock(&(DecoderMemory->rwmutex));
				space = (DecoderMemoryElement_t *) jobqueue_pull(DecoderMemory);
				pthread_mutex_unlock(&(DecoderMemory->rwmutex));
			} while (space == NULL);
			
			DecoderElement_t (* restrict outElement)[kDesiredHitCnt] = &(space->Anchors[0][0]);
			atomic_uint * const restrict Donecount = space->count;
			atomic_store(Donecount, 0U);
			
			//////////////////////////////////////////////////////////////////////////////////////////////
			// Loop on tags per decoder batch size : kDecoderTagBatch
			unsigned int count = 0U;
			const unsigned int TagCnt = task->TagCnt;
			if (TagCnt == 0U) break;
			const ReaderTag_t * restrict Tags = task->Tags;
			while( count < TagCnt ) {
				////////////////////////////////////////////////////////////////////////////////////////////
				// Aquire emtpy Decoder memory slot
				bsem_wait(DecoderPool->PoolMemory_q.has_items);
				pthread_mutex_lock(&(DecoderPool->PoolMemory_q.rwmutex));
				DecoderJobElement_t * const restrict work = (DecoderJobElement_t*) jobqueue_pull(&(DecoderPool->PoolMemory_q));	
				pthread_mutex_unlock(&(DecoderPool->PoolMemory_q.rwmutex));
				
				////////////////////////////////////////////////////////////////////////////////////////////
				// Fillin
				work->Decoder.Anchors = outElement;
				work->count = Donecount;
				work->nToBeProcessed = ((count + kDecoderTagBatch) < TagCnt ) ? kDecoderTagBatch : (TagCnt - count);
				work->Reader.Tags = Tags;
				
				////////////////////////////////////////////////////////////////////////////////////////////
				// Submitting
				//printf("%s adding task 0x%lx, count %u, in 0x%lx, out 0x%lx\n", __FUNCTION__, work, work->nToBeProcessed, Tags, outElement);
				pthread_mutex_lock(&(DecoderPool->PoolJob_q.rwmutex));
				jobqueue_push(&(DecoderPool->PoolJob_q), (job_t*) work);
				pthread_mutex_unlock(&(DecoderPool->PoolJob_q.rwmutex)); 
				count += kDecoderTagBatch;
				outElement += kDecoderTagBatch;
				Tags += kDecoderTagBatch;
			}
			
			//////////////////////////////////////////////////////////////////////////////////////////////
			// Create a dummy last job to have the last one sending the decoding data further
			// down the pipeline
			{
				////////////////////////////////////////////////////////////////////////////////////////////
				// Aquire emtpy Decoder memory slot
				bsem_wait(DecoderPool->PoolMemory_q.has_items);
				pthread_mutex_lock(&(DecoderPool->PoolMemory_q.rwmutex));
				DecoderJobElement_t * const restrict work = (DecoderJobElement_t*) jobqueue_pull(&(DecoderPool->PoolMemory_q));			
				pthread_mutex_unlock(&(DecoderPool->PoolMemory_q.rwmutex));
				
				////////////////////////////////////////////////////////////////////////////////////////////
				// Fillin
				work->Decoder.DecoderJob = space;
				work->count = Donecount;
				work->nToBeProcessed = 0U;
				work->Reader.ReaderJob = task;
				
				////////////////////////////////////////////////////////////////////////////////////////////
				// Submitting
				//printf("%s adding task 0x%lx, count %u, Decode space @ 0x%lx, ReaderJob @ 0x%lx\n", __FUNCTION__, work, work->nToBeProcessed, space, task);
				pthread_mutex_lock(&(DecoderPool->PoolJob_q.rwmutex));
				jobqueue_push(&(DecoderPool->PoolJob_q), (job_t*) work);
				pthread_mutex_unlock(&(DecoderPool->PoolJob_q.rwmutex));
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	// Wait for internal decoders to terminate and send all jobs to the aligners
	printf("%s: Terminating decoder pool\n", __FUNCTION__);
	waitForDecoderPool(DecoderPool);
	terminateDecoderPool(DecoderPool);
	
 	//////////////////////////////////////////////////////////////////////////////////////////////
 	// Push End of Sequence to aligners
// 	{
// 		bsem_wait(DecoderPool->PoolMemory_q.has_items);
// 		pthread_mutex_lock(&(DecoderPool->PoolMemory_q.rwmutex));
// 		DecoderJobElement_t * const restrict work = (DecoderJobElement_t*) jobqueue_pull(&(DecoderPool->PoolMemory_q));			
// 		pthread_mutex_unlock(&(DecoderPool->PoolMemory_q.rwmutex));
// 		
// 		work->Decoder.DecoderJob = NULL;
// 		work->count = NULL;
// 		work->nToBeProcessed = 0U;
// 		work->Reader.ReaderJob = NULL;
// 		
// 		pthread_mutex_lock(&(DecoderPool->ToBeAligned_q->rwmutex));
// 		jobqueue_push(DecoderPool->ToBeAligned_q, (job_t*) work);
// 		pthread_mutex_unlock(&(DecoderPool->ToBeAligned_q->rwmutex));
// 	}
	bsem_post(DecoderPool->ToBeAligned_q->has_items);	
	return (void*) 0;
}
//---------------------------------------------------------------

/* vim: tabstop=2 shiftwidth=2
 */
