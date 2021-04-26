/*
 * ------------------------------------------------------------------------------------------------------------------------
 *
 *                                * micmap *
 *      Mapping of short reads in fastq format onto a reference
 *
 *
 *  Copyright (C) SIB  - Swiss Institute of Bioinformatics,                2015-2019 Nicolas Guex, Thierry Schuepbach and Christian Iseli
 *  Copyright (C) UNIL - University of Lausanne, Switzerland                    2019 Nicolas Guex and Christian Iseli
 *  Copyright (C) EPFL - Ecole Polytechnique Fédérale de Lausanne, Switzerland  2021 Christian Iseli
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
 *      Contacts:   Nicolas.Guex@unil.ch and christian.iseli@epfl.ch
 *      Repository: https://github.com/sib-swiss/micmap
 *
 * ------------------------------------------------------------------------------------------------------------------------
 */
#include "constants.h"
#include "config.h"
#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <pthread.h>
#include <signal.h>
#include <sys/mman.h>
#include <alloca.h>
#include <unistd.h>

#include "global_var.h"
#include "Aligner.h"
#include "Decoder.h"
#include "cpuThreadAlignPool.h"

#define HUGE_2MB_PAGE_SIZE (2*1024*1024)
#define MAX_NANOSEC 999999999
#define CEIL(X) ((X-(int)(X)) > 0 ? (int)(X+1) : (int)(X))

extern cpuAlignerPool_t * CPU_Pools;
//static cpuSpace_t TriggerDone = { .OutBufferPtr=NULL, .length=0LU, .BelongsTo_q=NULL, .Used=0LU, .Stored=0LU };

void* cpuBlockAligners(AlignerBlockArgs_t * const restrict data)
{
	cpuAlignerPool_t Pool;
	char * restrict OutputMemoryBlocks;
	union Statistics * restrict worker_stat;
	int err;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Tell the global aligner pool we are working, there will be nInputBlock -1
	pthread_mutex_lock(&(data->Pool->thcount_lock));
	data->Pool->num_threads_alive++;
	pthread_mutex_unlock(&(data->Pool->thcount_lock));

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Initialization of cpu pool data
	initRevCompTable(gIUPACrevcomp);
		
	Pool.AlignerBlock = data;
	Pool.num_threads = data->nAligners;
	Pool.threads_keepalive = 1;
	atomic_store(&Pool.num_threads_alive, 0);
	atomic_store(&Pool.JobIndex, -1);
	atomic_store(&Pool.OutIndex, -1);
	Pool.job_q.front = NULL;
   	Pool.job_q.len = 0;
	Pool.job_q.rear = NULL;
	atomic_store(&(Pool.job_q.rwmutex), 0);
	
	jobqueue_init(&(Pool.space_q));
	
#ifdef BENCHMARKS
	memset(ProcessBenchmarks, 0, gRecvCnt*gSlaveCnt*sizeof(Benchmarks_t));
	memset(ReaderBenchmarks, 0, gRecvCnt*sizeof(Benchmarks_t));
#endif

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Allocating  memory
	/* Output buffers */
	const size_t OutBlockSize = ( data->OutputBlockSize + (size_t) (HUGE_2MB_PAGE_SIZE-1)) & ~((size_t) (HUGE_2MB_PAGE_SIZE-1));
	OutputMemoryBlocks = (char *) mmap(NULL, data->nOutputBlock*OutBlockSize, PROT_READ | PROT_WRITE,
	                                   MAP_PRIVATE | MAP_ANONYMOUS | MAP_POPULATE | MAP_HUGETLB,
	                                   -1, 0);
	if (OutputMemoryBlocks == MAP_FAILED) {
		fprintf(stderr, "Cannot allocate Huge TLB for %u input blocks of %lu Bytes, trying standard...\n", data->nOutputBlock, OutBlockSize);
		OutputMemoryBlocks = (char *) mmap(NULL, data->nOutputBlock*OutBlockSize, PROT_READ | PROT_WRITE,
	                                   MAP_PRIVATE | MAP_ANONYMOUS | MAP_POPULATE, -1, 0);
		if (OutputMemoryBlocks == MAP_FAILED) {
			fprintf(stderr, "Cannot allocate memory for %u input blocks of %lu Bytes\n", data->nOutputBlock, OutBlockSize);
			goto bail;
		}
		//if (madvise(OutputMemoryBlocks, data->nOutputBlock*OutBlockSize, MADV_HUGEPAGE | MADV_WILLNEED)) perror("madise");
	}
	else {
		printf("CPU block aligners allocated %u output block of %lu Bytes each.\n", data->nOutputBlock, OutBlockSize);
	}
	printf("CPU block aligners seeking within genome chunk size of %u\n", data->GenomeChunkSize);
	
	Pool.OutputBuffers = (cpuSpace_t*) malloc(data->nOutputBlock*sizeof(cpuSpace_t));
	pthread_mutex_lock(&Pool.space_q.rwmutex);
	for (unsigned int i=0; i<data->nOutputBlock; i++) {
		////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Set id and memory parameters
		Pool.OutputBuffers[i].OutBufferPtr = &OutputMemoryBlocks[i*OutBlockSize];
		Pool.OutputBuffers[i].index = i;
		Pool.OutputBuffers[i].prev = NULL;
		Pool.OutputBuffers[i].BelongsTo_q = NULL;
		atomic_store(&(Pool.OutputBuffers[i].Used), 0UL);
		atomic_store(&(Pool.OutputBuffers[i].Stored), 0UL);
		
		////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Push to queue
		jobqueue_push(&(Pool.space_q), (job_t*) &Pool.OutputBuffers[i]);
	}
	pthread_mutex_unlock(&Pool.space_q.rwmutex);
	
	/* 3. Statistics data */
	Pool.stats = (union Statistics*) _mm_malloc(Pool.num_threads*sizeof(Statistics_t), 64);
	if (Pool.stats == NULL) {
			puts("CPU block aligner cannot allocate statistics table\n");
			goto bail;
	}
	const __m128i __zero = _mm_setzero_si128();
	for (unsigned int i=0; i<Pool.num_threads; i++) {
		Pool.stats[i].Vector[0] = __zero;
		Pool.stats[i].Vector[1] = __zero;
		Pool.stats[i].Vector[2] = __zero;
		Pool.stats[i].Vector[3] = __zero;
	}
		
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Filling up the work jobs
	fputs("Filling up input buffers\n", stdout);
	Pool.Jobs = (cpuAlignerJob_t*) alloca(data->nInputBlock*sizeof(cpuAlignerJob_t));
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Set id and memory parameters
	Pool.Jobs[0].index = 0;
	atomic_store(&(Pool.Jobs[0].Treated), 0UL);
	atomic_store(&(Pool.Jobs[0].Stored), 0UL);

	// WARNING: Less block than nInputBlock and we are dead stalling here;
	DecoderJobElement_t * DecoderJob;
	do {
		sleep(1);
		bsem_wait(data->Pool->ToBeAligned_q->has_items);
		printf("CPU block aligner requesting input buffer %u loading...\n", 0); fflush(stdout);
		pthread_mutex_lock(&(data->Pool->ToBeAligned_q->rwmutex));
		DecoderJob = (DecoderJobElement_t *) jobqueue_pull(data->Pool->ToBeAligned_q);
		pthread_mutex_unlock(&(data->Pool->ToBeAligned_q->rwmutex));
	} while (DecoderJob == NULL);

	assert(DecoderJob != NULL);
	Pool.Jobs[0].TagsData = DecoderJob->Reader.ReaderJob->Tags;
	Pool.Jobs[0].Anchors = (const DecoderElement_t (*)[16]) DecoderJob->Decoder.DecoderJob->Anchors;
	Pool.Jobs[0].TagCnt = DecoderJob->Reader.ReaderJob->TagCnt;
	Pool.Jobs[0].DecoderJob = DecoderJob;

	pthread_mutex_lock(&(data->Pool->thcount_lock));
	data->Pool->num_threads_working++;
	pthread_mutex_unlock(&(data->Pool->thcount_lock));

	printf("CPU block got data with %u tags\n", Pool.Jobs[0].TagCnt);

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Push to queue
	pushToQueue(&(Pool.job_q), (Task_t*) &Pool.Jobs[0]);

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Terminal viewer can now process
	CPU_Pools = &Pool;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Starting worker threads
	fputs("Launching aligners\n", stdout);
	Pool.threads = (pthread_t*) alloca(Pool.num_threads*sizeof(pthread_t));
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, 1024UL*1024UL);
	for(unsigned int i=0; i<Pool.num_threads; i++) {
#ifdef USE_AFFINITY
		if (data->Affinities) {
			if(pthread_attr_setaffinity_np(&attr, sizeof(Affinity_Mask_t), (cpu_set_t*) &(data->Affinities[i]))) {
				fprintf(stderr, "Unable to set affinity to aligner worker thread %i\n", i);
				goto bail;
			}
		}
#endif
		if (pthread_create(&Pool.threads[i], &attr, (void* (*)(void*)) data->ExtraPtr, (void*) &Pool) != 0) {
			printf("CPU block %s unable to start thread %u\n", __FUNCTION__, i);
			goto bail;
		}
	}
                                                                                                                                        

	for (unsigned int i=1; i<data->nInputBlock; i++) {
		////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Set id and memory parameters
		Pool.Jobs[i].index = i;
		atomic_store(&(Pool.Jobs[i].Treated), 0UL);
		atomic_store(&(Pool.Jobs[i].Stored), 0UL);
		
		// WARNING: Incrementing working threads to make sure we correctly get off the semaphore 
		//          in the case threads go faster.

		pthread_mutex_lock(&(data->Pool->thcount_lock));
		data->Pool->num_threads_working++;
		pthread_mutex_unlock(&(data->Pool->thcount_lock));

		printf("CPU block aligner requesting input buffer %u (of %u) loading...\n", i, data->nInputBlock); fflush(stdout);
		bsem_wait(data->Pool->ToBeAligned_q->has_items);
		
		pthread_mutex_lock(&(data->Pool->ToBeAligned_q->rwmutex));
		DecoderJobElement_t * const DecoderJob = (DecoderJobElement_t *) jobqueue_pull(data->Pool->ToBeAligned_q);
		pthread_mutex_unlock(&(data->Pool->ToBeAligned_q->rwmutex));

		if (DecoderJob == NULL) {
			// FIXME - this scenario seems to cause deadlocks down the road...
			printf("%s: Not enough tags to satisfy %i input buffers\n", __FUNCTION__, data->nInputBlock);
			Pool.threads_keepalive = 0;
			bsem_post(data->Pool->ToBeAligned_q->has_items);
			break;
		}
		
		Pool.Jobs[i].TagsData = DecoderJob->Reader.ReaderJob->Tags;
		Pool.Jobs[i].Anchors = (const DecoderElement_t (*)[16]) DecoderJob->Decoder.DecoderJob->Anchors;
		Pool.Jobs[i].TagCnt = DecoderJob->Reader.ReaderJob->TagCnt;
		Pool.Jobs[i].DecoderJob = DecoderJob;
	
		printf("CPU block aligner got data with %u tags\n", Pool.Jobs[i].TagCnt);
		
		////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Push to queue
		int Test = 0;
		while (!atomic_compare_exchange_strong(&(Pool.job_q.rwmutex), &Test, -1)) {
			Test = 0;
			__asm__ __volatile__ ("pause" ::: "memory");
		}
		pushToQueue(&(Pool.job_q), (Task_t*) &Pool.Jobs[i]);
		atomic_store(&(Pool.job_q.rwmutex), 0);
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// wait for tag aligner threads completion
	fputs("CPU block aligner main thread awaiting jobs completion\n", stdout); fflush(stdout);
	for(unsigned int i=0; i<Pool.num_threads; i++) { 
		void *tres;
		pthread_join(Pool.threads[i], &tres);
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Flush non empty write buffers to host and signal done to the writers
	jobqueue_t * const restrict ToBeWritten_q = &(data->Pool->ToBeWritten_q);
	{
		fputs("CPU block aligner main thread flushing non empty write buffers\n", stdout); fflush(stdout);
		pthread_mutex_lock(&(ToBeWritten_q->rwmutex));
		for (unsigned int i=0; i<data->nOutputBlock; i++) {
			// If this is defined, chances are we already just pushed it out to be written and it is not back yet - no need to push again
			if (Pool.OutputBuffers[i].BelongsTo_q == NULL)
			{
				const register size_t Size = atomic_load(&(Pool.OutputBuffers[i].Stored));
				if (Size > 0) {
					assert(Pool.space_q.has_items != NULL);
					Pool.OutputBuffers[i].BelongsTo_q = &(Pool.space_q);
					Pool.OutputBuffers[i].length = Size;
					assert(Pool.OutputBuffers[i].prev == NULL);
					jobqueue_push(ToBeWritten_q, (job_t*) &Pool.OutputBuffers[i]);
				}
			}
		}
		cpuSpace_t *TriggerDone = calloc(1,sizeof(cpuSpace_t));
		jobqueue_push(ToBeWritten_q, (job_t*) TriggerDone);
		pthread_mutex_unlock(&(ToBeWritten_q->rwmutex));
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Prints out statistics
	{
		printf("buf    direct MM_ovrflw   LargeAlign(rescued)    Alignment(rescued) |     Quick HalfQuick HalfAlign FullAlign BstOfEach   HalfMap  Hopeless   ForSize | Alternate\n"); 
		union Statistics sum;
		__m128i __zero = _mm_setzero_si128();
		for (unsigned int cpu = 0; cpu < Pool.num_threads; cpu++) {
			sum.Vector[0] = __zero;
			sum.Vector[1] = __zero;
			sum.Vector[2] = __zero;
			sum.Vector[3] = __zero;
		}
		for (unsigned int cpu = 0; cpu < Pool.num_threads; cpu++) {
			printf("%3u %9u %9u %9u (%9u) %9u (%9u) | %9u %9u %9u %9u %9u %9u %9u %9u | %9u\n", cpu,
							Pool.stats[cpu].Counter[directhit],
							Pool.stats[cpu].Counter[mismatch_overflow],
							Pool.stats[cpu].Counter[LargeAlignment],
							Pool.stats[cpu].Counter[LargeAlignmentRescued],
							Pool.stats[cpu].Counter[Alignment],
							Pool.stats[cpu].Counter[AlignmentRescued],
							Pool.stats[cpu].Counter[quickalign],
							Pool.stats[cpu].Counter[HalfQuickAlign],
							Pool.stats[cpu].Counter[HalfAlignment],
							Pool.stats[cpu].Counter[FullAlignment],
							Pool.stats[cpu].Counter[BestOfEach],
							Pool.stats[cpu].Counter[HalfMapped],
							Pool.stats[cpu].Counter[hopeless],
							Pool.stats[cpu].Counter[ForCompression],
							Pool.stats[cpu].Counter[Alternate]
						);

			sum.Vector[0] = _mm_add_epi32(sum.Vector[0], Pool.stats[cpu].Vector[0]);
			sum.Vector[1] = _mm_add_epi32(sum.Vector[1], Pool.stats[cpu].Vector[1]);
			sum.Vector[2] = _mm_add_epi32(sum.Vector[2], Pool.stats[cpu].Vector[2]);
			sum.Vector[3] = _mm_add_epi32(sum.Vector[3], Pool.stats[cpu].Vector[3]);
		}
        printf("-----------------------------------------------------------------------------------------------------------------------------------------------------\n"
					 "    %9u %9u %9u (%9u) %9u (%9u) | %9u %9u %9u %9u %9u %9u %9u %9u | %9u\n",
						sum.Counter[directhit],
						sum.Counter[mismatch_overflow],
						sum.Counter[LargeAlignment],
						sum.Counter[LargeAlignmentRescued],
						sum.Counter[Alignment],
						sum.Counter[AlignmentRescued],
						sum.Counter[quickalign],
						sum.Counter[HalfQuickAlign],
						sum.Counter[HalfAlignment],
						sum.Counter[FullAlignment],
						sum.Counter[BestOfEach],
						sum.Counter[HalfMapped],
						sum.Counter[hopeless],
						sum.Counter[ForCompression],
						sum.Counter[Alternate]
					);
		printf("                                                                    ---------------------------------------------------------------------------------\n"
		       "                                     %70u\n",
					sum.Counter[quickalign] + sum.Counter[HalfQuickAlign] + sum.Counter[HalfAlignment] + 
					sum.Counter[FullAlignment] + sum.Counter[BestOfEach] + sum.Counter[hopeless] +
					sum.Counter[ForCompression] + sum.Counter[HalfMapped]);
	}
	err = 0;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Wait until all output memory buffers and terminate trigger are back in the space queue
	printf("Start waiting for output blocks to return\n");
	{
		/* Continuous polling */
		double timeout = 1.0;
		time_t start, end;
		double tpassed = 0.0;
		time (&start);
		while (tpassed < timeout && (Pool.space_q.len < data->nOutputBlock)) {
			time (&end);
			tpassed = difftime(end,start);
		}

		/* Exponential polling */
		long init_nano = 1; /* MUST be above 0 */
		long new_nano;
		double multiplier = 1.01;
		int  max_secs = 20;

		struct timespec polling_interval;
		polling_interval.tv_sec  = 0;
		polling_interval.tv_nsec = init_nano;

		while (Pool.space_q.len < data->nOutputBlock) {
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
		while (Pool.space_q.len < data->nOutputBlock)
		{
			sleep(max_secs);
		}
	}
	printf("Done waiting for output blocks to return\n");

	bail:;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Signal we are no longer active within the Aligner global pool
	pthread_mutex_lock(&(data->Pool->thcount_lock));
	data->Pool->num_threads_alive--;
	pthread_mutex_unlock(&(data->Pool->thcount_lock));
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Clean up and quit
	CPU_Pools = NULL;
	if (err != 0) { printf("CPU main block aligner finish with error!\n"); fflush(stdout);}
	if (OutputMemoryBlocks) munmap(OutputMemoryBlocks, OutBlockSize*data->nOutputBlock);
	
#ifdef BENCHMARKS
	char Line[256];
	fputs("-------------------------- BENCHMARKS ----------------------------------\n",stdout);
	
	for (int i=0; i<gRecvCnt; i++) {
		Printrusage("", &(ReaderBenchmarks[i].usage));
		printf("\nRECEIVER %d:\n", i);
		for (unsigned int j=0; j<gSlaveCnt; j++) {
			char * restrict cptr = Line;
			const tic_t * const restrict ltopics = ProcessBenchmarks[i*gSlaveCnt+j].topics;
			for (int j=0; j<MAX_PROCESS_TOPICS; j++) {
				cptr += snprintf(cptr, 256, "%s=%lu ", ProcessTopics[j], ltopics[j].counter );
			}
			printf("\n\tPROCESSER %d: %s\n", j, Line);
			float ftime[MAX_PROCESS_TOPICS];
			cptr = Line;
			
			*cptr++ = '\t'; *cptr++ = '\t';
			for (unsigned int j=0; j<MAX_PROCESS_TOPICS; j++) {
				cptr += snprintf(cptr, 256, "%s=%8.3f ", ProcessTopics[j], (float) ltopics[j].tics / (float) ltopics[j].counter );
			}
			Printrusage(Line, &(ProcessBenchmarks[i*gSlaveCnt+j].usage));
		}
	}	
	fputs("------------------------------------------------------------------------\n", stdout);
#endif
	printf("CPU main block aligner closing...\n"); fflush(stdout);	
	return (void*) ((size_t) err);
}
//---------------------------------------------------------------

/* vim: tabstop=2 shiftwidth=2
 */
