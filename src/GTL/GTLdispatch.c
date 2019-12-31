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
#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
 
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <pthread.h>
#include <signal.h>
#include <xmmintrin.h>
#include <assert.h>

#include "Genome.h"
#include "virt_chr.h"
#include "GTLdispatch.h"
//include "Aligner/align.h"

#define kMAXPAIRLEN (1024+10)
Genome_t Genome __attribute__((weak));
// VIRTUALCHR virtchr[kMAXCHRcnt] __attribute__((weak));
// const char * const genome_tbl __attribute__((weak));

static const char * const IndexFileSuffix[] = {
	"none",
	"idx",
	"idr"
};

typedef struct GTLFileIndex {
	off_t ReferenceSize;
	struct timespec ReferenceLastModification;
	int PositionType;
	unsigned int nBlocks;
} GTLFileIndex_t;

threadpool_t * createGTLThreadPool(void* (*Fct)(threadpool_t * const restrict), const cpu_set_t * const restrict affinities, const size_t jobvarsize, const int nThreads, const int onHold)
{
	threadpool_t * const restrict thpool = (threadpool_t*) malloc(sizeof(threadpool_t));
	if (thpool == NULL) return NULL;
	
	thpool->threads = (pthread_t*) malloc(nThreads*sizeof(pthread_t));
	if (thpool->threads == NULL) goto err1;
	
	thpool->threads_on_hold   = onHold;
	thpool->threads_keepalive = 1;
	
	thpool->num_threads         = nThreads;
	thpool->num_threads_alive   = 0;
	thpool->num_threads_working = 0;
	
	pthread_mutex_init(&(thpool->thcount_lock), NULL);
	
	if (jobqueue_init(&thpool->jobqueue_p) == -1) goto err2;

	if (jobvarsize > 0 ) {
		if (jobqueue_init(&thpool->donequeue_p) == -1) goto err3;
		
		/* Create free memory jobs and place them into donequeue, EXTRA could be NICE ?*/
		void * restrict Jobs = (void*) malloc((2*nThreads+1)*jobvarsize);
		if (Jobs == NULL) goto err4;
		
		thpool->jobs = (job_t*) Jobs;
		for(int iJob=0; iJob<=(2*nThreads); iJob++) {
			jobqueue_push(&thpool->donequeue_p, (job_t*) Jobs);
			Jobs += jobvarsize;
		}
	}
	else {
		thpool->jobs = NULL;
	}
	
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, (size_t) (2*1024*1024));
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
// 	pthread_attr_setguardsize(&attr, (size_t) (2*4096));
	if (affinities) {
		for (int iThread = 0; iThread<nThreads; iThread++) {
			pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &affinities[iThread]);
			if (pthread_create(&(thpool->threads[iThread]), &attr, (void* (*)(void*)) Fct, thpool) != 0) {
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
	return thpool;
	
	err5:
		jobqueue_destroy(&thpool->jobqueue_p);
	err4:
		jobqueue_destroy(&thpool->donequeue_p);
	err3:
		pthread_mutex_destroy(&(thpool->thcount_lock));
	err2:
		free(thpool->threads);
	err1:
		free(thpool);
	return NULL;
}
//---------------------------------------------------------------

void destroyGTLThreadPool(threadpool_t * restrict thpool)
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
	while (tpassed < TIMEOUT && thpool->num_threads_alive){
		bsem_post_all(thpool->jobqueue_p.has_items);
		time (&end);
		tpassed = difftime(end,start);
	}

	/* Poll remaining threads */
	while (thpool->num_threads_alive){
		bsem_post_all(thpool->jobqueue_p.has_items);
		sleep(1);
	}

	/* Job queue cleanup */
	jobqueue_destroy(&thpool->jobqueue_p);
	if (thpool->jobs) {
		jobqueue_destroy(&thpool->donequeue_p);
		free(thpool->jobs);
	}

	free(thpool->threads);
	free(thpool);
	thpool = NULL;
}
//---------------------------------------------------------------

GTLList_t * createGTLList(const char * const restrict * const Files, const unsigned int nFiles)
{
	GTLList_t * const restrict list = (GTLList_t*) malloc(nFiles*sizeof(GTLList_t));
	if (list == NULL) return NULL;
	
	list->fds = (int*) calloc(nFiles,sizeof(int));
	if (list->fds == NULL) goto err2;
	
	for (unsigned int iFile=0; iFile<nFiles; iFile++) {
		const int itmp = open(Files[iFile], O_RDONLY);
		if (itmp == -1) { 
			fprintf(stderr, "Unable to open file '%s'\n", Files[iFile]);
			goto err4;
		}
		list->fds[iFile] = itmp;
	}
	
	list->FileNames = Files;
	list->nFiles = nFiles;
	list->PositionType = None;
	list->nBlocks = 0U;
	list->Blocks = NULL;
	
	return list;
	err4:
		for (unsigned int iFile=0; iFile<nFiles; iFile++) {
			if (list->fds[iFile] > 0 ) close(list->fds[iFile]);
		}
		free(list->fds);
	err2:
		free(list);
		return NULL;
}
//---------------------------------------------------------------

void destroyGTLList(GTLList_t * GL) {
	if (GL) {
		if (GL->fds) free( GL->fds);
		if (GL->Blocks) {
			for (size_t iBlock=0UL; iBlock<GL->nBlocks; iBlock++) {
				pthread_mutex_destroy(&(GL->Blocks[iBlock].decompressing));
				if (GL->Blocks[iBlock].whatever) {
					freeBlock((TLBDATA *) GL->Blocks[iBlock].whatever);
					free(GL->Blocks[iBlock].whatever);
				}
			}
			free(GL->Blocks);
		}
		free(GL);
	}
	GL = NULL;
}
//---------------------------------------------------------------

static int CompareGTLBlockInfo(const void* A, const void* B)
{
	const GTLBlockInfo_t *const restrict ta = A, * const restrict tb = B;
	if (ta->thd.minPos < tb->thd.minPos)
		return -1;
	else if (ta->thd.minPos > tb->thd.minPos) 
		return 1;
	else
		return 0;
}
//---------------------------------------------------------------

static int testGTLFileIndex(const char * const restrict FileName, const enum GTLHeaderExtremaType Order)
{
	char IndexFileName[256];
	snprintf(IndexFileName, 256, "%s.%s", FileName, IndexFileSuffix[(int) Order]);
	fprintf(stderr, "Testing existance of index file %s, ", IndexFileName);
	const int IndexFile = open(IndexFileName, O_NOATIME | O_RDONLY);
	if (IndexFile > 0) {
		GTLFileIndex_t Src;
		again:;
		const ssize_t res = read(IndexFile, &Src, sizeof(GTLFileIndex_t));
		
		if ( res == sizeof(GTLFileIndex_t)) {
			struct stat st;
			if (stat(FileName, &st) == 0) {
				/* Check size and last modification time */
				if (st.st_size == Src.ReferenceSize && st.st_mtim.tv_nsec == Src.ReferenceLastModification.tv_nsec 
				    && st.st_mtim.tv_sec == Src.ReferenceLastModification.tv_sec && (int) Order == Src.PositionType) {
					close(IndexFile);
					if (Src.nBlocks < INT32_MAX) {
						fputs("OK\n", stderr);
						return Src.nBlocks;
					}
					else {
						fputs("Number of blocks exceed int 32\n", stderr);
						return -1;
					}
				}
			}
			else {
				perror("stat error");
			}
		}
		else if (res == EINTR) 
			goto again;
		fprintf(stderr, "bad or no header found\n");
		close(IndexFile);
		return -1;
	}
	fprintf(stderr, "not found\n");
	return -1;
}
//---------------------------------------------------------------

static int readGTLFileIndex(const char * const restrict FileName, const int fd, GTLBlockInfo_t * restrict BLINFO,
                            const enum GTLHeaderExtremaType Order)
{
	char IndexFileName[256];
	snprintf(IndexFileName, 256, "%s.%s", FileName, IndexFileSuffix[(int) Order]);
	fprintf(stderr, "Reading index for %s, ", IndexFileName);
	FILE* IndexFile = fopen(IndexFileName, "rb");
	if (IndexFile != NULL) {
		GTLFileIndex_t Src;
		again:;
		size_t res = fread(&Src, sizeof(GTLFileIndex_t), 1, IndexFile);
		if ( res == 1 ) {
			struct stat st;
			if (stat(FileName, &st) == 0) {
				/* Check size and last modification time */
				if (st.st_size == Src.ReferenceSize && st.st_mtim.tv_nsec == Src.ReferenceLastModification.tv_nsec 
				    && st.st_mtim.tv_sec == Src.ReferenceLastModification.tv_sec && (int) Order == Src.PositionType) {
					int OK = 1;
					for (unsigned int iBlock=0; iBlock<Src.nBlocks; iBlock++) {
						OK &= fread(&(BLINFO->thd), sizeof(TLBHEADER), 1, IndexFile) ==  1;
						OK &= fread(&(BLINFO->offset), sizeof(off_t), 1, IndexFile) ==  1;
						OK &= fread(&(BLINFO->PointToFraction[1]), sizeof(DirectLink_t), GTL_INDEX_SPLIT-1, IndexFile) == (GTL_INDEX_SPLIT-1);
						BLINFO->PointToFraction[0].ID       = 0U;
						BLINFO->PointToFraction[0].Left     = 0U;
						BLINFO->PointToFraction[0].Right    = 0U;
						BLINFO->PointToFraction[0].Header   = 0U;
						BLINFO->PointToFraction[0].Cigar    = 0U;
						BLINFO->PointToFraction[0].Mismatch = 0U;
						BLINFO->PointToFraction[0].qs       = 0U;
						BLINFO->whatever                    = NULL;
						BLINFO->fd                          = fd;
						BLINFO->usingIt                     = 0;
						pthread_mutex_init(&(BLINFO->decompressing), NULL);
						BLINFO++;
					}
					fclose(IndexFile);
					if (OK && Src.nBlocks < INT32_MAX) {
						fprintf(stderr, "read %u blocks\n", Src.nBlocks);
						return Src.nBlocks;
					}
					else {
						fputs("read error or number of blocks exceed int 32\n", stderr);
						return -1;
					}
				}
			}
		}
		else if (res == EINTR) 
			goto again;
		fprintf(stderr, "bad or no header found\n");
		fclose(IndexFile);
	}
	fprintf(stderr, "not found\n");
	return -1;
}
//---------------------------------------------------------------

int indexGTLList(GTLList_t * const restrict list, const enum GTLHeaderExtremaType Order, 
                 const _Bool DoNotUseExisting, const unsigned char provenance)
{
	TLBHEADER BlockHeader;
	struct stat st;
	int err=1;
	const unsigned int nFiles = list->nFiles;
	unsigned char IndexFound[nFiles];
	
	size_t TotalBlocks = 0UL;
	
	if (nFiles == 0) goto bail;

	//////////////////////////////////////////////////////////////////////////////////////////////////
	// Counting blocks
	int iFile = 0;
	unsigned char HasScanned = 0;

	while (iFile < nFiles) {
		/* Try to see if we have an index file */
		const int nBlocks = (DoNotUseExisting) ? -1 : testGTLFileIndex(list->FileNames[iFile], Order);
		if (nBlocks < 0) {
			/* Perform the scan */
			HasScanned = 1;
			IndexFound[iFile] = 0;
			const int fd = list->fds[iFile];
			if (fstat(fd, &st) != 0) {
				perror("fstat");
			}
			off_t FileOffset = (off_t) 0;
			while (1) {
				if (FileOffset >= st.st_size) break;
				
				//////////////////////////////////////////////////////////////////////////////////////////////
				// Read GTL Block headers
				const ssize_t res = pread(fd, &BlockHeader, sizeof(TLBHEADER), FileOffset);
				
				if (res == -1) {
					fprintf(stderr, "Failed to read header block in %s @ %lu\n", list->FileNames[iFile], FileOffset);
					perror("pread");
					goto bail;
				}
							
				/* Verify checksum */
				if (res != sizeof(TLBHEADER)) {
					fputs("Wrong header size returned from reading\n", stderr);
					goto bail;
				}
				const unsigned char hcs = simple8bitCS((unsigned char *) &BlockHeader, sizeof(TLBHEADER));
				if (hcs != 0) {
					fprintf(stderr,"Checksum mismatch: %02x %02x\n",BlockHeader.headerCS,hcs);
					goto bail;
				}
				
				/* Do we want this block ? */
				TotalBlocks++;
					
				/* Move forward in file */
				FileOffset += BlockHeader.blockLength;
			}
		}
		else {
			IndexFound[iFile] = 1;
			TotalBlocks += (size_t) nBlocks;
		}
		iFile++;
	}
	
	list->PositionType = (HasScanned) ? FirstTag : Range;
	list->nBlocks = TotalBlocks;
	GTLBlockInfo_t * restrict BLInfo = (GTLBlockInfo_t*) calloc(TotalBlocks, sizeof(GTLBlockInfo_t));
	if (BLInfo == NULL) goto bail;
	
	list->Blocks = BLInfo;
	
	//////////////////////////////////////////////////////////////////////////////////////////////////
	// Storing blocks
	iFile = 0;
	while (iFile < nFiles) {
		/* Try to see if we have an index file */
		const int fd = list->fds[iFile];
		const int nBlocks = (IndexFound[iFile]) ? readGTLFileIndex(list->FileNames[iFile], fd, BLInfo, Order) : -1;
		if (nBlocks < 0) {
			/* Perform the scan */
			if (fstat(fd, &st) != 0) {
				perror("fstat");
			}
			off_t FileOffset = (off_t) 0;
			while (1) {
				if (FileOffset >= st.st_size) break;
				
				//////////////////////////////////////////////////////////////////////////////////////////////
				// Read GTL Block headers
				const ssize_t res = pread(fd, &(BLInfo->thd), sizeof(TLBHEADER), FileOffset);
				
				if (res == -1) {
					fprintf(stderr, "Failed to read header block in %s @ %lu\n", list->FileNames[iFile], FileOffset);
					perror("pread");
					goto bail;
				}
				
				/* Do we want this block ? */
				BLInfo->whatever = NULL;
				BLInfo->offset = FileOffset;
				BLInfo->fd = fd;
				BLInfo->usingIt = 0;
				BLInfo->provenance = provenance;
				pthread_mutex_init(&(BLInfo->decompressing),0);
				FileOffset += BLInfo->thd.blockLength;
				BLInfo++;
			}
		}
		else {
			for (int i=0; i<nBlocks; i++) {
				BLInfo->provenance = provenance;
				BLInfo++; 
			}
		}
		iFile++;
	}
	
	err = 0;
	bail:;
	
	return err;
}
//---------------------------------------------------------------

int writeGTLFileIndex(const GTLList_t * const restrict list, const char * const restrict Directory) 
{
	char Buffer[512];
	FILE * out;
	char * FileNamePtr;
	GTLFileIndex_t Header;
	struct stat st;
	
	const unsigned int nFiles = list->nFiles;
	
	int DirLength; 
	if (Directory) {
		DirLength = snprintf(Buffer, 512, "%s/", Directory);
		if (DirLength >= 512) {
			fprintf(stderr, "The specified directory path is too long (512 max): %s\n", Directory);
			return -1;
		}
		FileNamePtr = Buffer + DirLength;
	}
	else { 
		FileNamePtr = Buffer;
		DirLength = 0;
	}
	Header.nBlocks = 0;
	Header.PositionType = list->PositionType;
	int err = 0;
	for (unsigned int iFile=0; iFile<nFiles; iFile++) {
		if (fstat(list->fds[iFile], &st) == 0) {
			Header.ReferenceLastModification = st.st_mtim;
			Header.ReferenceSize = st.st_size;
			const char * baseName;
			if (Directory) {
				const int FileNameLength = (int) strlen(list->FileNames[iFile]);
				baseName = &list->FileNames[iFile][FileNameLength-1];
				while ((uintptr_t) baseName > (uintptr_t) list->FileNames[iFile] && *baseName != '/') baseName--;
				if (*baseName == '/') baseName++;
			}
			else 
				baseName = list->FileNames[iFile];
			
			const int TotalLength = DirLength + snprintf(FileNamePtr, 512, "%s.%s", baseName, IndexFileSuffix[Header.PositionType]);
			if (TotalLength < 512) {
				out = fopen(Buffer, "rb");
				if ( out == NULL) {
					out = fopen(Buffer, "wb+");
					if (out) {
						int OK = fwrite(&Header, sizeof(GTLFileIndex_t), 1, out) == 1;
						const int WantedFD = list->fds[iFile];
						const GTLBlockInfo_t * const limit = list->Blocks + list->nBlocks;
						GTLBlockInfo_t * restrict CurrentBlock = list->Blocks;
						unsigned int count = 0U;
						while(CurrentBlock < limit) {
							if (CurrentBlock->fd == WantedFD) {
								OK &= fwrite(&(CurrentBlock->thd), sizeof(TLBHEADER), 1, out) ==  1;
								OK &= fwrite(&(CurrentBlock->offset), sizeof(off_t), 1, out) == 1;
								OK &= fwrite(&(CurrentBlock->PointToFraction[0]), sizeof(DirectLink_t), (GTL_INDEX_SPLIT-1), out) == (GTL_INDEX_SPLIT-1);
								count++;
							}
							++CurrentBlock;
						}
						rewind(out);
						Header.nBlocks = count;
						OK &= fwrite(&Header, sizeof(GTLFileIndex_t), 1, out) == 1;
						fclose(out);
						if (OK) {
							fprintf(stderr, "Creating index file %s\n", Buffer);
							continue;
						}
						else {
							fprintf(stderr, "Error writing to index file %s\n", Buffer);
						}
					}
					else {
						fprintf(stderr, "Cannot create index file %s\n", Buffer);
					}
				}
				else {
					fprintf(stderr, "Index file %s already exists.\n", Buffer);
					fclose(out);
				}
			}
			else {
				fprintf(stderr, "Cannot create index file for %s, overall path too long\n", list->FileNames[iFile]);
			}
		}
		else {
			fprintf(stderr, "Cannot stat file %s, indexing impossible for this file\n", list->FileNames[iFile]);
		}
		
		err = -1;
	}
	return err;
}
//---------------------------------------------------------------

/* Extract the Tags - COULD BE OPTIMIZED */  
void 
ExtractData(decodedPair_t * const restrict dpp, const unsigned char * restrict diff,
            const unsigned char * restrict cigar, GTLRawData_t * const restrict * data)
{
	unsigned char tmpStr[kMaxReadLen];
	
	GTLRawData_t * restrict GD = data[0];
	GD->TagLength = (unsigned char) dpp->taglen1;
	assert(dpp->taglen1<kMaxReadLen);
	memcpy(GD->Tag,&(Genome.table[dpp->genomepos]),dpp->taglen1);
	const unsigned char * restrict GenomePtr = &(Genome.table[dpp->genomepos]);
	GD->revNmer = dpp->reverseTAG1;
	int first = 1;
	
	SecondTag:
	
	{
		memcpy(GD->Tag, GenomePtr, GD->TagLength);
		GD->Tag[GD->TagLength] = '\0';
		GD->AlignmentRange[0] = 0;
		unsigned int cigarlength = 0U;
		if (cigar != NULL && *cigar != kCIGARTERMINATOR) {
			int r=0, a=0;
// 			printf("CIG: ");
			while (cigar[cigarlength] != kCIGARTERMINATOR)
			{
				register int pos = cigar[cigarlength];
				const char code  = cigar[cigarlength+1];
				GD->Cigar[cigarlength  ] = cigar[cigarlength  ];
				GD->Cigar[cigarlength+1] = cigar[cigarlength+1];
				
// 				printf("%u%c ", (unsigned int) cigar[cigarlength], (int) cigar[cigarlength+1]); 
				cigarlength += 2;

				switch(code)
				{
					case 'S': // same as 'M'
					case 'M': assert(a+pos<=kMaxReadLen); while(pos-- >0) { GD->Tag[a++] = GenomePtr[r]; r++; } break;
					case 'D':                             while(pos-- >0) { r++;                           } break;
					case 'I': assert(a+pos<=kMaxReadLen); while(pos-- >0) { GD->Tag[a++] = 'N';            } break;
				}
			}
// 			printf("\n");
			GD->Tag[a] = '\0';
			assert(r >= 1);
			GD->AlignmentRange[1] = r - 1;
		}
		else {
			GD->AlignmentRange[1] = GD->TagLength - 1;
		}
		if (cigar != NULL && cigar[cigarlength] == kCIGARTERMINATOR) GD->Cigar[cigarlength++] = kCIGARTERMINATOR;
		GD->CigarLength = (unsigned char) cigarlength;
		cigar += cigarlength;
		
		unsigned int count = 0;
		if (diff) {
			while (diff[count] != kMMTERMINATOR)
			{
				if (diff[count] == kMMSOFTCLIP)
				{
					// grab the nucleotides clipped in a tmp string
					int i = 0;
					count++;
					while (diff[count] != kMMTERMINATOR) { tmpStr[i++] = diff[count++] ; }
					// copy the clipped nucleotides at the end of the tag.
					char * const copyto = (char*) ((GD->revNmer) ? &(GD->Tag[0]) : &(GD->Tag[GD->TagLength-i]));
					memcpy(copyto,tmpStr,i);
					break; // kMMSOFTCLIP is the last thing in MMstring. only one occurence possible. we are done.
				}
				const int pos = diff[count];
				assert(pos<GD->TagLength);
				const unsigned char nt = diff[count+1];
				GD->Tag[pos] = nt;
				count += 2;
			}
	// 		printf("\n");
			count++; // skip the kMMTERMINATOR
		}
		
		GD->MismatchLength = (unsigned char) count;
		diff += count;
	}

	if (first) {
		first = 0;
		GD = data[1];
		GD->TagLength = (unsigned char) dpp->taglen2;
		assert(dpp->taglen2<kMaxReadLen);
		memcpy(GD->Tag,&(Genome.table[dpp->genomepos]),dpp->taglen2);
		GenomePtr = &(Genome.table[dpp->genomepos + dpp->delta]);
		GD->revNmer = dpp->reverseTAG2;
		goto SecondTag;
	}
	
}
//---------------------------------------------------------------

void 
ExtractDataNasSoftClip(decodedPair_t * const restrict dpp, const unsigned char * restrict diff,
	               const unsigned char * restrict cigar, GTLRawData_t * const restrict * data)
{
	GTLRawData_t * restrict GD = data[0];
	GD->TagLength = (unsigned char) dpp->taglen1;
	assert(dpp->taglen1<kMaxReadLen);
	memcpy(GD->Tag,&(Genome.table[dpp->genomepos]),dpp->taglen1);
	const unsigned char * restrict GenomePtr = &(Genome.table[dpp->genomepos]);
	GD->revNmer = dpp->reverseTAG1;
	int first = 1;
	
	SecondTag:
	
	{
		memcpy(GD->Tag, GenomePtr, GD->TagLength);
		GD->Tag[GD->TagLength] = '\0';
		GD->AlignmentRange[0] = 0;
		unsigned int cigarlength = 0U;
		if (cigar != NULL && *cigar != kCIGARTERMINATOR) {
			int r=0, a=0;
// 			printf("CIG: ");
			while (cigar[cigarlength] != kCIGARTERMINATOR)
			{
				register int pos = cigar[cigarlength];
				const char code  = cigar[cigarlength+1];
				GD->Cigar[cigarlength  ] = cigar[cigarlength  ];
				GD->Cigar[cigarlength+1] = cigar[cigarlength+1];
				
// 				printf("%u%c ", (unsigned int) cigar[cigarlength], (int) cigar[cigarlength+1]); 
				cigarlength += 2;

				switch(code)
				{
					case 'S': assert(a+pos<=kMaxReadLen); while(pos-- >0) { GD->Tag[a++] = 'N'; r++; } break;
					case 'M': assert(a+pos<=kMaxReadLen); while(pos-- >0) { GD->Tag[a++] = GenomePtr[r]; r++; } break;
					case 'D':                             while(pos-- >0) { r++;                           } break;
					case 'I': assert(a+pos<=kMaxReadLen); while(pos-- >0) { GD->Tag[a++] = 'N';            } break;
				}
			}
// 			printf("\n");
			GD->Tag[a] = '\0';
			assert(r >= 1);
			GD->AlignmentRange[1] = r - 1;
		}
		else {
			GD->AlignmentRange[1] = GD->TagLength - 1;
		}
		if (cigar != NULL && cigar[cigarlength] == kCIGARTERMINATOR) GD->Cigar[cigarlength++] = kCIGARTERMINATOR;
		GD->CigarLength = (unsigned char) cigarlength;
		cigar += cigarlength;
		
		unsigned int count = 0;
		if (diff) {
			while (diff[count] != kMMTERMINATOR)
			{
				if (diff[count] == kMMSOFTCLIP)
				{
					// grab the nucleotides clipped in a tmp string
					count++;
					while (diff[count] != kMMTERMINATOR) count++;
					break; // kMMSOFTCLIP is the last thing in MMstring. only one occurence possible. we are done.
				}
				const int pos = diff[count];
				assert(pos<GD->TagLength);
				const unsigned char nt = diff[count+1];
				GD->Tag[pos] = nt;
				count += 2;
			}
	// 		printf("\n");
			count++; // skip the kMMTERMINATOR
		}
		
		GD->MismatchLength = (unsigned char) count;
		diff += count;
	}

	if (first) {
		first = 0;
		GD = data[1];
		GD->TagLength = (unsigned char) dpp->taglen2;
		assert(dpp->taglen2<kMaxReadLen);
		memcpy(GD->Tag,&(Genome.table[dpp->genomepos]),dpp->taglen2);
		GenomePtr = &(Genome.table[dpp->genomepos + dpp->delta]);
		GD->revNmer = dpp->reverseTAG2;
		goto SecondTag;
	}
	
}
//---------------------------------------------------------------

// This is not really recommended, but much easier
int LinksCmp(const void * a, const void * b)
{
	DirectLink_t *const restrict ta = (DirectLink_t *) a, * const restrict tb = (DirectLink_t*) b;
	if (ta->Right < tb->Right) {
		return -1;
	}
	else if (ta->Right > tb->Right) {
		tb->ID = ta->ID;
		const __m128i tmp = _mm_load_si128((__m128i*) &(ta->Header));
		_mm_store_si128((__m128i*) &(tb->Header), tmp);
		return 1;
	}
	else
		return 0;
}
//---------------------------------------------------------------

static void* getTagRange(threadpool_t * const restrict thpool)
{
	TLBDATA td;
	GTLRawData_t lND[2]; 
	GTLRawData_t * ND[2] = { &lND[0], &lND[1] };
	DirectLink_t * restrict Links;
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Allocate memory
	allocBlock(&td, true);
	Links = (DirectLink_t*) malloc(kBlockMaxCount*sizeof(DirectLink_t));
	if (Links == NULL) return (void*) 1;
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Mark thread as alive (initialized)
	pthread_mutex_lock(&thpool->thcount_lock);
	thpool->num_threads_alive++;
	pthread_mutex_unlock(&thpool->thcount_lock);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Handle jobqueue
	while(thpool->threads_keepalive){
		bsem_wait(thpool->jobqueue_p.has_items);

		if (thpool->threads_keepalive){
			/* Signal I am working */
			pthread_mutex_lock(&thpool->thcount_lock);
			thpool->num_threads_working++;
			pthread_mutex_unlock(&thpool->thcount_lock);
			
			/* Read job from queue and execute it */
			pthread_mutex_lock(&thpool->jobqueue_p.rwmutex);
			GTLBlockInfo_t * const restrict task = (GTLBlockInfo_t*) jobqueue_pull(&thpool->jobqueue_p);
			pthread_mutex_unlock(&thpool->jobqueue_p.rwmutex);
			
			/*************************************************************************/
			/*                             START THE TASK                            */
			/*************************************************************************/
			if (task) {
				unsigned int MinimumGenomePos = -1;
				unsigned int MaximumGenomePos = 0U;
				
				/* Read and decompress the block */
				if (preadDecompressBlock(task->fd, task->offset + sizeof(TLBHEADER), &(task->thd), &td, 0)) {
						fputs("Error decompressing block\n", stderr);
						exit(1);
				}

				if ((td.header.flags_7_0 & kTLBHflagPairedReads)) {
					unsigned char * restrict diff = ((task->thd.flags_7_0 & 0xf0) == kTLBHflagNMismatchBlock) ? NULL : td.diff[0];
					decodedPair_t dp;
					unsigned int mmpos = 0U;
					unsigned int cigarpos = 0U;
					const uintptr_t Header0 = (uintptr_t) td.hdr;
					const uintptr_t qs0 = (uintptr_t) td.qs;
					
					for (int i=0; i<td.cnt; i++) {
						/* Basic Decode pair */
						{
							if (td.header.flags_7_0 & kTLBHflagFixedLength) {
								dp.taglen1 = td.lengthInfo;
								dp.taglen2 = td.lengthInfo;
							}
							else {
								dp.taglen1 = td.len[2*i];
								dp.taglen2 = td.len[2*i+1];
							}
							Links[i].Header   = (unsigned int) ( (uintptr_t) td.hdr - Header0);
							Links[i].Mismatch = mmpos;
							Links[i].Cigar    = cigarpos;
							Links[i].qs       = (unsigned int) ( (uintptr_t) td.qs - qs0);
							
							char *next = strchr(td.hdr, '\t');
							td.hdr = next + 1;
							next = strchr(td.hdr, '\n');
							td.hdr = next + 1;
							dp.ordinal = 0;
							for (unsigned int j = 0; j < 8; j++)
								dp.ordinal = (dp.ordinal << 8) | (unsigned long) *((unsigned char *) td.hdr++);

							dp.reverseTAG1 = (td.ppd[0][i].tag2offset & k32bREVCOMP_TAG1_BIT) ? 1 : 0;
							dp.reverseTAG2 = (td.ppd[0][i].tag2offset & k32bREVCOMP_TAG2_BIT) ? 1 : 0;

							dp.delta = (td.ppd[0][i].tag2offset & k32bOFFSET_VALUE_MASK);
							if (td.ppd[0][i].tag2offset & k32bNEGATIVE_OFFSET_BIT)
								dp.delta = -dp.delta;

							dp.genomepos = (Genome.virtchr[td.header.chr].chr << 28) + Genome.virtchr[td.header.chr].offset + td.ppd[0][i].tag1pos;
						}
						
						ExtractData(&dp, (diff != NULL) ? &diff[mmpos] : NULL, td.cigar[0] + cigarpos, ND );
						lND[0].Location = td.ppd[0][i].tag1pos;
						lND[1].Location = td.ppd[0][i].tag1pos + dp.delta; 
						unsigned int cmin = -1;
						unsigned int cmax = 0U;
						for (int iTag=0; iTag<2; iTag++) {
							GTLRawData_t * restrict GD = ND[iTag];
							
							unsigned int TagGenLen   = GD->AlignmentRange[1] - GD->AlignmentRange[0];
							unsigned int TagGenBegin = GD->Location;
							unsigned int TagGenEnd   = TagGenBegin + TagGenLen;
							
							if (TagGenBegin < MinimumGenomePos) MinimumGenomePos = TagGenBegin;
							if (TagGenEnd > MaximumGenomePos) MaximumGenomePos = TagGenEnd;
							
							if (TagGenBegin < cmin) cmin = TagGenBegin;
							if (TagGenEnd > cmax) cmax = TagGenEnd;
							
							/* Move pointers */
							td.qs    += GD->TagLength; /* might not be necessary */
							cigarpos += GD->CigarLength;
							mmpos    += GD->MismatchLength;
						}
						Links[i].ID    = (unsigned int) i;
						Links[i].Left  = cmin;
						Links[i].Right = cmax;
					}
					
					qsort(Links, (size_t) td.cnt, sizeof(DirectLink_t), LinksCmp);
// 					{
// 						FILE *tmp = fopen("sort0", "w");
// 						for (int i=0; i<td.cnt; i++) {
// 							fprintf(tmp, "%u\t%u\n", Links[i].Left, Links[i].Right);
// 						}
// 						fclose(tmp);
// 					}
					{
						unsigned int MostLeft = MaximumGenomePos;
						assert(td.cnt > 0);
						for (int i=td.cnt-1; i>=0; i--) {
							if (Links[i].Left < MostLeft) {
								MostLeft = Links[i].Left;
							}
							else {
								Links[i].Left = MostLeft;
							}
						}
					}
// 					{
// 						FILE *tmp = fopen("sort1", "w");
// 						for (int i=0; i<td.cnt; i++) {
// 							fprintf(tmp, "%u\t%u\n", Links[i].Left, Links[i].Right);
// 						}
// 						fclose(tmp);
// 					}
					
					const unsigned int range = MaximumGenomePos - MinimumGenomePos;
					const unsigned int step = 1U + (range/GTL_INDEX_SPLIT);
					unsigned int limit = MinimumGenomePos + step;
					const DirectLink_t * restrict LinksPtr = Links;
					for (size_t j=0; j<(GTL_INDEX_SPLIT-1);j++) { 
						while (limit < MaximumGenomePos && LinksPtr->Right < limit) LinksPtr++;
						task->PointToFraction[j] = *LinksPtr;
						assert(LinksPtr->ID <= kBlockMaxCount);
						limit += step;
					}
					
					/*{
						FILE *tmp = fopen("sort2", "w");
						for (int i=0; i<(GTL_INDEX_SPLIT-1); i++) {
							fprintf(tmp, "%u\t%u\t%u\t%u\t%u\t%u\t%u\n",
							        task->PointToFraction[i].Left,
							        task->PointToFraction[i].Right,
											task->PointToFraction[i].Header,
											task->PointToFraction[i].Mismatch,
											task->PointToFraction[i].Cigar,
											task->PointToFraction[i].qs,
							        task->PointToFraction[i].ID
 										);
						}
						fclose(tmp);
					}*/
					
				}
				
				task->thd.minPos = MinimumGenomePos;
				task->thd.maxPos = MaximumGenomePos;
				task->whatever = NULL;
				memset(&(task->decompressing), 0, sizeof(pthread_mutex_t));
			}
			
			/* Signal I am no longer working */
			pthread_mutex_lock(&thpool->thcount_lock);
			thpool->num_threads_working--;
			pthread_mutex_unlock(&thpool->thcount_lock);
		}
	}

	/* Signal I am done */
	pthread_mutex_lock(&thpool->thcount_lock);
	thpool->num_threads_alive--;
	pthread_mutex_unlock(&thpool->thcount_lock);
 
	freeBlock(&td);
	free(Links);
	
	return (void*) 0;
}
//---------------------------------------------------------------

int getRealBlockMinMaxPosition(GTLList_t * const restrict list, const int nThreads,
															 const unsigned char DoSort)
{
	//////////////////////////////////////////////////////////////////////////////////////////////////
	// Create a thread pool 
	threadpool_t * const restrict thpool = createGTLThreadPool( getTagRange, NULL, 0, nThreads, 0);
	if (thpool == NULL) return -1;
	
	//////////////////////////////////////////////////////////////////////////////////////////////////
	// dispatch blocks
	const size_t nBlocks = list->nBlocks;
	for (size_t iBlock=0; iBlock<nBlocks; iBlock++) {
		jobqueue_push(&(thpool->jobqueue_p), (job_t*) &(list->Blocks[iBlock]));
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////
	// Wait for pool to terminate
	thpool_wait(thpool);
	
	//////////////////////////////////////////////////////////////////////////////////////////////////
	// Destroy a thread pool 
	destroyGTLThreadPool(thpool);
	
	if (DoSort) {
		qsort(list->Blocks, list->nBlocks, sizeof(GTLBlockInfo_t), CompareGTLBlockInfo);
	}
	
	list->PositionType = Range;
	
	return 0;
	
}
//---------------------------------------------------------------

void sortGTLlist(GTLList_t * const restrict list)
{
	qsort(list->Blocks, list->nBlocks, sizeof(GTLBlockInfo_t), CompareGTLBlockInfo);
}
//---------------------------------------------------------------

unsigned int dispatchGTLBlockRange(const GTLList_t * const restrict list, threadpool_t * const restrict thpool,
                                   const unsigned int chr, const unsigned int GenomeStart, unsigned int GenomeEnd,
                                   const unsigned char flags, const _Bool EndingJob,
                                   void * const restrict params)
{
	//////////////////////////////////////////////////////////////////////////////////////////////////
	// Loop on blocks
	unsigned int TotalJobs = 0U;
	const GTLBlockInfo_t * const restrict Limit = list->Blocks + list->nBlocks;
	GTLBlockInfo_t * restrict GTLBlock = list->Blocks;
	
	GTLBlockDispatchRange_t * restrict work;
	
	{
		//////////////////////////////////////////////////////////////////////////////////////////////
		// Aquire empty memory slot
		do {
			bsem_wait(thpool->donequeue_p.has_items);
			pthread_mutex_lock(&(thpool->donequeue_p.rwmutex));
			work = (GTLBlockDispatchRange_t *) jobqueue_pull(&(thpool->donequeue_p));
			pthread_mutex_unlock(&(thpool->donequeue_p.rwmutex));
		} while (work == NULL);
		
		if (GenomeEnd <= GenomeStart) GenomeEnd = GenomeStart+1;
 	
		KeepOn: ;
		/* Do we want this block ? */
		if (GTLBlock->thd.chr == chr && GTLBlock->thd.minPos <= GenomeEnd && GTLBlock->thd.maxPos >= GenomeStart) {
			/*if (BlockHeader->flags_7_0 & flags)*/ {
				work->Block = GTLBlock;
				work->GenomeStart = GenomeStart;
				work->GenomeEnd = GenomeEnd;
				work->TotalJobs = 0U;
				work->params = params;

				pthread_mutex_lock(&(thpool->jobqueue_p.rwmutex));
				jobqueue_push(&(thpool->jobqueue_p), (job_t*) work);
				pthread_mutex_unlock(&(thpool->jobqueue_p.rwmutex));
				TotalJobs++;
				
				do {
					bsem_wait(thpool->donequeue_p.has_items);
					pthread_mutex_lock(&(thpool->donequeue_p.rwmutex));
					work = (GTLBlockDispatchRange_t *) jobqueue_pull(&(thpool->donequeue_p));
					pthread_mutex_unlock(&(thpool->donequeue_p.rwmutex));
				} while (work == NULL);
			}
		}
		if (++GTLBlock < Limit) goto KeepOn;
	}
	
	jobqueue_t * restrict queue; 
	
	if (TotalJobs && EndingJob) {
		/* Close this call by adding an ending job */
// 		printf("Closing search for %u--%u with %u tasks\n", GenomeStart, GenomeEnd, TotalJobs);
		work->TotalJobs = TotalJobs;
		work->params = params;
		work->GenomeStart = GenomeStart;
		work->GenomeEnd = GenomeEnd;
		work->Block = NULL;
		queue  = &(thpool->jobqueue_p);
	}
	else {/* Hey Dude. do not forget to put back the storage into the queue, dummy !!! */
		queue = &(thpool->donequeue_p);
	}
	pthread_mutex_lock(&(queue->rwmutex));
	jobqueue_push(queue, (job_t*) work);
	pthread_mutex_unlock(&(queue->rwmutex));
		
	return TotalJobs;
}
//---------------------------------------------------------------

int dispatchGTLBlockFile(const GTLList_t * const restrict list, threadpool_t * const restrict thpool,
												 const unsigned char flags, const _Bool EndingJob, const unsigned int EndingJobEvery,
												 void * const restrict params)
{
	GTLBlockDispatchFile_t * work;
	struct stat st;
	int err = 1;
	
	/* Check there is enough jobs available in threadpool to reach EndingJobEvery */
	if (EndingJob && (thpool->donequeue_p.len < EndingJobEvery)) {
		fprintf(stderr, "Stall in perspective, there is not enough available jobs to reach the EndingJobEvery argument!\n");
		goto bail;
	}
	
	unsigned int iFile = 0U;
	const unsigned int nFiles = list->nFiles;
	GTLBlockDispatchFile_t * * BlockList = NULL;
	sem_t * restrict semaphore = NULL;
	
	
	while (iFile < nFiles) {
		/* Try to see if we have an index file */
		const int fd = list->fds[iFile];
		
		/* Perform the scan */
		if (fstat(fd, &st) != 0) {
			perror("fstat");
		}
		off_t FileOffset = (off_t) 0;
		unsigned int activeBlocks = 0U;
		const char * const restrict FileName = list->FileNames[iFile];
		
		if (EndingJob) {
			BlockList = (GTLBlockDispatchFile_t * *) malloc(EndingJobEvery*sizeof(GTLBlockDispatchFile_t *));
			if (BlockList == NULL) {
				fprintf(stderr, "Unable to allocate memory for %u GTLBlockDispatchFile_t pointers\n", EndingJobEvery);
				goto bail;
			}
			semaphore = (sem_t*) malloc(sizeof(sem_t));
			if (semaphore == NULL) {
				fputs("Unable to allocate memory for semaphore\n", stderr);
				goto bail;
			}
			sem_init(semaphore, 0, 0);
		}
		
		while (1) {
			if (FileOffset >= st.st_size) break;
			
			//////////////////////////////////////////////////////////////////////////////////////////////
			// Aquire empty memory slot
			do {
				bsem_wait(thpool->donequeue_p.has_items);
				pthread_mutex_lock(&(thpool->donequeue_p.rwmutex));
				work = (GTLBlockDispatchFile_t *) jobqueue_pull(&(thpool->donequeue_p));
				pthread_mutex_unlock(&(thpool->donequeue_p.rwmutex));
			} while (work == NULL);
			
			if ((activeBlocks < EndingJobEvery) || (!EndingJob)) {
				//////////////////////////////////////////////////////////////////////////////////////////////
				// Read GTL Block headers
				const ssize_t res = pread(fd, &(work->Block.thd), sizeof(TLBHEADER), FileOffset);
				
				if (res == -1) {
					fprintf(stderr, "Failed to read header block in %s @ %lu\n", list->FileNames[iFile], FileOffset);
					perror("pread");
					goto bail;
				}
				
				work->Block.offset = FileOffset;
				work->Block.fd = fd;
				work->Block.usingIt = 0;
				work->TotalJobs = 0U;
				work->nPreviousBlocks = NULL;
				if (EndingJob) {
					BlockList[activeBlocks] = work;
					work->semaphore = semaphore;
				}
				activeBlocks++;
				FileOffset += work->Block.thd.blockLength;
			}
			else {
				work->TotalJobs = EndingJobEvery;
				work->nPreviousBlocks = BlockList;
				work->semaphore = semaphore;
				
				BlockList = (GTLBlockDispatchFile_t * *) malloc(EndingJobEvery*sizeof(GTLBlockDispatchFile_t *));
				if (BlockList == NULL) {
					fprintf(stderr, "Unable to allocate memory for %u GTLBlockDispatchFile_t pointers\n", EndingJobEvery);
					goto bail;
				}
				semaphore = (sem_t*) malloc(sizeof(sem_t));
				if (semaphore == NULL) {
					fputs("Unable to allocate memory for semaphore\n", stderr);
					goto bail;
				}
				sem_init(semaphore, 0, 0);
				activeBlocks = 0U;
			}
			
			work->params = params;
			work->FileName = FileName;
			
			pthread_mutex_lock(&(thpool->jobqueue_p.rwmutex));
			jobqueue_push(&(thpool->jobqueue_p), (job_t*) work);
			pthread_mutex_unlock(&(thpool->jobqueue_p.rwmutex));
			
		}
		
		//////////////////////////////////////////////////////////////////////////////////////////////
		// Send last Ending Job
		if ((activeBlocks > 0) && EndingJob) {
			do {
				bsem_wait(thpool->donequeue_p.has_items);
				pthread_mutex_lock(&(thpool->donequeue_p.rwmutex));
				work = (GTLBlockDispatchFile_t *) jobqueue_pull(&(thpool->donequeue_p));
				pthread_mutex_unlock(&(thpool->donequeue_p.rwmutex));
			} while (work == NULL);
			
			work->TotalJobs = activeBlocks;
			work->nPreviousBlocks = BlockList;
			work->semaphore = semaphore;
			work->FileName = FileName;
			work->params = params;
			
			semaphore = NULL;
			pthread_mutex_lock(&(thpool->jobqueue_p.rwmutex));
			jobqueue_push(&(thpool->jobqueue_p), (job_t*) work);
			pthread_mutex_unlock(&(thpool->jobqueue_p.rwmutex));
		}
		else {
			sem_destroy(semaphore);
			free(semaphore);
			free((void*) BlockList);
		}
		iFile++;
	}
	
	err = 0;
	
	bail:;
		return err;
}
//---------------------------------------------------------------
/* vim: tabstop=2 shiftwidth=2
 */
