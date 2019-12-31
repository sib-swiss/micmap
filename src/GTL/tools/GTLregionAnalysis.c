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

//=================================================================================================

#define _GNU_SOURCE     /* Expose declaration of tdestroy() */
#include "config.h"
#include <search.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <sys/time.h>

#include <sys/io.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>

#include <pthread.h>
#include <semaphore.h>
#include <assert.h>

#include <errno.h>
#include <stdatomic.h>

// basic compression
#include <zlib.h>

#include "virt_chr.h"
#include "GTL.h"
#include "gtlVersion.h"
#include "GTLdispatch.h"

//=================================================================================================
#define kMAXNBINDELS 128
#define kMAXINDELLEN 128
#define kMAXPAIRLEN (1024+10)
#define kMaxLineBuf 1024
#define kMAXPWMPOS 1024
#define kMAXDELTA 15
#define kMINCOVWIGGLEREPORT 4
#define RANGE 20

#define HasMISMATCH 1
#define HasCIGAR    2
// #define kNOT_ACGT 0xff

//=================================================================================================
typedef struct IntervalResults {
	job_t * jobs;                        	 /* pointer to all jobs memory, only for freeing */
	pthread_mutex_t thcount_lock;         /* used for thread count etc           */
	jobqueue_t jobqueue_p;                 /* pointer to the job queue            */
	jobqueue_t donequeue_p;                /* pointer to the free jobs queue      */
} IntervalResults_t;

#define kMAXPWMPOS 1024
typedef struct POS_struct {
	unsigned int chr;
	unsigned int pos;
} POS;

typedef struct OPTIONS_struct {
	int chr;
	int noCloneFilter;
	unsigned int reportVCF;
	unsigned int minorAllelePct;
	unsigned int filterQualValue;
	unsigned int PWMcnt;
	unsigned int selectAdd;
	unsigned int Realign;
	POS PWMpos[kMAXPWMPOS];
	char outbase[512];
	char selectRegions[512];
	char wiggleBed[512];
	char expressionBase[512];
	char gencodeGFF3[512];
	char tromerTranscripts[512];
	char sampleName1[64];
	char sampleName2[64];
	char refName[64];
} OPTIONS;

typedef struct Interval {
	struct Interval * prev;
	struct Counter_s {
		atomic_uint perfect;
		atomic_uint mismatch;
		atomic_uint gap;
	} Counters;
	volatile size_t size;
	volatile size_t count;
	unsigned int GenomeStart;
	unsigned int GenomeEnd;
	sem_t semaphore;
} IntervalJob_t;

typedef struct Type {
	int type;
	unsigned int Location;
	int AlignmentRange[2];
} AlignmentType_t;

//=================================================================================================
static int verbose = 0;
static int debug = 0;
static OPTIONS options;
static IntervalResults_t * restrict IntervalResultsPool = NULL;
static int noShortCut = 0;

//=================================================================================================

static void __attribute__((noreturn)) usage() 
{
		printf("usage:\n\n");
		printf("GTLregionAnalysis [options]\n\n");
//		printf("           -g genomefile         : full path of the genome encoded by encode_genome\n");
//		printf("                                   default: '%s' \n",genomefile);
		printf("           -C chromosome         : select chromosome to analyze\n");
		printf("           -v level              : verbose level\n");
		printf("           -T <uint>             : number of threads\n");
		printf("           -w <uint>             : dispatcher window width\n");
		printf("           -s <file>             : selected regions to compute coverage (bed-style)\n");
		printf("\n");
		printf(GTL_VERSION_FULL_STRING "\n");
		printf("(c) Thierry Schuepbach & Cie 2014-2018\n");
		exit(1);
}
//---------------------------------------------------------------

static IntervalResults_t * createIntervalResults(const size_t N)
{
	IntervalResults_t * const restrict pool = (IntervalResults_t*) malloc(sizeof(IntervalResults_t));
	if (pool == NULL) return NULL;
	
	
	pthread_mutex_init(&(pool->thcount_lock), NULL);
	
	if (jobqueue_init(&pool->jobqueue_p) == -1) goto err2;

	if (jobqueue_init(&pool->donequeue_p) == -1) goto err3;
	
	/* Create free memory jobs and place them into donequeue, EXTRA could be NICE ?*/
	IntervalJob_t * const restrict Jobs = (IntervalJob_t*) malloc(N*sizeof(IntervalJob_t));
	if (Jobs == NULL) goto err4;
	
	pool->jobs = (job_t*) Jobs;
	for(size_t iJob=0; iJob<N; iJob++) {
		if (sem_init(&(Jobs[iJob].semaphore), 0, 0)) {
			perror("Semaphore init interval jobs");
			goto err5;
		}
		Jobs[iJob].size = 400UL;
		Jobs[iJob].count = 0UL;
		atomic_store(&(Jobs[iJob].Counters.perfect), 0U);
		atomic_store(&(Jobs[iJob].Counters.mismatch), 0U);
		atomic_store(&(Jobs[iJob].Counters.gap), 0U);
		jobqueue_push(&pool->donequeue_p, (job_t*) &Jobs[iJob]);
	}
	
	return pool;

	err5:
		jobqueue_destroy(&pool->jobqueue_p);
	err4:
		jobqueue_destroy(&pool->donequeue_p);
	err3:
		pthread_mutex_destroy(&(pool->thcount_lock));
	err2:
		free(pool);
	return NULL;
}
//---------------------------------------------------------------

static void destroyIntervalResults(IntervalResults_t * restrict pool)
{
	/* No need to destory if it's NULL */
	if (pool == NULL) return ;

	/* Job queue cleanup */
	jobqueue_destroy(&pool->jobqueue_p);
	jobqueue_destroy(&pool->donequeue_p);
	
	if (pool->jobs) free(pool->jobs);
	
	free(pool);
	pool = NULL;
}
//---------------------------------------------------------------

static inline void __ALWAYS_INLINE
ExtractType(decodedPair_t * const restrict dpp, const unsigned char * restrict * diff,
						const unsigned char * restrict * cigar, AlignmentType_t * const restrict * ND)
{
	AlignmentType_t * restrict GD = ND[0];
	assert(dpp->taglen1<kMaxReadLen);
	int first = 1;	
	const unsigned char * restrict lcigar = *cigar;
	const unsigned char * restrict ldiff = *diff;
	
	SecondTag:	
	{
		GD->type = 0;
		GD->AlignmentRange[0] = 0;
		unsigned int cigarlength = 0U;
		if (lcigar != NULL && *lcigar != kCIGARTERMINATOR) {
			int r=0, a=0;
			while (lcigar[cigarlength] != kCIGARTERMINATOR)
			{
				register int pos = lcigar[cigarlength];
				const char code  = lcigar[cigarlength+1];
				cigarlength += 2;

				switch(code)
				{
					case 'S': // same as 'M'
					case 'M': assert(a+pos<=kMaxReadLen); { a+=pos; r+=pos; } break;
					case 'D':                             { r+=pos;         } break;
					case 'I': assert(a+pos<=kMaxReadLen); { a+=pos;         } break;
				}
			}
			assert(r >= 1);
			GD->AlignmentRange[1] = r - 1;
		}
		else {
			GD->AlignmentRange[1] = dpp->taglen1 - 1;
		}
		if (lcigar != NULL && lcigar[cigarlength] == kCIGARTERMINATOR) cigarlength++;
		if (cigarlength) GD->type |= HasCIGAR;
		lcigar += cigarlength;
		
		unsigned int count = 0;
		if (ldiff) {
			while (ldiff[count] != kMMTERMINATOR)
			{
				if (ldiff[count] == kMMSOFTCLIP)
				{
					// grab the nucleotides clipped in a tmp string
					count++;
					while (ldiff[count] != kMMTERMINATOR) { count++; }
					break; // kMMSOFTCLIP is the last thing in MMstring. only one occurence possible. we are done.
				}
				count += 2;
			}
			count++; // skip the kMMTERMINATOR
		}
		if (count) GD->type |= HasMISMATCH;
		ldiff += count;
	}

	if (first) {
		first = 0;
		GD = ND[1];
		assert(dpp->taglen2<kMaxReadLen);
		goto SecondTag;
	}
	
	*cigar = lcigar;
	*diff = ldiff;
}
//---------------------------------------------------------------

static void* FetchTag1(threadpool_t * const restrict thpool)
{
	AlignmentType_t ND_mem[2];
	AlignmentType_t * const restrict ND[2] = { &ND_mem[0], &ND_mem[1] };
	decodedPair_t dp;
	
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
			GTLBlockDispatchRange_t * const restrict task = (GTLBlockDispatchRange_t*) jobqueue_pull(&thpool->jobqueue_p);
			pthread_mutex_unlock(&thpool->jobqueue_p.rwmutex);
			
			/*************************************************************************/
			/*                             START THE TASK                            */
			/*************************************************************************/
			if (task) {
				GTLBlockInfo_t * const restrict Block = task->Block;
				IntervalJob_t * const restrict IntervalJob = (IntervalJob_t*) (task->params);
				
				/*
				 * TotalJobs triggers either:
				 * - TotalJobs = zero for fetchers that will analyze and export tags
				 * - TotalJobs > 0 for last job to collect and send further to caller
				 */
				
				if (task->TotalJobs == 0 ) {
					pthread_mutex_lock(&(Block->decompressing));
					TLBDATA* restrict td = Block->whatever;
					/* Has the block already been decompressed, if not do it */
					if (td == NULL) {
						td = (TLBDATA*) malloc(sizeof(TLBDATA)); 
						if(td == NULL) {
							fputs("Unable to allocate memory for TLBDATA\n", stderr);
							exit(1);
						}
						allocBlock(td, true);
						
						/* Save it for other to use */
						Block->whatever = td;
						Block->usingIt = 0;
						
						/* Read and decompress the block */
						if (preadDecompressBlock(Block->fd, Block->offset + sizeof(TLBHEADER), &(Block->thd), td, 0)) {
								fputs("Error decompressing block\n", stderr);
								exit(1);
						}
					}
					pthread_mutex_unlock(&(Block->decompressing));
 					
					if ((td->header.flags_7_0 & kTLBHflagPairedReads)) {
						/* Tell anyone else we are using that block */
						__sync_fetch_and_add(&(Block->usingIt), 1);

						const unsigned int GenomeWindowStart = task->GenomeStart;
						const unsigned int GenomeWindowStop  = task->GenomeEnd;
						
						const char * restrict hdr = (char *) &(td->hdrs[0]);
						const unsigned char * restrict cigar = &(td->cigar[0][0]);
						const unsigned char * restrict diff = ((td->header.flags_7_0 & 0xf0) == kTLBHflagNMismatchBlock) ? NULL : &(td->diff[0][0]);
						
						for (int i=0; i<td->cnt; i++) {
							/* Basic Decode pair */
							{
								if (td->header.flags_7_0 & kTLBHflagFixedLength) {
									dp.taglen1 = td->lengthInfo;
									dp.taglen2 = td->lengthInfo;
								}
								else {
									dp.taglen1 = td->len[2*i];
									dp.taglen2 = td->len[2*i+1];
								}
								char *next = strchr(hdr, '\t');
								hdr = next + 1;
								next = strchr(hdr, '\n');
								hdr = next + 1;
								dp.ordinal = 0;
								for (unsigned int j = 0; j < 8; j++)
									dp.ordinal = (dp.ordinal << 8) | (unsigned long) *((unsigned char *) hdr++);

								dp.reverseTAG1 = (td->ppd[0][i].tag2offset & k32bREVCOMP_TAG1_BIT) ? 1 : 0;
								dp.reverseTAG2 = (td->ppd[0][i].tag2offset & k32bREVCOMP_TAG2_BIT) ? 1 : 0;

								dp.delta = (td->ppd[0][i].tag2offset & k32bOFFSET_VALUE_MASK);
								if (td->ppd[0][i].tag2offset & k32bNEGATIVE_OFFSET_BIT)
									dp.delta = -dp.delta;

								//dp.genomepos = (virtchr[td->header.chr].chr << 28) + virtchr[td->header.chr].offset + td->ppd[0][i].tag1pos;
							}
							
							ExtractType(&dp, &diff, &cigar, ND );
							
							ND[0]->Location = td->ppd[0][i].tag1pos;
							ND[1]->Location = td->ppd[0][i].tag1pos + dp.delta; 
						
							for (int iTag=0; iTag<2; iTag++) {
								AlignmentType_t * const restrict GD = ND[iTag];
								const unsigned int TagGenLen   = GD->AlignmentRange[1] - GD->AlignmentRange[0];
								const unsigned int TagGenBegin = GD->Location;
								const unsigned int TagGenEnd   = GD->Location + TagGenLen;
								
								if (TagGenBegin <= GenomeWindowStop && TagGenEnd >= GenomeWindowStart) {
									atomic_uint * Ptr = &(IntervalJob->Counters.perfect);
									if (GD->type & HasCIGAR) {
										Ptr = &(IntervalJob->Counters.gap);
									}
									else if (GD->type & HasMISMATCH) {
										Ptr = &(IntervalJob->Counters.mismatch);
									}
									atomic_fetch_add(Ptr, 1U);
								}								
							}
						}
					
						/* Tell anyone else we no longer use that block */
						__sync_fetch_and_add(&(Block->usingIt), -1);
					}
					
					/* I have done my part , signal it */
					sem_post(&(IntervalJob->semaphore));
// 					printf("Signaling fetcher search done for file %i @ %lu\n", Block->fd, Block->offset);
				}
				else {
					/* Wait for the working threads to perform their task */
// 					printf("Waiting on %u fetchers to terminate\n", task->TotalJobs); 
					for (unsigned int k=0; k<task->TotalJobs; k++) {
						int res;
						again:;
							res = sem_wait(&(IntervalJob->semaphore));
							if (res == EINTR) goto again;
							else if (res != 0) {
								perror("sem_wait");
								exit(1);
							}
					}
					
					/* Print out interval counters */
					const unsigned int sum = atomic_load(&(IntervalJob->Counters.perfect)) \
					                       + atomic_load(&(IntervalJob->Counters.mismatch)) \
					                       + atomic_load(&(IntervalJob->Counters.gap));
					fprintf(stdout, "%i\t%u\t%u\t%u\t%u\t%u\t%u\n", options.chr, task->GenomeStart, task->GenomeEnd,
									atomic_load(&(IntervalJob->Counters.perfect)),
									atomic_load(&(IntervalJob->Counters.mismatch)),
									atomic_load(&(IntervalJob->Counters.gap)), sum);
					
					/* Recycle the interval result memory */
					pthread_mutex_lock(&(IntervalResultsPool->donequeue_p.rwmutex));
					jobqueue_push(&(IntervalResultsPool->donequeue_p), (job_t*) IntervalJob);
					pthread_mutex_unlock(&(IntervalResultsPool->donequeue_p.rwmutex));
				}
				/* Give back the memory task to the available queue */
				pthread_mutex_lock(&(thpool->donequeue_p.rwmutex));
				jobqueue_push(&(thpool->donequeue_p), (job_t*) task);
				pthread_mutex_unlock(&(thpool->donequeue_p.rwmutex));
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
	 
	return (void*) 0;
}
//---------------------------------------------------------------

static void* FetchTag2(threadpool_t * const restrict thpool)
{
	AlignmentType_t ND_mem[2];
	AlignmentType_t * const restrict ND[2] = { &ND_mem[0], &ND_mem[1] };
	decodedPair_t dp;
	
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
			GTLBlockDispatchRange_t * const restrict task = (GTLBlockDispatchRange_t*) jobqueue_pull(&thpool->jobqueue_p);
			pthread_mutex_unlock(&thpool->jobqueue_p.rwmutex);
			
			/*************************************************************************/
			/*                             START THE TASK                            */
			/*************************************************************************/
			if (task) {
				GTLBlockInfo_t * const restrict Block = task->Block;
				IntervalJob_t * const restrict IntervalJob = (IntervalJob_t*) (task->params);
				
				/*
				 * TotalJobs triggers either:
				 * - TotalJobs = zero for fetchers that will analyze and export tags
				 * - TotalJobs > 0 for last job to collect and send further to caller
				 */
				
				if (task->TotalJobs == 0 ) {
					pthread_mutex_lock(&(Block->decompressing));
					TLBDATA* restrict td = Block->whatever;
					/* Has the block already been decompressed, if not do it */
					if (td == NULL) {
						td = (TLBDATA*) malloc(sizeof(TLBDATA)); 
						if(td == NULL) {
							fputs("Unable to allocate memory for TLBDATA\n", stderr);
							exit(1);
						}
						allocBlock(td, true);
						
						/* Save it for other to use */
						Block->whatever = td;
						Block->usingIt = 0;
						
						/* Read and decompress the block */
						if (preadDecompressBlock(Block->fd, Block->offset + sizeof(TLBHEADER), &(Block->thd), td, 0)) {
								fputs("Error decompressing block\n", stderr);
								exit(1);
						}
					}
					pthread_mutex_unlock(&(Block->decompressing));
 					
					if ((td->header.flags_7_0 & kTLBHflagPairedReads)) {
						/* Tell anyone else we are using that block */
						__sync_fetch_and_add(&(Block->usingIt), 1);

						const unsigned int GenomeWindowStart = task->GenomeStart;
						const unsigned int GenomeWindowStop  = task->GenomeEnd;
						int j;						
						for (j=0; j<GTL_INDEX_SPLIT; j++) {
							//printf("%u < %u\n", GenomeWindowStart, Block->PointToFraction[j+1].Right);
							if (Block->PointToFraction[j+1].Right >= GenomeWindowStart) break;
						}
						
						/*printf("Header: %u\nCigar: %u\nMismatch: %u\nQS:%u\nID: %u\n", 
									 Block->PointToFraction[j].Header,
						       Block->PointToFraction[j].Cigar,
						       Block->PointToFraction[j].Mismatch,
						       Block->PointToFraction[j].qs,
							     Block->PointToFraction[j].ID);*/

						const int start_pair = Block->PointToFraction[j].ID;
						int stop_pair;
						for (int k=j+1; k<GTL_INDEX_SPLIT;k++) {
							if (Block->PointToFraction[k].Left > GenomeWindowStop) {
								stop_pair = Block->PointToFraction[k].ID;
								goto ShortCut;
							}
						}
						stop_pair = td->cnt;
						ShortCut:;
						
// 						const char * restrict hdr = &(td->hdrs[Block->PointToFraction[j].Header]);
						const unsigned char * restrict cigar = &(td->cigar[0][Block->PointToFraction[j].Cigar]);
						const unsigned char * restrict diff = ((td->header.flags_7_0 & 0xf0) == kTLBHflagNMismatchBlock) ? NULL : &(td->diff[0][Block->PointToFraction[j].Mismatch]);
						
						assert(start_pair <= td->cnt);
						assert(stop_pair <= td->cnt);
						for ( int i=start_pair; i<stop_pair; i++) {
							/* Basic Decode pair */
							{
								if (td->header.flags_7_0 & kTLBHflagFixedLength) {
									dp.taglen1 = td->lengthInfo;
									dp.taglen2 = td->lengthInfo;
								}
								else {
									dp.taglen1 = td->len[2*i];
									dp.taglen2 = td->len[2*i+1];
								}
// 								char *next = strchr(hdr, '\t');
// 								hdr = next + 1;
// 								next = strchr(hdr, '\n');
// 								hdr = next + 1;
// 								dp.ordinal = 0;
// 								for (unsigned int k = 0; k < 8; k++)
// 									dp.ordinal = (dp.ordinal << 8) | (unsigned long) *((unsigned char *) hdr++);

								dp.reverseTAG1 = (td->ppd[0][i].tag2offset & k32bREVCOMP_TAG1_BIT) ? 1 : 0;
								dp.reverseTAG2 = (td->ppd[0][i].tag2offset & k32bREVCOMP_TAG2_BIT) ? 1 : 0;

								dp.delta = (td->ppd[0][i].tag2offset & k32bOFFSET_VALUE_MASK);
								if (td->ppd[0][i].tag2offset & k32bNEGATIVE_OFFSET_BIT)
									dp.delta = -dp.delta;

								//dp.genomepos = (virtchr[td->header.chr].chr << 28) + virtchr[td->header.chr].offset + td->ppd[0][i].tag1pos;
							}
							
							ExtractType(&dp, &diff, &cigar, ND );
							
							ND[0]->Location = td->ppd[0][i].tag1pos;
							ND[1]->Location = td->ppd[0][i].tag1pos + dp.delta; 
						
							for (int iTag=0; iTag<2; iTag++) {
								AlignmentType_t * const restrict GD = ND[iTag];
								const unsigned int TagGenLen   = GD->AlignmentRange[1] - GD->AlignmentRange[0];
								const unsigned int TagGenBegin = GD->Location;
								const unsigned int TagGenEnd   = GD->Location + TagGenLen;
								
								if (TagGenBegin <= GenomeWindowStop && TagGenEnd >= GenomeWindowStart) {
									atomic_uint * Ptr = &(IntervalJob->Counters.perfect);
									if (i<start_pair || i>stop_pair) printf("Tag %u is missed! not int [%u, %u]\n", i, start_pair, stop_pair);
									if (GD->type & HasCIGAR) {
										Ptr = &(IntervalJob->Counters.gap);
									}
									else if (GD->type & HasMISMATCH) {
										Ptr = &(IntervalJob->Counters.mismatch);
									}
									atomic_fetch_add(Ptr, 1U);
								}								
							}
						}
					
						/* Tell anyone else we no longer use that block */
						__sync_fetch_and_add(&(Block->usingIt), -1);
					}
					
					/* I have done my part , signal it */
					sem_post(&(IntervalJob->semaphore));
// 					printf("Signaling fetcher search done for file %i @ %lu\n", Block->fd, Block->offset);
				}
				else {
					
					/* Wait for the working threads to perform their task */
// 					printf("Waiting on %u fetchers to terminate\n", task->TotalJobs); 
					for (unsigned int k=0; k<task->TotalJobs; k++) {
						int res;
						again:;
							res = sem_wait(&(IntervalJob->semaphore));
							if (res == EINTR) goto again;
							else if (res != 0) {
								perror("sem_wait");
								exit(1);
							}
					}
					
					/* Print out interval counters */
					const unsigned int sum = atomic_load(&(IntervalJob->Counters.perfect)) \
					                       + atomic_load(&(IntervalJob->Counters.mismatch)) \
					                       + atomic_load(&(IntervalJob->Counters.gap));
					fprintf(stdout, "%i\t%u\t%u\t%u\t%u\t%u\t%u\n", options.chr, task->GenomeStart, task->GenomeEnd,
									atomic_load(&(IntervalJob->Counters.perfect)),
									atomic_load(&(IntervalJob->Counters.mismatch)),
									atomic_load(&(IntervalJob->Counters.gap)), sum);
					
					/* Recycle the interval result memory */
					pthread_mutex_lock(&(IntervalResultsPool->donequeue_p.rwmutex));
					jobqueue_push(&(IntervalResultsPool->donequeue_p), (job_t*) IntervalJob);
					pthread_mutex_unlock(&(IntervalResultsPool->donequeue_p.rwmutex));
				}
				/* Give back the memory task to the available queue */
				pthread_mutex_lock(&(thpool->donequeue_p.rwmutex));
				jobqueue_push(&(thpool->donequeue_p), (job_t*) task);
				pthread_mutex_unlock(&(thpool->donequeue_p.rwmutex));
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
	 
	return (void*) 0;
}
//---------------------------------------------------------------

int main (int argc, char **argv)
{
	const char * restrict rsltfile = NULL;
	const char * restrict VCFFile = NULL;
	const char * * InputGTLFiles = NULL;
	GTLList_t * restrict GTLInputFileList = NULL;
	threadpool_t * restrict thpool = NULL;
	cpu_set_t * affinities = NULL;
	int c, err;
	unsigned int dispatcherWindow[2] = { 0U, 0U};
	unsigned int WindowsWidth = 0U;
	unsigned int nFetchThreads = 4U;
	unsigned int UseAffinities = 0U;
	
	unsigned int nInputGTLFiles = 0U;
	
	memset(&options,0,sizeof(OPTIONS));
	strcpy(options.outbase,"/tmp/");

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Process arguments

	opterr = 0;
	while ((c = getopt (argc, argv, "g:C:v:T:w:s:ahz")) != -1)
	switch (c)
	{
		case 'C':
			sscanf(optarg,"%d",&options.chr);
			break;

//		case 'g':
//			genomefile = optarg;
//			break;

		case 'v':
			sscanf(optarg,"%d",&verbose);
			break;

		case 'T':
			{
				int dummy = atoi(optarg);
				if (dummy > 0) { 
					nFetchThreads = (unsigned int) dummy;
				}
				else {
					fputs("Thread number must be positive\n", stderr);
					exit(1);
				}
			}
			break;
		case 'w':
			{
				const int itmp = atoi(optarg); 
				if (itmp <= 0) {
					fprintf(stderr, "Window width is wrong (%s)\n", optarg);
					exit(1);
				}
				WindowsWidth = (unsigned int) itmp;
				if (verbose) fprintf(stderr, "Genome window width set to %u.\n", WindowsWidth);
				
			}
			break;
		case 's':
			strncpy(options.selectRegions, optarg, 512);
			break;
		case 'a':
			UseAffinities = 1;
			break;
		case 'z':
			noShortCut = 1;
			break;
		case 'h':
			usage();
			break;
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Check requirements
	if (WindowsWidth == 0) {
		fputs("Please provide a genome floating window width!\n", stderr);
		exit(1);
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Get input file names
	if (optind < argc) {
		nInputGTLFiles = argc - optind;
		if (verbose) fprintf(stderr,"%i files provided\n", argc - optind);
		InputGTLFiles = alloca(nInputGTLFiles*sizeof(char*));
		for (unsigned int iFile=0; iFile<nInputGTLFiles; iFile++) {
			InputGTLFiles[iFile] = argv[iFile+optind];
			if (verbose) fprintf(stderr,"Input file %2i is %s\n", iFile, InputGTLFiles[iFile]);
		}
	}
	else {
		fputs("Expected input files\n", stderr);
		usage();
	}
	
	if (verbose) fprintf(stderr, "Will use %u threads for fetching data and report types\n", nFetchThreads);
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Prepare affinities
	if (UseAffinities) {
		const unsigned int nTotalThreads = nFetchThreads;
		affinities = (cpu_set_t*) malloc(nTotalThreads*sizeof(cpu_set_t));
		if (affinities != NULL) {
			unsigned int iThread;
			const unsigned int maxN = nTotalThreads < (unsigned int) sysconf(_SC_NPROCESSORS_ONLN) ? nTotalThreads :  (unsigned int) sysconf(_SC_NPROCESSORS_ONLN);
			for (iThread=0; iThread<maxN; iThread++) {
					CPU_ZERO(&affinities[iThread]);
					CPU_SET(iThread, &affinities[iThread]);
			}
			
			while (iThread <nTotalThreads) {
				CPU_ZERO(&affinities[iThread]);
				for (int j = 0; j < CPU_SETSIZE; j++) CPU_SET(j, &affinities[iThread]);
				iThread++;
			}
		}
		else {
			fprintf(stderr, "Affinities could not be allocated, keeping on without...\n");
		}
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Prepare dispatcher, indexing file if necessary and sorting
	GTLInputFileList = createGTLList(InputGTLFiles, nInputGTLFiles);
	if (GTLInputFileList == NULL) {
		fputs("Error while preparing GTL input files\n", stderr);
		goto bail;
	}
	
	if (indexGTLList(GTLInputFileList, Range, false, 0)) {
		fputs("Error indexing GTL list\n", stderr);
		goto bail;
	}
	fprintf(stderr, "Found %lu blocks\n", GTLInputFileList->nBlocks);
	if ( GTLInputFileList->PositionType != Range ) {
		fputs("Getting real Block position extrema\n", stderr);
		getRealBlockMinMaxPosition( GTLInputFileList, nFetchThreads, 1 /*DoSort*/);
		fprintf(stderr, "Writing those indices to file...\n");
		if(writeGTLFileIndex(GTLInputFileList, NULL) != 0) {
			fputs("Error in writing indices...\n", stderr);
		}
	}
	else {
		sortGTLlist(GTLInputFileList);
	}

	thpool = createGTLThreadPool((noShortCut) ? FetchTag1 : FetchTag2, affinities, sizeof(GTLBlockDispatchRange_t), nFetchThreads, 0);
	
	if (thpool == NULL) {
		fputs("Thread pool creation error\n", stderr);
		goto bail;
	}
		
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Prepare Interval Results pool
	IntervalResultsPool = createIntervalResults(2*nFetchThreads);
	if (IntervalResultsPool == NULL) {
		fprintf(stderr, "Failed to create the interval queue results' pool\n");
		exit(1);
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Maximum range available in input files
	unsigned int InputChrRangeBegin = -1;
	unsigned int InputChrRangeEnd = 0U;
	{
		const size_t nBlocks = GTLInputFileList->nBlocks;
		for (size_t iBlock=0; iBlock<nBlocks; iBlock++) {
			if (GTLInputFileList->Blocks[iBlock].thd.minPos < InputChrRangeBegin) 
				InputChrRangeBegin = GTLInputFileList->Blocks[iBlock].thd.minPos;
			if (GTLInputFileList->Blocks[iBlock].thd.maxPos > InputChrRangeEnd) 
				InputChrRangeEnd = GTLInputFileList->Blocks[iBlock].thd.maxPos;
		}
		if (verbose) fprintf(stderr, "Available chromosome %u range is [%u,%u]\n", options.chr, InputChrRangeBegin, InputChrRangeEnd);
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Output Header
	printf("#CHROM\tSTART POS\tEND POS\tnPERFECT\tnMISMATCH\tnGAP\tsum\n");

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// NO BED file
	if (options.selectRegions[0] == 0) {
		////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Get available chromosome length
		unsigned int ChrRangeBegin = -1;
		unsigned int ChrRangeEnd = 0U;
		{
			const size_t nBlocks = GTLInputFileList->nBlocks;
			for (size_t iBlock=0; iBlock<nBlocks; iBlock++) {
				if (GTLInputFileList->Blocks[iBlock].thd.minPos < ChrRangeBegin) 
					ChrRangeBegin = GTLInputFileList->Blocks[iBlock].thd.minPos;
				if (GTLInputFileList->Blocks[iBlock].thd.maxPos > ChrRangeEnd) 
					ChrRangeEnd = GTLInputFileList->Blocks[iBlock].thd.maxPos;
			}
			if (verbose) fprintf(stderr, "Available chromosome %u range is [%u,%u]\n", options.chr, ChrRangeBegin, ChrRangeEnd);
		}
		
		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Loop on genome regions
		err = 0;
		unsigned int GarbageCollect = 0U;
		IntervalJob_t * restrict work;
		const unsigned int InputChrRangeBegin = ChrRangeBegin;
		const unsigned int InputChrRangeEnd = ChrRangeEnd;
		const float invScale = 100.0f/((float) InputChrRangeEnd);
		do {
			///////////////////////////////////////////////////////////////////////////////////////////////////////
			// Acquire a free Interval memory slot
			if (err == 0) {
				bsem_wait(IntervalResultsPool->donequeue_p.has_items);
				pthread_mutex_lock(&(IntervalResultsPool->donequeue_p.rwmutex));
				work = (IntervalJob_t *) jobqueue_pull(&(IntervalResultsPool->donequeue_p));
				pthread_mutex_unlock(&(IntervalResultsPool->donequeue_p.rwmutex));
			}

			///////////////////////////////////////////////////////////////////////////////////////////////////////
			// Dispatch GTL blocks to threads
			work->count = 0U;
			work->GenomeStart = ChrRangeBegin;
			ChrRangeBegin += WindowsWidth;
			work->GenomeEnd = ChrRangeBegin - 1;
			atomic_store(&(work->Counters.perfect), 0U);
			atomic_store(&(work->Counters.mismatch), 0U);
			atomic_store(&(work->Counters.gap), 0U);
							
			if (work->GenomeEnd > ChrRangeEnd) work->GenomeEnd = ChrRangeEnd;
			//fprintf(stderr, "Submitting window %u--%u\n", work->GenomeStart, work->GenomeEnd);
			err = dispatchGTLBlockRange(GTLInputFileList, thpool, options.chr, work->GenomeStart, work->GenomeEnd, -1, 1, work);
			if (err < 0) {
				fprintf(stderr, "dispatcher returned error code %i\n", err);
				goto bail1;
			}
			
			if ((++GarbageCollect & 0xFFF) == 0) {
				/* Free potential blocks */ 
				for (size_t iBlock=0; iBlock<GTLInputFileList->nBlocks; iBlock++) {
					if (GTLInputFileList->Blocks[iBlock].thd.maxPos < dispatcherWindow[0]) {
						if (GTLInputFileList->Blocks[iBlock].whatever != NULL && GTLInputFileList->Blocks[iBlock].usingIt == 0) {
							pthread_mutex_lock(&(GTLInputFileList->Blocks[iBlock].decompressing));
							freeBlock(GTLInputFileList->Blocks[iBlock].whatever);
							GTLInputFileList->Blocks[iBlock].whatever = NULL;
							pthread_mutex_unlock(&(GTLInputFileList->Blocks[iBlock].decompressing));
						}
					}
				}
				fprintf(stderr, "%15u\t%15u\t%15u\t%5.1f\r", InputChrRangeBegin, ChrRangeBegin, InputChrRangeEnd, ((float) ChrRangeBegin)*invScale);
			}
		} while (ChrRangeBegin < ChrRangeEnd);
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// BED file
	else {
		FILE * bedFile = fopen(options.selectRegions, "r");
		if (bedFile) {
			unsigned int GarbageCollect = 0U;
			IntervalJob_t * restrict work;
			err = 0;
			const float invScale = 100.0f/((float) InputChrRangeEnd);
			unsigned int line = 0U;
			char ChrString[16];
			while (!feof(bedFile)) {
				unsigned int ChrRangeBegin, ChrRangeEnd, Chr;
				const int num = fscanf(bedFile, "%s\t%u\t%u\n", ChrString, &ChrRangeBegin, &ChrRangeEnd);
				if (ChrString[0] == 'X' || ChrString[0] == 'x') {
					Chr = 23;
				}
				else if (ChrString[0] == 'Y' || ChrString[0] == 'y') {
					Chr = 24;
				}
				else {
					Chr = (unsigned int) atoi(ChrString);
				}
				
				line++;
				if (num == 3 && Chr == options.chr) {
// 					if (err == 0) {
						do {
							//////////////////////////////////////////////////////////////////////////////////////////////
							// Acquire a free Interval memory slot
							if (err == 0) {
								bsem_wait(IntervalResultsPool->donequeue_p.has_items);
								pthread_mutex_lock(&(IntervalResultsPool->donequeue_p.rwmutex));
								work = (IntervalJob_t *) jobqueue_pull(&(IntervalResultsPool->donequeue_p));
								pthread_mutex_unlock(&(IntervalResultsPool->donequeue_p.rwmutex));
							}

							//////////////////////////////////////////////////////////////////////////////////////////////
							// Dispatch GTL blocks to threads
							work->count = 0U;
							work->GenomeStart = ChrRangeBegin;
							ChrRangeBegin += WindowsWidth;
							work->GenomeEnd = ChrRangeBegin - 1;
							atomic_store(&(work->Counters.perfect), 0U);
							atomic_store(&(work->Counters.mismatch), 0U);
							atomic_store(&(work->Counters.gap), 0U);
							if (work->GenomeEnd > ChrRangeEnd) work->GenomeEnd = ChrRangeEnd;
							err = dispatchGTLBlockRange(GTLInputFileList, thpool, options.chr, work->GenomeStart, work->GenomeEnd, -1, 1, work);
							if (err < 0) {
								fprintf(stderr, "dispatcher returned error code %i\n", err);
								goto bail1;
							}
							
							if ((++GarbageCollect & 0xFFF) == 0) {
								/* Free potential blocks */ 
								for (size_t iBlock=0; iBlock<GTLInputFileList->nBlocks; iBlock++) {
									if (GTLInputFileList->Blocks[iBlock].thd.maxPos < work->GenomeStart) {
										if (GTLInputFileList->Blocks[iBlock].whatever != NULL 
										 && GTLInputFileList->Blocks[iBlock].usingIt == 0) {
											pthread_mutex_lock(&(GTLInputFileList->Blocks[iBlock].decompressing));
											freeBlock(GTLInputFileList->Blocks[iBlock].whatever);
											free(GTLInputFileList->Blocks[iBlock].whatever);
											GTLInputFileList->Blocks[iBlock].whatever = NULL;
											pthread_mutex_unlock(&(GTLInputFileList->Blocks[iBlock].decompressing));
										}
									}
								}
							}
							
//							fprintf(stderr, "%15u\t%15u\t%15u\t%5.1f\r", InputChrRangeBegin, ChrRangeBegin, InputChrRangeEnd, ((float) ChrRangeBegin)*invScale);
							
						} while (ChrRangeBegin < ChrRangeEnd);
// 					}
				}
				else if (num != 3) {
					fprintf(stderr, "Error reading region coordinate, got %i values out of 3 on line %u\n", num, line);
					fclose(bedFile);
					goto bail1;
				}
			}
			fclose(bedFile);
		}
		else {
			fprintf(stderr, "bed file %s is not available.\n", options.selectRegions);
			goto bail1;
		}
	}
	
bail1:;
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Wait for jobs to terminate and empty queues
	thpool_wait(thpool);
	
	destroyIntervalResults(IntervalResultsPool);
	destroyGTLThreadPool(thpool);
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Clean up and quit

bail:;
	fflush(stdout);
	if (GTLInputFileList) destroyGTLList(GTLInputFileList);

	if (affinities) free(affinities);
	
	return err;

} /* main */
//---------------------------------------------------------------

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */

