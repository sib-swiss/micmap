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
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>

#include <pthread.h>
#include <signal.h>
#include <semaphore.h>
#include <assert.h>
#include <errno.h>
#include <stdatomic.h>

#include "Genome.h"
#include "virt_chr.h"
#include "GTL.h"
#include "gtlVersion.h"
#include "GTLdispatch.h"
#include "htslib/vcf.h"
#include "htslib/vcfutils.h"

//=================================================================================================
#define kMAXNBINDELS 128
#define kMAXINDELLEN 128
#define kMAXPAIRLEN (1024+10)
#define kMaxLineBuf 1024
#define kMAXPWMPOS 1024
#define kMAXDELTA 15
#define kMINCOVWIGGLEREPORT 4
#define RANGE 20
// #define kNOT_ACGT 0xff

//=================================================================================================
#include "AnalyzeTags.h"

//=================================================================================================
unsigned int longuestAllowedDeletion = 2000; 

static int verbose = 0;
static int debug = 0;
OPTIONS options = {	
	.chr = 0,
	.noCloneFilter = 0,
	.reportVCF = 0,
	.minorAllelePct = 15,
	.filterQualValue = 0,
	.PWMcnt = 0,
	.selectAdd = 0,
	.Realign = 0,
	.PWMpos = { [0 ... kMAXPWMPOS-1] = {0,0}},
	.outbase = "/tmp/",
	.selectRegions = { [0 ... 511] = '\0'},
	.wiggleBed = { [0 ... 511] = '\0'},
	.expressionBase = { [0 ... 511] = '\0'},
	.gencodeGFF3 = { [0 ... 511] = '\0'},
	.tromerTranscripts = { [0 ... 511] = '\0'},
	.sampleName1 = { [0 ... 63] = '\0'},
	.sampleName2 = { [0 ... 63] = '\0'},
	.refName = { [0 ... 63] = '\0'}
};
static char default_genome[] = "/tmp/nodelete/nguex/data/hg19.bin";
// const char * restrict genomefile = default_genome;

static threadpool_t * restrict IntervalResultsPool = NULL;
bcf_hdr_t * hdr = NULL;

Genome_t Genome = {
	.FileName = &default_genome[0],
	.table = NULL,
	.configBuffer = NULL,
	.tableSize = 0UL,
	.configBufferSize = 0UL
};

atomic_uint nVCFRecord;
atomic_uint nRecordsFound;
atomic_uint nZygoteError;
atomic_uint nLocationNotFound;
atomic_uint nAllele0NotFound;
atomic_uint nAllele1NotFound;
atomic_uint nAlignmentFailure;
atomic_uint nNotEnoughCoverage;
atomic_uint nRemovedHasFork;
atomic_uint nNoVariant;

extern FILE * outfile;

//=================================================================================================
static void __attribute__((noreturn)) usage() 
{
		printf("usage:\n\n");
		printf("Nirvana [options]\n\n");
		printf("           -g genomefile         : full path of the genome encoded by encode_genome\n");
		printf("                                   default: '%s' \n",default_genome);
		printf("           -C chromosome         : select chromosome to analyze\n");
		printf("           -f <char>             : filter out nucleotides where the quality value is *below* the specified character [%c]\n", options.filterQualValue + '#');

		printf("           -m <int>              : minor allele minimum percent coverage\n");
		printf("           -v level              : verbose level\n");
		printf("           -T <uint>/<uint>      : number of thread to fetch tags / number of thread to analyze\n");
		printf("           -A <vcf file>         : analyze and report position given in vcf file\n");
		printf("           -R <uint>             : margin left and right to define window\n");
		printf("           -s <bed file>         : region list to filter\n");
		printf("           -d                    : deprecated version of analysis\n");
		printf("\n");
		printf(GTL_VERSION_FULL_STRING "\n");
		printf("(c) Thierry Schuepbach & Cie 2014-2019\n");
		exit(1);
}
//---------------------------------------------------------------

static threadpool_t* createIntervalResultsPool(void* (*Fct)(threadpool_t * const restrict),
																							 const cpu_set_t * const restrict affinities,
                                               const int nThreads, const int onHold)
{
	threadpool_t * const restrict thpool = (threadpool_t*) malloc(sizeof(threadpool_t));
	if (thpool == NULL) return NULL;
	
	thpool->threads = (pthread_t*) malloc(nThreads*sizeof(pthread_t));
	if (thpool->threads == NULL) goto err1;
	
	thpool->common = (void*) &options;
	thpool->threads_on_hold   = onHold;
	thpool->threads_keepalive = 1;
	
	thpool->num_threads         = nThreads;
	thpool->num_threads_alive   = 0;
	thpool->num_threads_working = 0;
	
	pthread_mutex_init(&(thpool->thcount_lock), NULL);
	
	if (jobqueue_init(&thpool->jobqueue_p) == -1) goto err2;

	if (jobqueue_init(&thpool->donequeue_p) == -1) goto err3;
	
	/* Create free memory jobs and place them into donequeue, EXTRA could be NICE ?*/
	IntervalJob_t * const restrict Jobs = (IntervalJob_t*) calloc((3*nThreads), sizeof(IntervalJob_t));
	if (Jobs == NULL) goto err4;
	
	thpool->jobs = (job_t*) Jobs;
	for(int iJob=0; iJob<(3*nThreads); iJob++) {
		if (sem_init(&(Jobs[iJob].semaphore), 0, 0)) {
			perror("Semaphore init interval jobs");
			goto err5;
		}
		Jobs[iJob].size = 400UL;
		if (pthread_mutex_init(&(Jobs[iJob].mutex), NULL) ) {
			perror("Mutex init interval jobs");
			goto err5;
		}
		Jobs[iJob].data = (void**) calloc(Jobs[iJob].size, sizeof(void*));
		if (Jobs[iJob].data == NULL) {
			while (--iJob) {
				free(Jobs[iJob].data);
			}
			fprintf(stderr, "Unable to allocate memory for Results pointers\n");
			goto err5;
		}
		
		Jobs[iJob].VCFRecord = bcf_init();
		if (Jobs[iJob].VCFRecord == NULL) {
			while (--iJob) {
				free(Jobs[iJob].data);
				bcf_destroy(Jobs[iJob].VCFRecord);
			}
			fputs("Unable to initialize record\n", stderr);
			goto err5;
		}
		jobqueue_push(&thpool->donequeue_p, (job_t*) &Jobs[iJob]);
	}
	
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, (size_t) (2*1024*1024));
// 	pthread_attr_setguardsize(&attr, (size_t) (100*4096));
// 	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
	if (affinities) {
		for (int iThread = 0; iThread<nThreads; iThread++) {
			pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &affinities[iThread]);
			if (pthread_create(&thpool->threads[iThread], &attr, (void* (*)(void*)) Fct, thpool) != 0) {
				while (--iThread > 0) {
					pthread_kill(thpool->threads[iThread], SIGKILL);
				}
				goto err5;
			}
		}
	}
	else {
		for (int iThread = 0; iThread<nThreads; iThread++) {
			if (pthread_create(&thpool->threads[iThread], &attr, (void* (*)(void*)) Fct, thpool) != 0) {
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

static void destroyIntervalResultsPool(threadpool_t * restrict thpool)
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
	jobqueue_destroy(&thpool->donequeue_p);
	
	{
		IntervalJob_t * const restrict Jobs = (IntervalJob_t *) thpool->jobs;
		if (Jobs) {
			for(int iJob=0; iJob<(3*thpool->num_threads); iJob++) {
				if (Jobs[iJob].data) free(Jobs[iJob].data);
				pthread_mutex_destroy(&(Jobs[iJob].mutex));
			}
			free(Jobs);
		}		
	}
	
	free(thpool->threads);
	free(thpool);
	thpool = NULL;
}
//---------------------------------------------------------------

static void* FetchTag(threadpool_t * const restrict thpool)
{
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
						
						unsigned char * restrict diff = ((td->header.flags_7_0 & 0xf0) == kTLBHflagNMismatchBlock) ? NULL : td->diff[0];
						decodedPair_t dp;
						unsigned int mmpos = 0U;
						unsigned int cigarpos = 0U;
						GTLRawData_t *ND[2] = { NULL, NULL };
						const unsigned int GenomeWindowStart = task->GenomeStart;
						const unsigned int GenomeWindowStop  = task->GenomeEnd;
						
						const char * restrict hdr = (char *) td->hdrs;
						const unsigned char * restrict qs = td->qual;
						
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
								const char *next = strchr(hdr, '\t');
								hdr = next + 1;
								next = strchr(hdr, '\n');
								hdr = next + 1;
								dp.ordinal = 0;
								for (unsigned int j = 0; j < 8; j++)
									dp.ordinal = (dp.ordinal << 8) | (unsigned long) *((unsigned char *) hdr++);

								if (td->ppd[0][i].tag2offset & k32bREVCOMP_TAG1_BIT)
									dp.reverseTAG1 = 1;
								else
									dp.reverseTAG1 = 0;

								if (td->ppd[0][i].tag2offset & k32bREVCOMP_TAG2_BIT)
									dp.reverseTAG2 = 1;
								else
									dp.reverseTAG2 = 0;

								dp.delta = (td->ppd[0][i].tag2offset & k32bOFFSET_VALUE_MASK);
								if (td->ppd[0][i].tag2offset & k32bNEGATIVE_OFFSET_BIT)
									dp.delta = -dp.delta;

								dp.genomepos = (Genome.virtchr[td->header.chr].chr << 28) + Genome.virtchr[td->header.chr].offset + td->ppd[0][i].tag1pos;
							}
							
							for (unsigned int l=0; l<2; l++) {
								if (ND[l] == NULL) {
									ND[l] = (GTLRawData_t*) calloc(1,sizeof(GTLRawData_t));
									if (ND[l] == NULL) {
										fputs("Unable to allocate memory for GTLRawData_t\n", stderr);
										exit(1);
									}
								}
								else 
									memset(ND[l], 0, sizeof(GTLRawData_t));
							}
							
							ExtractDataNasSoftClip(&dp, (diff != NULL) ? &diff[mmpos] : NULL, td->cigar[0] + cigarpos, ND );
							
							ND[0]->Location = td->ppd[0][i].tag1pos;
							ND[1]->Location = td->ppd[0][i].tag1pos + dp.delta; 
						
							for (int iTag=0; iTag<2; iTag++) {
								GTLRawData_t * const restrict GD = ND[iTag];
								const unsigned int TagGenLen   = GD->AlignmentRange[1] - GD->AlignmentRange[0];
								const unsigned int TagGenBegin = GD->Location;
								const unsigned int TagGenEnd   = GD->Location + TagGenLen;
								
								if (TagGenBegin <= GenomeWindowStop && TagGenEnd >= GenomeWindowStart) {
									GD->Ordinal = dp.ordinal;
									if (!GD->revNmer) {
										memcpy(GD->qs, qs, GD->TagLength);
									}
									else {
										unsigned char * restrict ptr = &(GD->qs[GD->TagLength-1]);
										for(int l=0; l<GD->TagLength; l++ ) {
											*ptr-- = qs[l];
										}
									}
									pthread_mutex_lock(&(IntervalJob->mutex));
									if (IntervalJob->count >= IntervalJob->size) {
										IntervalJob->size = IntervalJob->count + 100UL;
										IntervalJob->data = (void**) realloc(IntervalJob->data, IntervalJob->size*sizeof(void*));
										if (IntervalJob->data == NULL) {
											fputs("Unable to reallocate more memory for Results pointer array\n", stderr);
											exit(1);
										}
									}
									IntervalJob->data[IntervalJob->count++] = GD;
									pthread_mutex_unlock(&(IntervalJob->mutex));
									
									/* Tag is useful, next use shall reallocate */
									ND[iTag] = NULL; 
								}
								
								/* Move pointers */
								qs       += GD->TagLength;
								cigarpos += GD->CigarLength;
								mmpos    += GD->MismatchLength;
							}
						}
					
						/* Tell anyone else we no longer use that block */
						__sync_fetch_and_add(&(Block->usingIt), -1);
						
						/* Free pot ential ND memory */
						if (ND[0]) free(ND[0]);
						if (ND[1]) free(ND[1]);
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
// 					printf("There were %u jobs started\n", task->TotalJobs);
					jobqueue_t * queue;
					if (IntervalJob->count >= 4 ) {
// 						printf("Sending search result to interval queue\n");
						queue = &(IntervalResultsPool->jobqueue_p);
					}
					else {
// 						fprintf(stderr, "Interval %u--%u does not have enough coverage: %lu tags found\n", task->GenomeStart, task->GenomeEnd, IntervalJob->count);
// 						fprintf(stderr, "Interval %u--%u : %lu tags found\n", task->GenomeStart, task->GenomeEnd, IntervalJob->count);
						
						/* Free the memory */
						for (unsigned int k=0; k<IntervalJob->count; k++) {
							free(IntervalJob->data[k]);
						}
						queue = &(IntervalResultsPool->donequeue_p);
					}
					/* Recycle the interval result memory */
					pthread_mutex_lock(&(queue->rwmutex));
					jobqueue_push(queue, (job_t*) IntervalJob);
					pthread_mutex_unlock(&(queue->rwmutex));
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
						decodedPair_t dp;
						GTLRawData_t *ND[2] = { NULL, NULL };
						const unsigned int GenomeWindowStart = task->GenomeStart;
						const unsigned int GenomeWindowStop  = task->GenomeEnd;
						
						int StartingBlock, start_pair, stop_pair;
						for (StartingBlock=0; StartingBlock<(GTL_INDEX_SPLIT-1); StartingBlock++) {
							//printf("%u < %u\n", GenomeWindowStart, Block->PointToFraction[j+1].Right);
							if (Block->PointToFraction[StartingBlock+1].Right >= GenomeWindowStart) break;
						}
						
						if (StartingBlock ==(GTL_INDEX_SPLIT-1)) {
							start_pair = Block->PointToFraction[--StartingBlock].ID;
							stop_pair = td->cnt;
						}
						else {
							start_pair = Block->PointToFraction[StartingBlock].ID;
							stop_pair = td->cnt;
							for (int k=StartingBlock+1; k<GTL_INDEX_SPLIT;k++) {
								if (Block->PointToFraction[k].Left > GenomeWindowStop) {
									stop_pair = Block->PointToFraction[k].ID;
									break;
								}
							}
						}
						
						unsigned char * restrict diff = ((td->header.flags_7_0 & 0xf0) == kTLBHflagNMismatchBlock) ? NULL : td->diff[0];
						const char * restrict hdr = (char *) &(td->hdrs[Block->PointToFraction[StartingBlock].Header]);
						unsigned int cigarpos = Block->PointToFraction[StartingBlock].Cigar;
						unsigned int mmpos = Block->PointToFraction[StartingBlock].Mismatch;
						const unsigned char * restrict qs = &(td->qual[Block->PointToFraction[StartingBlock].qs]);
						const char * const restrict hdrMax = (char *) &(td->hdrs[td->readHdr.origLen]);
						
						assert(start_pair <= td->cnt);
						assert(stop_pair <= td->cnt);
// 						if (stop_pair > td->cnt) {
// 							printf("Block %lu is bad [%i,%i], %i > %i\n", StartingBlock, start_pair, stop_pair, stop_pair, td->cnt);
// 							printf("Genome: %u - %u\n", GenomeWindowStart, GenomeWindowStop);
// 							for (int j=0; j<32; j++) {
// 									printf("%i\t%u\t%u\n", Block->PointToFraction[j].ID, Block->PointToFraction[j].Left, Block->PointToFraction[j].Right);
// 							}
// 							assert(false);
// 						}
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
								assert(hdr);
								assert((uintptr_t) hdr < (uintptr_t) &td->hdrs[td->readHdr.origLen]);
								while (*hdr != '\t' && hdr < hdrMax) {
									assert((uintptr_t) hdr < (uintptr_t) &td->hdrs[td->readHdr.origLen]);
									++hdr;
								}
								while (*hdr != '\n' && hdr < hdrMax) { 
									assert((uintptr_t) hdr < (uintptr_t) &td->hdrs[td->readHdr.origLen]);
									++hdr;
								}
								if (hdr >= hdrMax) {
									fputs("Error reading block, end of header block reached!\n", stderr);
									exit(1);
								}
								++hdr;
								assert((uintptr_t) hdr < (uintptr_t) &td->hdrs[td->readHdr.origLen]);

								dp.ordinal = 0UL;
								for (unsigned int k = 0; k < 8; k++)
									dp.ordinal = (dp.ordinal << 8) | (unsigned long) *((unsigned char *) hdr++);
								
								dp.reverseTAG1 = (td->ppd[0][i].tag2offset & k32bREVCOMP_TAG1_BIT) ? 1 : 0;
								dp.reverseTAG2 = (td->ppd[0][i].tag2offset & k32bREVCOMP_TAG2_BIT) ? 1 : 0;

								dp.delta = (td->ppd[0][i].tag2offset & k32bOFFSET_VALUE_MASK);
								if (td->ppd[0][i].tag2offset & k32bNEGATIVE_OFFSET_BIT)
									dp.delta = -dp.delta;

								dp.genomepos = (Genome.virtchr[td->header.chr].chr << 28) + Genome.virtchr[td->header.chr].offset + td->ppd[0][i].tag1pos;
							}
							
							for (unsigned int l=0; l<2; l++) {
								if (ND[l] == NULL) {
									ND[l] = (GTLRawData_t*) calloc(1,sizeof(GTLRawData_t));
									if (ND[l] == NULL) {
										fputs("Unable to allocate memory for GTLRawData_t\n", stderr);
										exit(1);
									}
								}
								else 
									memset(ND[l], 0, sizeof(GTLRawData_t));
							}
							
							ExtractDataNasSoftClip(&dp, (diff != NULL) ? &diff[mmpos] : NULL, td->cigar[0] + cigarpos, ND );
							
							ND[0]->Location = td->ppd[0][i].tag1pos;
							ND[1]->Location = td->ppd[0][i].tag1pos + dp.delta; 
						
							for (int iTag=0; iTag<2; iTag++) {
								GTLRawData_t * const restrict GD = ND[iTag];
								const unsigned int TagGenLen   = GD->AlignmentRange[1] - GD->AlignmentRange[0];
								const unsigned int TagGenBegin = GD->Location;
								const unsigned int TagGenEnd   = GD->Location + TagGenLen;
								
								if (TagGenBegin <= GenomeWindowStop && TagGenEnd >= GenomeWindowStart) {
									GD->Ordinal = dp.ordinal;
									if (!GD->revNmer) {
										memcpy(GD->qs, qs, GD->TagLength);
									}
									else {
										unsigned char * restrict ptr = &(GD->qs[GD->TagLength-1]);
										for(int l=0; l<GD->TagLength; l++ ) {
											*ptr-- = qs[l];
										}
									}
									pthread_mutex_lock(&(IntervalJob->mutex));
									if (IntervalJob->count >= IntervalJob->size) {
										IntervalJob->size = IntervalJob->count + 100UL;
										IntervalJob->data = (void**) realloc(IntervalJob->data, IntervalJob->size*sizeof(void*));
										if (IntervalJob->data == NULL) {
											fputs("Unable to reallocate more memory for Results pointer array\n", stderr);
											exit(1);
										}
									}
									IntervalJob->data[IntervalJob->count++] = GD;
									pthread_mutex_unlock(&(IntervalJob->mutex));
									
									/* Tag is useful, next use shall reallocate */
									ND[iTag] = NULL; 
								}
								
								/* Move pointers */
								qs       += GD->TagLength;
								cigarpos += GD->CigarLength;
								mmpos    += GD->MismatchLength;
							}
						}
					
						/* Tell anyone else we no longer use that block */
						__sync_fetch_and_add(&(Block->usingIt), -1);
						
						/* Free potential ND memory */
						if (ND[0]) free(ND[0]);
						if (ND[1]) free(ND[1]);
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
// 					printf("There were %u jobs started\n", task->TotalJobs);
					jobqueue_t * queue;
					if (IntervalJob->count >= 4 ) {
// 						printf("Sending search result to interval queue\n");
						queue = &(IntervalResultsPool->jobqueue_p);
					}
					else {
// 						fprintf(stderr, "Interval %u--%u does not have enough coverage: %lu tags found\n", task->GenomeStart, task->GenomeEnd, IntervalJob->count);
// 						fprintf(stderr, "Interval %u--%u : %lu tags found\n", task->GenomeStart, task->GenomeEnd, IntervalJob->count);
						
						/* Free the memory */
						for (unsigned int k=0; k<IntervalJob->count; k++) {
							free(IntervalJob->data[k]);
						}
						queue = &(IntervalResultsPool->donequeue_p);
					}
					/* Recycle the interval result memory */
					pthread_mutex_lock(&(queue->rwmutex));
					jobqueue_push(queue, (job_t*) IntervalJob);
					pthread_mutex_unlock(&(queue->rwmutex));
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

// #include "AnalyzeTags2.c"
void* AnalyzeTags(threadpool_t * restrict const thpool);
void* AnalyzeTags3(threadpool_t * restrict const thpool);

int main (int argc, char **argv)
{
	const char * restrict rsltfile = NULL;
	const char * restrict VCFFile = NULL;
	const char * * InputGTLFiles = NULL;
	const char * OutFileName = NULL;
	const char * restrict BedFile = NULL;
	GTLList_t * restrict GTLInputFileList = NULL;
	threadpool_t * restrict thpool = NULL;
	unsigned int (*restrict Filtering)[2] = NULL;
	void* (*AnalysisFct)(threadpool_t * restrict const) = AnalyzeTags3;
	cpu_set_t * affinities = NULL;
	
	int c, err=1, KeepWindow=0;
	unsigned int dispatcherWindow[2] = { 0U, 0U};
	unsigned int nFetchThreads = 4U;
	unsigned int nAnalizeThreads = (unsigned int) sysconf(_SC_NPROCESSORS_ONLN) - 4U;
	unsigned int UseAffinities = 0U;
	
	unsigned int nInputGTLFiles = 0U;
	unsigned int noShortCut = 0U;
	int Margin = RANGE;
	int nBedRegions = 0;
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Process arguments
	opterr = 0;
	while ((c = getopt (argc, argv, "g:C:v:f:m:T:A:hakR:o:s:zd")) != -1)
	switch (c)
	{
		case 'C':
			options.chr = atoi(optarg);
			break;

		case 'g':
			Genome.FileName = optarg;
			break;

		case 'f':
			options.filterQualValue = (unsigned int) *optarg;
			break;

		case 'm':
			options.minorAllelePct = atoi(optarg);
			break;

		case 'v':
			sscanf(optarg,"%d",&verbose);
			break;

		case 'T':
			{
				int dummy[2];
				if (sscanf(optarg, "%i/%i", &dummy[0], &dummy[1]) == 2 ) {
					if (dummy[0] > 0 && dummy[1] > 0) { 
						nFetchThreads = (unsigned int) dummy[0];
						nAnalizeThreads = (unsigned int) dummy[1];
					}
					else {
						fputs("Thread number must be positive\n", stderr);
						exit(1);
					}
				}
				else {
					fprintf(stderr, "Unable to parse the thread numbers in '%s'\n", optarg);
					exit(1);
				}
			}
			break;
		case 'A':
			VCFFile = optarg;
			break;
		case 'k':
			KeepWindow = 1;
			break;
		case 'R':
				Margin = atoi(optarg);
				break;
		case 'o':
			OutFileName = optarg;
			break;
		case 's':
			BedFile = optarg;
			break;
		case 'a':
			UseAffinities = 1U;
			break;
		case 'z':
			noShortCut = 1U;
			break;
		case 'd':
			AnalysisFct = AnalyzeTags;
			break;
		case 'h':
		default:
			usage();
			break;
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Initialize
	atomic_store(&nVCFRecord, 0U);
	atomic_store(&nRecordsFound, 0U);
	atomic_store(&nZygoteError, 0U);
	atomic_store(&nLocationNotFound, 0U);
	atomic_store(&nAllele0NotFound, 0U);
	atomic_store(&nAllele1NotFound, 0U);
	atomic_store(&nAlignmentFailure, 0U);
	atomic_store(&nNotEnoughCoverage, 0U);
	atomic_store(&nRemovedHasFork, 0U);
	atomic_store(&nNoVariant, 0U);
	
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
	if (VCFFile == NULL) {
		fputs("Please provid ea VCF file to check\n", stderr);
		usage();
	}
	
	if (options.sampleName1[0] == 0) snprintf(options.sampleName1, 64, "%s", rsltfile);

	if (verbose) fprintf(stderr, "Will use %u threads for fetching data, %u for analyzing\n", nFetchThreads, nAnalizeThreads);
	
	if (OutFileName) {
		outfile = fopen(OutFileName, "w");
		if (outfile == NULL) fprintf(stderr, "Unable to open file %s\n", OutFileName);
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Potential Bed file
	if (BedFile) {
		FILE * bed = fopen(BedFile, "r");
		if (bed) {
			unsigned int lines = 0UL;
			char ChrString[32];
			unsigned int range[2];
			while (!feof(bed)) {
				if (fscanf(bed, "%31s\t%u\t%u\n", ChrString, &range[0], &range[1]) != 3) {
					fprintf(stderr, "Unable to read bed file at line %u\n", lines+1);
					fclose(bed);
					exit(1);
				}
				const char * ChrPtr = ChrString;
				if (ChrString[0] == 'c' || ChrString[0] == 'C') ChrPtr += 3;
				int chr = 0;
				if (ChrPtr[0] == 'x' || ChrPtr[0] == 'X') chr = 23;
				else if (ChrPtr[0] == 'y' || ChrPtr[0] == 'Y') chr = 24;
				else if (ChrPtr[0] == 'm' || ChrPtr[0] == 'M') chr = 25;
				if (chr == 0) chr = atoi(ChrPtr);
				if (chr == options.chr) lines++;
			}
			if (verbose) printf("Bed file %s holds %u lines\n", BedFile, lines);
			Filtering = (unsigned int (*)[2]) malloc(2*lines*sizeof(unsigned int));
			if (Filtering == NULL) {
				fputs("Unable to allocate memory for bed filtering\n", stderr);
				fclose(bed);
				exit(1);
			}
			rewind(bed);
			nBedRegions = lines;
			lines = 0UL;
			while (!feof(bed)) {
				if (fscanf(bed, "%31s\t%u\t%u\n", ChrString, &range[0], &range[1]) != 3) {
					fprintf(stderr, "Unable to read bed file at line %u\n", lines+1);
					fclose(bed);
					exit(1);
				}
				const char * ChrPtr = ChrString;
				if (ChrString[0] == 'c' || ChrString[0] == 'C') ChrPtr += 3;
				int chr = 0;
				if (ChrPtr[0] == 'x' || ChrPtr[0] == 'X') chr = 23;
				else if (ChrPtr[0] == 'y' || ChrPtr[0] == 'Y') chr = 24;
				else if (ChrPtr[0] == 'm' || ChrPtr[0] == 'M') chr = 25;
				if (chr == 0) chr = atoi(ChrPtr);		
				if (chr == options.chr) { 
					if (lines >= nBedRegions) { fprintf(stderr, "going beyond control at line %u\n" , lines);
						exit(1);
					}
					Filtering[lines][0] = range[0]; Filtering[lines][1] = range[1];
					lines++;
				}
			}
			fclose(bed);
		}
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Prepare affinities
	if (UseAffinities) {
		const unsigned int nTotalThreads = nAnalizeThreads + nFetchThreads;
		affinities = (cpu_set_t*) malloc(nTotalThreads*sizeof(cpu_set_t));
		if (affinities) {
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
	// Load genome
	if (verbose) fprintf(stderr, "Loading genome file %s\n", Genome.FileName);
	if (LoadGenome(&Genome) != NULL) {
		fprintf(stderr, "failed to load genome file %s\n", Genome.FileName); 
		err = 1;
		goto bail;
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
	
	thpool = createGTLThreadPool((noShortCut) ? FetchTag : FetchTag2, affinities, sizeof(GTLBlockDispatchRange_t), nFetchThreads, 0);
	
	if (thpool == NULL) {
		fputs("Thread pool creation error\n", stderr);
		goto bail;
	}
	/* Assign the window pointer to threadpool common data place holder*/
	thpool->common = (volatile void * volatile) &dispatcherWindow[0]; // Check this out, I do not think it is necessary anymore
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Prepare Interval Results pool
	IntervalResultsPool = createIntervalResultsPool(AnalysisFct, NULL, nAnalizeThreads, 0);
	if (IntervalResultsPool == NULL) {
		fprintf(stderr, "Failed to create the interval queue results' pool\n");
		exit(1);
	}
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Open VCF file
	htsFile * inf = bcf_open(VCFFile, "r");
	if (inf == NULL) {
		exit(1);
	}
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Read header and sequences' name
	hdr = bcf_hdr_read(inf);
	int nseq;
	const char **seqnames = bcf_hdr_seqnames(hdr, &nseq);
	if (seqnames == NULL) {
		fputs("Unable to get sequence names in VCF file\n", stderr);
		exit(1);
	}
	int WantedSequenceId = -1;
	{
		char tID[3], tIDext[6];
		snprintf(tID,3, "%i", options.chr);
		snprintf(tIDext, 6, "chr%i", options.chr);
		for (int i = 0; i < nseq; i++) {
			if ((strcmp(tID, seqnames[i]) == 0) || (strcmp(tIDext, seqnames[i]) == 0) ) {
				WantedSequenceId = i;
				break;
			}
		}
	}
	
	if (WantedSequenceId == -1) {
		fprintf(stderr, "Unable to get chromosome %i sequence in VCF file\n", options.chr);
		exit(1);
	}
	
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Loop on VCF file
	unsigned int GarbageCollect = 0U;
	unsigned int AllocatedBlocks = 0U;
	err = 0;
	IntervalJob_t * restrict work;
	unsigned int Count = 0U;
	while (1) {
		
		///////////////////////////////////////////////////////////////////////////////////////////////////////
		// Acquire a free Interval memory slot
		if ( err == 0 ) {
			bsem_wait(IntervalResultsPool->donequeue_p.has_items);
			pthread_mutex_lock(&(IntervalResultsPool->donequeue_p.rwmutex));
			work = (IntervalJob_t *) jobqueue_pull(&(IntervalResultsPool->donequeue_p));
			pthread_mutex_unlock(&(IntervalResultsPool->donequeue_p.rwmutex));
		}
		
		///////////////////////////////////////////////////////////////////////////////////////////////////////
		// Read VCF record
		bcf1_t * const restrict rec = work->VCFRecord;
	Reread:
		if (bcf_read(inf, hdr, rec) < 0) break;
		
		if (rec->rid == WantedSequenceId) {
			///////////////////////////////////////////////////////////////////////////////////////////////////////
			// We try to have all vcf loaded here from file 
			bcf_get_variant_types(rec);
		
			///////////////////////////////////////////////////////////////////////////////////////////////////////
			// Filter on Bed file if given
			if (Filtering) {
				_Bool NotFound = true;
				for (int i=0; i<nBedRegions;i++) {
					if ((rec->pos >= Filtering[i][0]) && (rec->pos <= Filtering[i][1])) {
						NotFound = false;
						break;
					}
				}
				if (NotFound) goto Reread;
			}

			///////////////////////////////////////////////////////////////////////////////////////////////////////
			// Dispatch GTL blocks to threads
			if (! KeepWindow) {
				dispatcherWindow[0] = (rec->pos > Margin) ? rec->pos - Margin : 0;
				dispatcherWindow[1] = rec->pos + Margin;
			}
			else {
				fprintf(stderr, "This part is to be written again!!!\n");
				exit(1);
// 				bcf_info_t * info = bcf_get_info(hdr, rec, "AlleleOK");
// 				if (info == NULL) {
// 						fprintf(stderr, "Unable to retrieve INFO AlleleOK\n");
// 						goto bail;
// 				}
// 				const char * ptr = info->vptr;
// 				unsigned char RangeFound = 0;
// 				while (*ptr != '\0') {
// 						if (ptr[0] == 'W' && ptr[1] == 'i') {
// 								const int cnt = sscanf(ptr, "Window=%u-%u", &dispatcherWindow[0], &dispatcherWindow[1]);
// 								if (cnt == 2) {
// 									RangeFound = 1;
// 									break;
// 								}
// 						}
// 						++ptr;
// 				}
// 				if (RangeFound == 0) {
// 						fprintf(stderr, "Unable to get window in '%s'\n", info->vptr);
// 						goto bail;
// 				}
			}
// 			printf("Dispatching window %u--%u\n", dispatcherWindow[0], dispatcherWindow[1]);
			work->count = 0U;
			work->GenomeStart = dispatcherWindow[0];
			work->GenomeEnd = dispatcherWindow[1];
			const int dispatchCount = dispatchGTLBlockRange(GTLInputFileList, thpool, options.chr,
																											dispatcherWindow[0], dispatcherWindow[1], -1, true, work);
			if (dispatchCount == 0) {
				fprintf(stdout, "\x1b[2K\x1b[GNo coverage at all for location %u\n", rec->pos);
			}
			Count++;
			fprintf(stderr, " %9i/%9i (%lf) | Zygote %9i Location %9i Allele0 %9i Allele1 %9i No variant %9i  Alignment %9i Coverage %9i Fork %9i Allocated Blocks %12u\r",
			        atomic_load(&nRecordsFound), atomic_load(&nVCFRecord), 100.0f*((float) atomic_load(&nRecordsFound))/((float) atomic_load(&nVCFRecord)), atomic_load(&nZygoteError), atomic_load(&nLocationNotFound), atomic_load(&nAllele0NotFound),
			        atomic_load(&nAllele1NotFound), atomic_load(&nNoVariant), atomic_load(&nAlignmentFailure), atomic_load(&nNotEnoughCoverage), atomic_load(&nRemovedHasFork), AllocatedBlocks);
			if ((++GarbageCollect & 0xFFF) == 0) {
				AllocatedBlocks = 0U;
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
					
					if (GTLInputFileList->Blocks[iBlock].whatever != NULL) AllocatedBlocks++;
				}
			}
			if (dispatchCount == 0) goto Reread;
			
		}
		else 
			goto Reread;		
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Wait for jobs to terminate and empty queues
	thpool_wait(thpool);
	thpool_wait(IntervalResultsPool);
	
	destroyGTLThreadPool(thpool);
	destroyIntervalResultsPool(IntervalResultsPool);
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Print last statistics in stderr and stdout 
	fprintf(stderr, " %9i/%9i (%lf) | Zygote %9i Location %9i Allele0 %9i Allele1 %9i No variant %9i Alignment %9i Coverage %9i Fork %9i\n\n",
			        atomic_load(&nRecordsFound), atomic_load(&nVCFRecord), 100.0f*((float) atomic_load(&nRecordsFound))/((float) atomic_load(&nVCFRecord)), atomic_load(&nZygoteError), atomic_load(&nLocationNotFound), atomic_load(&nAllele0NotFound), atomic_load(&nAllele1NotFound), atomic_load(&nNoVariant), atomic_load(&nAlignmentFailure), atomic_load(&nNotEnoughCoverage), atomic_load(&nRemovedHasFork));
	fprintf(stdout, "STATISTICS %9i/%9i (%lf) | Zygote %9i Location %9i Allele0 %9i Allele1 %9i No variant %9i  Alignment %9i Coverage %9i Fork %9i\n\n",
					atomic_load(&nRecordsFound), atomic_load(&nVCFRecord), 100.0f*((float) atomic_load(&nRecordsFound))/((float) atomic_load(&nVCFRecord)), atomic_load(&nZygoteError), atomic_load(&nLocationNotFound), atomic_load(&nAllele0NotFound), atomic_load(&nAllele1NotFound), atomic_load(&nNoVariant), atomic_load(&nAlignmentFailure), atomic_load(&nNotEnoughCoverage), atomic_load(&nRemovedHasFork));

	fprintf(stdout, "Overall read count is %u\n", Count);
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Clean up and quit

	if (seqnames) free(seqnames);

	bcf_hdr_destroy(hdr);


	////////////////////////////////////////////////////////////////////////////////////////////
	// Close VCF file
	bcf_close(inf);

bail:;
	if (GTLInputFileList) destroyGTLList(GTLInputFileList);
	FreeGenome(&Genome);
	
	if (outfile) fclose(outfile);
	if (affinities) free(affinities);
	
	return err;

} /* main */
//---------------------------------------------------------------

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
