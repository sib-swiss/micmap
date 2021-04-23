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
#include "htslib/vcf.h"
#include "htslib/vcfutils.h"

#include "Genome.h"
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
// #define kNOT_ACGT 0xff

//=================================================================================================
#include "AnalyzeTags.h"

#ifdef DEBUG_INFO
unsigned int longuestAllowedDeletion = 2000;
#endif

//=================================================================================================
static threadpool_t * restrict IntervalResultsPool = NULL;
static int verbose = 0;
static int debug = 0;
static char default_genome[] = "/tmp/nodelete/nguex/data/hg19.bin";

//=================================================================================================
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

const char * restrict genomefile = default_genome;

Genome_t Genome = {
	.FileName = &default_genome[0],
	.table = NULL,
	.configBuffer = NULL,
	.tableSize = 0UL,
	.configBufferSize = 0UL
};

//=================================================================================================
void* AnalyzeTags(threadpool_t * restrict const thpool);
void* AnalyzeTags3(threadpool_t * restrict const thpool);
void* AnalyzeTags3_FullOutput(threadpool_t * restrict const thpool);

//=================================================================================================
static void __attribute__((noreturn)) usage()
{
		printf("GTLcaller3 [options] GTL files\n\n");
		printf("           -g genomefile         : full path of the genome encoded by encode_genome\n");
		printf("                                   default: '%s' \n",genomefile);
		printf("           -C chromosome         : select chromosome to analyze\n");
		printf("           -c                    : do not filter clones (default is to filter them out)\n");
		printf("           -f <char>             : filter out nucleotides where the quality value is *below* the specified character [%c]\n", options.filterQualValue + '#');
		printf("           -m <int>              : minor allele minimum percent coverage\n");
		printf("           -v level              : verbose level\n");
		printf("           -T <uint>/<uint>      : number of thread to fetch tags/analyze tags\n");
		printf("           -w <uint>             : dispatcher window width\n");
		printf("           -s <file>             : selected regions to compute coverage (bed-style)\n");
		printf("           -V <file>             : VCF input file to recompute\n");
		printf("           -2                    : perform 2 pass with profile generation\n");
		printf("           -F                    : extended output information in VCF\n");
		printf("           -W                    : use WND information in VCF file\n");
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
					unsigned int prevPos1 = 0;
					unsigned int prevPos2 = 0;
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
							// Filter duplicates...
							int keep = 1;
							if (!options.noCloneFilter && ((ND[0]->Location == prevPos1 && ND[1]->Location == prevPos2) || (ND[0]->Location == prevPos2 && ND[1]->Location == prevPos1)))
							{
								keep = 0;
								//fprintf(stderr, "Filtering %u %u %u %u %u %u\n", ND[0]->Location, ND[0]->AlignmentRange[0], ND[0]->AlignmentRange[1], ND[1]->Location, ND[1]->AlignmentRange[0], ND[1]->AlignmentRange[1]);
							}
							prevPos1 = ND[0]->Location;
							prevPos2 = ND[1]->Location;

							for (int iTag=0; iTag<2; iTag++) {
								GTLRawData_t * const restrict GD = ND[iTag];
								const unsigned int TagGenLen   = GD->AlignmentRange[1] - GD->AlignmentRange[0];
								const unsigned int TagGenBegin = GD->Location;
								const unsigned int TagGenEnd   = GD->Location + TagGenLen;

								if (keep && TagGenBegin <= GenomeWindowStop && TagGenEnd >= GenomeWindowStart) {
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
					unsigned int prevPos1 = 0;
					unsigned int prevPos2 = 0;
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
							// Filter duplicates...
							int keep = 1;
							if (!options.noCloneFilter && ((ND[0]->Location == prevPos1 && ND[1]->Location == prevPos2) || (ND[0]->Location == prevPos2 && ND[1]->Location == prevPos1)))
							{
								keep = 0;
								//fprintf(stderr, "Filtering %u %u %u %u %u %u\n", ND[0]->Location, ND[0]->AlignmentRange[0], ND[0]->AlignmentRange[1], ND[1]->Location, ND[1]->AlignmentRange[0], ND[1]->AlignmentRange[1]);
							}
							prevPos1 = ND[0]->Location;
							prevPos2 = ND[1]->Location;

							for (int iTag=0; iTag<2; iTag++) {
								GTLRawData_t * const restrict GD = ND[iTag];
								const unsigned int TagGenLen   = GD->AlignmentRange[1] - GD->AlignmentRange[0];
								const unsigned int TagGenBegin = GD->Location;
								const unsigned int TagGenEnd   = GD->Location + TagGenLen;

								if (keep && TagGenBegin <= GenomeWindowStop && TagGenEnd >= GenomeWindowStart) {
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

static size_t LoadBedFile(const char * const restrict BedFileName, unsigned int (** const Filtering)[2])
{
	FILE * const restrict bed = fopen(BedFileName, "r");
	size_t nBedRegions = 0UL;
	if (bed) {
		char *Buffer = NULL;
		size_t Buffer_size = 0UL;
		size_t lines = 0UL;
		size_t GoodLines = 0UL;
		char ChrString[32];
		unsigned int range[2];
		while (!feof(bed)) {
			const ssize_t nRead = getline(&Buffer, &Buffer_size, bed);
			if (nRead == -1) {
				if (feof(bed)) break;
				perror("Getline");
				return 0UL;
			}
			int count = sscanf(Buffer, "%31s\t%u\t%u", ChrString, &range[0], &range[1]);
			if ( count != 3) {
				if ( GoodLines == 0UL) {
					lines++;
					continue;
				}
				else {
					fprintf(stderr, "Unable to read bed file at line %lu, read only %i/3\n Line is : %s", lines+1, count, Buffer);
					fclose(bed);
					return 0UL;
				}
			}
			const char * ChrPtr = ChrString;
			if (ChrString[0] == 'c' || ChrString[0] == 'C') ChrPtr += 3;
			int chr = 0;
			if (ChrPtr[0] == 'x' || ChrPtr[0] == 'X') chr = 23;
			else if (ChrPtr[0] == 'y' || ChrPtr[0] == 'Y') chr = 24;
			else if (ChrPtr[0] == 'm' || ChrPtr[0] == 'M') chr = 25;
			if (chr == 0) chr = atoi(ChrPtr);
			if (chr == options.chr) GoodLines++;
			lines++;
		}
		if (verbose) printf("Bed file %s holds %lu lines\n", BedFileName, lines);
		*Filtering = (unsigned int (*)[2]) malloc(2*GoodLines*sizeof(unsigned int));
		if (*Filtering == NULL) {
			fputs("Unable to allocate memory for bed filtering\n", stderr);
			fclose(bed);
			exit(1);
		}
		rewind(bed);
		nBedRegions = GoodLines;
		lines = 0UL;
		GoodLines = 0UL;
		while (!feof(bed)) {
			const ssize_t nRead = getline(&Buffer, &Buffer_size, bed);
			if (nRead == -1) {
				if (feof(bed)) break;
				perror("Getline");
				return 0UL;
			}
			int count = sscanf(Buffer, "%31s\t%u\t%u", ChrString, &range[0], &range[1]);
			if ( count != 3) {
				if ( GoodLines == 0UL) {
					lines++;
					continue;
				}
				else {
					fprintf(stderr, "Unable to read bed file at line %lu, read only %i/3\n Line is : %s", lines+1, count, Buffer);
					fclose(bed);
					return 0UL;
				}
			}
			const char * ChrPtr = ChrString;
			if (ChrString[0] == 'c' || ChrString[0] == 'C') ChrPtr += 3;
			int chr = 0;
			if (ChrPtr[0] == 'x' || ChrPtr[0] == 'X') chr = 23;
			else if (ChrPtr[0] == 'y' || ChrPtr[0] == 'Y') chr = 24;
			else if (ChrPtr[0] == 'm' || ChrPtr[0] == 'M') chr = 25;
			if (chr == 0) chr = atoi(ChrPtr);
			if (chr == options.chr) {
				if (GoodLines>0 && range[0] > 20U) {
					const unsigned int uitmp = range[0] - 20U;
					if (uitmp >= (*Filtering)[GoodLines-1][1]) range[0] = uitmp;
				}
				(*Filtering)[GoodLines][0] = range[0]; (*Filtering)[GoodLines][1] = range[1];
				if (++GoodLines > nBedRegions) { fprintf(stderr, "going beyond control at line %lu\n" , lines);
					exit(1);
				}
				lines++;
			}
		}
		fclose(bed);
		if (Buffer) free(Buffer);
	}
	return nBedRegions;
}
//---------------------------------------------------------------

//=================================================================================================
int main (int argc, char **argv)
{
	const char * restrict rsltfile = NULL;
	const char * restrict VCFFile = NULL;
	const char * * InputGTLFiles = NULL;
	GTLList_t * restrict GTLInputFileList = NULL;
	threadpool_t * restrict thpool = NULL;
	void* (*AnalysisFct)(threadpool_t * restrict const) = AnalyzeTags3;
	cpu_set_t * affinities = NULL;

	int c, err = 1;
	unsigned int dispatcherWindow[2] = { 0U, 0U};
	unsigned int WindowsWidth = 0U;
	unsigned int nFetchThreads = 4U;
	unsigned int nAnalizeThreads = (unsigned int) sysconf(_SC_NPROCESSORS_ONLN) - 4U;
	unsigned int UseAffinities = 0U;

	unsigned int nInputGTLFiles = 0U;
	unsigned int noShortCut = 0U;
	_Bool DeprecatedAnalysis = false;
	_Bool FullOutputMode = false;
	_Bool KeepVCFWindow = false;

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Process arguments

	opterr = 0;
	while ((c = getopt (argc, argv, "g:cC:v:f:m:T:w:s:V:ah2zdFW")) != -1)
	switch (c)
	{
		case 'c':
			options.noCloneFilter = 1;
			break;

		case 'C':
			sscanf(optarg,"%d",&options.chr);
			break;

		case 'g':
			Genome.FileName = optarg;
			break;

		case 'f':
			options.filterQualValue = (unsigned int) *optarg;
			break;

		case 'm':
			options.minorAllelePct = abs(atoi(optarg));
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

		case 'V':
			VCFFile = optarg;
			break;

		case 'a':
			UseAffinities = 1U;
			break;

		case '2':
		    PerformTwoPass = true;
		    break;

		case 'z':
				noShortCut = 1U;
				break;

		case 'd':
				DeprecatedAnalysis = true;
				break;

		case 'F':
				FullOutputMode = true;
				break;

		case 'W':
			KeepVCFWindow = true;
			break;

		case 'h':
			usage();
			break;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Check Consistency
	if (argc == 1) usage();

	if (FullOutputMode && DeprecatedAnalysis) {
		fprintf(stderr, "Extended format with deprecated analysis method is not consistent\n");
		exit(1);
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Check requirements
	if (!(VCFFile && KeepVCFWindow) && (WindowsWidth == 0U)) {
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

	if (verbose) fprintf(stderr, "Will use %u threads for fetching data, %u for analyzing\n", nFetchThreads, nAnalizeThreads);

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
// 	thpool->common = (volatile void * volatile) &dispatcherWindow[0]; // Check this out, I do not think it is necessary anymore

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Prepare Interval Results pool
	if (FullOutputMode)
		AnalysisFct = &AnalyzeTags3_FullOutput;
	else if (DeprecatedAnalysis)
		AnalysisFct = &AnalyzeTags;

	IntervalResultsPool = createIntervalResultsPool(AnalysisFct, affinities ? &affinities[nFetchThreads] : NULL, nAnalizeThreads, 0);
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
	// VCF Header
	printf("##fileformat=VCFv4.2\n");
	{
		struct tm t;
		const time_t now = time(NULL);
		struct tm * dummy = gmtime_r(&now, &t);
		printf("##fileDate=%4.4i%2.2i%2.2i\n", 1900+t.tm_year, 1+t.tm_mon, t.tm_mday);
	}
	// 	if (opt->selectRegions[0] != 0)
	printf("##SAMPLE=<ID=CALLER3,CONTIG=1,IQR25=0,MEDIAN=0,IQR75=0,CNT=0,SEL=0>\n"
	       "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth\">\n"
	       "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
	       "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Per-sample read depths for each allele\">\n"

	       "##INFO=<ID=WA,Number=1,Type=Integer,Description=\"Window Allele Number\">\n"
	       "##INFO=<ID=WAF,Number=.,Type=Float,Description=\"Window Allele Frequencies\">\n"
	       "##INFO=<ID=WREG,Number=1,Type=Integer,Description=\"Number of regions within window\">\n"
	       "##INFO=<ID=WND,Number=1,Type=String,Description=\"Caller window\">\n"
	       "##INFO=<ID=WTRASH,Number=1,Type=Float,Description=\"Window Trash, Not Assigned\">\n"
	       "##INFO=<ID=REGWRN,Number=0,Type=Flag,Description=\"Region not covered\">\n");
	if (FullOutputMode) {
		printf("##INFO=<ID=ORI,Number=.,Type=Integer,Description=\"Allele orientation coverage\">\n"
		       "##INFO=<ID=GEN,Number=1,Type=String,Description=\"Genome region\">\n"
		       "##INFO=<ID=PHRED,Number=.,Type=Float,Description=\"PHRED accounting for values 0-5,6-9,10-19,20-max\">\n"
		       "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth\">\n"
	         "##INFO=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
	         "##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Per-sample read depths for each allele\">\n");
	}
	for (unsigned int c = 1; c <= chrMax; c += 1)
		printf("##contig=<ID=%s,length=%u,assembly=%s>\n",Genome.virtchr[c].AC,Genome.virtchr[c].len,options.refName);
	printf("##phasing=partial\n"
	       "##source=GTL caller 3 - " GTL_VERSION_FULL_STRING "\n"
	       "##variants_justified=left\n"
	       "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCALLER3\n");

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// NO BED file
	if (options.selectRegions[0] == '\0') {
		////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Do we have a VCF file
		if (VCFFile) {
			/////////////////////////////////////////////////////////////////////////////////////////////////////////
			// Open VCF file
			htsFile * inf = bcf_open(VCFFile, "r");
			if (inf == NULL) {
				exit(1);
			}

			/////////////////////////////////////////////////////////////////////////////////////////////////////////
			// Read header and sequences' name
			bcf_hdr_t * hdr = bcf_hdr_read(inf);
			int nseq;
			const char **seqnames = bcf_hdr_seqnames(hdr, &nseq);
			if (seqnames == NULL) {
				fputs("Unable to get sequence names in VCF file\n", stderr);
				exit(1);
			}
			int WantedSequenceId = -1;
			{
				char tID[3], tIDext[6];
				if (options.chr <= 22) {
					snprintf(tID,3, "%i", options.chr);
					snprintf(tIDext, 6, "chr%i", options.chr);
				}
				else if (options.chr == 23) {
					tID[0] = 'X'; tID[1] = '\0';
					snprintf(tIDext, 6, "chr%c", 'X');
				}
				else if (options.chr == 24) {
					tID[0] = 'Y'; tID[1] = '\0';
					snprintf(tIDext, 6, "chr%c", 'Y');
				}
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

			bcf1_t * const restrict rec = bcf_init();
			unsigned int GarbageCollect = 0U;
			unsigned int dispatchCount = 1;
			IntervalJob_t * restrict work;
			const int Margin = WindowsWidth/2;
			unsigned int PreviousWindow[2] = { 0U, 0U };
			while (1) {
				///////////////////////////////////////////////////////////////////////////////////////////////////////
				// Acquire a free Interval memory slot
				if ( dispatchCount > 0 ) {
					bsem_wait(IntervalResultsPool->donequeue_p.has_items);
					pthread_mutex_lock(&(IntervalResultsPool->donequeue_p.rwmutex));
					work = (IntervalJob_t *) jobqueue_pull(&(IntervalResultsPool->donequeue_p));
					pthread_mutex_unlock(&(IntervalResultsPool->donequeue_p.rwmutex));
				}

				///////////////////////////////////////////////////////////////////////////////////////////////////////
				// Read VCF record
				if (bcf_read(inf, hdr, rec) < 0) break;
				if (rec->rid == WantedSequenceId) {
					///////////////////////////////////////////////////////////////////////////////////////////////////////
					// We try to have all vcf loaded here from file
					bcf_get_variant_types(rec);

					///////////////////////////////////////////////////////////////////////////////////////////////////////
					// Dispatch GTL blocks to threads
					if (! KeepVCFWindow) {
						dispatcherWindow[0] = (rec->pos > Margin) ? rec->pos - Margin : 0;
						dispatcherWindow[1] = rec->pos + Margin;
					}
					else {
						bcf_info_t * info = bcf_get_info(hdr, rec, "WND");
						if (info == NULL) {
								fprintf(stderr, "Unable to retrieve in INFO field WND\n");
								goto bail;
						}
						const int cnt = sscanf((char *) info->vptr, "%u-%u", &dispatcherWindow[0], &dispatcherWindow[1]);
						if (cnt != 2) {
							fprintf(stderr, "Unable to get window in '%s'\n", info->vptr);
							goto bail;
						}
					}
// 		 			printf("Dispatching window %u-%u\n", dispatcherWindow[0], dispatcherWindow[1]);
					if (PreviousWindow[0] != dispatcherWindow[0] && PreviousWindow[1] != dispatcherWindow[1]) {
						work->count = 0U;
						work->GenomeStart = dispatcherWindow[0];
						work->GenomeEnd = dispatcherWindow[1];
						dispatchCount = dispatchGTLBlockRange(GTLInputFileList, thpool, options.chr,
																									dispatcherWindow[0], dispatcherWindow[1], -1, true, work);
						if (dispatchCount == 0) {
							fprintf(stderr, "No coverage at all for location %u\n", rec->pos);
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
							GarbageCollect = 0U;
						}
						PreviousWindow[0] = dispatcherWindow[0]; PreviousWindow[1] = dispatcherWindow[1];
					}
					else
						dispatchCount = 0;
				}
			}
			bcf_destroy(rec);
		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////
		// No VCF file
		else {
			//////////////////////////////////////////////////////////////////////////////////////////////////////
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

			///////////////////////////////////////////////////////////////////////////////////////////////////////
			// Loop on genome regions
			unsigned int DispatchCount = 1U;
			unsigned int GarbageCollect = 0U;
			IntervalJob_t * restrict work;
			const unsigned int InputChrRangeBegin = ChrRangeBegin;
			const unsigned int InputChrRangeEnd = ChrRangeEnd;
			const float invScale = 100.0f/((float) InputChrRangeEnd);
			do {
				///////////////////////////////////////////////////////////////////////////////////////////////////////
				// Acquire a free Interval memory slot
				if (DispatchCount > 0) {
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
				work->GenomeEnd = ChrRangeBegin;
				if (work->GenomeEnd > ChrRangeEnd) work->GenomeEnd = ChrRangeEnd;
				//fprintf(stderr, "Submitting window %u--%u\n", work->GenomeStart, work->GenomeEnd);
				DispatchCount = dispatchGTLBlockRange(GTLInputFileList, thpool, options.chr, work->GenomeStart, work->GenomeEnd, -1, true, work);

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
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// BED file
	else {
		unsigned int (*Filtering)[2] = NULL;
		IntervalJob_t * restrict work;
		const size_t nBedFileLines = LoadBedFile(options.selectRegions, &Filtering);
		unsigned int DispatchCount = 1U;
		unsigned int GarbageCollect = 0U;

		for (size_t iLine=0UL; iLine < nBedFileLines; iLine++) {
			unsigned int ChrRangeBegin = Filtering[iLine][0];
			do {
				//////////////////////////////////////////////////////////////////////////////////////////////
				// Acquire a free Interval memory slot
				if (DispatchCount > 0) {
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
				work->GenomeEnd = ChrRangeBegin;
				if (work->GenomeEnd > Filtering[iLine][1]) work->GenomeEnd = Filtering[iLine][1];
				DispatchCount = dispatchGTLBlockRange(GTLInputFileList, thpool, options.chr, work->GenomeStart, work->GenomeEnd, -1, true, work);

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
					GarbageCollect = 0U;
				}
			} while (ChrRangeBegin < Filtering[iLine][1]);
		}
	}
	err = 0;
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Wait for jobs to terminate and empty queues
	thpool_wait(thpool);
	thpool_wait(IntervalResultsPool);

	destroyGTLThreadPool(thpool);
	destroyIntervalResultsPool(IntervalResultsPool);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Clean up and quit

bail:;
	fflush(stdout);
	if (GTLInputFileList) destroyGTLList(GTLInputFileList);

	FreeGenome(&Genome);

	if (affinities) free(affinities);

	return err;

} /* main */
//---------------------------------------------------------------

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */

