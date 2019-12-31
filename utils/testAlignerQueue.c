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
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdatomic.h>
#include <pthread.h>
#include <string.h>
#include <assert.h>
#include <getopt.h>
#include <unistd.h>

#ifdef __GNUC__
#define align aligned
#define __ALWAYS_INLINE __attribute__((always_inline))
#else
#define __ALWAYS_INLINE __attribute__((always_inline, gnu_inline))
#endif

// #define WORKLOAD 512
#define LENGTH_MASK 0xF

/* Generic job */
typedef struct job {
	struct job * prev;
} job_t;

/* Job queue */
typedef struct jobqueue{
	atomic_ulong rwmutex;            /* used for queue r/w access */
	job_t  *front;                   /* pointer to front of queue */
	job_t  *rear;                    /* pointer to rear  of queue */
	int   len;                       /* number of jobs in queue   */
} jobqueue_t;

typedef struct MyJob {
	struct AlignerJob * prev;
	atomic_ulong Treaded;
	int * restrict DataPtr;
	int index;
	atomic_ulong Stored;
} MyJob_t;

typedef struct MySpace {
	struct AlignerJob * prev;
	int * restrict DataPtr;
	atomic_ulong Used;
	int index;
} MySpace_t;

typedef struct AlignerPool {
	pthread_t * threads;                   /* pointer to threads                  */
	jobqueue_t jobqueue_p;                 /* pointer to the job queue            */
	jobqueue_t spacequeue_p;               /* pointer to th espace queue          */
	int num_threads;                       /* total number of threads in the pool */
	atomic_uint num_threads_alive;         /* threads currently alive             */
	volatile int threads_keepalive;
} AlignerPool_t;

atomic_int JobIndex;
atomic_int OutIndex;
atomic_ulong CopiedIndex;
MyJob_t * restrict Jobs = NULL;
MySpace_t * OutputBuffers = NULL;
unsigned long int JobSize;
unsigned long int OutputBufferSize = 0UL;
int nOutputBuffer = 0;
int * CopiedOutData = NULL;

#define MAX_NANOSEC 999999999
#define CEIL(X) ((X-(int)(X)) > 0 ? (int)(X+1) : (int)(X))

static __inline__  void __ALWAYS_INLINE
jobqueue_push(jobqueue_t * const restrict jobqueue_p, job_t * const restrict newjob)
{
	newjob->prev = NULL;
	switch(jobqueue_p->len){
		case 0:  /* if no jobs in queue */
					jobqueue_p->front = newjob;
					jobqueue_p->rear  = newjob;
					break;

		default: /* if jobs in queue */
					jobqueue_p->rear->prev = newjob;
					jobqueue_p->rear = newjob;
	}
	jobqueue_p->len++;
}

static __inline__ job_t* __ALWAYS_INLINE
jobqueue_pull(jobqueue_t * const restrict jobqueue_p)
{
	job_t* job_p = jobqueue_p->front;

	switch(jobqueue_p->len){
		case 0:  /* if no jobs in queue */
					return NULL;
		
		case 1:  /* if one job in queue */
					jobqueue_p->front = NULL;
					jobqueue_p->rear  = NULL;
					break;
		
		default: /* if >1 jobs in queue */
					jobqueue_p->front = job_p->prev;
					
	}
	jobqueue_p->len--;
	return job_p;
}

static __inline__  void __ALWAYS_INLINE
jobqueue_init(jobqueue_t * const restrict jobqueue_p)
{
	memset(jobqueue_p, 0, sizeof(jobqueue_t));
}


/* ========================== THREADPOOL ============================ */

void* thread(AlignerPool_t * const restrict data) 
{
	int CurrentInputIndex;
	const size_t Me = pthread_self();
	unsigned long int cnt = 0UL;
	const int Busy = -2;
	
	////////////////////////////////////////////////////////////////////////////////
	// thread alive
	printf("Activating thread 0x%lx\n", Me); fflush(stdout);
	atomic_fetch_add(&data->num_threads_alive,1);
	
	////////////////////////////////////////////////////////////////////////////////
	// Loop until queue done
	int IGotTheLastJob = 0;
	while (data->threads_keepalive) {
		//////////////////////////////////////////////////////////////////////////////
		// Keep working on that batch until done
		CurrentInputIndex = atomic_load(&JobIndex);
		if (CurrentInputIndex >= 0) {
// 			printf("Thread 0x%lx working on job %i\n", Me, CurrentInputIndex);fflush(stdout);
		Compute:;
			do {
				/* Randomly choose a WorkLoad */
				unsigned long int WorkLoad =  1UL + ( (unsigned long int) rand() & LENGTH_MASK);
				unsigned long int TreatedSoFar = atomic_fetch_add(&(Jobs[CurrentInputIndex].Treaded), WorkLoad);
					
				////////////////////////////////////////////////////////////////////////////
				// Work
				if ((TreatedSoFar) < JobSize) {
					/* Advise other to fetch new batch */
					const int IamLast = ((TreatedSoFar+WorkLoad) >= JobSize);
					if (IamLast) {
						int itmp = CurrentInputIndex;
						WorkLoad = (JobSize-TreatedSoFar);
						// We only change to Empty if JobIndex is still ours 
 						atomic_compare_exchange_strong(&JobIndex, &itmp, -1);
					}
					
					/* Get output buffer */
					unsigned long int SpaceUsedSoFar;
					int CurrentOutputIndex;
					while (1) {
						while ((CurrentOutputIndex = atomic_load(&OutIndex)) < 0) {
							CurrentOutputIndex = -1;
							if (atomic_compare_exchange_strong(&OutIndex,  &CurrentOutputIndex, Busy)) {
								printf("Thread 0x%lx seeking new space\n", Me);fflush(stdout);
								while (1) {
									unsigned long int Test;
									do { Test = 0UL; } while (!atomic_compare_exchange_strong(&(data->spacequeue_p.rwmutex), &Test, -1));
									MySpace_t * const restrict newSpace = (MySpace_t*) jobqueue_pull(&(data->spacequeue_p));
									atomic_store(&(data->spacequeue_p.rwmutex), 0UL);
									if ( newSpace != NULL ) {
										atomic_store(&OutIndex, newSpace->index);
										CurrentOutputIndex = newSpace->index;
										printf("Thread 0x%lx pulled space %i\n", Me, newSpace->index);fflush(stdout);
										break;
									}
								}
							}
						}
						
						assert(WorkLoad > 0 && WorkLoad <= OutputBufferSize);
						
						SpaceUsedSoFar = atomic_fetch_add(&(OutputBuffers[CurrentOutputIndex].Used), WorkLoad);
						if ((SpaceUsedSoFar + WorkLoad) >= OutputBufferSize) {
							/* WARNING:
							 * If still using the current buffer and we are the one just passing over the limit,
							 * advise other to fetch new buffer, then dump data 
							 */
							if (SpaceUsedSoFar < OutputBufferSize) {
								int itmp = CurrentOutputIndex; assert(CurrentOutputIndex >= 0);
								if (atomic_compare_exchange_strong(&OutIndex, &itmp, -1)) {
									printf("Thread 0x%lx not enough space in %i, dumping %lu\n", Me, CurrentOutputIndex, SpaceUsedSoFar);
									unsigned long int dst = atomic_fetch_add(&CopiedIndex, SpaceUsedSoFar);
									
									/* Dump data */
									memcpy(&CopiedOutData[dst], OutputBuffers[CurrentOutputIndex].DataPtr, SpaceUsedSoFar*sizeof(int));
									
									/* Replace the space in queue for later use */
									atomic_store(&OutputBuffers[CurrentOutputIndex].Used, 0UL);
									unsigned long int Test;
									do { Test = 0UL; } while (!atomic_compare_exchange_strong(&(data->spacequeue_p.rwmutex), &Test, -1));
									jobqueue_push(&(data->spacequeue_p), (job_t*) &(OutputBuffers[CurrentOutputIndex]));
									atomic_store(&(data->spacequeue_p.rwmutex), 0UL);
								}
							}
						}
						else 
							break;
					}
					
					/* Perform your partial part */
					assert(CurrentInputIndex >= 0);
					int * const restrict InPtr = Jobs[CurrentInputIndex].DataPtr + TreatedSoFar;
					int * const restrict OutPtr = OutputBuffers[CurrentOutputIndex].DataPtr + SpaceUsedSoFar;
					for (unsigned long int i=0; i<WorkLoad; i++) OutPtr[i] = -InPtr[i]; 
					atomic_fetch_add(&(Jobs[CurrentInputIndex].Stored), WorkLoad);
					cnt += WorkLoad;
					
					if (IamLast) {
						/* Wait for other to terminate */
						while (atomic_load(&Jobs[CurrentInputIndex].Stored) < JobSize) ;
// 						printf("Thread 0x%lx job %i fully stored\n", Me, CurrentInputIndex); fflush(stdout);
						
						/* Ask for new job */
					}
				}
			} while (CurrentInputIndex == atomic_load(&JobIndex));
// 			printf("Thread 0x%lx done with job %i\n", Me, CurrentInputIndex);fflush(stdout);
		}
		
		//////////////////////////////////////////////////////////////////////////////
		// Load next
		int Empty = -1;
		if (atomic_compare_exchange_strong(&JobIndex,  &Empty, Busy)) {
// 			printf("Thread 0x%lx seeking new job\n", Me);fflush(stdout);
			while (1) {
				unsigned long int Test;
				do { Test = 0UL; } while (!atomic_compare_exchange_strong(&(data->jobqueue_p.rwmutex), &Test, -1));
				MyJob_t * const restrict newJob = (MyJob_t*) jobqueue_pull(&(data->jobqueue_p));
				atomic_store(&(data->jobqueue_p.rwmutex), 0UL);
				
				if ( newJob != NULL ) {
					if ( newJob->DataPtr != NULL) {
						atomic_store(&JobIndex, newJob->index);
						CurrentInputIndex = newJob->index;
						printf("Thread 0x%lx pulled job %i\n", Me, newJob->index);fflush(stdout);
	// 					printf("Thread 0x%lx working on job %i\n", Me, CurrentInputIndex);fflush(stdout);
						goto Compute;
					}
					else {
	// 					printf("Thread 0x%lx got terminate job\n", Me);fflush(stdout);
						data->threads_keepalive = 0;
						IGotTheLastJob = 1;
					}
					break;
				}
			}
		}
	}
	
	atomic_fetch_sub(&data->num_threads_alive,1);
	
	/* Dump the partially used output buffers */
	if (IGotTheLastJob) {
		while (atomic_load(&data->num_threads_alive) > 0) {};
		int CurrentOutputIndex = atomic_load(&OutIndex);
		if (CurrentOutputIndex >= 0) {
			unsigned long int SpaceUsedSoFar = atomic_load(&(OutputBuffers[CurrentOutputIndex].Used));
			printf("Thread 0x%lx dumping last used space %i, dumping %lu\n", Me, CurrentOutputIndex, SpaceUsedSoFar);
			unsigned long int dst = atomic_fetch_add(&CopiedIndex, SpaceUsedSoFar);
			memcpy(&CopiedOutData[dst], OutputBuffers[CurrentOutputIndex].DataPtr, SpaceUsedSoFar*sizeof(int));
			atomic_store(&OutputBuffers[CurrentOutputIndex].Used, 0UL);
		}
	}
	
	printf("Ending thread 0x%lx with %lu computations\n", Me, cnt); fflush(stdout);
	return (void*) (uintptr_t) cnt;
}

int compare(const void *a, const void *b)
{
	if (*((int*) a) < *((int*) b))
		return -1;
	else 
		return 1;
}

int main (int argc, char * argv[]) 
{
	AlignerPool_t Pool;
	int nthreads = (int) sysconf(_SC_NPROCESSORS_ONLN);
	int nBatch = 2;
	size_t BatchSize = 0UL;
	int * restrict InputData = NULL;
	
	////////////////////////////////////////////////////////////////////////////////
	// Process arguments
	{
		int c;
		opterr = 0;
		while ((c = getopt (argc, argv, "n:t:s:N:S:h")) != -1) {
			switch (c)
			{
				case 'n':
					nBatch = atoi(optarg);
					break;
				case 's':
					BatchSize = (size_t) atoi(optarg);
					break;
				case 't':
					nthreads = atoi(optarg);
					break;
				case 'N':
					nOutputBuffer = atoi(optarg);
					break;
				case 'S':
					OutputBufferSize = (unsigned long int) atoi(optarg);
					break;
				case 'h':
				default:
					fprintf(stderr, "Usage: %s [-t threads] -n batches -s batch size -N # output buffer -S output buffer size\n", argv[0]);
					goto bail;
			}
		}
		
		if (BatchSize == 0) {
			fputs("Please provide a batch size!\n", stderr);
			goto bail;
		}
		
		if (nOutputBuffer == 0) {
			fputs("Please provide some output buffer!\n", stderr);
			goto bail;
		}
		
		if (OutputBufferSize == 0) {
			fputs("Please provide output buffer size!\n", stderr);
			goto bail;
		}
		
		JobSize = BatchSize;
	}
		
	////////////////////////////////////////////////////////////////////////////////
	// Prepare the data
	InputData = (int*) malloc(2*nBatch*BatchSize*sizeof(int));
	CopiedOutData = (int*) malloc(nBatch*BatchSize*sizeof(int));
	Jobs = (MyJob_t*) malloc(nBatch*sizeof(MyJob_t));
	OutputBuffers = (MySpace_t*) malloc(nOutputBuffer*sizeof(MySpace_t));
	int * const restrict OutputData = (int*) malloc(nOutputBuffer*OutputBufferSize*sizeof(int));
	
	if ( InputData == NULL || Jobs == NULL || OutputBuffers == NULL || OutputData == NULL || CopiedOutData == NULL)  {
		fputs("Memory allocation error\n", stderr);
		goto bail;
	}
	const int TotalSize = (int) (nBatch*BatchSize);
	for (int i=0; i<TotalSize; i++) InputData[i] = i;
		
	////////////////////////////////////////////////////////////////////////////////
	// Set the jobs
	jobqueue_init(&Pool.jobqueue_p);
	for (int i=0; i<nBatch; i++) {
			atomic_store(&(Jobs[i].Stored), 0UL);
			atomic_store(&(Jobs[i].Treaded), 0UL);
			Jobs[i].index = i;
			Jobs[i].DataPtr = InputData + i*JobSize;
	}
	
	////////////////////////////////////////////////////////////////////////////////
	// Set the Space	
	jobqueue_init(&Pool.spacequeue_p);
	for (int i=0; i<nOutputBuffer; i++) {
		OutputBuffers[i].DataPtr = OutputData + i*OutputBufferSize;
		atomic_store(&(OutputBuffers[i].Used), 0UL);
		OutputBuffers[i].index = i;
		jobqueue_push(&Pool.spacequeue_p, (job_t*) &OutputBuffers[i]);
	}
	
	////////////////////////////////////////////////////////////////////////////////
	// Set the pool
  atomic_store(&JobIndex, (int) -1);
	atomic_store(&OutIndex, (int) -1);
	atomic_store(&CopiedIndex, 0UL);
	Pool.threads = (pthread_t*) malloc(nthreads*sizeof(pthread_t));
	if (Pool.threads == NULL) {
		fputs("Memory allocation error\n", stderr);
		goto bail;
	}
	Pool.num_threads = nthreads;
	atomic_store(&Pool.num_threads_alive, 0);
	Pool.threads_keepalive = 1;
	
	const int half = nBatch/2;
	for (int i=0; i<half; i++) jobqueue_push(&(Pool.jobqueue_p), (job_t*) &Jobs[i]);
	
	////////////////////////////////////////////////////////////////////////////////
	// Start threads
	for (int i=0; i<nthreads; i++) {
		pthread_create(&(Pool.threads[i]), NULL, (void* (*)(void*)) thread, (void*) &Pool);
	}
	
	////////////////////////////////////////////////////////////////////////////////
	// In the mean time compute the real solution
	int * const restrict CorrectValue = InputData + nBatch*BatchSize;
	for (int i=0; i<TotalSize; i++) CorrectValue[i] = -(TotalSize-1 - i);
	
	////////////////////////////////////////////////////////////////////////////////
	// Submit the remaining ghalf jobs
	for (int i=half; i<nBatch; i++) jobqueue_push(&(Pool.jobqueue_p), (job_t*) &Jobs[i]);
	MyJob_t EndJob = { .DataPtr = NULL, .prev = NULL, .index = -1 };
	jobqueue_push(&(Pool.jobqueue_p), (job_t*) &EndJob);
	
	////////////////////////////////////////////////////////////////////////////////
	// Wait for threads to terminate
	size_t * const Counts = alloca(sizeof(size_t)*nthreads);
	for (int i=0; i<nthreads; i++) pthread_join(Pool.threads[i], (void**) &Counts[i]);
	printf("Statistics: %lu", Counts[0]);
	size_t sum = Counts[0];
	for (int i=1; i<nthreads; i++) { sum += Counts[i]; printf(" | %lu", Counts[i]);}
	printf(" = %lu\n", sum);
	
	////////////////////////////////////////////////////////////////////////////////
	// Sort results 
	qsort(CopiedOutData, (size_t) TotalSize, sizeof(int), compare);
	
	////////////////////////////////////////////////////////////////////////////////
	// Compare results
	fflush(stdout);
	int good = 0, bad = 0;
	for (int i=0; i<TotalSize; i++) {
		if (CopiedOutData[i] != CorrectValue[i]) {
			fprintf(stderr, "Pos %i: %i <> %i\n", i, CopiedOutData[i], CorrectValue[i]);
			bad++;
		}
		else
			good++;
	}
	printf("%i/%i goods, %i/%i bads\n", good, TotalSize, bad, TotalSize);
	
	////////////////////////////////////////////////////////////////////////////////
	// Clean

	
bail:;
	if (InputData) free(InputData);
	if (Jobs) free(Jobs);
	if (OutputBuffers) free(OutputBuffers);
	if (OutputData) free(OutputData);
	if (CopiedOutData) free(CopiedOutData);
		
	
	return 0;
}
/* vim: tabstop=2 shiftwidth=2
 */
