/*
 * ------------------------------------------------------------------------------------------------------------------------
 *
 *                                * micmap *
 *      Mapping of short reads in fastq format onto a reference
 *
 *
 *  Copyright (C) SIB  - Swiss Institute of Bioinformatics,  2015-2019 Nicolas Guex, Thierry Schuepbach and Christian Iseli
 *  Copyright (C) UNIL - University of Lausanne, Switzerland      2019 Nicolas Guex and Christian Iseli
 *  Copyright (C) EPFL - EPFL, Lausanne, Switzerland              2022 Christian Iseli
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
#ifndef _CPU_THREAD_ALIGNER_POOL_H
#define _CPU_THREAD_ALIGNER_POOL_H
#include <stdatomic.h>
#include <smmintrin.h>
#include "jobqueue.h"
#include "Reader.h"
#include "Decoder.h"
#include "Aligner.h"

//--------------------------------------------------------------- 
// DEFINITIONS 
//---------------------------------------------------------------

//--------------------------------------------------------------- 
// STRUCTURE DEFINITIONS 
//--------------------------------------------------------------- 
/* Generic job */
typedef struct Task {
	struct Task * prev;
} Task_t;

/* Job queue */
typedef struct Queue {
	atomic_int rwmutex;    /* used for queue r/w access */
	Task_t * front;        /* pointer to front of queue */
	Task_t * rear;         /* pointer to rear  of queue */
	int   len;             /* number of jobs in queue   */
} Queue_t;

typedef struct cpuAlignerJob {
	struct cpuAlignerJob * prev;
	atomic_ulong Treated;
	const ReaderTag_t * restrict TagsData;
	const DecoderElement_t (* restrict Anchors)[kDesiredHitCnt];
	unsigned int TagCnt;
	int index;
	atomic_ulong Stored;
	DecoderJobElement_t * DecoderJob; 
} cpuAlignerJob_t;

/* 
 * The following structure must be compatible with the job received 
 * directly by the writers, that is AlignerOutputBlock_t
 */
typedef struct cpuSpace {
	struct cpuSpace * prev;
	char * restrict OutBufferPtr;
	memqueue_t * BelongsTo_q;
	size_t length;
	atomic_ulong Used;
	atomic_ulong Stored;
	unsigned int index;
} cpuSpace_t;

enum PROCESS_STATS_TYPE {
	LargeAlignment=0,
	LargeAlignmentRescued=1,
	Alignment=2,
	AlignmentRescued=3,
	HalfQuickAlign=4,
	HalfAlignment=5,
	FullAlignment=6,
	BestOfEach=7,
	HalfMapped=8,
	ForCompression=9,
	directhit=10,
	hopeless=11,
	quickalign=12,
	mismatch_overflow=13,
	Alternate=14
};

typedef union Statistics {
	unsigned int Counter[16];
	__m128i Vector[4];
} Statistics_t;

typedef struct cpuAlignerPool {
	AlignerBlockArgs_t * restrict AlignerBlock; /* pointer to parent aligner block                         */
	pthread_t * threads;                        /* pointer to threads                                      */
	union Statistics * restrict stats;          /* pointer to threads statistics                           */
	Queue_t job_q;                              /* pointer to the job queue                                */
	int num_threads;                            /* total number of threads in the pool                     */
	jobqueue_t space_q;                         /* pointer to th espace queue                              */
	volatile int threads_keepalive;             /* workers use this to trigger self ending                 */
	atomic_int num_threads_alive;               /* workers use this to get their ID in statistics          */
	atomic_int JobIndex;
	atomic_int OutIndex;
	cpuAlignerJob_t * restrict Jobs;
	cpuSpace_t * OutputBuffers;
} cpuAlignerPool_t;

//--------------------------------------------------------------- 
// GLOBALS 
//---------------------------------------------------------------


//--------------------------------------------------------------- 
// FUNCTIONS
//--------------------------------------------------------------- 
void* avx2_pair(cpuAlignerPool_t * const restrict data);
void* sse41_pair(cpuAlignerPool_t * const restrict data);
void* avx512f_pair(cpuAlignerPool_t * const restrict data);

//--------------------------------------------------------------- 
// INLINE FUNCTIONS
//--------------------------------------------------------------- 
static __inline__  void __ALWAYS_INLINE
pushToQueue(Queue_t * const restrict queue, Task_t * const restrict newjob)
{
	newjob->prev = NULL;
	switch(queue->len){
		case 0:  /* if no jobs in queue */
					queue->front = newjob;
					queue->rear  = newjob;
					break;

		default: /* if jobs in queue */
					queue->rear->prev = newjob;
					queue->rear = newjob;
	}
	queue->len++;
}

static __inline__ Task_t* __ALWAYS_INLINE
pullFromQueue(Queue_t * const restrict queue)
{
	Task_t* job_p = queue->front;

	switch(queue->len){
		case 0:  /* if no jobs in queue */
					return NULL;
		
		case 1:  /* if one job in queue */
					queue->front = NULL;
					queue->rear  = NULL;
					break;
		
		default: /* if >1 jobs in queue */
					queue->front = job_p->prev;
					
	}
	queue->len--;
	return job_p;
}

#endif /* _CPU_THREAD_ALIGNER_POOL_H */
/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
