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
#ifndef _JOB_QUEUE_H
#define _JOB_QUEUE_H
#include <assert.h>
#include <pthread.h>

#ifdef NO_JOBQUEUE_INLINE
#define INLINE_JOBQUEUE(type) type
#else
#define INLINE_JOBQUEUE(type) __inline__  type __ALWAYS_INLINE
#endif

/********************************** 
 * @author      Johan Hanssen Seferidis
 * License:     MIT
 * 
 **********************************/

/* ========================== STRUCTURES ============================ */

/* Binary semaphore */
typedef struct bsem {
	pthread_mutex_t mutex;
	pthread_cond_t   cond;
	int v;
} bsem;

/* Generic job */
typedef struct job {
	struct job * prev;
} job_t;

/* Job queue */
typedef struct jobqueue{
	pthread_mutex_t rwmutex;            /* used for queue r/w access */
	job_t  *front;                      /* pointer to front of queue */
	job_t  *rear;                       /* pointer to rear  of queue */
	bsem *has_items;                    /* flag as binary semaphore  */
	int   len;                          /* number of jobs in queue   */
} jobqueue_t;

/* ======================== SYNCHRONISATION ========================= */

/* Init semaphore to 1 or 0 */
static INLINE_JOBQUEUE(void) bsem_init(bsem *bsem_p, int value) {
	if (value < 0 || value > 1) {
		fprintf(stderr, "bsem_init(): Binary semaphore can take only values 1 or 0");
		exit(1);
	}
	pthread_mutex_init(&(bsem_p->mutex), NULL);
	pthread_cond_init(&(bsem_p->cond), NULL);
	bsem_p->v = value;
}


/* Reset semaphore to 0 */
/* FIXME - not sure it is a good idea to re-init the structures here */
static INLINE_JOBQUEUE(void) bsem_reset(bsem *bsem_p) {
	bsem_init(bsem_p, 0);
}


/* Post to at least one thread */
static INLINE_JOBQUEUE(void) bsem_post(bsem *bsem_p) {
	assert(bsem_p != NULL);
	pthread_mutex_lock(&bsem_p->mutex);
	bsem_p->v = 1;
	pthread_cond_signal(&bsem_p->cond);
	pthread_mutex_unlock(&bsem_p->mutex);
}


/* Post to all threads */
static INLINE_JOBQUEUE(void) bsem_post_all(bsem *bsem_p) {
	pthread_mutex_lock(&bsem_p->mutex);
	bsem_p->v = 1;
	pthread_cond_broadcast(&bsem_p->cond);
	pthread_mutex_unlock(&bsem_p->mutex);
}


/* Wait on semaphore until semaphore has value 0 */
static INLINE_JOBQUEUE(void) bsem_wait(bsem* bsem_p) {
	pthread_mutex_lock(&bsem_p->mutex);
	while (bsem_p->v != 1) {
		pthread_cond_wait(&bsem_p->cond, &bsem_p->mutex);
	}
	bsem_p->v = 0;
	pthread_mutex_unlock(&bsem_p->mutex);
}

/* ============================ JOB QUEUE =========================== */

/* Add (allocated) job to queue
 *
 * Notice: Caller MUST hold a mutex
 */
static INLINE_JOBQUEUE(void)
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
	assert(jobqueue_p->has_items != NULL);
	bsem_post(jobqueue_p->has_items);
}

#if 0
/*
 * Notice: No NEED for Caller to hold a mutex
 */
static INLINE_JOBQUEUE(void)
jobqueue_push_batch(jobqueue_t * const restrict jobqueue_p, job_t * const restrict newjobs, const size_t N)
{
//	job_t * ljob         = &newjobs[N-1];
	job_t * const newjob = &newjobs[N-1];
// 	while ( (uintptr_t) ljob > (uintptr_t) newjobs) {
// 		ljob->prev = --ljob;
// 	}
// 	ljob->prev = NULL;
	for (int i=N-1; i>0; i--) newjobs[i].prev = &newjobs[i-1];
	newjobs[0].prev = NULL;
	
	pthread_mutex_lock(&jobqueue_p->rwmutex);
	switch(jobqueue_p->len){

		case 0:  /* if no jobs in queue */
					jobqueue_p->front = newjob;
					jobqueue_p->rear  = &newjobs[0];
					break;

		default: /* if jobs in queue */
					jobqueue_p->rear->prev = newjob;
					jobqueue_p->rear = &newjobs[0];
	}
	jobqueue_p->len += (int) N;
	
	assert(jobqueue_p->has_items != NULL);
	bsem_post(jobqueue_p->has_items);
	pthread_mutex_unlock(&jobqueue_p->rwmutex);
}
#endif

/* Get first job from queue(removes it from queue)
 * 
 * Notice: Caller MUST hold a mutex
 */
static INLINE_JOBQUEUE(job_t*)
jobqueue_pull(jobqueue_t * const restrict jobqueue_p)
{
	assert(jobqueue_p);
	job_t* job_p = jobqueue_p->front;
	assert(jobqueue_p->len >= 0);
	switch(jobqueue_p->len){
		case 0:  /* if no jobs in queue */
					return NULL;
		
		case 1:  /* if one job in queue */
					jobqueue_p->front = NULL;
					jobqueue_p->rear  = NULL;
					break;
		
		default: /* if >1 jobs in queue */
#ifndef NDEBUG
					if (job_p == NULL || job_p->prev == NULL)
					  fprintf(stderr,"len : %d; front : %p; rear : %p; job_p->prev : %p\n",jobqueue_p->len,jobqueue_p->front,jobqueue_p->rear,job_p->prev);
#endif
					assert(job_p != NULL);
					assert(job_p->prev != NULL);
					jobqueue_p->front = job_p->prev;
					job_p->prev = NULL;
	}
	jobqueue_p->len--;
	
	/* Make sure has_items has right value */
	if (jobqueue_p->len > 0) {
		bsem_post(jobqueue_p->has_items);
	}

	return job_p;
}

/* Clear the queue */
static INLINE_JOBQUEUE(void)
jobqueue_clear(jobqueue_t * const restrict jobqueue_p)
{
// 	while(jobqueue_p->len){
// 		free(jobqueue_pull(jobqueue_p));
// 	}

	jobqueue_p->front = NULL;
	jobqueue_p->rear  = NULL;
	bsem_reset(jobqueue_p->has_items);
	jobqueue_p->len = 0;
}

/* Initialize queue */
static INLINE_JOBQUEUE(int)
jobqueue_init(jobqueue_t * const restrict jobqueue_p)
{
	if (jobqueue_p == NULL){
		return -1;
	}
	pthread_mutex_init(&(jobqueue_p->rwmutex), NULL);
	
	jobqueue_p->has_items = (struct bsem*)malloc(sizeof(struct bsem));
	if (jobqueue_p->has_items == NULL) {
		return -1;
	}
	bsem_init(jobqueue_p->has_items, 0);
	jobqueue_p->len = 0;
	jobqueue_clear(jobqueue_p);
	return 0;
}


/* Free all queue resources back to the system */
static INLINE_JOBQUEUE(void) jobqueue_destroy(jobqueue_t * const restrict jobqueue_p){
	jobqueue_clear(jobqueue_p);
	free(jobqueue_p->has_items);
}


#endif
/* vim: tabstop=2 shiftwidth=2
 */
