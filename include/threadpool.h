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
#ifndef _THREAD_POOL_H
#define _THREAD_POOL_H
#include <pthread.h>
#include "jobqueue.h"

/********************************** 
 * @author      Johan Hanssen Seferidis
 * License:     MIT
 * 
 **********************************/

/* ========================== STRUCTURES ============================ */

/* Threadpool */
typedef struct threadpool {
	pthread_t * threads;                   /* pointer to threads                  */
	job_t * jobs;                        	 /* pointer to all jobs memory, only for freeing */
	const volatile void * volatile common; /* place holder for common thread data */
	pthread_mutex_t  thcount_lock;         /* used for thread count etc           */
	jobqueue_t jobqueue_p;                 /* pointer to the job queue            */
	jobqueue_t donequeue_p;                /* pointer to the free jobs queue      */
	int num_threads;                       /* total number of threads in the pool */
	volatile int num_threads_alive;        /* threads currently alive             */
	volatile int num_threads_working;      /* threads currently working           */
	volatile int threads_keepalive;        /* Used to stop threads                */
	volatile int threads_on_hold;
} threadpool_t;


/* =================================== API ======================================= */

void thpool_wait(threadpool_t * const);
void thpool_destroy(threadpool_t * const);
void thpool_pause(threadpool_t * const thpool_p);

static __inline__ void __ALWAYS_INLINE thpool_resume(threadpool_t * const  thpool_p) {
	thpool_p->threads_on_hold = 0;
}

/* ============================ THREAD ============================== */


/* Sets the calling thread on hold */
void thread_hold (threadpool_t * const thpool_p); 

#endif /* _THREAD_POOL_H */

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
