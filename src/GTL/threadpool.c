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
#include <errno.h>
#include <unistd.h>
 
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <pthread.h>

#include "config.h"
#include "threadpool.h" 

/* *****************************************************************************************************
 * Author:       Johan Hanssen Seferidis
 * License:          MIT
 * Description:  Library providing a threading pool where you can add
 *               work. For usage, check the thpool.h file or README.md
 *
 *//** @file threadpool.h *//*
 * 
 ********************************/

#ifdef THREADPOOL_DEBUG
#define THREADPOOL_DEBUG 1
#else
#define THREADPOOL_DEBUG 0
#endif

#define MAX_NANOSEC 999999999
#define CEIL(X) ((X-(int)(X)) > 0 ? (int)(X+1) : (int)(X))

/* ========================== THREADPOOL ============================ */

/* Wait until all jobs have finished */
void thpool_wait(threadpool_t * const thpool_p){

	/* Continuous polling */
	double timeout = 1.0;
	time_t start, end;
	double tpassed = 0.0;
	time (&start);
	while (tpassed < timeout && (thpool_p->jobqueue_p.len || thpool_p->num_threads_working))
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

	while (thpool_p->jobqueue_p.len || thpool_p->num_threads_working)
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
	while (thpool_p->jobqueue_p.len || thpool_p->num_threads_working) sleep(max_secs);
}
//---------------------------------------------------------------

/* Destroy the threadpool */
void thpool_destroy(threadpool_t * const thpool_p){

	//volatile int threads_total = thpool_p->num_threads_alive;

	/* End each thread 's infinite loop */
	thpool_p->threads_keepalive = 0;

	/* Give one second to kill idle threads */
	double TIMEOUT = 1.0;
	time_t start, end;
	double tpassed = 0.0;
	time (&start);
	while (tpassed < TIMEOUT && thpool_p->num_threads_alive){
		bsem_post_all(thpool_p->jobqueue_p.has_items);
		time (&end);
		tpassed = difftime(end,start);
	}

	/* Poll remaining threads */
	while (thpool_p->num_threads_alive){
		bsem_post_all(thpool_p->jobqueue_p.has_items);
		sleep(1);
	}

	/* Job queue cleanup */
	jobqueue_destroy(&(thpool_p->jobqueue_p));
	jobqueue_destroy(&(thpool_p->donequeue_p));
}
//---------------------------------------------------------------

/* Pause threads */
// void thpool_pause(threadpool_t * const thpool_p) {
//         int n;
//         ThreadPool_Arg_t * const restrict threads = thpool_p->threads;
//         for (n=0; n < thpool_p->num_threads_alive; n++){
//                 pthread_kill(threads[n].pthread, SIGUSR1);
//         }
// }

/* ============================ THREAD ============================== */

/* Sets the calling thread on hold */
void thread_hold (threadpool_t * const thpool_p) {
	thpool_p->threads_on_hold = 1;
	while (thpool_p->threads_on_hold) sleep(1);
}
//---------------------------------------------------------------

/* ============================ JOBS ============================== */
