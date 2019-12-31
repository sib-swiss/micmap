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
#include "constants.h"
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <unistd.h>

#include "global_var.h"
#include "Aligner.h"

//--------------------------------------------------------------- 
// DEFINITIONS 
//--------------------------------------------------------------- 
#define MAX_NANOSEC 999999999
#define CEIL(X) ((X-(int)(X)) > 0 ? (int)(X+1) : (int)(X))

//--------------------------------------------------------------- 
// FCTS 
//--------------------------------------------------------------- 
AlignerPool_t* createAlignerPool(jobqueue_t * const restrict ToBeAligned_q)
{
	AlignerPool_t * const restrict data = (AlignerPool_t*) malloc(sizeof(AlignerPool_t));
	if (data == NULL) goto err1;
	
	if (jobqueue_init(&(data->threads_q)) == -1) goto err2;
	if (jobqueue_init(&(data->ToBeWritten_q)) == -1) goto err3;
	if (pthread_mutex_init(&(data->thcount_lock), NULL) != 0) goto err4;
	
	data->num_threads_alive = 0;
	data->num_threads_working = 0;
	data->threads_keepalive = 1;
	data->threads_on_hold = 0;
	data->ToBeAligned_q = ToBeAligned_q;
	
	initRevCompTable(gIUPACrevcomp);

	return data;
	err4: ;
		jobqueue_destroy(&(data->ToBeWritten_q));
	err3:;
		jobqueue_destroy(&(data->threads_q));
	err2:;
		free(data);
	err1: ;
		return NULL;
}
//---------------------------------------------------------------

void waitForAlignerPool(AlignerPool_t * const restrict thpool)
{
	/* Continuous polling */
	double timeout = 1.0;
	time_t start, end;
	double tpassed = 0.0;
	time (&start);
	while ((tpassed < timeout) && ((thpool->ToBeAligned_q->len> 0) || (thpool->num_threads_working > 0))) {
		time (&end);
		tpassed = difftime(end,start);
		//printf("Aligner queue : len=%i, workers=%i\n", thpool->ToBeAligned_q->len, thpool->num_threads_working);
	}

	/* Exponential polling */
	long init_nano = 1; /* MUST be above 0 */
	long new_nano;
	double multiplier = 1.01;
	int  max_secs = 20;

	struct timespec polling_interval;
	polling_interval.tv_sec  = 0;
	polling_interval.tv_nsec = init_nano;

	while ((thpool->ToBeAligned_q->len > 0) || (thpool->num_threads_working > 0)) {
		//printf("Aligner queue : len=%i, workers=%i\n", thpool->ToBeAligned_q->len, thpool->num_threads_working);
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
	while ((thpool->ToBeAligned_q->len > 0) || (thpool->num_threads_working > 0)) {
		sleep(max_secs);
		//printf("Aligner queue : len=%i, workers=%i\n", thpool->ToBeAligned_q->len, thpool->num_threads_working);
	}
}
//---------------------------------------------------------------

void terminateAlignerPool(AlignerPool_t * const thpool)
{
	/* End each thread 's infinite loop */
	thpool->threads_keepalive = 0;

	/* Give one second to kill idle threads */
	time_t start, end;
	double tpassed = 0.0;
	time (&start);
	while (tpassed < 1.0 && (thpool->num_threads_alive > 0)) {
		bsem_post_all(thpool->ToBeAligned_q->has_items);
		time (&end);
		tpassed = difftime(end,start);
	}

	/* Poll remaining threads */
	struct timespec polling_interval = { .tv_sec = 1, .tv_nsec = 0 };
	while (thpool->num_threads_alive > 0) {
		bsem_post_all(thpool->ToBeAligned_q->has_items);
		nanosleep(&polling_interval, NULL);
	}
}
//---------------------------------------------------------------

/* vim: tabstop=2 shiftwidth=2
 */
