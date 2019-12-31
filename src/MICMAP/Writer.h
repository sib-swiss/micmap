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
#ifndef _WRITER_H
#define _WRITER_H
#include <pthread.h>
#include <stdatomic.h>
#include <stdbool.h>
#include "Genome.h"
#include "statistics.h"
#include "topology.h"

//--------------------------------------------------------------- 
// DEFINITIONS 
//--------------------------------------------------------------- 
#define kNOT_ACGT 0xFF

//--------------------------------------------------------------- 
// STRUCTURE DEFINITIONS 
//---------------------------------------------------------------
typedef struct WriterPool {
	pthread_t * threads;                   /* queue of pointer to threads               */
	jobqueue_t * ToBeProcessed_q;          /* pointer to the job queue to be processed  */
	const char * seedfn;                   /* seed file name                            */
	COUNTS * cnts;                         /* pointer to the statistics counts          */
	const Genome_t * Genome;               /* pointer to the genome                     */
#ifdef USE_AFFINITY
	const cpu_pool_t * Affinities;          /* affinities for flush buffer threads       */
#endif
	atomic_ulong WrittenAlignments;        /* keep an eye on the total number treated   */
	unsigned int FlushThreadPerWriter;     /* number of flush thread per writer thread  */
	bool PairedEnd;                        /* expecting paired ends                     */
	int num_threads;                       /* total number of threads in the pool       */
	atomic_int num_threads_alive;          /* threads currently alive                   */
	atomic_int num_threads_working;        /* threads currently working                 */
	volatile int threads_keepalive;        /* Used to stop threads                      */
} WriterPool_t;

//--------------------------------------------------------------- 
// FUNCTIONS 
//---------------------------------------------------------------
WriterPool_t* createWriterPool(const Genome_t * const restrict Genome,
                               jobqueue_t * const restrict AlignerProcessedQueue,
                               const char * const restrict seedfn, COUNTS * const restrict cnts,
                               const cpu_pool_t * const restrict Affinities,
                               const bool PairedEnd, const unsigned int N, const unsigned int nFlushBuffers);
void freeWriterPool(WriterPool_t * restrict Pool);
void terminateWriterPool(WriterPool_t * restrict Pool);

#endif /* _WRITER_H */

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
