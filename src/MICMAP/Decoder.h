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
#ifndef _DECODER_H
#define _DECODER_H
#include <stdatomic.h>
#include <sched.h>
#include <stdio.h>
#include <sys/types.h>
#include "jobqueue.h"
#include "topology.h"
#include "Reader.h"

//--------------------------------------------------------------- 
// DEFINITIONS 
//---------------------------------------------------------------
#define WHICH16MER_chr(w)      ((w) >> 28)
#define WHICH16MER_pos(w)      ((w) & 0x0FFFFFFF)
#define WHICH16MER_chr_part(w) ((w) & 0xF0000000)

#define INFO_MATCHPOS(I)       (((I) >> 8) & 0xFF)

//--------------------------------------------------------------- 
// STRUCTURE DEFINITIONS 
//--------------------------------------------------------------- 
typedef struct DecoderElement
{
	// identifies which unique exact match 16mers was matched [ can be used to retrieve the genomic pos with host ]
	unsigned int which16mer;
	/* Structure of which16mer is as follow
	 * 31 28                            0 bits
	 *  xxxx|xxxxxxxxxxxxxxxxxxxxxxxxxxxx
	 *     |                            |
	 *     |                            +--  position on chromosome (28 bits)
	 *     +-------------------------------  chromosome number      (4  bits)
	 */
	// identifies which of the 2nt to use, the position of the match, if tag is reverse, and number of matches 
	unsigned int info;
	/* Structure of info is as follow: 
	 * 31                       8        0 bits
	 *  xxxxxxxxxxxxxxxx|xxxxxxxx|xxxxxxxx
	 *                 |        |      |||
	 *                 |        |      ||+-- revNmer (1 bit)
	 *                 |        |      |+--- rev (1 bit) - tested tag was reversed
	 *                 |        |      +---- pos (6-bit - which of the kMAX16mersToTest was successful)
	 *                 |        +----------- match position (8 bits) (tagpos)
	 *                 +-------------------- which4nt or whichNnt (2 or 4 bits)
	 */
} DecoderElement_t;

typedef struct DecoderMemoryElement {
	struct DecoderMemoryElement * prev;
	DecoderElement_t (* restrict Anchors)[kMaxTagsPerBuffer][kDesiredHitCnt];
	atomic_uint * restrict count;
} DecoderMemoryElement_t;

typedef struct DecoderMemoryBulk {
	DecoderElement_t (* restrict Anchors)[kMaxTagsPerBuffer][kDesiredHitCnt];
	atomic_uint * restrict counters;
	jobqueue_t * restrict DecoderMemory_q;
	DecoderMemoryElement_t * restrict DecoderMemoryElements;
	jobqueue_t ToBeAligned_q;
} DecoderMemoryBulk_t;

// Following job_t, this must start with a pointer prev to its class def.
typedef struct DecoderJobElement {
	struct DecoderJobElement * prev;
	union DecoderData {
		DecoderElement_t (* restrict Anchors)[kDesiredHitCnt]; /* pointer to where to store decoding */
		DecoderMemoryElement_t * DecoderJob;                   /* if last, points to Decoder memory element */
	} Decoder;
	union ReaderData {
		const ReaderTag_t * Tags;    /* pointer to input tags */
		ReaderJob_t * ReaderJob;     /* if last job within the reader job, this is its address */
	} Reader;
	unsigned int nToBeProcessed;   /* Number of tags to decode in this job, 0 means last job */
	atomic_uint * count;           /* Overall count of all the decoding jobs working on this Reader Job */
} DecoderJobElement_t;

typedef struct ValidData_struct {
	unsigned char * restrict valid_tbl; // contain table of which entries of the Nmers table are valid (contain actual data)
	unsigned int * restrict following_nt_tbl;	// contain table of unique exact match Nmers.
	off_t valid_tbl_size;
	off_t following_nt_tbl_size;
} VALIDDATA;

typedef struct DecoderPool {
	pthread_t * threads;                   /* pointer to threads                           */
	job_t * jobs;                          /* pointer to all jobs memory, only for freeing */
	jobqueue_t * ToBeAligned_q;	           /* Pool of aligner masters                      */
	VALIDDATA * ValidDataTable;            /* Valid data table                             */
	const unsigned int * chrpos_tbl[kFollowNtCount];
	const unsigned char * strand_tbl;
	const unsigned int * following_nt_tbl;
	
	jobqueue_t PoolJob_q;                  /* pointer to the job queue                     */
	jobqueue_t PoolMemory_q;               /* empty memory slot queue                      */
	int num_threads;                       /* total number of threads in the pool          */
	int SharedMemoryFileDecriptor;
	atomic_uint num_threads_alive;         /* threads currently alive                      */
	atomic_uint num_threads_working;       /* threads currently working                    */
	volatile int threads_keepalive;        /* Used to stop threads                         */
} DecoderPool_t;

typedef struct ReaderToDecoderArgs {
	FASTQ * restrict ReaderData;
	DecoderPool_t * restrict DecoderPool;
	jobqueue_t * restrict DecoderMemory_q;
} ReaderToDecoderArgs_t;

typedef void* (*DecodeTagsPtr)(DecoderPool_t * const restrict);

//--------------------------------------------------------------- 
// FUNCTIONS
//--------------------------------------------------------------- 
void freeValidTable(VALIDDATA * restrict ValidDataTable);
VALIDDATA * loadValidTable(const char * const restrict prefix, const char NT);

int loadDecodingTable(DecoderPool_t * const restrict thpool, const char * const restrict prefix,
                      const char NT, const _Bool useHugeTLB);
int DecodingTableToRAM(const char * const restrict path, const char * const restrict prefix,
                       const char NT);
void freeDecodingTable(DecoderPool_t * const restrict thpool);

int createDecoderPool(DecoderPool_t * const restrict thpool, const DecodeTagsPtr Fct,
                      const Affinity_Mask_t * const restrict affinities,
                      DecoderMemoryBulk_t * const restrict InOut, const unsigned int nThreads);
void destroyThreadPool(DecoderPool_t * restrict thpool);
DecoderMemoryBulk_t * allocateDecoderMemory(const size_t nBlocks);
void freeDecoderMemory(DecoderMemoryBulk_t * data);

void* ReaderToDecoderThread(ReaderToDecoderArgs_t * const restrict data);

void* DecodeTags(DecoderPool_t * const restrict data);
void* DecodeTags_avx2(DecoderPool_t * const restrict data);

#endif /*_DECODER_H*/
/* vim: tabstop=2 shiftwidth=2
 */
