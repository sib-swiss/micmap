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
#ifndef _ALIGNER_H
#define _ALIGNER_H
#include <pthread.h>
#include "jobqueue.h"
#include "topology.h" 
#include "Genome.h"
#include "Reader.h"
#include "Decoder.h"

typedef jobqueue_t memqueue_t;

typedef struct AlignerPool {
	memqueue_t threads_q;                  /* queue of pointer to threads               */
	pthread_mutex_t thcount_lock;          /* used for thread count etc                 */
	jobqueue_t * ToBeAligned_q;            /* pointer to the job queue to be processed  */
	jobqueue_t ToBeWritten_q;              /* job queue of aligned data                 */
	volatile int num_threads_alive;        /* threads currently alive                   */
	volatile int num_threads_working;      /* threads currently working                 */
	volatile int threads_keepalive;        /* Used to stop threads                      */
	volatile int threads_on_hold;
} AlignerPool_t;

typedef struct ThreadItem {
	struct ThreadItem * prev;
	pthread_t thread;
} ThreadItem_t;

typedef struct AlignerOutputBlock {
	struct AlignerOutputBlock * prev;
	void * data;
	memqueue_t * BelongsTo_q;
	size_t length;
	atomic_ulong Used;
	atomic_ulong Stored;
	unsigned int index;
} AlignerOutputBlock_t;

/* Structure of the output is as follow, that is 7 Bytes for the header information,
 * then we will have the block consisting of those informations
 */
#if kHdrSize >= 256
#error HeaderLenght is not sufficient here!
#endif
typedef struct outputPairBlockHeader {
	unsigned long Ordinal1;
	unsigned long Ordinal2;
	unsigned short TagLength1;
	unsigned short TagLength2;
	unsigned char HeaderLength1;
	unsigned char HeaderLength2;
	unsigned char nAlignmentResults;
	unsigned char StreamedNext; // bitfield to know who is streamed after AlignmentResults
	                            // if any. Bit 0 is for Tag1, Bit 1 is for Tag2, order must
	                            // kept.
	                            // 0x1111 0011
	                            //   ^^^^   ^^
	                            //   ||||   |+- Tag 1
	                            //   ||||   +-- Tag 2
	                            //   |||+------ SoftClip Tag 1
	                            //   ||+------- SoftClip Tag 2
	                            //   |+-------- SoftClip alternate Tag 1
	                            //   +--------- SoftClip alternate Tag 2
	                            
} outputPairBlockHeader_t ;

/* Each Alignment is then a stream of the following */
typedef struct Alignment_struct {
	int score;
	unsigned int DecisionTreeHistory;   // record decision tree
	unsigned int chr_pos;               // chromosome position
#ifndef NDEBUG
	unsigned int GenomeSeekStart;       // starting position to seek tag
	unsigned int GenomSeekStop;         // ending position to seek tag
#endif
	unsigned short chr;                 // chromosome number
	unsigned short MismatchCount;
	unsigned char revNmer;              // is the Nmer reported issued from a reverse
	                                    // complement of the genome
#ifdef EXPORT_STATES
	char States[3*kTagSize/2];          // State string for debugging
#endif
#ifdef DEBUG_TAG_MAPPING_PROCESS
	unsigned char AnchorPositionInTag;
#endif 
	unsigned char Mismatch[kMAXmismatchSize]; // enough to encode up to kMaxMismatch
	                                          // mismatches in a tag up to 255nt.
	unsigned char MismatchLength;
	unsigned char Cigar[kMAXcigarSize];
	unsigned char CigarLength;
#ifndef NDEBUG
	unsigned char AnchorID;
#endif
} AlignmentResults_t;

typedef struct AlignerBlockArgs {
	AlignerPool_t * Pool;
	memqueue_t * ReaderMemoryQueue;    // Required in sender thread to place back memory slot in the queue 
	memqueue_t * DecoderMemoryQueue;   // Required in sender thread to place back memory slot in the queue
	jobqueue_t * DecoderInternalJob_q; // Required in sender thread to place back internal decoder job
	Genome_t * Genome;
	void* ExtraPtr;                    // cpud void* (*ThreadFct)(void*), MIC: scif_epd_t*
#ifdef USE_AFFINITY
	Affinity_Mask_t * Affinities;      // Affinity array for CPU thread worker threads
#endif
	off_t OutputBlockSize;
	int ID;
	int SingleEnd;
	unsigned int GenomeChunkSize;
	unsigned int nAligners;
	unsigned int nInputBlock;
	unsigned int nOutputBlock;
} AlignerBlockArgs_t;

typedef void* (*AlignerFct)(AlignerBlockArgs_t * const restrict data);

AlignerPool_t* createAlignerPool(jobqueue_t * const restrict ToBeAligned_q);
void waitForAlignerPool(AlignerPool_t * const restrict thpool);
void terminateAlignerPool(AlignerPool_t * const thpool);
void* MicAligner(AlignerBlockArgs_t * const restrict data);
void* cpuBlockAligners(AlignerBlockArgs_t * const restrict data);

#endif

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
