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
#ifndef _GTL_DISTPATCH_H
#define _GTL_DISTPATCH_H
#include <stdbool.h>
#include <sys/types.h>
#include <sys/time.h>
#include <semaphore.h>
#include <pthread.h> 

#ifdef __cplusplus
#define restrict 
extern "C" {
#endif

#include "GTL.h"
#include "MATCHkonst.h"
#include "threadpool.h"

#define GTL_INDEX_SPLIT 32

/* Do not forget to alter the Suffix array in code file accordingly */
enum GTLHeaderExtremaType {
	None = 0,
	FirstTag = 1,
	Range = 2
};

typedef struct DirectLink {
	unsigned int ID;
	unsigned int Left;
	unsigned int Right;
	unsigned int padding;
	unsigned int Header;
	unsigned int Mismatch;
	unsigned int Cigar;
	unsigned int qs;
} DirectLink_t;

typedef struct GTLBlockInfo {
	void * volatile whatever; 	/* use to hold either the prev in the jobqueue, or the pointer to the decompressed block */
	TLBHEADER thd;
	off_t offset;
	int fd;
	volatile int usingIt;
	DirectLink_t PointToFraction[GTL_INDEX_SPLIT];
	pthread_mutex_t decompressing; /* Avoid race conditions when one is decompressing and the other needs it */
	unsigned char provenance;
} GTLBlockInfo_t;

typedef struct GTLBlockDispatchRange {
	struct GTLBlockDispatchRange * prev; /* Needs to be first as we dereference it to job_t */
	void * params;
	GTLBlockInfo_t * restrict Block;
	unsigned int GenomeStart;
	unsigned int GenomeEnd;
	unsigned int TotalJobs;
} GTLBlockDispatchRange_t;

typedef struct GTLBlockDispatchFile {
	GTLBlockInfo_t Block;
	void * params;
	const char * FileName;
	struct GTLBlockDispatchFile * const restrict * nPreviousBlocks;
	sem_t * semaphore;
	unsigned int TotalJobs;
} GTLBlockDispatchFile_t;

typedef struct GTLlist {
	unsigned int nFiles;
	int * restrict fds;
	const char * const restrict * FileNames;
	GTLBlockInfo_t * restrict Blocks;
	size_t nBlocks;
	enum GTLHeaderExtremaType PositionType;
} GTLList_t;

typedef struct GTLRawData {
 	unsigned long Ordinal;
	unsigned int Location;
	int AlignmentRange[2];            // That is genome offset and genome length, thus using cigar to compute
	int ProfileRange[2];
	unsigned char Tag[kMaxReadLen+1]; // carriage return
	unsigned char qs[kMaxReadLen];
	unsigned char Cigar[kMAXcigarSize];
	unsigned char Mismatch[kMAXmismatchSize]; 
 	unsigned char revNmer;
	unsigned char CigarLength;
 	unsigned char MismatchLength;
	unsigned char TagLength;
	unsigned char Provenance;
} GTLRawData_t;

typedef struct _decodedPair {
	unsigned long ordinal;
	unsigned int genomepos;
	int taglen1;
	int taglen2;
	int reverseTAG1;
	int reverseTAG2;
	int delta;
} decodedPair_t;

threadpool_t * createGTLThreadPool(void* (*Fct)(threadpool_t * const restrict),
                                   const cpu_set_t * const restrict affinities, const size_t jobvarsize, const int nThreads, const int onHold);
GTLList_t * createGTLList(const char * const restrict * const Files, const unsigned int nFiles);
void destroyGTLThreadPool(threadpool_t * const restrict thpool);
void destroyGTLList(GTLList_t * GL);
int indexGTLList(GTLList_t * const restrict list, const enum GTLHeaderExtremaType Order,
                 const _Bool DoNotUseExisting, const unsigned char provenance);
int writeGTLFileIndex(const GTLList_t * const restrict list, const char * const restrict Directory);
void sortGTLlist(GTLList_t * const restrict list);
int getRealBlockMinMaxPosition( GTLList_t * const restrict list, const int nThreads, 
																const unsigned char DoSort);
unsigned int dispatchGTLBlockRange(const GTLList_t * const restrict list,
                                   threadpool_t * const restrict thpool,
                                   const unsigned int chr, const unsigned int GenomeStart,
                                   unsigned int GenomeEnd,
                                   const unsigned char flags, const _Bool EndingJob,
													         void * const restrict params);
int dispatchGTLBlockFile(const GTLList_t * const restrict list, threadpool_t * const restrict thpool,
												 const unsigned char flags, const _Bool EndingJob, const unsigned int EndingJobEvery,
												 void * const restrict params);
void ExtractData(decodedPair_t * const restrict dpp, const unsigned char * restrict diff,
								 const unsigned char * restrict cigar, GTLRawData_t * const restrict * data);
void ExtractDataNasSoftClip(decodedPair_t * const restrict dpp, const unsigned char * restrict diff,
                            const unsigned char * restrict cigar, GTLRawData_t * const restrict * data);
#ifdef __cplusplus
};
#endif

#endif /* _GTL_DISTPATCH_H */

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
