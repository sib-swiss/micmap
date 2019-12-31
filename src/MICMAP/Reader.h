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
#ifndef _READER_H
#define _READER_H
#include <stdbool.h>
#include <linux/limits.h>
#include "jobqueue.h"

typedef struct ReaderTag {
	unsigned char Tag[kTagSize] __attribute__((aligned(16)));
	unsigned char Quality[kTagSize] __attribute__((aligned(16)));
	char Header[kHdrSize];
	unsigned long ordinal; 
	unsigned short AlignLen;
	unsigned short TagLen;
	unsigned short ValidQualityLen;
	unsigned short HeaderLen;
} ReaderTag_t;

// Folowwing job_t, this must start with a pointer prev to its class def.
typedef struct ReaderJob {
	struct ReaderJob * prev;
	ReaderTag_t * Tags;
	unsigned int TagCnt;
} ReaderJob_t;

typedef struct ReaderMemory {
	char * Ptr;
	size_t Size;
	ReaderJob_t * restrict Jobs;
	unsigned char IsHugeTLB;
} ReaderMemoryBulk_t;

typedef	struct FASTQ_struct {
	jobqueue_t ReaderJobQueue;
	jobqueue_t ReaderMemorySlot;
	unsigned long startOrdinal;
	const char * restrict tagfile1;
	const char * restrict tagfile2;
	ReaderMemoryBulk_t Memory;
	bool PairedEnd;
	char MinimalQualityScore;
} FASTQ;

void* readFASTQ(FASTQ * const fastq);

int allocateReaderMemory(FASTQ * const restrict ReaderInputs, const unsigned int nBlocks);
void freeReaderMemory(FASTQ * const restrict ReaderInputs);

#endif

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
