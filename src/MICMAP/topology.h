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
#ifndef _TOPOLOGY_H
#define _TOPOLOGY_H
#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <sched.h>
#include "system.h"

//---------------------------------------------------------------
// DEFINITIONS 
//---------------------------------------------------------------
// #define kDECODETHREADCNT 1
// #define kWRITE_TO_DISK_THREAD_CNT 3 /* must be power of 2 ?  check relation to gFilledBufCnt */
// #define kCPU_ASSIGN_WRITING 10
// #define kMaxALIGNTHREADcnt 32


//---------------------------------------------------------------
// STRUCTURE DEFINITIONS 
//---------------------------------------------------------------
typedef uint64_t MicCoreMask;
typedef uint8_t MicThreadMask;

struct MicTopologyId {
	MicCoreMask Cores;
	MicThreadMask Threads;
};

typedef struct cpu_pool_struct {
	Affinity_Mask_t Master;
	Affinity_Mask_t * Slaves;
} cpu_pool_t;

typedef struct Topology_struct {
#ifdef USE_AFFINITY
	Affinity_Mask_t * Reader;
	Affinity_Mask_t * MicAligners;
	cpu_pool_t * Aligners; 
	cpu_pool_t * Writers;
	cpu_pool_t * Decoders;
#endif

	unsigned int nDecoders;
	unsigned int nAligners;
	unsigned int nAlignerThreads;
	unsigned int nWriters;
	unsigned int nWriterThreads;
		
	unsigned int nReaderBlocks;
	unsigned int nDecoderBlocks;
	unsigned int nAlignerInBlocks;
	unsigned int nAlignerOutBlocks;
	off_t nAlignerOutBlockSize;
} Topology_t;

//---------------------------------------------------------------
// Globals 
//---------------------------------------------------------------
extern Topology_t gTopology;

//---------------------------------------------------------------
// FUNCTIONS 
//---------------------------------------------------------------
void printJSON(const SystemInfo * const restrict info, FILE * const restrict fd);
int loadTopology(const SystemInfo * const restrict info, const char * const restrict FileName);
void printParameters(FILE* const out);

#endif /* _TOPOLOGY_H */

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
