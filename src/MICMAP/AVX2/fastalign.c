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
#include "config.h"
#include "constants.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <immintrin.h>
#include <assert.h>
#include <stdio.h> 
#include <sys/mman.h>
#define USE_SHORT_INT
#include "ALIGNkonst.h"
#include "GTLkonst.h"
#include "MATCHkonst.h"
#include "../FastAlign/fastalign.h"

// #define STORE_FUNCTION(where,what) _mm512_store_si512(where, what)
#define PREFETCHL1(address) _mm_prefetch((char*)(address), _MM_HINT_T0);
#define EVICTL1(address)

#define __DECLARE_STATIC_VARIABLES__
#include "pftools_systolic.h"
#undef __DECLARE_STATIC_VARIABLES__


#define __EXTRA_FUNCTIONS__
#include "pftools_systolic.h"
#undef __EXTRA_FUNCTIONS__

static int freeStorage(GlobalData_t * const restrict data)
{
	munmap(data->MatrixStorage, HUGE_PAGE_SIZE);
	return 0;
}
//---------------------------------------------------------------

static int allocStorage(GlobalData_t * const restrict data)
{
		return 0;
}
//---------------------------------------------------------------

static int CreateMatrix(GlobalData_t * const restrict data)
{
	#define __DECLARE_VARIABLES__
	#include "pftools_systolic.h"
	#undef __DECLARE_VARIABLES__
	
	const int ReadLength = data->TagLength;
	const int GenomeLength = data->GenomeLength;
	const int MatrixLD = data->MatrixLD;
	unsigned char * restrict ptrState = &(data->StateSequence[STATE_MEMORY_SIZE-1]);
	const unsigned char * const restrict ptrTag = data->Tag;
	const unsigned char * const restrict genome = data->Genome;
	short int * const restrict IDM_Matrix = (short int*) data->MatrixStorage;
	
	#define __FILL_MATRIX__
	#include "pftools_systolic.h"
	#undef __FILL_MATRIX__
	
	const size_t num = (ReadLength-1) & (size_t) 0xF;
	data->EndOnMatchOrDeletion = MorD[num];
	data->BestAlignmentIndex = Id[num];
	return 0;
}
//---------------------------------------------------------------

static int GetBestAlignment(GlobalData_t * const restrict data)
{
	int ret = 0;
	size_t counter = 0;
	unsigned int MismatchCount = 0;
	unsigned int InsertionCount = 0;
	unsigned StateChange = 0;
	
	const unsigned char * const restrict ptrTag = data->Tag;
	unsigned char * restrict ptrCigar = data->Cigar + kMAXcigarSize;
	unsigned char * restrict ptrMismatch = data->Mismatch + kMAXmismatchSize;
	const unsigned char * const restrict cigarStorageBottomLimit = data->Cigar;
	const unsigned char * const restrict mismatchStorageBottomLimit = data->Mismatch;
	const unsigned char * const restrict genome = data->Genome;
	const int ReadLength = data->TagLength;
	int Align[2];
	const unsigned int EndOnMatchOrDeletion = data->EndOnMatchOrDeletion;
	int iseq = data->BestAlignmentIndex;
	short int * const restrict IDM_Matrix = data->MatrixStorage;
	const int MatrixLD = data->MatrixLD;
	
	#define __GET_ALIGNMENT__
	#include "pftools_systolic.h"
	#undef __GET_ALIGNMENT__
	if (ret) return ret;
	return 0;
}
//---------------------------------------------------------------

static int GetStateSequence(GlobalData_t * const restrict data)
{
	int ret;
	size_t counter = 0;
	const int ReadLength = data->TagLength;
	int iseq = data->BestAlignmentIndex;
	int iprf = ReadLength-1;
	const unsigned char * const restrict ptrTag = data->Tag;
	const unsigned char * const restrict genome = data->Genome;
	const unsigned int EndOnMatchOrDeletion = data->EndOnMatchOrDeletion;
	const int MatrixLD = data->MatrixLD;
	short int * const restrict IDM_Matrix = (short int*) data->MatrixStorage;
	unsigned char * restrict ptrState = &data->StateSequence[STATE_MEMORY_SIZE-1];
	const unsigned char * const restrict StateLimit = data->StateSequence;
	int Align[2];
	
	int SoftClipBoundary = data->SoftClipBoundary;
	/* Limit SoftClip Size according to kMinTagLenNotSoftClipped */
	assert(ReadLength >= kMinTagLenNotSoftClipped);
	if (!(data->revNmer)) {
		if(SoftClipBoundary<kMinTagLenNotSoftClipped) SoftClipBoundary = kMinTagLenNotSoftClipped;
	}
	else {
		register const int itmp = (ReadLength-kMinTagLenNotSoftClipped);
		if (SoftClipBoundary > itmp) SoftClipBoundary = itmp;
	}
						
	unsigned int SoftClipMismatchCount = 0U;
	unsigned int MismatchCount  = 0U;
	unsigned char ContainsInDel = 0;
	unsigned int RequiredMMSize = 0U;
	unsigned int SoftClipEvictedMMSize = 0U;
	unsigned int RequiredCigSize = 1U;
	unsigned int SoftClipEvictedCigSize = 0U;
		
	*ptrState = '\0';
	
	#define __GET_ALIGNMENT_STATE__
	#include "pftools_systolic.h"
	#undef __GET_ALIGNMENT_STATE__
	
	data->States = ptrState;
	data->SoftClipMismatchCount = SoftClipMismatchCount;
	data->MismatchCount  = MismatchCount;
	data->RequiredMMSize = RequiredMMSize;
	data->SoftClipEvictedMMSize = SoftClipEvictedMMSize;
	data->RequiredCigSize = RequiredCigSize;
	data->SoftClipEvictedCigSize = SoftClipEvictedCigSize;
	data->AlignmentRange[0] = Align[0];
	data->AlignmentRange[1] = Align[1];
	const enum StateTypeOffset Offset = (EndOnMatchOrDeletion) ? DeletionOffset : MatchOffset;
	unsigned int uitmp = *(Location(IDM_Matrix, MatrixLD, ReadLength, data->BestAlignmentIndex, Offset));
	data->score = ToScore((int) uitmp);
	
	if (ret) return ret;
	return 0;
}
//---------------------------------------------------------------

extern int EncodeInternalAlignment(GlobalData_t * const restrict data, const unsigned int SoftClipMismatchToRemove );

static int ComputeCigarAndMismatch(GlobalData_t * const restrict data)
{
	unsigned int SoftClipMismatchToRemove = 0U;
	// POTENTIAL SOFTCLIP RESCUE
	if (data->MismatchCount > kMaxMismatch) {
		SoftClipMismatchToRemove = data->MismatchCount - kMaxMismatch;
		if (!(data->revNmer)) {
			// check we will be able to encode, keep min mind that SoftClipEvictedCigSize in this orientation
			// cannot be totally removed, only  SoftClipEvictedCigSize-1. In addition, we need to keep space 
			// for the Sofclip. Consequently, the maximum allowed size is kMaxEncodedCigar-1 with 
			//     kMaxEncodedCigar-1 >= RequiredCigSize - (SoftClipEvictedCigSize-1)
			if (   (data->RequiredMMSize-data->SoftClipEvictedMMSize) <= kMaxEncodedMismatch
				  && (data->RequiredCigSize-data->SoftClipEvictedCigSize) < (kMaxEncodedCigar-1) )
				data->MismatchCount -= data->SoftClipMismatchCount;
			else
				return -1;
		}
		else {
			// check we will be able to encode
			if (data->SoftClipEvictedMMSize <= kMaxEncodedMismatch && data->SoftClipEvictedCigSize <= kMaxEncodedCigar)
				data->MismatchCount = data->SoftClipMismatchCount;
			else 
				return -1;
		}
		if (data->MismatchCount < kMaxMismatch) data->MismatchCount = kMaxMismatch;
	}
	else if (data->RequiredMMSize > kMaxEncodedMismatch || data->RequiredCigSize > kMaxEncodedCigar) {
		return -1;
	}

	if (data->MismatchCount <= kMaxMismatch) {
		return EncodeInternalAlignment(data, SoftClipMismatchToRemove);
	}
	
	return -1;
}
//---------------------------------------------------------------

ComputeFunctions_t avx2_systolic = {
	.createMatrix = CreateMatrix,
	.getBestAlignment = GetBestAlignment,
	.getStateSequence = GetStateSequence,
	.computeCigarAndMismatch = ComputeCigarAndMismatch,
	.allocStorage = allocStorage,
	.freeStorage = freeStorage
};

/* vim: tabstop=2 shiftwidth=2
 */
