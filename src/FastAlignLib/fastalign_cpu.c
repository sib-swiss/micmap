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
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <smmintrin.h>
#include <assert.h>
#include "fastalign.h"

////////////////////////////// FOR SSE to WORK //////////////////////////////////////////////
// WARNING: MATCH and INSERTION should be next to each other (in this ordering)
//          and either in the upper part of xmm or the lower part.
//          Do no mix them and correct the storing function according to your above choice.
static inline void __attribute__((always_inline))
StoreMatchInsertion( __m64 * const _address, const __m128 _reg)
{
	//_mm_storel_pi(_address, _reg);
	_mm_storeh_pi(_address, _reg);
}
//---------------------------------------------------------------

static inline void __attribute__((always_inline))
operation(__m128i x, __m128i* inout, const size_t row, const size_t column, const size_t lda)
{
	_mm_store_si128(&inout[row*lda+column], x);
}
//---------------------------------------------------------------

static int freeStorage(GlobalData_t * const restrict data)
{
	if (data->MatrixStorage) { _mm_free(data->MatrixStorage); data->MatrixStorage = NULL; }
  if (data->WorkStorage) { _mm_free(data->WorkStorage); data->WorkStorage = NULL;}
  return 0;
}
//---------------------------------------------------------------

static int allocStorage(GlobalData_t * const restrict data)
{
	data->MatrixLD = (data->TagLength+1);
	data->WorkStorage = (void*) _mm_malloc(63+(2*data->MatrixLD)*sizeof(sIOP),64);
	if (data->WorkStorage) {
		data->MatrixLD = (data->MatrixLD + 3) & ~(3);
		data->MatrixStorage = (void*) _mm_malloc(data->MatrixLD*(1+data->GenomeLength)*sizeof(lScores), 64);
		if (data->MatrixStorage) return 0;
	}
	freeStorage(data);
	return -1;
}
//---------------------------------------------------------------

static int CreateMatrix(GlobalData_t * const restrict data)
/*
 * WARNING: for SSE version, WORK should be aligned on cache size (64b) and 4 times the (profile size + 1)*sizeof(int)
 *          + 63 to align to cache line
 *          Matrix should be of size at least (profile size + 1)*(sequence size) and aligned on 16b.
 */
{
	__m128i KOPD;
	
	const __m128i __ClearMask = _mm_set1_epi32(CLEAR_MASK);
	const __m128i __N         = _mm_set1_epi8(N);
	
	register const __m128i __Mx          = prf.Match.Vector;
	register const __m128i __Ix          = prf.Insertion.Vector;
	register const __m128i __Dx          = prf.Deletion.Vector;
	register const __m128i __FirstColumn = _mm_load_si128((__m128i*)RowInit);
	register const __m128i __KOPDInit    = _mm_set1_epi32 (RowInit[DELETION]);
	
	register __m128i __MatchScore     = _mm_set1_epi32 (prf.M);
	register __m128i __MismatchScore  = _mm_set1_epi32 (prf.m);
	register __m128i __InsertionScore = _mm_set1_epi32 (prf.I);
	register __m128i __DeletionScore  = _mm_set1_epi32 (prf.D);
	
	__MatchScore     = _mm_add_epi32(__MatchScore,     __Mx);
	__MismatchScore  = _mm_add_epi32(__MismatchScore,  __Mx);
	__DeletionScore  = _mm_add_epi32(__DeletionScore,  __Dx);
	__InsertionScore = _mm_add_epi32(__InsertionScore, __Ix);
	
	const union sIOP * restrict IOP_R;
	union sIOP * restrict IOP_W = (union sIOP *) data->WorkStorage;
	
	const __m128i __MatchMask         = _mm_set1_epi32(PRIORITY_MATCH     << STATE_SHIFT);
	const __m128i __InsertionMask     = _mm_set1_epi32(PRIORITY_INSERTION << STATE_SHIFT);
	const __m128i __DeletionMask      = _mm_set1_epi32(PRIORITY_DELETION  << STATE_SHIFT);
	const __m128i __ExtraMask         = _mm_set1_epi32(PRIORITY_EXTRA     << STATE_SHIFT);
	const __m128i __UpperNegativeMask = _mm_set1_epi32(LEFT_NEGATIVE_SCORE_MASK); /* three major bits on */
	const __m128i __RealMatchMask     = _mm_set1_epi32(MISMATCH_MASK);
	
	const size_t ReadLength = data->TagLength;
	const size_t MatrixLD = data->MatrixLD;
	unsigned char * const restrict genome = data->Genome;
	const unsigned char * restrict ptrTag = data->Tag; 
	/*
	 * Initialize Insertion and Match Entrance Line using FirstSequenceProtein
	 */  
	{
		__m128i * const restrict MatrixPtr = (__m128i*) data->MatrixStorage;
		/*
		 * PROFILE COLUMN 0 entrance
		 */
		IOP_W[0].Element.Match     = RowInit[MATCH];
		IOP_W[0].Element.Insertion = RowInit[INSERTION];
		KOPD                       = _mm_set1_epi32(RowInit[DELETION]);
		
		// Store all scores to Matrix
		operation(__FirstColumn, MatrixPtr, 0, 0, MatrixLD);
		
		
		register union sIOP * restrict pIOP = &IOP_W[1];
		
		/*
		 * LOOP THROUGH THE REST OF THE PROFILE
		 */
		for (size_t iprf=1; iprf<=ReadLength; ++iprf ) {
			KOPD = _mm_and_si128(KOPD, __ClearMask);
			
			// Add KD to Transitions
			__m128i __TransitionsD = _mm_add_epi32(__DeletionScore, KOPD);
			
			KOPD = _mm_shuffle_epi32(__TransitionsD, (DELETION << 6) | (DELETION << 4) | (DELETION << 2) | DELETION);
			
			// Store IOPI and IOPM
			StoreMatchInsertion( &(pIOP->mm), (__m128) __TransitionsD);
			pIOP++;
			
			// Store all scores to Matrix
			operation(__TransitionsD, MatrixPtr, 0, iprf, MatrixLD);     
		}
	}
		
	// Swap and assign Read and write pointers
	IOP_R = IOP_W;
	IOP_W = (union sIOP*) (((uintptr_t) &((union sIOP*)data->WorkStorage)[MatrixLD] + 63) & ~63);
	
	__m128i * const restrict MatrixPtr = ((__m128i*) data->MatrixStorage) + MatrixLD;
	/*
	 * LOOP THROUGH THE SEQUENCE STRING
	 */
	
	for ( size_t iseq=0; iseq < data->GenomeLength; ++iseq) {
		const __m128i __SequenceDNA = _mm_set1_epi8(genome[iseq]);
		__m128i KOPM = _mm_set1_epi32(IOP_R[0].Element.Match);
		/*
		 * PROFILE COLUMN 0 entrance
		 */
		{
			KOPD = __KOPDInit;
			
			// Store IOPI and IOPM
			StoreMatchInsertion( &(IOP_W[0].mm), (__m128) __FirstColumn);
			
			// Store all scores to Matrix
			operation(__FirstColumn, MatrixPtr, iseq, 0, MatrixLD);
		}
				
		/*
		 * LOOP THROUGH THE REST OF THE PROFILE
		 */
		for (size_t iprf=1; iprf<=ReadLength; ++iprf ) {
			const __m128i __ReadSequence = _mm_set1_epi8(ptrTag[iprf-1]);
			const __m128i ReadIsN        = _mm_cmpeq_epi8(__ReadSequence, __N);
			__m128i DoWeMatch            = _mm_cmpeq_epi8(__SequenceDNA, __ReadSequence);
			DoWeMatch                    = _mm_or_si128(DoWeMatch, ReadIsN);
			const __m128i __AddtoM       = _mm_blendv_epi8(__MismatchScore, __MatchScore, DoWeMatch);
			
			__m128i __KI = _mm_set1_epi32(IOP_R[iprf].Element.Insertion);
			
			KOPD = _mm_and_si128(KOPD, __ClearMask);
			KOPM = _mm_and_si128(KOPM, __ClearMask);
			__KI = _mm_and_si128(__KI, __ClearMask);
			
			__m128i __KM = _mm_add_epi32(KOPM, __AddtoM);
			__KI = _mm_add_epi32(__KI, __InsertionScore);
			
			KOPM = _mm_set1_epi32(IOP_R[iprf].Element.Match);
			
			__m128i __max = _mm_max_epi32(__KM, __KI);
			
			__m128i __KD = _mm_add_epi32(KOPD, __DeletionScore);
			
			__max = _mm_max_epi32(__max, __KD);
			
			
			StoreMatchInsertion( &IOP_W[iprf].mm, (__m128) __max);
			operation(__max, MatrixPtr, iseq, iprf, MatrixLD);
			
			KOPD = _mm_shuffle_epi32(__max, (DELETION << 6) | (DELETION << 4) | (DELETION << 2) | DELETION);
		}

		// Swap Read and Write pointers
		const register union sIOP * const ptr = IOP_W;
		IOP_W = (union sIOP *) IOP_R;
		IOP_R = ptr;
	}
	return 0;
}
//---------------------------------------------------------------

static int GetBestAlignment(GlobalData_t * const restrict data)
{
	const int (* restrict Scores)[4]  = (const int (* restrict )[4]) data->MatrixStorage;
	int index = 0;
	int extrema = Scores[data->TagLength][EXTRA];
	const size_t ld = data->MatrixLD;
	const size_t prfLength = data->TagLength;
	const size_t SeqLength = data->GenomeLength;
	for (size_t i=1; i<=SeqLength; ++i) {
		if (Scores[i*ld+prfLength][EXTRA] > extrema)
		{
			extrema = Scores[i*ld+prfLength][EXTRA];
			index = i;
		}
	}
	
	data->score = ToScore(Scores[index*ld+prfLength][EXTRA]);
  data->AlignmentRange[1] = index;
	size_t counter     = 0;
	int iprf           = prfLength;
	size_t State       = EXTRA;
	while ( iprf >= 0 && index >= 0)
	{
		const unsigned int Move = (Scores[index*ld+iprf][State] & STATE_MASK) >> STATE_SHIFT;
		switch (Move)
		{
			case PRIORITY_MATCH:
				--index;
				--iprf;
				State = MATCH;
				break;
			case PRIORITY_INSERTION:
				--index;
				State = INSERTION;
				break;
			case PRIORITY_DELETION:
				--iprf;
				State = DELETION;
				break;
			case PRIORITY_EXTRA:
				State = EXTRA;
				goto OUT;
				break;
			default:
				fprintf(stderr,"\nUnknown move (%u) encountered from cell (%i,%i)\n", Move, index, iprf);
				--iprf;
				return -1;
				break;
		}
		++counter;
	}
OUT:
	data->AlignmentRange[0] = index+1;
	return counter;
}
//---------------------------------------------------------------

static int GetAlignmentSequence(GlobalData_t * const restrict data)
{
	const int (*Scores)[4]  = (const int (*)[4]) data->MatrixStorage;
	const size_t prfLength = data->TagLength;
	const size_t SeqLength = data->GenomeLength;
	
	int index          = data->AlignmentRange[1];
	size_t counter     = 0;
	int iprf           = prfLength;
	size_t State       = EXTRA;
	const size_t ld = data->MatrixLD;
	char * const restrict AlignmentSequence = data->AlignmentSequence;    
	const unsigned char * const restrict Sequence = data->Genome; 
	while ( iprf >= 0 && index >= 0)
	{
		if (counter >= data->AlignmentSequenceLength) {
			AlignmentSequence[-1] = '\0';
					return -1;
		}
		const unsigned int Move = Scores[index*ld+iprf][State] & 0x3;
		const unsigned char C = Sequence[index];
		switch (Move)
		{
			case PRIORITY_MATCH:
				--index;
				--iprf;
				State = MATCH;
				AlignmentSequence[counter] = C;
				break;
			case PRIORITY_INSERTION:
				--index;
				State = INSERTION;
				AlignmentSequence[counter] = (char) (  C + ( 'a' - 'A'));
				break;
			case PRIORITY_DELETION:
				--iprf;
				State = DELETION;
				AlignmentSequence[counter] = '-';
					    break;
			case PRIORITY_EXTRA:
				State = EXTRA;
				goto OUT;
				break;
			default:
				fprintf(stderr,"\nUnknown move (%u) encountered from cell (%i,%i)\n", Move, index, iprf);
				--iprf;
				goto OUT;
				break;
		}
		++counter;
	}
OUT:;
	/* Reverse the string */
	char * BackPtr = &AlignmentSequence[counter-1];
	 
	for (size_t i=0; i<counter/2; ++i) {
		const char c = AlignmentSequence[i];
		AlignmentSequence[i] = *BackPtr;
		*BackPtr-- = c;
	}
	if (counter & 0x1) {
		const char c = AlignmentSequence[counter/2];
		AlignmentSequence[counter/2] = *BackPtr;
		*BackPtr = c;
	 }
	AlignmentSequence[counter] = '\0';
	return 0;
}
//---------------------------------------------------------------

static int GetStateSequence(GlobalData_t * const restrict data)
{
	const int (* restrict Scores)[4]  = (const int (* restrict )[4]) data->MatrixStorage;
	const size_t ld = data->MatrixLD;
	int index = 0;
	{
		const size_t prfLength = data->TagLength;
		int extrema = NLOW;
		const size_t SeqLength = data->GenomeLength;
		for (size_t i=1; i<=SeqLength; ++i) {
			const int val = ToScore(Scores[i*ld+prfLength][EXTRA]);
			if (val > extrema)
			{
				extrema = val;
				index = i;
			}
		}
		
		data->score = extrema;
	  data->AlignmentRange[1] = index;
	}
	unsigned int State = EXTRA;
	
	const register int SoftClipBoundary = data->SoftClipBoundary;
	unsigned int SoftClipMismatchCount = 0U;
	unsigned int MismatchCount  = 0U;
	// 						unsigned int InsertionCount = 0U;
	unsigned int RequiredMMSize = 0U;
	unsigned int SoftClipEvictedMMSize = 0U;
	unsigned int RequiredCigSize = 0U;
	unsigned int SoftClipEvictedCigSize = 0U;
	
	unsigned char * restrict ptrState = data->StateSequence + STATE_MEMORY_SIZE - 1;
	const unsigned char * const restrict StateLimit = data->StateSequence; 
	const unsigned char * const restrict genome = data->Genome; 
	const unsigned char * restrict ptrTag = data->Tag;
	int iprf = data->TagLength;
	unsigned int previous_State = ((unsigned int) Scores[index*ld+iprf][EXTRA] & STATE_MASK) >> STATE_SHIFT;
	*ptrState = '\0';
	
	int ret = -1;

	while ( iprf >= 0 && index >= 0)
	{
		if (--ptrState < StateLimit) goto OUT;
		const unsigned int lScore = (unsigned int) Scores[index*ld+iprf][State];
		const unsigned int current_State = (lScore & STATE_MASK) >> STATE_SHIFT;
		if (current_State != previous_State) {
			RequiredCigSize++;
			previous_State = current_State;
		}
		
		switch(current_State) {
			case PRIORITY_MATCH:
				if (lScore & MISMATCH_MASK) {
					if (genome[index] == 'N' && ptrTag[iprf] != 'N') {
						*ptrState = 'n';
						RequiredMMSize++;
					}
					else if (ptrTag[iprf] == 'N') {
						*ptrState = 'N';
						RequiredMMSize++;
					}
					else {
						*ptrState = 'M';
					}
				}
				else {
					*ptrState = 'm';
					MismatchCount++;
					RequiredMMSize++;
				}
				State = MATCH;
				--index;
				--iprf;
				break;
			case PRIORITY_INSERTION:
				--index;
				*ptrState = 'D';
				State = INSERTION;
				break;
			case PRIORITY_DELETION:
				*ptrState = 'I';
				--iprf;
				RequiredMMSize++;
				State = DELETION;
				break;
			case PRIORITY_EXTRA:
				goto DONE;
				break;
			default:
				fprintf(stderr,"\nUnknown state (%u) encountered from cell (%i,%i)\n", current_State, index, iprf);
				--iprf;
				goto OUT;
				break;
		}
		
		if (iprf == SoftClipBoundary) {
			SoftClipMismatchCount  = MismatchCount;
			SoftClipEvictedMMSize  = RequiredMMSize;
			SoftClipEvictedCigSize = RequiredCigSize;
		}
	}
	DONE:;
	ret = 0;
	
	data->MismatchCount = MismatchCount;
	data->SoftClipMismatchCount = SoftClipMismatchCount;
	data->RequiredMMSize = RequiredMMSize;
	data->RequiredCigSize = RequiredCigSize;
	data->SoftClipEvictedMMSize = SoftClipEvictedMMSize;
	data->SoftClipEvictedCigSize = SoftClipEvictedCigSize;
	data->States = ptrState + 1; // avoid the last move
	data->AlignmentRange[0] = index;
	data->ProfileRange[0] = iprf;
	data->ProfileRange[1] = data->TagLength;
	
	OUT:
	return ret;
}
//---------------------------------------------------------------

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

ComputeFunctions_t cpu_std = {
	.createMatrix = CreateMatrix,
	.getBestAlignment = GetBestAlignment,
	.getStateSequence = GetStateSequence,
	.computeCigarAndMismatch = ComputeCigarAndMismatch,
	.allocStorage = allocStorage,
	.freeStorage = freeStorage
};

static int CreateMatrixBorder(GlobalData_t * const restrict data)
/*
 * WARNING: for SSE version, WORK should be aligned on cache size (64b) and 4 times the (profile size + 1)*sizeof(int)
 *          + 63 to align to cache line
 *          Matrix should be of size at least (profile size + 1)*(sequence size) and aligned on 16b.
 */
{
	__m128i KOPD;
	
	const __m128i __ClearMask = _mm_set1_epi32(CLEAR_MASK);
	const __m128i __N         = _mm_set1_epi8(N);
	
	register const __m128i __FirstColumn = _mm_load_si128((__m128i*)RowInit);
	register const __m128i __KOPDInit    = _mm_set1_epi32 (RowInit[DELETION]);
	
	register __m128i __MatchScore     = _mm_set1_epi32 (prf.M);
	register __m128i __InsertionScore = _mm_set1_epi32 (prf.I);
	register __m128i __DeletionScore  = _mm_set1_epi32 (prf.D);
	
	__MatchScore     = _mm_add_epi32(__MatchScore,     prf.Match.Vector);
	__DeletionScore  = _mm_add_epi32(__DeletionScore,  prf.Deletion.Vector);
	__InsertionScore = _mm_add_epi32(__InsertionScore, prf.Insertion.Vector);
	
	const union sIOP * restrict IOP_R;
	union sIOP * restrict IOP_W = (union sIOP *) data->WorkStorage;
	
	const size_t ReadLength = data->TagLength;
	const size_t MatrixLD = data->MatrixLD;
	unsigned char * const restrict genome = data->Genome;
	const unsigned char * restrict ptrTag = data->Tag;
	assert(ReadLength > prf.AllowedJumpOut);
	const size_t BorderEnd = ReadLength - prf.AllowedJumpOut;
	/*
	 * Initialize Insertion and Match Entrance Line using FirstSequenceProtein
	 */  
	{
		__m128i * const restrict MatrixPtr = (__m128i*) data->MatrixStorage;
		/*
		 * PROFILE COLUMN 0 entrance
		 */
		IOP_W[0].Element.Match     = RowInit[MATCH];
		IOP_W[0].Element.Insertion = RowInit[INSERTION];
		KOPD                       = _mm_set1_epi32(RowInit[DELETION]);
		
		// Store all scores to Matrix
		operation(__FirstColumn, MatrixPtr, 0, 0, MatrixLD);
		
		
		register union sIOP * restrict pIOP = &IOP_W[1];
		
		/*
		 * LOOP THROUGH THE REST OF THE PROFILE
		 */
		size_t iprf = 1;
		
		while(iprf<=prf.AllowedJumpIn ) {
			KOPD = _mm_and_si128(KOPD, __ClearMask);
			
			// Add KD to Transitions
			__m128i __TransitionsD = _mm_add_epi32(__DeletionScore, KOPD);
			__TransitionsD         = _mm_max_epi32(__TransitionsD, __FirstColumn);
			
			KOPD = _mm_shuffle_epi32(__TransitionsD, (DELETION << 6) | (DELETION << 4) | (DELETION << 2) | DELETION);
			
			// Store IOPI and IOPM
			StoreMatchInsertion( &(pIOP->mm), (__m128) __TransitionsD);
			pIOP++;
			
			// Store all scores to Matrix
			operation(__TransitionsD, MatrixPtr, 0, iprf++, MatrixLD);     
		}
		
		while(iprf<=ReadLength ) {
			KOPD = _mm_and_si128(KOPD, __ClearMask);
			
			// Add KD to Transitions
			__m128i __TransitionsD = _mm_add_epi32(__DeletionScore, KOPD);
			
			KOPD = _mm_shuffle_epi32(__TransitionsD, (DELETION << 6) | (DELETION << 4) | (DELETION << 2) | DELETION);
			
			// Store IOPI and IOPM
			StoreMatchInsertion( &(pIOP->mm), (__m128) __TransitionsD);
			pIOP++;
			
			// Store all scores to Matrix
			operation(__TransitionsD, MatrixPtr, 0, iprf++, MatrixLD);     
		}
	}
		
	// Swap and assign Read and write pointers
	IOP_R = IOP_W;
	IOP_W = (union sIOP*) (((uintptr_t) &((union sIOP*)data->WorkStorage)[MatrixLD] + 63) & ~63);
	
	__m128i * const restrict MatrixPtr = ((__m128i*) data->MatrixStorage) + MatrixLD;
	/*
	 * LOOP THROUGH THE SEQUENCE STRING
	 */
	
	for ( size_t iseq=0; iseq < data->GenomeLength; ++iseq) {
		const __m128i __SequenceDNA = _mm_set1_epi8(genome[iseq]);
		__m128i KOPM = _mm_set1_epi32(IOP_R[0].Element.Match);
		/*
		 * PROFILE COLUMN 0 entrance
		 */
		{
			KOPD = __KOPDInit;
			
			// Store IOPI and IOPM
			StoreMatchInsertion( &(IOP_W[0].mm), (__m128) __FirstColumn);
			
			// Store all scores to Matrix
			operation(__FirstColumn, MatrixPtr, iseq, 0, MatrixLD);
		}
		
		/*
		 * LOOP THROUGH THE REST OF THE PROFILE
		 */
		size_t iprf = 1;
		
		/* Allowed Jump In range */
		__m128i __MismatchScore = _mm_add_epi32(_mm_set1_epi32(prf.mB), prf.Match.Vector);
		while ( iprf<=prf.AllowedJumpIn) {
			const __m128i __ReadSequence = _mm_set1_epi8(ptrTag[iprf-1]);
			const __m128i ReadIsN        = _mm_cmpeq_epi8(__ReadSequence, __N);
			__m128i DoWeMatch            = _mm_cmpeq_epi8(__SequenceDNA, __ReadSequence);
			DoWeMatch                    = _mm_or_si128(DoWeMatch, ReadIsN);
			const __m128i __AddtoM       = _mm_blendv_epi8(__MismatchScore, __MatchScore, DoWeMatch);
			
			__m128i __KI = _mm_set1_epi32(IOP_R[iprf].Element.Insertion);
			
			KOPD = _mm_and_si128(KOPD, __ClearMask);
			KOPM = _mm_and_si128(KOPM, __ClearMask);
			__KI = _mm_and_si128(__KI, __ClearMask);
			
			__m128i __KM = _mm_add_epi32(KOPM, __AddtoM);
			__KI = _mm_add_epi32(__KI, __InsertionScore);
			
			KOPM = _mm_set1_epi32(IOP_R[iprf].Element.Match);
			
			__m128i __max = _mm_max_epi32(__KM, __KI);
			
			__m128i __KD = _mm_add_epi32(KOPD, __DeletionScore);
			__KD = _mm_max_epi32(__KD, __FirstColumn);
			
			__max = _mm_max_epi32(__max, __KD);
			
			
			StoreMatchInsertion( &IOP_W[iprf].mm, (__m128) __max);
			operation(__max, MatrixPtr, iseq, iprf++, MatrixLD);
			
			KOPD = _mm_shuffle_epi32(__max, (DELETION << 6) | (DELETION << 4) | (DELETION << 2) | DELETION);
		}
		
		/* Core range */
		__MismatchScore = _mm_add_epi32(_mm_set1_epi32(prf.m), prf.Match.Vector);
		while ( iprf<BorderEnd) {
			const __m128i __ReadSequence = _mm_set1_epi8(ptrTag[iprf-1]);
			const __m128i ReadIsN        = _mm_cmpeq_epi8(__ReadSequence, __N);
			__m128i DoWeMatch            = _mm_cmpeq_epi8(__SequenceDNA, __ReadSequence);
			DoWeMatch                    = _mm_or_si128(DoWeMatch, ReadIsN);
			const __m128i __AddtoM       = _mm_blendv_epi8(__MismatchScore, __MatchScore, DoWeMatch);
			
			__m128i __KI = _mm_set1_epi32(IOP_R[iprf].Element.Insertion);
			
			KOPD = _mm_and_si128(KOPD, __ClearMask);
			KOPM = _mm_and_si128(KOPM, __ClearMask);
			__KI = _mm_and_si128(__KI, __ClearMask);
			
			__m128i __KM = _mm_add_epi32(KOPM, __AddtoM);
			__KI = _mm_add_epi32(__KI, __InsertionScore);
			
			KOPM = _mm_set1_epi32(IOP_R[iprf].Element.Match);
			
			__m128i __max = _mm_max_epi32(__KM, __KI);
			
			__m128i __KD = _mm_add_epi32(KOPD, __DeletionScore);
			
			__max = _mm_max_epi32(__max, __KD);
			
			
			StoreMatchInsertion( &IOP_W[iprf].mm, (__m128) __max);
			operation(__max, MatrixPtr, iseq, iprf++, MatrixLD);
			
			KOPD = _mm_shuffle_epi32(__max, (DELETION << 6) | (DELETION << 4) | (DELETION << 2) | DELETION);
		}
		
		/* Allowed Jump Out range */
		__MismatchScore = _mm_add_epi32(_mm_set1_epi32(prf.mE), prf.Match.Vector);
		while ( iprf<=ReadLength) {
			const __m128i __ReadSequence = _mm_set1_epi8(ptrTag[iprf-1]);
			const __m128i ReadIsN        = _mm_cmpeq_epi8(__ReadSequence, __N);
			__m128i DoWeMatch            = _mm_cmpeq_epi8(__SequenceDNA, __ReadSequence);
			DoWeMatch                    = _mm_or_si128(DoWeMatch, ReadIsN);
			const __m128i __AddtoM       = _mm_blendv_epi8(__MismatchScore, __MatchScore, DoWeMatch);
			
			__m128i __KI = _mm_set1_epi32(IOP_R[iprf].Element.Insertion);
			
			KOPD = _mm_and_si128(KOPD, __ClearMask);
			KOPM = _mm_and_si128(KOPM, __ClearMask);
			__KI = _mm_and_si128(__KI, __ClearMask);
			
			__m128i __KM = _mm_add_epi32(KOPM, __AddtoM);
			__KI = _mm_add_epi32(__KI, __InsertionScore);
			
			KOPM = _mm_set1_epi32(IOP_R[iprf].Element.Match);
			
			__m128i __max = _mm_max_epi32(__KM, __KI);
			
			__m128i __KD = _mm_add_epi32(KOPD, __DeletionScore);
			
			__max = _mm_max_epi32(__max, __KD);
			
			
			StoreMatchInsertion( &IOP_W[iprf].mm, (__m128) __max);
			operation(__max, MatrixPtr, iseq, iprf++, MatrixLD);
			
			KOPD = _mm_shuffle_epi32(__max, (DELETION << 6) | (DELETION << 4) | (DELETION << 2) | DELETION);
		}
		

		// Swap Read and Write pointers
		const register union sIOP * const ptr = IOP_W;
		IOP_W = (union sIOP *) IOP_R;
		IOP_R = ptr;
	}
	return 0;
}
//---------------------------------------------------------------

static int GetAlignmentSequenceBorder(GlobalData_t * const restrict data)
{
	fputs("GetAlignmentSequenceBorder no yet implemented\n", stderr);
	return 1;
}
//---------------------------------------------------------------

static int GetStateSequenceBorder(GlobalData_t * const restrict data)
{
	const int (* restrict Scores)[4]  = (const int (* restrict )[4]) data->MatrixStorage;
	const size_t ld = data->MatrixLD;
	int index_r = 0, index_c = 0;
	{
		const size_t prfLength = data->TagLength;
		int extrema = NLOW;
		const size_t SeqLength = data->GenomeLength;
		assert(prfLength > prf.AllowedJumpOut);
		const size_t BorderEnd = prfLength - prf.AllowedJumpOut;
		for (size_t i=1; i<=SeqLength; ++i) {
			for (size_t j=BorderEnd; j<=prfLength; j++) {
				const int val = ToScore(Scores[i*ld+j][EXTRA]);
				if (val > extrema)
				{
					extrema = val;
					index_r = i;
					index_c = j;
				}
			}
		}
		
		data->score = extrema;
	  data->AlignmentRange[1] = index_r;
		data->ProfileRange[1] = index_c;
	}
	unsigned int State = EXTRA;
	
	const register int SoftClipBoundary = data->SoftClipBoundary;
	unsigned int SoftClipMismatchCount = 0U;
	unsigned int MismatchCount  = 0U;
	// 						unsigned int InsertionCount = 0U;
	unsigned int RequiredMMSize = 0U;
	unsigned int SoftClipEvictedMMSize = 0U;
	unsigned int RequiredCigSize = 0U;
	unsigned int SoftClipEvictedCigSize = 0U;
	
	unsigned char * restrict ptrState = data->StateSequence + STATE_MEMORY_SIZE - 1;
	const unsigned char * const restrict StateLimit = data->StateSequence; 
	const unsigned char * const restrict genome = data->Genome; 
	const unsigned char * restrict ptrTag = data->Tag;
	unsigned int previous_State = ((unsigned int) Scores[index_r*ld+index_c][EXTRA] & STATE_MASK) >> STATE_SHIFT;
	*ptrState = '\0';
	
	int ret = -1;

	while ( index_c >= 0 && index_r >= 0)
	{
		if (--ptrState < StateLimit) {
			fputs("Allocated state sequence length is not sufficient\n", stderr);
			goto OUT;
		}
		const unsigned int lScore = (unsigned int) Scores[index_r*ld+index_c][State];
		const unsigned int current_State = (lScore & STATE_MASK) >> STATE_SHIFT;
		if (current_State != previous_State) {
			RequiredCigSize++;
			previous_State = current_State;
		}
		
		switch(current_State) {
			case PRIORITY_MATCH:
				if (lScore & MISMATCH_MASK) {
					if (genome[index_r] == 'N' && ptrTag[index_c] != 'N') {
						*ptrState = 'n';
						RequiredMMSize++;
					}
					else if (ptrTag[index_c] == 'N') {
						*ptrState = 'N';
						RequiredMMSize++;
					}
					else {
						*ptrState = 'M';
					}
				}
				else {
					*ptrState = 'm';
					MismatchCount++;
					RequiredMMSize++;
				}
				State = MATCH;
				--index_r;
				--index_c;
				break;
			case PRIORITY_INSERTION:
				--index_r;
				*ptrState = 'D';
				State = INSERTION;
				break;
			case PRIORITY_DELETION:
				*ptrState = 'I';
				--index_c;
				RequiredMMSize++;
				State = DELETION;
				break;
			case PRIORITY_EXTRA:
				goto DONE;
				break;
			default:
				fprintf(stderr,"\nUnknown state (%u) encountered from cell (%i,%i)\n", current_State, index_r, index_c);
				//--index_c;
				goto OUT;
				break;
		}
		
		if (index_c == SoftClipBoundary) {
			SoftClipMismatchCount  = MismatchCount;
			SoftClipEvictedMMSize  = RequiredMMSize;
			SoftClipEvictedCigSize = RequiredCigSize;
		}
	}
	DONE:;
	ret = 0;
	
	data->MismatchCount = MismatchCount;
	data->SoftClipMismatchCount = SoftClipMismatchCount;
	data->RequiredMMSize = RequiredMMSize;
	data->RequiredCigSize = RequiredCigSize;
	data->SoftClipEvictedMMSize = SoftClipEvictedMMSize;
	data->SoftClipEvictedCigSize = SoftClipEvictedCigSize;
	data->States = ptrState + 1; // avoid the last move
	data->AlignmentRange[0] = index_r;
	data->ProfileRange[0] = index_c;
	
	OUT:
	return ret;
}
//---------------------------------------------------------------

ComputeFunctions_t cpu_std_border = {
	.createMatrix = CreateMatrixBorder,
	.getBestAlignment = GetBestAlignment,
	.getStateSequence = GetStateSequenceBorder,
	.computeCigarAndMismatch = ComputeCigarAndMismatch,
	.allocStorage = allocStorage,
	.freeStorage = freeStorage
};


static int allocStorageVector(GlobalData_t * const restrict data)
{
	data->MatrixLD = (data->TagLength + ((ALIGNMENT/sizeof(int))-1) ) & ~( (ALIGNMENT/sizeof(int)) - 1 );
	
	data->WorkStorage = (void*) _mm_malloc(2*data->MatrixLD*sizeof(int), ALIGNMENT);
	if (data->WorkStorage) {
		data->MatrixStorage = (void*) _mm_malloc(3*data->MatrixLD*data->GenomeLength*sizeof(int), ALIGNMENT);
		if (data->MatrixStorage) return 0;
	}
	return -1;
}
//---------------------------------------------------------------

static int CreateMatrixVector(GlobalData_t * const restrict data)
{
	static const unsigned int BestAlignmentMasks[4][4] __attribute__((aligned(ALIGNMENT))) = {
		{ -1,  0,  0,  0},
		{  0, -1,  0,  0},
		{  0,  0, -1,  0},
		{  0,  0,  0, -1}
	};
	const __m128i __ClearMask   = _mm_set1_epi32(CLEAR_MASK);
	const __m128i __Broadcast_D = _mm_set1_epi32(DeletionIndex);
	const __m128i __N           = _mm_set1_epi32((int)N);
	
	const __m128i __MI       = _mm_set1_epi32(prf.Match.To[INSERTION]);
	const __m128i __MD       = _mm_set1_epi32(prf.Match.To[DELETION]);
	const __m128i __MM       = _mm_set1_epi32(prf.Match.To[MATCH]);
	const __m128i __MX       = _mm_set1_epi32(prf.Match.To[EXTRA]);
	
	const __m128i __II       = _mm_set1_epi32(prf.I + prf.Insertion.To[INSERTION]);
	const __m128i __IM       = _mm_set1_epi32(prf.I + prf.Insertion.To[MATCH]);
	const __m128i __IX       = _mm_set1_epi32(prf.I + prf.Insertion.To[EXTRA]);
	const __m128i __DD       = _mm_set1_epi32(prf.Deletion.To[DELETION]);
	
	const __m128i __DM_minus_DD = _mm_sub_epi32(_mm_set1_epi32(prf.Deletion.To[MATCH]), __DD);
	const __m128i __DX_minus_DD = _mm_sub_epi32(_mm_set1_epi32(prf.Deletion.To[EXTRA]), __DD);
	
	const __m128i __NLOW          = _mm_set1_epi32(NLOW << SCORE_SHIFT);
	const __m128i __KOPMInit      = _mm_set1_epi32(RowInit[MATCH]);
	const __m128i __MismatchScore = _mm_set1_epi32(prf.m);
	const __m128i  __MatchScore   = _mm_set1_epi32(prf.M);
	const int KOPDInit = RowInit[DELETION];
	
	const size_t MatrixLD = data->MatrixLD;
	const size_t ReadLength = data->TagLength;
	const unsigned char * const restrict genome = data->Genome;
	const unsigned char * restrict ptrTag = data->Tag; 
	
	int * const restrict WORK_match = (int*) data->WorkStorage;
	int * const restrict WORK_insertion = ((int*)data->WorkStorage) + MatrixLD;
	
	int * const MatrixMatch = data->MatrixStorage;
	int * const MatrixInsertion = ((int*) data->MatrixStorage) + data->GenomeLength*MatrixLD;
	int * const MatrixDeletion = ((int*) data->MatrixStorage) + 2*(data->GenomeLength*MatrixLD);
	
	size_t BestAlignment = 0;
	__m128i __BestAlignmentValue = __NLOW;
	/* Fills in the Matrix */
	{
		/*
		 * Initialize Insertion and Match Entrance Line
		 */  
		{
			/*
			 * PROFILE COLUMN 0 entrance
			 */
			WORK_match[0]     = RowInit[MATCH];
			WORK_insertion[0] = RowInit[INSERTION];
			int KOPD          = 0;
			
			/*
			 * LOOP THROUGH THE REST OF THE PROFILE
			 */
			for (size_t iprf=1; iprf<ReadLength; ++iprf ) {
				KOPD                += prf.D;
				WORK_match[iprf]     = KOPD + prf.Deletion.To[MATCH];
				WORK_insertion[iprf] = KOPD + prf.Deletion.To[INSERTION];
				KOPD                += prf.Deletion.To[DELETION];
			}
		}

		const size_t lReadLength = ReadLength - 4;
		const __m128i __BestAlignmentMask = _mm_load_si128((const __m128i*) BestAlignmentMasks[(ReadLength-1) & 0x3]);
		/*
		 * LOOP THROUGH THE SEQUENCE STRING
		 */
		int * restrict MatrixMatchPtr     =  MatrixMatch;
		int * restrict MatrixInsertionPtr =  MatrixInsertion;
		int * restrict MatrixDeletionPtr  =  MatrixDeletion;
		for ( size_t iseq=0; iseq < kGAPPED_ALIGN_GENOME_LENGTH; ++iseq) {
			const __m128i __SequenceDNA = _mm_set1_epi32((int)genome[iseq]);
			const __m128i SeqIsN = _mm_cmpeq_epi32(__SequenceDNA, __N);
			
			/*
			 * PROFILE COLUMN 0 entrance
			 */
			int KOPD = KOPDInit;
			__m128i __previousKOPM = __KOPMInit;
			/*
			 * LOOP THROUGH THE REST OF THE PROFILE
			 */
			size_t iprf;
			for (iprf=0; iprf<lReadLength; iprf+=4 ) {
				const __m128i __ReadSequence = _mm_cvtepi8_epi32(*((__m128i*)&ptrTag[iprf]));
				const __m128i __ReadIsN      = _mm_cmpeq_epi32(__ReadSequence, __N);
				__m128i   DoWeMatch          = _mm_cmpeq_epi32(__SequenceDNA, __ReadSequence);
				DoWeMatch                    = _mm_or_si128(DoWeMatch, _mm_or_si128(SeqIsN, __ReadIsN));
				const __m128i __AddtoM       = _mm_blendv_epi8 (__MismatchScore, __MatchScore, DoWeMatch);
				
				/* Load the previous line insertion vector */
				__m128i __KOPI = _mm_load_si128((__m128i*)&WORK_insertion[iprf]);
				
				/* Load the previous line match vector */
				__m128i __tempKOPM = _mm_load_si128((__m128i*)&WORK_match[iprf]);
				
				/* shift KOPM 1 position inserting column 0 data
				 *    |                   KOPM                        |          KOPMInit or previous KOPM            | 
				 *    |15 14 13 12 11 10 09 08 07 06 05 04 03 02 01 00|15 14 13 12 11 10 09 08 07 06 05 04 03 02 01 00|
				 * RowInit					* ->    |-----------------------------------------------|
				 *       shift 15                                     |-----------------------------------------------|
				 */
				__m128i __KOPM = _mm_alignr_epi8(__tempKOPM, __previousKOPM, 12);
				__previousKOPM = __tempKOPM;
				
				/* clean past states */
				__KOPI = _mm_and_si128(__KOPI, __ClearMask);
				__KOPM = _mm_and_si128(__KOPM, __ClearMask);
				
				
				/* treats mismatch and match */
				__KOPM = _mm_add_epi32(__KOPM, __AddtoM);
				
				/* compute new insertions: I + II, M + MI and store it*/
				const __m128i __tempI_II = _mm_add_epi32(__KOPI, __II);
				const __m128i __tempM_MI = _mm_add_epi32(__KOPM, __MI);
				const __m128i __tempwrite =  _mm_max_epi32(__tempI_II, __tempM_MI);
				_mm_store_si128((__m128i*)&WORK_insertion[iprf], __tempwrite);
				_mm_store_si128((__m128i*)&MatrixInsertionPtr[iprf], __tempwrite);
				
				
				/* compute partial new matches: I + IM, M + MM */
				const __m128i __tempI_IM = _mm_add_epi32(__KOPI, __IM);
				const __m128i __tempM_MM = _mm_add_epi32(__KOPM, __MM);
				const __m128i __Matches  = _mm_max_epi32(__tempI_IM, __tempM_MM);
				
				/* treats deletion */
				const __m128i __tempM_MD = _mm_add_epi32(__KOPM, __MD);
				
				__m128i __Deletions;
#if 0
			{
				__m128i __KOPD = _mm_set1_epi32(KOPD);
				__m128i __D_DD = _mm_add_epi32(_mm_set1_epi32(prf.D), __DD);
				__m128i __mask = _mm_set1_epi32(0xFFFFFFFF);
				for (int i=0; i<4; i++) {
					__KOPD           = _mm_and_si128(__KOPD, __ClearMask);
					__KOPD           = _mm_add_epi32(__KOPD, __D_DD);
					__m128i __bigger = _mm_cmpgt_epi32(__tempM_MD, __KOPD); 
					__bigger         = _mm_and_si128(__bigger, __mask);
					__m128i __temp   = _mm_blendv_epi8(__KOPD, __tempM_MD, __bigger);
					__KOPD           = _mm_alignr_epi8(__temp, __NLOW, 12);
					__Deletions      = _mm_blendv_epi8(__Deletions, __temp, __mask);
					__mask           = _mm_slli_si128(__mask, 4);
				}
				_mm_store_si128( (__m128i*) &MatrixDeletionPtr[iprf], __Deletions);
				KOPD = _mm_extract_epi32(__KOPD, 3);
			}
#else
			{
				int * const restrict MDPtr = &MatrixDeletionPtr[iprf];
				_mm_store_si128((__m128i*)MDPtr, __tempM_MD);
				const int DD = prf.D + prf.Deletion.To[DELETION];
				
				for (int i=0; i< 4; i++) {
					KOPD &= CLEAR_MASK;
					KOPD += DD;
					int tmp = MDPtr[i];
					tmp = (KOPD > tmp ) ? KOPD : tmp;
					KOPD = tmp;
					MDPtr[i] = tmp;
				}
				__Deletions = _mm_load_si128((__m128i*)MDPtr);
			}
// 			KOPD = MDPtr[3]; implicit 
#endif

				/* finalize matches */
				const __m128i __FinalMatch = _mm_max_epi32(__Matches, _mm_add_epi32(__Deletions, __DM_minus_DD));
				
				/*store matches */
				_mm_store_si128((__m128i*)&WORK_match[iprf], __FinalMatch);
				_mm_store_si128((__m128i*)&MatrixMatchPtr[iprf], __FinalMatch);
			}

			/* Perform last bulk of columns, computing best alignment found */
			{
				const __m128i __ReadSequence = _mm_cvtepi8_epi32(*((__m128i*)&ptrTag[iprf]));
				const __m128i __ReadIsN      = _mm_cmpeq_epi32(__ReadSequence, __N);
				__m128i   DoWeMatch          = _mm_cmpeq_epi32(__SequenceDNA, __ReadSequence);
				DoWeMatch                    = _mm_or_si128(DoWeMatch, _mm_or_si128(SeqIsN, __ReadIsN));
				const __m128i __AddtoM       = _mm_blendv_epi8 (__MismatchScore, __MatchScore, DoWeMatch);
				
				/* Load the previous line insertion vector */
				__m128i __KOPI = _mm_load_si128((__m128i*)&WORK_insertion[iprf]);
				
				/* Load the previous line match vector */
				__m128i __tempKOPM = _mm_load_si128((__m128i*)&WORK_match[iprf]);
				
				/* shift KOPM 1 position inserting column 0 data
				 *    |                   KOPM                        |          KOPMInit or previous KOPM            | 
				 *    |15 14 13 12 11 10 09 08 07 06 05 04 03 02 01 00|15 14 13 12 11 10 09 08 07 06 05 04 03 02 01 00|
				 * ->    |-----------------------------------------------|
				 *       shift 15                                     |-----------------------------------------------|
				 */
				__m128i __KOPM = _mm_alignr_epi8(__tempKOPM, __previousKOPM, 12);
				__previousKOPM = __tempKOPM;
				
				/* clean past states */
				__KOPI = _mm_and_si128(__KOPI, __ClearMask);
				__KOPM = _mm_and_si128(__KOPM, __ClearMask);
				
				
				/* treats mismatch and match */
				__KOPM = _mm_add_epi32(__KOPM, __AddtoM);
				
				/* compute new insertions: I + II, M + MI and store it*/
				const __m128i __tempI_II = _mm_add_epi32(__KOPI, __II);
				const __m128i __tempM_MI = _mm_add_epi32(__KOPM, __MI);
				const __m128i __tempwrite =  _mm_max_epi32(__tempI_II, __tempM_MI);
				const __m128i __tempI_IX = _mm_add_epi32(__KOPI, __IX);
				const __m128i __tempM_MX = _mm_add_epi32(__KOPM, __MX);
				_mm_store_si128((__m128i*)&WORK_insertion[iprf], __tempwrite);
				_mm_store_si128((__m128i*)&MatrixInsertionPtr[iprf], __tempwrite);
				__m128i __maxOut =  _mm_max_epi32(__tempI_IX, __tempM_MX);
				
				/* compute partial new matches: I + IM, M + MM */
				const __m128i __tempI_IM = _mm_add_epi32(__KOPI, __IM);
				const __m128i __tempM_MM = _mm_add_epi32(__KOPM, __MM);
				const __m128i __Matches  = _mm_max_epi32(__tempI_IM, __tempM_MM);
				
				/* treats deletion */
				const __m128i __tempM_MD = _mm_add_epi32(__KOPM, __MD);
				
				__m128i __Deletions;
#if 0
			{
				__m128i __KOPD = _mm_set1_epi32(KOPD);
				__m128i __D_DD = _mm_add_epi32(_mm_set1_epi32(prf.D), __DD);
				__m128i __mask = _mm_set1_epi32(0xFFFFFFFF);
				for (int i=0; i<4; i++) {
					__KOPD           = _mm_and_si128(__KOPD, __ClearMask);
					__KOPD           = _mm_add_epi32(__KOPD, __D_DD);
					__m128i __bigger = _mm_cmpgt_epi32(__tempM_MD, __KOPD); 
					__bigger         = _mm_and_si128(__bigger, __mask);
					__m128i __temp   = _mm_blendv_epi8(__KOPD, __tempM_MD, __bigger);
					__KOPD           = _mm_alignr_epi8(__temp, __NLOW, 12);
					__Deletions      = _mm_blendv_epi8(__Deletions, __temp, __mask);
					__mask           = _mm_slli_si128(__mask, 4);
				}
				_mm_store_si128( (__m128i*) &MatrixDeletionPtr[iprf], __Deletions);
				KOPD = _mm_extract_epi32(__KOPD, 3);
			}
#else
			{
				int * const restrict MDPtr = &MatrixDeletionPtr[iprf];
				_mm_store_si128((__m128i*)MDPtr, __tempM_MD);
				const int DD = prf.D + prf.Deletion.To[DELETION];
				
				for (int i=0; i< 4; i++) {
					KOPD &= CLEAR_MASK;
					KOPD += DD;
					int tmp = MDPtr[i];
					tmp = (KOPD > tmp ) ? KOPD : tmp;
					KOPD = tmp;
					MDPtr[i] = tmp;
				}
				__Deletions = _mm_load_si128((__m128i*)MDPtr);
			}
#endif
																					
			/* finalize matches */
			const __m128i __FinalMatch = _mm_max_epi32(__Matches, _mm_add_epi32(__Deletions, __DM_minus_DD));
			__maxOut                   = _mm_max_epi32(__maxOut, _mm_add_epi32(__Deletions, __DX_minus_DD));
			
			/*store matches */
			_mm_store_si128((__m128i*)&WORK_match[iprf], __FinalMatch);
			_mm_store_si128((__m128i*)&MatrixMatchPtr[iprf], __FinalMatch);
			
			/* update best alignement */
			{
				const __m128i __test = _mm_cmpgt_epi32(__BestAlignmentValue, __maxOut); 
				const int allzeros = _mm_test_all_zeros (__test, __BestAlignmentMask);
				__BestAlignmentValue = _mm_blendv_epi8(__maxOut, __BestAlignmentValue, __test);
				BestAlignment = ( allzeros ) ? iseq : BestAlignment;
			}
		}

			MatrixDeletionPtr  += MatrixLD;
			MatrixInsertionPtr += MatrixLD;
			MatrixMatchPtr     += MatrixLD;
		}
	}
		
	/* Place best state within best alignemnt */
	{
		unsigned int itmp[4] __attribute__((aligned(ALIGNMENT)));
		_mm_store_si128((__m128i*) itmp, __BestAlignmentValue);
		register size_t BestOut = ((size_t) (itmp[(ReadLength-1) & 0x3])) << 32;
		BestAlignment |= BestOut;
	}

	return BestAlignment;
}
//---------------------------------------------------------------

ComputeFunctions_t cpu_vector = {
	.createMatrix = CreateMatrixVector,
	.getBestAlignment = NULL,
	.getStateSequence = NULL,
	.computeCigarAndMismatch = NULL,
	.allocStorage = allocStorageVector,
	.freeStorage = freeStorage
};      

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
