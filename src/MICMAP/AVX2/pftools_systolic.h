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
/* Check that only one is active */
#if defined(__DECLARE_VARIABLES__)
#  if defined(__EXTRA_FUNCTIONS__) || defined(__ALLOCATE_MEMORY__) || defined(__FILL_MATRIX__) || defined(__GET_ALIGNMENT__) || defined(__FREE_MEMORY__) || defined (__GET_ALIGNMENT_STATE__)
#    error "Too many preprocessing variables defined"
#  endif
#elif defined(__EXTRA_FUNCTIONS__)
#  if defined(__DECLARE_VARIABLES__) || defined(__ALLOCATE_MEMORY__) || defined(__FILL_MATRIX__) || defined(__GET_ALIGNMENT__) || defined(__FREE_MEMORY__) || defined (__GET_ALIGNMENT_STATE__)
#    error "Too many preprocessing variables defined"
#  endif
#elif defined(__ALLOCATE_MEMORY__)
#  if defined(__EXTRA_FUNCTIONS__) || defined(__DECLARE_VARIABLES__) || defined(__FILL_MATRIX__) || defined(__GET_ALIGNMENT__) || defined(__FREE_MEMORY__) || defined (__GET_ALIGNMENT_STATE__)
#    error "Too many preprocessing variables defined"
#  endif
#elif defined(__FILL_MATRIX__)
#  if defined(__EXTRA_FUNCTIONS__) || defined(__ALLOCATE_MEMORY__) || defined(__DECLARE_VARIABLES__) || defined(__GET_ALIGNMENT__) || defined(__FREE_MEMORY__) || defined (__GET_ALIGNMENT_STATE__)
#    error "Too many preprocessing variables defined"
#  endif
#elif defined(__GET_ALIGNMENT__)
#  if defined(__EXTRA_FUNCTIONS__) || defined(__ALLOCATE_MEMORY__) || defined(__FILL_MATRIX__) || defined(__DECLARE_VARIABLES__) || defined(__FREE_MEMORY__) || defined (__GET_ALIGNMENT_STATE__)
#    error "Too many preprocessing variables defined"
#  endif
#elif defined(__GET_ALIGNMENT_STATE__)
#  if defined(__EXTRA_FUNCTIONS__) || defined(__ALLOCATE_MEMORY__) || defined(__FILL_MATRIX__) || defined(__DECLARE_VARIABLES__) || defined(__FREE_MEMORY__) || defined (__GET_ALIGNMENT__)
#    error "Too many preprocessing variables defined"
#  endif
#elif defined(__FREE_MEMORY__)
#  if defined(__EXTRA_FUNCTIONS__) || defined(__ALLOCATE_MEMORY__) || defined(__FILL_MATRIX__) || defined(__GET_ALIGNMENT__) || defined(__DECLARE_VARIABLES__) || defined (__GET_ALIGNMENT_STATE__)
#    error "Too many preprocessing variables defined"
#  endif
#endif
//
// -------------------------------------------------------
// DECLARE VARIABLE TO BE GLOBALLY USED WITH ALIGNMENT
// -------------------------------------------------------
//

// WARNING: With the addition of a third state not to account N in genome, the tag has to be set to 0 when MATCH
//          and 1 when mismatch, otherwise we must add instruction and cannot use masks.

#ifdef __DECLARE_STATIC_VARIABLES__
	static const short int At[16] __attribute__((aligned(64))) = {16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1};
	static const union { char i[32]; __m256i ymm; } reverse = { .i = {14,15,12,13,10,11,8,9,6,7,4,5,2,3,0,1, 
																																    14,15,12,13,10,11,8,9,6,7,4,5,2,3,0,1}};
	static const union { short int i[16]; __m256i ymm; } OnlySecond = { .i = {0,0xFFFF,0,0,0,0,0,0,0,0,0,0,0,0,0,0} };
	static const short int __N[16] __attribute__((aligned(32))) = { (short int)'N',(short int)'N',(short int)'N',(short int)'N',
	                                                                (short int)'N',(short int)'N',(short int)'N',(short int)'N',
	                                                                (short int)'N',(short int)'N',(short int)'N',(short int)'N',
	                                                                (short int)'N',(short int)'N',(short int)'N',(short int)'N' };

	static const short int Globals[16] __attribute__((aligned(32))) = {
		(_M<<SCORE_SHIFT)+(PRIORITY_MATCH    <<STATE_SHIFT)+(0<<MISMATCH_SHIFT),
		(_m<<SCORE_SHIFT)+(PRIORITY_MATCH    <<STATE_SHIFT)+(1<<MISMATCH_SHIFT),
		(_D<<SCORE_SHIFT)+(PRIORITY_DELETION <<STATE_SHIFT)+(1<<MISMATCH_SHIFT),
		(_I<<SCORE_SHIFT)+(PRIORITY_INSERTION<<STATE_SHIFT)+(1<<MISMATCH_SHIFT),
		(_MD<<SCORE_SHIFT),
		(_MX<<SCORE_SHIFT),
		(_MM<<SCORE_SHIFT),
		(_MI<<SCORE_SHIFT),
		((_D+_DD)<<SCORE_SHIFT)+(PRIORITY_DELETION <<STATE_SHIFT)+(1<<MISMATCH_SHIFT),
		((_D+_DX)<<SCORE_SHIFT)+(PRIORITY_DELETION <<STATE_SHIFT)+(1<<MISMATCH_SHIFT),
		((_D+_DM)<<SCORE_SHIFT)+(PRIORITY_DELETION <<STATE_SHIFT)+(1<<MISMATCH_SHIFT),
		((_D+_DI)<<SCORE_SHIFT)+(PRIORITY_DELETION <<STATE_SHIFT)+(1<<MISMATCH_SHIFT),
		((_I+_ID)<<SCORE_SHIFT)+(PRIORITY_INSERTION<<STATE_SHIFT)+(1<<MISMATCH_SHIFT),
		((_I+_IX)<<SCORE_SHIFT)+(PRIORITY_INSERTION<<STATE_SHIFT)+(1<<MISMATCH_SHIFT),
		((_I+_IM)<<SCORE_SHIFT)+(PRIORITY_INSERTION<<STATE_SHIFT)+(1<<MISMATCH_SHIFT),
		((_I+_II)<<SCORE_SHIFT)+(PRIORITY_INSERTION<<STATE_SHIFT)+(1<<MISMATCH_SHIFT)
	};
#define M_INDEX 0
#define m_INDEX 1
#define D_INDEX 2
#define I_INDEX 3
#define MD_INDEX 4
#define MM_INDEX 6
#define MI_INDEX 7
#define D_DD_INDEX 8
#define D_DM_INDEX 10
#define I_IM_INDEX 14
#define I_II_INDEX 15
#endif

#ifdef __DECLARE_VARIABLES__
	short int results[16]  __attribute__((aligned(32)));
	short int Id[16]   __attribute__((aligned(32)));
	unsigned short int MorD[16] __attribute__((aligned(32)));

	const __m256i __ClearMask    = _mm256_set1_epi16(CLEAR_MASK);
	const __m256i __DeletionMask = _mm256_set1_epi16(PRIORITY_DELETION <<STATE_SHIFT);
	const __m256i __MatchMask    = _mm256_set1_epi16(PRIORITY_MATCH <<STATE_SHIFT);

	const __m256i __MI = _mm256_set1_epi16 (Globals[MI_INDEX]);
	const __m256i __MD = _mm256_set1_epi16 (Globals[MD_INDEX]);
	const __m256i __MM = _mm256_set1_epi16 (Globals[MM_INDEX]);

	const __m256i __DD = _mm256_set1_epi16 (Globals[D_DD_INDEX]);
	const __m256i __II = _mm256_set1_epi16 (Globals[I_II_INDEX]);

	const __m256i __NLOW          = _mm256_set1_epi16(NLOW << SCORE_SHIFT);
	const __m256i __MismatchScore = _mm256_set1_epi16 (Globals[m_INDEX]);
	const __m256i __MatchScore    = _mm256_set1_epi16 (Globals[M_INDEX]);
#endif

//
// -------------------------------------------------------
// EXTRA INLINED FUNCTIONS
// -------------------------------------------------------
//
#ifdef __EXTRA_FUNCTIONS__
/* Get the memory location within IDM_Matrix*/
enum StateTypeOffset { InsertionOffset=0, DeletionOffset=16, MatchOffset=32};
static inline unsigned short int*  __attribute__((always_inline))
Location(short int * const restrict _Matrix, const size_t MatrixLD, const int iprf, const int iseq, const enum StateTypeOffset Offset)
{
	size_t BlockY       = iseq >> 0x4;
	const size_t BlockX = iprf >> 0x4;

	size_t Y = iseq & 0xF;
	size_t X = iprf & 0xF;
	const size_t Sum = X+Y; 
	if ( Sum > 15) {
		BlockY++;
		Y = Sum - 16;
	} else {
		Y = Sum;
	}
//  	fprintf(stderr, "%i %i %i -> %i\n", iseq, iprf, Offset, BlockX*MatrixLD*48 + BlockY*48 + Y*48 + X + Offset);
	return (unsigned short int*) &_Matrix[BlockX*MatrixLD*48 + BlockY*16*48 + Y*48 + X + Offset];
}
#endif

//
// -------------------------------------------------------
// ALLOCATE GLOBAL MEMORY
// -------------------------------------------------------
//
#ifdef __ALLOCATE_MEMORY__
{
  const size_t maxTag = (kTagSize + 15UL) & ~(15UL);
	IDM_Matrix = malloc_huge_pages(3*maxTag*maxMatrixLD*sizeof(int));
	if (IDM_Matrix == NULL) {
		fputs("Unable to allocate enough memory\n", stderr);
		err = -1 ; goto zeend;
	}
}
#endif

//
// -------------------------------------------------------
// FILLS IN THE MATRIX
// -------------------------------------------------------
//
#ifdef __FILL_MATRIX__
#define CONVERT_BROADCAST_UINT8(ptr) _mm256_cvtepu8_epi16(_mm_load_si128( (__m128i*) ptr ))
#define REVERSE_SEQUENCE(a) {\
		a = _mm256_shuffle_epi8( a, reverse.ymm);\
		a = _mm256_permute2x128_si256(a, a, 0x01);\
}

#define _mm256_alignr_epi16(a, b, shift)  _mm256_alignr_epi8(a,_mm256_permute2f128_si256(a,b,0x21),(const int) (2*shift))

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TwoDiagsLoop1
{
#define TwoDiagsLoop1(iloop, MatrixOffset, PastSequence, MatchLoadRegister, MatchStoreRegister) \
{\
	/* Clean past states */\
	__m256i __KOPD     = _mm256_and_si256(__Deletions,  __ClearMask);\
	__m256i __KOPI     = _mm256_and_si256(__Insertions, __ClearMask);\
	\
	MatchStoreRegister = _mm256_and_si256(__Match,      __ClearMask);\
	\
	/* Get the previous diagonal for deletion, shift it */\
	__KOPD = _mm256_alignr_epi16(__KOPD, __KOPDInit, 15);\
	\
	/* Coming from Deletion */\
	const __m256i __DtoM = _mm256_add_epi16(__KOPD, __DM);\
	const __m256i __DtoD = _mm256_add_epi16(__KOPD, __DD);\
	\
	/* Coming from Insertion */\
	const __m256i __ItoM = _mm256_add_epi16(__KOPI, __IM);\
	const __m256i __ItoI = _mm256_add_epi16(__KOPI, __II);\
	\
	__Match = _mm256_max_epi16(__DtoM, __ItoM);\
	\
	/* Get the correct Genome string, for dummy NLOW is pefect */\
	{\
		const __m256i __NewSeq = _mm256_alignr_epi16(PastSequence, __SequenceDNA, 15-iloop);\
		SeqIsN    = _mm256_cmpeq_epi16(__NewSeq, *(__m256i*)__N);\
		DoWeMatch = _mm256_or_si256(_mm256_cmpeq_epi16(__NewSeq, __ReadSequence), ReadIsN);\
		__AddtoM  = _mm256_blendv_epi8(__MismatchScore, __MatchScore, DoWeMatch);\
	}\
	/* Coming from Match */\
	__m256i __M;\
	__M = _mm256_alignr_epi16(MatchLoadRegister, __KOPMInit, 15);\
	/*__M = _mm512_mask_add_epi32(__M, SeqIsNotN, __M, __AddtoM);*/\
	__M = _mm256_add_epi16(__M, _mm256_andnot_si256(SeqIsN, __AddtoM));\
	/*__M = _mm512_mask_add_epi32(__M, _mm512_knot(SeqIsNotN), __M, __MatchMask);*/\
	__M = _mm256_add_epi16(__M, _mm256_and_si256(__MatchMask, SeqIsN));\
	const __m256i __MtoI  = _mm256_add_epi16(__M, __MI);\
	const __m256i __MtoD  = _mm256_add_epi16(__M, __MD);\
	const __m256i __MtoM  = _mm256_add_epi16(__M, __MM);\
	\
	__Insertions = _mm256_max_epi16(__ItoI, __MtoI);\
	__Deletions  = _mm256_max_epi16(__DtoD, __MtoD);\
	__Match      = _mm256_max_epi16(__Match, __MtoM);\
	\
	/* Store results in order I,D,M */\
	_mm256_store_si256((__m256i*) &IDM_Matrix[MatrixOffset+iloop*48],    __Insertions);\
	_mm256_store_si256((__m256i*) &IDM_Matrix[MatrixOffset+16+iloop*48], __Deletions);\
	_mm256_store_si256((__m256i*) &IDM_Matrix[MatrixOffset+32+iloop*48], __Match);\
}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TwoDiagsLoop2
{
#define TwoDiagsLoop2(iloop, MatrixOffset, PastSequence, MatchLoadRegister, MatchStoreRegister) \
{\
	/* Get and clean past column needs */\
	const __m256i __PreviousDeletion = _mm256_load_si256((__m256i*) (PreviousColumn+MatrixOffset+iloop*48+48+16));\
	const __m256i __PreviousMatch    = _mm256_load_si256((__m256i*) (PreviousColumn+MatrixOffset+iloop*48+32));\
	\
	/* Evict from L2 cache */\
	EVICTL1(PreviousColumn+MatrixOffset+iloop*48+48+16);\
	EVICTL1(PreviousColumn+MatrixOffset+iloop*48+32+16);\
	\
	/* Prefetch to L1 from L2*/\
	PREFETCHL1(PreviousColumn+MatrixOffset+(2+iloop)*48+48+16);\
	PREFETCHL1(PreviousColumn+MatrixOffset+(2+iloop)*48+32);\
	\
	/* Clean past states */\
	const __m256i __KOPI = _mm256_and_si256(__Insertions, __ClearMask);\
	MatchStoreRegister   = __Match;\
	\
	/* Get the previous diagonal for deletion, shift it */\
	__m256i __KOPD = _mm256_and_si256(_mm256_alignr_epi16(__Deletions, __PreviousDeletion, 15), __ClearMask);\
	\
	/* Coming from Deletion */\
	const __m256i __DtoM = _mm256_add_epi16(__KOPD, __DM);\
	const __m256i __DtoD = _mm256_add_epi16(__KOPD, __DD);\
	\
	/* Coming from Insertion */\
	const __m256i __ItoM = _mm256_add_epi16(__KOPI, __IM);\
	const __m256i __ItoI = _mm256_add_epi16(__KOPI, __II);\
	\
	__Match = _mm256_max_epi16(__DtoM, __ItoM);\
	\
	/* Get the correct Genome string, for dummy NLOW is pefect */\
	__m256i __AddtoM, SeqIsN;\
	{\
		const __m256i __NewSeq = _mm256_alignr_epi16(PastSequence, __SequenceDNA, 15-iloop);\
		SeqIsN = _mm256_cmpeq_epi16(__NewSeq,*((__m256i*)__N));\
		const __m256i DoWeMatch = _mm256_or_si256(_mm256_cmpeq_epi16(__NewSeq, __ReadSequence), ReadIsN);\
		__AddtoM  = _mm256_blendv_epi8(__MismatchScore, __MatchScore, DoWeMatch);\
	}\
	/* Coming from Match */\
	__m256i __M;\
	__M = _mm256_and_si256(_mm256_alignr_epi16(MatchLoadRegister, __PreviousMatch, 15), __ClearMask);\
	__M = _mm256_add_epi16(__M, _mm256_andnot_si256(SeqIsN, __AddtoM));\
	__M = _mm256_add_epi16(__M, _mm256_and_si256(SeqIsN, __MatchMask));\
	const __m256i __MtoI  = _mm256_add_epi16(__M, __MI);\
	const __m256i __MtoD  = _mm256_add_epi16(__M, __MD);\
	const __m256i __MtoM  = _mm256_add_epi16(__M, __MM);\
	\
	__Insertions = _mm256_max_epi16(__ItoI, __MtoI);\
	__Deletions  = _mm256_max_epi16(__DtoD, __MtoD);\
	__Match      = _mm256_max_epi16(__Match, __MtoM);\
	\
	/* Store results in order I,D,M */\
	_mm256_store_si256((__m256i*) &CurrentColumn[MatrixOffset+iloop*48],    __Insertions);\
	_mm256_store_si256((__m256i*) &CurrentColumn[MatrixOffset+16+iloop*48], __Deletions);\
	_mm256_store_si256((__m256i*) &CurrentColumn[MatrixOffset+32+iloop*48], __Match);\
}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TwoDiagsLoop3
{
#define ASM_CODE(iloop, MatrixOffset) {\
	_mm256_store_si256((__m256i*) &CurrentColumn[MatrixOffset+iloop*48], __Insertions);\
	_mm256_store_si256((__m256i*) &CurrentColumn[MatrixOffset+iloop*48+16], __Deletions);\
	_mm256_store_si256((__m256i*) &CurrentColumn[MatrixOffset+iloop*48+32], __Match);\
	const __m256i KD = _mm256_or_si256(_mm256_cmpgt_epi16(__Deletions, __BestScore), _mm256_cmpeq_epi16(__Deletions, __BestScore));\
	__BestScore = _mm256_blendv_epi8(__BestScore, __Deletions, KD);\
	__BestAt    = _mm256_blendv_epi8(__BestAt, __At, KD);\
	__MorD      = _mm256_blendv_epi8(__MorD, __One, KD);\
	const __m256i KM = _mm256_or_si256(_mm256_cmpgt_epi16(__Match, __BestScore), _mm256_cmpeq_epi16(__Match, __BestScore));\
	__BestScore = _mm256_blendv_epi8(__BestScore, __Match, KM);\
	__BestAt    = _mm256_blendv_epi8(__BestAt, __At, KM);\
	__MorD      = _mm256_blendv_epi8(__MorD, __Zero, KM);\
}
#define TwoDiagsLoop3(iloop, MatrixOffset, PastSequence, MatchLoadRegister, MatchStoreRegister) \
{\
	/* Get and clean past column needs */\
	const __m256i __PreviousDeletion = _mm256_load_si256((__m256i*) (PreviousColumn+MatrixOffset+iloop*48+48+16));\
	const __m256i __PreviousMatch    = _mm256_load_si256((__m256i*) (PreviousColumn+MatrixOffset+iloop*48+32));\
	\
	/* Evict from L2 cache */\
	EVICTL1(PreviousColumn+MatrixOffset+iloop*48+48+16);\
	EVICTL1(PreviousColumn+MatrixOffset+iloop*48+32+16);\
	\
	/* Prefetch to L1 from L2*/\
	PREFETCHL1(PreviousColumn+MatrixOffset+(2+iloop)*48+48+16);\
	PREFETCHL1(PreviousColumn+MatrixOffset+(2+iloop)*48+32);\
	\
	/* Clean past states */\
	const __m256i __KOPI = _mm256_and_si256(__Insertions, __ClearMask);\
	MatchStoreRegister   = __Match;\
	\
	/* Get the previous diagonal for deletion, shift it */\
	__m256i __KOPD = _mm256_and_si256(_mm256_alignr_epi16(__Deletions, __PreviousDeletion, 15), __ClearMask);\
	\
	/* Coming from Deletion */\
	/*__m256i __DtoM;\
	__asm__ ( "vpaddd %2{{1to16}},%1,%0;\n        " : "=v"(__DtoM) : "v"(__KOPD), "m"(Globals[D_DM_INDEX]) );*/\
	const __m256i __DtoM = _mm256_add_epi16(__KOPD, _mm256_set1_epi16(Globals[D_DM_INDEX]));\
	const __m256i __DtoD = _mm256_add_epi16(__KOPD, __DD);\
	\
	/* Coming from Insertion */\
	/*const __m256i __ItoM = _mm256_add_epi16(__KOPI, __IM);*/\
	/*__m256i __ItoM;\
	__asm__ ( "vpaddd %2{{1to16}},%1,%0;\n        " : "=v"(__ItoM) : "v"(__KOPI), "m"(Globals[I_IM_INDEX]) );*/\
	const __m256i __ItoM = _mm256_add_epi16(__KOPI, _mm256_set1_epi16(Globals[I_IM_INDEX]));\
	const __m256i __ItoI = _mm256_add_epi16(__KOPI, __II);\
	\
	__Match = _mm256_max_epi16(__DtoM, __ItoM);\
	\
	/* Get the correct Genome string, for dummy NLOW is pefect */\
	__m256i __AddtoM, SeqIsN;\
	{\
		const __m256i __NewSeq = _mm256_alignr_epi16(PastSequence, __SequenceDNA, 15-iloop);\
		SeqIsN = _mm256_cmpeq_epi16(__NewSeq,*((__m256i*)__N));\
		__m256i DoWeMatch = _mm256_or_si256(_mm256_cmpeq_epi16(__NewSeq, __ReadSequence), ReadIsN);\
		__AddtoM  = _mm256_blendv_epi8(__MismatchScore, __MatchScore, DoWeMatch);\
	}\
	/* Coming from Match */\
	__m256i __M;\
	__M = _mm256_and_si256(_mm256_alignr_epi16(MatchLoadRegister, __PreviousMatch, 15), __ClearMask);\
	__M = _mm256_add_epi16(__M, _mm256_andnot_si256(SeqIsN, __AddtoM));\
	__M = _mm256_add_epi16(__M, _mm256_and_si256(SeqIsN, __MatchMask));\
	const __m256i __MtoI  = _mm256_add_epi16(__M, __MI);\
	const __m256i __MtoD  = _mm256_add_epi16(__M, __MD);\
	const __m256i __MtoM  = _mm256_add_epi16(__M, __MM);\
	\
	__Insertions = _mm256_max_epi16(__ItoI, __MtoI);\
	__Deletions  = _mm256_max_epi16(__DtoD, __MtoD);\
	__Match      = _mm256_max_epi16(__Match, __MtoM);\
	\
	/* Store results in order I,D,M */\
	ASM_CODE(iloop, MatrixOffset);\
	__At = _mm256_add_epi16(__At, __One);\
}
}

{
	assert(ReadLength >= 32);
	/*
	 * First 16 wide column
	 */
	{
		/* Get the read, check for N's*/

		const __m256i __ReadSequence = CONVERT_BROADCAST_UINT8(&ptrTag[0]);
		const __m256i ReadIsN        = _mm256_cmpeq_epi16(__ReadSequence, *((__m256i*)__N));
		
		__m256i __SequenceDNA = CONVERT_BROADCAST_UINT8(&genome[0]);
		__m256i SeqIsN        = _mm256_cmpeq_epi16(__SequenceDNA, *((__m256i*)__N));
		__m256i __AddtoM      = _mm256_andnot_si256(SeqIsN, __MismatchScore);
		__m256i DoWeMatch     = _mm256_or_si256(_mm256_cmpeq_epi16(__SequenceDNA, __ReadSequence), ReadIsN);
		__AddtoM              = _mm256_blendv_epi8( __AddtoM, __MatchScore, DoWeMatch);
		short int vtmp[16] __attribute__((aligned(32)));
		/* Prepare Match diagonal:
		 *  /15 14 ... 0/ = / 0 14 NLOW RowInit[MATCH] / 
		 */
		vtmp[0] = RowInit[MATCH] & CLEAR_MASK;
		for (int k=1; k<15; k++) vtmp[k] = (NLOW << SCORE_SHIFT);
		// Set a correct RowInit in the last position for alignr later on, this should not be a problem
		// for the initial use
		vtmp[15] = RowInit[MATCH] & CLEAR_MASK;
		const __m256i __KOPMInit = _mm256_load_si256((__m256i*) &vtmp[0]);

		/* Prepare Deletion diagonal:
		 * /15 14 ... 0/ = / 14 NLOW RowInit[DELETION]+DD  RowInit[DELETION] / 
		 */
		vtmp[0] = RowInit[DELETION] & CLEAR_MASK;
		vtmp[1] = (vtmp[0] + Globals[D_DD_INDEX]) & CLEAR_MASK;
		for (int k=2; k<16; k++) vtmp[k] = (NLOW << SCORE_SHIFT); 
		const  __m256i __KOPDInit = _mm256_load_si256((__m256i*) &vtmp[0]);
		
		/* Reversing Genome string */
		REVERSE_SEQUENCE(__SequenceDNA);

		/*  Compute first block of diags (16) */
		__m256i __KOPM, __KOPM2, __Insertions, __Match, __Deletions;
		const __m256i __DM = _mm256_set1_epi16(Globals[D_DM_INDEX]);
		const __m256i __IM = _mm256_set1_epi16(Globals[I_IM_INDEX]);
		{
			/************************************************ LINE 0 *********************************************/
			{
				/* Coming from Match */
				__m256i __M    = _mm256_add_epi16(__KOPMInit, __AddtoM);
				__m256i __MtoM = _mm256_add_epi16(__M, __MM);
				__m256i __MtoD = _mm256_add_epi16(__M, __MD);
				__Insertions   = _mm256_add_epi16(__M, __MI);
				
				/* Coming from Deletion */
				const __m256i __DtoM = _mm256_add_epi16(__KOPDInit, __DM);
				__Match              = _mm256_max_epi16(__MtoM, __DtoM);
				const __m256i __DtoD = _mm256_add_epi16(__KOPDInit, __DD);
				__Deletions          = _mm256_max_epi16(__MtoD, __DtoD);
				
				/* Coming from Insertion */ // Should not happen !!!
// 					const __m256i __ItoM  = _mm256_add_epi16(__KOPIInit, __IM);
// 					__Match               = _mm256_max_epi16(__Match, __ItoM);

				/* Store results in order I,D,M */
				_mm256_store_si256( (__m256i*) &IDM_Matrix[ 0], __Insertions);
				_mm256_store_si256( (__m256i*) &IDM_Matrix[16], __Deletions);
				_mm256_store_si256( (__m256i*) &IDM_Matrix[32], __Match);
			}
			
			/************************************************ LINE 1 *********************************************/
			{
				/* Clean past states */
				__m256i __KOPI  = _mm256_add_epi16(__Insertions, __ClearMask);
				__m256i __KOPD  = _mm256_add_epi16(__Deletions,  __ClearMask);
				__KOPM2         = _mm256_add_epi16(__Match,      __ClearMask);
				
				/* Get the previous diagonal for deletion, shift it */
				__KOPD = _mm256_alignr_epi16(__KOPD, __KOPDInit, 15);
				
				/* Coming from Deletion */
				const __m256i __DtoM = _mm256_add_epi16(__KOPD, __DM);
				const __m256i __DtoD = _mm256_add_epi16(__KOPD, __DD);
				
				/* Coming from Insertion */
				const __m256i __ItoM = _mm256_add_epi16(__KOPI, __IM);
				const __m256i __ItoI = _mm256_add_epi16(__KOPI, __II);
				
				__Match = _mm256_max_epi16(__DtoM, __ItoM);
				
				/* Get the correct Genome string, for dummy NLOW is pefect */
				{
					const __m256i __NewSeq = _mm256_alignr_epi16(__NLOW, __SequenceDNA, 14);
					SeqIsN                 = _mm256_cmpeq_epi16(__NewSeq,*((__m256i*)__N));
					DoWeMatch              = _mm256_or_si256(_mm256_cmpeq_epi16(__NewSeq, __ReadSequence), ReadIsN);
					__AddtoM               = _mm256_blendv_epi8(__MismatchScore, __MatchScore, DoWeMatch);
				}
				/* Coming from Match */
				__m256i __M;
// 				__M = _mm256_blendv_epi8(_mm256_and_si256(__KOPDInit, __ClearMask), __KOPMInit, 0/*FAUX*/);
				/*__M = _mm512_mask_add_epi32(__M, SeqIsNotN, __M, __AddtoM);*/
				__M = _mm256_add_epi16(__M, _mm256_andnot_si256(SeqIsN, __AddtoM));
				
				const __m256i __MtoI  = _mm256_add_epi16(__M, __MI);
				const __m256i __MtoD  = _mm256_add_epi16(__M, __MD);
				const __m256i __MtoM  = _mm256_add_epi16(__M, __MM);
				
				__Insertions = _mm256_max_epi16(__ItoI, __MtoI);
				__Deletions  = _mm256_max_epi16(__DtoD, __MtoD);
				__Match      = _mm256_max_epi16(__Match, __MtoM);
				
				/* Store results in order I,D,M */
				_mm256_store_si256( (__m256i*) &IDM_Matrix[48], __Insertions);
				_mm256_store_si256( (__m256i*) &IDM_Matrix[64], __Deletions);
				_mm256_store_si256( (__m256i*) &IDM_Matrix[80], __Match);
			}
			
			/************************************************ LINE 2-15 *********************************************/

			TwoDiagsLoop1( 2 , 0, __NLOW, __KOPM2, __KOPM);
			TwoDiagsLoop1( 3 , 0, __NLOW, __KOPM,  __KOPM2);
			TwoDiagsLoop1( 4 , 0, __NLOW, __KOPM2, __KOPM);
			TwoDiagsLoop1( 5 , 0, __NLOW, __KOPM,  __KOPM2);
			TwoDiagsLoop1( 6 , 0, __NLOW, __KOPM2, __KOPM);
			TwoDiagsLoop1( 7 , 0, __NLOW, __KOPM,  __KOPM2);
			TwoDiagsLoop1( 8 , 0, __NLOW, __KOPM2, __KOPM);
			TwoDiagsLoop1( 9 , 0, __NLOW, __KOPM,  __KOPM2);
			TwoDiagsLoop1(10 , 0, __NLOW, __KOPM2, __KOPM);
			TwoDiagsLoop1(11 , 0, __NLOW, __KOPM,  __KOPM2);
			TwoDiagsLoop1(12 , 0, __NLOW, __KOPM2, __KOPM);
			TwoDiagsLoop1(13 , 0, __NLOW, __KOPM,  __KOPM2);
			TwoDiagsLoop1(14 , 0, __NLOW, __KOPM2, __KOPM);
			TwoDiagsLoop1(15 , 0, __NLOW, __KOPM,  __KOPM2);
		}
		
		/* Compute the following blocks */
		__m256i __PastReversedSequenceDNA = __SequenceDNA;
		size_t iseq = 0;
		#pragma noprefetch
		while (iseq < (GenomeLength)) {
// 				printf("Block 16 : %lu\n", iseq);
			/* Get the sequence */
			__SequenceDNA = CONVERT_BROADCAST_UINT8(&genome[16+iseq]);
			/* Reversing Genome string */
			REVERSE_SEQUENCE(__SequenceDNA);
			const size_t itmp = (16+iseq)*48;
			
			TwoDiagsLoop1( 0 , itmp, __PastReversedSequenceDNA, __KOPM2, __KOPM);
			TwoDiagsLoop1( 1 , itmp, __PastReversedSequenceDNA, __KOPM,  __KOPM2);
			TwoDiagsLoop1( 2 , itmp, __PastReversedSequenceDNA, __KOPM2, __KOPM);
			TwoDiagsLoop1( 3 , itmp, __PastReversedSequenceDNA, __KOPM,  __KOPM2);
			TwoDiagsLoop1( 4 , itmp, __PastReversedSequenceDNA, __KOPM2, __KOPM);
			TwoDiagsLoop1( 5 , itmp, __PastReversedSequenceDNA, __KOPM,  __KOPM2);
			TwoDiagsLoop1( 6 , itmp, __PastReversedSequenceDNA, __KOPM2, __KOPM);
			TwoDiagsLoop1( 7 , itmp, __PastReversedSequenceDNA, __KOPM,  __KOPM2);
			TwoDiagsLoop1( 8 , itmp, __PastReversedSequenceDNA, __KOPM2, __KOPM);
			TwoDiagsLoop1( 9 , itmp, __PastReversedSequenceDNA, __KOPM,  __KOPM2);
			TwoDiagsLoop1(10 , itmp, __PastReversedSequenceDNA, __KOPM2, __KOPM);
			TwoDiagsLoop1(11 , itmp, __PastReversedSequenceDNA, __KOPM,  __KOPM2);
			TwoDiagsLoop1(12 , itmp, __PastReversedSequenceDNA, __KOPM2, __KOPM);
			TwoDiagsLoop1(13 , itmp, __PastReversedSequenceDNA, __KOPM,  __KOPM2);
			TwoDiagsLoop1(14 , itmp, __PastReversedSequenceDNA, __KOPM2, __KOPM);
			TwoDiagsLoop1(15 , itmp, __PastReversedSequenceDNA, __KOPM,  __KOPM2);
			
			__PastReversedSequenceDNA = __SequenceDNA;
			iseq += 16;
		}
	}
	/*
	 * Following columns
	 */
	const size_t lReadLength = ReadLength - 16;
	size_t iprf = 16;
	const short int * restrict PreviousColumn = &IDM_Matrix[48*14];
	short int * restrict CurrentColumn = &IDM_Matrix[48*MatrixLD];
	while ( iprf < lReadLength ) {
		/* Get the read, check for N's*/
		const __m256i __ReadSequence = CONVERT_BROADCAST_UINT8(&ptrTag[iprf]);
		const __m256i ReadIsN = _mm256_cmpeq_epi16(__ReadSequence, *((__m256i*)__N));
		__m256i __SequenceDNA = CONVERT_BROADCAST_UINT8(&genome[0]);
		
		short int vtmp[16] __attribute__((aligned(64)));
		/* Prepare Deletion diagonal:
		 * /15 14 ... 0/ = / 14 NLOW RowInit[DELETION]+DD  RowInit[DELETION] / 
		 */
		vtmp[0] = RowInit[DELETION] & CLEAR_MASK;
		vtmp[1] = (vtmp[0] + Globals[D_DD_INDEX]) & CLEAR_MASK;
		for (int k=2; k<16; k++) vtmp[k] = (NLOW << SCORE_SHIFT); 
		const __m256i __KOPDInit = _mm256_load_si256((__m256i*) &vtmp[0]);

		/*  Compute first block of diags (16) */
		__m256i __KOPM, __KOPM2, __Insertions, __Match, __Deletions;
		const __m256i __DM = _mm256_set1_epi16(Globals[D_DM_INDEX]);
		const __m256i __IM = _mm256_set1_epi16(Globals[I_IM_INDEX]);
		{
			/************************************************ LINE 0 *********************************************/
			{
				/* For computing first diagonal we need:
				 * 1: upper diagonal MATCH 
				 *      /15 14 ... 0/ = / 15 NLOW
				 *                        Element 15 MATCH of previous column row 14
				 *                      /
				 */
				__m256i __PreviousM1 = _mm256_and_si256(_mm256_load_si256((__m256i*) (PreviousColumn+32)), __ClearMask);
				EVICTL1(PreviousColumn+32);
				__PreviousM1         = _mm256_alignr_epi16(__NLOW, __PreviousM1, 15);
				/*
				 * 2: Diagonal shifted left DELETION 
				 *      /15 14 ... 0/ = / 14 NLOW  
				 *                        Element 15 DELETION of previous column row 14 + D + DD
				 *                        Element 15 DELETION of previous column row 15
				 *                      /
				 */
				__m256i __PreviousD1;
				{
					__m256i __PreviousD1_0 = _mm256_and_si256(_mm256_load_si256((__m256i*) (PreviousColumn+16))   , __ClearMask);
					__m256i __PreviousD1_1 = _mm256_and_si256(_mm256_load_si256((__m256i*) (PreviousColumn+16+48)), __ClearMask);
					
	// 					EVICTL1(PreviousColumn+16);
					EVICTL1(PreviousColumn+16+48);
					PREFETCHL1(PreviousColumn+96+48+16);
					PREFETCHL1(PreviousColumn+96+32);
					
					__PreviousD1_0         = _mm256_alignr_epi16(__NLOW, __PreviousD1_0, 14);
					__PreviousD1_1         = _mm256_alignr_epi16(__NLOW, __PreviousD1_1, 15);
					__PreviousD1_0         = _mm256_add_epi16(__PreviousD1_0, _mm256_and_si256(__DD, __ClearMask));
					__PreviousD1           = _mm256_blendv_epi8(__PreviousD1_1, __PreviousD1_0, OnlySecond.ymm);
				}
			
				const __m256i SeqIsN    = _mm256_cmpeq_epi16(__SequenceDNA, *((__m256i*)__N));
				const __m256i DoWeMatch = _mm256_or_si256(_mm256_cmpeq_epi16(__SequenceDNA, __ReadSequence), ReadIsN);
				const __m256i __AddtoM  = _mm256_blendv_epi8(__MismatchScore, __MatchScore, DoWeMatch);
		
				/* Coming from Match */
				__m256i __M    = _mm256_add_epi16(__PreviousM1, _mm256_andnot_si256(SeqIsN, __AddtoM));
				__m256i __MtoD = _mm256_add_epi16(__M, __MD);
				__Insertions   = _mm256_add_epi16(__M, __MI);
				__Match        = _mm256_add_epi16(__M, __MM);   
				
				/* Coming from Deletion */
				
				const __m256i __DtoM = _mm256_add_epi16(__PreviousD1, __DM);
				__Match              = _mm256_max_epi16(__Match, __DtoM);
				const __m256i __DtoD = _mm256_add_epi16(__PreviousD1, __DD);
				__Deletions          = _mm256_max_epi16(__MtoD, __DtoD);
				
				/* Coming from Insertion */ // Should not happen
// 				const __m256i __ItoM  = _mm256_add_epi16(__KOPIInit, __IM);
// 				__Match               = _mm256_max_epi16(__Match, __ItoM);
				
				/* Store results in order I,D,M */
				_mm256_store_si256((__m256i*) &CurrentColumn[ 0], __Insertions);
				_mm256_store_si256((__m256i*) &CurrentColumn[16], __Deletions);
				_mm256_store_si256((__m256i*) &CurrentColumn[32], __Match);
			}
			
			/* Reversing Genome string */
			REVERSE_SEQUENCE(__SequenceDNA);
			
			/************************************************ LINE 1 *********************************************/
			{
				/* For computing second diagonal we need:
				 * 1: upper diagonal MATCH 
				 *      /15 14 ... 0/ = / 14 NLOW
				 *                        Element 15 DELETION of previous column row 14 + D + DM
				 *                        Element 15 MATCH of previous column row 15
				 *                      /
				 */
				__m256i __PreviousM2;
				{
					__m256i __PreviousM2_0 = _mm256_and_si256(_mm256_load_si256((__m256i*) (PreviousColumn+16))   , __ClearMask);
					__m256i __PreviousM2_1 = _mm256_and_si256(_mm256_load_si256((__m256i*) (PreviousColumn+32+48)), __ClearMask);
					
					EVICTL1(PreviousColumn+16);
					EVICTL1(PreviousColumn+32+48);
					PREFETCHL1(PreviousColumn+144+48+16);
					PREFETCHL1(PreviousColumn+144+32);
					
					__PreviousM2_0 = _mm256_alignr_epi16(__NLOW, __PreviousM2_0, 14);
					__PreviousM2_1 = _mm256_alignr_epi16(__NLOW, __PreviousM2_1, 15);
					
					__PreviousM2_0 = _mm256_add_epi16(__PreviousM2_0, _mm256_and_si256(__DM, __ClearMask));
					__PreviousM2   = _mm256_blendv_epi8(__PreviousM2_1, __PreviousM2_0, OnlySecond.ymm);
				}
			
				/* Get and clean past column needs */
				const __m256i __PreviousDeletion = _mm256_load_si256((__m256i*) (PreviousColumn+96+16));
				EVICTL1(PreviousColumn+96+16);
// 					const __m256i __PreviousMatch    = _mm256_load_si256(PreviousColumn+48+32);
				
				/* Clean past states */
				__m256i __KOPI  = _mm256_and_si256(__Insertions, __ClearMask);
				__KOPM2         = __Match;
				
				/* Get the previous diagonal for deletion, shift it */
				__m256i __KOPD = _mm256_and_si256(_mm256_alignr_epi16(__Deletions, __PreviousDeletion, 15), __ClearMask);
				
				/* Coming from Deletion */
				const __m256i __DtoM = _mm256_add_epi16(__KOPD, __DM);
				const __m256i __DtoD = _mm256_add_epi16(__KOPD, __DD);
				
				/* Coming from Insertion */
				const __m256i __ItoM = _mm256_add_epi16(__KOPI, __IM);
				const __m256i __ItoI = _mm256_add_epi16(__KOPI, __II);
				
				__Match = _mm256_max_epi16(__DtoM, __ItoM);
				
				/* Get the correct Genome string, for dummy NLOW is pefect */
				__m256i __AddtoM, SeqIsN;
				{
					const __m256i __NewSeq  = _mm256_alignr_epi16(__NLOW, __SequenceDNA, 14);
					SeqIsN                  = _mm256_cmpeq_epi16(__NewSeq,*((__m256i*)__N));
					const __m256i DoWeMatch = _mm256_or_si256(_mm256_cmpeq_epi16(__NewSeq, __ReadSequence), ReadIsN);
					__AddtoM                = _mm256_blendv_epi8( __MismatchScore, __MatchScore, DoWeMatch);
				}
				/* Coming from Match */
				__m256i __M;
// 					__M = _mm256_alignr_epi16(__PreviousM2, __PreviousMatch, 15);
				/*__M = _mm512_mask_add_epi32(__PreviousM2, SeqIsNotN, __PreviousM2, __AddtoM);*/
				__M = _mm256_add_epi16(__PreviousM2, _mm256_andnot_si256(SeqIsN, __AddtoM));
				const __m256i __MtoI  = _mm256_add_epi16(__M, __MI);
				const __m256i __MtoD  = _mm256_add_epi16(__M, __MD);
				const __m256i __MtoM  = _mm256_add_epi16(__M, __MM);
				
				__Insertions = _mm256_max_epi16(__ItoI, __MtoI);
				__Deletions  = _mm256_max_epi16(__DtoD, __MtoD);
				__Match      = _mm256_max_epi16(__Match, __MtoM);
				
				/* Store results in order I,D,M */
				_mm256_store_si256( (__m256i*) &CurrentColumn[48], __Insertions);
				_mm256_store_si256( (__m256i*) &CurrentColumn[64], __Deletions);
				_mm256_store_si256( (__m256i*) &CurrentColumn[80], __Match);
			}
			
			/************************************************ LINE 2-15 *********************************************/

			TwoDiagsLoop2( 2, 0, __NLOW, __KOPM2, __KOPM);
			TwoDiagsLoop2( 3, 0, __NLOW, __KOPM,  __KOPM2);
			TwoDiagsLoop2( 4, 0, __NLOW, __KOPM2, __KOPM);
			TwoDiagsLoop2( 5, 0, __NLOW, __KOPM,  __KOPM2);
			TwoDiagsLoop2( 6, 0, __NLOW, __KOPM2, __KOPM);
			TwoDiagsLoop2( 7, 0, __NLOW, __KOPM,  __KOPM2);
			TwoDiagsLoop2( 8, 0, __NLOW, __KOPM2, __KOPM);
			TwoDiagsLoop2( 9, 0, __NLOW, __KOPM,  __KOPM2);
			TwoDiagsLoop2(10, 0, __NLOW, __KOPM2, __KOPM);
			TwoDiagsLoop2(11, 0, __NLOW, __KOPM,  __KOPM2);
			TwoDiagsLoop2(12, 0, __NLOW, __KOPM2, __KOPM);
			TwoDiagsLoop2(13, 0, __NLOW, __KOPM,  __KOPM2);
			TwoDiagsLoop2(14, 0, __NLOW, __KOPM2, __KOPM);
			TwoDiagsLoop2(15, 0, __NLOW, __KOPM,  __KOPM2);
		}
		
		/* Compute the following blocks */
		__m256i __PastReversedSequenceDNA = __SequenceDNA;
		size_t iseq = 0;

		while ( iseq < (GenomeLength) ) {
// 				printf("Block 16 : %lu\n", iseq);
			/* Get the sequence */
			__SequenceDNA = CONVERT_BROADCAST_UINT8(&genome[16+iseq]);
			/* Reversing Genome string */
			REVERSE_SEQUENCE(__SequenceDNA);
			const size_t itmp = 48*(iseq+16);
			TwoDiagsLoop2( 0 , itmp, __PastReversedSequenceDNA, __KOPM2, __KOPM);
			TwoDiagsLoop2( 1 , itmp, __PastReversedSequenceDNA, __KOPM,  __KOPM2);
			TwoDiagsLoop2( 2 , itmp, __PastReversedSequenceDNA, __KOPM2, __KOPM);
			TwoDiagsLoop2( 3 , itmp, __PastReversedSequenceDNA, __KOPM,  __KOPM2);
			TwoDiagsLoop2( 4 , itmp, __PastReversedSequenceDNA, __KOPM2, __KOPM);
			TwoDiagsLoop2( 5 , itmp, __PastReversedSequenceDNA, __KOPM,  __KOPM2);
			TwoDiagsLoop2( 6 , itmp, __PastReversedSequenceDNA, __KOPM2, __KOPM);
			TwoDiagsLoop2( 7 , itmp, __PastReversedSequenceDNA, __KOPM,  __KOPM2);
			TwoDiagsLoop2( 8 , itmp, __PastReversedSequenceDNA, __KOPM2, __KOPM);
			TwoDiagsLoop2( 9 , itmp, __PastReversedSequenceDNA, __KOPM,  __KOPM2);
			TwoDiagsLoop2(10 , itmp, __PastReversedSequenceDNA, __KOPM2, __KOPM);
			TwoDiagsLoop2(11 , itmp, __PastReversedSequenceDNA, __KOPM,  __KOPM2);
			TwoDiagsLoop2(12 , itmp, __PastReversedSequenceDNA, __KOPM2, __KOPM);
			TwoDiagsLoop2(13 , itmp, __PastReversedSequenceDNA, __KOPM,  __KOPM2);
			TwoDiagsLoop2(14 , itmp, __PastReversedSequenceDNA, __KOPM2, __KOPM);
			TwoDiagsLoop2(15 , itmp, __PastReversedSequenceDNA, __KOPM,  __KOPM2);
			
			__PastReversedSequenceDNA = __SequenceDNA;
			iseq += 16;
		}
		
		iprf += 16;
		CurrentColumn  += 48*MatrixLD;
		PreviousColumn += 48*MatrixLD;
	}
	
	/* Last Column */
	{
		/* Get the read, check for N's*/
		const __m256i __ReadSequence   = CONVERT_BROADCAST_UINT8(&ptrTag[iprf]);
		const __m256i ReadIsN          = _mm256_cmpeq_epi16(__ReadSequence, *((__m256i*)__N));
		register __m256i __SequenceDNA = CONVERT_BROADCAST_UINT8(&genome[0]);
		/*  Compute first block of diags (16) */
		register __m256i __KOPM, __KOPM2, __Insertions, __Match, __Deletions;
		{
			/************************************************ LINE 0 *********************************************/
			const __m256i __DM = _mm256_set1_epi16(Globals[D_DM_INDEX]);
			const __m256i __IM = _mm256_set1_epi16(Globals[I_IM_INDEX]);
			{
				/* For computing first diagonal we need:
				 * 1: upper diagonal MATCH 
				 *      /15 14 ... 0/ = / 15 NLOW
				 *                        Element 15 MATCH of previous column row 14
				 *                      /
				 */
				__m256i __PreviousM1 = _mm256_and_si256(_mm256_load_si256( (__m256i*) PreviousColumn+32), __ClearMask);
				EVICTL1(PreviousColumn+32);
				__PreviousM1         = _mm256_alignr_epi16(__NLOW, __PreviousM1, 15);
				/*
				 * 2: Diagonal shifted left DELETION 
				 *      /15 14 ... 0/ = / 14 NLOW  
				 *                        Element 15 DELETION of previous column row 14 + D + DD
				 *                        Element 15 DELETION of previous column row 15
				 *                      /
				 */
				__m256i __PreviousD1;
				{
					__m256i __PreviousD1_0 = _mm256_and_si256(_mm256_load_si256((__m256i*) PreviousColumn+16)   , __ClearMask);
					__m256i __PreviousD1_1 = _mm256_and_si256(_mm256_load_si256((__m256i*) PreviousColumn+16+48), __ClearMask);
					
	// 					EVICTL1(PreviousColumn+16);
					EVICTL1(PreviousColumn+16+48);
					PREFETCHL1(PreviousColumn+96+48+16);
					PREFETCHL1(PreviousColumn+96+32);
					
					__PreviousD1_0 = _mm256_alignr_epi16(__NLOW, __PreviousD1_0, 14);
					__PreviousD1_1 = _mm256_alignr_epi16(__NLOW, __PreviousD1_1, 15);
					__PreviousD1_0 = _mm256_add_epi16(__PreviousD1_0, _mm256_and_si256(__DD, __ClearMask));
					__PreviousD1   = _mm256_blendv_epi8(__PreviousD1_1, __PreviousD1_0, OnlySecond.ymm);
				}
				
				const __m256i SeqIsN    = _mm256_cmpeq_epi16(__SequenceDNA, *((__m256i*)__N));
				const __m256i DoWeMatch = _mm256_or_si256(_mm256_cmpeq_epi16(__SequenceDNA, __ReadSequence), ReadIsN);
				const __m256i __AddtoM  = _mm256_blendv_epi8(__MismatchScore, __MatchScore, DoWeMatch);
				
				/* Coming from Match */
				__m256i __M          = _mm256_add_epi16(__PreviousM1, _mm256_andnot_si256(SeqIsN, __AddtoM));
				__m256i __MtoD       = _mm256_add_epi16(__M, __MD);
				__Insertions         = _mm256_add_epi16(__M, __MI);
				__Match              = _mm256_add_epi16(__M, __MM);
				
				/* Coming from Deletion */
				
				const __m256i __DtoM = _mm256_add_epi16(__PreviousD1, __DM);
				__Match              = _mm256_max_epi16(__Match, __DtoM);
				const __m256i __DtoD = _mm256_add_epi16(__PreviousD1, __DD);
				__Deletions          = _mm256_max_epi16(__MtoD, __DtoD);
				
				/* Coming from Insertion */ // Should not happen !!!
// 					const __m256i __ItoM  = _mm256_add_epi16(__KOPIInit, __IM);
// 					__Match               = _mm256_max_epi16(__Match, __ItoM);
				
				/* Store results in order I,D,M */
				_mm256_store_si256((__m256i*) &CurrentColumn[ 0], __Insertions);
				_mm256_store_si256((__m256i*) &CurrentColumn[16], __Deletions);
				_mm256_store_si256((__m256i*) &CurrentColumn[32], __Match);
			}
			
			/* Reversing Genome string */
			REVERSE_SEQUENCE(__SequenceDNA);
			
			/************************************************ LINE 1 *********************************************/
			{
				/* For computing second diagonal we need:
				 * 1: upper diagonal MATCH 
				 *      /15 14 ... 0/ = / 14 NLOW
				 *                        Element 15 DELETION of previous column row 14 + D + DM
				 *                        Element 15 MATCH of previous column row 15
				 *                      /
				 */
				__m256i __PreviousM2;
				{
					__m256i __PreviousM2_0 = _mm256_and_si256(_mm256_load_si256((__m256i*) (PreviousColumn+16))   , __ClearMask);
					__m256i __PreviousM2_1 = _mm256_and_si256(_mm256_load_si256((__m256i*) (PreviousColumn+32+48)), __ClearMask);
					
					EVICTL1(PreviousColumn+16);
					EVICTL1(PreviousColumn+32+48);
					PREFETCHL1(PreviousColumn+144+48+16);
					PREFETCHL1(PreviousColumn+144+32);
					
					__PreviousM2_0         = _mm256_alignr_epi16(__NLOW, __PreviousM2_0, 14);
					__PreviousM2_1         = _mm256_alignr_epi16(__NLOW, __PreviousM2_1, 15);
					
					__PreviousM2_0         = _mm256_add_epi16(__PreviousM2_0, _mm256_and_si256(__DM, __ClearMask));
					__PreviousM2           = _mm256_blendv_epi8(__PreviousM2_1, __PreviousM2_0, OnlySecond.ymm);
				}
				/* Get and clean past column needs */
				const __m256i __PreviousDeletion = _mm256_load_si256((__m256i*) (PreviousColumn+96+16));
				EVICTL1(PreviousColumn+96+16);
// 					const __m256i __PreviousMatch    = _mm256_load_si256(PreviousColumn+48+32);
				
				/* Clean past states */
				__m256i __KOPI  = _mm256_and_si256(__Insertions, __ClearMask);
				__KOPM2         = __Match;
				
				/* Get the previous diagonal for deletion, shift it */
				__m256i __KOPD = _mm256_and_si256(_mm256_alignr_epi16(__Deletions, __PreviousDeletion, 15), __ClearMask);
				
				/* Coming from Deletion */
				const __m256i __DtoM = _mm256_add_epi16(__KOPD, __DM);
				const __m256i __DtoD = _mm256_add_epi16(__KOPD, __DD);
				
				/* Coming from Insertion */
				const __m256i __ItoM = _mm256_add_epi16(__KOPI, __IM);
				const __m256i __ItoI = _mm256_add_epi16(__KOPI, __II);
				
				__Match = _mm256_max_epi16(__DtoM, __ItoM);
				
				/* Get the correct Genome string, for dummy NLOW is pefect */
				__m256i __AddtoM, SeqIsN; 
				{
					const __m256i __NewSeq  = _mm256_alignr_epi16(__NLOW, __SequenceDNA, 14);
					SeqIsN                  = _mm256_cmpeq_epi16(__NewSeq,*((__m256i*)__N));
					const __m256i DoWeMatch = _mm256_or_si256(_mm256_cmpeq_epi16(__NewSeq, __ReadSequence), ReadIsN);
					__AddtoM                = _mm256_blendv_epi8(__MismatchScore, __MatchScore, DoWeMatch);
				}
				/* Coming from Match */
				__m256i __M;
// 					__M = _mm256_alignr_epi16(__PreviousM2, __PreviousMatch, 15);
				/*__M = _mm512_mask_add_epi32(__PreviousM2, SeqIsNotN, __PreviousM2, __AddtoM);*/
				__M = _mm256_add_epi16(__PreviousM2, _mm256_andnot_si256(SeqIsN, __AddtoM));
				const __m256i __MtoI  = _mm256_add_epi16(__M, __MI);
				const __m256i __MtoD  = _mm256_add_epi16(__M, __MD);
				const __m256i __MtoM  = _mm256_add_epi16(__M, __MM);
				
				__Insertions = _mm256_max_epi16(__ItoI, __MtoI);
				__Deletions  = _mm256_max_epi16(__DtoD, __MtoD);
				__Match      = _mm256_max_epi16(__Match, __MtoM);
				
				/* Store results in order I,D,M */
				_mm256_store_si256((__m256i*) &CurrentColumn[48], __Insertions);
				_mm256_store_si256((__m256i*) &CurrentColumn[64], __Deletions);
				_mm256_store_si256((__m256i*) &CurrentColumn[80], __Match);
			}
			
			/************************************************ LINE 2-15 *********************************************/

			TwoDiagsLoop2( 2 , 0, __NLOW, __KOPM2, __KOPM);
			TwoDiagsLoop2( 3 , 0, __NLOW, __KOPM,  __KOPM2);
			TwoDiagsLoop2( 4 , 0, __NLOW, __KOPM2, __KOPM);
			TwoDiagsLoop2( 5 , 0, __NLOW, __KOPM,  __KOPM2);
			TwoDiagsLoop2( 6 , 0, __NLOW, __KOPM2, __KOPM);
			TwoDiagsLoop2( 7 , 0, __NLOW, __KOPM,  __KOPM2);
			TwoDiagsLoop2( 8 , 0, __NLOW, __KOPM2, __KOPM);
			TwoDiagsLoop2( 9 , 0, __NLOW, __KOPM,  __KOPM2);
			TwoDiagsLoop2(10 , 0, __NLOW, __KOPM2, __KOPM);
			TwoDiagsLoop2(11 , 0, __NLOW, __KOPM,  __KOPM2);
			TwoDiagsLoop2(12 , 0, __NLOW, __KOPM2, __KOPM);
			TwoDiagsLoop2(13 , 0, __NLOW, __KOPM,  __KOPM2);
			TwoDiagsLoop2(14 , 0, __NLOW, __KOPM2, __KOPM);
			TwoDiagsLoop2(15 , 0, __NLOW, __KOPM,  __KOPM2);
		}
		
		/* Compute the following blocks */
		register __m256i __PastReversedSequenceDNA = __SequenceDNA;
		size_t iseq = 0;
		static int One = 1;
		const __m256i __Zero = _mm256_setzero_si256();
		__m256i __BestScore  = _mm256_setzero_si256();
		__m256i __BestAt     = _mm256_setzero_si256();
		__m256i __MorD       = _mm256_setzero_si256();
		__m256i __At         = _mm256_load_si256((__m256i*) &At[0]);
		const __m256i __One  = _mm256_set1_epi16 (One);

		while ( iseq < (GenomeLength) ) {
			/* Get the sequence */
			__SequenceDNA = CONVERT_BROADCAST_UINT8(&genome[iseq+16]);

			/* Reversing Genome string */
			REVERSE_SEQUENCE(__SequenceDNA);
			const size_t itmp = (16+iseq)*48;
			TwoDiagsLoop3( 0 , itmp, __PastReversedSequenceDNA, __KOPM2, __KOPM);
			TwoDiagsLoop3( 1 , itmp, __PastReversedSequenceDNA, __KOPM,  __KOPM2);
			TwoDiagsLoop3( 2 , itmp, __PastReversedSequenceDNA, __KOPM2, __KOPM);
			TwoDiagsLoop3( 3 , itmp, __PastReversedSequenceDNA, __KOPM,  __KOPM2);
			TwoDiagsLoop3( 4 , itmp, __PastReversedSequenceDNA, __KOPM2, __KOPM);
			TwoDiagsLoop3( 5 , itmp, __PastReversedSequenceDNA, __KOPM,  __KOPM2);
			TwoDiagsLoop3( 6 , itmp, __PastReversedSequenceDNA, __KOPM2, __KOPM);
			TwoDiagsLoop3( 7 , itmp, __PastReversedSequenceDNA, __KOPM,  __KOPM2);
			TwoDiagsLoop3( 8 , itmp, __PastReversedSequenceDNA, __KOPM2, __KOPM);
			TwoDiagsLoop3( 9 , itmp, __PastReversedSequenceDNA, __KOPM,  __KOPM2);
			TwoDiagsLoop3(10 , itmp, __PastReversedSequenceDNA, __KOPM2, __KOPM);
			TwoDiagsLoop3(11 , itmp, __PastReversedSequenceDNA, __KOPM,  __KOPM2);
			TwoDiagsLoop3(12 , itmp, __PastReversedSequenceDNA, __KOPM2, __KOPM);
			TwoDiagsLoop3(13 , itmp, __PastReversedSequenceDNA, __KOPM,  __KOPM2);
			TwoDiagsLoop3(14 , itmp, __PastReversedSequenceDNA, __KOPM2, __KOPM);
			TwoDiagsLoop3(15 , itmp, __PastReversedSequenceDNA, __KOPM,  __KOPM2);
			
			__PastReversedSequenceDNA = __SequenceDNA;
			iseq += 16;
		}
		_mm256_store_si256((__m256i*)results, __BestScore);
		_mm256_store_si256((__m256i*)Id, __BestAt);
		_mm256_store_si256((__m256i*)MorD, __MorD);
	}
	/* Memory fence */
	// no non temporal streaming
}
#undef TwoDiagsLoop1
#undef TwoDiagsLoop2
#undef TwoDiagsLoop3
#undef ASM_CODE
#undef CONVERT_BROADCAST_UINT8
#endif

//
// -------------------------------------------------------
// GET BEST ALIGNMENT 
// -------------------------------------------------------
//
#ifdef __GET_ALIGNMENT__
{
	int iprf = ReadLength-1;
	enum StateTypeOffset StateOffset = (EndOnMatchOrDeletion & 0x1) ? DeletionOffset : MatchOffset;	
	unsigned int previous_State = (EndOnMatchOrDeletion & 0x1) ? PRIORITY_DELETION : PRIORITY_MATCH;
	unsigned int StateCount = 0;
	unsigned int current_State;
	Align[1] = iseq;
	
	while ( iprf >= 0 && iseq >= 0)
	{
		unsigned short int * ptr        = Location(IDM_Matrix, MatrixLD, iprf, iseq, StateOffset);
		const unsigned short int lScore = *(ptr);
		current_State                   = (lScore & STATE_MASK) >> STATE_SHIFT;
		
		if (current_State != previous_State) {
			ptrCigar -= 2;
			if (ptrCigar < cigarStorageBottomLimit) goto OUT;
			ptrCigar[0] = (unsigned char) StateCount;
			ptrCigar[1] = StateLetter[previous_State];
			previous_State = current_State;
			StateCount = 1;
			StateChange += StateChangeAdd[current_State];
		}
		else {
			StateCount++;
		}
		
		switch (current_State)
		{
			case PRIORITY_MATCH:
// 				fprintf(stderr, "%i %i %i %c\n", iseq, iprf, ToScore(lScore), 'M');
				
				StateOffset = MatchOffset;
				if (lScore & MISMATCH_MASK){
					ptrMismatch -= 2;
					if (ptrMismatch <  mismatchStorageBottomLimit) goto OUT;
					ptrMismatch[0] = (unsigned char) iprf;
#ifdef __WITH_SOFTCLIP_EXTENSION__
					ptrMismatch[1] = ptrTag[iprf] | 0x80;
#else
					ptrMismatch[1] = ptrTag[iprf];
#endif
					MismatchCount++;
// 					fputc('m', stderr);
				}
				else {
					if (genome[iseq] == 'N' && ptrTag[iprf] != 'N') {
						ptrMismatch -= 2;
						if (ptrMismatch <  mismatchStorageBottomLimit) goto OUT;
						ptrMismatch[0] = (unsigned char) iprf;
						ptrMismatch[1] = ptrTag[iprf];
					}
					else if (ptrTag[iprf] == 'N') {
						ptrMismatch -= 2;
						if (ptrMismatch <  mismatchStorageBottomLimit) goto OUT;
						ptrMismatch[0] = (unsigned char) iprf;
						ptrMismatch[1] = 'N';
					}
// 					fputc('M', stderr);
				}
				--iseq;
				--iprf;
				break;
			case PRIORITY_INSERTION:
// 				fprintf(stderr, "%i %i %i %c\n", iseq, iprf, ToScore(lScore), 'I');
				--iseq;
				StateOffset = InsertionOffset;
// 				fputc('D', stderr);
				break;
			case PRIORITY_DELETION:
// 				fprintf(stderr, "%i %i %i %c\n", iseq, iprf, ToScore(lScore), 'D');
				StateOffset = DeletionOffset;
				ptrMismatch -= 2;
				if (ptrMismatch <  mismatchStorageBottomLimit) goto OUT;
				ptrMismatch[0] = (unsigned char) iprf;
#ifdef __WITH_SOFTCLIP_EXTENSION__
				ptrMismatch[1] = ptrTag[iprf] | 0x80;
#else				
				ptrMismatch[1] = ptrTag[iprf];
#endif
				--iprf;
				InsertionCount++;
// 				fputc('I', stderr);
				break;
			default:
				fprintf(stderr,"\nUnknown state (%u) encountered from cell (%i,%i)\n", current_State, iseq, iprf);
				--iprf;
				ret = -1;
				goto OUT;
				break;
		}
#ifdef _TEST_
		++counter;
#endif
	}
	
	ptrCigar -= 2;
	if (ptrCigar < cigarStorageBottomLimit) goto OUT;
	ptrCigar[0] = (unsigned char) StateCount;
	ptrCigar[1] = StateLetter[previous_State];
	StateChange += StateChangeAdd[previous_State];
	
	 
	if (iprf >= 0 ) {
		ptrCigar -= 2;
		if (ptrCigar < cigarStorageBottomLimit) goto OUT;
		ptrCigar[0] = (unsigned char) 1+iprf;
		ptrCigar[1] = 'I';
		StateChange += 1; // Adding last deletion state coming from match, that is 1 here.
		ptrMismatch -= 2*(1+iprf);
		if (ptrMismatch <  mismatchStorageBottomLimit) {
			ret = -1;
			goto OUT;
		}
		else {
			int i=0;
			while (i <= iprf) {
				ptrMismatch[2*i]   = (unsigned char) i;
				ptrMismatch[2*i+1] = ptrTag[i]; 
				i++;
			}
			InsertionCount += 1+iprf;
		}
	}
	Align[0] = iseq < 0 ? 0 : (current_State == PRIORITY_MATCH) ? iseq+1 : iseq+1;
	ret = 0;
	OUT:;
// 	fputc('\n', stderr);
// 	printf("CUrrent state: %d iseq %i\n", (int) current_State, iseq);
}
#endif

#ifdef __GET_ALIGNMENT_STATE__
{
	enum StateTypeOffset StateOffset = (EndOnMatchOrDeletion & 0x1) ? DeletionOffset : MatchOffset;
	unsigned int previous_State = (EndOnMatchOrDeletion & 0x1) ? PRIORITY_DELETION : PRIORITY_MATCH;
	unsigned int current_State;
	Align[1] = iseq;
	ret = -1;

	while ( iprf >= 0 && iseq >= 0)
	{
		if (--ptrState < StateLimit) goto OUT;
		
		unsigned short int* ptr         = Location(IDM_Matrix, MatrixLD, iprf, iseq, StateOffset);
		const unsigned short int lScore = *(ptr);
		current_State                   = (lScore & STATE_MASK) >> STATE_SHIFT;
		if (current_State != previous_State) {
			RequiredCigSize++;
			previous_State = current_State;
		}
		
		switch(current_State) {
			case PRIORITY_MATCH:
				StateOffset = MatchOffset;
				if (lScore & MISMATCH_MASK) {
					*ptrState = 'm';
					MismatchCount++;
					RequiredMMSize++;
				}
				else {
					if (genome[iseq] == 'N' && ptrTag[iprf] != 'N') {
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
				--iseq;
				--iprf;
				break;
			case PRIORITY_INSERTION:
				--iseq;
				*ptrState = 'D';
				StateOffset = InsertionOffset;
				ContainsInDel = 1;
				break;
			case PRIORITY_DELETION:
				StateOffset = DeletionOffset;
				*ptrState = 'I';
				--iprf;
				ContainsInDel = 1;
// 				InsertionCount++;
				RequiredMMSize++;
				break;
			default:
				fprintf(stderr,"\nUnknown state (%u) encountered from cell (%i,%i)\n", current_State, iseq, iprf);
				--iprf;
				goto OUT;
				break;
		}
		
		if (iprf == SoftClipBoundary) {
			SoftClipMismatchCount  = MismatchCount;
			SoftClipEvictedMMSize  = RequiredMMSize;
			SoftClipEvictedCigSize = RequiredCigSize;
		}
		
#ifdef _TEST_
		++counter;
#endif
	}
	
	if (iprf >= 0) {
		ptrState -= (1+iprf);
		if (ptrState < StateLimit) {
			goto OUT;
		}
		else {
			int i=0;
			while (i <= iprf) {
				ptrState[i] = 'I';
				i++;
			}
			
			if (iprf >= SoftClipBoundary) {
				const int count = iprf - SoftClipBoundary;
				SoftClipEvictedMMSize  = RequiredMMSize + count;
				SoftClipEvictedCigSize = RequiredCigSize + (current_State == PRIORITY_DELETION) ? 0 : 1;
			}
			
// 			InsertionCount += 1+iprf;
			RequiredMMSize  += 1+iprf;
			RequiredCigSize += (current_State == PRIORITY_DELETION) ? 0 : 1; 
			ContainsInDel = 1;
		}
	}
	Align[0] = iseq < 0 ? 0 : iseq+1;
	ret = 0;
	OUT:;
}
#endif


//
// -------------------------------------------------------
// FREE GLOBAL MEMORY
// -------------------------------------------------------
//
#ifdef __FREE_MEMORY__
free_huge_pages(IDM_Matrix);
#endif

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
