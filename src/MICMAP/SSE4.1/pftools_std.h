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
#  if defined(__EXTRA_FUNCTIONS__) || defined(__ALLOCATE_MEMORY__) || defined(__FILL_MATRIX__) || defined(__GET_ALIGNMENT__) || defined(__FREE_MEMORY__) || defined(__DECLARE_DEFINITIONS)
#    error "Too many preprocessing variables defined"
#  endif
#elif defined(__EXTRA_FUNCTIONS__)
#  if defined(__DECLARE_VARIABLES__) || defined(__ALLOCATE_MEMORY__) || defined(__FILL_MATRIX__) || defined(__GET_ALIGNMENT__) || defined(__FREE_MEMORY__) || defined(__DECLARE_DEFINITIONS)
#    error "Too many preprocessing variables defined"
#  endif
#elif defined(__ALLOCATE_MEMORY__)
#  if defined(__EXTRA_FUNCTIONS__) || defined(__DECLARE_VARIABLES__) || defined(__FILL_MATRIX__) || defined(__GET_ALIGNMENT__) || defined(__FREE_MEMORY__) || defined(__DECLARE_DEFINITIONS)
#    error "Too many preprocessing variables defined"
#  endif
#elif defined(__FILL_MATRIX__)
#  if defined(__EXTRA_FUNCTIONS__) || defined(__ALLOCATE_MEMORY__) || defined(__DECLARE_VARIABLES__) || defined(__GET_ALIGNMENT__) || defined(__FREE_MEMORY__) || defined(__DECLARE_DEFINITIONS)
#    error "Too many preprocessing variables defined"
#  endif
#elif defined(__GET_ALIGNMENT__)
#  if defined(__EXTRA_FUNCTIONS__) || defined(__ALLOCATE_MEMORY__) || defined(__FILL_MATRIX__) || defined(__DECLARE_VARIABLES__) || defined(__FREE_MEMORY__) || defined(__DECLARE_DEFINITIONS)
#    error "Too many preprocessing variables defined"
#  endif
#elif defined(__FREE_MEMORY__)
#  if defined(__EXTRA_FUNCTIONS__) || defined(__ALLOCATE_MEMORY__) || defined(__FILL_MATRIX__) || defined(__GET_ALIGNMENT__) || defined(__DECLARE_VARIABLES__) || defined(__DECLARE_DEFINITIONS)
#    error "Too many preprocessing variables defined"
#  endif
#elif defined(__DECLARE_DEFINITIONS__)
#  if defined(__EXTRA_FUNCTIONS__) || defined(__ALLOCATE_MEMORY__) || defined(__FILL_MATRIX__) || defined(__GET_ALIGNMENT__) || defined(__DECLARE_VARIABLES__) || defined(__FREE_MEMORY__)
#    error "Too many preprocessing variables defined"
#  endif
#endif

#ifdef __DECLARE_DEFINITIONS__
	typedef union sIOP {
		struct {
			int Match;
			int Insertion;
		} Element;
		__m64 mm;
	} sIOP;
#endif
	
/*
 * -------------------------------------------------------
 * DECLARE VARIABLE TO BE GLOBALLY USED WITH ALIGNMENT   
 * -------------------------------------------------------
 */

#ifdef __DECLARE_STATIC_VARIABLES__
/*                  D                 X                 M                    I (see vector position) */
	static const union { int i[4]; __m128i xmm; } MatchScore = {
		.i = {
			((_M+_MD) << SCORE_SHIFT) + (PRIORITY_MATCH << STATE_SHIFT) + (0 << MISMATCH_SHIFT),
			((_M+_MX) << SCORE_SHIFT) + (PRIORITY_MATCH << STATE_SHIFT) + (0 << MISMATCH_SHIFT),
			((_M+_MM) << SCORE_SHIFT) + (PRIORITY_MATCH << STATE_SHIFT) + (0 << MISMATCH_SHIFT),
			((_M+_MI) << SCORE_SHIFT) + (PRIORITY_MATCH << STATE_SHIFT) + (0 << MISMATCH_SHIFT)
		}
	};
	static const union { int i[4]; __m128i xmm; } MismatchScore = {
		.i = {
			((_m+_MD) << SCORE_SHIFT) + (PRIORITY_MATCH << STATE_SHIFT) + (1 << MISMATCH_SHIFT),
			((_m+_MX) << SCORE_SHIFT) + (PRIORITY_MATCH << STATE_SHIFT) + (1 << MISMATCH_SHIFT),
			((_m+_MM) << SCORE_SHIFT) + (PRIORITY_MATCH << STATE_SHIFT) + (1 << MISMATCH_SHIFT),
			((_m+_MI) << SCORE_SHIFT) + (PRIORITY_MATCH << STATE_SHIFT) + (1 << MISMATCH_SHIFT)
		}
	};
	static const union { int i[4]; __m128i xmm; } SequenceNMatchScore = {
		.i = {
			((_MD) << SCORE_SHIFT) + (PRIORITY_MATCH << STATE_SHIFT) + (0 << MISMATCH_SHIFT),
			((_MX) << SCORE_SHIFT) + (PRIORITY_MATCH << STATE_SHIFT) + (0 << MISMATCH_SHIFT),
			((_MM) << SCORE_SHIFT) + (PRIORITY_MATCH << STATE_SHIFT) + (0 << MISMATCH_SHIFT),
			((_MI) << SCORE_SHIFT) + (PRIORITY_MATCH << STATE_SHIFT) + (0 << MISMATCH_SHIFT)
		}
	};
	
	static const union { int i[4]; __m128i xmm; } InsertionScore = {
		.i = {
			((NLOW  ) << SCORE_SHIFT) + (PRIORITY_INSERTION << STATE_SHIFT) + (1 << MISMATCH_SHIFT),
			((_I+_IX) << SCORE_SHIFT) + (PRIORITY_INSERTION << STATE_SHIFT) + (1 << MISMATCH_SHIFT),
			((_I+_IM) << SCORE_SHIFT) + (PRIORITY_INSERTION << STATE_SHIFT) + (1 << MISMATCH_SHIFT),
			((_I+_II) << SCORE_SHIFT) + (PRIORITY_INSERTION << STATE_SHIFT) + (1 << MISMATCH_SHIFT)
		}
	};
	static const union { int i[4]; __m128i xmm; } DeletionScore = {
		.i = {
			((_D+_DD) << SCORE_SHIFT) + (PRIORITY_DELETION << STATE_SHIFT) + (1 << MISMATCH_SHIFT),
			((_D+_DX) << SCORE_SHIFT) + (PRIORITY_DELETION << STATE_SHIFT) + (1 << MISMATCH_SHIFT),
			((_D+_DM) << SCORE_SHIFT) + (PRIORITY_DELETION << STATE_SHIFT) + (1 << MISMATCH_SHIFT),
			((NLOW  ) << SCORE_SHIFT) + (PRIORITY_DELETION << STATE_SHIFT) + (1 << MISMATCH_SHIFT)
		}
	};

	static const union { int i[4]; __m128i xmm; } ClearMask = { .i = { CLEAR_MASK, CLEAR_MASK, CLEAR_MASK, CLEAR_MASK} };

	static const union { int i[4]; __m128i xmm; } nlow = { .i = { NLOW<<SCORE_SHIFT, NLOW<<SCORE_SHIFT, NLOW<<SCORE_SHIFT, NLOW<<SCORE_SHIFT} };
	static const unsigned int BitOnes[4] __attribute__((aligned(16))) = { -1, -1, -1, -1 };
	static const union { int i[4]; __m128i xmm; } Ones = { .i = { 1, 1, 1, 1 } };
#endif

#ifdef __DECLARE_VARIABLES__
	const __m128i __FirstColumn = _mm_load_si128((__m128i*)RowInit);
	const __m128i __KOPDInit    = _mm_set1_epi32 (RowInit[DELETION]);
	const __m128i __NLOW        = _mm_set1_epi32(NLOW << SCORE_SHIFT);
#endif 
/*
 * -------------------------------------------------------
 * ALLOCATE GLOBAL MEMORY   
 * -------------------------------------------------------
 */
#ifdef __ALLOCATE_MEMORY__
	// WARNING: We have used the fact that the largest tag is kTagSize - 1 , so that plus 1 is not required !
	const size_t maxMatrixLD = (kTagSize + 3) & ~(3);
	sIOP * restrict WORK = _mm_malloc(63+(2*maxMatrixLD)*sizeof(union sIOP),64);
	lScores * restrict Matrix = _mm_malloc(maxMatrixLD*(1+kGAPPED_ALIGN_GENOME_LENGTH)*sizeof(lScores),64);
	if (WORK == NULL || Matrix == NULL ) {
		fputs("Cannot allocate matrix and work array\n", stderr);
		err = -1 ; goto zeend;
	}
#endif

/*
 * -------------------------------------------------------
 * FILLS IN THE MATRIX 
 * -------------------------------------------------------
 */
#ifdef __FILL_MATRIX__
{
	__m128i KOPD;
	const union sIOP * restrict IOP_R;
	union sIOP * restrict IOP_W = (union sIOP *) WorkSpace;
	
	/*
	 * Initialize Insertion and Match Entrance Line using FirstSequenceProtein
	 */  
	{
		/*
		 * PROFILE COLUMN 0 entrance
		 */
		IOP_W[0].Element.Match     = RowInit[MATCH];
		IOP_W[0].Element.Insertion = RowInit[INSERTION];
		KOPD                       = _mm_set1_epi32(RowInit[DELETION]);
		
		// Store all scores to Matrix
		//operation(__FirstColumn, MatrixPtr, 0, 0, MatrixLD);
		_mm_store_si128(&IDM_Matrix[0], __FirstColumn);
		register union sIOP * restrict pIOP = &IOP_W[1];
		
		/*
		 * LOOP THROUGH THE REST OF THE PROFILE
		 */
		for (size_t iprf=1; iprf<=ReadLength; ++iprf ) {
			KOPD = _mm_and_si128(KOPD, ClearMask.xmm);
			
			// Add KD to Transitions
			__m128i __TransitionsD = _mm_add_epi32(DeletionScore.xmm, KOPD);
			
			KOPD = _mm_shuffle_epi32(__TransitionsD, (DELETION << 6) | (DELETION << 4) | (DELETION << 2) | DELETION);
			
			// Store IOPI and IOPM
			_mm_storeh_pi( &(pIOP->mm), (__m128) __TransitionsD);

			pIOP++;
			
			// Store all scores to Matrix
			//operation(__TransitionsD, MatrixPtr, 0, iprf, MatrixLD);
			_mm_store_si128(&IDM_Matrix[iprf], __TransitionsD);
		}
	}
		
	// Swap and assign Read and write pointers
	IOP_R = IOP_W;
	IOP_W = (union sIOP*) (((uintptr_t) &((union sIOP*) WorkSpace)[MatrixLD] + 63) & ~63);
	
	/*
	 * LOOP THROUGH THE SEQUENCE STRING
	 */
	
	const register __m128i __BitOnes = _mm_load_si128((__m128i*) BitOnes);
	const register __m128i __Zeros   = _mm_setzero_si128();
	__BestScore = nlow.xmm;
	__location = __Zeros;
	__m128i __Currentloc = Ones.xmm;
	for ( size_t iseq=0; iseq < GenomeLength; ++iseq) {
		const register char SequenceDNA = genome[iseq];
		__m128i KOPM = _mm_set1_epi32(IOP_R[0].Element.Match);
		/*
		 * PROFILE COLUMN 0 entrance
		 */
		{
			KOPD = __KOPDInit;
			
			// Store IOPI and IOPM
// 			StoreMatchInsertion( &(IOP_W[0].mm), (__m128) __FirstColumn);
			_mm_storeh_pi(&(IOP_W[0].mm), (__m128) __FirstColumn);
			
			// Store all scores to Matrix
// 			operation(__FirstColumn, MatrixPtr, iseq, 0, MatrixLD);
			_mm_store_si128(&IDM_Matrix[iseq*MatrixLD], __FirstColumn);
		}
				
		/*
		 * LOOP THROUGH THE REST OF THE PROFILE
		 */
		size_t iprf;
		if (SequenceDNA != 'N') {
			for (iprf=1; iprf<ReadLength; ++iprf ) {
				const int ReadIsN        = (ptrTag[iprf-1] == 'N') ? 1 : 0;
				const int DoWeMatch      = (ptrTag[iprf-1] == SequenceDNA) ? 1 : 0;
				const __m128i __movemask = (DoWeMatch || ReadIsN) ? __BitOnes : __Zeros;
				const __m128i __AddtoM   = _mm_blendv_epi8(MismatchScore.xmm, MatchScore.xmm, __movemask);
				
				__m128i __KI = _mm_set1_epi32(IOP_R[iprf].Element.Insertion);
				
				KOPD = _mm_and_si128(KOPD, ClearMask.xmm);
				KOPM = _mm_and_si128(KOPM, ClearMask.xmm);
				__KI = _mm_and_si128(__KI, ClearMask.xmm);
				
				__m128i __KM = _mm_add_epi32(KOPM, __AddtoM);
				__KI = _mm_add_epi32(__KI, InsertionScore.xmm);
				
				KOPM = _mm_set1_epi32(IOP_R[iprf].Element.Match);
				
				__m128i __max = _mm_max_epi32(__KM, __KI);
				
				__m128i __KD = _mm_add_epi32(KOPD, DeletionScore.xmm);
				
				__max = _mm_max_epi32(__max, __KD);
				
				
	// 			StoreMatchInsertion( &IOP_W[iprf].mm, (__m128) __max);
				_mm_storeh_pi(&IOP_W[iprf].mm, (__m128) __max);
	// 			operation(__max, MatrixPtr, iseq, iprf, MatrixLD);
				_mm_store_si128(&IDM_Matrix[(iseq+1)*MatrixLD+iprf], __max);
				
				KOPD = _mm_shuffle_epi32(__max, (DELETION << 6) | (DELETION << 4) | (DELETION << 2) | DELETION);
			}
			
			{
				const int ReadIsN        = (ptrTag[iprf-1] == 'N') ? 1 : 0;
				const int DoWeMatch      = (ptrTag[iprf-1] == SequenceDNA) ? 1 : 0;
				const __m128i __movemask = ( DoWeMatch || ReadIsN) ? __BitOnes : __Zeros;
				const __m128i __AddtoM   = _mm_blendv_epi8(MismatchScore.xmm, MatchScore.xmm, __movemask);
				
				__m128i __KI = _mm_set1_epi32(IOP_R[iprf].Element.Insertion);
				
				KOPD = _mm_and_si128(KOPD, ClearMask.xmm);
				KOPM = _mm_and_si128(KOPM, ClearMask.xmm);
				__KI = _mm_and_si128(__KI, ClearMask.xmm);
				
				__m128i __KM = _mm_add_epi32(KOPM, __AddtoM);
				__KI = _mm_add_epi32(__KI, InsertionScore.xmm);
				
				KOPM = _mm_set1_epi32(IOP_R[iprf].Element.Match);
				
				__m128i __max = _mm_max_epi32(__KM, __KI);
				
				__m128i __KD = _mm_add_epi32(KOPD, DeletionScore.xmm);
				
				__max = _mm_max_epi32(__max, __KD);
				
				const __m128i __mask = _mm_cmpgt_epi32(__max, __BestScore);
				__BestScore = _mm_blendv_epi8(__BestScore, __max, __mask);
				__location  = _mm_blendv_epi8(__location, __Currentloc, __mask);
				
	// 			StoreMatchInsertion( &IOP_W[iprf].mm, (__m128) __max);
				_mm_storeh_pi(&IOP_W[iprf].mm, (__m128) __max);
	// 			operation(__max, MatrixPtr, iseq, iprf, MatrixLD);
				_mm_store_si128(&IDM_Matrix[(1+iseq)*MatrixLD+iprf], __max);
				
				__Currentloc = _mm_add_epi32(__Currentloc, Ones.xmm);
			}
		}
		else {
			for (iprf=1; iprf<ReadLength; ++iprf ) {
				__m128i __KI = _mm_set1_epi32(IOP_R[iprf].Element.Insertion);
				
				KOPD = _mm_and_si128(KOPD, ClearMask.xmm);
				KOPM = _mm_and_si128(KOPM, ClearMask.xmm);
				__KI = _mm_and_si128(__KI, ClearMask.xmm);
				
				__m128i __KM = _mm_add_epi32(KOPM, SequenceNMatchScore.xmm);
				__KI = _mm_add_epi32(__KI, InsertionScore.xmm);
				
				KOPM = _mm_set1_epi32(IOP_R[iprf].Element.Match);
				
				__m128i __max = _mm_max_epi32(__KM, __KI);
				
				__m128i __KD = _mm_add_epi32(KOPD, DeletionScore.xmm);
				
				__max = _mm_max_epi32(__max, __KD);
				
				
	// 			StoreMatchInsertion( &IOP_W[iprf].mm, (__m128) __max);
				_mm_storeh_pi(&IOP_W[iprf].mm, (__m128) __max);
	// 			operation(__max, MatrixPtr, iseq, iprf, MatrixLD);
				_mm_store_si128(&IDM_Matrix[(iseq+1)*MatrixLD+iprf], __max);
				
				KOPD = _mm_shuffle_epi32(__max, (DELETION << 6) | (DELETION << 4) | (DELETION << 2) | DELETION);
			}
			
			{
				__m128i __KI = _mm_set1_epi32(IOP_R[iprf].Element.Insertion);
				
				KOPD = _mm_and_si128(KOPD, ClearMask.xmm);
				KOPM = _mm_and_si128(KOPM, ClearMask.xmm);
				__KI = _mm_and_si128(__KI, ClearMask.xmm);
				
				__m128i __KM = _mm_add_epi32(KOPM, SequenceNMatchScore.xmm);
				__KI = _mm_add_epi32(__KI, InsertionScore.xmm);
				
				KOPM = _mm_set1_epi32(IOP_R[iprf].Element.Match);
				
				__m128i __max = _mm_max_epi32(__KM, __KI);
				
				__m128i __KD = _mm_add_epi32(KOPD, DeletionScore.xmm);
				
				__max = _mm_max_epi32(__max, __KD);
				
				const __m128i __mask = _mm_cmpgt_epi32(__max, __BestScore);
				__BestScore = _mm_blendv_epi8(__BestScore, __max, __mask);
				__location  = _mm_blendv_epi8(__location, __Currentloc, __mask);
				
	// 			StoreMatchInsertion( &IOP_W[iprf].mm, (__m128) __max);
				_mm_storeh_pi(&IOP_W[iprf].mm, (__m128) __max);
	// 			operation(__max, MatrixPtr, iseq, iprf, MatrixLD);
				_mm_store_si128(&IDM_Matrix[(1+iseq)*MatrixLD+iprf], __max);
				
				__Currentloc = _mm_add_epi32(__Currentloc, Ones.xmm);
			}
		}
		// Swap Read and Write pointers
		const register union sIOP * const ptr = IOP_W;
		IOP_W = (union sIOP *) IOP_R;
		IOP_R = ptr;
	}
}
#endif
/*
 * -------------------------------------------------------
 * GET BEST ALIGNMENT 
 * -------------------------------------------------------
 */
#ifdef __GET_ALIGNMENT__
{
	int iprf = ReadLength;
	int StateOffset = (EndOnMatchOrDeletion & 0x1) ? DELETION : MATCH;	
	unsigned int previous_State = (EndOnMatchOrDeletion & 0x1) ? PRIORITY_DELETION : PRIORITY_MATCH;
	unsigned int StateCount = 0;
	unsigned int current_State;
	Align[1] = iseq;
	
	while ( iprf >= 0 && iseq >= 0)
	{
		unsigned int* ptr         = (unsigned int*) &IDM_Matrix[MatrixLD*iseq+iprf];
		const unsigned int lScore = ptr[StateOffset];
		            current_State = (lScore & STATE_MASK) >> STATE_SHIFT;
		
		if (current_State != previous_State) {
			ptrCigar -= 2;
			if (ptrCigar < cigarStorageBottomLimit) { ret = -1; goto OUT;}
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
				
				StateOffset = MATCH;
				if (lScore & MISMATCH_MASK) {
					ptrMismatch -= 2;
					if (ptrMismatch <  mismatchStorageBottomLimit) { ret = -1; goto OUT;}
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
						if (ptrMismatch <  mismatchStorageBottomLimit) { ret = -1; goto OUT;}
						ptrMismatch[0] = (unsigned char) iprf;
						ptrMismatch[1] = ptrTag[iprf];
					}
					else if (ptrTag[iprf] == 'N') {
						ptrMismatch -= 2;
						if (ptrMismatch <  mismatchStorageBottomLimit) { ret = -1; goto OUT;}
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
				StateOffset = INSERTION;
// 				fputc('D', stderr);
				break;
			case PRIORITY_DELETION:
// 				fprintf(stderr, "%i %i %i %c\n", iseq, iprf, ToScore(lScore), 'D');
				StateOffset = DELETION;
				ptrMismatch -= 2;
				if (ptrMismatch <  mismatchStorageBottomLimit) { ret = -1; goto OUT;}
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
	if (ptrCigar < cigarStorageBottomLimit) { ret = -1; goto OUT;}
	ptrCigar[0] = (unsigned char) StateCount;
	ptrCigar[1] = StateLetter[previous_State];
	StateChange += StateChangeAdd[previous_State];
	
	 
	if (iprf >= 0 ) {
		ptrCigar -= 2;
		if (ptrCigar < cigarStorageBottomLimit) { ret = -1; goto OUT;}
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
	OUT:;
// 	fputc('\n', stderr);
// 	printf("CUrrent state: %d iseq %i\n", (int) current_State, iseq);
}
#endif

#ifdef __GET_ALIGNMENT_STATE__
{
	int StateOffset = (EndOnMatchOrDeletion & 0x1) ? DELETION : MATCH;
	unsigned int previous_State = (EndOnMatchOrDeletion & 0x1) ? PRIORITY_DELETION : PRIORITY_MATCH;
	unsigned int current_State;
	unsigned int StateCount = 0U;
	Align[1] = iseq-1;
	ret = -1;
	if (--ptrState < StateLimit) goto OUT;

	while ( iprf >= 0 && iseq >= 0)
	{
		unsigned int* ptr         = (unsigned int*) &IDM_Matrix[MatrixLD*iseq+iprf];
		const unsigned int lScore = ptr[StateOffset];
		            current_State = (lScore & STATE_MASK) >> STATE_SHIFT;
		if (current_State != previous_State) {
			RequiredCigSize++;
			previous_State = current_State;
			StateCount = 1U;
		}
		else {
			if ((++StateCount) >= 256) goto OUT;
		}
		
		switch(current_State) {
			case PRIORITY_MATCH:
				StateOffset = MATCH;
				if (!(lScore & MISMATCH_MASK)) {
					if (genome[iseq-1] == 'N' && ptrTag[iprf-1] != 'N') {
						*ptrState = 'n';
						RequiredMMSize++;
					}
					else if (ptrTag[iprf-1] == 'N') {
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
				--iseq;
				--iprf;
				break;
			case PRIORITY_INSERTION:
				--iseq;
				*ptrState = 'D';
				StateOffset = INSERTION;
				ContainsInDel = 1;
				break;
			case PRIORITY_DELETION:
				StateOffset = DELETION;
				*ptrState = 'I';
				--iprf;
				ContainsInDel = 1;
// 				InsertionCount++;
				RequiredMMSize++;
				break;
			case PRIORITY_EXTRA:
				goto DONE;
				break;
			default:
				fprintf(stderr,"\nUnknown state (%u) encountered from cell (%i,%i)\n", current_State, iseq, iprf);
				assert(0);
				--iprf;
				goto OUT;
				break;
		}
		
		if (iprf == SoftClipBoundary) {
			SoftClipMismatchCount  = MismatchCount;
			SoftClipEvictedMMSize  = RequiredMMSize;
			SoftClipEvictedCigSize = RequiredCigSize;
		}
		if (ptrState-- < StateLimit) goto OUT;	
#ifdef _TEST_
		++counter;
#endif
	}
	DONE:;
	ptrState++;
	assert(iprf == 0);	
	if (iprf >/*=*/ 0) {
		assert(0);
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
				SoftClipMismatchCount  - MismatchCount;
			}
			
// 			InsertionCount += 1+iprf;
			RequiredMMSize  += 1+iprf;
			RequiredCigSize += (current_State == PRIORITY_DELETION) ? 0 : 1; 
			ContainsInDel = 1;
		}
	}
	Align[0] = iseq < 0 ? 0 : (current_State == PRIORITY_DELETION) ? iseq-1 : iseq /*-1+1*/;
	ret = 0;
	OUT:;
}
#endif

/*
 * -------------------------------------------------------
 * FREE GLOBAL MEMORY
 * -------------------------------------------------------
 */
#ifdef __FREE_MEMORY__
	if (WORK) _mm_free(WORK);
	if (Matrix) _mm_free(Matrix);
#endif

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
