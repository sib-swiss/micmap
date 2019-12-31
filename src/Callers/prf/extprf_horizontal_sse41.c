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
/*******************************************************
                        PFTOOLS
 *******************************************************
  Oct 12, 2012 xali1_print_sse41.c
 *******************************************************
 (C) 2012 Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@isb-sib.ch)
 *******************************************************/
#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <smmintrin.h>
#include <assert.h>
#include "pfCompute.h"
#include "sse41_inline_fcts.h"

#define TAG

/* Priority */
enum StatePriority {
  PRIORITY_MATCH     = 0,
  PRIORITY_INSERTION = 1,
  PRIORITY_DELETION  = 2,
  PRIORITY_EXTRA     = 3
};

/*
  * INTEGER FORMAT  -------------------------------------
  *
  * |xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx|xx|
  *                               |  |
  *                               |  +-- State
  *                               +----- Score
  */
# define SCORE_SHIFT               2
# define SCORE_MASK                0xFFFFFFFC
# define STATE_SHIFT               0
# define STATE_MASK                0x00000003
# define LEFT_NEGATIVE_SCORE_MASK  0xC0000000
# define CLEAN_STATE               0xFFFFFFF3
// # define CLEAN_PHASE               0xFFFFFFFC
// # define CLEAN_STATE_AND_PHASE     0xFFFFFFF0


#define TO_SCORE(x) (int)(((x)<0) ? ((((unsigned int)(x)) >> SCORE_SHIFT) | LEFT_NEGATIVE_SCORE_MASK ) : (((unsigned int)x) >> SCORE_SHIFT))
#define MAX(a,b) ((a>b) ? a : b)

union sIOP {
    struct {
      int Match;
      int Insertion;
    } Element;
    __m64 mm;
};

#define operation(x, inout, row, column, lda) _mm_store_si128(&((inout)[(row)*(lda)+(column)]), x)

static void xali1_print_sse41(const struct Profile * const restrict prf,
                              const unsigned char * const restrict Sequence,
                              union lScores * const restrict matrix, int * const restrict WORK,
                              unsigned short int * const restrict NumberOfInsertions,
                              const size_t BSEQ, const size_t LSEQ)
/*
 * WARNING: for SSE version, WORK should be aligned on cache size (64b) and 4 times the (profile size + 1)*sizeof(int)
 *          + 63 to align to cache line
 *          matrix should be of size at least (profile size + 1)*(sequence size+1) and aligned on 16b.
 */
{
  int KOPD;
  const union sIOP * restrict IOP_R;
  union sIOP * restrict IOP_W = (union sIOP*) WORK;
  const __m128i __MatchMask     = _mm_set1_epi32(PRIORITY_MATCH);
  const __m128i __InsertionMask = _mm_set1_epi32(PRIORITY_INSERTION);
  const __m128i __DeletionMask  = _mm_set1_epi32(PRIORITY_DELETION);
  const __m128i __ExtraMask     = _mm_set1_epi32(PRIORITY_EXTRA);
  
  register const TransitionScores * const restrict Transitions = prf->Scores.Insertion.Transitions;
  const StoredIntegerFormat * const restrict Match             = prf->Scores.Match.Alphabet;
  const StoredIntegerFormat * const restrict Insertion         = prf->Scores.Insertion.Alphabet;
  const size_t AlignStep                                       = prf->Scores.Match.AlignStep;
  const size_t prfLength = prf->Length;
  
  /* Set matrix ptr according to BSEQ */
  __m128i * const restrict MatrixPtr = (BSEQ == 0) ? &matrix[1+prfLength].xmm : &matrix[BSEQ].xmm;
  
  /* NOTE: The following part could be replaced and performed only once for a profile as it
   *       is profile dependent. Nevertheless it does a good job loading Match and Transition
   *       matrices into the cache hierarchy.
   */

  /*
   * Initialize Insertion and Match Entrance Line using FirstSequenceProtein
   */  
  {
    register const StoredIntegerFormat * restrict lMatch = (const StoredIntegerFormat *) &Match[_D];
    register const ScoreTuple * restrict FirstSequenceProtein = prf->Scores.Insertion.FirstSequenceProtein;
    /*
     * PROFILE COLUMN 0 entrance
     */
    IOP_W[0].Element.Match     = (int) FirstSequenceProtein[0].To[MATCH];
    IOP_W[0].Element.Insertion = (int) FirstSequenceProtein[0].To[INSERTION];
    KOPD                       = (int) FirstSequenceProtein[0].To[DELETION];


    {
			__m128i __FirstSequenceProtein = LoadStoredIntegerVector(&(FirstSequenceProtein[0]));
			__FirstSequenceProtein = _mm_slli_epi32(__FirstSequenceProtein, 2);
			__FirstSequenceProtein = _mm_or_si128(__FirstSequenceProtein, __ExtraMask);
			// Store all scores to matrix
			operation(__FirstSequenceProtein, MatrixPtr - (1+prfLength), 0, 0, 1+prfLength);
    }
    
    FirstSequenceProtein++;
    register const TransitionScores (* restrict pTransitions) = &Transitions[1];
    register union sIOP * restrict pIOP = &IOP_W[1];
    
    /*
     * LOOP THROUGH THE REST OF THE PROFILE
     */
    for (size_t iprf=1; iprf<=prfLength; ++iprf ) {
      register const int KD = KOPD + (int) *lMatch;
      lMatch += AlignStep;

      // Transform KD into a vector
      __m128i __KD = _mm_set1_epi32(KD);
      // Load Transitions and  Convert signed WORD into signed DWORD
      __m128i __TransitionsD = LoadStoredIntegerVector(&(pTransitions->From[DELETION]));

      // Add KD to Transitions
      __TransitionsD = _mm_add_epi32(__TransitionsD, __KD);

      // Move to next profile transitions
      pTransitions++;

      // Load FirstSequenceProtein
      __m128i __FirstSequenceProtein = LoadStoredIntegerVector(&(FirstSequenceProtein[0]));

      // Move to next profile First Sequence
      FirstSequenceProtein++;
      
      // Paste index in lowest 2 bits
      __TransitionsD = _mm_slli_epi32(__TransitionsD, SCORE_SHIFT);
      __FirstSequenceProtein = _mm_slli_epi32(__FirstSequenceProtein, SCORE_SHIFT);
      __TransitionsD = _mm_or_si128(__TransitionsD, __DeletionMask);
      __FirstSequenceProtein = _mm_or_si128(__FirstSequenceProtein, __ExtraMask);
      
      // Get maximum ( this is SSE 4.1 )
      __m128i __max = _mm_max_epi32(__TransitionsD, __FirstSequenceProtein);

      // Store all scores to matrix
      operation(__max, MatrixPtr - (1+prfLength), 0, iprf, 1+prfLength);

      // Clean extra bits
      __max = _mm_srai_epi32(__max, SCORE_SHIFT);
      
      // Store IOPI and IOPM
      StoreMatchInsertion( &(pIOP->mm), (__m128) __max);
      pIOP++;

      // Set KOPD ( this is SSE 4.1 )
      KOPD = _mm_extract_epi32(__max, DELETION);
    }
  }

  // Swap and assign Read and write pointers
  IOP_R = IOP_W;
  IOP_W = (union sIOP*) (((uintptr_t) &WORK[2*(prf->Length+1)] + 63) & ~63);

  /*
   * LOOP THROUGH THE SEQUENCE STRING
   */
  for ( size_t iseq=BSEQ; iseq < LSEQ-1; ++iseq) {
    register const size_t j1 = (size_t) Sequence[iseq];
    register const StoredIntegerFormat * restrict lInsertion = Insertion;
    /*
     * PROFILE COLUMN 0 entrance
     */
    {
      register const int KI = NLOW;

      // Transform KI into a vector
      __m128i __KI = _mm_set1_epi32(KI);
      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsI = LoadStoredIntegerVector(&(Transitions[0].From[INSERTION]));
     // Add KI to Transition
      __TransitionsI = _mm_add_epi32(__TransitionsI, __KI);

       // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsX = LoadStoredIntegerVector(&(Transitions[0].From[EXTRA]));

      // Insert lScore into __TransitionsX not necessary as profile loading store NLOW there automatically
      //__TransitionsX = _mm_insert_epi32(__TransitionsX, NLOW, DUMMY);
      
      // Paste index in lowest 2 bits
      __TransitionsI = _mm_slli_epi32(__TransitionsI, SCORE_SHIFT);
      __TransitionsX = _mm_slli_epi32(__TransitionsX, SCORE_SHIFT);
      __TransitionsI = _mm_or_si128(__TransitionsI, __InsertionMask);
      __TransitionsX = _mm_or_si128(__TransitionsX, __ExtraMask);

      // Get maximum ( this is SSE 4.1 )
      __m128i __max = _mm_max_epi32(__TransitionsI, __TransitionsX);
      
      // Store all scores to matrix
      operation(__max, MatrixPtr, iseq, 0, 1+prfLength); //_mm_store_si128(&pmatrix[iprf-1], __max1);

      // Clean extra bits
      __max = _mm_srai_epi32(__max, SCORE_SHIFT);

			// Store IOPI and IOPM
      StoreMatchInsertion( &(IOP_W[0].mm), (__m128) __max);

      // Store KOPD
      KOPD = _mm_extract_epi32(__max, DELETION);

      // Backup new score to xmm register
      //lScore = _mm_extract_epi32(__max, DUMMY);
    }

//     lInsertion += AlignStep;
    register const StoredIntegerFormat * restrict lMatch = Match;

    /*
     * LOOP THROUGH THE REST OF THE PROFILE
     */
		int KOPM = IOP_R[0].Element.Match;
		int KOPI = IOP_R[0].Element.Insertion;
    for (size_t iprf=1; iprf<=prfLength; ++iprf ) {
      const int KM = KOPM + (int) lMatch[j1];
      const int KI = KOPI + (int) lInsertion[j1];
      const int KD = KOPD + (int) lMatch[_D];

      lMatch     += AlignStep;
      lInsertion += AlignStep;

      KOPM = IOP_R[iprf].Element.Match;
			KOPI = IOP_R[iprf].Element.Insertion;

      // Transform KM into a vector
      __m128i __KM = _mm_set1_epi32(KM);
      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsM = LoadStoredIntegerVector(&(Transitions[iprf].From[MATCH]));
      // Add KM to Transition
      __TransitionsM = _mm_add_epi32(__TransitionsM, __KM);


      // Transform KI into a vector
      __m128i __KI = _mm_set1_epi32(KI);
      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsI = LoadStoredIntegerVector(&(Transitions[iprf].From[INSERTION]));
      // Add KI to Transition
      __TransitionsI = _mm_add_epi32(__TransitionsI, __KI);

      // Paste index in lowest 2 bits
      __TransitionsM = _mm_slli_epi32(__TransitionsM, SCORE_SHIFT);
      __TransitionsI = _mm_slli_epi32(__TransitionsI, SCORE_SHIFT);
      __TransitionsM = _mm_or_si128(__TransitionsM, __MatchMask);
      __TransitionsI = _mm_or_si128(__TransitionsI, __InsertionMask);

      // Get maximum ( this is SSE 4.1 )
      __m128i __max1 = _mm_max_epi32(__TransitionsM, __TransitionsI);      

      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsX = LoadStoredIntegerVector(&(Transitions[iprf].From[EXTRA]));
       // Insert lscore into TransitionX not necessary as profile loading should have already done it
      // __TransitionsX = _mm_insert_epi32(__TransitionsX, lScore, DUMMY);

      // Transform KD into a vector
      __m128i __KD = _mm_set1_epi32(KD);
      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsD = LoadStoredIntegerVector(&(Transitions[iprf].From[DELETION]));
      // Add KD to Transition
      __TransitionsD = _mm_add_epi32(__TransitionsD, __KD);

      // Paste index in lowest 2 bits
      __TransitionsX = _mm_slli_epi32(__TransitionsX, SCORE_SHIFT);
      __TransitionsD = _mm_slli_epi32(__TransitionsD, SCORE_SHIFT);
      __TransitionsX = _mm_or_si128(__TransitionsX, __ExtraMask);
      __TransitionsD = _mm_or_si128(__TransitionsD, __DeletionMask);

      // Get maximum ( this is SSE 4.1 )
      __m128i __max2 = _mm_max_epi32(__TransitionsD, __TransitionsX);
      __max1 = _mm_max_epi32(__max1, __max2);

      // Store all scores to matrix
      operation(__max1, MatrixPtr, iseq, iprf, 1+prfLength); //_mm_store_si128(&pmatrix[iprf-1], __max1);

      // Clean extra bits
      __max1 = _mm_srai_epi32(__max1, SCORE_SHIFT);

      // Store IOPI and IOPM
      StoreMatchInsertion( &IOP_W[iprf].mm, (__m128) __max1);


      // Set KOPD ( this is SSE 4.1 )
      KOPD = _mm_extract_epi32(__max1, DELETION);

//       lScore = _mm_extract_epi32(__max1, DUMMY);

    }

    // Swap Read and Write pointers
    const register union sIOP * const ptr = IOP_W;
    IOP_W = (union sIOP *) IOP_R;
    IOP_R = ptr;
  }

  /*
   * Last position on the Sequence using LastSequenceProtein
   */
  {
    register const StoredIntegerFormat * restrict lInsertion = Insertion;
    const int j1 = (int) Sequence[LSEQ-1];
    /*
     * PROFILE COLUMN 0 entrance
     */

    int KI = NLOW;

    KOPD = MAX( KI + (int) Transitions[0].Element[_ID], (int) Transitions[0].Element[_XD] );
    register const ScoreTuple * const restrict LastSequenceProtein = prf->Scores.Insertion.LastSequenceProtein;

    register const StoredIntegerFormat * restrict lMatch = Match;
    lInsertion += AlignStep;

		/*
     * LOOP THROUGH THE REST OF THE PROFILE
     */
		int KOPM = IOP_R[0].Element.Match;
		int KOPI = IOP_R[0].Element.Insertion;
    for (size_t iprf=1; iprf<=prfLength; ++iprf) {
      const int KM = KOPM + lMatch[j1];
      KI           = KOPI + lInsertion[j1];
      const int KD = KOPD + lMatch[_D];

      lMatch     += AlignStep;
      lInsertion += AlignStep;

      KOPM = IOP_R[iprf].Element.Match;

      // Transform KM into a vector
      __m128i __KM = _mm_set1_epi32(KM);
      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsM = LoadStoredIntegerVector(&(Transitions[iprf].From[MATCH]));
      // Insert LastProteinSequence
      __TransitionsM = _mm_insert_epi32(__TransitionsM, (int) LastSequenceProtein[iprf].From[MATCH], DUMMY);
      // Add KM to Transition
      __TransitionsM = _mm_add_epi32(__TransitionsM, __KM);


      // Transform KI into a vector
      __m128i __KI = _mm_set1_epi32(KI);
      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsI = LoadStoredIntegerVector(&(Transitions[iprf].From[INSERTION]));
      // Insert LastProteinSequence
      __TransitionsI = _mm_insert_epi32(__TransitionsI, (int) LastSequenceProtein[iprf].From[INSERTION], DUMMY);
      // Add KI to Transition
      __TransitionsI = _mm_add_epi32(__TransitionsI, __KI);

      // Paste index in lowest 2 bits
      __TransitionsM = _mm_slli_epi32(__TransitionsM, SCORE_SHIFT);
      __TransitionsI = _mm_slli_epi32(__TransitionsI, SCORE_SHIFT);
      __TransitionsM = _mm_or_si128(__TransitionsM, __MatchMask);
      __TransitionsI = _mm_or_si128(__TransitionsI, __InsertionMask);

      // Get maximum ( this is SSE 4.1 )
      __m128i __max1 = _mm_max_epi32(__TransitionsM, __TransitionsI);

      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsX = LoadStoredIntegerVector(&(Transitions[iprf].From[EXTRA]));
      // Insert lscore into TransitionX not necessary as profile loading should have already done it
      // __TransitionsX = _mm_insert_epi32(__TransitionsX, lScore, DUMMY);

      // Transform KD into a vector
      __m128i __KD = _mm_set1_epi32(KD);
      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsD = LoadStoredIntegerVector(&(Transitions[iprf].From[DELETION]));
      // Insert LastProteinSequence
      __TransitionsD = _mm_insert_epi32(__TransitionsD, (int) LastSequenceProtein[iprf].From[DELETION], DUMMY);
      // Add KD to Transition
      __TransitionsD = _mm_add_epi32(__TransitionsD, __KD);

      // Paste index in lowest 2 bits
      __TransitionsX = _mm_slli_epi32(__TransitionsX, SCORE_SHIFT);
      __TransitionsD = _mm_slli_epi32(__TransitionsD, SCORE_SHIFT);
      __TransitionsX = _mm_or_si128(__TransitionsX, __ExtraMask);
      __TransitionsD = _mm_or_si128(__TransitionsD, __DeletionMask);

      // Get maximum ( this is SSE 4.1 )
      __m128i __max2 = _mm_max_epi32(__TransitionsD, __TransitionsX);
      __max1 = _mm_max_epi32(__max1, __max2);

      // Store all scores to matrix
      operation(__max1, MatrixPtr, LSEQ-1, iprf, prfLength+1); //_mm_store_si128(&pmatrix[iprf-1], __max1);

      // Clean extra bits
      __max1 = _mm_srai_epi32(__max1, SCORE_SHIFT);

      // Set KOPD ( this is SSE 4.1 )
      KOPD = _mm_extract_epi32(__max1, DELETION);
    }
  }
}
//---------------------------------------------------------------

static size_t GetNextBestAlignment(const union lScores * const matrix, const _Bool * const HasAPath, Alignment_t * const Alignment,
                                   const size_t SeqLength, const size_t prfLength)
{
	const size_t matrix_ld = prfLength + 1;
	const int (* restrict Scores)[4]  = (const int (* restrict)[4]) &matrix[matrix_ld+1];
	int index_r = -1, index_c = -1;
	int extrema = NLOW;
	const int lSeqLength = (int) SeqLength;
	for (int i=0; i<lSeqLength; ++i) {
		if (!HasAPath[1+i]) {
			for (int j=0; j<prfLength; ++j) {
				const int Value = TO_SCORE(Scores[j][EXTRA]);
				if (Value >= extrema) {
					extrema = Value;
					index_r = i;
					index_c = j;
				}
			}
		}
		Scores += matrix_ld;
	}
	Alignment->Score = extrema;
	if (extrema == NLOW) return 0;
	
	/* Correct for initial entrance column and row */
	index_c += 1;
	index_r += 1;
	
	Alignment->Matrix.row.End      = index_r;
	Alignment->Matrix.column.End   = index_c;
	/* Correcting indexing to account extra top row and left column */
	Alignment->Region.profile.End  = index_c - 1;
	Alignment->Region.sequence.End = index_r - 1;
	
	size_t counter     = 1; /* take into account initial X */
	int iprf           = index_c;
	size_t State       = EXTRA;
	Scores  = (const int (* restrict)[4]) matrix;
	while ( iprf >= 0 && index_r >= 0) {
		const unsigned int MoveToState = (Scores[index_r*matrix_ld+iprf][State] & STATE_MASK);
		switch (MoveToState)
		{
			case PRIORITY_MATCH << STATE_SHIFT:
				iprf    -= 1;
				index_r -= 1;
				State = MATCH;
				break;
			case PRIORITY_INSERTION << STATE_SHIFT:
				iprf    -= 1;
				index_r -= 1;
				State = INSERTION;
				break;
			case PRIORITY_DELETION << STATE_SHIFT:
					--iprf;
				State = DELETION;
				break;
			case PRIORITY_EXTRA << STATE_SHIFT:
				++counter;
				goto OUT;
				break;
		}
		++counter;
// 		index_r = (index_r < 0 ) ? 0 : index_r;
	}
	return 0;
	
	OUT:;
	Alignment->Matrix.column.Begin = iprf;
	Alignment->Matrix.row.Begin    = index_r;
	assert(Alignment->Matrix.column.Begin >= 0);
	assert(Alignment->Matrix.row.Begin  >= 0);
	
	/* Correcting indexing to account extra top row and left column */
	Alignment->Region.profile.Begin  = (State == INSERTION) ? iprf - 1: iprf /*-1+1*/;
// 	Alignment->Region.sequence.Begin = (State == DELETION) ? index_r - 1 : index_r /*-1+1*/;
	Alignment->Region.sequence.Begin = index_r;
	assert(Alignment->Region.sequence.Begin >= 0);
	assert(Alignment->Region.profile.Begin >= 0);
	
	Alignment->Cycles = 0;
	
	return ++counter; // Account for ternimal C code character \0
}
//---------------------------------------------------------------

static void GetAlignmentSequence(const union lScores * const restrict matrix,
			  const unsigned char * const restrict Sequence, unsigned char * const restrict AlignmentSequence,
			  const Alignment_t * const restrict Alignment, const size_t SeqLength, const size_t prfLength)
{
    const int (*Scores)[4]  = (const int (*)[4]) matrix;

    int index       = Alignment->Matrix.row.End;
    size_t counter  = 0;
    int iprf        = Alignment->Matrix.column.End;
    size_t State    = EXTRA;
    const size_t ld = prfLength+1;
    
    while ( iprf >= 0 && index >= 0)
    {
			const unsigned int Move = Scores[index*ld+iprf][State] & 0x3;
			const unsigned char C = (Sequence[index-1] > 'Z') ? Sequence[index-1] - ('a' - 'A') : Sequence[index-1];
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
					--iprf;
			    State = INSERTION;
			    AlignmentSequence[counter] = (char) (  C + ( 'a' - 'A'));
			    break;
			  case PRIORITY_DELETION:
			    --iprf;
			    State = DELETION;
			    AlignmentSequence[counter] = '-';
			    break;
			  case PRIORITY_EXTRA:
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
    fputs("Potential issue as alignment does not start with an ENTRY point in std GetAlignmentSequence!\n", stderr);

    OUT:;

    /* Reverse the string */
    unsigned char * BackPtr = &AlignmentSequence[counter-1];

    for (size_t i=0; i<counter/2; ++i) {
			const unsigned char c = AlignmentSequence[i];
			AlignmentSequence[i] = *BackPtr;
			*BackPtr-- = c;
    }
    if (counter & 0x1) {
			const unsigned char c = AlignmentSequence[counter/2];
			AlignmentSequence[counter/2] = *BackPtr;
			*BackPtr = c;
    }
    AlignmentSequence[counter] = '\0';
}
//---------------------------------------------------------------

const Compute_t ExtendedProfile_sse41 = {
	.BuildMatrix = xali1_print_sse41,
	.GetBestScore = NULL,
	.GetBestAlignment = NULL,
	.GetNextBestAlignment = GetNextBestAlignment,
	.GetStateSequence = NULL,
	.GetAlignmentSequence = GetAlignmentSequence,
	.GetAlignments = NULL,
	.ScoreMask = SCORE_MASK,
	.LeftNegativeScoreMask = LEFT_NEGATIVE_SCORE_MASK,
	.ScoreShift = SCORE_SHIFT,
	.StateMask = STATE_MASK,
	.StateShift = STATE_SHIFT,
	.MatrixColumnMultiplier = 1U,
	.MatrixExtraColumn = 1U,
	.MatrixExtraRow = 1U,
	.WorkCellSize = 4*sizeof(int)
};

#undef MAX

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
