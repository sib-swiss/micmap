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
/***************************************************************************************************
                        PFTOOLS
 ***************************************************************************************************
  Oct 3, 2011 pfProfileInline.h
 ***************************************************************************************************
 (C) 2011 Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@isb-sib.ch)
 ***************************************************************************************************/
#ifdef BUILD_LIBRARY
#define PF_STATIC_INLINE(type) PFEXPORT type
#define PF_EXTERN_INLINE(type) PFEXPORT type
#else
#define PF_STATIC_INLINE(type) static inline type __ALWAYS_INLINE
#define PF_EXTERN_INLINE(type) extern inline type __ALWAYS_INLINE
#endif

PF_STATIC_INLINE(void) InitializeDefault(union Scores * const matrices,
                                         char * const MatchSymbol,
                                         char * const InsertionSymbol )
{
  /* MATCHES */
  memset(matrices->Match.Alphabet, 0, (ALPHABET_MEMORY_SIZE)*sizeof(StoredIntegerFormat));
  /* INSERTIONS */
  memset(matrices->Insertion.Alphabet,    0, (ALPHABET_MEMORY_SIZE)*sizeof(StoredIntegerFormat));
  memset(matrices->Insertion.Boundaries,  0, (INSERTION_BOUNDARIES_SIZE)*sizeof(StoredIntegerFormat));
  memset(matrices->Insertion.Transitions, 0, sizeof(TransitionScores));

  matrices->Insertion.Boundaries[_BI] = NLOW;
  matrices->Insertion.Boundaries[_BD] = NLOW;
  matrices->Insertion.Boundaries[_BE] = NLOW;
  matrices->Insertion.Boundaries[_DE] = NLOW;
  matrices->Insertion.Boundaries[_IE] = NLOW;

  matrices->Insertion.Transitions->Element[_MI] = NLOW;
  matrices->Insertion.Transitions->Element[_MD] = NLOW;
  matrices->Insertion.Transitions->Element[_ID] = NLOW;
  matrices->Insertion.Transitions->Element[_IM] = NLOW;
  matrices->Insertion.Transitions->Element[_DM] = NLOW;
  matrices->Insertion.Transitions->Element[_DI] = NLOW;

  matrices->Insertion.Alphabet[_STOP] = NLOW_16;
  matrices->Match.Alphabet[_STOP]     = NLOW_16;
#ifdef _BEST_IS_NEGATIVE_
  matrices->Insertion.Alphabet[0] = NLOW_16;
  matrices->Match.Alphabet[0]     = NLOW_16;
#endif
  
  *MatchSymbol = 'X';
  *InsertionSymbol = '-';
}

PF_STATIC_INLINE(void) FreeScores(union Scores * const matrices)
{
  if (matrices->Insertion.Alphabet != NULL) _mm_free(matrices->Insertion.Alphabet);
  if (matrices->Insertion.Boundaries != NULL) _mm_free(matrices->Insertion.Boundaries);
  if (matrices->Insertion.Transitions != NULL) _mm_free(matrices->Insertion.Transitions);
  if (matrices->Insertion.FirstSequenceProtein != NULL) _mm_free(matrices->Insertion.FirstSequenceProtein);
  if (matrices->Insertion.LastSequenceProtein != NULL) _mm_free(matrices->Insertion.LastSequenceProtein);
  matrices->Insertion.AlignStep = 0;
  if (matrices->Match.Alphabet != NULL) _mm_free(matrices->Match.Alphabet);
}

PF_STATIC_INLINE(int) AllocateScores(union Scores * const matrices, const size_t Alphabet_Length, const size_t Profile_Length)
{
  register void * ptr;
  const size_t Aligned_Alphabet_Length = (Alphabet_Length+ALPHABET_EXTRA_LETTERS + 15) & ~15 ;
  matrices->Insertion.AlignStep = Aligned_Alphabet_Length;
//   matrices->Match.AlignStep = Aligned_Alphabet_Length;

  ptr =  _mm_malloc(Aligned_Alphabet_Length*(Profile_Length+1)*sizeof(StoredIntegerFormat), 64);
  if (ptr == NULL) goto FIN;
  matrices->Match.Alphabet = (StoredIntegerFormat *) ptr;

  ptr =  _mm_malloc(Aligned_Alphabet_Length*(Profile_Length+1)*sizeof(StoredIntegerFormat), 64);
  if (ptr == NULL) goto FIN;
  matrices->Insertion.Alphabet = (StoredIntegerFormat *) ptr;

  ptr =  _mm_malloc(INSERTION_BOUNDARIES_SIZE*(Profile_Length+1)*sizeof(StoredIntegerFormat), 64);
  if (ptr == NULL) goto FIN;
  matrices->Insertion.Boundaries  = (StoredIntegerFormat *) ptr;

  ptr =  _mm_malloc((Profile_Length+1)*sizeof(TransitionScores), 64);
   if (ptr == NULL) goto FIN;
  matrices->Insertion.Transitions = (TransitionScores *) ptr;

  ptr =  _mm_malloc((Profile_Length+1)*sizeof(ScoreTuple), 64);
  if (ptr == NULL) goto FIN;
  matrices->Insertion.FirstSequenceProtein = (ScoreTuple *) ptr;

  ptr =  _mm_malloc((Profile_Length+1)*sizeof(ScoreTuple), 64);
  if (ptr == NULL) goto FIN;
  matrices->Insertion.LastSequenceProtein = (ScoreTuple *) ptr;

  return 0;
FIN:
   FreeScores(matrices);
   return 1;
}

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
