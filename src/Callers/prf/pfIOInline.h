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
#ifndef PF_IO_INLINE_H
#define PF_IO_INLINE_H
/*
 ************************************************************************************************
 *                                   INLINE FUNCTIONS                                           *
 ************************************************************************************************
 */

extern inline void __ALWAYS_INLINE NextInsertionProfile( struct SInsertion * Insertion)
{
  const size_t step       = Insertion->AlignStep;
  Insertion->Alphabet    += step;
  Insertion->Boundaries  += INSERTION_BOUNDARIES_SIZE;
  Insertion->Transitions ++;
}

extern inline void __ALWAYS_INLINE PreviousInsertionProfile( struct SInsertion * Insertion)
{
  const size_t step       = Insertion->AlignStep;
  Insertion->Alphabet    -= step;
  Insertion->Boundaries  -= INSERTION_BOUNDARIES_SIZE;
  Insertion->Transitions --;
}

extern inline void __ALWAYS_INLINE CopyPreviousInsertionProfile(struct SInsertion * Insertion)
{
  const size_t step       = Insertion->AlignStep;
  memcpy(Insertion->Alphabet, Insertion->Alphabet - step, step*sizeof(StoredIntegerFormat));
  memcpy(Insertion->Boundaries, Insertion->Boundaries - INSERTION_BOUNDARIES_SIZE, INSERTION_BOUNDARIES_SIZE*sizeof(StoredIntegerFormat));
  memcpy(Insertion->Transitions, Insertion->Transitions - 1, sizeof(TransitionScores));
}

extern inline void __ALWAYS_INLINE NextMatchProfile( struct SMatch * Match)
{
  const size_t step = Match->AlignStep;
  Match->Alphabet  += step;
}

extern inline void __ALWAYS_INLINE PreviousMatchProfile( struct SMatch * Match)
{
  const size_t step = Match->AlignStep;
  Match->Alphabet  -= step;
}

extern inline size_t __ALWAYS_INLINE GetInsertionMemory(const char * const key, struct SInsertion * const Insertion, StoredIntegerFormat ** pointer)
{
  /* return 0 if single value, 1 if vector and 2 if error */

  StoredIntegerFormat * Alphabet    = Insertion->Alphabet;
  StoredIntegerFormat * Boundaries  = Insertion->Boundaries;
  StoredIntegerFormat * Transitions = Insertion->Transitions->Element;

  switch(key[0]) {
    case 'I':
      switch(key[1]) {
        case '\0': *pointer = Alphabet    +   1; return 1; break;
        case '0' : *pointer = Alphabet    +   0; return 0; break;
        case 'M' : *pointer = Transitions + _IM; return 0; break;
        case 'I' : *pointer = Transitions + _II; return 0; break;
        case 'D' : *pointer = Transitions + _ID; return 0; break;
        case 'E' : *pointer = Boundaries  + _IE; return 0; break;
        default  : return 2;
      };
      break;
    case 'B':
      switch(key[1]) {
        case '0' : *pointer = Boundaries + _B0; return 0; break;
        case '1' : *pointer = Boundaries + _B1; return 0; break;
        case 'M' : *pointer = Boundaries + _BM; return 0; break;
        case 'I' : *pointer = Boundaries + _BI; return 0; break;
        case 'D' : *pointer = Boundaries + _BD; return 0; break;
        case 'E' : *pointer = Boundaries + _BE; return 0; break;
        default  : return 2;
      };
      break;
    case 'E':
      switch(key[1]) {
        case '0' : *pointer = Boundaries + _E0; return 0; break;
        case '1' : *pointer = Boundaries + _E1; return 0; break;
        default  : return 2;
      };
      break;
    case 'M':
      switch(key[1]) {
        case 'M' : *pointer = Transitions + _MM; return 0; break;
        case 'I' : *pointer = Transitions + _MI; return 0; break;
        case 'D' : *pointer = Transitions + _MD; return 0; break;
	case 'E' : *pointer = Boundaries  + _ME; return 0; break;
        default  : return 2;
      };
      break;
    case 'D':
      switch(key[1]) {
        case 'M' : *pointer = Transitions + _DM; return 0; break;
        case 'I' : *pointer = Transitions + _DI; return 0; break;
        case 'D' : *pointer = Transitions + _DD; return 0; break;
        case 'E' : *pointer = Boundaries  + _DE; return 0; break;
        default  : return 2;
      };
      break;
    default:
      return 2;
  };
}

static inline void __ALWAYS_INLINE FreeAverage(SAverage * const Average)
{
    if (Average->Weights) _mm_free(Average->Weights);
}
#endif /* PF_IO_INLINE_H */

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
