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
#ifndef _STATISTICS_H 
#define _STATISTICS_H 
#include "virt_chr.h"

//--------------------------------------------------------------- 
// STRUCTURE DEFINITIONS 
//--------------------------------------------------------------- 
typedef struct COUNTS_TYPE_s { 
  unsigned int pppf;    // Perfect matches no N 
  unsigned int mNppf;   // Only mismatches that are Ns 
  unsigned int mppf;    // Only mismatches, no Ns 
  unsigned int mppfN;   // Only mismatches with some Ns 
  unsigned int gppf;    // gap alignment no Ns 
  unsigned int gppfN;   // gap alignment with some Ns 
  unsigned int smppf;   // Soft clipped with only mismatches no Ns 
  unsigned int smppfN;  // Soft clipped with only mismatches with some Ns 
  unsigned int sgppf;   // Soft clipped with gap alignment no Ns 
  unsigned int sgppfN;  // Soft clipped with gap alignment with some Ns 
  unsigned int genomeN; // Ns within genome 
} COUNTS_TYPE_t; 

typedef struct COUNTS_struct {
  COUNTS_TYPE_t Pair[kMAXCHRcnt]; 
  COUNTS_TYPE_t TooFar[kMAXCHRcnt]; 
  COUNTS_TYPE_t HalfMap[kMAXCHRcnt]; 
  COUNTS_TYPE_t Chimera;
  unsigned int upf;         // unmapped  
  unsigned int upf_encoded; // unmapped due to criteria but correclty encoded 
} COUNTS;

//--------------------------------------------------------------- 
// FUNCTIONS 
//---------------------------------------------------------------
void reportStatistics(const COUNTS * const restrict cnts, FILE* const out, unsigned long elapsed);

//--------------------------------------------------------------- 
// INLINED FUNCTIONS 
//---------------------------------------------------------------
static inline void __attribute__((always_inline))
cumulateStatistics(COUNTS * const restrict cnts, const size_t N) {
	unsigned int (*const restrict CountPtr)[sizeof(COUNTS)/sizeof(unsigned int)] = (unsigned int (*const restrict)[sizeof(COUNTS)/sizeof(unsigned int)]) cnts;
	for (size_t iCnt=1; iCnt<N; iCnt++) {
		for (size_t iMember=0; iMember < (sizeof(COUNTS)/sizeof(unsigned int)); iMember++) {
			CountPtr[0][iMember] += CountPtr[iCnt][iMember];
		}
	}
}

#endif /*_STATISTICS_H*/

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
