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
  Apr 29, 2014 sse41_inline_fcts.h
 *******************************************************
 (C) 2011-14 Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@isb-sib.ch)
 *******************************************************/
#ifndef SSE41_INLINE_FCTS_H_
#define SSE41_INLINE_FCTS_H_
#ifdef __SSE4_1__
#ifndef __USE_32BIT_INTEGER__
#define LoadStoredIntegerVector(address) _mm_cvtepi16_epi32(_mm_loadl_epi64((__m128i*) &((address)->vector)));
#else  /*__USE_32BIT_INTEGER__*/
#define LoadStoredIntegerVector(address) _mm_load_si128(&((address)->vector));
#endif /*__USE_32BIT_INTEGER__*/
#else  /*__SSE4_1__*/
#error "sse41_inline_fcts.h should not be included in non sse 4.1 files"
#endif /*__SSE4_1__*/

#endif

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
