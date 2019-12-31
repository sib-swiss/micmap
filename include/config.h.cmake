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
#ifndef _CONFIG_H
#define _CONFIG_H

/*
 * Hardware specific
 */
#define MAP_HUGE_SHIFT	26
#define MAP_HUGE_2MB    (21 << MAP_HUGE_SHIFT)
#define MAP_HUGE_1GB    (30 << MAP_HUGE_SHIFT)
#define HUGETLBFS_MAGIC 0x958458f6
/*
 * COMPILER SPECIFIC CONFIGURATION 
 */
#ifndef _LARGEFILE64_SOURCE
#define _LARGEFILE64_SOURCE
#endif

// GCC
#ifdef __GNUC__
#define align aligned
#define __ALWAYS_INLINE __attribute__((always_inline))
#else
#define __ALWAYS_INLINE __attribute__((always_inline, gnu_inline))
#endif

#ifndef __INTEL_COMPILER
#define __assume_aligned(a,b) __builtin_assume_aligned(a,b)
#endif

#define HAS_AVX512_INT2MASK
#ifndef __MIC__
#ifndef HAS_AVX512_INT2MASK
#define _mm512_int2mask(x) (x)
#endif
#endif

/*
 * CONSTANT TO TEXT 
 */
#define TEXT(x) TEXT2(x)
#define TEXT2(x) #x

/*
 * AVAILABLE INSTRUCTION SET
 */
#cmakedefine HANDLE_SSE
#cmakedefine HANDLE_SSE2
#cmakedefine HANDLE_SSE41
#cmakedefine HANDLE_AVX
#cmakedefine HANDLE_AVX2
#cmakedefine HANDLE_AVX512F
#cmakedefine HANDLE_AVX512CD
#cmakedefine HANDLE_AVX512BW

#endif /* _CONFIG_H */

/* vim: tabstop=2 shiftwidth=2 filetype=c
 */
/* ------------------------------------------------------------------------------------ */
