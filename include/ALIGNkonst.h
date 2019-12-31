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
//---------------------------------------------------------------
// GLOBAL DEFINITIONS
//---------------------------------------------------------------
#define kMaxGenomeChunkSize (1024-kTagSize)	/* same as above crashes... ?  */
#define MIC_ALIGNER_LISTENING_PORT 2049
//(SCIF_PORT_RSVD + 3)
#ifdef __MIC__
#define ALIGNMENT 64
#else
#define ALIGNMENT 16
#endif

//#define STORE_FUNCTION(where,what) _mm512_store_epi32(where, what)
#ifdef __MIC__
#define STORE_FUNCTION(where,what) _mm512_storenrngo_ps(where, _mm512_castsi512_ps(what))
#else
#define STORE_FUNCTION(where,what) _mm512_stream_si512((__m512i*)where, what)
#endif
//#define PREFETCHL1(address) 
#define PREFETCHL1(address) _mm_prefetch((char*)(address), _MM_HINT_T0);
#define EVICTL1(address)
//#define EVICTL1(address) _mm_clevict(address, _MM_HINT_T1);

/*
  * INTEGER FORMAT  -------------------------------------
  *
  * |xxxxxxxxxxxxxxxxxxxxxxxxxxxxx|xx|x
  *                              |  | |
  *                              |  | + Mismatch or true match
  *                              |  +-- State
  *                              +----- Score
  */
#define SCORE_SHIFT               3
#define STATE_SHIFT               1
#define MISMATCH_SHIFT            0

/* 
 * CONSTANTS USED IN ALIGNMENT
 */
// Division by 8 is required to prevent tagging to underflow
// #  define NLOW             (((unsigned int) -536870912) >> SCORE_SHIFT )
#ifndef USE_SHORT_INT
#  define NLOW            -536870912/8
#  define STORED_INT_MIN  INT_MIN
#  define STORED_INT_MAX	INT_MAX
#  define SCORE_MASK                0xFFFFFFF8
#  define STATE_MASK                0x00000006
#  define MISMATCH_MASK             0x00000001
#  define LEFT_NEGATIVE_SCORE_MASK  0xE0000000
#  define CLEAR_MASK                0xFFFFFFF8
#else
#  define NLOW            (-((1<<12)>>3))
// #  define STORED_INT_MIN  INT_MIN
// #  define STORED_INT_MAX	INT_MAX
#  define SCORE_MASK                0xFFF8
#  define STATE_MASK                0x0006
#  define MISMATCH_MASK             0x0001
#  define LEFT_NEGATIVE_SCORE_MASK  0xE000
#  define CLEAR_MASK                0xFFF8
#endif

#define _I		 0
#define _IM		 0 
#define _II		 0 
#define _ID		 NLOW
#define _IX		 0

#define _D		-1
#define _DM		 0
#define _DI		 NLOW
#define _DD		 0
#define _DX		 0

#define _m		-1
#define _M		 3
#define _MM		 1
#define _MI		-15
#define _MD		-3
#define _MX		 0 

#define _BM		 0
#define _BD		-2


# define TO_SCORE(x) ((int) ((x < 0) ? (((unsigned int) x)>>SCORE_SHIFT) | LEFT_NEGATIVE_SCORE_MASK : ((unsigned int) x)>>SCORE_SHIFT))
//---------------------------------------------------------------
// ENUMERATION DEFINITIONS
//---------------------------------------------------------------
// This is the position within the SSE register tuple
// Pay attention that Storing functions needs to be aligned on that.
enum VectorPosition {
  /* Positions within both 4-tuple and array of 4-tuple */
  MATCH=2,
  INSERTION=3,
  DELETION=0,
  EXTRA=1
};

// This is the mask used to define the priority when performing comparisons
enum StatePriority {
  PRIORITY_MATCH     = 2,
  PRIORITY_INSERTION = 0,
  PRIORITY_DELETION  = 1,
  PRIORITY_EXTRA     = 3
};

//---------------------------------------------------------------
// STRUCTURE DEFINITIONS
//---------------------------------------------------------------

//---------------------------------------------------------------
// GLOBAL PRIVATE VARIABLES
//---------------------------------------------------------------
/* 
 *WARNING: Remember we inverse profile and sequence compared to pftools 
 */
static const int RowInit[4] __attribute__((aligned(ALIGNMENT))) = {
  ( _BD << SCORE_SHIFT) | (PRIORITY_EXTRA << STATE_SHIFT),
  (NLOW << SCORE_SHIFT) | (PRIORITY_EXTRA << STATE_SHIFT),
  ( _BM << SCORE_SHIFT) | (PRIORITY_EXTRA << STATE_SHIFT),
  (NLOW << SCORE_SHIFT) | (PRIORITY_EXTRA << STATE_SHIFT) 
};
static const int Masks[4] = { 0xF, 0xF0, 0xF00, 0xF000 };
static const int DeletionIndex = DELETION;
static const char N = 'N';
// <- according to state priority, but due to reverse compare to pftools I <-> D
static const char StateLetter[4] = "DIMX";
static const unsigned int StateChangeAdd[4] = { 1, 1, 0, 0};

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
