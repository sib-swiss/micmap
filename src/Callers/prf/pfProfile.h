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
  Oct 3, 2011 pfProfile.h
 ***************************************************************************************************
 (C) 2011 Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@isb-sib.ch)
 ***************************************************************************************************/
#ifndef _PROFILE_H
#define _PROFILE_H
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <sys/types.h>
#include <string.h>
#include <stdio.h>
#include <emmintrin.h>

#define ALPHABET                "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
#define ALPHABET_SIZE 		26
/* Extra space for :
 *	undefined letter -> index 0
 * 	deletion cost    -> ALPHABET_SIZE+1
 *	stop codon       -> ALPHABET_SIZE+2
 */
#define ALPHABET_EXTRA_LETTERS  3
#define ALPHABET_MEMORY_SIZE	(ALPHABET_SIZE+ALPHABET_EXTRA_LETTERS)
#define TRANSITION_SIZE 	9
#define PROFILE_MAX_LINE_SIZE   512

#define KDIS                    2
#define KDPM                    2
/* Maximum number of normalization modes, half is for score and the other for heuristic score */ 
#define MAXN                    8
/* Maximum number of function type in normalization (LINEAR,...) */
#define KNOR                    3
/* Maximum number of coefficients used in normalization function */
#define KNPM                    5
/* Maximum number of cutoffs */
#define MAXC                    8
/* Maximum number of alignments per sequence */
#define NALI 1024

#ifndef _BEST_IS_NEGATIVE_
/* 
 * Lowest score value used when forbidding paths in the computation 
 * Remember that we need to be able to go underneath without underflow
 */
# ifndef __USE_32BIT_INTEGER__
#  define NLOW                    -16383
#  define STORED_INT_MIN	SHRT_MIN
#  define STORED_INT_MAX	SHRT_MAX
# else
// Division by 4 is required to prevent the pfplot tagging to underflow
#  define NLOW             -536870912/4 
#  define STORED_INT_MIN        INT_MIN
#  define STORED_INT_MAX	INT_MAX
# endif
# define NLOW_16                 -16383
#else
/* 
 * Highest score value used when forbidding paths in the computation 
 * Remember that we need to be able to go over without underflow
 */
# ifndef __USE_32BIT_INTEGER__
#  define NLOW                     16383
#  define STORED_INT_MIN	SHRT_MIN
#  define STORED_INT_MAX	SHRT_MAX
# else
// Division by 4 is required to prevent the pfplot tagging to underflow
#  define NLOW              536870912/4 
#  define STORED_INT_MIN	INT_MIN
#  define STORED_INT_MAX	INT_MAX
# endif
# define NLOW_16                 16383
#endif

/************************** Insertion Matrices ************************
 * There are 3 different matrices for:
 *    - alphabet
 *    - boundaries
 *    - transitions
 */
#define _UNKNOWN 0
#define _D       ALPHABET_SIZE+1
#define _STOP    ALPHABET_SIZE+2
/*
 ************************************************************************************************
 *                                    BOUNDARIES                                                *
 ************************************************************************************************
 */
/* Price to pay to enter the profile alignment in First Sequence Protein */
#define _B0  0
/* Price to pay to enter the profile alignment elsewhere in the sequence */
#define _B1  1
/* Price to pay to leave the profile alignment in Last Sequence Protein */
#define _E0  2
/* Price to pay to leave the profile alignment elsewhere in the sequence */
#define _E1  3

/* Price to pay to enter the profile alignment used in both First Sequence Protein and Extra*/
/* From entrance to Match */
#define _BM  4
/* From entrance to Insertion */
#define _BI  5
/* From entrance to Deletion */
#define _BD  6
/* From entrance to ??? THIS SEEMS TO BE UNUSED*/
#define _BE  7

/* Extensions to pay to end the profile alignment used in both LastSequence and Extra*/
/* From Match to exit */
#define _ME  8
/* From Insertion to exit */
#define _IE  9
/* From Deletion to exit */
#define _DE 10

#define INSERTION_BOUNDARIES_SIZE 11

/*
 ************************************************************************************************
 *                                   TRANSITIONS                                                *
 ************************************************************************************************
 */
 /* TABLE OF TRANSITIONS
  *
  * CODE    | x < NDIP1  | NDIP1 <= i < NDIP2  | NDIP2 <= i  | part of
  * =============================================================================================
  * _XM     | B1   + BM  |       B1 + BM       |  NLOW + BM  |
  * _XI     | B1   + BI  |       B1 + BI       |  NLOW + BI  | Insertion score EXTRA
  * _XD     | B1   + BD  |       B1 + BD       |  NLOW + BD  |
  *         |            |                     |             |
  * _MX     | NLOW + ME  |       E1 + ME       |   E1  + ME  | Insertion score MATCH
  * _IX     | NLOW + IE  |       E1 + IE       |   E1  + IE  | Insertion score INSERTION
  * _DX     | NLOW + DE  |       E1 + DE       |   E1  + DE  | Insertion score DELETION
  *---------+------------+---------------------+-------------+-----------------------------------
  * _YM     | B0 + BM    |       B0 + BM       |  NLOW + BM  |
  * _YI     | B0 + BI    |       B0 + BI       |  NLOW + BM  | FIRST SEQUENCE PROTEIN
  * _YD     | B0 + BD    |       B0 + BD       |  NLOW + BM  |
  *---------+------------+---------------------+-------------+-----------------------------------
  * _MY     | NLOW + ME  |       E0 + ME       |   E0  + ME  |
  * _IY     | NLOW + IE  |       E0 + IE       |   E0  + IE  | LAST SEQUENCE PROTEIN
  * _DY     | NLOW + DE  |       E0 + DE       |   E0  + DE  |
  *---------+------------+---------------------+-------------+-----------------------------------
  */

////////////////////////////// FOR SSE2 to WORK //////////////////////////////////////////////
// WARNING: MATCH=2,INSERTION=3,DELETION=0,DUMMY=1
#define CPP_MATCH 2
#define CPP_INSERTION 3
#define CPP_DELETION 0
#define CPP_EXTRA 1
enum VectorPosition {
  /* Positions within both 4-tuple and array of 4-tuple */
  MATCH = CPP_MATCH,
  INSERTION = CPP_INSERTION,
  DELETION = CPP_DELETION,
  EXTRA = CPP_EXTRA,
  /* Position of empty space within 4-tuple */
  DUMMY = CPP_EXTRA,
  /* Positions of transition from ? to ? within array of 4-tuples */
  /* MATCH VECTOR */
  _MM = 4*MATCH+MATCH,     _MI = 4*MATCH+INSERTION,     _MD = 4*MATCH+DELETION,     _MX    = 4*MATCH+EXTRA,
  /* INSERTION VECTOR */
  _IM = 4*INSERTION+MATCH, _II = 4*INSERTION+INSERTION, _ID = 4*INSERTION+DELETION, _IX    = 4*INSERTION+EXTRA,
  /* DELETION VECTOR */
  _DM = 4*DELETION+MATCH,  _DI = 4*DELETION+INSERTION,  _DD = 4*DELETION+DELETION,  _DX    = 4*DELETION+EXTRA,
  /* EXTRA VECTOR */
  _XM = 4*EXTRA+MATCH,     _XI = 4*EXTRA+INSERTION,     _XD = 4*EXTRA+DELETION,     _DUMMY = 4*EXTRA+EXTRA,
  /* Overall Size of transtion structure */
  INSERTION_TRANSITIONS_SIZE = 16 /* 4 times sizeof(ScoreTuple) */
};

enum ProfileType { PF_MATRIX=0b1, PF_PATTERN=0b10};

#define TRANSITION_INDEX_FROM_TO(A,B) (4*(A)+(B))

////////////////////////////// FOR SSE to WORK //////////////////////////////////////////////
// WARNING: MATCH and INSERTION should be next to each other (in this ordering)
//          and either in the upper part of xmm or the lower part.
//          Do no mix them and correct the storing function according to your above choice.
//
// StoreMatchInsertion( __m64 * const _address, const __m128 _reg)
#define StoreMatchInsertion(_address, _reg) { _mm_storeh_pi(_address, _reg); }
// For vertical xali1 we have to define a store match and deletion, hence shuffle will be
// necessary to accomodate the above.
// StoreMatchDeletion( __m64 * const _address, const __m128i _reg)
#define StoreMatchDeletion(_address, _reg) { _mm_storel_epi64((__m128i*) _address, _mm_shuffle_epi32 (_reg, _MM_SHUFFLE(DUMMY,DUMMY,DELETION,MATCH))); }

// COVERAGE FUNCTION TO STORE SCORE AND LENGTH in 64 bits
union ScoreLength {
  struct {
    int Score;
    int Length;
  } Element;
  __m64 mm;
};
//StoreScoreLength( union ScoreLength * const _address, __m128i _score, __m128i _length)
#define StoreScoreLength(_address, _score, _length) {  _mm_storeh_pi(&_address->mm, (__m128) _mm_unpacklo_epi32(_score, _length);}

/*
 ************************************************************************************************
 *                              FUNCTION POINTERS PART 1                                        *
 ************************************************************************************************
 */
typedef int (*NormalizedToRawFunctionPtr)(const float, const float * const restrict, const float, const size_t);
typedef float (*RawToNormalizedFunctionPtr)(const int, const float * const restrict, const float, const size_t);

/*
 ************************************************************************************************
 *                                    DEFINITIONS                                               *
 ************************************************************************************************
 */
#ifdef __USE_32BIT_INTEGER__
typedef int StoredIntegerFormat;
typedef union ScoreTuple {
  StoredIntegerFormat To[4];
  StoredIntegerFormat From[4];
  __m128i vector;
} ScoreTuple;
#else
typedef short int StoredIntegerFormat;
typedef union ScoreTuple {
  StoredIntegerFormat To[4];
  StoredIntegerFormat From[4];
  __m64 vector;
} ScoreTuple;
#endif

typedef union {
  ScoreTuple From[4];
  StoredIntegerFormat Element[INSERTION_TRANSITIONS_SIZE];
} TransitionScores;

typedef struct InsertionScores {
  StoredIntegerFormat * Alphabet;
  StoredIntegerFormat * Boundaries;
  TransitionScores * Transitions;
  ScoreTuple * FirstSequenceProtein;
  ScoreTuple * LastSequenceProtein;
} struct_InsertionScores;

/************************** Match Matrix   ************************/
typedef struct MatchScores {
  StoredIntegerFormat * alphabet;
} struct_MatchScores;

typedef struct Disjoint {
  char CDIS[KDIS][16];
  int  JDIP[KDIS];
  int  NDIP[KDPM];
  int  MDIS;
} SDisjoint;

typedef struct NormalizationItem {
  float RNOP[KNPM];		// Normalization coefficients
  char  CNTX[32];		// description as read in profile TEXT
  int   MNOR;			// index of the function type within CNOR (MNOR<KNOR)
  int   NNOR;			// Mode as read in profile MODE
  int   NNPR;			// Priority as read in profile PRIORITY
} SNormalizationItem;

typedef struct Normalization {
  SNormalizationItem Values[MAXN];
  char  CNOR[KNOR][16];		// description of the mode function
  int   JNOP[KNOR];             // number of coefficients 
  int   JNOR;
} SNormalization;

enum NormalizationMode { LINEAR=0, GLE_ZSCORE=1, GLE_ZSCAVE=2 };
extern const char NormalizationModeName[3][16];

typedef struct Average {
  float * Weights;
  size_t size;
} SAverage;

typedef struct CutOffItem {
  char  CCUT[32];		// description as read in profile TEXT
  float RCUT[MAXN];		// Mode normalized cutoff as read in profile N_SCORES
  int   MCUT[MAXN];		// Modes as read in profile MODE
  int   MCLE;			// Level as read in profile LEVEL
  int   ICUT;			// Filter cutoff SCORE
  unsigned int HCUT;		// Heuristic cutoff H_SCORE
  int   JCNM;			// Number of modes defined, used dimension of mode arrays
} SCutOffItem;

typedef struct CutOff {
  SCutOffItem Values[MAXC];
  int JCUT;
} SCutOff;

struct Profile {
  char Identification[64] __attribute__((aligned(16)));
  char AC_Number[64] __attribute__((aligned(16)));
  char Date[128] __attribute__((aligned(16)));
  char Description[256] __attribute__((aligned(16)));
  unsigned char Alphabet_Mapping[ALPHABET_SIZE+2] __attribute__((aligned(16)));
  char CABC[ALPHABET_SIZE+2];
  enum ProfileType Type;
  char * Pattern;
  char * Sequence;
  size_t Length;
  size_t Alphabet_Length;

  union Scores {
    struct SInsertion {
      StoredIntegerFormat * Alphabet;
      StoredIntegerFormat * Boundaries;
      TransitionScores * Transitions;
      ScoreTuple * FirstSequenceProtein;
      ScoreTuple * LastSequenceProtein;
      size_t AlignStep;
      struct_MatchScores _dummy;
    } Insertion;
    struct SMatch {
      struct_InsertionScores _dummy;
      size_t AlignStep;
      StoredIntegerFormat * Alphabet;
    } Match;
  } Scores;

  SNormalization NormalizationData;
  SDisjoint DisjointData;
  SCutOff CutOffData;

  NormalizedToRawFunctionPtr NormalizedToRaw;
  RawToNormalizedFunctionPtr RawToNormalized;
  SAverage Average;

  short int LevelIndex;						/* WARNING: This is not the real level but the index corresponding within the cutoff array */
  short int ModeIndex;          	/* WARNING: This is not the real mode but the index corresponding within the normalization array */
  short int HeuristicModeIndex;		/* WARNING: This is not the real mode but the index corresponding within the normalization array */
  short int ModeIndexWithinLevel;	/* WARNING: This is not the real Mode but the index within the cutoff array of modes. */

  unsigned int HeuristicCutOff;
  int CutOff;
  float NormalizedCutOff;
  enum NormalizationMode NormalizationType;
  float * restrict NormalizationCoefs;

  _Bool isCircular;
	_Bool CompleteCycleOnly;
  _Bool isReversed;
	_Bool ReverseSequence;

  struct Profile * next;
  struct Profile * previous;
};

#ifndef PFSEQ
#define PFSEQ
typedef struct PFSequence {
  unsigned char * ProfileIndex;
  size_t Length;
} PFSequence;
#endif

struct UniProtMatch {
  char (*UniqueIdentifier)[16];
  char (*EntryName)[16];
  char *State;
  char (*DB)[2];
  int *DBIndex;
  unsigned int * HeuristicScore;
  int * FilterScore;
  size_t Count;
  unsigned int True_positive;
  unsigned int Unknown;
  unsigned int Partial;
  unsigned int False_negative;
  unsigned int False_posistive;
};

/******************* Function choice enumerator *****************/
enum Version { SSE2=0, SSE41=1};

/************** Sequence generated from profile options *********/
enum GeneratedSequenceOptions { 
  GENERATE_MATCH=1, 
  GENERATE_INSERTION=2,
  GENERATE_DELETION=4
};

/********************** Complementar mode  **********************/
enum ComplementaryMode { NONE=0, DNA, IUPAC };

/*
 ************************************************************************************************
 *                              PROFILE FUNCTION DECLARATIONS                                   *
 ************************************************************************************************
 */
int PrepareExtraTable(struct Profile * const prf);
struct Profile * ReverseProfile(const struct Profile * const restrict inprf);
int ReadProfile(const char * const restrict FileName, struct Profile * const prf,
                const _Bool SetExtraTable, const _Bool CompleteCycleOnly);
void FreeProfile(struct Profile * const prf, const _Bool IsPointer);

#include "pfProfileInline.h"

#endif

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
