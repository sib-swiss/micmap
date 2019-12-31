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
#ifndef _FAST_ALIGN_H
#define _FAST_ALIGN_H
#include <stdio.h> 
#include "ALIGNkonst.h"
#include "GTLkonst.h"
#include "MATCHkonst.h"
// #include "../align.h"
#define kGAPPED_ALIGN_GENOME_LENGTH 256
#define STATE_MEMORY_SIZE  (kGAPPED_ALIGN_GENOME_LENGTH+kTagSize+1) 

#ifndef __MIC__
#include <xmmintrin.h>
#endif

typedef struct GlobalData {
	void * restrict MatrixStorage;
	void * restrict WorkStorage;
	unsigned char * restrict Genome;
	unsigned char * restrict Tag;
	char * restrict AlignmentSequence;
	unsigned char * States; 
	size_t AlignmentSequenceLength;
	size_t MatrixLD;
	size_t GenomeLength;
	size_t TagLength;
	int score;
	int AlignmentRange[2];
	int ProfileRange[2];
	unsigned int MismatchCount;
	int SoftClipBoundary;
	unsigned int SoftClipMismatchCount;
	unsigned int RequiredMMSize;
	unsigned int SoftClipEvictedMMSize;
	unsigned int RequiredCigSize;
	unsigned int SoftClipEvictedCigSize;
	int BestAlignmentIndex;
	unsigned int CigarLength;		/* Needed for gtl caller for the time being */
	unsigned int MismatchLength;
	enum VectorPosition LastState;
	unsigned char revNmer;
	unsigned char StateSequence[STATE_MEMORY_SIZE];
	unsigned char Cigar[kMAXcigarSize];
	unsigned char Mismatch[kMAXmismatchSize];  
	unsigned char EndOnMatchOrDeletion; 
} GlobalData_t; 

typedef int (*GlobalFctPtr)(GlobalData_t * const restrict);

typedef struct ComputeFunctions {
	GlobalFctPtr createMatrix;
	GlobalFctPtr getBestAlignment;
	GlobalFctPtr getStateSequence;
	GlobalFctPtr computeCigarAndMismatch;
	GlobalFctPtr allocStorage;
	GlobalFctPtr freeStorage;
} ComputeFunctions_t;

#ifdef __MIC__
extern ComputeFunctions_t mic_std;
extern ComputeFunctions_t mic_vector;
extern ComputeFunctions_t mic_systolic;
extern ComputeFunctions_t mic_systolic_simplified;
#else
extern ComputeFunctions_t cpu_std;
ComputeFunctions_t cpu_std_border;
#endif

typedef int StoredIntegerFormat;
typedef union ScoreTuple {
	StoredIntegerFormat To[4] __attribute__((aligned(ALIGNMENT)));
	StoredIntegerFormat From[4] __attribute__((aligned(ALIGNMENT)));
#if !defined(__MIC__)
	__m128i Vector;
#endif
} ScoreTuple;

#if !defined(__MIC__)
typedef union sIOP {
	struct {
		int Match;
		int Insertion;
	} Element;
	__m64 mm;
} sIOP;
typedef union lScores {
  int Element[4];
	__m128i Vector;
} lScores;
#else
typedef struct {
	int dumm_deletion;
	int dummy_extra;
	int Match;
	int Insertion;
} sIOP;
typedef struct {
  int Element[4];
} lScores;
#endif

/* 
 *WARNING: Remember we inverse profile and sequence compared to pftools 
 */
struct Profile {
	/* Enlarge genome range */
	unsigned int ExtraLeftGenomeShift;
	/* Starting range for allowed jump in */
	unsigned int AllowedJumpIn;
	/* Ending range for jump out */
	unsigned int AllowedJumpOut;
	/* Match */
	int M;
	/* Mismatch */
	int m;
	int mB; // in allowed starting region
	int mE; // in allowed ending region
	/* Deletion */
	int D;
	/* Insertion */
	int I;  
	/* Transitions */
	ScoreTuple Match;
	ScoreTuple Insertion;
	ScoreTuple Deletion;
};

extern struct Profile prf; 

void dumpAlignerScores(const char * FileName);
int loadAlignerScores(const char * FileName);
void printCigar(const GlobalData_t * const restrict data);
void printMismatch(const GlobalData_t * const restrict data);
void printCigarStr(const unsigned char * const restrict data, FILE* );
void printMismatchStr(const unsigned char * const restrict data, FILE*);
int EncodeInternalAlignment(GlobalData_t * const restrict data, const unsigned int SoftClipMismatchToRemove );

static inline int __attribute__((always_inline))  ToScore(const int value)
{
	unsigned int res = ((unsigned int) value) >> SCORE_SHIFT;
	
	if (value < 0) {
		res |=  LEFT_NEGATIVE_SCORE_MASK;
	}
		return (int) res;
}

#endif
/* vim: tabstop=2 shiftwidth=2
 */
