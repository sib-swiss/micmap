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
#ifndef _GTLcaller_common_H
#define _GTLcaller_common_H
/*

	common definitions for various callers

	(c) N.Guex and C.Iseli 2017

*/
#include <stdio.h>
#include "fastalign.h"

//=================================================================================================
#define kMAXNBINDELS 128
#define kMAXINDELLEN 128
#define kMAXPAIRLEN (1024+10)
#define kMaxLineBuf 1024
#define kMAXPWMPOS 1024
#define kMAXDELTA 15
#define kMINCOVWIGGLEREPORT 4
// #define kNOT_ACGT 0xff

//=================================================================================================
typedef struct POS_struct {
	unsigned int chr;
	unsigned int pos;
} POS;

typedef struct OPTIONS_struct {
	int chr;
	int noCloneFilter;
	int additionalFilter;
	unsigned int reportVCF;
	unsigned int minorAllelePct;
	unsigned int filterQualValue;
	unsigned int PWMcnt;
	unsigned int selectAdd;
	unsigned int maxPairLen;
	unsigned int Realign;
	POS PWMpos[kMAXPWMPOS];
	char outbase[512];
	char selectRegions[512];
	char wiggleBed[512];
	char expressionBase[512];
	char gencodeGFF3[512];
	char tromerTranscripts[512];
	char sampleName1[64];
	char sampleName2[64];
	char refName[64];
	FILE *ref1;
	FILE *ref2;
} OPTIONS;

typedef struct _decodedSingle {
	unsigned long ordinal;
	unsigned int genomepos;
	int taglen1;
	int reverseTAG1;
} decodedSingle_t, *decodedSingle_p_t;

typedef struct _decodedPair {
	unsigned long ordinal;
	unsigned int genomepos;
	int taglen1;
	int taglen2;
	int reverseTAG1;
	int reverseTAG2;
	int delta;
} decodedPair_t, *decodedPair_p_t;

typedef struct _bt_node_t {
	unsigned long refNo:61;
	unsigned long var:3;
} bt_node_t, *bt_node_p_t;

typedef struct _variant {
	void *ref;
	unsigned short ACGTP[4];
	unsigned short ACGTM[4];
	struct InternalDel {
		unsigned short cnt;
		unsigned short len;
	} del[kMAXNBINDELS];
	struct InternalIns {
		unsigned short cnt;
		char nt[38]; // FIXME - need something better, and check is wrong in the code
	} ins[kMAXNBINDELS];
} gtl_variant_t, *gtl_variant_p_t;

typedef struct _coverage {
	gtl_variant_p_t var;
	unsigned short cntP;
	unsigned short cntM;
	unsigned short cov[3]; // 0 -> perfect; 1 -> mismatches; 2-> gapped
	unsigned short flags;
} coverage_t, *coverage_p_t;

typedef struct _refRead {
	unsigned long offset:8;
	unsigned long len:8;
	unsigned long reverse:1;
	unsigned long second:1;
	unsigned long ordinal:46;
} refRead_t, *refRead_p_t;

typedef struct _reference {
	refRead_p_t rr;
	unsigned int max;
	unsigned int cnt;
} reference_t, *reference_p_t;

typedef struct _covStat {
	unsigned int max[4];
	unsigned int mean[4];
	unsigned int median[4];
} covStat_t, *covStat_p_t;

//=================================================================================================
extern VIRTUALCHR virtchr[kMAXCHRcnt];
extern unsigned char *genome_tbl; // contain genome in 16 virtual chromosomes.
extern unsigned char gNT2bits[128];
extern unsigned int nThreads;

//=================================================================================================

int loadGenome(const char * const restrict , OPTIONS *);
void reportRef(FILE *, gtl_variant_p_t, char, reference_p_t, unsigned int, OPTIONS *);
void loadSelected(coverage_p_t, OPTIONS *);
void computeCov(coverage_p_t, OPTIONS *, coverage_p_t, unsigned int *, unsigned int *, unsigned int *, unsigned int *);
int decompress(const char * const restrict , OPTIONS *, coverage_p_t, reference_p_t);
int decompress_realign(const char * const restrict , OPTIONS *, coverage_p_t, reference_p_t);
int decompress_PairedThreaded(const char * const restrict fn, OPTIONS *opt, coverage_p_t cov, reference_p_t ref);

void 
ExtractStateSequence(decodedPair_p_t dpp, unsigned char * qual,
                     const unsigned char * restrict diff, const unsigned char * restrict cigar,
                     GlobalData_t * const restrict GD);

#endif /* _GTLcaller_common_H */
/* vim: tabstop=2 shiftwidth=2
 */
