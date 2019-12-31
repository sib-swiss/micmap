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
#ifndef _ANALYZE_TAGS_H
#define _ANALYZE_TAGS_H
#define NIRVANA 1
#define CALLER  2

#ifndef TAG_ANALYSIS_MODE
#error "MODE IS NOT SET, SPECIFY TAG_ANALYSIS_MODE"
#endif

#ifdef DEBUG_INFO
#  define IF_DEBUG_INFO(s) s
#else 
#  define IF_DEBUG_INFO(s)  
#endif

#define MAX_ALLELE_SIZE (2*kMaxReadLen)

#define kMAXPWMPOS 1024
typedef struct POS_struct {
	unsigned int chr;
	unsigned int pos;
} POS;

typedef struct OPTIONS_struct {
	int chr;
	int noCloneFilter;
	unsigned int reportVCF;
	unsigned int minorAllelePct;
	unsigned int filterQualValue;
	unsigned int PWMcnt;
	unsigned int selectAdd;
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
} OPTIONS;

typedef struct BAMRawData {
	unsigned long Ordinal;
	unsigned char * Tag;
	unsigned char * qs;
	unsigned char * Cigar;
	unsigned int Location;
	int AlignmentRange[2];
	unsigned int TagLength;
	unsigned short int CigarLength;
 	unsigned char revNmer;
	char Header[17]; 
} BAMRawData_t ;

typedef struct Interval {
	struct Interval * prev;
	void ** volatile data;
	volatile size_t size;
	volatile size_t count;
	unsigned int GenomeStart;
	unsigned int GenomeEnd;
#if TAG_ANALYSIS_MODE == NIRVANA
	bcf1_t * VCFRecord;
#endif
	pthread_mutex_t mutex;
	sem_t semaphore;
} IntervalJob_t;

extern OPTIONS options;
extern VIRTUALCHR virtchr[kMAXCHRcnt];
extern Genome_t Genome;
extern unsigned char *genome_tbl;
extern char MinQuality;

#if TAG_ANALYSIS_MODE == NIRVANA
extern atomic_uint nVCFRecord;
extern atomic_uint nRecordsFound;
extern atomic_uint nZygoteError;
extern atomic_uint nLocationNotFound;
extern atomic_uint nAllele0NotFound;
extern atomic_uint nAllele1NotFound;
extern atomic_uint nAlignmentFailure;
extern atomic_uint nNotEnoughCoverage;
extern atomic_uint nRemovedHasFork;
extern atomic_uint nNoVariant;
extern bcf_hdr_t * hdr;
#endif

#ifdef DEBUG_INFO
extern unsigned int longuestAllowedDeletion;
#endif

extern _Bool PerformTwoPass;

#endif

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
