/*
 * ------------------------------------------------------------------------------------------------------------------------
 *
 *                                * micmap *
 *      Mapping of short reads in fastq format onto a reference
 *
 *
 *  Copyright (C) SIB  - Swiss Institute of Bioinformatics,                2015-2019 Nicolas Guex, Thierry Schuepbach and Christian Iseli
 *  Copyright (C) UNIL - University of Lausanne, Switzerland               2019-2020 Nicolas Guex and Christian Iseli
 *  Copyright (C) EPFL - Ecole Polytechnique Fédérale de Lausanne, Switzerland  2020 Christian Iseli
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
 *      Contacts:   Nicolas.Guex@unil.ch and Christian.Iseli@epfl.ch
 *      Repository: https://github.com/sib-swiss/micmap
 *
 * ------------------------------------------------------------------------------------------------------------------------
 */
#ifndef _GTL_H
#define _GTL_H
#include <stdio.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdatomic.h>
#include "GTLkonst.h"

typedef	struct PAIRCHRPOSDATA_struct  PAIRCHRPOSDATA;
struct PAIRCHRPOSDATA_struct
{
	unsigned int	chr;     /* (2x16 bits)  if bit #15 is set, it means reverse comp */
	unsigned int	tag1pos; /* mapped position */
	unsigned int	tag2pos; /* mapped position */
};


typedef	struct PAIRPOSDATA_struct  PAIRPOSDATA;
struct PAIRPOSDATA_struct
{
	unsigned int	tag1pos;    /* mapped position */
	unsigned int	tag2offset; /* offset to start of second tag (also encodes the two reversecomp flags) */
	/* [note could be stored on 16bits in majority of cases - thus need to make two distinct data structures for maximal compression ] */
};

// These all need to stay together for sorting purposes
typedef struct PAIRDATA_struct {
	PAIRCHRPOSDATA pcpd[kMAXnbMatches];
	unsigned long ordinal;
	unsigned int nbMatches;
	unsigned int len1;
	unsigned int len2;
	unsigned int hlen1;
	unsigned int hlen2;
	unsigned int ml[kMAXnbMatches];
	unsigned int cl[kMAXnbMatches];
	unsigned int pl;
	char hdr1[kMaxHdrSize];
	char hdr2[kMaxHdrSize];
	char qs1[kMaxReadLen];
	char qs2[kMaxReadLen];
	//unsigned char MMstring[2*((2*kMaxEncodedMismatch)+1)]; this was for things other than unmapped
	unsigned char MMstring[kMAXnbMatches][4*kMaxReadLen+2]; // worse case is pair of tag, with pos + nt + terminators (when dealing with unmapped)
	unsigned char cigar[kMAXnbMatches][2*kMAXcigarSize];

	unsigned char packed4nt[kMaxReadLen/2];
} PAIRDATA;

typedef struct TLBHEADER_struct {
	char magic[4];
	unsigned int blockLength;
	unsigned int chr;
	unsigned int minPos;
	unsigned int maxPos;
	unsigned char headerCS;
	unsigned char flags_23_16;
	unsigned char flags_15_8;
	unsigned char flags_7_0;
} TLBHEADER;

typedef struct BASICBLOCK_struct {
	unsigned int compressor;
	unsigned int origLen;
	unsigned int packedLen;
	unsigned char *data;
} BASICBLOCK;
// we have
// max readLen  = kBlockMaxCount * 2 * sizeof(unsigned int) (= kBlockMaxCount *    8)
// max pairPos  = kBlockMaxCount * sizeof(PAIRPOSDATA)      (= kBlockMaxCount *    8) * kMAXnbMatches
// max pairPos  = kBlockMaxCount * sizeof(PAIRCHRPOSDATA)   (= kBlockMaxCount *   12) * kMAXnbMatches
// max unmapped = kBlockMaxCount * kMaxReadLen/2            (= kBlockMaxCount *  128)
// max mismatch = kBlockMaxCount * 4*kMaxReadLen+2          (= kBlockMaxCount * 1026) * kMAXnbMatches
// max cigar    = kBlockMaxCount * 2*kMAXcigarSize          (= kBlockMaxCount *  129) * kMAXnbMatches
// max header   = kBlockMaxCount * (kMaxHdrSize * 2 + 10)   (= kBlockMaxCount *  266)
// max qual     = kBlockMaxCount * kMaxReadLen * 2          (= kBlockMaxCount *  512)
// total                                                       kBlockMaxCount * 3264

typedef struct TLBDATA_struct {
	TLBHEADER header;
	unsigned int fletcher32CS;
	unsigned int cnt;
	unsigned int lengthInfo;
	unsigned int nbMatches;
	unsigned int hdrSize;												// Not populated when reading
	unsigned int qualSize;											// Not populated when reading
	unsigned int mmSize[kMAXnbMatches];
	unsigned int cigarSize[kMAXnbMatches];
	unsigned int p4Size;
	PAIRDATA *pd;
	PAIRDATA **pdp; // sorted index into pd
	BASICBLOCK readLen;
	BASICBLOCK pairPos[kMAXnbMatches];
	BASICBLOCK mismatch[kMAXnbMatches];
	BASICBLOCK cigarBlock[kMAXnbMatches];
	BASICBLOCK unmapped;
	BASICBLOCK readHdr;
	BASICBLOCK readQual;
	// the decoded stuff
	unsigned int *len;
	PAIRPOSDATA *ppd[kMAXnbMatches];
	PAIRCHRPOSDATA *pcpd[kMAXnbMatches];
	unsigned char *p4nt; // stay at buffer head
	unsigned char *hdrs; // stay at buffer head
	unsigned char *qual; // stay at buffer head
	char *hdr; // can move around
	unsigned char *qs; // can move around
	unsigned char *diff[kMAXnbMatches];
	unsigned char *cigar[kMAXnbMatches];
	unsigned char *packed4nt;
	unsigned char *diskBuffer;
} TLBDATA;


typedef struct flush_queue
{
	unsigned int head;
	unsigned int tail;
	unsigned int cnt;
	unsigned int nbFreeBlocks;
	TLBDATA *block[kFLUSHQUEUEBUFFERSIZE];
	TLBDATA **freeBlock; // pool of free blocks; put it here under control of the mutex
	pthread_mutex_t mutex;
	pthread_cond_t not_empty_cond;
	pthread_cond_t not_full_cond;
} FlushQueue;

typedef struct GTLBUFFERS_struct {
// should be enough room to get all threads busy, queue full, and all buffers being filled (+1)
  unsigned int OUTPUTFILEScnt;
  unsigned int chrMax;
  unsigned int totalBlocks;
  TLBDATA *allBlock; // all the available data blocks
  TLBDATA **perfectBlock; // blocks set aside for perfect matches
  TLBDATA **mismatchNBlock; // blocks set aside for mismatches with only Ns
  TLBDATA **mismatchBlock; // blocks set aside for mismatches
  TLBDATA **alignBlock; // blocks set aside for mismatches + cigar
  TLBDATA **halfmapBlock; // blocks set aside for halfmap (only one tag of the pair mapped)
  TLBDATA **unmappedBlock; // blocks set aside for unmapped
  TLBDATA **chimeraBlock; // blocks set aside for chimera (both tags of the pair map on different chromosomes)
  TLBDATA **MperfectBlock; // blocks set aside for multiple perfect matches (max 2 at the moment - kMAXnbMatches)
  TLBDATA **MmismatchNBlock; // blocks set aside for multiple matches with mismatches with only Ns
  TLBDATA **MmismatchBlock; // blocks set aside for multiple matches with mismatches
  TLBDATA **MalignBlock; // blocks set aside for multiple matches with mismatches + cigar
  FlushQueue fq; // ng moved out of globals
  char	seedfn[512];
	atomic_ulong * WrittenCnt;
} GTLBUFFERS;

// Need to define this here, since we need the size of the header
#define kMaxBlockHeadSize (sizeof(TLBHEADER) + (4 + ((4 + 3 * kMAXnbMatches) * 3)) * sizeof(unsigned int))
// look into allocBlock of GTL.c for a breakdown of the 3264 number
#define kMaxBlockDataSize (kMaxBlockHeadSize + kBlockMaxCount * 3264)
// we love round numbers...  - go for a disk block of 512 bytes...
#define kMaxBlockLength (kMaxBlockDataSize + (~((kMaxBlockDataSize - 1) & 0x1ff) & 0x1ff))

// compression stuff
typedef struct 
{
	enum 
	{
		CMP_RAW = 0,
		CMP_ZIP = 1,
		CMP_BZ2 = 2,
		CMP_REDT6 = 3, // Range encoder from Daniel and Thierry (legacy version 6)
		CMP_REDT = 4,  // Range encoder from Daniel and Thierry (current version 7)
		CMP_LZ4 = 5   // Very fast compressor
	}
	def,    // default compression
	scoredef;       // default compression for score blocks
}
compress_t;

// global UMI settings
extern int gHasUMI; // used for sorting in comparePAIRDATA()

// global compression settings
extern compress_t gCompress;

typedef const struct 
{
	int id;
	const char* name;
	const char* description;
}
compress_dict_t;

extern compress_dict_t compress_dict[];


//=================================================================================================

int compress_find(const char* name);
const char* compress_name(int id);
void initBuffers(GTLBUFFERS * const, const unsigned int, const unsigned int,
                 const unsigned int * const restrict, const int);
void flushAllBuffers(GTLBUFFERS *G);
TLBDATA * get_block(FlushQueue *fq);
int flushBlock(FlushQueue *fq,TLBDATA **);
unsigned int countPat1InUnmappedBlock(void *, unsigned int, TLBHEADER *, unsigned int *);
void countPat2InUnmappedBlock(void *, unsigned int, TLBHEADER *, unsigned int *, unsigned int, unsigned int *);
int readBlock(int, TLBHEADER*, TLBDATA*);
int readDecompressBlock(int, TLBHEADER *, TLBDATA *, int);
int preadDecompressBlock(int, off_t, TLBHEADER*, TLBDATA *, int);
int decompressBlock(void *, const unsigned int, TLBHEADER*, TLBDATA*, int);
int decodeBlock(void *, unsigned int, TLBHEADER *, PAIRDATA *, TLBDATA *);
int readDecodeBlock(int, TLBHEADER *, PAIRDATA *, TLBDATA *);
void allocBlock(TLBDATA *, const _Bool);
void zeroBlock(TLBDATA *);
void freeBlock(TLBDATA *);
unsigned char simple8bitCS(unsigned char const *, size_t);
int comparePAIRCHRPOSDATA(const void *, const void *);
int comparePAIRDATAptr(const void *, const void *);
int comparePAIRDATA(const void *, const void *);
int comparePAIRUDATAptr(const void *, const void *);
int comparePAIRUDATA(const void *, const void *);
int compareUNMAPDATAptr(const void *, const void *);
int compareUNMAPDATA(const void *, const void *);
unsigned char * packBlock(TLBDATA *, TLBDATA *);
unsigned char * compressBlock(TLBDATA *);
void report_compress_stat(FILE* const);

#endif /* _GTL_H */
/* vim: tabstop=2 shiftwidth=2
 */
