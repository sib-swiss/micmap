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
/*
 *
 *	The Genomic Transport Layer block routines
 *
 *	(c) N.Guex C.Iseli I.Topolsky D.Zerzion T.Schuepbach 2015
 *
 */

//=================================================================================================
#define _GNU_SOURCE
#include "config.h"
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>

#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <pthread.h>

#include <semaphore.h>

// basic compression
#include <zlib.h>
#include <bzlib.h>
#include "Compresser/lz4frame.h"

// keep timing stats
#include <sys/time.h>

#include "GTL.h"
#include "Compresser/rangecod.h"

#define FROMTO_SIZE 44

//=================================================================================================
//GTLBUFFERS G;

// global compression settings
compress_t gCompress =
{
	.def = CMP_BZ2,
	.scoredef = CMP_REDT
};

compress_dict_t compress_dict[] =
{
	{ CMP_RAW, "raw", "do not compress" },
	{ CMP_ZIP, "zip", "deflate from libz - average" },
	{ CMP_BZ2, "bz2", "bzip2" },
	{ CMP_REDT6, "re6", "range encoder (legacy version 6)" },
	{ CMP_REDT, "re", "range encoder (current version 7)" },
	{ CMP_LZ4, "lz4", "lz4 - lightning fast" },
	{ -1, NULL, NULL }
};
typedef struct struct_compress_stat_t {
	int compressor;
	unsigned long raw_bytes;
	unsigned long packed_bytes;
	struct timeval tv;
} compress_stat_t;
static compress_stat_t len_stat_c, len_stat_d;
static compress_stat_t pairs_stat_c, pairs_stat_d;
static compress_stat_t mismatch_stat_c, mismatch_stat_d;
static compress_stat_t cigar_stat_c, cigar_stat_d;
static compress_stat_t unmapped_stat_c, unmapped_stat_d;
static compress_stat_t header_stat_c, header_stat_d;
static compress_stat_t qs_stat_c, qs_stat_d;

//=================================================================================================

//FlushQueue fq;

//---------------------------------------------------------------
const size_t BlockMemorySize = kBlockMaxCount * sizeof(PAIRDATA) + kBlockMaxCount * sizeof(PAIRDATA *) + kBlockMaxCount * 8\
                             + kMAXnbMatches*kBlockMaxCount*(12+770+129) + kBlockMaxCount*(96+266+512) + kMaxBlockLength;
														 
														 
														 
int compress_find(const char* name)
{
	unsigned i;
	for(i = 0; compress_dict[i].name != NULL; i++)
		if (0 == strcmp(name, compress_dict[i].name))
			break;
	return compress_dict[i].id;
}
//---------------------------------------------------------------

const char* compress_name(int id)
{
	unsigned i;
	for(i = 0; compress_dict[i].name != NULL; i++)
		if (id == compress_dict[i].id) 
			return compress_dict[i].name;
	return "error!";
}
//---------------------------------------------------------------

void
allocBlock(TLBDATA *tdp, const _Bool WithDiskBuffer)
{
	tdp->pd = (PAIRDATA *) malloc(kBlockMaxCount * sizeof(PAIRDATA)); assert(tdp->pd);
	tdp->pdp = (PAIRDATA **) malloc(kBlockMaxCount * sizeof(PAIRDATA *)); assert(tdp->pdp);
	tdp->readLen.data = malloc(kBlockMaxCount * 8); assert(tdp->readLen.data);
	for (unsigned int i = 0; i < kMAXnbMatches; i++)
	{
		tdp->pairPos[i].data = malloc(kBlockMaxCount * 12); assert(tdp->pairPos[i].data);
		tdp->mismatch[i].data = malloc(kBlockMaxCount * 770); assert(tdp->mismatch[i].data);
		tdp->cigarBlock[i].data = malloc(kBlockMaxCount * 129); assert(tdp->cigarBlock[i].data);
	}
	tdp->unmapped.data = malloc(kBlockMaxCount * 96); assert(tdp->unmapped.data);
	tdp->readHdr.data = malloc(kBlockMaxCount * 266); assert(tdp->readHdr.data);
	tdp->readQual.data = malloc(kBlockMaxCount * 512); assert(tdp->readQual.data);
	if (WithDiskBuffer) {
		tdp->diskBuffer = malloc(kMaxBlockLength); assert(tdp->diskBuffer);
	}
	else {
		tdp->diskBuffer = NULL;
	}
}
//---------------------------------------------------------------

void
zeroBlock(TLBDATA *tdp)
{
	memset(&tdp->header,0,sizeof(TLBHEADER));
	tdp->header.minPos = 0xffffffffU;
	tdp->fletcher32CS = 0;
	tdp->cnt = 0;
	tdp->lengthInfo = 0;
	tdp->nbMatches = 0;
	tdp->hdrSize = 0;
	tdp->qualSize = 0;
	tdp->p4Size = 0;
	memset(&tdp->readLen,0,3 * sizeof(unsigned int));
	memset(&tdp->unmapped,0,3 * sizeof(unsigned int));
	memset(&tdp->readHdr,0,3 * sizeof(unsigned int));
	memset(&tdp->readQual,0,3 * sizeof(unsigned int));
	tdp->len = NULL;
	tdp->p4nt = NULL;
	tdp->hdrs = NULL;
	tdp->qual = NULL;
	tdp->hdr = NULL;
	tdp->qs = NULL;
	tdp->packed4nt = NULL;
	for (unsigned int i = 0; i < kMAXnbMatches; i++)
	{
		tdp->mmSize[i] = 0;
		tdp->cigarSize[i] = 0;
		memset(&tdp->pairPos[i],0,3 * sizeof(unsigned int));
		memset(&tdp->mismatch[i],0,3 * sizeof(unsigned int));
		memset(&tdp->cigarBlock[i],0,3 * sizeof(unsigned int));
		tdp->ppd[i] = NULL;
		tdp->pcpd[i] = NULL;
		tdp->diff[i] = NULL;
		tdp->cigar[i] = NULL;
	}
}
//---------------------------------------------------------------

void
initBuffers(GTLBUFFERS * const G, const unsigned int cnt, const unsigned int chrMax,
						const unsigned int * const restrict chrPartOffset, const int singleEnd)
{
	unsigned int i;
	G->OUTPUTFILEScnt = cnt;
	G->chrMax = chrMax;
	// we need 8 times cnt for p,N,m,g (normal + multiple matches) + chrMax for halfmapped and chimera + 16 for unmapped + 8 for chimera + buffer reserve
	G->totalBlocks = G->OUTPUTFILEScnt * 8 + chrMax * 2 + 16 + 8 + kFLUSHQUEUEBUFFERSIZE + kFLUSHTHREADCOUNT;
	G->allBlock = (TLBDATA *) calloc(G->totalBlocks, sizeof(TLBDATA));
	G->fq.freeBlock = (TLBDATA **) calloc(G->totalBlocks, sizeof(TLBDATA *));
	for (i = 0; i < G->totalBlocks; i++)
	{
		G->allBlock[i].header.magic[0] = 'G';
		G->allBlock[i].header.magic[1] = 'T';
		G->allBlock[i].header.magic[2] = 'L';
		G->allBlock[i].header.minPos = 0xffffffffU;
		allocBlock(G->allBlock + i, false);
	}
	printf("InitBuffers allocated %lu bytes RAM, %u blocks of %lu, kMaxBlockLength is %lu\n", G->totalBlocks*BlockMemorySize,G->totalBlocks, BlockMemorySize, kMaxBlockLength); 
	// we assume index 0 will be reserved for unmapped
	unsigned int cb = 0;
	G->perfectBlock = malloc(G->OUTPUTFILEScnt * sizeof(TLBDATA *));
	for (i = 0; i < G->OUTPUTFILEScnt; i++)
	{
		G->perfectBlock[i] = G->allBlock + cb++;
		if (!singleEnd)
			G->perfectBlock[i]->header.flags_7_0 |= kTLBHflagPairedReads;
		G->perfectBlock[i]->header.flags_7_0 |= kTLBHflagPairPosBlock;
	}
	G->mismatchNBlock = malloc(G->OUTPUTFILEScnt * sizeof(TLBDATA *));
	for (i = 0; i < G->OUTPUTFILEScnt; i++)
	{
		G->mismatchNBlock[i] = G->allBlock + cb++;
		if (!singleEnd)
			G->mismatchNBlock[i]->header.flags_7_0 |= kTLBHflagPairedReads;
		G->mismatchNBlock[i]->header.flags_7_0 |= kTLBHflagPairPosBlock;
		G->mismatchNBlock[i]->header.flags_7_0 |= kTLBHflagNMismatchBlock;
	}
	G->mismatchBlock = malloc(G->OUTPUTFILEScnt * sizeof(TLBDATA *));
	for (i = 0; i < G->OUTPUTFILEScnt; i++)
	{
		G->mismatchBlock[i] = G->allBlock + cb++;
		if (!singleEnd)
			G->mismatchBlock[i]->header.flags_7_0 |= kTLBHflagPairedReads;
		G->mismatchBlock[i]->header.flags_7_0 |= kTLBHflagPairPosBlock;
		G->mismatchBlock[i]->header.flags_7_0 |= kTLBHflagMismatchBlock;
	}
	G->alignBlock = malloc(G->OUTPUTFILEScnt * sizeof(TLBDATA *));
	for (i = 0; i < G->OUTPUTFILEScnt; i++)
	{
		G->alignBlock[i] = G->allBlock + cb++;
		if (!singleEnd)
			G->alignBlock[i]->header.flags_7_0 |= kTLBHflagPairedReads;
		G->alignBlock[i]->header.flags_7_0 |= kTLBHflagPairPosBlock;
		G->alignBlock[i]->header.flags_7_0 |= kTLBHflagMismatchBlock;
		G->alignBlock[i]->header.flags_7_0 |= kTLBHflagCigarBlock;
	}
	G->MperfectBlock = malloc(G->OUTPUTFILEScnt * sizeof(TLBDATA *));
	for (i = 0; i < G->OUTPUTFILEScnt; i++)
	{
		G->MperfectBlock[i] = G->allBlock + cb++;
		if (!singleEnd)
			G->MperfectBlock[i]->header.flags_7_0 |= kTLBHflagPairedReads;
		G->MperfectBlock[i]->header.flags_7_0 |= kTLBHflagPairPosBlock;
		G->MperfectBlock[i]->header.flags_15_8 |= kTLBHflagMultipleMatches;
	}
	G->MmismatchNBlock = malloc(G->OUTPUTFILEScnt * sizeof(TLBDATA *));
	for (i = 0; i < G->OUTPUTFILEScnt; i++)
	{
		G->MmismatchNBlock[i] = G->allBlock + cb++;
		if (!singleEnd)
			G->MmismatchNBlock[i]->header.flags_7_0 |= kTLBHflagPairedReads;
		G->MmismatchNBlock[i]->header.flags_7_0 |= kTLBHflagPairPosBlock;
		G->MmismatchNBlock[i]->header.flags_7_0 |= kTLBHflagNMismatchBlock;
		G->MmismatchNBlock[i]->header.flags_15_8 |= kTLBHflagMultipleMatches;
	}
	G->MmismatchBlock = malloc(G->OUTPUTFILEScnt * sizeof(TLBDATA *));
	for (i = 0; i < G->OUTPUTFILEScnt; i++)
	{
		G->MmismatchBlock[i] = G->allBlock + cb++;
		if (!singleEnd)
			G->MmismatchBlock[i]->header.flags_7_0 |= kTLBHflagPairedReads;
		G->MmismatchBlock[i]->header.flags_7_0 |= kTLBHflagPairPosBlock;
		G->MmismatchBlock[i]->header.flags_7_0 |= kTLBHflagMismatchBlock;
		G->MmismatchBlock[i]->header.flags_15_8 |= kTLBHflagMultipleMatches;
	}
	G->MalignBlock = malloc(G->OUTPUTFILEScnt * sizeof(TLBDATA *));
	for (i = 0; i < G->OUTPUTFILEScnt; i++)
	{
		G->MalignBlock[i] = G->allBlock + cb++;
		if (!singleEnd)
			G->MalignBlock[i]->header.flags_7_0 |= kTLBHflagPairedReads;
		G->MalignBlock[i]->header.flags_7_0 |= kTLBHflagPairPosBlock;
		G->MalignBlock[i]->header.flags_7_0 |= kTLBHflagMismatchBlock;
		G->MalignBlock[i]->header.flags_7_0 |= kTLBHflagCigarBlock;
		G->MalignBlock[i]->header.flags_15_8 |= kTLBHflagMultipleMatches;
	}
	unsigned int chr = 1;
	for (i = 0; i < G->OUTPUTFILEScnt; i++)
	{
		// FIXME - this seems a bit shaky... Already bit me when there is only 1
		// chromosome so added chr < chrMax part, but not sure whether this is
		// enough
		if (chr < chrMax && chrPartOffset[chr+1] == i)
			chr += 1;
		G->perfectBlock[i]->header.chr = chr;
		G->mismatchNBlock[i]->header.chr = chr;
		G->mismatchBlock[i]->header.chr = chr;
		G->alignBlock[i]->header.chr = chr;
		G->MperfectBlock[i]->header.chr = chr;
		G->MmismatchNBlock[i]->header.chr = chr;
		G->MmismatchBlock[i]->header.chr = chr;
		G->MalignBlock[i]->header.chr = chr;
	}
	G->halfmapBlock = malloc(chrMax * sizeof(TLBDATA *));
	for (i = 0; i < chrMax; i++)
	{
		G->halfmapBlock[i] = G->allBlock + cb++;
		if (!singleEnd)
			G->halfmapBlock[i]->header.flags_7_0 |= kTLBHflagPairedReads;
		G->halfmapBlock[i]->header.flags_15_8 |= kTLBHflagHalfmapBlock;
		G->halfmapBlock[i]->header.flags_7_0 |= kTLBHflagPairChrPosBlock;
		G->halfmapBlock[i]->header.flags_7_0 |= kTLBHflagMismatchBlock;
		G->halfmapBlock[i]->header.chr = ((i + 1) << 16) | kHalfMapChrID;
	}
	if (singleEnd)
	{
		G->unmappedBlock = malloc(4 * sizeof(TLBDATA *));
		for (i = 0; i < 4; i++)
		{
			G->unmappedBlock[i] = G->allBlock + cb++;
			G->unmappedBlock[i]->header.flags_15_8 |= kTLBHflagUnmappedBlock;
			G->unmappedBlock[i]->header.flags_7_0 |= kTLBHflagMismatchBlock;
			G->unmappedBlock[i]->header.chr = (i << 16) | kUnmappedChrID;
		}
	}
	else
	{
		G->unmappedBlock = malloc(16 * sizeof(TLBDATA *));
		for (i = 0; i < 16; i++)
		{
			G->unmappedBlock[i] = G->allBlock + cb++;
			G->unmappedBlock[i]->header.flags_7_0 |= kTLBHflagPairedReads;
			G->unmappedBlock[i]->header.flags_15_8 |= kTLBHflagUnmappedBlock;
			G->unmappedBlock[i]->header.flags_7_0 |= kTLBHflagMismatchBlock;
			G->unmappedBlock[i]->header.chr = (i << 16) | kUnmappedChrID;
		}
	}
	G->chimeraBlock = malloc((chrMax + 8) * sizeof(TLBDATA *));
	for (i = 0; i < (chrMax + 8); i++)
	{
		G->chimeraBlock[i] = G->allBlock + cb++;
		if (!singleEnd)
			G->chimeraBlock[i]->header.flags_7_0 |= kTLBHflagPairedReads;
		G->chimeraBlock[i]->header.flags_15_8 |= kTLBHflagChimeraBlock;
		G->chimeraBlock[i]->header.flags_7_0 |= kTLBHflagPairChrPosBlock;
		G->chimeraBlock[i]->header.flags_7_0 |= kTLBHflagMismatchBlock;
		G->chimeraBlock[i]->header.chr = ((i + 1) << 16) | kChimeraChrID;
	}
	if (cb >= G->totalBlocks)
	{
		fprintf(stderr, "HWHAP: not enough buffers... %d vs %d\n", cb, G->totalBlocks);
		exit(1);
	}
	i=0;
	while (cb < G->totalBlocks)
	{
		G->fq.freeBlock[i++] = G->allBlock + cb++;
		G->fq.nbFreeBlocks += 1;
	}
}
//---------------------------------------------------------------

/* snarfed from wikipedia...  */
static uint32_t
fletcher32(uint16_t const *data, size_t words)
{
	uint32_t sum1 = 0xffff, sum2 = 0xffff;

	while (words) {
		unsigned tlen = words > 359 ? 359 : words;
		words -= tlen;
		do {
			sum2 += sum1 += *data++;
		} while (--tlen);
		sum1 = (sum1 & 0xffff) + (sum1 >> 16);
		sum2 = (sum2 & 0xffff) + (sum2 >> 16);
	}
	/* Second reduction step to reduce sums to 16 bits */
	sum1 = (sum1 & 0xffff) + (sum1 >> 16);
	sum2 = (sum2 & 0xffff) + (sum2 >> 16);
	return sum2 << 16 | sum1;
}
//---------------------------------------------------------------

unsigned char
simple8bitCS(unsigned char const *data, size_t len)
{
	unsigned int acc = 0xff;
	unsigned int i;
	for (i = 0; i < len; i++)
	{
		acc += data[i];
	}
	while (acc>>8)
		acc = (acc & 0xff) + (acc >> 8);
	return ~acc;
}
//---------------------------------------------------------------

int
comparePAIRCHRPOSDATA(const void *a, const void *b)
{
	PAIRCHRPOSDATA *pcpda = ((PAIRCHRPOSDATA *) a);
	PAIRCHRPOSDATA *pcpdb = ((PAIRCHRPOSDATA *) b);
	unsigned short ac1 = (pcpda->chr >> 16) & 0x7fff;
	unsigned short ac2 = pcpda->chr & 0x7fff;
	// try to deal reasonably with halfmaped and single-end
	// WARNING : this relies on the fact that tag1pos and tag2pos cannot be 0 for real matches
	if (pcpda->tag1pos == 0)
		ac1 = ac2;
	if (pcpda->tag2pos == 0)
		ac2 = ac1;
	if (ac2 < ac1)
	{
		ac2 = ac1;
		ac1 = pcpda->chr & 0x7fff;
	}
	unsigned short bc1 = (pcpdb->chr >> 16) & 0x7fff;
	unsigned short bc2 = pcpdb->chr & 0x7fff;
	if (pcpdb->tag1pos == 0)
		bc1 = bc2;
	if (pcpdb->tag2pos == 0)
		bc2 = bc1;
	if (bc2 < bc1)
	{
		bc2 = bc1;
		bc1 = pcpdb->chr & 0x7fff;
	}
	if (ac1 < bc1)
		return -1;
	if (ac1 > bc1)
		return 1;
	if (ac2 < bc2)
		return -1;
	if (ac2 > bc2)
		return 1;
	unsigned int mina,maxa;
	if (pcpda->tag1pos == 0)
	{
		mina = pcpda->tag2pos;
		maxa = pcpda->tag2pos;
	}
	else
	{
		if (pcpda->tag2pos == 0)
		{
			mina = pcpda->tag1pos;
			maxa = pcpda->tag1pos;
		}
		else
		{
			if (pcpda->tag1pos <= pcpda->tag2pos)
			{
				mina = pcpda->tag1pos;
				maxa = pcpda->tag2pos;
			}
			else
			{
				mina = pcpda->tag2pos;
				maxa = pcpda->tag1pos;
			}
		}
	}
	unsigned int minb,maxb;
	if (pcpdb->tag1pos == 0)
	{
		minb = pcpdb->tag2pos;
		maxb = pcpdb->tag2pos;
	}
	else
	{
		if (pcpdb->tag2pos == 0)
		{
			minb = pcpdb->tag1pos;
			maxb = pcpdb->tag1pos;
		}
		else
		{
			if (pcpdb->tag1pos <= pcpdb->tag2pos)
			{
				minb = pcpdb->tag1pos;
				maxb = pcpdb->tag2pos;
			}
			else
			{
				minb = pcpdb->tag2pos;
				maxb = pcpdb->tag1pos;
			}
		}
	}
	if (mina < minb)
		return -1;
	if (mina > minb)
		return 1;
	if (maxa < maxb)
		return -1;
	if (maxa > maxb)
		return 1;
	return 0;
}
//---------------------------------------------------------------

int
comparePAIRDATA(const void *a, const void *b)
{
	PAIRDATA *pda = ((PAIRDATA *) a);
	PAIRDATA *pdb = ((PAIRDATA *) b);
	int res = comparePAIRCHRPOSDATA(pda->pcpd,pdb->pcpd);
	if (res != 0)
		return res;
	if (pda->ml[0] < pdb->ml[0])
		return -1;
	if (pda->ml[0] > pdb->ml[0])
		return 1;
	if (pda->cl[0] < pdb->cl[0])
		return -1;
	if (pda->cl[0] > pdb->cl[0])
		return 1;
	return strcmp(pda->hdr1,pdb->hdr1);
}
//---------------------------------------------------------------

int
comparePAIRDATAptr(const void *a, const void *b)
{
	PAIRDATA *pda = *((PAIRDATA **) a);
	PAIRDATA *pdb = *((PAIRDATA **) b);
	return comparePAIRDATA(pda,pdb);
}
//---------------------------------------------------------------

int
compareUNMAPDATA(const void *a, const void *b)
{
	PAIRDATA *pda = ((PAIRDATA *) a);
	PAIRDATA *pdb = ((PAIRDATA *) b);
	if (pda->pl > 0 && pdb->pl > 0)
	{
		unsigned int pl = pda->pl;
		if (pdb->pl < pl)
			pl = pdb->pl;
		int res = memcmp(pda->packed4nt,pdb->packed4nt,pl);
		if (res != 0)
			return res;
	}
	if (pda->pl < pdb->pl)
		return -1;
	if (pda->pl > pdb->pl)
		return 1;
	return strcmp(pda->hdr1,pdb->hdr1);
}
//---------------------------------------------------------------

int
comparePAIRUDATA(const void *a, const void *b)
{
	PAIRDATA *pda = ((PAIRDATA *) a);
	PAIRDATA *pdb = ((PAIRDATA *) b);
	int res = comparePAIRDATA(pda,pdb);
	if (res != 0)
		return res;
	// let's see if we have unmapped that we can compare
	return compareUNMAPDATA(a,b);
}
//---------------------------------------------------------------

int
comparePAIRUDATAptr(const void *a, const void *b)
{
	PAIRDATA *pda = *((PAIRDATA **) a);
	PAIRDATA *pdb = *((PAIRDATA **) b);
	return comparePAIRUDATA(pda,pdb);
}
//---------------------------------------------------------------

int
compareUNMAPDATAptr(const void *a, const void *b)
{
	PAIRDATA *pda = *((PAIRDATA **) a);
	PAIRDATA *pdb = *((PAIRDATA **) b);
	return compareUNMAPDATA(pda,pdb);
}
//---------------------------------------------------------------

static void
simpleCompress(BASICBLOCK *bb, int compressor, compress_stat_t *stat, unsigned char *packed, unsigned int pMaxSize)
{
	size_t uSize = bb->origLen;
	size_t pSize = 0;
	struct timeval ts, te, td;

	gettimeofday(&ts,NULL);
	switch(compressor)
	{
		case CMP_ZIP:
		{
			pSize = compressBound(uSize);
			assert(pSize <= pMaxSize);
			stat->compressor = CMP_ZIP;
			int res = compress(packed,&pSize,bb->data,uSize);
			if (res == Z_OK && pSize < uSize)
			{
				// packing worked and saved some space
				bb->compressor = CMP_ZIP;
				bb->packedLen = pSize;
			}
			else
			{
				memcpy(packed,bb->data,uSize);
				bb->packedLen = bb->origLen;
			}
		}
		break;
		case CMP_BZ2:
		case CMP_REDT6: // we get here if the real RE could not be applied for some reason
		case CMP_REDT: // we get here if the real RE could not be applied for some reason
		{
			unsigned int dSize = uSize + uSize / 100 + 600; // from manual
			assert(dSize <= pMaxSize);
			stat->compressor = CMP_BZ2;
			int res = BZ2_bzBuffToBuffCompress((char*)packed,&dSize,(char*)bb->data,uSize,9,0,0); // Workfactor 0 = use default (30, according to doc)
			if (res == BZ_OK && dSize < uSize)
			{
				// packing worked and saved some space
				bb->compressor = CMP_BZ2;
				bb->packedLen = dSize;
			}
			else
			{
				memcpy(packed,bb->data,uSize);
				bb->packedLen = bb->origLen;
			}
		}
		break;
		case CMP_LZ4:
		{
			// NOTE LZ4 has two possible use. Raw compressed bit stream (as kernel does) or frames (as lz4 command does)
			// we use the frame format
			LZ4F_preferences_t pref  = 
			{
				.frameInfo = 
				{
					.blockSizeID =	max4MB,
					.blockMode =	blockLinked, // otherwise possible to decompress blocks in threads for even more performance
					.contentChecksumFlag =	contentChecksumEnabled,
					.frameType =	LZ4F_frame
				},
				.compressionLevel = 0,
					// lh4 :	compresses bloody fast, but has not so good performance (dict search is not exhaustive)
					// lh4hc :	a bit slower but better performance (dictionnary is the best achievable with lz4 format)
					// usually lz4 is best for in-ram targets and high speed SSD (like PCI-Express or NXE based)
					// HDD spinning media can be slower than lh4hc and thus its worth
				.autoFlush = 0
			};
			pSize = LZ4F_compressFrameBound(uSize, &pref);
			assert(pSize <= pMaxSize);
			stat->compressor = CMP_LZ4;
			size_t res = LZ4F_compressFrame(packed,pSize,bb->data,uSize,&pref);
			if (! LZ4F_isError(res) &&  res < uSize)
			{
				// packing worked and saved some space
				bb->compressor = CMP_LZ4;
				bb->packedLen = res;
			}
			else
			{
				memcpy(packed,bb->data,uSize);
				bb->packedLen = bb->origLen;
			}
		}
		break;
		case CMP_RAW:
			// just chain to worst case handler
			memcpy(packed,bb->data,uSize);
			bb->packedLen = bb->origLen;
		break;
		default:
			fprintf(stderr,"Bug: unimplemented compressor... %d ???\n", compressor);
			exit(1);
	}
	gettimeofday(&te,NULL);
	timersub(&te, &ts, &td);
	timeradd(&td, &(stat->tv), &(stat->tv));
	stat->raw_bytes += bb->origLen;
	stat->packed_bytes += bb->packedLen;
}
//---------------------------------------------------------------

extern inline unsigned int __attribute__((always_inline))
TailNullFrom(const unsigned char * const restrict LineBuffer, const int length)
{
	for (int i=(length-1); i>=0; i--) {
		if (LineBuffer[i] != '#' )
			return 1 + i;
	}
	return 0;
}
//---------------------------------------------------------------

static void
FromToFrequenciesPosition(const unsigned char * const restrict streamin, freq * const restrict counters, const unsigned int nlines, const freq length)
{
	/* first zero the counters */
	memset(counters, 0, FROMTO_SIZE*FROMTO_SIZE*length*sizeof(freq));

	/* then count the number of occurances of each byte */
	unsigned int iline = 0;
	const freq limit = length - 1;
	const unsigned char * restrict LineBuffer = streamin;

	int previous = 0;

	while (iline < nlines) {
		freq upto = TailNullFrom(LineBuffer, (int) length);
		freq * restrict colcounters = counters;

		if (upto < limit) {
			for (freq i=0; i<upto; ++i) {
				const int ch = (int) (LineBuffer[i] - 33 + 1);
				colcounters[previous*FROMTO_SIZE + ch] += 1;
				previous = ch;
				colcounters += FROMTO_SIZE*FROMTO_SIZE;
			}
			colcounters[previous*FROMTO_SIZE + 0 ] += 1;
			previous = 0;
		}
		else {
			for (freq i=0; i<length; ++i) {
				int ch = ((int) LineBuffer[i]) - 33 + 1;
				ch = ( ch < 0 ) ? 0 : ch;
				colcounters[previous*FROMTO_SIZE + ch] += 1;
				previous = ch;
				colcounters += FROMTO_SIZE*FROMTO_SIZE;
			}
		}
		iline++;
		LineBuffer  += length;
	}
}
//---------------------------------------------------------------

#ifndef __APPLE__
static int
EncodeQS7(char * * streamout, const unsigned char * restrict const streamin, const unsigned int nlines, const unsigned int nchars)
{
	rangecoder rc;

	/* Open a stream for writing */
	size_t streamSize;
	FILE * out = open_memstream(streamout, &streamSize);
	if (out == NULL) {
		return -1;
	}

	/* initialize the range coder, first byte 0, no header */
	start_encoding(&rc,0,0, out);

	encode_freq(&rc,1,1,2); /* a stupid way to code a bit */

	freq * const counters = (freq*) malloc(FROMTO_SIZE*FROMTO_SIZE*nchars*sizeof(freq));
	if (counters == NULL) {
		fprintf(stderr,"Unable to allocate memory\n");
		return -2;
	}

	FromToFrequenciesPosition(streamin, counters, nlines, nchars);

	/*
	 * write the statistics of position dependent from on a per
	 * position 2D array of FROMTO_SIZE * FROMTO_SIZE with range
	 * encoder adjusted to maximum value
	 */

	/* compute the maximum value for the FromTo table counter */
	freq maxValue = 0;
	for(freq i=0; i<(FROMTO_SIZE*FROMTO_SIZE)*nchars; i++) maxValue = (maxValue > counters[i]) ? maxValue : counters[i];

	/* compute the maximum value for the run length encoder */
	freq max = 0;
	{
		freq i =0;
		while (i<(FROMTO_SIZE*FROMTO_SIZE)*nchars) {
			const freq Value = counters[i];
			if (Value == 0) {
				freq count = 1;
				freq j = i + 1;
				while ((j < (FROMTO_SIZE*FROMTO_SIZE)*nchars) && (counters[j] == Value)) {
					count++;
					j++;
				}
				i += count;
				max = (max > count) ? max : count;
			}
			else {
				i++;
			}
		}
	}

	/* Cant use putchar or other since we are after start of the rangecoder */
	/* as you can see the rangecoder doesn't care where probabilities come */
	/* from, it uses a flat distribution of 0..0xffff in encode_short. */
	encode_short(&rc, (short) nchars);
	encode_shift(&rc, (freq)1, nlines, (freq)24);
	encode_shift(&rc, (freq)1, max, (freq) 24);
	max += maxValue;
	encode_shift(&rc, (freq)1, max, (freq) 24);

	{
		freq i =0;
		while (i<(FROMTO_SIZE*FROMTO_SIZE)*nchars) {
			freq Value = counters[i];
			if (Value == 0 ) {
				freq count = 1;
				freq j = i + 1;
				while ((j < (FROMTO_SIZE*FROMTO_SIZE)*nchars) && (counters[j] == Value)) {
					count++;
					j++;
				}
				Value = maxValue + count;
				i += count;
			}
			else {
				i++;
			}
			encode_freq(&rc, (freq)1, Value, max+1);
		}
	}

	/* compute CDF */
	{
		freq * restrict colcounters = counters;
		for (freq k=0; k<FROMTO_SIZE*nchars; k++) {
			freq sum = 0;
			for (freq j=0; j<FROMTO_SIZE; j++) {
				const freq tmp = colcounters[j];
				colcounters[j] = sum;
				sum += tmp;
			}
			colcounters += FROMTO_SIZE;
		}
	}


	/* output the encoded symbols */
	unsigned int iline = 0;
	const unsigned char * restrict LineBuffer = streamin;
	const freq limit = (freq) (nchars-1);

	int previous = 0;
	while(iline < nlines) {
		const freq upto = TailNullFrom(LineBuffer, (int) nchars);
		freq * restrict colcounters = counters;
		for (freq i=0; i<((upto < limit) ? upto : nchars); ++i) {
			int ch = (int) (LineBuffer[i]-33+1);
			ch = (ch < 0) ? 0 : ch;
			const freq count = colcounters[previous*FROMTO_SIZE+ch+1]-colcounters[previous*FROMTO_SIZE+ch];
			if (count == 0) {
				fprintf(stderr,"%.*s\nCharacter '%c'=%i is not in alphabet @ %lu\n->\n",(int) nchars, LineBuffer, LineBuffer[i], ch,i);
				return -3;
			}
			encode_freq(&rc,count,colcounters[previous*FROMTO_SIZE+ch],colcounters[previous*FROMTO_SIZE+FROMTO_SIZE-1]);
			for (freq k=ch+1; k<FROMTO_SIZE;k++) colcounters[previous*FROMTO_SIZE + k]--;
			colcounters += FROMTO_SIZE*FROMTO_SIZE;
			previous = ch;
		}
		if (upto < limit ) {
			int ch = (int) 0;
			const freq count = colcounters[previous*FROMTO_SIZE+ch+1]-colcounters[previous*FROMTO_SIZE+ch];
			if (count == 0) {
				fprintf(stderr,"%.*s\nCharacter '%c' is not in alphabet @ %lu\n->\n",(int) nchars,LineBuffer, ch, upto);
				return -4;
			}
			encode_freq(&rc,count,colcounters[previous*FROMTO_SIZE+ch],colcounters[previous*FROMTO_SIZE+FROMTO_SIZE-1]);
			for (freq k=ch+1; k<FROMTO_SIZE;k++) colcounters[previous*FROMTO_SIZE + k]--;
			previous = ch;
		}
		LineBuffer += nchars;
		iline++;
	}

	/* flag absence of next block by a bit */
	encode_freq(&rc,1,0,2);

	/* close the encoder */
	done_encoding(&rc);

	free(counters);

	fclose(out);
	if (streamSize != (size_t) rc.bytecount) return -9;

	return (int) rc.bytecount;
}
//---------------------------------------------------------------

static void
simpleCompressQS(BASICBLOCK *bb, const unsigned int nlines, const unsigned int nchars, compress_stat_t *stat, unsigned char *packedRes, unsigned int pMaxSize)
{
	unsigned long uSize = bb->origLen;
	unsigned char *packed = NULL;
	struct timeval ts, te, td;

	gettimeofday(&ts,NULL);
	stat->compressor = CMP_REDT;
	int res = EncodeQS7((char **)&packed,bb->data,nlines,nchars);
	if (res > 0 && res < uSize)
	{
		// packing worked and saved some space
		bb->compressor = CMP_REDT; // Daniel & Thierry's range encoder
		bb->packedLen = (unsigned int) res;
		assert(bb->packedLen <= pMaxSize);
		memcpy(packedRes,packed,bb->packedLen);
	}
	else
	{
		// packing did not work or did not save space
		// FIXME - maybe we still want to try with simpleCompress ?
		bb->packedLen = bb->origLen;
		assert(bb->packedLen <= pMaxSize);
		memcpy(packedRes,bb->data,bb->packedLen);
	}
	free(packed);
	gettimeofday(&te,NULL);
	timersub(&te, &ts, &td);
	timeradd(&td, &(stat->tv), &(stat->tv));
	stat->raw_bytes += bb->origLen;
	stat->packed_bytes += bb->packedLen;
}
//---------------------------------------------------------------
#endif

static void
report_compress_stat_1(FILE* const out, compress_stat_t *stat, char *name)
{
	if (stat->raw_bytes == 0)
		return;
	fprintf(
		out,
		"%s(%s): %ld vs %ld (%.2f%% ; CR %.2f) %ds %.2f MB/s\n",
		name,
		compress_dict[stat->compressor].name,
		stat->packed_bytes,
		stat->raw_bytes,
		(stat->raw_bytes == 0) ? 0.0 : ((double) stat->packed_bytes) * 100.0 / ((double) stat->raw_bytes),
		(stat->packed_bytes == 0) ? 0.0 : ((double) stat->raw_bytes) / ((double) stat->packed_bytes),
		(int)stat->tv.tv_sec,
		(stat->tv.tv_sec == 0) ? ((double) (stat->raw_bytes / (1024 * 1024))) : ((double) stat->raw_bytes) / ((double) (stat->tv.tv_sec * 1024 * 1024)));
}
//---------------------------------------------------------------

void
report_compress_stat(FILE * const out)
{
	report_compress_stat_1(out, &len_stat_c,"C len");
	report_compress_stat_1(out, &pairs_stat_c,"C pairs");
	report_compress_stat_1(out, &mismatch_stat_c,"C mismatch");
	report_compress_stat_1(out, &cigar_stat_c,"C cigar");
	report_compress_stat_1(out, &unmapped_stat_c,"C unmapped");
	report_compress_stat_1(out, &header_stat_c,"C header");
	report_compress_stat_1(out, &qs_stat_c,"C qs");
	report_compress_stat_1(out, &len_stat_d,"D len");
	report_compress_stat_1(out, &pairs_stat_d,"D pairs");
	report_compress_stat_1(out, &mismatch_stat_d,"D mismatch");
	report_compress_stat_1(out, &cigar_stat_d,"D cigar");
	report_compress_stat_1(out, &unmapped_stat_d,"D unmapped");
	report_compress_stat_1(out, &header_stat_d,"D header");
	report_compress_stat_1(out, &qs_stat_d,"D qs");
}
//---------------------------------------------------------------

/**
 * This function obtains a queued task from the pool and returns it.
 * If no such task is available the operation blocks.
 */
TLBDATA *
get_block(FlushQueue *fq)
{
	//printf("Requesting one task (of %d)\n",fq.cnt);
	if (pthread_mutex_lock(&fq->mutex) != 0) {
		perror("pthread_mutex_lock: ");
		return NULL;
	}
	while (fq->cnt == 0) {
		/* Block until a new task arrives. */
		//printf("Waiting for task\n");
		if (pthread_cond_wait(&fq->not_empty_cond,&fq->mutex)) {
			perror("pthread_cond_wait: ");
			if (pthread_mutex_unlock(&fq->mutex)) {
				perror("pthread_mutex_unlock: ");
			}
			return NULL;
		}
	}
	//printf("Got 1 out of %d tasks\n",fq->cnt);
	TLBDATA *block = fq->block[fq->tail++];
	if (fq->tail >= kFLUSHQUEUEBUFFERSIZE)
		fq->tail -= kFLUSHQUEUEBUFFERSIZE;
	fq->cnt -= 1;
	pthread_cond_signal(&fq->not_full_cond);
	if (fq->cnt > 0)
		pthread_cond_signal(&fq->not_empty_cond);
	if (pthread_mutex_unlock(&fq->mutex)) {
		perror("pthread_mutex_unlock: ");
		return NULL;
	}
	//printf("Returning 1 tasks\n");
	return block;
}
//---------------------------------------------------------------

/**
 * This function adds a block to be flushed in the queue and replaces it with an empty block
 * If the queue is full the operation blocks.
 */
int
flushBlock(FlushQueue *fq,TLBDATA **bp)
{
	TLBDATA *block = *bp;
#ifndef NDEBUG
	TLBHEADER keptH;
	memcpy(&keptH,&block->header,sizeof(TLBHEADER));
#endif
	//printf("entering flushBlock G.fq->cnt=%d\n",G.fq->cnt); fflush(stdout);
	if (pthread_mutex_lock(&fq->mutex) != 0) {
		perror("pthread_mutex_lock: ");
		return -1;
	}
	while (fq->cnt + 1 > kFLUSHQUEUEBUFFERSIZE) {
		/* Block until there is enough space. */
		//printf("Waiting for space fq->cnt=%d cnt=%d\n",fq->cnt,cnt); fflush(stdout);
		if (pthread_cond_wait(&fq->not_full_cond,&fq->mutex)) {
			perror("pthread_cond_wait: ");
			if (pthread_mutex_unlock(&fq->mutex)) {
				perror("pthread_mutex_unlock: ");
			}
			return -1;
		}
	}
	//printf("Got enough space for %d (%d)\n", cnt, fq->cnt); fflush(stdout);
	fq->block[fq->head++] = block;
	if (fq->head >= kFLUSHQUEUEBUFFERSIZE)
		fq->head = 0;
	fq->cnt += 1;
	if (fq->cnt < kFLUSHQUEUEBUFFERSIZE)
		pthread_cond_signal(&fq->not_full_cond);
	pthread_cond_signal(&fq->not_empty_cond);
	// replace the full block with an empty one from the free list
	if (fq->nbFreeBlocks == 0) {
		fprintf(stderr,"HWHAP: we have no more free blocks... but we should...\n");
		return -1;
	}
	*bp = fq->freeBlock[--fq->nbFreeBlocks];
	// FIXME - probably safer to keep a local copy of the header ?
	memcpy(&(*bp)->header,&block->header,sizeof(TLBHEADER));
#ifndef NDEBUG
	// haven't seen it fail yet...  keep this as a debug conditional
	if (memcmp(&keptH, &(*bp)->header, sizeof(TLBHEADER)))
	{
		fprintf(stderr,"HWHAP: I'm afraid it changed, Jim!\n");
	}
#endif
	(*bp)->header.minPos = 0xffffffffU;
	(*bp)->header.maxPos = 0;
	if (pthread_mutex_unlock(&fq->mutex)) {
		perror("pthread_mutex_unlock: ");
		return -1;
	}
	//printf("stored tasks\n"); fflush(stdout);
	return 0;
}
//---------------------------------------------------------------

void
flushAllBuffers(GTLBUFFERS *G)
{
	unsigned int i;
	unsigned long cnt = 0UL;
	for (i = 0; i < G->OUTPUTFILEScnt; i++)
	{
		if (G->perfectBlock[i]->cnt != 0) {
			flushBlock(&G->fq,G->perfectBlock + i); cnt += G->perfectBlock[i]->cnt;
		}
		if (G->mismatchNBlock[i]->cnt != 0) { 
			flushBlock(&G->fq,G->mismatchNBlock + i); cnt += G->mismatchNBlock[i]->cnt;
		}
		if (G->mismatchBlock[i]->cnt != 0) {
			flushBlock(&G->fq,G->mismatchBlock + i); cnt += G->mismatchBlock[i]->cnt;
		}
		if (G->alignBlock[i]->cnt != 0) {
			flushBlock(&G->fq,G->alignBlock + i); cnt += G->alignBlock[i]->cnt;
		}
		if (G->MperfectBlock[i]->cnt != 0) {
			flushBlock(&G->fq,G->MperfectBlock + i); cnt += G->MperfectBlock[i]->cnt;
		}
		if (G->MmismatchNBlock[i]->cnt != 0) {
			flushBlock(&G->fq,G->MmismatchNBlock + i); cnt += G->MmismatchNBlock[i]->cnt;
		}
		if (G->MmismatchBlock[i]->cnt != 0) {
			flushBlock(&G->fq,G->MmismatchBlock + i); cnt += G->MmismatchBlock[i]->cnt;
		}
		if (G->MalignBlock[i]->cnt != 0) {
			flushBlock(&G->fq,G->MalignBlock + i); cnt += G->MalignBlock[i]->cnt;
		}
	}
	for (i = 0; i < G->chrMax; i++)
		if (G->halfmapBlock[i]->cnt != 0) {
			flushBlock(&G->fq,G->halfmapBlock + i); cnt += G->halfmapBlock[i]->cnt;
		}
	unsigned int nbUnmapped = 16;
	if ((G->unmappedBlock[0]->header.flags_7_0 & kTLBHflagPairedReads) == 0)
		nbUnmapped = 4;
	for (i = 0; i < nbUnmapped; i++)
		if (G->unmappedBlock[i]->cnt != 0) {
			flushBlock(&G->fq,G->unmappedBlock + i); cnt += G->unmappedBlock[i]->cnt;
		}
	for (i = 0; i < G->chrMax + 8; i++)
		if (G->chimeraBlock[i]->cnt != 0) {
			flushBlock(&G->fq,G->chimeraBlock + i); cnt += G->chimeraBlock[i]->cnt;
		}
	// they now should all be empty.  Just enqueue enough empty blocks to reap all the threads
	for (i = 0; i < kFLUSHTHREADCOUNT;i++)
		flushBlock(&G->fq,G->perfectBlock);
	
	atomic_fetch_add(G->WrittenCnt, cnt);
}
//---------------------------------------------------------------

unsigned char *
compressBlock(TLBDATA *wb)
{
	unsigned char *dataBlock = wb->diskBuffer; assert(dataBlock);
	const int multi = (wb->header.flags_15_8 & kTLBHflagMultipleMatches) != 0;
	const size_t headerSize = sizeof(TLBHEADER) + (multi ? 4 : 3) * sizeof(unsigned int);
	// leave room for header and set it to zeroes
	memset(dataBlock,0,headerSize);
	unsigned char *wpos = dataBlock + headerSize;
	unsigned char *wmax = dataBlock + kMaxBlockLength;
	// leave room for each block triplet if present
	if (wb->readLen.origLen != 0)
	{
		wpos += 3 * sizeof(unsigned int);
	}
	if (wb->pairPos[0].origLen != 0)
		for (unsigned int i = 0; i < wb->nbMatches; i++)
		{
			wpos += 3 * sizeof(unsigned int);
		}
	if (wb->mismatch[0].origLen != 0)
		for (unsigned int i = 0; i < wb->nbMatches; i++)
		{
			wpos += 3 * sizeof(unsigned int);
		}
	if (wb->cigarBlock[0].origLen != 0)
		for (unsigned int i = 0; i < wb->nbMatches; i++)
		{
			wpos += 3 * sizeof(unsigned int);
		}
	if (wb->unmapped.origLen != 0)
	{
		wpos += 3 * sizeof(unsigned int);
	}
	wpos += 3 * sizeof(unsigned int);
	wpos += 3 * sizeof(unsigned int);
	// compress the blocks directly in place
	if ((wb->header.flags_7_0 & kTLBHflagFixedLength) == 0)
	{
		simpleCompress(&wb->readLen,gCompress.def,&len_stat_c,wpos,wmax - wpos);
		if (wb->readLen.packedLen == 0)
		{
			fprintf(stderr,"HWHAP: should have readLen but there are none...\n");
			exit(1);
		}
		wpos += wb->readLen.packedLen;
	}
	if ((wb->header.flags_7_0 & kTLBHflagPairPosBlock) != 0
			|| (wb->header.flags_7_0 & kTLBHflagPairChrPosBlock) != 0)
		for (unsigned int i = 0; i < wb->nbMatches; i++)
		{
			simpleCompress(&wb->pairPos[i],gCompress.def,&pairs_stat_c,wpos,wmax - wpos);
			if (wb->pairPos[i].packedLen == 0)
			{
				fprintf(stderr,"HWHAP: should have pairPos %u but there are none...\n",i);
				exit(1);
			}
			wpos += wb->pairPos[i].packedLen;
		}
	if ((wb->header.flags_7_0 & kTLBHflagNMismatchBlock) != 0
			|| (wb->header.flags_7_0 & kTLBHflagMismatchBlock) != 0)
		for (unsigned int i = 0; i < wb->nbMatches; i++)
		{
			simpleCompress(&wb->mismatch[i],gCompress.def,&mismatch_stat_c,wpos,wmax - wpos);
			if (wb->mismatch[i].packedLen == 0)
			{
				fprintf(stderr,"HWHAP: should have mismatches %u but there are none...\n",i);
				exit(1);
			}
			wpos += wb->mismatch[i].packedLen;
		}
	if ((wb->header.flags_7_0 & kTLBHflagCigarBlock) != 0)
		for (unsigned int i = 0; i < wb->nbMatches; i++)
		{
			simpleCompress(&wb->cigarBlock[i],gCompress.def,&cigar_stat_c,wpos,wmax - wpos);
			if (wb->cigarBlock[i].packedLen == 0)
			{
				fprintf(stderr,"HWHAP: should have cigars %u but there are none...\n",i);
				exit(1);
			}
			wpos += wb->cigarBlock[i].packedLen;
		}
	if ((wb->header.flags_15_8 & kTLBHflagUnmappedBlock) == kTLBHflagUnmappedBlock
			|| (wb->header.flags_15_8 & kTLBHflagHalfmapBlock) == kTLBHflagHalfmapBlock)
	{
		simpleCompress(&wb->unmapped,gCompress.def,&unmapped_stat_c,wpos,wmax - wpos);
		if (wb->unmapped.packedLen == 0)
		{
			fprintf(stderr,"HWHAP: should have unmapped but there are none...\n");
			exit(1);
		}
		wpos += wb->unmapped.packedLen;
	}
	simpleCompress(&wb->readHdr,gCompress.def,&header_stat_c,wpos,wmax - wpos);
	wpos += wb->readHdr.packedLen;
	// FIXME offload compression of QS to MIC0
	// move it at beginning of function so that it can work while host compress the rest of the block data
#ifndef __APPLE__
	if ((wb->header.flags_7_0 & kTLBHflagFixedLength) != 0
			&& (gCompress.scoredef == CMP_REDT6 || gCompress.scoredef == CMP_REDT))
	{
		if (gCompress.scoredef == CMP_REDT6)
			fprintf(stderr,"forcing use of new RE instead of selected legacy RE\n");
		const unsigned int tcnt = (wb->header.flags_7_0 & kTLBHflagPairedReads) ? 2*wb->cnt : wb->cnt;
		simpleCompressQS(&wb->readQual,tcnt,wb->lengthInfo,&qs_stat_c,wpos,wmax - wpos);
	}
	else
#endif
		simpleCompress(&wb->readQual,gCompress.scoredef,&qs_stat_c,wpos,wmax - wpos);
	wpos += wb->readQual.packedLen;
	wb->header.blockLength =
		headerSize // header + data block header
		+ (wb->readLen.packedLen != 0 ? 3 * sizeof(unsigned int) + wb->readLen.packedLen : 0) // the read length block if it exists
		+ (wb->unmapped.packedLen != 0 ? 3 * sizeof(unsigned int) + wb->unmapped.packedLen : 0) // the unmapped block if it exists
		+ 3 * sizeof(unsigned int) + wb->readHdr.packedLen // the read headers block
		+ 3 * sizeof(unsigned int) + wb->readQual.packedLen // the read quality scores block
		;
	for (unsigned int i = 0; i < wb->nbMatches; i++)
	{
		wb->header.blockLength +=
			(wb->pairPos[i].packedLen != 0 ? 3 * sizeof(unsigned int) + wb->pairPos[i].packedLen : 0) // the pair positions block if it exists
			+ (wb->mismatch[i].packedLen != 0 ? 3 * sizeof(unsigned int) + wb->mismatch[i].packedLen : 0) // the mismatches if there are some
			+ (wb->cigarBlock[i].packedLen != 0 ? 3 * sizeof(unsigned int) + wb->cigarBlock[i].packedLen : 0) // the cigars if there are some
			;
	}
	assert(wb->header.blockLength == wpos - dataBlock);
	// we love round numbers...  - go for a disk block of 512 bytes...
	unsigned int padding = ~((wb->header.blockLength - 1) & 0x1ff) & 0x1ff;
	memset(wpos,0,padding); // set the end of the block to 0s...
	wb->header.blockLength += padding;
	wb->header.headerCS = 0;
	memcpy(dataBlock,&wb->header,sizeof(TLBHEADER));
	TLBDATA *td = (TLBDATA *) dataBlock;
	// must not alter header after this point
	td->header.headerCS = simple8bitCS(dataBlock,sizeof(TLBHEADER));
	td->cnt = wb->cnt;
	td->lengthInfo = wb->lengthInfo;
	if (multi)
		td->nbMatches = wb->nbMatches;
	// now write block triplets
	wpos = dataBlock + headerSize;
	if (wb->readLen.packedLen != 0)
	{
		memcpy(wpos,&wb->readLen,3 * sizeof(unsigned int));
		wpos += 3 * sizeof(unsigned int);
	}
	if (wb->pairPos[0].packedLen != 0)
		for (unsigned int i = 0; i < wb->nbMatches; i++)
		{
			memcpy(wpos,&wb->pairPos[i],3 * sizeof(unsigned int));
			wpos += 3 * sizeof(unsigned int);
		}
	if (wb->mismatch[0].packedLen != 0)
		for (unsigned int i = 0; i < wb->nbMatches; i++)
		{
			memcpy(wpos,&wb->mismatch[i],3 * sizeof(unsigned int));
			wpos += 3 * sizeof(unsigned int);
		}
	if (wb->cigarBlock[0].packedLen != 0)
		for (unsigned int i = 0; i < wb->nbMatches; i++)
		{
			memcpy(wpos,&wb->cigarBlock[i],3 * sizeof(unsigned int));
			wpos += 3 * sizeof(unsigned int);
		}
	if (wb->unmapped.packedLen != 0)
	{
		memcpy(wpos,&wb->unmapped,3 * sizeof(unsigned int));
		wpos += 3 * sizeof(unsigned int);
	}
	memcpy(wpos,&wb->readHdr,3 * sizeof(unsigned int));
	wpos += 3 * sizeof(unsigned int);
	memcpy(wpos,&wb->readQual,3 * sizeof(unsigned int));
	td->fletcher32CS = fletcher32((uint16_t *)(dataBlock + sizeof(TLBHEADER)), (wb->header.blockLength - sizeof(TLBHEADER)) / 2);
	return dataBlock;
}
//---------------------------------------------------------------

unsigned char *
packBlock(TLBDATA *block, TLBDATA *wb)
{
	unsigned int i;
	if (block != wb)
	{
		zeroBlock(wb);
		memcpy(&wb->header,&block->header,sizeof(TLBHEADER));
		wb->cnt = block->cnt;
		wb->nbMatches = block->nbMatches;
	}
	// If we do not have constant length, we need to keep a block with the length of each read
	if (block->lengthInfo == 0)
	{
		// we do not have constant length
		wb->lengthInfo = 4; // all int for now
		if ((wb->header.flags_7_0 & kTLBHflagPairedReads) == 0)
			wb->readLen.origLen = block->cnt * sizeof(unsigned int);
		else
			wb->readLen.origLen = block->cnt * 2 * sizeof(unsigned int);
		assert(wb->readLen.origLen <= kBlockMaxCount * 8);
		unsigned int *rlp = (unsigned int *) wb->readLen.data;
		// need to set all the length info
		if ((wb->header.flags_7_0 & kTLBHflagPairedReads) == 0)
		{
			for (i = 0; i < block->cnt; i++)
			{
				rlp[i] = block->pdp[i]->len1;
			}
		}
		else
		{
			for (i = 0; i < block->cnt; i++)
			{
				rlp[2 * i] = block->pdp[i]->len1;
				rlp[2 * i + 1] = block->pdp[i]->len2;
			}
		}
	}
	else
	{
		// we do have constant length
		wb->lengthInfo = block->lengthInfo; // the constant length
		wb->header.flags_7_0 |= kTLBHflagFixedLength;
	}
	// for the time being, we have paired_positions or paired_chr_positions with or without unmapped, mismatches, cigars, headers, and qualities
	// Let's prepare a block for each and compress
	wb->header.flags_15_8 |= kTLBHflagHeaderBlock | kTLBHflagQualityBlock;
	PAIRPOSDATA *ppd[kMAXnbMatches] = { NULL };
	PAIRCHRPOSDATA *pcpd[kMAXnbMatches] = { NULL };
	if ((wb->header.flags_7_0 & kTLBHflagPairPosBlock) != 0)
		for (unsigned int j = 0; j < wb->nbMatches; j++)
		{
			wb->pairPos[j].origLen = block->cnt * sizeof(PAIRPOSDATA); // choose encoding with just  pos+offset
			assert(wb->pairPos[j].origLen <= kBlockMaxCount * 12);
			ppd[j] = (PAIRPOSDATA *) wb->pairPos[j].data;
		}
	if ((wb->header.flags_7_0 & kTLBHflagPairChrPosBlock) != 0)
		for (unsigned int j = 0; j < wb->nbMatches; j++)
		{
			wb->pairPos[j].origLen = block->cnt * sizeof(PAIRCHRPOSDATA); // choose encoding with    chr1 chr2 pos1 pos2
			assert(wb->pairPos[j].origLen <= kBlockMaxCount * 12);
			pcpd[j] = (PAIRCHRPOSDATA *) wb->pairPos[j].data;
		}
	unsigned char *wposu = NULL;
	if ((wb->header.flags_15_8 & kTLBHflagUnmappedBlock) == kTLBHflagUnmappedBlock
			|| (wb->header.flags_15_8 & kTLBHflagHalfmapBlock) == kTLBHflagHalfmapBlock)
	{
		wb->unmapped.origLen = block->p4Size;
		assert(wb->unmapped.origLen <= kBlockMaxCount * 96);
		wposu = wb->unmapped.data;
	}

	unsigned char *wposm[kMAXnbMatches] = { NULL };
	unsigned char *wposc[kMAXnbMatches] = { NULL };
	for (unsigned int j = 0; j < wb->nbMatches; j++)
	{
		wb->mismatch[j].origLen = block->mmSize[j];
		assert(wb->mismatch[j].origLen <= kBlockMaxCount * 770);
		wposm[j] = wb->mismatch[j].data;
		wb->cigarBlock[j].origLen = block->cigarSize[j];
		assert(wb->cigarBlock[j].origLen <= kBlockMaxCount * 129);
		wposc[j] = wb->cigarBlock[j].data;
	}
	if ((wb->header.flags_7_0 & kTLBHflagPairedReads) == 0)
		wb->readHdr.origLen = block->hdrSize + block->cnt + 8 * block->cnt; // the headers + a newline at the end + 8 bytes for ordinal
	else
		wb->readHdr.origLen = block->hdrSize + 2 * block->cnt + 8 * block->cnt; // the headers + a tab between the 2 headers and a newline at the end + 8 bytes for ordinal
	assert(wb->readHdr.origLen <= kBlockMaxCount * 266);
	unsigned char *wposh = wb->readHdr.data;
	wb->readQual.origLen = block->qualSize;
	assert(wb->readQual.origLen <= kBlockMaxCount * 512);
	unsigned char *wposq = wb->readQual.data;
	for (i = 0; i < block->cnt; i++)
	{
		PAIRDATA *pi = block->pdp[i];
		for (unsigned int j = 0; j < wb->nbMatches; j++)
		{
			if ((wb->header.flags_7_0 & kTLBHflagPairPosBlock) != 0)
			{
				if ((wb->header.flags_7_0 & kTLBHflagPairedReads) == 0)
				{
					assert(((pi->pcpd[0].chr >> 16) & 0x7FFF) == wb->header.chr);
					assert(pi->pcpd[j].tag2pos == 0);
					ppd[j][i].tag1pos = pi->pcpd[j].tag1pos;
					ppd[j][i].tag2offset = ((pi->pcpd[j].chr & 0x80000000) != 0) ? k32bREVCOMP_TAG1_BIT : 0;
					if (((pi->pcpd[j].chr >> 16) & 0x7FFF) != wb->header.chr)
					{
						// need to encode chromosome
						unsigned int chr2 = (pi->pcpd[j].chr >> 16) & 0x7FFF;
						if (chr2 > (k32b2nd_MATCH_CHR_MASK >> 19))
						{
							fprintf(stderr,"HWHAP: chr of match %u is too large: %u\n",j,chr2);
							exit(1);
						}
						ppd[j][i].tag2offset |= chr2 << 19;
					}
				}
				else
				{
					assert((pi->pcpd[j].chr & 0x7FFF) == ((pi->pcpd[j].chr >> 16) & 0x7FFF));
					assert((pi->pcpd[0].chr & 0x7FFF) == wb->header.chr);
					// pack format for disk block
					int delta = (int)(pi->pcpd[j].tag2pos - pi->pcpd[j].tag1pos);
					if (delta < 0)
					{
						delta = (0-delta);
						if (delta > k32bOFFSET_VALUE_MASK)
						{
							fprintf(stderr,"HWHAP: delta is too big (-) %u: %u %u %u %u %u\n",j,delta,(pi->pcpd[j].chr >> 16) & 0x7FFF,pi->pcpd[j].tag1pos,pi->pcpd[j].chr & 0x7FFF,pi->pcpd[j].tag2pos);
							exit(1);
						}
						delta |= k32bNEGATIVE_OFFSET_BIT;
					}
					else
					{
						if (delta > k32bOFFSET_VALUE_MASK)
						{
							fprintf(stderr,"HWHAP: delta is too big (+) %u: %u %u %u %u %u\n",j,delta,(pi->pcpd[j].chr >> 16) & 0x7FFF,pi->pcpd[j].tag1pos,pi->pcpd[j].chr & 0x7FFF,pi->pcpd[j].tag2pos);
							exit(1);
						}
					}
					if ((pi->pcpd[j].chr & 0x80000000) != 0)
						delta |= k32bREVCOMP_TAG1_BIT;
					if ((pi->pcpd[j].chr & 0x00008000) != 0)
						delta |= k32bREVCOMP_TAG2_BIT;
					if ((pi->pcpd[j].chr & 0x7FFF) != wb->header.chr)
					{
						// need to encode chromosome
						unsigned int chr2 = pi->pcpd[j].chr & 0x7FFF;
						if (chr2 > (k32b2nd_MATCH_CHR_MASK >> 19))
						{
							fprintf(stderr,"HWHAP: chr of match %u is too large: %u\n",j,chr2);
							exit(1);
						}
						delta |= chr2 << 19;
					}
					ppd[j][i].tag1pos = pi->pcpd[j].tag1pos;
					ppd[j][i].tag2offset = delta;
				}
			}
			if ((wb->header.flags_7_0 & kTLBHflagPairChrPosBlock) != 0)
				pcpd[j][i] = pi->pcpd[j];
			memcpy(wposm[j],pi->MMstring[j],pi->ml[j]);
			wposm[j] += pi->ml[j];
			memcpy(wposc[j],pi->cigar[j],pi->cl[j]);
			wposc[j] += pi->cl[j];
		}
		if ((wb->header.flags_15_8 & kTLBHflagUnmappedBlock) == kTLBHflagUnmappedBlock
				|| (wb->header.flags_15_8 & kTLBHflagHalfmapBlock) == kTLBHflagHalfmapBlock)
		{
			memcpy(wposu,pi->packed4nt,pi->pl);
			wposu += pi->pl;
		}
		memcpy(wposh,pi->hdr1,pi->hlen1);
		wposh += pi->hlen1;
		if ((wb->header.flags_7_0 & kTLBHflagPairedReads) != 0)
		{
			*wposh++ = '\t';
			memcpy(wposh,pi->hdr2,pi->hlen2);
			wposh += pi->hlen2;
		}
		*wposh++ = '\n';
		{
			unsigned int j;
			for (j = 0; j < 8; j++)
				*wposh++ = (pi->ordinal >> (56 - j * 8)) & 0xff;
		}
		memcpy(wposq,pi->qs1,pi->len1);
		wposq += pi->len1;
		if ((wb->header.flags_7_0 & kTLBHflagPairedReads) != 0)
		{
			memcpy(wposq,pi->qs2,pi->len2);
			wposq += pi->len2;
		}
	}
	return compressBlock(wb);
}
//---------------------------------------------------------------

static int
DecodeQS6m(char * streamout, const unsigned char * restrict const streamin)
{
	rangecoder rc;
	int res = -1;

	if (start_decoding(&rc, streamin) != 0)
	{
		fprintf(stderr,"could not suceessfully open input data\n");
		return -1;
	}

	if (decode_culfreq(&rc,2))
	{
		freq j, i, blocksize;

		decode_update(&rc,1,1,2);

		freq LineLength;
		short dummy = decode_short(&rc);
		LineLength = (freq) dummy;
		blocksize =  decode_culshift(&rc, 24);
		decode_update(&rc, (freq) 1, blocksize, (freq) 1 << 24);

		//fprintf(stderr, "%lu reads of length %lu\n", blocksize, LineLength);
		res = blocksize * LineLength;
		char * restrict out = streamout;

		freq* counters = (freq*) malloc(FROMTO_SIZE*FROMTO_SIZE*LineLength*sizeof(freq));
		if (counters == NULL) {
			fprintf(stderr,"Unable to allocate memory\n");
			return -2;
		}
		freq avant = 0;
		for (i=0; i<FROMTO_SIZE*FROMTO_SIZE*LineLength; i++) {
			freq tmp = decode_culshift(&rc, 18);
			counters[i] = tmp;
			if (tmp > ((freq) (1 << 17))) {
				const freq length = tmp - 1 - ((freq) (1 << 17));
				for (freq k=0; k<length; k++) counters[k+i] = avant;
				i += length-1;
			}
			else {
				counters[i] = tmp;
				avant = tmp;
			}
			decode_update(&rc, (freq) 1, tmp, (freq) 1 << 18);
		}

		{
			freq * restrict colcounters = counters;
			for (freq k=0; k<FROMTO_SIZE*LineLength; k++) {
				freq sum = 0;
				for (freq j=0; j<FROMTO_SIZE; j++) {
					const freq tmp = colcounters[j];
					colcounters[j] = sum;
					sum += tmp;
				}
				colcounters += FROMTO_SIZE;
			}
		}

		int previous = 0;
		freq cf, symbol;
		for (i=0; i<blocksize; i++) {
			freq * colcounters = counters;
			for (j=0; j<LineLength; j++) {
				cf = decode_culfreq(&rc,colcounters[previous*FROMTO_SIZE+FROMTO_SIZE-1]);
				for (symbol=0; colcounters[previous*FROMTO_SIZE+symbol+1] <= cf; symbol++)
					;
				decode_update(&rc,
					colcounters[previous*FROMTO_SIZE+symbol+1]-colcounters[previous*FROMTO_SIZE+symbol],
					colcounters[previous*FROMTO_SIZE+symbol],
					colcounters[previous*FROMTO_SIZE+FROMTO_SIZE-1]);
				for(freq k=symbol+1; k<FROMTO_SIZE; k++)
					colcounters[previous*FROMTO_SIZE+k]--;
				previous = (int) symbol;

				if (symbol == 0 && j<(LineLength-1)) {
					for (;j<LineLength;j++) { *out++ = '#'; }
					break;
				}
				else {
					*out++ = 33-1+symbol;
				}
				colcounters += FROMTO_SIZE*FROMTO_SIZE;
			}
		}
		free(counters);
	}

	done_decoding(&rc);

	return res;
}
//---------------------------------------------------------------

static int
DecodeQS7m(char * const streamout, const unsigned char * restrict const streamin)
{
	rangecoder rc;
	int res = -1;

	/* Open a stream for writing */
	if (start_decoding(&rc, streamin) != 0)
	{   fprintf(stderr,"could not suceessfully open input data\n");
		return -1;
	}

	if (decode_culfreq(&rc,2))
	{   freq j, i, blocksize, maxRunLength, maxValue;

		decode_update(&rc,1,1,2);

		const short dummy = decode_short(&rc);
		const freq LineLength = (freq) dummy;
		blocksize =  decode_culshift(&rc, 24);
		decode_update(&rc, (freq) 1, blocksize, (freq) 1 << 24);
		maxRunLength =  decode_culshift(&rc, 24);
		decode_update(&rc, (freq) 1, maxRunLength, (freq) 1 << 24);
		maxValue = decode_culshift(&rc, 24);
		decode_update(&rc, (freq) 1, maxValue, (freq) 1 << 24);

		res = blocksize * LineLength;
		char * restrict out = streamout;

		freq* counters = (freq*) malloc(FROMTO_SIZE*FROMTO_SIZE*LineLength*sizeof(freq));
		if (counters == NULL) {
			fprintf(stderr,"Unable to allocate memory\n");
			return -2;
		}

		const freq threshold = maxValue - maxRunLength;
		maxValue++;
		{
			freq i=0;
			while (i<FROMTO_SIZE*FROMTO_SIZE*LineLength) {
				freq tmp = decode_culfreq(&rc, maxValue);

				if (tmp > threshold) {
					const freq length = tmp - threshold;
					for (freq k=0; k<length; k++) counters[k+i] = 0;
					i += length;
				}
				else {
					counters[i] = tmp;
					i++;
				}
				decode_update(&rc, (freq) 1, tmp, maxValue);
			}
		}

		/* compute CDF */
		{
			freq * restrict colcounters = counters;
			for (freq k=0; k<FROMTO_SIZE*LineLength; k++) {
				freq sum = 0;
				for (freq j=0; j<FROMTO_SIZE; j++) {
					const freq tmp = colcounters[j];
					colcounters[j] = sum;
					sum += tmp;
				}
				colcounters += FROMTO_SIZE;
			}
		}

		int previous = 0;
		freq cf, symbol;
		for (i=0; i<blocksize; i++) {
			freq * colcounters = counters;
			for (j=0; j<LineLength; j++) {
				cf = decode_culfreq(&rc,colcounters[previous*FROMTO_SIZE+FROMTO_SIZE-1]);

				/* Faster binary search */
				if (colcounters[previous*FROMTO_SIZE+FROMTO_SIZE/2+1] <= cf) {
					if (colcounters[previous*FROMTO_SIZE+3*FROMTO_SIZE/4+1] <= cf) {
						symbol = 3*FROMTO_SIZE/4+1;
					}
					else {
						symbol = FROMTO_SIZE/2+1;
					}
				}
				else {
					if (colcounters[previous*FROMTO_SIZE+FROMTO_SIZE/4+1] <= cf) {
						symbol = FROMTO_SIZE/4+1;
					}
					else {
						symbol = 0;
					}
				}
				// 	    if (colcounters[previous*FROMTO_SIZE+symbol+1+5] <= cf) symbol += 5;
				for (; colcounters[previous*FROMTO_SIZE+symbol+1] <= cf; symbol++);

				decode_update(&rc,
					colcounters[previous*FROMTO_SIZE+symbol+1]-colcounters[previous*FROMTO_SIZE+symbol],
					colcounters[previous*FROMTO_SIZE+symbol],
					colcounters[previous*FROMTO_SIZE+FROMTO_SIZE-1]);
				for(freq k=symbol+1; k<FROMTO_SIZE; k++) colcounters[previous*FROMTO_SIZE+k]--;
				previous = (int) symbol;

				if (symbol == 0 && j<(LineLength-1)) {
					for (;j<LineLength;j++) { *out++ = '#'; }
					break;
				}
				else {
					*out++ = 33-1+symbol;
				}
				colcounters += FROMTO_SIZE*FROMTO_SIZE;
			}
		}
		free(counters);
	}

	done_decoding(&rc);

	return res;
}
//---------------------------------------------------------------

static void
simpleUncompress(BASICBLOCK *bb, const unsigned char *src, compress_stat_t *stat)
{
	struct timeval ts, te, td;

	gettimeofday(&ts,NULL);
	stat->compressor = bb->compressor;
	switch(bb->compressor)
	{
		case CMP_RAW:
			memcpy(bb->data,src,bb->origLen);
			break;
		case CMP_ZIP:
		{
			unsigned long size = bb->origLen;
			int res = uncompress(bb->data,&size,src,bb->packedLen);
			if (res != Z_OK || size != bb->origLen)
			{
				fprintf(stderr, "Zlib's unpack failed: %d or %ld != %d\n", res, size, bb->origLen);
				exit(1);
			}
			break;
		}
		case CMP_BZ2:
		{
			unsigned size = bb->origLen;
			int res = BZ2_bzBuffToBuffDecompress((char*)bb->data,&size,(char*)src,bb->packedLen,0,0);
			if (res != BZ_OK || size != bb->origLen)
			{
				fprintf(stderr, "libBZ2's unpack failed: %d or %u != %u\n", res, size, bb->origLen);
				exit(1);
			}
			break;
		}
		case CMP_REDT6:
		{
			unsigned long size = bb->origLen;
			int res = DecodeQS6m((char *) bb->data,src);
			if (res != bb->origLen)
			{
				fprintf(stderr, "Range Encoder6's unpack failed: %d or %ld != %d\n", res, size, bb->origLen);
				exit(1);
			}
			break;
		}
		case CMP_REDT:
		{
			unsigned long size = bb->origLen;
			int res = DecodeQS7m((char *) bb->data,src);
			if (res != bb->origLen)
			{
				fprintf(stderr, "Range Encoder's unpack failed: %d or %ld != %d\n", res, size, bb->origLen);
				exit(1);
			}
			break;
		}
		case CMP_LZ4:
		{
			LZ4F_decompressOptions_t opt =
			{
				.stableDst = 1,
			};
			size_t dSize = bb->origLen, pSize = bb->packedLen;
			LZ4F_decompressionContext_t dctx;
			LZ4F_createDecompressionContext(&dctx, LZ4F_VERSION);
			size_t res = LZ4F_decompress(dctx,bb->data,&dSize,src,&pSize,&opt);
			// HACK in theory decompress expects a loop, but we try all in one go
			if (LZ4F_isError(res) || dSize < bb->origLen) {
				fprintf(stderr, "LZ4 decompress failed (%ld) : %s\n", res, LZ4F_getErrorName(res));
				exit(1);
			}
			if (dSize < bb->origLen) {
				fprintf(stderr, "LZ4 wrong size %lu != %u\nmaybe we need a loop, after all...", dSize, bb->origLen);
				exit(1);
			}
			LZ4F_freeDecompressionContext(dctx);
			break;
		}
		default:
			fprintf(stderr, "Unknown compressor %d\n",bb->compressor);
			exit(1);
	}
	gettimeofday(&te,NULL);
	timersub(&te, &ts, &td);
	timeradd(&td, &(stat->tv), &(stat->tv));
	stat->raw_bytes += bb->origLen;
	stat->packed_bytes += bb->packedLen;
}
//---------------------------------------------------------------

static void
decodeHeaderSingle(TLBDATA *tpp, PAIRDATA *pdp)
{
	unsigned int i;
	char *next = strchr(tpp->hdr, '\n');

	pdp->hlen1 = next - tpp->hdr;
	memcpy(pdp->hdr1,tpp->hdr,pdp->hlen1);
	//pdp->hdr1[pdp->hlen1] = 0;
	tpp->hdr = next + 1;
	pdp->hlen2 = 0;
	pdp->ordinal = 0;
	for (i = 0; i < 8; i++)
		pdp->ordinal = (pdp->ordinal << 8) | (unsigned long) *((unsigned char *) tpp->hdr++);
}
//---------------------------------------------------------------

static void
decodeHeaderPair(TLBDATA *tpp, PAIRDATA *pdp)
{
	unsigned int i;
	char *next = strchr(tpp->hdr, '\t');

	pdp->hlen1 = next - tpp->hdr;
	memcpy(pdp->hdr1,tpp->hdr,pdp->hlen1);
	//pdp->hdr1[pdp->hlen1] = 0;
	tpp->hdr = next + 1;
	next = strchr(tpp->hdr, '\n');
	pdp->hlen2 = next - tpp->hdr;
	memcpy(pdp->hdr2,tpp->hdr,pdp->hlen2);
	//pdp->hdr2[pdp->hlen2] = 0;
	tpp->hdr = next + 1;
	pdp->ordinal = 0;
	for (i = 0; i < 8; i++)
		pdp->ordinal = (pdp->ordinal << 8) | (unsigned long) *((unsigned char *) tpp->hdr++);
}
//---------------------------------------------------------------

static void
decodeSingle(TLBDATA *tpp, unsigned int i, PAIRDATA *pdp)
{
	for (unsigned int i = 0; i < kMAXnbMatches; i++)
	{
		pdp->ml[i] = 0;
		pdp->cl[i] = 0;
	}
	pdp->pl = 0;
	if (tpp->header.flags_7_0 & kTLBHflagFixedLength)
	{
		pdp->len1 = tpp->lengthInfo;
	}
	else
	{
		pdp->len1 = tpp->len[i];
	}
	pdp->nbMatches = tpp->nbMatches;
	for (unsigned int j = 0; j < tpp->nbMatches; j++)
	{
		pdp->pcpd[j].chr = 0;
		pdp->pcpd[j].tag2pos = 0;
		if (tpp->ppd[j] != NULL)
		{
			unsigned int chr = (tpp->ppd[j][i].tag2offset & k32b2nd_MATCH_CHR_MASK) >> 19;
			if (chr == 0)
				chr = tpp->header.chr;
			pdp->pcpd[j].tag1pos = tpp->ppd[j][i].tag1pos;
			pdp->pcpd[j].chr = chr << 16;
			if (tpp->ppd[j][i].tag2offset & k32bREVCOMP_TAG1_BIT)
				pdp->pcpd[j].chr |= 0x80000000;
		}
	}
	memcpy(pdp->qs1,tpp->qs,pdp->len1);
	pdp->qs1[pdp->len1] = 0;
	tpp->qs += pdp->len1;
	decodeHeaderSingle(tpp,pdp);
	if (tpp->p4nt != NULL)
	{
		unsigned int pl = pdp->len1 >> 2;
		if ((pdp->len1 & 3) != 0)
			pl += 1;
		memcpy(pdp->packed4nt,tpp->packed4nt,pl);
		tpp->packed4nt += pl;
		pdp->pl = pl;
	}
}
//---------------------------------------------------------------

static void
decodePair(TLBDATA *tpp, unsigned int i, PAIRDATA *pdp)
{
	for (unsigned int i = 0; i < kMAXnbMatches; i++)
	{
		pdp->ml[i] = 0;
		pdp->cl[i] = 0;
	}
	pdp->pl = 0;
	if (tpp->header.flags_7_0 & kTLBHflagFixedLength)
	{
		pdp->len1 = tpp->lengthInfo;
		pdp->len2 = tpp->lengthInfo;
	}
	else
	{
		pdp->len1 = tpp->len[2*i];
		pdp->len2 = tpp->len[2*i+1];
	}
	pdp->nbMatches = tpp->nbMatches;
	for (unsigned int j = 0; j < tpp->nbMatches; j++)
	{
		if (tpp->ppd[j] != NULL)
		{
			unsigned int chr = (tpp->ppd[j][i].tag2offset & k32b2nd_MATCH_CHR_MASK) >> 19;
			if (chr == 0)
				chr = tpp->header.chr;
			pdp->pcpd[j].tag1pos = tpp->ppd[j][i].tag1pos;
			pdp->pcpd[j].chr = (chr << 16) | chr;
			if (tpp->ppd[j][i].tag2offset & k32bREVCOMP_TAG1_BIT)
				pdp->pcpd[j].chr |= 0x80000000;
			if (tpp->ppd[j][i].tag2offset & k32bREVCOMP_TAG2_BIT)
				pdp->pcpd[j].chr |= 0x00008000;
			int delta = (tpp->ppd[j][i].tag2offset & k32bOFFSET_VALUE_MASK);
			if (tpp->ppd[j][i].tag2offset & k32bNEGATIVE_OFFSET_BIT)
				delta = -delta;
			pdp->pcpd[j].tag2pos = pdp->pcpd[j].tag1pos + delta;
		}
		else if (tpp->pcpd[j] != NULL)
			pdp->pcpd[j] = tpp->pcpd[j][i];
		else
		{
			pdp->pcpd[j].chr = 0;
			pdp->pcpd[j].tag1pos = 0;
			pdp->pcpd[j].tag2pos = 0;
		}
	}
	memcpy(pdp->qs1,tpp->qs,pdp->len1);
	pdp->qs1[pdp->len1] = 0;
	tpp->qs += pdp->len1;
	memcpy(pdp->qs2,tpp->qs,pdp->len2);
	pdp->qs2[pdp->len2] = 0;
	tpp->qs += pdp->len2;
	decodeHeaderPair(tpp,pdp);
	if (tpp->p4nt != NULL)
	{
		unsigned int chr1 = (pdp->pcpd[0].chr >> 16) & 0x7fff;
		unsigned int chr2 = pdp->pcpd[0].chr & 0x7fff;
		if (chr1 == 0)
		{
			unsigned int pl = pdp->len1 >> 2;
			if ((pdp->len1 & 3) != 0)
				pl += 1;
			memcpy(pdp->packed4nt,tpp->packed4nt,pl);
			tpp->packed4nt += pl;
			pdp->pl = pl;
		}
		if (chr2 == 0)
		{
			unsigned int pl = pdp->len2 >> 2;
			if ((pdp->len2 & 3) != 0)
				pl += 1;
			memcpy(pdp->packed4nt + pdp->pl,tpp->packed4nt,pl);
			tpp->packed4nt += pl;
			pdp->pl += pl;
		}
	}
}
//---------------------------------------------------------------

static int
addDiffSingle(PAIRDATA *pdp, unsigned char *diff, int mmpos, unsigned int i)
{
	pdp->ml[i] = 0;
	while (diff[mmpos] != kMMTERMINATOR)
	{
		pdp->MMstring[i][pdp->ml[i]++] = diff[mmpos++];
	}
	pdp->MMstring[i][pdp->ml[i]++] = diff[mmpos++]; // copy the kMMTERMINATOR
	return mmpos;
}
//---------------------------------------------------------------

static int
addDiffPair(PAIRDATA *pdp, unsigned char *diff, int mmpos, unsigned int i)
{
	pdp->ml[i] = 0;
	while (diff[mmpos] != kMMTERMINATOR)
	{
		pdp->MMstring[i][pdp->ml[i]++] = diff[mmpos++];
	}
	pdp->MMstring[i][pdp->ml[i]++] = diff[mmpos++]; // copy the kMMTERMINATOR
	while (diff[mmpos] != kMMTERMINATOR)
	{
		pdp->MMstring[i][pdp->ml[i]++] = diff[mmpos++];
	}
	pdp->MMstring[i][pdp->ml[i]++] = diff[mmpos++]; // copy the kMMTERMINATOR
	return mmpos;
}
//---------------------------------------------------------------

static int
addCigarSingle(PAIRDATA *pdp, unsigned char *cigar, int cigarpos, unsigned int i)
{
	pdp->cl[i] = 0;
	while (cigar[cigarpos] != kCIGARTERMINATOR)
	{
		pdp->cigar[i][pdp->cl[i]++] = cigar[cigarpos++];
	}
	pdp->cigar[i][pdp->cl[i]++] = cigar[cigarpos++]; // copy the kCIGARTERMINATOR
	return cigarpos; // skip the kCIGARTERMINATOR
}
//---------------------------------------------------------------

static int
addCigarPair(PAIRDATA *pdp, unsigned char *cigar, int cigarpos, unsigned int i)
{
	pdp->cl[i] = 0;
	while (cigar[cigarpos] != kCIGARTERMINATOR)
	{
		pdp->cigar[i][pdp->cl[i]++] = cigar[cigarpos++];
	}
	pdp->cigar[i][pdp->cl[i]++] = cigar[cigarpos++]; // copy the kCIGARTERMINATOR
	while (cigar[cigarpos] != kCIGARTERMINATOR)
	{
		pdp->cigar[i][pdp->cl[i]++] = cigar[cigarpos++];
	}
	pdp->cigar[i][pdp->cl[i]++] = cigar[cigarpos++]; // copy the kCIGARTERMINATOR
	return cigarpos; // skip the kCIGARTERMINATOR
}
//---------------------------------------------------------------

int
decompressBlock(void *dataBlock, const unsigned int size, TLBHEADER *thp, TLBDATA *tdp, int seqOnly)
{
	zeroBlock(tdp);
	memcpy(&(tdp->header),thp,sizeof(TLBHEADER));
	unsigned char *rp;
	const int multi = (thp->flags_15_8 & kTLBHflagMultipleMatches) != 0;
	if (multi)
	{
		memcpy(&tdp->fletcher32CS,dataBlock,4 * sizeof(unsigned int)); // get CS, cnt, lengthInfo and nbMatches
		rp = dataBlock + 4 * sizeof(unsigned int);
	}
	else
	{
		memcpy(&tdp->fletcher32CS,dataBlock,3 * sizeof(unsigned int)); // get CS, cnt and lengthInfo
		rp = dataBlock + 3 * sizeof(unsigned int);
		tdp->nbMatches = 1;
	}
	
	if (tdp->cnt > kBlockMaxCount) {
		fprintf(stderr, "Maxiumum block size exceeded, please recompile with kBlockMaxCount > %u\n", tdp->cnt);
		return 1;
	}
	
	memset(dataBlock,0,sizeof(unsigned int));
	unsigned int cs = fletcher32((uint16_t *) dataBlock, size / 2);
	if (tdp->fletcher32CS != cs)
	{
		fprintf(stderr,"Bad checksum : %08x != %08x\n", tdp->fletcher32CS, cs);
		return 1;
	}
	// need to unpack the various blocks
	if ((thp->flags_7_0 & kTLBHflagFixedLength) == 0)
	{
		memcpy(&tdp->readLen,rp,3 * sizeof(unsigned int));
		rp += 3 * sizeof(unsigned int);
		//fprintf(stderr, "readLen %d %d %d\n", tdp->readLen.compressor, tdp->readLen.origLen, tdp->readLen.packedLen);
	}
	if ((thp->flags_7_0 & kTLBHflagPairPosBlock) != 0)
		for (unsigned int i = 0; i < tdp->nbMatches; i++)
		{
			memcpy(&tdp->pairPos[i],rp,3 * sizeof(unsigned int));
			rp += 3 * sizeof(unsigned int);
			//fprintf(stderr, "pairPos[i] %d %d %d\n", tdp->pairPos[i].compressor, tdp->pairPos[i].origLen, tdp->pairPos[i].packedLen);
		}
	if ((thp->flags_7_0 & kTLBHflagPairChrPosBlock) != 0)
		for (unsigned int i = 0; i < tdp->nbMatches; i++)
		{
			memcpy(&tdp->pairPos[i],rp,3 * sizeof(unsigned int));
			rp += 3 * sizeof(unsigned int);
			//fprintf(stderr, "pairPos[i] %d %d %d\n", tdp->pairPos[i].compressor, tdp->pairPos[i].origLen, tdp->pairPos[i].packedLen);
		}
	if ((thp->flags_7_0 & kTLBHflagNMismatchBlock) != 0)
		for (unsigned int i = 0; i < tdp->nbMatches; i++)
		{
			memcpy(&tdp->mismatch[i],rp,3 * sizeof(unsigned int));
			rp += 3 * sizeof(unsigned int);
			//fprintf(stderr, "N mismatches %d %d %d\n", tdp->mismatch[i].compressor, tdp->mimatch[i].origLen, tdp->mismatch[i].packedLen);
		}
	if ((thp->flags_7_0 & kTLBHflagMismatchBlock) != 0)
	{
		if ((thp->flags_7_0 & kTLBHflagNMismatchBlock) != 0)
		{
			fprintf(stderr, "HWHAP: can't have both N-only and normal mismatches");
			return 1;
		}
		for (unsigned int i = 0; i < tdp->nbMatches; i++)
		{
			memcpy(&tdp->mismatch[i],rp,3 * sizeof(unsigned int));
			rp += 3 * sizeof(unsigned int);
			//fprintf(stderr, "mismatches %d %d %d\n", tdp->mismatch[i].compressor, tdp->mimatch[i].origLen, tdp->mismatch[i].packedLen);
		}
	}
	if ((thp->flags_7_0 & kTLBHflagCigarBlock) != 0)
		for (unsigned int i = 0; i < tdp->nbMatches; i++)
		{
			memcpy(&tdp->cigarBlock[i],rp,3 * sizeof(unsigned int));
			rp += 3 * sizeof(unsigned int);
			//fprintf(stderr, "Cigar %d %d %d\n", tdp->cigarBlock[i].compressor, tdp->cigarBlock[i].origLen, tdp->cigarBlock[i].packedLen);
		}
	if ((thp->flags_15_8 & kTLBHflagUnmappedBlock) != 0 || (thp->flags_15_8 & kTLBHflagHalfmapBlock) != 0)
	{
		memcpy(&tdp->unmapped,rp,3 * sizeof(unsigned int));
		rp += 3 * sizeof(unsigned int);
		//fprintf(stderr, "Unmapped %d %d %d\n", tdp->unmapped.compressor, tdp->unmapped.origLen, tdp->unmapped.packedLen);
	}
	if ((thp->flags_15_8 & kTLBHflagHeaderBlock) != 0)
	{
		memcpy(&tdp->readHdr,rp,3 * sizeof(unsigned int));
		rp += 3 * sizeof(unsigned int);
		//fprintf(stderr, "readHdr %d %d %d\n", tdp->readHdr.compressor, tdp->readHdr.origLen, tdp->readHdr.packedLen);
	}
	if ((thp->flags_15_8 & kTLBHflagQualityBlock) != 0)
	{
		memcpy(&tdp->readQual,rp,3 * sizeof(unsigned int));
		rp += 3 * sizeof(unsigned int);
		//fprintf(stderr, "readQual %d %d %d\n", tdp->readQual.compressor, tdp->readQual.origLen, tdp->readQual.packedLen);
	}
	// get and uncompress the bits
	if ((thp->flags_7_0 & kTLBHflagFixedLength) == 0)
	{
		assert(tdp->readLen.origLen <= kBlockMaxCount * 8);
		simpleUncompress(&tdp->readLen,rp,&len_stat_d);
		rp += tdp->readLen.packedLen;
		tdp->len = (unsigned int *) tdp->readLen.data;
	}
	if ((thp->flags_7_0 & kTLBHflagPairPosBlock) != 0)
		for (unsigned int i = 0; i < tdp->nbMatches; i++)
		{
			assert(tdp->pairPos[i].origLen <= kBlockMaxCount * 12);
			simpleUncompress(&tdp->pairPos[i],rp,&pairs_stat_d);
			rp += tdp->pairPos[i].packedLen;
			tdp->ppd[i] = (PAIRPOSDATA *) tdp->pairPos[i].data;
		}
	if ((thp->flags_7_0 & kTLBHflagPairChrPosBlock) != 0)
		for (unsigned int i = 0; i < tdp->nbMatches; i++)
		{
			assert(tdp->pairPos[i].origLen <= kBlockMaxCount * 12);
			simpleUncompress(&tdp->pairPos[i],rp,&pairs_stat_d);
			rp += tdp->pairPos[i].packedLen;
			tdp->pcpd[i] = (PAIRCHRPOSDATA *) tdp->pairPos[i].data;
		}
	if ((thp->flags_7_0 & kTLBHflagNMismatchBlock) != 0)
		for (unsigned int i = 0; i < tdp->nbMatches; i++)
		{
			assert(tdp->mismatch[i].origLen <= kBlockMaxCount * 770);
			simpleUncompress(&tdp->mismatch[i],rp,&mismatch_stat_d);
			rp += tdp->mismatch[i].packedLen;
			tdp->diff[i] = tdp->mismatch[i].data;
		}
	if ((thp->flags_7_0 & kTLBHflagMismatchBlock) != 0)
		for (unsigned int i = 0; i < tdp->nbMatches; i++)
		{
			assert(tdp->mismatch[i].origLen <= kBlockMaxCount * 770);
			simpleUncompress(&tdp->mismatch[i],rp,&mismatch_stat_d);
			rp += tdp->mismatch[i].packedLen;
			tdp->diff[i] = tdp->mismatch[i].data;
		}
	if ((thp->flags_7_0 & kTLBHflagCigarBlock) != 0)
		for (unsigned int i = 0; i < tdp->nbMatches; i++)
		{
			assert(tdp->cigarBlock[i].origLen <= kBlockMaxCount * 129);
			simpleUncompress(&tdp->cigarBlock[i],rp,&cigar_stat_d);
			rp += tdp->cigarBlock[i].packedLen;
			tdp->cigar[i] = tdp->cigarBlock[i].data;
		}
	if ((thp->flags_15_8 & kTLBHflagUnmappedBlock) != 0 || (thp->flags_15_8 & kTLBHflagHalfmapBlock) != 0)
	{
		assert(tdp->unmapped.origLen <= kBlockMaxCount * 96);
		simpleUncompress(&tdp->unmapped,rp,&unmapped_stat_d);
		rp += tdp->unmapped.packedLen;
		tdp->p4nt = tdp->unmapped.data;
		tdp->packed4nt = tdp->p4nt;
	}
	if ((thp->flags_15_8 & kTLBHflagHeaderBlock) != 0)
	{
		if (!seqOnly)
		{
			assert(tdp->readHdr.origLen <= kBlockMaxCount * 266);
			simpleUncompress(&tdp->readHdr,rp,&header_stat_d);
			tdp->hdrs = tdp->readHdr.data;
			tdp->hdr = (char *) tdp->hdrs;
		}
		rp += tdp->readHdr.packedLen;
	}
	if ((thp->flags_15_8 & kTLBHflagQualityBlock) != 0)
	{
		if (!seqOnly)
		{
			assert(tdp->readQual.origLen <= kBlockMaxCount * 512);
			simpleUncompress(&tdp->readQual,rp,&qs_stat_d);
			tdp->qual = tdp->readQual.data;
			tdp->qs = tdp->qual;
		}
		rp += tdp->readQual.packedLen;
	}
	return 0;
}
//---------------------------------------------------------------

// Basically same as decompressBlock but only count patterns
unsigned int
countPatInUnmappedBlock(void *dataBlock, unsigned int size, TLBHEADER *thp, unsigned int *patCnt1, unsigned int *patCnt2)
{
	TLBDATA td;
	TLBDATA *tdp = &td;
	memset(tdp,0,sizeof(TLBDATA));
	unsigned char *rp;
	const int multi = (thp->flags_15_8 & kTLBHflagMultipleMatches) != 0;
	if (multi)
	{
		memcpy(&tdp->fletcher32CS,dataBlock,4 * sizeof(unsigned int)); // get CS, cnt, lengthInfo and nbMatches
		rp = dataBlock + 4 * sizeof(unsigned int);
	}
	else
	{
		memcpy(&tdp->fletcher32CS,dataBlock,3 * sizeof(unsigned int)); // get CS, cnt and lengthInfo
		rp = dataBlock + 3 * sizeof(unsigned int);
		tdp->nbMatches = 1;
	}
	// need to unpack the various blocks
	if ((thp->flags_7_0 & kTLBHflagFixedLength) == 0)
	{
		memcpy(&tdp->readLen,rp,3 * sizeof(unsigned int));
		rp += 3 * sizeof(unsigned int);
	}
	assert((thp->flags_7_0 & kTLBHflagPairPosBlock) == 0);
	assert((thp->flags_7_0 & kTLBHflagPairChrPosBlock) == 0);
	assert((thp->flags_7_0 & kTLBHflagNMismatchBlock) == 0);
	if ((thp->flags_7_0 & kTLBHflagMismatchBlock) != 0)
		for (unsigned int i = 0; i < tdp->nbMatches; i++)
		{
			memcpy(&tdp->mismatch[i],rp,3 * sizeof(unsigned int));
			rp += 3 * sizeof(unsigned int);
		}
	assert((thp->flags_7_0 & kTLBHflagCigarBlock) == 0);
	assert((thp->flags_15_8 & kTLBHflagHalfmapBlock) == 0);
	if ((thp->flags_15_8 & kTLBHflagUnmappedBlock) != 0)
	{
		memcpy(&tdp->unmapped,rp,3 * sizeof(unsigned int));
		rp += 3 * sizeof(unsigned int);
	}
	if ((thp->flags_15_8 & kTLBHflagHeaderBlock) != 0)
	{
		memcpy(&tdp->readHdr,rp,3 * sizeof(unsigned int));
		rp += 3 * sizeof(unsigned int);
	}
	if ((thp->flags_15_8 & kTLBHflagQualityBlock) != 0)
	{
		memcpy(&tdp->readQual,rp,3 * sizeof(unsigned int));
		rp += 3 * sizeof(unsigned int);
	}
	// get and uncompress the bits
	if ((thp->flags_7_0 & kTLBHflagFixedLength) == 0)
	{
		tdp->len = malloc(tdp->readLen.origLen);
		tdp->readLen.data = (unsigned char *)tdp->len;
		simpleUncompress(&tdp->readLen,rp,&len_stat_d);
		rp += tdp->readLen.packedLen;
	}
	if ((thp->flags_7_0 & kTLBHflagMismatchBlock) != 0)
		for (unsigned int i = 0; i < tdp->nbMatches; i++)
		{
			tdp->mismatch[i].data = rp;
			rp += tdp->mismatch[i].packedLen;
		}
	if ((thp->flags_15_8 & kTLBHflagUnmappedBlock) != 0)
	{
		tdp->p4nt = malloc(tdp->unmapped.origLen);
		tdp->unmapped.data = tdp->p4nt;
		simpleUncompress(&tdp->unmapped,rp,&unmapped_stat_d);
		rp += tdp->unmapped.packedLen;
		tdp->packed4nt = tdp->p4nt;
	}
	if (tdp->p4nt != NULL)
	{
		unsigned int i;
		if ((thp->flags_7_0 & kTLBHflagPairedReads) == 0)
		{
			for (i = 0; i < tdp->cnt; i++)
			{
				unsigned int len1 = tdp->lengthInfo;
				if ((thp->flags_7_0 & kTLBHflagFixedLength) == 0)
				{
					len1 = tdp->len[i];
				}
				unsigned int pl = len1 >> 2;
				if ((len1 & 3) != 0)
					pl += 1;
				unsigned int pat = ((tdp->packed4nt[0] & 0x3f) << 8) | tdp->packed4nt[1];
				patCnt1[pat] += 1;
				tdp->packed4nt += pl;
			}
		}
		else
		{
			for (i = 0; i < tdp->cnt; i++)
			{
				unsigned int len1 = tdp->lengthInfo;
				unsigned int len2 = tdp->lengthInfo;
				if ((thp->flags_7_0 & kTLBHflagFixedLength) == 0)
				{
					len1 = tdp->len[2*i];
					len2 = tdp->len[2*i+1];
				}
				unsigned int pl = len1 >> 2;
				if ((len1 & 3) != 0)
					pl += 1;
				unsigned int pat = ((tdp->packed4nt[0] & 0x3f) << 8) | tdp->packed4nt[1];
				patCnt1[pat] += 1;
				tdp->packed4nt += pl;
				pl = len2 >> 2;
				if ((len2 & 3) != 0)
					pl += 1;
				pat = ((tdp->packed4nt[0] & 0x3f) << 8) | tdp->packed4nt[1];
				patCnt2[pat] += 1;
				tdp->packed4nt += pl;
			}
		}
	}
	free(tdp->len);
	free(tdp->p4nt);
	return tdp->cnt;
}
//---------------------------------------------------------------

int
readBlock(int fd, TLBHEADER *thp, TLBDATA *tdp)
{
	unsigned int size = thp->blockLength - sizeof(TLBHEADER);
	size_t rs = size;
	assert(size <= kMaxBlockLength);
	void *rp = tdp->diskBuffer; assert(rp);
	while (1)
	{
		ssize_t res = read(fd,rp,rs);
		if (res == rs)
			break;
		if (res > 0)
		{
			rp += res;
			rs -= res;
			continue;
		}
		if (res == 0)
		{
			fprintf(stderr,"unexpected EOF while reading block in readBlock\n");
			return 1;
		}
		if (res == -1 && (errno == EINTR || errno == EAGAIN))
			continue;
		perror("readBlock read:");
		return 1;
	}
	return 0;
}
//---------------------------------------------------------------

int
readDecompressBlock(int fd, TLBHEADER *thp, TLBDATA *tdp, int seqOnly)
{
	unsigned int size = thp->blockLength - sizeof(TLBHEADER);
	size_t rs = size;
	assert(size <= kMaxBlockLength);
	void *dataBlock = tdp->diskBuffer; assert(dataBlock);
	void *rp = dataBlock;
	while (1)
	{
		ssize_t res = read(fd,rp,rs);
		if (res == rs)
			break;
		if (res > 0)
		{
			rp += res;
			rs -= res;
			continue;
		}
		if (res == 0)
		{
			fprintf(stderr,"unexpected EOF while reading block in readDecompressBlock\n");
			return 1;
		}
		if (res == -1 && (errno == EINTR || errno == EAGAIN))
			continue;
		perror("readDecompressBlock read:");
		return 1;
	}
	int err = decompressBlock(dataBlock,size,thp,tdp,seqOnly);
	return err;
}
//---------------------------------------------------------------

int
preadDecompressBlock(int fd, off_t offset, TLBHEADER *thp, TLBDATA *tdp, int seqOnly)
{
	unsigned int size = thp->blockLength - sizeof(TLBHEADER);
	size_t rs = size;
	assert(size <= kMaxBlockLength);
	void *dataBlock = tdp->diskBuffer; assert(dataBlock);
	void *rp = dataBlock;
	while (1)
	{
		ssize_t res = pread(fd,rp,rs,offset);
		if (res == rs)
			break;
		if (res > 0)
		{
			rp += res;
			rs -= res;
			offset += (off_t) res;
			continue;
		}
		if (res == 0)
		{
			fprintf(stderr,"unexpected EOF while reading block in readDecompressBlock\n");
			return 1;
		}
		if (res == -1 && (errno == EINTR || errno == EAGAIN))
			continue;
		perror("readDecompressBlock read:");
		return 1;
	}
	return decompressBlock(dataBlock,size,thp,tdp,seqOnly);
}
//---------------------------------------------------------------

void
freeBlock(TLBDATA *tdp)
{
	free(tdp->pd);
	free(tdp->pdp);
	free(tdp->readLen.data);
	for (unsigned int i = 0; i < kMAXnbMatches; i++)
	{
		free(tdp->pairPos[i].data);
		free(tdp->mismatch[i].data);
		free(tdp->cigarBlock[i].data);
	}
	free(tdp->unmapped.data);
	free(tdp->readHdr.data);
	free(tdp->readQual.data);
	if(tdp->diskBuffer) free(tdp->diskBuffer);
}
//---------------------------------------------------------------

int
decodeBlock(void *dataBlock, unsigned int size, TLBHEADER *thp, PAIRDATA *dest, TLBDATA *tdp)
{
	void (*decode)(TLBDATA *, unsigned int, PAIRDATA *) = decodePair;
	int (*addDiff)(PAIRDATA *, unsigned char *, int, unsigned int) = addDiffPair;
	int (*addCigar)(PAIRDATA *, unsigned char *, int, unsigned int) = addCigarPair;
	if (decompressBlock(dataBlock,size,thp,tdp,0) != 0)
		return -1;
	if ((thp->flags_7_0 & kTLBHflagPairedReads) == 0)
	{
		decode = decodeSingle;
		addDiff = addDiffSingle;
		addCigar = addCigarSingle;
	}
	if (tdp->cigar[0] != NULL)
	{
		int mmpos[kMAXnbMatches] = { 0 };
		int cigarpos[kMAXnbMatches] = { 0 };
		unsigned int i;
		for (i=0; i<tdp->cnt; i++)
		{
			decode(tdp,i,dest + i);
			for (unsigned int j = 0; j < tdp->nbMatches; j++)
			{
				cigarpos[j] = addCigar(dest + i,tdp->cigar[j],cigarpos[j],j);
				mmpos[j] = addDiff(dest + i,tdp->diff[j],mmpos[j],j);
			}
		}
	}
	else if (tdp->diff[0] != NULL)
	{
		int mmpos[kMAXnbMatches] = { 0 };
		unsigned int i;
		for (i=0; i<tdp->cnt; i++)
		{
			decode(tdp,i,dest + i);
			for (unsigned int j = 0; j < tdp->nbMatches; j++)
				mmpos[j] = addDiff(dest + i,tdp->diff[j],mmpos[j],j);
		}
	}
	else
	{
		unsigned int i;
		for (i=0; i<tdp->cnt; i++)
		{
			decode(tdp,i,dest + i);
		}
	}
	return tdp->cnt;
}
//---------------------------------------------------------------

int
readDecodeBlock(int fd, TLBHEADER *thp, PAIRDATA *dest, TLBDATA *tdp)
{
	unsigned int size = thp->blockLength - sizeof(TLBHEADER);
	assert(size <= kMaxBlockLength);
	void *dataBlock = tdp->diskBuffer; assert(dataBlock);
	ssize_t res = read(fd,dataBlock,size);
	if (res != size)
		return -1;
	int cnt = decodeBlock(dataBlock,size,thp,dest,tdp);
	return cnt;
}
//---------------------------------------------------------------

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
