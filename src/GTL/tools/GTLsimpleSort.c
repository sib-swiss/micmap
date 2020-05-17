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
/*

  This software will sort a GTL file in RAM and write sorted pieces

	(c) N.Guex and C.Iseli 2020

*/

//=================================================================================================
#define _GNU_SOURCE
#include "config.h"
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <sys/time.h>

#include <sys/io.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <errno.h>
#include <assert.h>

#include "GTL.h"
#include "gtlVersion.h"

// number of possibilities in the 8nt head of each unmapped read, given that the 1st nt is known
#define NPOSSIBILITIES 0x4000
// factor to multiply the range between 1% and 99% of the sorted 8 nt head of each unmapped read
#define IQR99_MULT_FACTOR 100

#define ONE_GIGA (1024*1024*1024UL)

//=================================================================================================

static unsigned int chrMax = 25;

static void
writeBuffer(TLBDATA *block, unsigned int cnt, PAIRDATA **pdp, const char* ofprefix, const char *partName, TLBDATA *wbp)
{
	unsigned int i;
	if (cnt == 0)
		return;
	zeroBlock(wbp);
	memcpy(&wbp->header,&block->header,sizeof(TLBHEADER));
	wbp->cnt = cnt;
	// Get the various info from the block components
	for (i = 0; i < cnt; i++)
	{
		if (i == 0)
		{
			wbp->lengthInfo = pdp[i]->len1;
			wbp->nbMatches = pdp[i]->nbMatches;
			assert(pdp[i]->nbMatches != 0);
		}
		assert(wbp->nbMatches == pdp[i]->nbMatches);
		if (pdp[i]->len1 != wbp->lengthInfo || pdp[i]->len2 != wbp->lengthInfo)
			wbp->lengthInfo = 0;
		if ((block->header.flags_7_0 & kTLBHflagPairPosBlock) != 0
				|| (block->header.flags_7_0 & kTLBHflagPairChrPosBlock) != 0)
		{
			unsigned int minp;
			if ((block->header.flags_7_0 & kTLBHflagPairedReads) == 0)
			{
				minp = pdp[i]->pcpd[0].tag1pos;
			}
			else
			{
				if (pdp[i]->pcpd[0].tag1pos != 0)
				{
					if (pdp[i]->pcpd[0].tag2pos != 0)
					{
						if (pdp[i]->pcpd[0].tag1pos <= pdp[i]->pcpd[0].tag2pos)
							minp = pdp[i]->pcpd[0].tag1pos;
						else
							minp = pdp[i]->pcpd[0].tag2pos;
					}
					else
						minp = pdp[i]->pcpd[0].tag1pos;
				}
				else
					minp = pdp[i]->pcpd[0].tag2pos;
			}
			if (wbp->header.minPos > minp)
				wbp->header.minPos = minp;
			if (wbp->header.maxPos < minp)
				wbp->header.maxPos = minp;
		}
		wbp->hdrSize += pdp[i]->hlen1;
		wbp->hdrSize += pdp[i]->hlen2;
		wbp->qualSize += pdp[i]->len1;
		wbp->qualSize += pdp[i]->len2;
		for (unsigned int j = 0; j < kMAXnbMatches; j++)
		{
			wbp->mmSize[j] += pdp[i]->ml[j];
			wbp->cigarSize[j] += pdp[i]->cl[j];
		}
		wbp->p4Size += pdp[i]->pl;
	}
	wbp->pdp = pdp;
	unsigned char *dataBlock = packBlock(wbp, wbp);
	char fname[1024];
	if ((wbp->header.chr & 0xffff) == kUnmappedChrID)
	{
		if ((wbp->header.flags_7_0 & kTLBHflagPairedReads) == 0)
		{
			unsigned int c = (wbp->header.chr >> 16) & 3;
			snprintf(fname,1024,"%s/UM_%c%s.gtl",ofprefix,"ACGT"[c],partName);
		}
		else
		{
			unsigned int c1 = (wbp->header.chr >> 18) & 3;
			unsigned int c2 = (wbp->header.chr >> 16) & 3;
			snprintf(fname,1024,"%s/UM_%c%c%s.gtl",ofprefix,"ACGT"[c1],"ACGT"[c2],partName);
		}
	}
	else if ((wbp->header.chr & 0xffff) == kHalfMapChrID)
	{
		unsigned int chr = wbp->header.chr >> 16;
		snprintf(fname,1024,"%s/HM_chr%d.gtl",ofprefix,chr);
	}
	else if ((wbp->header.chr & 0xffff) == kChimeraChrID)
	{
		unsigned int chr = wbp->header.chr >> 16;
		if (chr > chrMax)
			snprintf(fname,1024,"%s/C_d%d.gtl",ofprefix,chr - chrMax);
		else
			snprintf(fname,1024,"%s/C_chr%d.gtl",ofprefix,chr);
	}
	else
	{
		if ((wbp->header.flags_15_8 & kTLBHflagMultipleMatches) != 0)
			snprintf(fname,1024,"%s/M_chr%d%c_%s.gtl",ofprefix,wbp->header.chr,'a' + wbp->header.minPos / kCountsPerChrPart,partName);
		else
			snprintf(fname,1024,"%s/chr%d%c_%s.gtl",ofprefix,wbp->header.chr,'a' + wbp->header.minPos / kCountsPerChrPart,partName);
	}
	int fd = open(fname, O_WRONLY|O_CREAT|O_APPEND, 0664);
	if (fd == -1)
	{
		perror("open:");
		exit(1);
	}
	if (write(fd,dataBlock,wbp->header.blockLength) != wbp->header.blockLength)
	{
		perror("write failed");
		exit(1);
	}
	else
	{
		snprintf(fname,1024,"%s/chr%d_i.txt",ofprefix,wbp->header.chr & 0xffff);
		FILE *f = fopen(fname,"a");
		fprintf(f,"GTLINDEX\t%d\t%s\t%u\t%u\t%u\n",wbp->header.chr & 0xffff,partName,wbp->header.minPos,wbp->header.maxPos,wbp->header.blockLength);
		fclose(f);
	}
	close(fd);
}
//---------------------------------------------------------------

int
intcompare(const void *a, const void *b)
{
	unsigned int ia = *((unsigned int *)a);
	unsigned int ib = *((unsigned int *)b);
	if (ia < ib)
		return -1;
	if (ia > ib)
		return 1;
	return 0;
}
//---------------------------------------------------------------

static int
sort(char *fn, unsigned int chr, const char* ofprefix, int multi)
{
	int fd = open(fn, O_RDONLY);
	struct stat sb;
	if (fstat(fd,&sb) != 0)
	{
		perror("fstat");
		return 1;
	}
	unsigned char *readBuf = NULL;
	size_t fSize = 0;
	size_t rCnt;
	unsigned char *rPos;
	if (S_ISREG(sb.st_mode))
	{
		readBuf = malloc(sb.st_size);
		if (readBuf == NULL)
			return 1;
		fSize = sb.st_size;
		rCnt = sb.st_size;
		rPos = readBuf;
		ssize_t res;
		while ((res = read(fd,rPos,rCnt)) != rCnt)
		{
			if (res == 0)
			{
				fprintf(stderr,"we miss %ld bits (got 0)\n", rCnt);
				return 1;
			}
			if (res == -1)
			{
				if (errno != EINTR)
				{
					perror("read");
					return 1;
				}
			}
			else
			{
				rPos += res;
				rCnt -= res;
			}
		}
	}
	else
	{
		unsigned int nb_giga = 1;
		readBuf = malloc(ONE_GIGA);
		if (readBuf == NULL)
		{
			fprintf(stderr, "Could not malloc %ld: %s(%d)\n", ONE_GIGA, strerror(errno), errno);
			return 1;
		}
		rCnt = ONE_GIGA;
		rPos = readBuf;
		ssize_t res;
		while ((res = read(fd,rPos,rCnt)) != 0)
		{
			if (res == -1)
			{
				if (errno != EINTR)
				{
					perror("read");
					return 1;
				}
			}
			else
			{
				rPos += res;
				fSize += res;
				rCnt -= res;
			}
			if (rCnt == 0)
			{
				// need more space
				readBuf = realloc(readBuf, (nb_giga + 1) * ONE_GIGA);
				if (readBuf == NULL)
				{
					fprintf(stderr, "Could not realloc %ld: %s(%d)\n", nb_giga * ONE_GIGA, strerror(errno), errno);
					return 1;
				}
				rPos = readBuf + nb_giga * ONE_GIGA;
				nb_giga += 1;
				rCnt = ONE_GIGA;
			}
		}
	}
	close(fd);
	// count the blocks
	unsigned int allBlocks = 0;
	unsigned int perfectBlocks = 0;
	unsigned int mismatchNBlocks = 0;
	unsigned int mismatchBlocks = 0;
	unsigned int alignBlocks = 0;
	unsigned int unmapBlocks = 0;
	unsigned int halfmapBlocks = 0;
	unsigned int chimeraBlocks = 0;
	// keep track if we deal with paired reads
	unsigned int pairedFlag = 0;
	rPos = readBuf;
	rCnt = 0;
	while (rCnt < fSize)
	{
		unsigned char hcs = simple8bitCS(rPos, sizeof(TLBHEADER));
		TLBHEADER *th = (TLBHEADER *) rPos;
		if (hcs != 0)
		{
			fprintf(stderr,"Checksum mismatch: %02x %02x\n",th->headerCS,hcs);
			return 1;
		}
		if (th->chr == chr && (th->flags_15_8 & kTLBHflagMultipleMatches) == multi)
		{
			// we want this block
			allBlocks += 1;
			if ((chr & 0xffff) == kUnmappedChrID)
				unmapBlocks += 1;
			else if ((chr & 0xffff) == kHalfMapChrID)
				halfmapBlocks += 1;
			else if ((chr & 0xffff) == kChimeraChrID)
				chimeraBlocks += 1;
			else if ((th->flags_7_0 & kTLBHflagPairPosBlock) != 0)
			{
				if ((th->flags_7_0 & 0xf0) == 0)
					perfectBlocks += 1;
				if ((th->flags_7_0 & 0xf0) == kTLBHflagNMismatchBlock)
					mismatchNBlocks += 1;
				if ((th->flags_7_0 & 0xf0) == kTLBHflagMismatchBlock)
					mismatchBlocks += 1;
				if ((th->flags_7_0 & 0xf0) == (kTLBHflagMismatchBlock | kTLBHflagCigarBlock))
					alignBlocks += 1;
			}
			if ((th->flags_7_0 & kTLBHflagPairedReads) != 0)
				pairedFlag = kTLBHflagPairedReads;
		}
		rPos += th->blockLength;
		rCnt += th->blockLength;
	}
	unsigned int *patCnt1 = NULL, *patCnt2 = NULL;
	unsigned char *patOut1 = NULL, *patOut2 = NULL;
	char patSuffix[16 * 0x100];
	unsigned int cutoff1 = 0, cutoff2 = 0;
	unsigned int nbFiltered = 0;
	unsigned int filterBlocks = 0;
	if ((chr & 0xffff) == kUnmappedChrID)
	{
		// flag reads probably issued from vector
		unsigned int i;
		unsigned int totalPairs = 0;
		patCnt1 = calloc(NPOSSIBILITIES,sizeof(unsigned int));
		patOut1 = malloc(NPOSSIBILITIES * sizeof(unsigned char));
		if (pairedFlag)
		{
			patCnt2 = calloc(NPOSSIBILITIES,sizeof(unsigned int));
			patOut2 = malloc(NPOSSIBILITIES * sizeof(unsigned char));
		}
		rPos = readBuf;
		rCnt = 0;
		while (rCnt < fSize)
		{
			// we already checked the header CS above, so just proceed...
			TLBHEADER *th = (TLBHEADER *) rPos;
			if (th->chr == chr && (th->flags_15_8 & kTLBHflagMultipleMatches) == multi)
			{
				// we want this block
				unsigned int size = th->blockLength - sizeof(TLBHEADER);
				totalPairs += countPat1InUnmappedBlock(rPos + sizeof(TLBHEADER),size,th,patCnt1);
				fprintf(stderr,".");
			}
			rPos += th->blockLength;
			rCnt += th->blockLength;
		}
		fprintf(stderr,"\n");
		unsigned int *sa = malloc(NPOSSIBILITIES * sizeof(unsigned int));
		memcpy(sa,patCnt1,NPOSSIBILITIES * sizeof(unsigned int));
		qsort(sa,NPOSSIBILITIES,sizeof(unsigned int),intcompare);
		cutoff1 = sa[NPOSSIBILITIES/2] + IQR99_MULT_FACTOR * (sa[NPOSSIBILITIES*99/100] - sa[NPOSSIBILITIES/100]);
		if (cutoff1 == 0)
			cutoff1 = sa[NPOSSIBILITIES-1] + 1;
		for (i = 0; i < NPOSSIBILITIES; i++)
		{
			if (patCnt1[i] > cutoff1)
			{
				nbFiltered += patCnt1[i];
				if (filterBlocks >= 0x100)
				{
					fprintf(stderr,"HWHAP: too many patterns to filter at cutoff1 %u sa[0] %u sa[NPOSSIBILITIES/2] %u sa[NPOSSIBILITIES-1] %u\n",cutoff1,sa[0],sa[NPOSSIBILITIES/2],sa[NPOSSIBILITIES-1]);
					return 1;
				}
				snprintf(patSuffix + 16 * filterBlocks,16,"_1%c%c%c%c%c%c%c%c",
					"ACGT"[(chr >> (pairedFlag ? 18 : 16)) & 3],
					"ACGT"[(i >> 12) & 3],
					"ACGT"[(i >> 10) & 3],
					"ACGT"[(i >> 8) & 3],
					"ACGT"[(i >> 6) & 3],
					"ACGT"[(i >> 4) & 3],
					"ACGT"[(i >> 2) & 3],
					"ACGT"[i & 3]);
				patOut1[i] = filterBlocks++;
			}
		}
		if (pairedFlag)
		{
			rPos = readBuf;
			rCnt = 0;
			while (rCnt < fSize)
			{
				// we already checked the header CS above, so just proceed...
				TLBHEADER *th = (TLBHEADER *) rPos;
				if (th->chr == chr && (th->flags_15_8 & kTLBHflagMultipleMatches) == multi)
				{
					// we want this block
					unsigned int size = th->blockLength - sizeof(TLBHEADER);
					countPat2InUnmappedBlock(rPos + sizeof(TLBHEADER),size,th,patCnt1,cutoff1,patCnt2);
					fprintf(stderr,".");
				}
				rPos += th->blockLength;
				rCnt += th->blockLength;
			}
			fprintf(stderr,"\n");
			memcpy(sa,patCnt2,NPOSSIBILITIES * sizeof(unsigned int));
			qsort(sa,NPOSSIBILITIES,sizeof(unsigned int),intcompare);
			cutoff2 = sa[NPOSSIBILITIES/2] + IQR99_MULT_FACTOR * (sa[NPOSSIBILITIES*99/100] - sa[NPOSSIBILITIES/100]);
			if (cutoff2 == 0)
				cutoff2 = sa[NPOSSIBILITIES-1] + 1;
			for (i = 0; i < NPOSSIBILITIES; i++)
			{
				if (patCnt2[i] > cutoff2)
				{
					nbFiltered += patCnt2[i];
					if (filterBlocks >= 0x100)
					{
						fprintf(stderr,"HWHAP: too many patterns to filter at cutoff2 %u sa[0] %u sa[NPOSSIBILITIES/2] %u sa[NPOSSIBILITIES-1] %u\n",cutoff2,sa[0],sa[NPOSSIBILITIES/2],sa[NPOSSIBILITIES-1]);
						return 1;
					}
					snprintf(patSuffix + 16 * filterBlocks,16,"_2%c%c%c%c%c%c%c%c",
						"ACGT"[(chr >> 16) & 3],
						"ACGT"[(i >> 12) & 3],
						"ACGT"[(i >> 10) & 3],
						"ACGT"[(i >> 8) & 3],
						"ACGT"[(i >> 6) & 3],
						"ACGT"[(i >> 4) & 3],
						"ACGT"[(i >> 2) & 3],
						"ACGT"[i & 3]);
					patOut2[i] = filterBlocks++;
				}
			}
		}
		free(sa);
		printf("we will filter %u out of %u (%.2f %%) and use %u filterBlocks\n",nbFiltered,totalPairs,100.0*nbFiltered/totalPairs,filterBlocks);
		// we will need less blocks to sort, but add a few blocks to output the filtered stuff
		// tricky to get right, might temporarily overflow -  allBlocks -= nbFiltered / kBlockMaxCount;
		unmapBlocks -= nbFiltered / kBlockMaxCount;
	}
	PAIRDATA *allData = malloc(kBlockMaxCount * allBlocks * sizeof(PAIRDATA));
	if (allData == NULL)
		return 1;
	unsigned int allCnt = 0;
	PAIRDATA **perfect = malloc(kBlockMaxCount * perfectBlocks * sizeof(PAIRDATA *));
	if (perfect == NULL)
		return 1;
	unsigned int perfectCnt = 0;
	PAIRDATA **mismatchN = malloc(kBlockMaxCount * mismatchNBlocks * sizeof(PAIRDATA *));
	if (mismatchN == NULL)
		return 1;
	unsigned int mismatchNCnt = 0;
	PAIRDATA **mismatch = malloc(kBlockMaxCount * mismatchBlocks * sizeof(PAIRDATA *));
	if (mismatch == NULL)
		return 1;
	unsigned int mismatchCnt = 0;
	PAIRDATA **align = malloc(kBlockMaxCount * alignBlocks * sizeof(PAIRDATA *));
	if (align == NULL)
		return 1;
	unsigned int alignCnt = 0;
	PAIRDATA **unmap = malloc(kBlockMaxCount * unmapBlocks * sizeof(PAIRDATA *));
	if (unmap == NULL)
		return 1;
	unsigned int unmapCnt = 0;
	PAIRDATA **halfmap = malloc(kBlockMaxCount * halfmapBlocks * sizeof(PAIRDATA *));
	if (halfmap == NULL)
		return 1;
	unsigned int halfmapCnt = 0;
	PAIRDATA **chimera = malloc(kBlockMaxCount * chimeraBlocks * sizeof(PAIRDATA *));
	if (chimera == NULL)
		return 1;
	unsigned int chimeraCnt = 0;
	unsigned int i;
	const unsigned int nbBlocks = 7;
	TLBDATA *allBlock = (TLBDATA *) calloc(nbBlocks, sizeof(TLBDATA));
	for (i = 0; i < nbBlocks; i++)
	{
		allBlock[i].header.magic[0] = 'G';
		allBlock[i].header.magic[1] = 'T';
		allBlock[i].header.magic[2] = 'L';
		allBlock[i].header.chr = chr;
		allBlock[i].header.minPos = 0xffffffffU;
	}
	int cb = 0;
	TLBDATA *perfectBlock = allBlock + cb++;
	perfectBlock->header.flags_7_0 |= pairedFlag;
	perfectBlock->header.flags_7_0 |= kTLBHflagPairPosBlock;
	perfectBlock->header.flags_15_8 |= multi;
	TLBDATA *mismatchNBlock = allBlock + cb++;
	mismatchNBlock->header.flags_7_0 |= pairedFlag;
	mismatchNBlock->header.flags_7_0 |= kTLBHflagPairPosBlock;
	mismatchNBlock->header.flags_7_0 |= kTLBHflagNMismatchBlock;
	mismatchNBlock->header.flags_15_8 |= multi;
	TLBDATA *mismatchBlock = allBlock + cb++;
	mismatchBlock->header.flags_7_0 |= pairedFlag;
	mismatchBlock->header.flags_7_0 |= kTLBHflagPairPosBlock;
	mismatchBlock->header.flags_7_0 |= kTLBHflagMismatchBlock;
	mismatchBlock->header.flags_15_8 |= multi;
	TLBDATA *alignBlock = allBlock + cb++;
	alignBlock->header.flags_7_0 |= pairedFlag;
	alignBlock->header.flags_7_0 |= kTLBHflagPairPosBlock;
	alignBlock->header.flags_7_0 |= kTLBHflagMismatchBlock;
	alignBlock->header.flags_7_0 |= kTLBHflagCigarBlock;
	alignBlock->header.flags_15_8 |= multi;
	TLBDATA *unmapBlock = allBlock + cb++;
	unmapBlock->header.flags_7_0 |= pairedFlag;
	unmapBlock->header.flags_7_0 |= kTLBHflagMismatchBlock;
	unmapBlock->header.flags_15_8 |= kTLBHflagUnmappedBlock;
	unmapBlock->header.minPos = 0;
	TLBDATA *halfmapBlock = allBlock + cb++;
	halfmapBlock->header.flags_7_0 |= pairedFlag;
	halfmapBlock->header.flags_15_8 |= kTLBHflagHalfmapBlock;
	halfmapBlock->header.flags_7_0 |= kTLBHflagPairChrPosBlock;
	halfmapBlock->header.flags_7_0 |= kTLBHflagMismatchBlock;
	halfmapBlock->header.flags_15_8 |= multi;
	TLBDATA *chimeraBlock = allBlock + cb++;
	chimeraBlock->header.flags_7_0 |= pairedFlag;
	chimeraBlock->header.flags_15_8 |= kTLBHflagChimeraBlock;
	chimeraBlock->header.flags_7_0 |= kTLBHflagPairChrPosBlock;
	chimeraBlock->header.flags_7_0 |= kTLBHflagMismatchBlock;
	chimeraBlock->header.flags_15_8 |= multi;
	TLBDATA *filterBlock = (TLBDATA *) calloc(filterBlocks, sizeof(TLBDATA));
	for (i = 0; i < filterBlocks; i++)
	{
		filterBlock[i].header.magic[0] = 'G';
		filterBlock[i].header.magic[1] = 'T';
		filterBlock[i].header.magic[2] = 'L';
		filterBlock[i].header.chr = chr;
		filterBlock[i].header.flags_7_0 |= pairedFlag;
		filterBlock[i].header.flags_7_0 |= kTLBHflagMismatchBlock;
		filterBlock[i].header.flags_15_8 |= kTLBHflagUnmappedBlock;
		filterBlock[i].pd = malloc(kBlockMaxCount * sizeof(PAIRDATA));
		filterBlock[i].pdp = malloc(kBlockMaxCount * sizeof(PAIRDATA *));
	}
	rPos = readBuf;
	rCnt = 0;
	TLBDATA td;
	allocBlock(&td, true);
	free(td.pdp); // won't be used
	td.pdp = NULL;
	while (rCnt < fSize)
	{
		// we already checked the header CS above, so just proceed...
		TLBHEADER *th = (TLBHEADER *) rPos;
		if (th->chr == chr && (th->flags_15_8 & kTLBHflagMultipleMatches) == multi)
		{
			// we want this block
			unsigned int size = th->blockLength - sizeof(TLBHEADER);
			int cnt = decodeBlock(rPos + sizeof(TLBHEADER),size,th,allData + allCnt,&td);
			if (cnt < 0)
				return 1;
			PAIRDATA **selected = NULL;
			if ((chr & 0xffff) == kUnmappedChrID)
			{
				selected = unmap + unmapCnt;
				// filter reads
				unsigned int good = 0;
				for (i = 0; i < cnt; i++)
				{
					PAIRDATA *dbp = allData + allCnt + i;
					unsigned int pat1 = ((dbp->packed4nt[0] & 0x3f) << 8) | dbp->packed4nt[1];
					int filter = patCnt1[pat1] > cutoff1;
					unsigned int outIdx;
					if (filter)
						outIdx = patOut1[pat1];
					else if (pairedFlag)
					{
						unsigned int pl = dbp->len1 >> 2;
						if ((dbp->len1 & 3) != 0)
							pl += 1;
						unsigned int pat2 = ((dbp->packed4nt[pl] & 0x3f) << 8) | dbp->packed4nt[pl + 1];
						filter = patCnt2[pat2] > cutoff2;
						if (filter)
							outIdx = patOut2[pat2];
					}
					if (filter)
					{
						filterBlock[outIdx].pd[filterBlock[outIdx].cnt] = allData[allCnt + i];
						filterBlock[outIdx].pdp[filterBlock[outIdx].cnt] = filterBlock[outIdx].pd + filterBlock[outIdx].cnt;
						filterBlock[outIdx].cnt += 1;
						if (filterBlock[outIdx].cnt == kBlockMaxCount)
						{
							qsort(filterBlock[outIdx].pdp,kBlockMaxCount,sizeof(PAIRDATA *),compareUNMAPDATAptr);
							writeBuffer(filterBlock + outIdx,kBlockMaxCount,filterBlock[outIdx].pdp,ofprefix,patSuffix + 16 * outIdx,&td);
							fprintf(stderr,"f");
							filterBlock[outIdx].cnt = 0;
						}
					}
					else
					{
						if (good != i)
							// need to move this
							allData[allCnt + good] = allData[allCnt + i];
						good += 1;
					}
				}
				cnt = good;
				unmapCnt += cnt;
				fprintf(stderr,"u");
			}
			else if ((chr & 0xffff) == kHalfMapChrID)
			{
				selected = halfmap + halfmapCnt;
				halfmapCnt += cnt;
				fprintf(stderr,"h");
			}
			else if ((chr & 0xffff) == kChimeraChrID)
			{
				selected = chimera + chimeraCnt;
				chimeraCnt += cnt;
				fprintf(stderr,"c");
			}
			else if ((th->flags_7_0 & kTLBHflagPairPosBlock) != 0)
			{
				if ((th->flags_7_0 & 0xf0) == 0)
				{
					selected = perfect + perfectCnt;
					perfectCnt += cnt;
					fprintf(stderr,"p");
				}
				if ((th->flags_7_0 & 0xf0) == kTLBHflagNMismatchBlock)
				{
					selected = mismatchN + mismatchNCnt;
					mismatchNCnt += cnt;
					fprintf(stderr,"N");
				}
				if ((th->flags_7_0 & 0xf0) == kTLBHflagMismatchBlock)
				{
					selected = mismatch + mismatchCnt;
					mismatchCnt += cnt;
					fprintf(stderr,"m");
				}
				if ((th->flags_7_0 & 0xf0) == (kTLBHflagMismatchBlock | kTLBHflagCigarBlock))
				{
					selected = align + alignCnt;
					alignCnt += cnt;
					fprintf(stderr,"a");
				}
			}
			int i;
			for (i = 0; i < cnt; i++)
			{
				selected[i] = allData + allCnt + i;
			}
			allCnt += cnt;
		}
		rPos += th->blockLength;
		rCnt += th->blockLength;
	}
	for (i = 0; i < filterBlocks; i++)
	{
		if (filterBlock[i].cnt != 0)
		{
			qsort(filterBlock[i].pdp,filterBlock[i].cnt,sizeof(PAIRDATA *),compareUNMAPDATAptr);
			writeBuffer(filterBlock + i,filterBlock[i].cnt,filterBlock[i].pdp,ofprefix,patSuffix + 16 * i,&td);
			fprintf(stderr,"f");
			free(filterBlock[i].pd);
			free(filterBlock[i].pdp);
		}
	}
	free(filterBlock);
	filterBlock = NULL;
	if ((chr & 0xffff) == kUnmappedChrID)
	{
		fprintf(stderr,"\nsorting %u unmapped\n",unmapCnt);
		qsort(unmap,unmapCnt,sizeof(PAIRDATA *),compareUNMAPDATAptr);
	}
	else if ((chr & 0xffff) == kHalfMapChrID)
	{
		fprintf(stderr,"\nsorting %u halfmapped\n",halfmapCnt);
		qsort(halfmap,halfmapCnt,sizeof(PAIRDATA *),comparePAIRUDATAptr);
	}
	else if ((chr & 0xffff) == kChimeraChrID)
	{
		fprintf(stderr,"\nsorting %u chimera\n",chimeraCnt);
		qsort(chimera,chimeraCnt,sizeof(PAIRDATA *),comparePAIRDATAptr);
	}
	else
	{
		fprintf(stderr,"\nsorting %u perfect %u mismatchN %u mismatch %u align\n",perfectCnt,mismatchNCnt,mismatchCnt,alignCnt);
		qsort(perfect,perfectCnt,sizeof(PAIRDATA *),comparePAIRDATAptr);
		qsort(mismatchN,mismatchNCnt,sizeof(PAIRDATA *),comparePAIRDATAptr);
		qsort(mismatch,mismatchCnt,sizeof(PAIRDATA *),comparePAIRDATAptr);
		qsort(align,alignCnt,sizeof(PAIRDATA *),comparePAIRDATAptr);
	}
	unsigned int cur = 0;
	unsigned int chunk = kBlockMaxCount;
	while (cur + chunk <= perfectCnt)
	{
		writeBuffer(perfectBlock,chunk,perfect + cur,ofprefix,"p",&td);
		fprintf(stderr,"p");
		cur += chunk;
	}
	chunk = perfectCnt - cur;
	if (chunk > 0)
	{
		writeBuffer(perfectBlock,chunk,perfect + cur,ofprefix,"p",&td);
		fprintf(stderr,"p");
	}
	cur = 0;
	chunk = kBlockMaxCount;
	while (cur + chunk <= mismatchNCnt)
	{
		writeBuffer(mismatchNBlock,chunk,mismatchN + cur,ofprefix,"N",&td);
		fprintf(stderr,"N");
		cur += chunk;
	}
	chunk = mismatchNCnt - cur;
	if (chunk > 0)
	{
		writeBuffer(mismatchNBlock,chunk,mismatchN + cur,ofprefix,"N",&td);
		fprintf(stderr,"N");
	}
	cur = 0;
	chunk = kBlockMaxCount;
	while (cur + chunk <= mismatchCnt)
	{
		writeBuffer(mismatchBlock,chunk,mismatch + cur,ofprefix,"m",&td);
		fprintf(stderr,"m");
		cur += chunk;
	}
	chunk = mismatchCnt - cur;
	if (chunk > 0)
	{
		writeBuffer(mismatchBlock,chunk,mismatch + cur,ofprefix,"m",&td);
		fprintf(stderr,"m");
	}
	cur = 0;
	chunk = kBlockMaxCount;
	while (cur + chunk <= alignCnt)
	{
		writeBuffer(alignBlock,chunk,align + cur,ofprefix,"g",&td);
		fprintf(stderr,"a");
		cur += chunk;
	}
	chunk = alignCnt - cur;
	if (chunk > 0)
	{
		writeBuffer(alignBlock,chunk,align + cur,ofprefix,"g",&td);
		fprintf(stderr,"a");
	}
	cur = 0;
	chunk = kBlockMaxCount;
	while (cur + chunk <= unmapCnt)
	{
		writeBuffer(unmapBlock,chunk,unmap + cur,ofprefix,"",&td);
		fprintf(stderr,"u");
		cur += chunk;
	}
	chunk = unmapCnt - cur;
	if (chunk > 0)
	{
		writeBuffer(unmapBlock,chunk,unmap + cur,ofprefix,"",&td);
		fprintf(stderr,"u");
	}
	cur = 0;
	chunk = kBlockMaxCount;
	while (cur + chunk <= halfmapCnt)
	{
		writeBuffer(halfmapBlock,chunk,halfmap + cur,ofprefix,"h",&td);
		fprintf(stderr,"h");
		cur += chunk;
	}
	chunk = halfmapCnt - cur;
	if (chunk > 0)
	{
		writeBuffer(halfmapBlock,chunk,halfmap + cur,ofprefix,"h",&td);
		fprintf(stderr,"h");
	}
	cur = 0;
	chunk = kBlockMaxCount;
	while (cur + chunk <= chimeraCnt)
	{
		writeBuffer(chimeraBlock,chunk,chimera + cur,ofprefix,"c",&td);
		fprintf(stderr,"c");
		cur += chunk;
	}
	chunk = chimeraCnt - cur;
	if (chunk > 0)
	{
		writeBuffer(chimeraBlock,chunk,chimera + cur,ofprefix,"c",&td);
		fprintf(stderr,"c");
	}
	fprintf(stderr,"\n");

	td.pdp = NULL; // was freed already
	freeBlock(&td);

	return 0;
} // sort
//---------------------------------------------------------------

int
main (int argc, char **argv)
{
	int c;
	char rsltfile[512];
	char ofprefix[512] = "/tmp/sorted_chr";
	int err;
	unsigned int chr = 0xffffffffU;
	int multi = 0;

	/* --------- process arguments ----------------------*/

	rsltfile[0] = 0;
	opterr = 0;
	while ((c = getopt (argc, argv, "c:m:o:r:s:u:z:M")) != -1)
	switch (c)
	{
	  case 'c':
			sscanf(optarg,"%u",&chr);
			break;

	  case 'm':
			sscanf(optarg,"%u",&chrMax);
			break;

		case 'o':
			sscanf(optarg,"%511s",ofprefix);
			break;

	  case 'r':
			strcpy(rsltfile,optarg);
			break;

		case 's':
			gCompress.scoredef = compress_find(optarg);
			break;

	  case 'u':
			sscanf(optarg,"%u",&gHasUMI);
			break;

		case 'z':
			gCompress.def = compress_find(optarg);
			break;

	  case 'M':
			multi = kTLBHflagMultipleMatches;
			break;
	}

	if (rsltfile[0] == 0 || gCompress.scoredef < 0 || gCompress.def < 0)
	{
		unsigned int i;
		printf("usage:\n\n");
        printf("simpleSort -c chromosome -r resultfile [-o prefix]\n\n");
		printf("           -c chromosome         : select chromosome to sort (1..25 + special cases for chimera, half- and un-mapped)\n");
		printf("           -m chrMax             : define max chr number of reference (used for delta chimera grouping), default is %u\n",chrMax);
		printf("           -o prefix             : prefix to prepend to all output files, default is '%s'\n", ofprefix);
		printf("           -z compression        : compression to use for data, default is '%s'\n", compress_name(gCompress.def));
		printf("           -s compression        : compression to use for scores, default is '%s'\n", compress_name(gCompress.scoredef));
		// TODO add an option either generic "level" or per compressor flag (like "hc" for lz4)
		printf("                    available compression :\n");
		for (i = 0; compress_dict[i].name != NULL; i++)
			printf("                     %11.11s : %s\n", compress_dict[i].name, compress_dict[i].description);
		printf("           -r resultfile         : name of the result file (produced by match)\n");
		printf("           -M                    : select multiple matches\n");
		printf("           -u <int>              : reads have UMI of that length in their header, 0 means no UMI, default is %d\n",gHasUMI);
		printf("\n");
		printf(GTL_VERSION_FULL_STRING "\n");
		printf("(c) Nicolas Guex & Christian Iseli & Thierry Schuepbach & Ivan Topolsky 2015-2020\n");
		return(1);
	}

	err = sort(rsltfile,chr,ofprefix,multi);
	report_compress_stat(stderr);

	return(err);

} /* main */
//---------------------------------------------------------------

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
