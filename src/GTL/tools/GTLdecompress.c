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
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <assert.h>
// basic compression
#include <zlib.h>

#include "virt_chr.h"
#include "GTL.h"
#include "gtlVersion.h"

#ifndef __APPLE__
#include "fastalign.h"
#endif

//=================================================================================================

#define kIsPairedEnd 0xFFFF /* for ADNIview */
#define kMaxTrimPairs 600000000
// this loads a binary blob output containing the trim data
#define BIN_OUTPUT

//=================================================================================================

//static struct timeval ts;
//static struct timeval te;
static int verbose = 0;
static int debug = 0;
static VIRTUALCHR virtchr[kMAXCHRcnt];
static unsigned int statMMpos[256];
static unsigned int maxMMpos;
static int gFilterSAM;

typedef struct OPTIONS_struct {
	int fromchr;
	int tochr;
	int fromPos;
	int toPos;
	int decompressPerfectMatches;
	int decompressMismatches;
	int decompressMismatchesN;
	int decompressAligned;
	int decompressUnmapped;
	int decompressHalfmapped;
	int decompressChimeric;
	int filterSAM;
} OPTIONS;

typedef struct _decodedSingle {
	unsigned long ordinal;
	char tag1[kMaxReadLen * 2];
	char hdr1[kHdrSize];
	char qual1[kMaxReadLen];
	unsigned int genomepos[kMAXnbMatches];
	unsigned int mmposStart[kMAXnbMatches];
	unsigned int cigarposStart[kMAXnbMatches];
	int taglen1;
	int hlen1;
	int hhlen1; // the non-blank first part of the header, used in SAM and other formats
	int reverseTAG1[kMAXnbMatches];
	unsigned int chr[kMAXnbMatches];
	int flag1[kMAXnbMatches];
	int bamStyle;
	int nbMismatch1;
} decodedSingle_t, *decodedSingle_p_t;

typedef struct _decodedPair {
	unsigned long ordinal;
	//char *hdrs;
	//char *qual;
	char tag1[kMaxReadLen];
	char tag2[kMaxReadLen];
	char hdr1[kHdrSize];
	char hdr2[kHdrSize];
	char qual1[kMaxReadLen];
	char qual2[kMaxReadLen];
	unsigned int genomepos[kMAXnbMatches];
	unsigned int mmposStart[kMAXnbMatches];
	unsigned int cigarposStart[kMAXnbMatches];
	int taglen1;
	int taglen2;
	int hlen1;
	int hlen2;
	int hhlen1; // the non-blank first part of the header, used in SAM and other formats
	int hhlen2; // the non-blank first part of the header, used in SAM and other formats
	int reverseTAG1[kMAXnbMatches];
	int reverseTAG2[kMAXnbMatches];
	int delta[kMAXnbMatches];
	unsigned int chr[kMAXnbMatches];
	int flag1[kMAXnbMatches];
	int flag2[kMAXnbMatches];
	int ilen[kMAXnbMatches];
	int bamStyle;
	int nbMismatch1;
	int nbMismatch2;
} decodedPair_t, *decodedPair_p_t;

static void print_single_fastq(TLBDATA *, unsigned int, decodedSingle_p_t);
static void print_pair_fastq(TLBDATA *, unsigned int, decodedPair_p_t);
static void print_single_fetchGWI_p(TLBDATA *, unsigned int, decodedSingle_p_t);
static void print_single_fetchGWI_m(TLBDATA *, unsigned int, decodedSingle_p_t);
static void print_single_fetchGWI_mN(TLBDATA *, unsigned int, decodedSingle_p_t);
static void print_single_fetchGWI_g(TLBDATA *, unsigned int, decodedSingle_p_t);
static void print_single_fetchGWI_u(TLBDATA *, unsigned int, decodedSingle_p_t);
static void print_pair_fetchGWI_p(TLBDATA *, unsigned int, decodedPair_p_t);
static void print_pair_fetchGWI_m(TLBDATA *, unsigned int, decodedPair_p_t);
static void print_pair_fetchGWI_mN(TLBDATA *, unsigned int, decodedPair_p_t);
static void print_pair_fetchGWI_g(TLBDATA *, unsigned int, decodedPair_p_t);
static void print_pair_fetchGWI_u(TLBDATA *, unsigned int, decodedPair_p_t);
static void print_single_SAM_p(TLBDATA *, unsigned int, decodedSingle_p_t);
static void print_single_SAM_m(TLBDATA *, unsigned int, decodedSingle_p_t);
static void print_single_SAM_mN(TLBDATA *, unsigned int, decodedSingle_p_t);
static void print_single_SAM_g(TLBDATA *, unsigned int, decodedSingle_p_t);
static void print_single_SAM_u(TLBDATA *, unsigned int, decodedSingle_p_t);
static void print_single_contig_pmNg(TLBDATA *, unsigned int, decodedSingle_p_t);
static void print_single_contig_u(TLBDATA *, unsigned int, decodedSingle_p_t);
static void print_pair_SAM_p(TLBDATA *, unsigned int, decodedPair_p_t);
static void print_pair_SAM_m(TLBDATA *, unsigned int, decodedPair_p_t);
static void print_pair_SAM_mN(TLBDATA *, unsigned int, decodedPair_p_t);
static void print_pair_SAM_g(TLBDATA *, unsigned int, decodedPair_p_t);
static void print_pair_SAM_u(TLBDATA *, unsigned int, decodedPair_p_t);
static void print_single_ADNIview_p(TLBDATA *, unsigned int, decodedSingle_p_t);
static void print_single_ADNIview_m(TLBDATA *, unsigned int, decodedSingle_p_t);
static void print_single_ADNIview_mN(TLBDATA *, unsigned int, decodedSingle_p_t);
static void print_single_ADNIview_g(TLBDATA *, unsigned int, decodedSingle_p_t);
static void print_single_ADNIview_u(TLBDATA *, unsigned int, decodedSingle_p_t);
static void print_pair_ADNIview_p(TLBDATA *, unsigned int, decodedPair_p_t);
static void print_pair_ADNIview_m(TLBDATA *, unsigned int, decodedPair_p_t);
static void print_pair_ADNIview_mN(TLBDATA *, unsigned int, decodedPair_p_t);
static void print_pair_ADNIview_g(TLBDATA *, unsigned int, decodedPair_p_t);
static void print_pair_ADNIview_u(TLBDATA *, unsigned int, decodedPair_p_t);
static void print_single_refless_p(TLBDATA *, unsigned int, decodedSingle_p_t);
static void print_single_refless_m(TLBDATA *, unsigned int, decodedSingle_p_t);
static void print_single_refless_mN(TLBDATA *, unsigned int, decodedSingle_p_t);
static void print_single_refless_g(TLBDATA *, unsigned int, decodedSingle_p_t);
static void print_single_refless_u(TLBDATA *, unsigned int, decodedSingle_p_t);
static void print_pair_refless_p(TLBDATA *, unsigned int, decodedPair_p_t);
static void print_pair_refless_m(TLBDATA *, unsigned int, decodedPair_p_t);
static void print_pair_refless_mN(TLBDATA *, unsigned int, decodedPair_p_t);
static void print_pair_refless_g(TLBDATA *, unsigned int, decodedPair_p_t);
static void print_pair_refless_u(TLBDATA *, unsigned int, decodedPair_p_t);
static void print_pair_contig_pmNg(TLBDATA *, unsigned int, decodedPair_p_t);
static void print_pair_contig_u(TLBDATA *, unsigned int, decodedPair_p_t);
static void print_single_tromer_p(TLBDATA *, unsigned int, decodedSingle_p_t);
static void print_single_tromer_m(TLBDATA *, unsigned int, decodedSingle_p_t);
static void print_single_tromer_mN(TLBDATA *, unsigned int, decodedSingle_p_t);
static void print_single_tromer_g(TLBDATA *, unsigned int, decodedSingle_p_t);
static void print_single_tromer_u(TLBDATA *, unsigned int, decodedSingle_p_t);
static void print_pair_tromer_p(TLBDATA *, unsigned int, decodedPair_p_t);
static void print_pair_tromer_m(TLBDATA *, unsigned int, decodedPair_p_t);
static void print_pair_tromer_mN(TLBDATA *, unsigned int, decodedPair_p_t);
static void print_pair_tromer_g(TLBDATA *, unsigned int, decodedPair_p_t);
static void print_pair_tromer_u(TLBDATA *, unsigned int, decodedPair_p_t);
#ifndef __APPLE__
static void print_single_realign(TLBDATA *, unsigned int, decodedSingle_p_t);
static void print_pair_realign(TLBDATA *, unsigned int, decodedPair_p_t);
#endif
static void print_single_raw(TLBDATA *, unsigned int, decodedSingle_p_t);
static void print_pair_raw(TLBDATA *, unsigned int, decodedPair_p_t);

typedef struct _printFunc {
	char *name;
	int bamStyle;
	int lowercase;
	int rangeSelect;
	void (*print_single_p) (TLBDATA *, unsigned int, decodedSingle_p_t);
	void (*print_single_m) (TLBDATA *, unsigned int, decodedSingle_p_t);
	void (*print_single_mN) (TLBDATA *, unsigned int, decodedSingle_p_t);
	void (*print_single_g) (TLBDATA *, unsigned int, decodedSingle_p_t);
	void (*print_single_u) (TLBDATA *, unsigned int, decodedSingle_p_t);
	void (*print_pair_p) (TLBDATA *, unsigned int, decodedPair_p_t);
	void (*print_pair_m) (TLBDATA *, unsigned int, decodedPair_p_t);
	void (*print_pair_mN) (TLBDATA *, unsigned int, decodedPair_p_t);
	void (*print_pair_g) (TLBDATA *, unsigned int, decodedPair_p_t);
	void (*print_pair_u) (TLBDATA *, unsigned int, decodedPair_p_t);
} printFunc_t, *printFunc_p_t;

// Output types
static const printFunc_t PrintFuncs[] = {
	{"fastq",
		0,
		0,
		0,
		print_single_fastq, print_single_fastq, print_single_fastq, print_single_fastq, print_single_fastq,
		print_pair_fastq, print_pair_fastq, print_pair_fastq, print_pair_fastq, print_pair_fastq},
	{"fetchGWI",
		0,
		0,
		0,
		print_single_fetchGWI_p, print_single_fetchGWI_m, print_single_fetchGWI_mN, print_single_fetchGWI_g, print_single_fetchGWI_u,
		print_pair_fetchGWI_p, print_pair_fetchGWI_m, print_pair_fetchGWI_mN, print_pair_fetchGWI_g, print_pair_fetchGWI_u},
	{"SAM",
		1,
		0,
		0,
		print_single_SAM_p, print_single_SAM_m, print_single_SAM_mN, print_single_SAM_g, print_single_SAM_u,
		print_pair_SAM_p, print_pair_SAM_m, print_pair_SAM_mN, print_pair_SAM_g, print_pair_SAM_u},
	{"ADNIview",
		1,
		32,
		1,
		print_single_ADNIview_p, print_single_ADNIview_m, print_single_ADNIview_mN, print_single_ADNIview_g, print_single_ADNIview_u,
		print_pair_ADNIview_p, print_pair_ADNIview_m, print_pair_ADNIview_mN, print_pair_ADNIview_g, print_pair_ADNIview_u},
	{"refless",
		0,
		0,
		0,
		print_single_refless_p, print_single_refless_m, print_single_refless_mN, print_single_refless_g, print_single_refless_u,
		print_pair_refless_p, print_pair_refless_m, print_pair_refless_mN, print_pair_refless_g, print_pair_refless_u},
	{"contig",
		0,
		0,
		1,
		print_single_contig_pmNg, print_single_contig_pmNg, print_single_contig_pmNg, print_single_contig_pmNg, print_single_contig_u,
		print_pair_contig_pmNg, print_pair_contig_pmNg, print_pair_contig_pmNg, print_pair_contig_pmNg, print_pair_contig_u},
	{"tromer",
		0,
		0,
		0,
		print_single_tromer_p, print_single_tromer_m, print_single_tromer_mN, print_single_tromer_g, print_single_tromer_u,
		print_pair_tromer_p, print_pair_tromer_m, print_pair_tromer_mN, print_pair_tromer_g, print_pair_tromer_u},
#ifndef __APPLE__
	{"realign",
		1,
		32,
		1,
		print_single_realign, print_single_realign, print_single_realign, print_single_realign, print_single_realign,
		print_pair_realign, print_pair_realign, print_pair_realign, print_pair_realign, print_pair_ADNIview_u},
#endif
	{"raw",
		1,
		32,
		0,
		print_single_raw, print_single_raw, print_single_raw, print_single_raw, print_single_raw,
		print_pair_raw, print_pair_raw, print_pair_raw, print_pair_raw, print_pair_raw}
};
#define kOUTTYPESNB (sizeof(PrintFuncs) / sizeof(printFunc_t))
static int outSelect;

static int gFromPos;
static int gToPos;

typedef struct _trimData {
	unsigned char cnt[4];
	unsigned char len[4];
} trimData_t, *trimData_p_t;

static trimData_p_t trim;

//=================================================================================================



unsigned char *genome_tbl;	// contain genome in 16 virtual chromosomes.
int loadGenome(char *genomefile)
{
	int err = 1;
	int sf = 0;
	struct stat genomeStat;

		/* ----------------------- load genome data ----------------------*/

		sf = open ((char*)genomefile, O_RDONLY);
		if (!sf)
		{
			printf("can't open bin file\n");
			goto bail;
		}

		if (fstat (sf,&genomeStat) != 0)
		{
			printf("cannot stat file %s\n", (char*)genomefile);
			goto bail;
		}
		if (genomeStat.st_size != kGENOME_DATA_SIZE)
		{
			printf("genomeStat.st_size %tx != %lx\n",genomeStat.st_size,kGENOME_DATA_SIZE);
			goto bail;
		}

		genome_tbl = (unsigned char *)mmap (genome_tbl, kGENOME_DATA_SIZE, PROT_READ, MAP_PRIVATE, sf, 0);
		if (genome_tbl == MAP_FAILED)
		{
			printf("mapping failed.\n");
			goto bail;
		}

		err = 0;

	bail:;

		return(err);

} // loadGenome
//---------------------------------------------------------------

static void
decodeHeaderSingle(TLBDATA *tpp, decodedSingle_p_t dsp)
{
	char *next = strchr(tpp->hdr, '\n');

	dsp->hlen1 = next - tpp->hdr;
	memcpy(dsp->hdr1,tpp->hdr,dsp->hlen1);
	tpp->hdr = next + 1;
	next = strchr(dsp->hdr1, ' ');
	dsp->hhlen1 = dsp->hlen1;
	if (next != NULL)
		dsp->hhlen1 = next - dsp->hdr1;
	for (unsigned int i = 0; i < 8; i++)
		dsp->ordinal = (dsp->ordinal << 8) | (unsigned long) *((unsigned char *) tpp->hdr++);
}
//---------------------------------------------------------------

static void
decodeSingle(TLBDATA *tpp, unsigned int i, decodedSingle_p_t dsp, int bamStyle)
{
	memset(dsp,0,sizeof(decodedSingle_t));
	dsp->bamStyle = bamStyle;

	if (tpp->header.flags_7_0 & kTLBHflagFixedLength)
	{
		memcpy(dsp->qual1,tpp->qs,tpp->lengthInfo);
		tpp->qs += tpp->lengthInfo;
		dsp->taglen1 = tpp->lengthInfo;
	}
	else
	{
		memcpy(dsp->qual1,tpp->qs,tpp->len[i]);
		tpp->qs += tpp->len[i];
		dsp->taglen1 = tpp->len[i];
	}
	decodeHeaderSingle(tpp,dsp);

	for (unsigned int j = 0; j < tpp->nbMatches; j++)
	{
		if (tpp->ppd[j][i].tag2offset & k32bREVCOMP_TAG1_BIT)
			dsp->reverseTAG1[j] = 1;
		dsp->chr[j] = (tpp->ppd[j][i].tag2offset & k32b2nd_MATCH_CHR_MASK) >> 19;
		if (dsp->chr[j] == 0)
			dsp->chr[j] = tpp->header.chr;
		if (dsp->reverseTAG1[j]) // reverse
			dsp->flag1[j] |= 0x10;
	}

	dsp->genomepos[0] = (virtchr[tpp->header.chr].chr << 28) + virtchr[tpp->header.chr].offset + tpp->ppd[0][i].tag1pos;
	memcpy(dsp->tag1,&genome_tbl[dsp->genomepos[0]],dsp->taglen1);
	if (dsp->reverseTAG1[0]) // reverse
	{
		int r;
		int a = 0;
		if (bamStyle) {
			char qualr[kMaxReadLen];
			memcpy(qualr,dsp->qual1,dsp->taglen1);
			for (r = dsp->taglen1-1; r>=0; r--)
			{
				dsp->qual1[a++] = qualr[r];
			}
		}
		else
		{
			for (r = dsp->taglen1-1; r>=0; r--)
			{
				dsp->tag1[a++] = gIUPACrevcomp[ genome_tbl[dsp->genomepos[0]+r] ];
			}
		}
	}
}
//---------------------------------------------------------------

static void
decodeHeaderPair(TLBDATA *tpp, decodedPair_p_t dpp)
{
	char *next = strchr(tpp->hdr, '\t');

	dpp->hlen1 = next - tpp->hdr;
	memcpy(dpp->hdr1,tpp->hdr,dpp->hlen1);
	tpp->hdr = next + 1;
	next = strchr(tpp->hdr, '\n');
	dpp->hlen2 = next - tpp->hdr;
	memcpy(dpp->hdr2,tpp->hdr,dpp->hlen2);
	tpp->hdr = next + 1;
	// FIXME !!  we can do faster if necessary...
	if (dpp->hdr2[0] == 0)
	{
		if (dpp->bamStyle) // remove trailing dot or slash - maybe should check if really . or / ... ?
		{
			dpp->hdr1[--(dpp->hlen1)] = 0;
		}
		strcpy(dpp->hdr2,dpp->hdr1);
		if (!dpp->bamStyle)
		{
			strcat(dpp->hdr2,"2");
			strcat(dpp->hdr1,"1");
			dpp->hlen1 += 1;
		}
		dpp->hlen2 = dpp->hlen1;
	}
	next = strchr(dpp->hdr1, ' ');
	dpp->hhlen1 = dpp->hlen1;
	if (next != NULL)
		dpp->hhlen1 = next - dpp->hdr1;
	next = strchr(dpp->hdr2, ' ');
	dpp->hhlen2 = dpp->hlen2;
	if (next != NULL)
		dpp->hhlen2 = next - dpp->hdr2;
	for (unsigned int i = 0; i < 8; i++)
		dpp->ordinal = (dpp->ordinal << 8) | (unsigned long) *((unsigned char *) tpp->hdr++);
}
//---------------------------------------------------------------

static void
decodePair(TLBDATA *tpp, unsigned int i, decodedPair_p_t dpp, int bamStyle)
{
	memset(dpp,0,sizeof(decodedPair_t));
	dpp->bamStyle = bamStyle;

	if (tpp->header.flags_7_0 & kTLBHflagFixedLength)
	{
		memcpy(dpp->qual1,tpp->qs,tpp->lengthInfo);
		tpp->qs += tpp->lengthInfo;
		memcpy(dpp->qual2,tpp->qs,tpp->lengthInfo);
		tpp->qs += tpp->lengthInfo;
		dpp->taglen1 = tpp->lengthInfo;
		dpp->taglen2 = tpp->lengthInfo;
	}
	else
	{
		memcpy(dpp->qual1,tpp->qs,tpp->len[2*i]);
		tpp->qs += tpp->len[2*i];
		memcpy(dpp->qual2,tpp->qs,tpp->len[2*i+1]);
		tpp->qs += tpp->len[2*i+1];
		dpp->taglen1 = tpp->len[2*i];
		dpp->taglen2 = tpp->len[2*i+1];
	}
	decodeHeaderPair(tpp,dpp);

	for (unsigned int j = 0; j < tpp->nbMatches; j++)
	{
		dpp->flag1[j] = (j > 0) ? 0x143 : 0x43;
		dpp->flag2[j] = (j > 0) ? 0x183 : 0x83;
		if (tpp->ppd[j][i].tag2offset & k32bREVCOMP_TAG1_BIT)
			dpp->reverseTAG1[j] = 1;

		if (tpp->ppd[j][i].tag2offset & k32bREVCOMP_TAG2_BIT)
			dpp->reverseTAG2[j] = 1;
		dpp->delta[j] = (tpp->ppd[j][i].tag2offset & k32bOFFSET_VALUE_MASK);
		if (tpp->ppd[j][i].tag2offset & k32bNEGATIVE_OFFSET_BIT)
			dpp->delta[j] = -dpp->delta[j];
		dpp->chr[j] = (tpp->ppd[j][i].tag2offset & k32b2nd_MATCH_CHR_MASK) >> 19;
		if (dpp->chr[j] == 0)
			dpp->chr[j] = tpp->header.chr;

		if (dpp->delta[j] >= 0) {
			dpp->ilen[j] = dpp->delta[j] + dpp->taglen2 - 1;
		} else {
			dpp->ilen[j] = dpp->delta[j] - dpp->taglen1 + 1;
		}
		if (dpp->reverseTAG1[j]) // reverse
		{
			dpp->flag1[j] |= 0x10;
			dpp->flag2[j] |= 0x20;
		}
		if (dpp->reverseTAG2[j]) // reverse
		{
			dpp->flag1[j] |= 0x20;
			dpp->flag2[j] |= 0x10;
		}
	}

	dpp->genomepos[0] = (virtchr[tpp->header.chr].chr << 28) + virtchr[tpp->header.chr].offset + tpp->ppd[0][i].tag1pos;
	memcpy(dpp->tag1,&genome_tbl[dpp->genomepos[0]],dpp->taglen1);
	if (dpp->reverseTAG1[0]) // reverse
	{
		int r;
		int a = 0;
		if (bamStyle) {
			char qualr[kMaxReadLen];
			memcpy(qualr,dpp->qual1,dpp->taglen1);
			for (r = dpp->taglen1-1; r>=0; r--)
			{
				dpp->qual1[a++] = qualr[r];
			}
		}
		else
		{
			for (r = dpp->taglen1-1; r>=0; r--)
			{
				dpp->tag1[a++] = gIUPACrevcomp[ genome_tbl[dpp->genomepos[0]+r] ];
			}
		}
	}

	memcpy(dpp->tag2,&genome_tbl[dpp->genomepos[0] + dpp->delta[0]],dpp->taglen2);
	if (dpp->reverseTAG2[0]) // reverse
	{
		int r;
		int a = 0;
		if (bamStyle) {
			char qualr[kMaxReadLen];
			memcpy(qualr,dpp->qual2,dpp->taglen2);
			for (r = dpp->taglen2-1; r>=0; r--)
			{
				dpp->qual2[a++] = qualr[r];
			}
		}
		else
		{
			for (r = dpp->taglen2-1; r>=0; r--)
			{
				dpp->tag2[a++] = gIUPACrevcomp[ genome_tbl[dpp->genomepos[0]+dpp->delta[0]+r] ];
			}
		}
	}
}
//---------------------------------------------------------------

static int
applyNdiffSingle(decodedSingle_p_t dsp, unsigned char *diff, int mmpos)
{
	dsp->mmposStart[0] = mmpos;
	while (diff[mmpos] != kMMTERMINATOR)
	{
		int pos = diff[mmpos++];
		if (dsp->reverseTAG1[0] && !dsp->bamStyle) // reverse tag 1
			dsp->tag1[dsp->taglen1-1-pos] = 'N';
		else
			dsp->tag1[pos] = 'N';
	}
	return mmpos+1; // skip the kMMTERMINATOR
}
//---------------------------------------------------------------

static int
applyNdiffPair(decodedPair_p_t dpp, unsigned char *diff, int mmpos)
{
	dpp->mmposStart[0] = mmpos;
	// check if 1st tag has mismatches (warning no need to revcomp differences)
	while (diff[mmpos] != kMMTERMINATOR)
	{
		int pos = diff[mmpos++];
		if (dpp->reverseTAG1[0] && !dpp->bamStyle) // reverse tag 1
			dpp->tag1[dpp->taglen1-1-pos] = 'N';
		else
			dpp->tag1[pos] = 'N';
	}
	// check if 2nd tag has mismatches (warning no need to revcomp differences)
	mmpos++; // skip the kMMTERMINATOR
	while (diff[mmpos] != kMMTERMINATOR)
	{
		int pos = diff[mmpos++];
		if (dpp->reverseTAG2[0] && !dpp->bamStyle) // reverse tag 2
			dpp->tag2[dpp->taglen2-1-pos] = 'N';
		else
			dpp->tag2[pos] = 'N';
	}
	return mmpos+1; // skip the kMMTERMINATOR
}
//---------------------------------------------------------------

static int
applyDiffSingle(decodedSingle_p_t dsp, unsigned char *diff, int mmpos, int lowercase)
{
	dsp->mmposStart[0] = mmpos;
	while (diff[mmpos] != kMMTERMINATOR)
	{
		if (diff[mmpos] == kMMSOFTCLIP)
		{
			// grab the nucleotides clipped in a tmp string
			unsigned char tmpStr[kMaxReadLen];
			int i = 0;
			mmpos++;
			while (diff[mmpos] != kMMTERMINATOR) { tmpStr[i++] = diff[mmpos++] + lowercase; }
			// copy the clipped nucleotides at the end of the tag.
			if (dsp->reverseTAG1[0] && dsp->bamStyle) // reverse tag 1
			{
				unsigned int j;
				for (j = 0; j < i; j++)
					dsp->tag1[j] = gIUPACrevcomp[ tmpStr[i - j - 1] ];
			}
			else
				memcpy(&dsp->tag1[dsp->taglen1-i],tmpStr,i);
			break; // kMMSOFTCLIP is the last thing in MMstring. only one occurence possible. we are done.
		}
		int pos = diff[mmpos++];
		unsigned char nt = diff[mmpos++];
		if (nt != 'N')
		{
			dsp->nbMismatch1 += 1;
			if (pos > maxMMpos)
				maxMMpos = pos;
			if (dsp->reverseTAG1[0])
				statMMpos[dsp->taglen1-1-pos] += 1;
			else
				statMMpos[pos] += 1;
		}
		if (dsp->reverseTAG1[0] && !dsp->bamStyle) // reverse tag 1
			dsp->tag1[dsp->taglen1-1-pos] = gIUPACrevcomp[ nt ];
		else
			dsp->tag1[pos] = nt;
	}
	return mmpos+1; // skip the kMMTERMINATOR
}
//---------------------------------------------------------------

static int
applyDiffPair(decodedPair_p_t dpp, unsigned char *diff, int mmpos, int lowercase)
{
	dpp->mmposStart[0] = mmpos;
	// check if 1st tag has mismatches (warning no need to revcomp differences)
	while (diff[mmpos] != kMMTERMINATOR)
	{
		if (diff[mmpos] == kMMSOFTCLIP)
		{
			// grab the nucleotides clipped in a tmp string
			unsigned char tmpStr[kMaxReadLen];
			int i = 0;
			mmpos++;
			while (diff[mmpos] != kMMTERMINATOR) { tmpStr[i++] = diff[mmpos++] + lowercase; }
			// copy the clipped nucleotides at the end of the tag.
			if (dpp->reverseTAG1[0] && dpp->bamStyle) // reverse tag 1
			{
				unsigned int j;
				for (j = 0; j < i; j++)
					dpp->tag1[j] = gIUPACrevcomp[ tmpStr[i - j - 1] ];
			}
			else
				memcpy(&dpp->tag1[dpp->taglen1-i],tmpStr,i);
			break; // kMMSOFTCLIP is the last thing in MMstring. only one occurence possible. we are done.
		}
		int pos = diff[mmpos++];
		unsigned char nt = diff[mmpos++];
		if (nt != 'N')
		{
			dpp->nbMismatch1 += 1;
			if (pos > maxMMpos)
				maxMMpos = pos;
			if (dpp->reverseTAG1[0])
				statMMpos[dpp->taglen1-1-pos] += 1;
			else
				statMMpos[pos] += 1;
		}
		if (dpp->reverseTAG1[0] && !dpp->bamStyle) // reverse tag 1
			dpp->tag1[dpp->taglen1-1-pos] = gIUPACrevcomp[ nt ];
		else
			dpp->tag1[pos] = nt;
	}
	// check if 2nd tag has mismatches (warning no need to revcomp differences)
	mmpos++; // skip the kMMTERMINATOR
	//printf("entering tag2MM with mmpos=%d\n",mmpos);
	while (diff[mmpos] != kMMTERMINATOR)
	{
		//printf("diff[%d]=%d\n",mmpos,diff[mmpos]);
		if (diff[mmpos] == kMMSOFTCLIP)
		{
			// grab the nucleotides clipped in a tmp string
			unsigned char tmpStr[kMaxReadLen];
			int i = 0;
			mmpos++;
			while (diff[mmpos] != kMMTERMINATOR) { tmpStr[i++] = diff[mmpos++] + lowercase; }
			// copy the clipped nucleotides at the end of the tag.
			//printf("i=%d replacepos=%d\n",i,dpp->taglen2-1-i);
			if (dpp->reverseTAG2[0] && dpp->bamStyle) // reverse tag 2
			{
				unsigned int j;
				for (j = 0; j < i; j++)
					dpp->tag2[j] = gIUPACrevcomp[ tmpStr[i - j - 1] ];
			}
			else
				memcpy(&dpp->tag2[dpp->taglen2-i],tmpStr,i);
			break; // kMMSOFTCLIP is the last thing in MMstring. only one occurence possible. we are done.
		}
		int pos = diff[mmpos++];
		unsigned char nt = diff[mmpos++];
		if (nt != 'N')
		{
			dpp->nbMismatch2 += 1;
			if (pos > maxMMpos)
				maxMMpos = pos;
			if (dpp->reverseTAG2[0])
				statMMpos[dpp->taglen2-1-pos] += 1;
			else
				statMMpos[pos] += 1;
		}
		if (dpp->reverseTAG2[0] && !dpp->bamStyle) // reverse tag 2
			dpp->tag2[dpp->taglen2-1-pos] = gIUPACrevcomp[ nt ];
		else
			dpp->tag2[pos] = nt;
	}
	return mmpos+1; // skip the kMMTERMINATOR
}
//---------------------------------------------------------------

static void
AlignTagsOnGenomeUsingCigarSingle(char *tag1, char *qual1, decodedSingle_p_t dsp, unsigned char *cigar, unsigned int *off1)
{
	unsigned int cigarpos = dsp->cigarposStart[0];
	// check if 1st tag has a gapped alignment cigar
	if (cigar[cigarpos] != kCIGARTERMINATOR)	// if not, normal decoding.
	{
		int r = 0;
		int a = 0;

		while (cigar[cigarpos] != kCIGARTERMINATOR)
		{
			int pos = cigar[cigarpos++];
			char code = cigar[cigarpos++];
			switch(code)
			{
				case 'S': // same as 'M'
				case 'M':  while(pos-- >0) { tag1[a]   = dsp->tag1[r]; qual1[a]   = dsp->qual1[r]; a++; r++;                     }             break;
				case 'D':  while(pos-- >0) { tag1[a]   = '_';          qual1[a]   = '_';           a++;       dsp->taglen1++;    }             break;
				case 'I':
					if (a > 0)
					{
						tag1[a-1] = '*';
						qual1[a-1] = '*';
						r += pos;
						dsp->taglen1 -= pos;
					}
					else
					{
						tag1[a] = '*';
						qual1[a] = '*';
						r += pos;
						a += 1;
						dsp->taglen1 -= (pos - 1);
						*off1 = 1;
					}
					break;
			}
			assert(a < kMaxReadLen * 2);
		}
		tag1[a] = '\0';
		qual1[a] = '\0';
	}
	else
	{
		strcpy(tag1,dsp->tag1);
		strcpy(qual1,dsp->qual1);
	}
}
//---------------------------------------------------------------

static void
AlignTagsOnGenomeUsingCigar(char *tag1, char *tag2, char *qual1, char *qual2, decodedPair_p_t dpp, unsigned char *cigar, unsigned int *off1, unsigned int *off2)
{
	unsigned int cigarpos = dpp->cigarposStart[0];
	// check if 1st tag has a gapped alignment cigar
	if (cigar[cigarpos] != kCIGARTERMINATOR)	// if not, normal decoding.
	{
		int r = 0;
		int a = 0;

		while (cigar[cigarpos] != kCIGARTERMINATOR)
		{
			int pos = cigar[cigarpos++];
			char code = cigar[cigarpos++];
			switch(code)
			{
				case 'S': // same as 'M'
				case 'M':  while(pos-- >0) { tag1[a]   = dpp->tag1[r]; qual1[a]   = dpp->qual1[r]; a++; r++;                     }             break;
				case 'D':  while(pos-- >0) { tag1[a]   = '_';          qual1[a]   = '_';           a++;       dpp->taglen1++;    }             break;
				case 'I':
					if (a > 0)
					{
						tag1[a-1] = '*';
						qual1[a-1] = '*';
						r += pos;
						dpp->taglen1 -= pos;
					}
					else
					{
						tag1[a] = '*';
						qual1[a] = '*';
						r += pos;
						a += 1;
						dpp->taglen1 -= (pos - 1);
						*off1 = 1;
					}
					break;
			}
			assert(a < kMaxReadLen * 2);
		}
		tag1[a] = '\0';
		qual1[a] = '\0';
	}
	else
	{
		strcpy(tag1,dpp->tag1);
		strcpy(qual1,dpp->qual1);
	}
	cigarpos++; // skip the kCIGARTERMINATOR

	// check if 2nd tag has a gapped alignment cigar
	if (cigar[cigarpos] != kCIGARTERMINATOR)	// if not, normal decoding.
	{
		int r = 0;
		int a = 0;

		while (cigar[cigarpos] != kCIGARTERMINATOR)
		{
			int pos = cigar[cigarpos++];
			char code = cigar[cigarpos++];
			switch(code)
			{
				case 'S': // same as 'M'
				case 'M':  while(pos-- >0) { tag2[a]   = dpp->tag2[r]; qual2[a]   = dpp->qual2[r]; a++; r++;                     }             break;
				case 'D':  while(pos-- >0) { tag2[a]   = '_';          qual2[a]   = '_';           a++;       dpp->taglen2++;    }             break;
				case 'I':
					if (a > 0)
					{
						tag2[a-1] = '*';
						qual2[a-1] = '*';
						r += pos;
						dpp->taglen2 -= pos;
					}
					else
					{
						tag2[a] = '*';
						qual2[a] = '*';
						r += pos;
						a += 1;
						dpp->taglen2 -= (pos - 1);
						*off2 = 1;
					}
					break;
			}
			assert(a < kMaxReadLen * 2);
		}
		tag2[a] = '\0';
		qual2[a] = '\0';
	}
	else
	{
		strcpy(tag2,dpp->tag2);
		strcpy(qual2,dpp->qual2);
	}

} // AlignTagsOnGenomeUsingCigar
//---------------------------------------------------------------

static int
applyCigarSingle(decodedSingle_p_t dsp, unsigned char *cigar, int cigarpos)
{
	dsp->cigarposStart[0] = cigarpos;
	// check if 1st tag has a gapped alignment cigar
	if (cigar[cigarpos] != kCIGARTERMINATOR)	// if not, normal decoding.
	{
		int r = 0;
		int a = 0;
		char tmp[kMaxReadLen];

		while (cigar[cigarpos] != kCIGARTERMINATOR)
		{
			int pos = cigar[cigarpos++];
			char code = cigar[cigarpos++];
			switch(code)
			{
				case 'S': // same as 'M'
				case 'M': assert(a+pos<=kMaxReadLen); while(pos-- >0) { tmp[a++] = genome_tbl[dsp->genomepos[0]+r]; r++; } break;
				case 'D':                          while(pos-- >0) { r++;                                          } break;
				case 'I': assert(a+pos<=kMaxReadLen); while(pos-- >0) { tmp[a++] = 'N';                               } break;
			}
		}
		assert(dsp->taglen1<kMaxReadLen);
		tmp[dsp->taglen1] = '\0';
		if (dsp->reverseTAG1[0] && !dsp->bamStyle) // reverse tag 1
		{
			int r;
			int a = 0;
			for (r = dsp->taglen1-1; r>=0; r--)
				dsp->tag1[a++] = gIUPACrevcomp[ (unsigned char) tmp[r] ];
			dsp->tag1[a] = '\0';
		}
		else
			strcpy(dsp->tag1,tmp);
	}
	return cigarpos+1; // skip the kCIGARTERMINATOR
}
//---------------------------------------------------------------

static int
applyCigarPair(decodedPair_p_t dpp, unsigned char *cigar, int cigarpos)
{
	char tmp[kMaxReadLen];
	int a, r;
	dpp->cigarposStart[0] = cigarpos;
	// check if 1st tag has a gapped alignment cigar
	if (cigar[cigarpos] != kCIGARTERMINATOR)	// if not, normal decoding.
	{
		r = 0;
		a = 0;

		while (cigar[cigarpos] != kCIGARTERMINATOR)
		{
			int pos = cigar[cigarpos++];
			char code = cigar[cigarpos++];
			if (a >= dpp->taglen1 || pos > kMaxReadLen)
			{
				// FIXME - this should not have been encoded this way...
				fprintf(stderr,"FIXME : Ordinal %lu Read 1 r:%d a:%d pos:%d code:%c taglen:%d\n",dpp->ordinal,r,a,pos,code,dpp->taglen2);
			}
			else
			{
				switch(code)
				{
					case 'S': // same as 'M'
					case 'M': assert(a+pos<=kMaxReadLen); while(pos-- >0) { tmp[a++] = genome_tbl[dpp->genomepos[0]+r]; r++; } break;
					case 'D':                             while(pos-- >0) { r++;                                          } break;
					case 'I': assert(a+pos<=kMaxReadLen); while(pos-- >0) { tmp[a++] = 'N';                               } break;
				}
			}
		}
		assert(dpp->taglen1<kMaxReadLen);
		tmp[dpp->taglen1] = '\0';
		if (dpp->reverseTAG1[0] && !dpp->bamStyle) // reverse tag 1
		{
			a = 0;
			for (r = dpp->taglen1-1; r>=0; r--)
				dpp->tag1[a++] = gIUPACrevcomp[ (unsigned char) tmp[r] ];
			dpp->tag1[a] = '\0';
		}
		else
			strcpy(dpp->tag1,tmp);
	}
	cigarpos++; // skip the kCIGARTERMINATOR
	// check if 2nd tag has a gapped alignment cigar
	if (cigar[cigarpos] != kCIGARTERMINATOR)	// if not, normal decoding.
	{
		r = 0;
		a = 0;
		while (cigar[cigarpos] != kCIGARTERMINATOR)
		{
			int pos = cigar[cigarpos++];
			char code = cigar[cigarpos++];
			if (a >= dpp->taglen2 || pos > kMaxReadLen)
			{
				// FIXME - this should not have been encoded this way...
				fprintf(stderr,"FIXME : Ordinal %lu Read 2 r:%d a:%d pos:%d code:%c taglen:%d\n",dpp->ordinal, r,a,pos,code,dpp->taglen2);
			}
			else
			{
				switch(code)
				{
					case 'S': // same as 'M'
					case 'M': assert(a+pos<=kMaxReadLen); while(pos-- >0) { tmp[a++] = genome_tbl[dpp->genomepos[0]+dpp->delta[0]+r]; r++; } break;
					case 'D':                          while(pos-- >0) { r++;                                                     } break;
					case 'I': assert(a+pos<=kMaxReadLen); while(pos-- >0) { tmp[a++] = 'N';                                          } break;
				}
			}
		}
		assert(dpp->taglen2<kMaxReadLen);
		tmp[dpp->taglen2] = '\0';
		if (dpp->reverseTAG2[0] && !dpp->bamStyle) // reverse tag 2
		{
			a = 0;
			for (r = dpp->taglen2-1; r>=0; r--)
				dpp->tag2[a++] = gIUPACrevcomp[ (unsigned char) tmp[r] ];
			dpp->tag2[a] = '\0';
		}
		else
			strcpy(dpp->tag2,tmp);
	}
	return cigarpos+1; // skip the kCIGARTERMINATOR
}
//---------------------------------------------------------------

const static char NT2bits[4] = { 'A', 'C', 'G', 'T' };

static int
decodeSingleUnmapped(TLBDATA *tpp, unsigned int i, decodedSingle_p_t dsp, int mmpos, int bamStyle)
{
	memset(dsp,0,sizeof(decodedSingle_t));
	dsp->bamStyle = bamStyle;
	dsp->flag1[0] = 4;

	if (tpp->header.flags_7_0 & kTLBHflagFixedLength)
	{
		memcpy(dsp->qual1,tpp->qs,tpp->lengthInfo);
		tpp->qs += tpp->lengthInfo;
		dsp->taglen1 = tpp->lengthInfo;
	}
	else
	{
		memcpy(dsp->qual1,tpp->qs,tpp->len[i]);
		tpp->qs += tpp->len[i];
		dsp->taglen1 = tpp->len[i];
	}
	decodeHeaderSingle(tpp,dsp);
	char *tp = &dsp->tag1[0];
	unsigned int pos =0;
	unsigned int shift = 8;
	while(pos < dsp->taglen1)
	{
		if (shift == 0)
		{
			shift = 8;
			tpp->packed4nt += 1;
		}
		shift -= 2;
		*tp = NT2bits[(*tpp->packed4nt) >> shift  & 0x3];
		tp++;
		pos++;
	}
	tpp->packed4nt += 1;

	// check if 1st tag has mismatches (warning no need to revcomp differences)
	while (tpp->diff[0][mmpos] != kMMTERMINATOR)
	{
		int pos = tpp->diff[0][mmpos++];
		char nt = tpp->diff[0][mmpos++];
		dsp->tag1[pos] = nt;
	}

	return mmpos+1; // skip the kMMTERMINATOR
}
//---------------------------------------------------------------

static int
decodePairUnmapped(TLBDATA *tpp, unsigned int i, decodedPair_p_t dpp, int mmpos, int bamStyle)
{
	memset(dpp,0,sizeof(decodedPair_t));
	dpp->flag1[0] = 0x41;
	dpp->flag2[0] = 0x81;
	dpp->bamStyle = bamStyle;

	if (tpp->header.flags_7_0 & kTLBHflagFixedLength)
	{
		memcpy(dpp->qual1,tpp->qs,tpp->lengthInfo);
		tpp->qs += tpp->lengthInfo;
		memcpy(dpp->qual2,tpp->qs,tpp->lengthInfo);
		tpp->qs += tpp->lengthInfo;
		dpp->taglen1 = tpp->lengthInfo;
		dpp->taglen2 = tpp->lengthInfo;
	}
	else
	{
		memcpy(dpp->qual1,tpp->qs,tpp->len[2*i]);
		tpp->qs += tpp->len[2*i];
		memcpy(dpp->qual2,tpp->qs,tpp->len[2*i+1]);
		tpp->qs += tpp->len[2*i+1];
		dpp->taglen1 = tpp->len[2*i];
		dpp->taglen2 = tpp->len[2*i+1];
	}
	decodeHeaderPair(tpp,dpp);

	if ((tpp->header.flags_7_0 & kTLBHflagPairChrPosBlock) != 0)
	{
		if (IS_TAG1_REVERSED(tpp->pcpd[0][i].chr))
			dpp->reverseTAG1[0] = 1;

		if (IS_TAG2_REVERSED(tpp->pcpd[0][i].chr))
			dpp->reverseTAG2[0] = 1;

		unsigned int chr1 = GET_TAG1_CHR(tpp->pcpd[0][i].chr);
		unsigned int chr2 = GET_TAG2_CHR(tpp->pcpd[0][i].chr);

		if (chr1 != 0)
		{
			dpp->genomepos[0] = (virtchr[chr1].chr << 28) + virtchr[chr1].offset + tpp->pcpd[0][i].tag1pos;
			memcpy(dpp->tag1,&genome_tbl[dpp->genomepos[0]],dpp->taglen1);
			if (dpp->reverseTAG1[0]) // reverse
			{
				int r;
				int a = 0;
				if (bamStyle) {
					char qualr[kMaxReadLen];
					memcpy(qualr,dpp->qual1,dpp->taglen1);
					for (r = dpp->taglen1-1; r>=0; r--)
					{
						dpp->qual1[a++] = qualr[r];
					}
				}
				else
				{
					for (r = dpp->taglen1-1; r>=0; r--)
					{
						dpp->tag1[a++] = gIUPACrevcomp[ genome_tbl[dpp->genomepos[0]+r] ];
					}
				}
				dpp->flag1[0] |= 0x10;
				dpp->flag2[0] |= 0x20;
			}
		}
		else
		{
			char *tp = dpp->tag1;
			unsigned int pos =0;
			unsigned int shift = 8;
			dpp->flag1[0] |= 4;
			dpp->flag2[0] |= 8;
			while(pos < dpp->taglen1)
			{
				if (shift == 0)
				{
					shift = 8;
					tpp->packed4nt += 1;
				}
				shift -= 2;
				*tp = NT2bits[(*tpp->packed4nt) >> shift  & 0x3];
				tp++;
				pos++;
			}
			tpp->packed4nt += 1;
		}

		if (chr2 != 0)
		{
			dpp->genomepos[0] = (virtchr[chr2].chr << 28) + virtchr[chr2].offset + tpp->pcpd[0][i].tag2pos;
			memcpy(dpp->tag2,&genome_tbl[dpp->genomepos[0]],dpp->taglen2);
			if (dpp->reverseTAG2[0]) // reverse
			{
				int r;
				int a = 0;
				if (bamStyle) {
					char qualr[kMaxReadLen];
					memcpy(qualr,dpp->qual2,dpp->taglen2);
					for (r = dpp->taglen2-1; r>=0; r--)
					{
						dpp->qual2[a++] = qualr[r];
					}
				}
				else
				{
					for (r = dpp->taglen2-1; r>=0; r--)
					{
						dpp->tag2[a++] = gIUPACrevcomp[ genome_tbl[dpp->genomepos[0]+r] ];
					}
				}
				dpp->flag1[0] |= 0x20;
				dpp->flag2[0] |= 0x10;
			}
		}
		else
		{
			char *tp = dpp->tag2;
			unsigned int pos =0;
			unsigned int shift = 8;
			dpp->flag1[0] |= 8;
			dpp->flag2[0] |= 4;
			pos =0;
			while(pos < dpp->taglen2)
			{
				if (shift == 0)
				{
					shift = 8;
					tpp->packed4nt += 1;
				}
				shift -= 2;
				*tp = NT2bits[(*tpp->packed4nt) >> shift  & 0x3];
				tp++;
				pos++;
			}
			tpp->packed4nt += 1;
		}
	}
	else
	{
		char *tp = dpp->tag1;
		unsigned int pos =0;
		unsigned int shift = 8;
		dpp->flag1[0] |= 0xc;
		dpp->flag2[0] |= 0xc;

		while(pos < dpp->taglen1)
		{
			if (shift == 0)
			{
				shift = 8;
				tpp->packed4nt += 1;
			}
			shift -= 2;
			*tp = NT2bits[(*tpp->packed4nt) >> shift  & 0x3];
			tp++;
			pos++;
		}
		tpp->packed4nt += 1;

		tp = dpp->tag2;
		pos =0;
		shift = 8;
		while(pos < dpp->taglen2)
		{
			if (shift == 0)
			{
				shift = 8;
				tpp->packed4nt += 1;
			}
			shift -= 2;
			*tp = NT2bits[(*tpp->packed4nt) >> shift  & 0x3];
			tp++;
			pos++;
		}
		tpp->packed4nt += 1;
	}

	mmpos = applyDiffPair(dpp,tpp->diff[0],mmpos,0);

	return mmpos; // skip the kMMTERMINATOR
}
//---------------------------------------------------------------

static int
validTrimSingle(trimData_p_t td)
{
	if (td->cnt[0] == 0 && td->len[0] != 0)
		return 0;
	if (td->cnt[1] == 0 && td->len[1] != 0)
		return 0;
	if (td->cnt[0] > 20 && td->cnt[1] > 20)
		return 0;
	return 1;
}
//---------------------------------------------------------------

static int
validTrim(trimData_p_t td)
{
	if (td->cnt[0] == 0 && td->len[0] != 0)
		return 0;
	if (td->cnt[1] == 0 && td->len[1] != 0)
		return 0;
	if (td->cnt[2] == 0 && td->len[2] != 0)
		return 0;
	if (td->cnt[3] == 0 && td->len[3] != 0)
		return 0;
	if (td->cnt[0] > 20 && td->cnt[2] > 20)
		return 0;
	if (td->cnt[1] > 20 && td->cnt[3] > 20)
		return 0;
	return 1;
}
//---------------------------------------------------------------

static void
print_single_fastq(TLBDATA *tpp, unsigned int i, decodedSingle_p_t dsp)
{
	if (debug)
		printf("%lu0 %s\n%lu1 %s\n%lu2 +\n%lu3 %s\n",dsp->ordinal,dsp->hdr1,dsp->ordinal,dsp->tag1,dsp->ordinal,dsp->ordinal,dsp->qual1);
	else if (trim != NULL)
	{
		if (validTrimSingle(trim + dsp->ordinal))
		{
			unsigned int start = trim[dsp->ordinal].len[0] - 1;
			if (trim[dsp->ordinal].len[0] == 0)
				start = 0;
			unsigned int len = trim[dsp->ordinal].len[1] - start;
			if (trim[dsp->ordinal].len[1] == 0)
				len = dsp->taglen1 - start;
			printf("@%lu%s\n%.*s\n+\n%.*s\n",dsp->ordinal,dsp->hdr1,len,dsp->tag1 + start,len,dsp->qual1 + start);
		}
	}
	else
		printf("%s\n%s\n+\n%s\n",dsp->hdr1,dsp->tag1,dsp->qual1);
}
//---------------------------------------------------------------

static void
print_single_fetchGWI_p(TLBDATA *tpp, unsigned int i, decodedSingle_p_t dsp)
{
// AGTCCTGGGATTACAGGTGTGAGCCATTACGCCTGCCCAGTTCCTAATTT      1_1101_1228_2057        0;0     @@BDFDEFHHHHFJIIIGHIIJJJIJHIIIIIJJJIJJJIHGIIJIJIIJ      AGTCCTGGGATTACAGGTGTGAGCCATTACGCCTGCCCAGTTCCTAATTT      NC_000009.11[103202063..103202112]      +
	printf("%s\t%.*s\t0;0\t%s\t%s\t%s[%d..%d]\t%c\n", dsp->tag1, dsp->hhlen1, dsp->hdr1, dsp->qual1, dsp->tag1, virtchr[tpp->header.chr].AC, tpp->ppd[0][i].tag1pos, tpp->ppd[0][i].tag1pos + dsp->taglen1 - 1, dsp->reverseTAG1[0] ? '-' : '+');
}
//---------------------------------------------------------------

static void
print_single_fetchGWI_m(TLBDATA *tpp, unsigned int i, decodedSingle_p_t dsp)
{
	printf("%s\t%.*s\t0;0\t%s\t%s\t%s[%d..%d]\t%c\t%d\t%.*s\t%.*s\n",
				 dsp->tag1,
				 dsp->hhlen1, dsp->hdr1,
				 dsp->qual1,
				 dsp->tag1,
				 virtchr[tpp->header.chr].AC,
				 tpp->ppd[0][i].tag1pos,
				 tpp->ppd[0][i].tag1pos + dsp->taglen1 - 1,
				 dsp->reverseTAG1[0] ? '-' : '+',
				 /* FIXME - this is not yet correct... */
				 0,
				 dsp->taglen1, ":::",
				 dsp->taglen1, dsp->tag1);
}
//---------------------------------------------------------------

static void
print_single_fetchGWI_mN(TLBDATA *tpp, unsigned int i, decodedSingle_p_t dsp)
{
	printf("%s\t%.*s\t0;0\t%s\t%s\t%s[%d..%d]\t%c\n", dsp->tag1, dsp->hhlen1, dsp->hdr1, dsp->qual1, dsp->tag1, virtchr[tpp->header.chr].AC, tpp->ppd[0][i].tag1pos, tpp->ppd[0][i].tag1pos + dsp->taglen1 - 1, dsp->reverseTAG1[0] ? '-' : '+');
}
//---------------------------------------------------------------

static void
print_single_fetchGWI_g(TLBDATA *tpp, unsigned int i, decodedSingle_p_t dsp)
{
	// unimplemented
}
//---------------------------------------------------------------

static void
print_single_fetchGWI_u(TLBDATA *tpp, unsigned int i, decodedSingle_p_t dsp)
{
	// unimplemented
}
//---------------------------------------------------------------

static void
print_single_SAM_p(TLBDATA *tpp, unsigned int i, decodedSingle_p_t dsp)
{
	for (unsigned int j = 0; j < tpp->nbMatches; j++)
		printf("%.*s\t%d\t%s\t%d\t%d\t%dM\t=\t%d\t0\t%.*s\t%.*s\n",
					 dsp->hhlen1 - 1, dsp->hdr1 + 1,
					 dsp->flag1[j],
					 virtchr[dsp->chr[j]].SAMname,
					 tpp->ppd[j][i].tag1pos,
					 254,
					 dsp->taglen1,
					 tpp->ppd[j][i].tag1pos,
					 dsp->taglen1, dsp->tag1,
					 dsp->taglen1, dsp->qual1);
}
//---------------------------------------------------------------

static void
print_single_SAM_m(TLBDATA *tpp, unsigned int i, decodedSingle_p_t dsp)
{
	for (unsigned int j = 0; j < tpp->nbMatches; j++)
		printf("%.*s\t%d\t%s\t%d\t%d\t%dM\t=\t%d\t0\t%.*s\t%.*s\n",
					 dsp->hhlen1 - 1, dsp->hdr1 + 1,
					 dsp->flag1[j],
					 virtchr[dsp->chr[j]].SAMname,
					 tpp->ppd[j][i].tag1pos,
					 254 * (dsp->taglen1 - dsp->nbMismatch1) / dsp->taglen1,
					 dsp->taglen1,
					 tpp->ppd[j][i].tag1pos,
					 dsp->taglen1, dsp->tag1,
					 dsp->taglen1, dsp->qual1);
}
//---------------------------------------------------------------

static void
print_single_SAM_mN(TLBDATA *tpp, unsigned int i, decodedSingle_p_t dsp)
{
	for (unsigned int j = 0; j < tpp->nbMatches; j++)
		printf("%.*s\t%d\t%s\t%d\t%d\t%dM\t=\t%d\t0\t%.*s\t%.*s\n",
					 dsp->hhlen1 - 1, dsp->hdr1 + 1,
					 dsp->flag1[j],
					 virtchr[dsp->chr[j]].SAMname,
					 tpp->ppd[j][i].tag1pos,
					 254,
					 dsp->taglen1,
					 tpp->ppd[j][i].tag1pos,
					 dsp->taglen1, dsp->tag1,
					 dsp->taglen1, dsp->qual1);
}
//---------------------------------------------------------------

static void
print_single_SAM_g(TLBDATA *tpp, unsigned int i, decodedSingle_p_t dsp)
{
	for (unsigned int j = 0; j < tpp->nbMatches; j++)
	{
		// get the cigar strings if they exist
		unsigned int cigarpos_l = dsp->cigarposStart[j];
		if (tpp->cigar[j][cigarpos_l] != kCIGARTERMINATOR)	// if not, normal decoding.
		{
			char CigarString[2048];
			int cl = 0;
			// FIXME - probably not the right place for this, but for the time being... hack for correct output length
			int total = 0;
			while (tpp->cigar[j][cigarpos_l] != kCIGARTERMINATOR)
			{
				int pos = tpp->cigar[j][cigarpos_l++];
				char code = tpp->cigar[j][cigarpos_l++];
				if (code != 'D')
					total += pos;
				if (tpp->cigar[j][cigarpos_l] == kCIGARTERMINATOR && total < dsp->taglen1) {
					fprintf(stderr, "cigar is too short at %.*s : %s%d%c vs %d\n", dsp->hhlen1 - 1, dsp->hdr1 + 1, CigarString, pos,code, dsp->taglen1);
					pos += dsp->taglen1 - total;
				}
				cl += sprintf(&CigarString[cl],"%d%c",pos,code);
			}
			printf("%.*s\t%d\t%s\t%d\t%d\t%s \t=\t%d\t0\t%.*s\t%.*s\n",
						 dsp->hhlen1 - 1, dsp->hdr1 + 1,
						 dsp->flag1[j],
						 virtchr[dsp->chr[j]].SAMname,
						 tpp->ppd[j][i].tag1pos,
						 254 * (dsp->taglen1 - dsp->nbMismatch1) / dsp->taglen1,
						 CigarString,
						 tpp->ppd[j][i].tag1pos,
						 dsp->taglen1, dsp->tag1,
						 dsp->taglen1, dsp->qual1);
		}
		else
			printf("%.*s\t%d\t%s\t%d\t%d\t%dM\t=\t%d\t0\t%.*s\t%.*s\n",
						 dsp->hhlen1 - 1, dsp->hdr1 + 1,
						 dsp->flag1[j],
						 virtchr[dsp->chr[j]].SAMname,
						 tpp->ppd[j][i].tag1pos,
						 254 * (dsp->taglen1 - dsp->nbMismatch1) / dsp->taglen1,
						 dsp->taglen1,
						 tpp->ppd[j][i].tag1pos,
						 dsp->taglen1, dsp->tag1,
						 dsp->taglen1, dsp->qual1);
	}
}
//---------------------------------------------------------------

static void
print_single_SAM_u(TLBDATA *tpp, unsigned int i, decodedSingle_p_t dsp)
{
	printf("%.*s\t%d\t%s\t%d\t%d\t*\t=\t%d\t0\t%.*s\t%.*s\n",
				 dsp->hhlen1 - 1, dsp->hdr1 + 1,
				 dsp->flag1[0],
				 "*",
				 0,
				 255,
				 0,
				 dsp->taglen1, dsp->tag1,
				 dsp->taglen1, dsp->qual1);
}
//---------------------------------------------------------------

static void
print_single_ADNIview_p(TLBDATA *tpp, unsigned int i, decodedSingle_p_t dsp)
{
	printf("%.*s\t%s\t%s\t%s\t%s\t%u\t%d\t%d\t%d\t%lu\n", dsp->hhlen1 - 1, dsp->hdr1 + 1, dsp->tag1, dsp->qual1, dsp->tag1, dsp->qual1, tpp->ppd[0][i].tag1pos-1, dsp->reverseTAG1[0], dsp->taglen1,0,dsp->ordinal);
}
//---------------------------------------------------------------

static void
print_single_ADNIview_m(TLBDATA *tpp, unsigned int i, decodedSingle_p_t dsp)
{
	printf("%.*s\t%s\t%s\t%s\t%s\t%u\t%d\t%d\t%d\t%lu\n", dsp->hhlen1 - 1, dsp->hdr1 + 1, dsp->tag1, dsp->qual1, dsp->tag1, dsp->qual1, tpp->ppd[0][i].tag1pos-1, dsp->reverseTAG1[0], dsp->taglen1,0,dsp->ordinal);
}
//---------------------------------------------------------------

static void
print_single_ADNIview_mN(TLBDATA *tpp, unsigned int i, decodedSingle_p_t dsp)
{
	printf("%.*s\t%s\t%s\t%s\t%s\t%u\t%d\t%d\t%d\t%lu\n", dsp->hhlen1 - 1, dsp->hdr1 + 1, dsp->tag1, dsp->qual1, dsp->tag1, dsp->qual1, tpp->ppd[0][i].tag1pos-1, dsp->reverseTAG1[0], dsp->taglen1,0,dsp->ordinal);
}
//---------------------------------------------------------------

static void
print_single_ADNIview_g(TLBDATA *tpp, unsigned int i, decodedSingle_p_t dsp)
{
	char tag1[kMaxReadLen * 2]; // hold the tag1 revised with indels for display on genome
	char qual1[kMaxReadLen * 2]; // hold the tag1 quality revised with indels for display on genome
	unsigned int off1 = 0;
	AlignTagsOnGenomeUsingCigarSingle(tag1,qual1,dsp,tpp->cigar[0],&off1);
	printf("%.*s\t%s\t%s\t%s\t%s\t%u\t%d\t%d\t%d\t%lu\n", dsp->hhlen1 - 1, dsp->hdr1 + 1, tag1, qual1, dsp->tag1, dsp->qual1, tpp->ppd[0][i].tag1pos-1-off1, dsp->reverseTAG1[0], dsp->taglen1,0,dsp->ordinal);
}
//---------------------------------------------------------------

static void
print_single_ADNIview_u(TLBDATA *tpp, unsigned int i, decodedSingle_p_t dsp)
{
}
//---------------------------------------------------------------

static void
print_single_refless_p(TLBDATA *tpp, unsigned int i, decodedSingle_p_t dsp)
{
	if (trim != NULL)
	{
		if (validTrimSingle(trim + dsp->ordinal))
		{
			unsigned int start = trim[dsp->ordinal].len[0] - 1;
			if (trim[dsp->ordinal].len[0] == 0)
				start = 0;
			unsigned int len = trim[dsp->ordinal].len[1] - start;
			if (trim[dsp->ordinal].len[1] == 0)
				len = dsp->taglen1 - start;
			printf("%lu\t%up%u\t%.*s\t%.*s\n",dsp->ordinal,tpp->header.chr,tpp->ppd[0][i].tag1pos,len,dsp->tag1 + start,len,dsp->qual1 + start);
		}
	}
	else
		printf("%lu\t%up%u\t%s\t%s\t%u\t%u\t%u\t%u\n",dsp->ordinal,tpp->header.chr,tpp->ppd[0][i].tag1pos,dsp->tag1,dsp->qual1,tpp->header.chr,tpp->ppd[0][i].tag1pos,dsp->taglen1,dsp->reverseTAG1[0]);
}
//---------------------------------------------------------------

static void
print_single_refless_m(TLBDATA *tpp, unsigned int i, decodedSingle_p_t dsp)
{
	printf("%lu\t%um%u\t%s\t%s\t%u\t%u\t%u\t%u\t%u\n",dsp->ordinal,tpp->header.chr,tpp->ppd[0][i].tag1pos,dsp->tag1,dsp->qual1,tpp->header.chr,tpp->ppd[0][i].tag1pos,dsp->taglen1,dsp->reverseTAG1[0],dsp->nbMismatch1);
}
//---------------------------------------------------------------

static void
print_single_refless_mN(TLBDATA *tpp, unsigned int i, decodedSingle_p_t dsp)
{
	printf("%lu\t%uN%u\t%s\t%s\t%u\t%u\t%u\t%u\n",dsp->ordinal,tpp->header.chr,tpp->ppd[0][i].tag1pos,dsp->tag1,dsp->qual1,tpp->header.chr,tpp->ppd[0][i].tag1pos,dsp->taglen1,dsp->reverseTAG1[0]);
}
//---------------------------------------------------------------

static void
print_single_refless_g(TLBDATA *tpp, unsigned int i, decodedSingle_p_t dsp)
{
	printf("%lu\t%ug%u\t%s\t%s\n",dsp->ordinal,tpp->header.chr,tpp->ppd[0][i].tag1pos,dsp->tag1,dsp->qual1);
}
//---------------------------------------------------------------

static void
print_single_refless_u(TLBDATA *tpp, unsigned int i, decodedSingle_p_t dsp)
{
	printf("%lu\tU\t%s\t%s\n",dsp->ordinal,dsp->tag1,dsp->qual1);
}
//---------------------------------------------------------------

static void
print_single_contig_pmNg(TLBDATA *tpp, unsigned int i, decodedSingle_p_t dsp)
{
	if (!(tpp->ppd[0][i].tag1pos-200 > gFromPos && tpp->ppd[0][i].tag1pos < gToPos))
		return;
	printf("%lu\t%s\n",dsp->ordinal,dsp->reverseTAG1[0] ? "l" : "L");
}
//---------------------------------------------------------------

static void
print_single_contig_u(TLBDATA *tpp, unsigned int i, decodedSingle_p_t dsp)
{
}
//---------------------------------------------------------------

static void
print_single_tromer_p(TLBDATA *tpp, unsigned int i, decodedSingle_p_t dsp)
{
	for (unsigned int j = 0; j < tpp->nbMatches; j++)
	{
		if (j == 0)
			printf("%lu\t%u\tp\t%u\t%u\t%s",dsp->ordinal,dsp->chr[j],tpp->ppd[j][i].tag1pos,tpp->ppd[j][i].tag1pos+dsp->taglen1-1,dsp->tag1);
		else
			printf("\t%u\tp\t%u\t%u",dsp->chr[j],tpp->ppd[j][i].tag1pos,tpp->ppd[j][i].tag1pos+dsp->taglen1-1);
	}
	fputc('\n',stdout);
}
//---------------------------------------------------------------

static void
print_single_tromer_m(TLBDATA *tpp, unsigned int i, decodedSingle_p_t dsp)
{
	for (unsigned int j = 0; j < tpp->nbMatches; j++)
	{
		if (j == 0)
			printf("%lu\t%u\tm\t%u\t%u\t%s",dsp->ordinal,dsp->chr[j],tpp->ppd[j][i].tag1pos,tpp->ppd[j][i].tag1pos+dsp->taglen1-1,dsp->tag1);
		else
			printf("\t%u\tm\t%u\t%u",dsp->chr[j],tpp->ppd[j][i].tag1pos,tpp->ppd[j][i].tag1pos+dsp->taglen1-1);
	}
	fputc('\n',stdout);
}
//---------------------------------------------------------------

static void
print_single_tromer_mN(TLBDATA *tpp, unsigned int i, decodedSingle_p_t dsp)
{
	//printf("%lu\t%uN%u\t%s\t%s\n"    ,dsp->ordinal,tpp->header.chr,tpp->ppd[0][i].tag1pos,dsp->tag1,dsp->qual1);
	for (unsigned int j = 0; j < tpp->nbMatches; j++)
	{
		if (j == 0)
			printf("%lu\t%u\tN\t%u\t%u\t%s",dsp->ordinal,dsp->chr[j],tpp->ppd[j][i].tag1pos,tpp->ppd[j][i].tag1pos+dsp->taglen1-1,dsp->tag1);
		else
			printf("\t%u\tN\t%u\t%u",dsp->chr[j],tpp->ppd[j][i].tag1pos,tpp->ppd[j][i].tag1pos+dsp->taglen1-1);
	}
	fputc('\n',stdout);
}
//---------------------------------------------------------------

static void
print_single_tromer_g(TLBDATA *tpp, unsigned int i, decodedSingle_p_t dsp)
{
	//printf("%lu\t%ug%u\t%s\t%s\n"    ,dsp->ordinal,tpp->header.chr,tpp->ppd[0][i].tag1pos,dsp->tag1,dsp->qual1);
	for (unsigned int j = 0; j < tpp->nbMatches; j++)
	{
		if (j == 0)
			printf("%lu\t%u\tg\t%u\t%u\t%s",dsp->ordinal,dsp->chr[j],tpp->ppd[j][i].tag1pos,tpp->ppd[j][i].tag1pos+dsp->taglen1-1,dsp->tag1);
		else
			printf("\t%u\tg\t%u\t%u",dsp->chr[j],tpp->ppd[j][i].tag1pos,tpp->ppd[j][i].tag1pos+dsp->taglen1-1);
	}
	fputc('\n',stdout);
}
//---------------------------------------------------------------

static void
print_single_tromer_u(TLBDATA *tpp, unsigned int i, decodedSingle_p_t dsp)
{
	printf("%lu\t%u\tU\t%u\t%u\t%s\n",dsp->ordinal,0,0,0,dsp->tag1);
}
//---------------------------------------------------------------

#ifndef __APPLE__
static void
print_single_realign(TLBDATA *tpp, unsigned int i, decodedSingle_p_t dpp)
{
	fputs("Caramba olé, ca manque!!\n", stdout);
}
//---------------------------------------------------------------
#endif

static void
print_pair_fastq(TLBDATA *tpp, unsigned int i, decodedPair_p_t dpp)
{
	if (debug)
	{
		printf("%lu0 %s\n%lu1 %s\n%lu2 +\n%lu3 %s\n",dpp->ordinal,dpp->hdr1,dpp->ordinal,dpp->tag1,dpp->ordinal,dpp->ordinal,dpp->qual1);
		fprintf(stderr,"%lu0 %s\n%lu1 %s\n%lu2 +\n%lu3 %s\n",dpp->ordinal,dpp->hdr2,dpp->ordinal,dpp->tag2,dpp->ordinal,dpp->ordinal,dpp->qual2);
	}
	else if (trim != NULL)
	{
		if (validTrim(trim + dpp->ordinal))
		{
			unsigned int start = trim[dpp->ordinal].len[0] - 1;
			if (trim[dpp->ordinal].len[0] == 0)
				start = 0;
			unsigned int len = trim[dpp->ordinal].len[1] - start;
			if (trim[dpp->ordinal].len[1] == 0)
				len = dpp->taglen1 - start;
			printf("@%lu%s\n%.*s\n+\n%.*s\n",dpp->ordinal,dpp->hdr1,len,dpp->tag1 + start,len,dpp->qual1 + start);
			start = trim[dpp->ordinal].len[2] - 1;
			if (trim[dpp->ordinal].len[2] == 0)
				start = 0;
			len = trim[dpp->ordinal].len[3] - start;
			if (trim[dpp->ordinal].len[3] == 0)
				len = dpp->taglen2 - start;
			fprintf(stderr,"%s\n%.*s\n+\n%.*s\n",dpp->hdr2,len,dpp->tag2 + start,len,dpp->qual2 + start);
		}
	}
	else
	{
		printf("%s\n%s\n+\n%s\n",dpp->hdr1,dpp->tag1,dpp->qual1);
		printf("%s\n%s\n+\n%s\n",dpp->hdr2,dpp->tag2,dpp->qual2);
	}
}
//---------------------------------------------------------------

static void
print_pair_fetchGWI_p(TLBDATA *tpp, unsigned int i, decodedPair_p_t dpp)
{
	printf("%s\t%.*s\t0;0\t%s\t%s\t%s[%d..%d]\t%c\n", dpp->tag1, dpp->hhlen1, dpp->hdr1, dpp->qual1, dpp->tag1, virtchr[tpp->header.chr].AC, tpp->ppd[0][i].tag1pos, tpp->ppd[0][i].tag1pos + dpp->taglen1 - 1, dpp->reverseTAG1[0] ? '-' : '+');
	printf("%s\t%.*s\t0;0\t%s\t%s\t%s[%d..%d]\t%c\n", dpp->tag2, dpp->hhlen2, dpp->hdr2, dpp->qual2, dpp->tag2, virtchr[tpp->header.chr].AC, tpp->ppd[0][i].tag1pos + dpp->delta[0], tpp->ppd[0][i].tag1pos + dpp->delta[0] + dpp->taglen2 - 1, dpp->reverseTAG2[0] ? '-' : '+');
}
//---------------------------------------------------------------

static void
print_pair_fetchGWI_m(TLBDATA *tpp, unsigned int i, decodedPair_p_t dpp)
{
	printf("%s\t%.*s\t0;0\t%s\t%s\t%s[%d..%d]\t%c\n",
				 dpp->tag1,
				 dpp->hhlen1, dpp->hdr1,
				 dpp->qual1,
				 dpp->tag1,
				 virtchr[tpp->header.chr].AC,
				 tpp->ppd[0][i].tag1pos,
				 tpp->ppd[0][i].tag1pos + dpp->taglen1 - 1,
				 dpp->reverseTAG1[0] ? '-' : '+');
	printf("%s\t%.*s\t0;0\t%s\t%s\t%s[%d..%d]\t%c\n",
				 dpp->tag2,
				 dpp->hhlen2, dpp->hdr2,
				 dpp->qual2,
				 dpp->tag2,
				 virtchr[tpp->header.chr].AC,
				 tpp->ppd[0][i].tag1pos + dpp->delta[0],
				 tpp->ppd[0][i].tag1pos + dpp->delta[0] + dpp->taglen2 - 1,
				 dpp->reverseTAG2[0] ? '-' : '+');
//TTTAAAGCCTGTATTTGCAAGAGACGGTGCTTATTAAACTATAATCTAATATACTTATCAGTTTGGGCAGTGAAAATAAATTCTAAGTTTCCTATGCTTT	1_48_638_804.2	0;0;2	_^^cc^cc^eg`cfffhfhdhhhhhhhgadgbffhf__ef`fdgffhhhhfhhhgdfdghhaWb\aae_ggfbghggggbgdde_dbZ`bddb`dbbddb	TTTAAAGCCTGTATTTGCAAGAGACGGTGCTTATTAAACTATAATCTAATATACTTATCAGTTTGGACAGTGAAAATAAATTCTAAGTTTCCTTTGCTTT	NC_000021.8[15893561..15893660]	-	482	::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::x::::::::::::::::::::::::::x::::::	TTTAAAGCCTGTATTTGCAAGAGACGGTGCTTATTAAACTATAATCTAATATACTTATCAGTTTGGGCAGTGAAAATAAATTCTAAGTTTCCTATGCTTT
}
//---------------------------------------------------------------

static void
print_pair_fetchGWI_mN(TLBDATA *tpp, unsigned int i, decodedPair_p_t dpp)
{
	printf("%s\t%.*s\t0;0\t%s\t%s\t%s[%d..%d]\t%c\n", dpp->tag1, dpp->hhlen1, dpp->hdr1, dpp->qual1, dpp->tag1, virtchr[tpp->header.chr].AC, tpp->ppd[0][i].tag1pos, tpp->ppd[0][i].tag1pos + dpp->taglen1 - 1, dpp->reverseTAG1[0] ? '-' : '+');
	printf("%s\t%.*s\t0;0\t%s\t%s\t%s[%d..%d]\t%c\n", dpp->tag2, dpp->hhlen2, dpp->hdr2, dpp->qual2, dpp->tag2, virtchr[tpp->header.chr].AC, tpp->ppd[0][i].tag1pos + dpp->delta[0], tpp->ppd[0][i].tag1pos + dpp->delta[0] + dpp->taglen2 - 1, dpp->reverseTAG2[0] ? '-' : '+');
}
//---------------------------------------------------------------

static void
print_pair_fetchGWI_g(TLBDATA *tpp, unsigned int i, decodedPair_p_t dpp)
{
	// unimplemented
}
//---------------------------------------------------------------

static void
print_pair_fetchGWI_u(TLBDATA *tpp, unsigned int i, decodedPair_p_t dpp)
{
	// unimplemented
}
//---------------------------------------------------------------

static void
print_pair_SAM_p(TLBDATA *tpp, unsigned int i, decodedPair_p_t dpp)
{
	for (unsigned int j = 0; j < tpp->nbMatches; j++)
	{
		if (gFilterSAM
		    && (tpp->ppd[j][i].tag1pos + dpp->taglen1 >= virtchr[dpp->chr[j]].len
				    || tpp->ppd[j][i].tag1pos + dpp->delta[j] + dpp->taglen2 >= virtchr[dpp->chr[j]].len))
			continue;
		printf("%.*s\t%d\t%s\t%d\t%d\t%dM\t=\t%d\t%d\t%.*s\t%.*s\n",
					 dpp->hhlen1 - 1, dpp->hdr1 + 1,
					 dpp->flag1[j],
					 virtchr[dpp->chr[j]].SAMname,
					 tpp->ppd[j][i].tag1pos,
					 254,
					 dpp->taglen1,
					 tpp->ppd[j][i].tag1pos + dpp->delta[j],
					 dpp->ilen[j],
					 dpp->taglen1, dpp->tag1,
					 dpp->taglen1, dpp->qual1);
		printf("%.*s\t%d\t%s\t%d\t%d\t%dM\t=\t%d\t%d\t%.*s\t%.*s\n",
					 dpp->hhlen2 - 1, dpp->hdr2 + 1,
					 dpp->flag2[j],
					 virtchr[dpp->chr[j]].SAMname,
					 tpp->ppd[j][i].tag1pos + dpp->delta[j],
					 254,
					 dpp->taglen2,
					 tpp->ppd[j][i].tag1pos,
					 - dpp->ilen[j],
					 dpp->taglen2, dpp->tag2,
					 dpp->taglen2, dpp->qual2);
	}
}
//---------------------------------------------------------------

static void
print_pair_SAM_m(TLBDATA *tpp, unsigned int i, decodedPair_p_t dpp)
{
	for (unsigned int j = 0; j < tpp->nbMatches; j++)
	{
		// FIXME - could be better handled
		if (gFilterSAM
		    && (tpp->ppd[j][i].tag1pos == 0
						|| tpp->ppd[j][i].tag1pos + dpp->taglen1 >= virtchr[dpp->chr[j]].len
						|| tpp->ppd[j][i].tag1pos + dpp->delta[j] <= 0
				    || tpp->ppd[j][i].tag1pos + dpp->delta[j] + dpp->taglen2 >= virtchr[dpp->chr[j]].len))
			continue;
		printf("%.*s\t%d\t%s\t%d\t%d\t%dM\t=\t%d\t%d\t%.*s\t%.*s\n",
					 dpp->hhlen1 - 1, dpp->hdr1 + 1,
					 dpp->flag1[j],
					 virtchr[dpp->chr[j]].SAMname,
					 tpp->ppd[j][i].tag1pos,
					 254 * (dpp->taglen1 - dpp->nbMismatch1) / dpp->taglen1,
					 dpp->taglen1,
					 tpp->ppd[j][i].tag1pos + dpp->delta[j],
					 dpp->ilen[j],
					 dpp->taglen1, dpp->tag1,
					 dpp->taglen1, dpp->qual1);
		printf("%.*s\t%d\t%s\t%d\t%d\t%dM\t=\t%d\t%d\t%.*s\t%.*s\n",
					 dpp->hhlen2 - 1, dpp->hdr2 + 1,
					 dpp->flag2[j],
					 virtchr[dpp->chr[j]].SAMname,
					 tpp->ppd[j][i].tag1pos + dpp->delta[j],
					 254 * (dpp->taglen2 - dpp->nbMismatch2) / dpp->taglen2,
					 dpp->taglen2,
					 tpp->ppd[j][i].tag1pos,
					 - dpp->ilen[j],
					 dpp->taglen2, dpp->tag2,
					 dpp->taglen2, dpp->qual2);
	}
}
//---------------------------------------------------------------

static void
print_pair_SAM_mN(TLBDATA *tpp, unsigned int i, decodedPair_p_t dpp)
{
	for (unsigned int j = 0; j < tpp->nbMatches; j++)
	{
		printf("%.*s\t%d\t%s\t%d\t%d\t%dM\t=\t%d\t%d\t%.*s\t%.*s\n",
					 dpp->hhlen1 - 1, dpp->hdr1 + 1,
					 dpp->flag1[j],
					 virtchr[dpp->chr[j]].SAMname,
					 tpp->ppd[j][i].tag1pos,
					 254,
					 dpp->taglen1,
					 tpp->ppd[j][i].tag1pos + dpp->delta[j],
					 dpp->ilen[j],
					 dpp->taglen1, dpp->tag1,
					 dpp->taglen1, dpp->qual1);
		printf("%.*s\t%d\t%s\t%d\t%d\t%dM\t=\t%d\t%d\t%.*s\t%.*s\n",
					 dpp->hhlen2 - 1, dpp->hdr2 + 1,
					 dpp->flag2[j],
					 virtchr[dpp->chr[j]].SAMname,
					 tpp->ppd[j][i].tag1pos + dpp->delta[j],
					 254,
					 dpp->taglen2,
					 tpp->ppd[j][i].tag1pos,
					 - dpp->ilen[j],
					 dpp->taglen2, dpp->tag2,
					 dpp->taglen2, dpp->qual2);
	}
}
//---------------------------------------------------------------

static void
print_pair_SAM_g(TLBDATA *tpp, unsigned int i, decodedPair_p_t dpp)
{
	for (unsigned int j = 0; j < tpp->nbMatches; j++)
	{
		// FIXME - should compute total of deletions instead of just using 255 as a safety margin...
		if (gFilterSAM
		    && (tpp->ppd[j][i].tag1pos == 0
						|| tpp->ppd[j][i].tag1pos + dpp->taglen1 + 255 >= virtchr[dpp->chr[j]].len
						|| tpp->ppd[j][i].tag1pos + dpp->delta[j] <= 0
				    || tpp->ppd[j][i].tag1pos + dpp->delta[j] + dpp->taglen2 + 255 >= virtchr[dpp->chr[j]].len))
			continue;
		// get the cigar strings if they exist
		unsigned int cigarpos_l = dpp->cigarposStart[j];
		if (tpp->cigar[j][cigarpos_l] != kCIGARTERMINATOR)	// if not, normal decoding.
		{
			char CigarString[2048];
			int cl = 0;
			// FIXME - probably not the right place for this, but for the time being... hack for correct output length
			int total = 0;
			while (tpp->cigar[j][cigarpos_l] != kCIGARTERMINATOR)
			{
				int pos = tpp->cigar[j][cigarpos_l++];
				char code = tpp->cigar[j][cigarpos_l++];
				if (code != 'D')
					total += pos;
				if (tpp->cigar[j][cigarpos_l] == kCIGARTERMINATOR && total < dpp->taglen1) {
					fprintf(stderr, "cigar is too short at %.*s : %s%d%c vs %d\n", dpp->hhlen1 - 1, dpp->hdr1 + 1, CigarString, pos,code, dpp->taglen1);
					pos += dpp->taglen1 - total;
				}
				cl += sprintf(&CigarString[cl],"%d%c",pos,code);
			}
			printf("%.*s\t%d\t%s\t%d\t%d\t%s\t=\t%d\t%d\t%.*s\t%.*s\n",
						 dpp->hhlen1 - 1, dpp->hdr1 + 1,
						 dpp->flag1[j],
						 virtchr[dpp->chr[j]].SAMname,
						 tpp->ppd[j][i].tag1pos,
						 254 * (dpp->taglen1 - dpp->nbMismatch1) / dpp->taglen1,
						 CigarString,
						 tpp->ppd[j][i].tag1pos + dpp->delta[j],
						 dpp->ilen[j],
						 dpp->taglen1, dpp->tag1,
						 dpp->taglen1, dpp->qual1);
		}
		else
			printf("%.*s\t%d\t%s\t%d\t%d\t%dM\t=\t%d\t%d\t%.*s\t%.*s\n",
						 dpp->hhlen1 - 1, dpp->hdr1 + 1,
						 dpp->flag1[j],
						 virtchr[dpp->chr[j]].SAMname,
						 tpp->ppd[j][i].tag1pos,
						 254 * (dpp->taglen1 - dpp->nbMismatch1) / dpp->taglen1,
						 dpp->taglen1,
						 tpp->ppd[j][i].tag1pos + dpp->delta[j],
						 dpp->ilen[j],
						 dpp->taglen1, dpp->tag1,
						 dpp->taglen1, dpp->qual1);
		cigarpos_l += 1;
		if (tpp->cigar[j][cigarpos_l] != kCIGARTERMINATOR)	// if not, normal decoding.
		{
			char CigarString[2048];
			int cl = 0;
			// FIXME - probably not the right place for this, but for the time being... hack for correct output length
			int total = 0;
			while (tpp->cigar[j][cigarpos_l] != kCIGARTERMINATOR)
			{
				int pos = tpp->cigar[j][cigarpos_l++];
				char code = tpp->cigar[j][cigarpos_l++];
				if (code != 'D')
					total += pos;
				if (tpp->cigar[j][cigarpos_l] == kCIGARTERMINATOR && total < dpp->taglen2) {
					fprintf(stderr, "cigar is too short at %.*s : %s%d%c vs %d\n", dpp->hhlen2 - 1, dpp->hdr2 + 1, CigarString, pos,code, dpp->taglen2);
					pos += dpp->taglen2 - total;
				}
				cl += sprintf(&CigarString[cl],"%d%c",pos,code);
			}
			printf("%.*s\t%d\t%s\t%d\t%d\t%s\t=\t%d\t%d\t%.*s\t%.*s\n",
						 dpp->hhlen2 - 1, dpp->hdr2 + 1,
						 dpp->flag2[j],
						 virtchr[dpp->chr[j]].SAMname,
						 tpp->ppd[j][i].tag1pos + dpp->delta[j],
						 254 * (dpp->taglen2 - dpp->nbMismatch2) / dpp->taglen2,
						 CigarString,
						 tpp->ppd[j][i].tag1pos,
						 - dpp->ilen[j],
						 dpp->taglen2, dpp->tag2,
						 dpp->taglen2, dpp->qual2);
		}
		else
			printf("%.*s\t%d\t%s\t%d\t%d\t%dM\t=\t%d\t%d\t%.*s\t%.*s\n",
						 dpp->hhlen2 - 1, dpp->hdr2 + 1,
						 dpp->flag2[j],
						 virtchr[dpp->chr[j]].SAMname,
						 tpp->ppd[j][i].tag1pos + dpp->delta[j],
						 254 * (dpp->taglen2 - dpp->nbMismatch2) / dpp->taglen2,
						 dpp->taglen2,
						 tpp->ppd[j][i].tag1pos,
						 - dpp->ilen[j],
						 dpp->taglen2, dpp->tag2,
						 dpp->taglen2, dpp->qual2);
	}
}
//---------------------------------------------------------------

static void
print_pair_SAM_u(TLBDATA *tpp, unsigned int i, decodedPair_p_t dpp)
{
	if ((tpp->header.flags_7_0 & kTLBHflagPairChrPosBlock) == 0)
	{
		printf("%.*s\t%d\t%s\t%d\t%d\t*\t*\t%d\t%d\t%.*s\t%.*s\n",
					 dpp->hhlen1 - 1, dpp->hdr1 + 1,
					 dpp->flag1[0],
					 "*",
					 0,
					 255,
					 0,
					 0,
					 dpp->taglen1, dpp->tag1,
					 dpp->taglen1, dpp->qual1);
		printf("%.*s\t%d\t%s\t%d\t%d\t*\t*\t%d\t%d\t%.*s\t%.*s\n",
					 dpp->hhlen2 - 1, dpp->hdr2 + 1,
					 dpp->flag2[0],
					 "*",
					 0,
					 255,
					 0,
					 0,
					 dpp->taglen2, dpp->tag2,
					 dpp->taglen2, dpp->qual2);
	}
	else
		for (unsigned int j = 0; j < tpp->nbMatches; j++)
		{
			char cigar1[8],cigar2[8];
			unsigned int chr1 = GET_TAG1_CHR(tpp->pcpd[j][i].chr);
			unsigned int chr2 = GET_TAG2_CHR(tpp->pcpd[j][i].chr);
			// FIXME - do not output those 0 coordinates... could use a better fix
			if (gFilterSAM
					&& ((chr1 != 0 && (tpp->pcpd[j][i].tag1pos == 0 || tpp->pcpd[j][i].tag1pos + dpp->taglen1 >= virtchr[chr1].len))
							|| (chr2 != 0 && (tpp->pcpd[j][i].tag2pos == 0 || tpp->pcpd[j][i].tag2pos + dpp->taglen2 >= virtchr[chr2].len))))
				continue;
			if (chr1 != 0)
				snprintf(cigar1,8,"%dM",dpp->taglen1);
			else
				strcpy(cigar1,"*");
			if (chr2 != 0)
				snprintf(cigar2,8,"%dM",dpp->taglen2);
			else
				strcpy(cigar2,"*");
			printf("%.*s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%.*s\t%.*s\n",
						 dpp->hhlen1 - 1, dpp->hdr1 + 1,
						 dpp->flag1[j],
						 chr1 == 0 ? "*" : virtchr[chr1].SAMname,
						 tpp->pcpd[j][i].tag1pos,
						 chr1 == 0 ? 0 : 254 * (dpp->taglen1 - dpp->nbMismatch1) / dpp->taglen1,
						 cigar1,
						 chr2 == 0 ? "*" : virtchr[chr2].SAMname,
						 tpp->pcpd[j][i].tag2pos,
						 dpp->ilen[j],
						 dpp->taglen1, dpp->tag1,
						 dpp->taglen1, dpp->qual1);
			printf("%.*s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%.*s\t%.*s\n",
						 dpp->hhlen2 - 1, dpp->hdr2 + 1,
						 dpp->flag2[j],
						 chr2 == 0 ? "*" : virtchr[chr2].SAMname,
						 tpp->pcpd[j][i].tag2pos,
						 chr2 == 0 ? 0 : 254 * (dpp->taglen2 - dpp->nbMismatch2) / dpp->taglen2,
						 cigar2,
						 chr1 == 0 ? "*" : virtchr[chr1].SAMname,
						 tpp->pcpd[j][i].tag1pos,
						 dpp->ilen[j],
						 dpp->taglen2, dpp->tag2,
						 dpp->taglen2, dpp->qual2);
		}
}
//---------------------------------------------------------------

static void
print_pair_ADNIview_p(TLBDATA *tpp, unsigned int i, decodedPair_p_t dpp)
{
	int leftTagPos,rightTagPos;
	int delta = (tpp->ppd[0][i].tag2offset & k32bOFFSET_VALUE_MASK);
	if (tpp->ppd[0][i].tag2offset & k32bNEGATIVE_OFFSET_BIT)
		delta = -delta;

	if (delta >= 0)
	{
		leftTagPos = tpp->ppd[0][i].tag1pos;
		rightTagPos = tpp->ppd[0][i].tag1pos + delta;
	}
	else
	{
		rightTagPos = tpp->ppd[0][i].tag1pos;
		leftTagPos = tpp->ppd[0][i].tag1pos + delta;
	}

	if (!((leftTagPos-200 > gFromPos && leftTagPos < gToPos)
		|| (rightTagPos-200 > gFromPos && rightTagPos < gToPos)))
		return;

	// FIXME -- strange ??
	int orientation = 0;

	if (delta >= 0)
	{
		printf("%.*s\t%s\t%s\t%s\t%s\t%u\t%d\t%d\t%d\t%lu\n",
		dpp->hhlen1 - 1, dpp->hdr1 + 1, dpp->tag1, dpp->qual1, dpp->tag1, dpp->qual1, tpp->ppd[0][i].tag1pos-1, kIsPairedEnd, dpp->taglen1,orientation,dpp->ordinal);
	}
	printf("%.*s\t%s\t%s\t%s\t%s\t%u\t%d\t%d\t%d\t%lu\n",
	dpp->hhlen2 - 1, dpp->hdr2 + 1, dpp->tag2, dpp->qual2, dpp->tag2, dpp->qual2, tpp->ppd[0][i].tag1pos-1+delta, kIsPairedEnd, dpp->taglen2,orientation,dpp->ordinal);
	if (delta < 0)
	{
		printf("%.*s\t%s\t%s\t%s\t%s\t%u\t%d\t%d\t%d\t%lu\n",
		dpp->hhlen1 - 1, dpp->hdr1 + 1, dpp->tag1, dpp->qual1, dpp->tag1, dpp->qual1, tpp->ppd[0][i].tag1pos-1, kIsPairedEnd, dpp->taglen1,orientation,dpp->ordinal);
	}
}
//---------------------------------------------------------------

static void
print_pair_ADNIview_m(TLBDATA *tpp, unsigned int i, decodedPair_p_t dpp)
{
	int leftTagPos,rightTagPos;
	int delta = (tpp->ppd[0][i].tag2offset & k32bOFFSET_VALUE_MASK);
	if (tpp->ppd[0][i].tag2offset & k32bNEGATIVE_OFFSET_BIT)
		delta = -delta;

	if (delta >= 0)
	{
		leftTagPos = tpp->ppd[0][i].tag1pos;
		rightTagPos = tpp->ppd[0][i].tag1pos + delta;
	}
	else
	{
		rightTagPos = tpp->ppd[0][i].tag1pos;
		leftTagPos = tpp->ppd[0][i].tag1pos + delta;
	}

	if (!((leftTagPos-200 > gFromPos && leftTagPos < gToPos)
		|| (rightTagPos-200 > gFromPos && rightTagPos < gToPos)))
		return;

	//!!! FIXME: not necessary to remove /1 /2 should be done in decodeHeaderPair ???  to fix in all ADNIview.
	// FIXME  orientation is always 0.... does not affect display, but...
	int orientation = 0;

	if (delta >= 0)
	{
		printf("%.*s\t%s\t%s\t%s\t%s\t%u\t%d\t%d\t%d\t%lu\n",
		dpp->hhlen1 - 1, dpp->hdr1 + 1, dpp->tag1, dpp->qual1, dpp->tag1, dpp->qual1, tpp->ppd[0][i].tag1pos-1, kIsPairedEnd, dpp->taglen1,orientation,dpp->ordinal);
	}
	printf("%.*s\t%s\t%s\t%s\t%s\t%u\t%d\t%d\t%d\t%lu\n",
	dpp->hhlen2 - 1, dpp->hdr2 + 1, dpp->tag2, dpp->qual2, dpp->tag2, dpp->qual2, tpp->ppd[0][i].tag1pos-1+delta, kIsPairedEnd, dpp->taglen2,orientation,dpp->ordinal);
	if (delta < 0)
	{
		printf("%.*s\t%s\t%s\t%s\t%s\t%u\t%d\t%d\t%d\t%lu\n",
		dpp->hhlen1 - 1, dpp->hdr1 + 1, dpp->tag1, dpp->qual1, dpp->tag1, dpp->qual1, tpp->ppd[0][i].tag1pos-1, kIsPairedEnd, dpp->taglen1,orientation,dpp->ordinal);
	}
}
//---------------------------------------------------------------

static void
print_pair_ADNIview_mN(TLBDATA *tpp, unsigned int i, decodedPair_p_t dpp)
{
	int leftTagPos,rightTagPos;
	int delta = (tpp->ppd[0][i].tag2offset & k32bOFFSET_VALUE_MASK);
	if (tpp->ppd[0][i].tag2offset & k32bNEGATIVE_OFFSET_BIT)
		delta = -delta;

	if (delta >= 0)
	{
		leftTagPos = tpp->ppd[0][i].tag1pos;
		rightTagPos = tpp->ppd[0][i].tag1pos + delta;
	}
	else
	{
		rightTagPos = tpp->ppd[0][i].tag1pos;
		leftTagPos = tpp->ppd[0][i].tag1pos + delta;
	}

	if (!((leftTagPos-200 > gFromPos && leftTagPos < gToPos)
		|| (rightTagPos-200 > gFromPos && rightTagPos < gToPos)))
		return;

	int orientation = 0;

	if (delta >= 0)
	{
		printf("%.*s\t%s\t%s\t%s\t%s\t%u\t%d\t%d\t%d\t%lu\n",
		dpp->hhlen1 - 1, dpp->hdr1 + 1, dpp->tag1, dpp->qual1, dpp->tag1, dpp->qual1, tpp->ppd[0][i].tag1pos-1, kIsPairedEnd, dpp->taglen1,orientation,dpp->ordinal);
	}
	printf("%.*s\t%s\t%s\t%s\t%s\t%u\t%d\t%d\t%d\t%lu\n",
	dpp->hhlen2 - 1, dpp->hdr2 + 1, dpp->tag2, dpp->qual2, dpp->tag2, dpp->qual2, tpp->ppd[0][i].tag1pos-1+delta, kIsPairedEnd, dpp->taglen2,orientation,dpp->ordinal);
	if (delta < 0)
	{
		printf("%.*s\t%s\t%s\t%s\t%s\t%u\t%d\t%d\t%d\t%lu\n",
		dpp->hhlen1 - 1, dpp->hdr1 + 1, dpp->tag1, dpp->qual1, dpp->tag1, dpp->qual1, tpp->ppd[0][i].tag1pos-1, kIsPairedEnd, dpp->taglen1,orientation,dpp->ordinal);
	}
}
//---------------------------------------------------------------

static void
print_pair_ADNIview_g(TLBDATA *tpp, unsigned int i, decodedPair_p_t dpp)
{
	char tag1[kMaxReadLen * 2]; // hold the tag1 revised with indels for display on genome
	char tag2[kMaxReadLen * 2]; // hold the tag2 revised with indels for display on genome
	char qual1[kMaxReadLen * 2]; // hold the tag1 quality revised with indels for display on genome
	char qual2[kMaxReadLen * 2]; // hold the tag2 quality revised with indels for display on genome
	int leftTagPos,rightTagPos;
	int delta = (tpp->ppd[0][i].tag2offset & k32bOFFSET_VALUE_MASK);
	if (tpp->ppd[0][i].tag2offset & k32bNEGATIVE_OFFSET_BIT)
		delta = -delta;

	if (delta >= 0)
	{
		leftTagPos = tpp->ppd[0][i].tag1pos;
		rightTagPos = tpp->ppd[0][i].tag1pos + delta;
	}
	else
	{
		rightTagPos = tpp->ppd[0][i].tag1pos;
		leftTagPos = tpp->ppd[0][i].tag1pos + delta;
	}

	if (!((leftTagPos-200 > gFromPos && leftTagPos < gToPos)
		|| (rightTagPos-200 > gFromPos && rightTagPos < gToPos)))
		return;

	int orientation = 0;
	unsigned int off1 = 0;
	unsigned int off2 = 0;

	AlignTagsOnGenomeUsingCigar(tag1,tag2,qual1,qual2,dpp,tpp->cigar[0],&off1,&off2);
	if (delta >= 0)
	{
		printf("%.*s\t%s\t%s\t%s\t%s\t%u\t%d\t%d\t%d\t%lu\n",
		dpp->hhlen1 - 1, dpp->hdr1 + 1, tag1, qual1, dpp->tag1, dpp->qual1, tpp->ppd[0][i].tag1pos-1-off1, kIsPairedEnd, dpp->taglen1,orientation,dpp->ordinal);
	}
	printf("%.*s\t%s\t%s\t%s\t%s\t%u\t%d\t%d\t%d\t%lu\n",
	dpp->hhlen2 - 1, dpp->hdr2 + 1, tag2, qual2, dpp->tag2, dpp->qual2, tpp->ppd[0][i].tag1pos-1+delta-off2, kIsPairedEnd, dpp->taglen2,orientation,dpp->ordinal);
	if (delta < 0)
	{
		printf("%.*s\t%s\t%s\t%s\t%s\t%u\t%d\t%d\t%d\t%lu\n",
		dpp->hhlen1 - 1, dpp->hdr1 + 1, tag1, qual1, dpp->tag1, dpp->qual1, tpp->ppd[0][i].tag1pos-1-off1, kIsPairedEnd, dpp->taglen1,orientation,dpp->ordinal);
	}
}
//---------------------------------------------------------------

static void
print_pair_ADNIview_u(TLBDATA *tpp, unsigned int i, decodedPair_p_t dpp)
{
}
//---------------------------------------------------------------

static void
print_pair_refless_p(TLBDATA *tpp, unsigned int i, decodedPair_p_t dpp)
{
	if (trim != NULL)
	{
		if (validTrim(trim + dpp->ordinal))
		{
			unsigned int start1 = trim[dpp->ordinal].len[0] - 1;
			if (trim[dpp->ordinal].len[0] == 0)
				start1 = 0;
			unsigned int len1 = trim[dpp->ordinal].len[1] - start1;
			if (trim[dpp->ordinal].len[1] == 0)
				len1 = dpp->taglen1 - start1;
			unsigned int start2 = trim[dpp->ordinal].len[2] - 1;
			if (trim[dpp->ordinal].len[2] == 0)
				start2 = 0;
			unsigned int len2 = trim[dpp->ordinal].len[3] - start2;
			if (trim[dpp->ordinal].len[3] == 0)
				len2 = dpp->taglen2 - start2;
			printf("%lu\t%up%u\t%.*s\t%.*s\t%.*s\t%.*s\n",dpp->ordinal,tpp->header.chr,tpp->ppd[0][i].tag1pos,len1,dpp->tag1 + start1,len2,dpp->tag2 + start2,len1,dpp->qual1 + start1,len2,dpp->qual2 + start2);
		}
	}
	else
	{
		unsigned int start = tpp->ppd[0][i].tag1pos;
		unsigned int end = start + dpp->taglen1;
		unsigned int tag2pos = tpp->ppd[0][i].tag1pos + dpp->delta[0];
		if (tag2pos < start)
			start = tag2pos;
		if (tag2pos + dpp->taglen2 > end)
			end = tag2pos + dpp->taglen2;
		printf("%lu\t%up%u\t%s\t%s\t%s\t%s\t%u\t%u\t%u\t%u\t%u\n",dpp->ordinal,tpp->header.chr,start,dpp->tag1,dpp->tag2,dpp->qual1,dpp->qual2,tpp->header.chr,start,end - start,dpp->reverseTAG1[0],dpp->reverseTAG2[0]);
	}
}
//---------------------------------------------------------------

static void
print_pair_refless_m(TLBDATA *tpp, unsigned int i, decodedPair_p_t dpp)
{
	if (trim != NULL)
	{
		if (validTrim(trim + dpp->ordinal))
		{
			unsigned int start1 = trim[dpp->ordinal].len[0] - 1;
			if (trim[dpp->ordinal].len[0] == 0)
				start1 = 0;
			unsigned int len1 = trim[dpp->ordinal].len[1] - start1;
			if (trim[dpp->ordinal].len[1] == 0)
				len1 = dpp->taglen1 - start1;
			unsigned int start2 = trim[dpp->ordinal].len[2] - 1;
			if (trim[dpp->ordinal].len[2] == 0)
				start2 = 0;
			unsigned int len2 = trim[dpp->ordinal].len[3] - start2;
			if (trim[dpp->ordinal].len[3] == 0)
				len2 = dpp->taglen2 - start2;
			printf("%lu\t%um%u\t%.*s\t%.*s\t%.*s\t%.*s\n",dpp->ordinal,tpp->header.chr,tpp->ppd[0][i].tag1pos,len1,dpp->tag1 + start1,len2,dpp->tag2 + start2,len1,dpp->qual1 + start1,len2,dpp->qual2 + start2);
		}
	}
	else
	{
		unsigned int start = tpp->ppd[0][i].tag1pos;
		unsigned int end = start + dpp->taglen1;
		unsigned int tag2pos = tpp->ppd[0][i].tag1pos + dpp->delta[0];
		if (tag2pos < start)
			start = tag2pos;
		if (tag2pos + dpp->taglen2 > end)
			end = tag2pos + dpp->taglen2;
		printf("%lu\t%um%u\t%s\t%s\t%s\t%s\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n",dpp->ordinal,tpp->header.chr,start,dpp->tag1,dpp->tag2,dpp->qual1,dpp->qual2,tpp->header.chr,start,end - start,dpp->reverseTAG1[0],dpp->reverseTAG2[0],dpp->nbMismatch1,dpp->nbMismatch2);
	}
}
//---------------------------------------------------------------

static void
print_pair_refless_mN(TLBDATA *tpp, unsigned int i, decodedPair_p_t dpp)
{
	if (trim != NULL)
	{
		if (validTrim(trim + dpp->ordinal))
		{
			unsigned int start1 = trim[dpp->ordinal].len[0] - 1;
			if (trim[dpp->ordinal].len[0] == 0)
				start1 = 0;
			unsigned int len1 = trim[dpp->ordinal].len[1] - start1;
			if (trim[dpp->ordinal].len[1] == 0)
				len1 = dpp->taglen1 - start1;
			unsigned int start2 = trim[dpp->ordinal].len[2] - 1;
			if (trim[dpp->ordinal].len[2] == 0)
				start2 = 0;
			unsigned int len2 = trim[dpp->ordinal].len[3] - start2;
			if (trim[dpp->ordinal].len[3] == 0)
				len2 = dpp->taglen2 - start2;
			printf("%lu\t%uN%u\t%.*s\t%.*s\t%.*s\t%.*s\n",dpp->ordinal,tpp->header.chr,tpp->ppd[0][i].tag1pos,len1,dpp->tag1 + start1,len2,dpp->tag2 + start2,len1,dpp->qual1 + start1,len2,dpp->qual2 + start2);
		}
	}
	else
	{
		unsigned int start = tpp->ppd[0][i].tag1pos;
		unsigned int end = start + dpp->taglen1;
		unsigned int tag2pos = tpp->ppd[0][i].tag1pos + dpp->delta[0];
		if (tag2pos < start)
			start = tag2pos;
		if (tag2pos + dpp->taglen2 > end)
			end = tag2pos + dpp->taglen2;
		printf("%lu\t%uN%u\t%s\t%s\t%s\t%s\t%u\t%u\t%u\t%u\t%u\n",dpp->ordinal,tpp->header.chr,start,dpp->tag1,dpp->tag2,dpp->qual1,dpp->qual2,tpp->header.chr,start,end - start,dpp->reverseTAG1[0],dpp->reverseTAG2[0]);
	}
}
//---------------------------------------------------------------

static void
print_pair_refless_g(TLBDATA *tpp, unsigned int i, decodedPair_p_t dpp)
{
	if (trim != NULL)
	{
		if (validTrim(trim + dpp->ordinal))
		{
			unsigned int start1 = trim[dpp->ordinal].len[0] - 1;
			if (trim[dpp->ordinal].len[0] == 0)
				start1 = 0;
			unsigned int len1 = trim[dpp->ordinal].len[1] - start1;
			if (trim[dpp->ordinal].len[1] == 0)
				len1 = dpp->taglen1 - start1;
			unsigned int start2 = trim[dpp->ordinal].len[2] - 1;
			if (trim[dpp->ordinal].len[2] == 0)
				start2 = 0;
			unsigned int len2 = trim[dpp->ordinal].len[3] - start2;
			if (trim[dpp->ordinal].len[3] == 0)
				len2 = dpp->taglen2 - start2;
			printf("%lu\t%ug%u\t%.*s\t%.*s\t%.*s\t%.*s\n",dpp->ordinal,tpp->header.chr,tpp->ppd[0][i].tag1pos,len1,dpp->tag1 + start1,len2,dpp->tag2 + start2,len1,dpp->qual1 + start1,len2,dpp->qual2 + start2);
		}
	}
	else
	{
		unsigned int start = tpp->ppd[0][i].tag1pos;
		unsigned int end = start + dpp->taglen1;
		unsigned int tag2pos = tpp->ppd[0][i].tag1pos + dpp->delta[0];
		if (tag2pos < start)
			start = tag2pos;
		if (tag2pos + dpp->taglen2 > end)
			end = tag2pos + dpp->taglen2;
		// get the cigar strings if they exist
		unsigned int cigarpos_l = dpp->cigarposStart[0];
		char CigarString1[2048];
		char CigarString2[2048];
		if (tpp->cigar[0][cigarpos_l] != kCIGARTERMINATOR)	// if not, normal decoding.
		{
			int cl = 0;
			while (tpp->cigar[0][cigarpos_l] != kCIGARTERMINATOR)
			{
				int pos = tpp->cigar[0][cigarpos_l++];
				char code = tpp->cigar[0][cigarpos_l++];
				cl += sprintf(CigarString1 + cl,"%d%c",pos,code);
			}
		}
		else
			sprintf(CigarString1,"%dM",dpp->taglen1);
		cigarpos_l += 1;
		if (tpp->cigar[0][cigarpos_l] != kCIGARTERMINATOR)	// if not, normal decoding.
		{
			int cl = 0;
			while (tpp->cigar[0][cigarpos_l] != kCIGARTERMINATOR)
			{
				int pos = tpp->cigar[0][cigarpos_l++];
				char code = tpp->cigar[0][cigarpos_l++];
				cl += sprintf(CigarString2 + cl,"%d%c",pos,code);
			}
		}
		else
			sprintf(CigarString2,"%dM",dpp->taglen2);
		printf("%lu\t%ug%u\t%s\t%s\t%s\t%s\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%s\t%s\n",dpp->ordinal,tpp->header.chr,start,dpp->tag1,dpp->tag2,dpp->qual1,dpp->qual2,tpp->header.chr,start,end - start,dpp->reverseTAG1[0],dpp->reverseTAG2[0],dpp->nbMismatch1,dpp->nbMismatch2,CigarString1,CigarString2);
	}
}
//---------------------------------------------------------------

static void
print_pair_refless_u(TLBDATA *tpp, unsigned int i, decodedPair_p_t dpp)
{
	if (trim != NULL)
	{
		if (validTrim(trim + dpp->ordinal))
		{
			unsigned int start1 = trim[dpp->ordinal].len[0] - 1;
			if (trim[dpp->ordinal].len[0] == 0)
				start1 = 0;
			unsigned int len1 = trim[dpp->ordinal].len[1] - start1;
			if (trim[dpp->ordinal].len[1] == 0)
				len1 = dpp->taglen1 - start1;
			unsigned int start2 = trim[dpp->ordinal].len[2] - 1;
			if (trim[dpp->ordinal].len[2] == 0)
				start2 = 0;
			unsigned int len2 = trim[dpp->ordinal].len[3] - start2;
			if (trim[dpp->ordinal].len[3] == 0)
				len2 = dpp->taglen2 - start2;
			if ((tpp->header.flags_7_0 & kTLBHflagPairChrPosBlock) == 0)
				printf("%lu\tU\t%.*s\t%.*s\t%.*s\t%.*s\n",dpp->ordinal,len1,dpp->tag1 + start1,len2,dpp->tag2 + start2,len1,dpp->qual1 + start1,len2,dpp->qual2 + start2);
			else
			{
				unsigned int chr1 = GET_TAG1_CHR(tpp->pcpd[0][i].chr);
				unsigned int chr2 = GET_TAG2_CHR(tpp->pcpd[0][i].chr);
				if (chr1 != 0 && chr2 != 0)
					// chimera
					printf("%lu\t%uc%u\t%.*s\t%.*s\t%.*s\t%.*s\n",dpp->ordinal,chr1,chr2,len1,dpp->tag1 + start1,len2,dpp->tag2 + start2,len1,dpp->qual1 + start1,len2,dpp->qual2 + start2);
				else
				{
					if (chr1 != 0)
						printf("%lu\t%uh%u\t%.*s\t%.*s\t%.*s\t%.*s\n",dpp->ordinal,chr1,tpp->pcpd[0][i].tag1pos,len1,dpp->tag1 + start1,len2,dpp->tag2 + start2,len1,dpp->qual1 + start1,len2,dpp->qual2 + start2);
					else
						printf("%lu\t%uh%u\t%.*s\t%.*s\t%.*s\t%.*s\n",dpp->ordinal,chr2,tpp->pcpd[0][i].tag2pos,len1,dpp->tag1 + start1,len2,dpp->tag2 + start2,len1,dpp->qual1 + start1,len2,dpp->qual2 + start2);
				}
			}
		}
	}
	else
	{
		if ((tpp->header.flags_7_0 & kTLBHflagPairChrPosBlock) == 0)
			printf("%lu\tU\t%s\t%s\t%s\t%s\n",dpp->ordinal,dpp->tag1,dpp->tag2,dpp->qual1,dpp->qual2);
		else
		{
			unsigned int chr1 = GET_TAG1_CHR(tpp->pcpd[0][i].chr);
			unsigned int chr2 = GET_TAG2_CHR(tpp->pcpd[0][i].chr);
			if (chr1 != 0 && chr2 != 0)
				// chimera
				printf("%lu\t%uc%u\t%s\t%s\t%s\t%s\n",dpp->ordinal,chr1,chr2,dpp->tag1,dpp->tag2,dpp->qual1,dpp->qual2);
			else
			{
				if (chr1 != 0)
					printf("%lu\t%uh%u\t%s\t%s\t%s\t%s\n",dpp->ordinal,chr1,tpp->pcpd[0][i].tag1pos,dpp->tag1,dpp->tag2,dpp->qual1,dpp->qual2);
				else
					printf("%lu\t%uh%u\t%s\t%s\t%s\t%s\n",dpp->ordinal,chr2,tpp->pcpd[0][i].tag2pos,dpp->tag1,dpp->tag2,dpp->qual1,dpp->qual2);
			}
		}
	}
}
//---------------------------------------------------------------

static void
print_pair_contig_pmNg(TLBDATA *tpp, unsigned int i, decodedPair_p_t dpp)
{
	int leftTagPos,rightTagPos;
	int delta = (tpp->ppd[0][i].tag2offset & k32bOFFSET_VALUE_MASK);
	if (tpp->ppd[0][i].tag2offset & k32bNEGATIVE_OFFSET_BIT)
		delta = -delta;

	if (delta >= 0)
	{
		leftTagPos = tpp->ppd[0][i].tag1pos;
		rightTagPos = tpp->ppd[0][i].tag1pos + delta;
	}
	else
	{
		rightTagPos = tpp->ppd[0][i].tag1pos;
		leftTagPos = tpp->ppd[0][i].tag1pos + delta;
	}

	if (!((leftTagPos-200 > gFromPos && leftTagPos < gToPos)
		|| (rightTagPos-200 > gFromPos && rightTagPos < gToPos)))
		return;
	printf("%lu\t%s\n",dpp->ordinal,dpp->reverseTAG1[0] ? "l" : "L");
	printf("%lu\t%s\n",dpp->ordinal,dpp->reverseTAG1[0] ? "r" : "R");
}
//---------------------------------------------------------------

static void
print_pair_contig_u(TLBDATA *tpp, unsigned int i, decodedPair_p_t dpp)
{
	if ((tpp->header.flags_7_0 & kTLBHflagPairChrPosBlock) == 0)
		return;
	unsigned int chr1 = GET_TAG1_CHR(tpp->pcpd[0][i].chr);
	unsigned int chr2 = GET_TAG2_CHR(tpp->pcpd[0][i].chr);
	if (chr1 != 0 && chr2 != 0)
	{
		// chimera
	}
	else
	{
		if (chr1 != 0)
		{
		}
		else
		{
		}
	}
	int leftTagPos,rightTagPos;
	int delta = (tpp->ppd[0][i].tag2offset & k32bOFFSET_VALUE_MASK);
	if (tpp->ppd[0][i].tag2offset & k32bNEGATIVE_OFFSET_BIT)
		delta = -delta;

	if (delta >= 0)
	{
		leftTagPos = tpp->ppd[0][i].tag1pos;
		rightTagPos = tpp->ppd[0][i].tag1pos + delta;
	}
	else
	{
		rightTagPos = tpp->ppd[0][i].tag1pos;
		leftTagPos = tpp->ppd[0][i].tag1pos + delta;
	}

	if (!((leftTagPos-200 > gFromPos && leftTagPos < gToPos)
		|| (rightTagPos-200 > gFromPos && rightTagPos < gToPos)))
		return;
	printf("%lu\t%s\n",dpp->ordinal,dpp->reverseTAG1[0] ? "lr" : "LR");
}
//---------------------------------------------------------------

static void
print_pair_tromer_p(TLBDATA *tpp, unsigned int i, decodedPair_p_t dpp)
{
	for (unsigned int j = 0; j < tpp->nbMatches; j++)
	{
		unsigned int start = tpp->ppd[j][i].tag1pos;
		unsigned int tem = tpp->ppd[j][i].tag1pos + dpp->delta[j];
		if (tem < start)
			start = tem;
		unsigned int end = tpp->ppd[j][i].tag1pos + dpp->taglen1 - 1;
		tem = tpp->ppd[j][i].tag1pos + dpp->delta[j] + dpp->taglen2 - 1;
		if (tem > end)
			end = tem;
		if (j == 0)
			printf("%lu\t%u\tp\t%u\t%u\t%s\t%s",dpp->ordinal,dpp->chr[j],start,end,dpp->tag1,dpp->tag2);
		else
			printf("\t%u\tp\t%u\t%u",dpp->chr[j],start,end);
	}
	fputc('\n',stdout);
}
//---------------------------------------------------------------

static void
print_pair_tromer_m(TLBDATA *tpp, unsigned int i, decodedPair_p_t dpp)
{
	for (unsigned int j = 0; j < tpp->nbMatches; j++)
	{
		unsigned int start = tpp->ppd[j][i].tag1pos;
		unsigned int tem = tpp->ppd[j][i].tag1pos + dpp->delta[j];
		if (tem < start)
			start = tem;
		unsigned int end = tpp->ppd[j][i].tag1pos + dpp->taglen1 - 1;
		tem = tpp->ppd[j][i].tag1pos + dpp->delta[j] + dpp->taglen2 - 1;
		if (tem > end)
			end = tem;
		if (j == 0)
			printf("%lu\t%u\tm\t%u\t%u\t%s\t%s",dpp->ordinal,dpp->chr[j],start,end,dpp->tag1,dpp->tag2);
		else
			printf("\t%u\tm\t%u\t%u",dpp->chr[j],start,end);
	}
	fputc('\n',stdout);
}
//---------------------------------------------------------------

static void
print_pair_tromer_mN(TLBDATA *tpp, unsigned int i, decodedPair_p_t dpp)
{
	for (unsigned int j = 0; j < tpp->nbMatches; j++)
	{
		unsigned int start = tpp->ppd[j][i].tag1pos;
		unsigned int tem = tpp->ppd[j][i].tag1pos + dpp->delta[j];
		if (tem < start)
			start = tem;
		unsigned int end = tpp->ppd[j][i].tag1pos + dpp->taglen1 - 1;
		tem = tpp->ppd[j][i].tag1pos + dpp->delta[j] + dpp->taglen2 - 1;
		if (tem > end)
			end = tem;
		if (j == 0)
			printf("%lu\t%u\tN\t%u\t%u\t%s\t%s",dpp->ordinal,dpp->chr[j],start,end,dpp->tag1,dpp->tag2);
		else
			printf("\t%u\tN\t%u\t%u",dpp->chr[j],start,end);
	}
	fputc('\n',stdout);
}
//---------------------------------------------------------------

static void
print_pair_tromer_g(TLBDATA *tpp, unsigned int i, decodedPair_p_t dpp)
{
	for (unsigned int j = 0; j < tpp->nbMatches; j++)
	{
		unsigned int start = tpp->ppd[j][i].tag1pos;
		unsigned int tem = tpp->ppd[j][i].tag1pos + dpp->delta[j];
		if (tem < start)
			start = tem;
		unsigned int end = tpp->ppd[j][i].tag1pos + dpp->taglen1 - 1;
		tem = tpp->ppd[j][i].tag1pos + dpp->delta[j] + dpp->taglen2 - 1;
		if (tem > end)
			end = tem;
		if (j == 0)
			printf("%lu\t%u\tg\t%u\t%u\t%s\t%s",dpp->ordinal,dpp->chr[j],start,end,dpp->tag1,dpp->tag2);
		else
			printf("\t%u\tg\t%u\t%u",dpp->chr[j],start,end);
	}
	fputc('\n',stdout);
}
//---------------------------------------------------------------

static void
print_pair_tromer_u(TLBDATA *tpp, unsigned int i, decodedPair_p_t dpp)
{
	if ((tpp->header.flags_7_0 & kTLBHflagPairChrPosBlock) == 0)
		printf("%lu\t%u\tU\t%u\t%u\t%s\t%s\n",dpp->ordinal,0,0,0,dpp->tag1,dpp->tag2);
	else
	{
		unsigned int chr1 = (tpp->pcpd[0][i].chr >> 16) & 0x7fff;
		unsigned int chr2 = tpp->pcpd[0][i].chr & 0x7fff;
		if (chr1 != 0 && chr2 != 0)
			// chimera
			printf("%lu\t%u\tc\t%u\t%u\t%s\t%s\n",dpp->ordinal,0,0,0,dpp->tag1,dpp->tag2);
		else
		{
			if (chr1 != 0)
				printf("%lu\t%u\tH\t%u\t%u\t%s\t%s\n",dpp->ordinal,chr1,tpp->pcpd[0][i].tag1pos,tpp->pcpd[0][i].tag1pos+dpp->taglen1-1,dpp->tag1,dpp->tag2);
			else
				printf("%lu\t%u\th\t%u\t%u\t%s\t%s\n",dpp->ordinal,chr2,tpp->pcpd[0][i].tag2pos,tpp->pcpd[0][i].tag2pos+dpp->taglen2-1,dpp->tag1,dpp->tag2);
		}
	}
}
//---------------------------------------------------------------

#ifndef __APPLE__
static GlobalData_t ReAlignSettings[2];
static void
print_pair_realign(TLBDATA *tpp, unsigned int i, decodedPair_p_t dpp)
{
	void GenerateTagString(const GlobalData_t * const restrict AlignData, unsigned char * const restrict tag,
	                       unsigned char * const restrict qual, unsigned char * const restrict tagSource, char * const restrict qualSource,
	                       int * const restrict taglen)
	{
		int r = 0;
		int a = 0;
		const int TagLength = *taglen;

		for (int l=0; l<AlignData->ProfileRange[0]; l++) {
			tagSource[r] += 'a' - 'A'; tag[a] = tagSource[r]; qual[a] = qualSource[r]; a++; r++;
		}
		const unsigned char * const restrict Cigar = AlignData->Cigar;
		if (Cigar[0] != kCIGARTERMINATOR) { // if not, normal decoding.
			int cigarpos = 0;
			int ToBeSubtractedOnNextState = 0;
			while (Cigar[cigarpos] != kCIGARTERMINATOR) {
				int pos   = Cigar[cigarpos] - ToBeSubtractedOnNextState;
				ToBeSubtractedOnNextState = 0;
				char code = Cigar[cigarpos+1];
				cigarpos += 2;
				switch(code) {
					case 'S':  while(pos-- >0) { tagSource[r] += 'a' - 'A'; tag[a] = tagSource[r]; qual[a] = qualSource[r]; a++; r++; } break;
					case 'M':  while(pos-- >0) { tag[a] = tagSource[r]; qual[a] = qualSource[r]; a++; r++; } break;
					case 'D':  while(pos-- >0) { tag[a] = '_'; qual[a] = '_'; a++; *taglen += 1;}; break;
					case 'I':  while(pos-- >0) { tag[a] = '*'; qual[a] = '*'; r++; *taglen -= 1;}; ToBeSubtractedOnNextState=1; a++; r++; break;
				}
				assert(a < kMaxReadLen * 2);
			}
		}
		else {
			while (r<AlignData->ProfileRange[1]) {
				tag[a] = tagSource[r]; qual[a] = qualSource[r]; a++; r++;
			}
		}

		while (r<TagLength) {
			tagSource[r] += 'a' - 'A'; tag[a] = tagSource[r]; qual[a] = qualSource[r]; a++; r++;
		}

		tag[a] = '\0';
		qual[a] = '\0';
	}

	unsigned char tag1[kMaxReadLen], tag2[kMaxReadLen];

	/* Avoid too far to speed things up a bit */
	int leftTagPos,rightTagPos;
	int delta = (tpp->ppd[0][i].tag2offset & k32bOFFSET_VALUE_MASK);
	if (tpp->ppd[0][i].tag2offset & k32bNEGATIVE_OFFSET_BIT)
		delta = -delta;

	if (dpp->delta[0] >= 0)
	{
		leftTagPos = tpp->ppd[0][i].tag1pos;
		rightTagPos = tpp->ppd[0][i].tag1pos + delta;
	}
	else
	{
		rightTagPos = tpp->ppd[0][i].tag1pos;
		leftTagPos = tpp->ppd[0][i].tag1pos + delta;
	}

	if (!((leftTagPos-200 > gFromPos && leftTagPos < gToPos)
		|| (rightTagPos-200 > gFromPos && rightTagPos < gToPos)))
		return;

	/* Copy to UPPERCASE the tags !!! */
	for (size_t l=0; l<dpp->taglen1; l++) {
		const unsigned char c = (unsigned char) dpp->tag1[l];
		tag1[l] = (c >= 'a') ? c-'a'+'A': c;
	}
	tag1[dpp->taglen1] = '\0';
	for (size_t l=0; l<dpp->taglen2; l++) {
		const unsigned char c = (unsigned char) dpp->tag2[l];
		tag2[l] = (c >= 'a') ? c-'a'+'A': c;
	}
	tag2[dpp->taglen2] = '\0';
#ifdef PRINTOUT_DEBUG
	printf("REALIGN:\n");
#endif
	/* Compute new alignments */
	ReAlignSettings[0].TagLength = dpp->taglen1;
	ReAlignSettings[0].Tag       = tag1;
	ReAlignSettings[0].Genome    = ((unsigned char*) &genome_tbl[dpp->genomepos[0]]) - prf.ExtraLeftGenomeShift;
	ReAlignSettings[0].revNmer   = dpp->reverseTAG1[0];
	int k;
	if (!dpp->reverseTAG1[0]) {
		k=dpp->taglen1-1;
		while (k >= 0) if (dpp->qual1[k] != '#') break; else --k;
		k = (k < 0) ? 0 : k;
	}
	else {
		k = 0;
		while( k<dpp->taglen1) if (dpp->qual1[k] != '#') break; else ++k;
		k = ( k > dpp->taglen1) ? dpp->taglen1 - 1 : k;
	}
	ReAlignSettings[0].SoftClipBoundary = 1+k; // matrix border of 1

	ReAlignSettings[1].TagLength = dpp->taglen2;
	ReAlignSettings[1].Tag       = tag2;
	ReAlignSettings[1].Genome    = ((unsigned char*) &genome_tbl[dpp->genomepos[0]]) + dpp->delta[0] - prf.ExtraLeftGenomeShift;
	ReAlignSettings[1].revNmer   = dpp->reverseTAG2[0];
	if (!dpp->reverseTAG2[0]) {
		k=dpp->taglen2-1;
		while (k >= 0) if (dpp->qual2[k] != '#') break; else --k;
	}
	else {
		k = 0;
		while( k<dpp->taglen2) if (dpp->qual2[k] != '#') break; else ++k;
	}
	ReAlignSettings[1].SoftClipBoundary = 1+k; // matrix border of 1

	int offset[2];

	for (int iTag=0; iTag<2; iTag++) {
		int ret = cpu_std_border.createMatrix(&ReAlignSettings[iTag]);
		ret = cpu_std_border.getStateSequence(&ReAlignSettings[iTag]);
#ifdef PRINTOUT_DEBUG
		printf("TAG %i (len=%zu) @ %i-%i/%i-%i, score=%i, mismatch=%u, sofclip mismatch=%u, RequiredMM=%u, softclip evicted MM=%u,"
		       "RequireCigar=%u, softclip evicted cigar=%u, sofclip boundary @ %i\n",
					1+iTag, ReAlignSettings[iTag].TagLength,
					ReAlignSettings[iTag].AlignmentRange[0], ReAlignSettings[iTag].AlignmentRange[1],
					ReAlignSettings[iTag].ProfileRange[0], ReAlignSettings[iTag].ProfileRange[1],
					ReAlignSettings[iTag].score,
					ReAlignSettings[iTag].MismatchCount,
					ReAlignSettings[iTag].SoftClipMismatchCount,
					ReAlignSettings[iTag].RequiredMMSize,
					ReAlignSettings[iTag].SoftClipEvictedMMSize,
					ReAlignSettings[iTag].RequiredCigSize,
					ReAlignSettings[iTag].SoftClipEvictedCigSize,
					ReAlignSettings[iTag].SoftClipBoundary);
		printf("STA: ");
		for ( int space=0; space<ReAlignSettings[iTag].ProfileRange[0]; space++) fputc(' ', stdout);
		printf("%s\n", ReAlignSettings[iTag].States);
		printf("TAG: %.*s\n",(int) ReAlignSettings[iTag].TagLength, ReAlignSettings[iTag].Tag);
		printf("GEN: ");
		for ( int space=0; space<ReAlignSettings[iTag].ProfileRange[0]; space++) fputc(' ', stdout);
		printf("%.*s\n", (int) (ReAlignSettings[iTag].AlignmentRange[1] - ReAlignSettings[iTag].AlignmentRange[0]), &(ReAlignSettings[iTag].Genome[ReAlignSettings[iTag].AlignmentRange[0]]));
#endif
		if (ret != 0) goto TOO_SMALL;
		ret = cpu_std_border.computeCigarAndMismatch(&ReAlignSettings[iTag]);
		if (ret != 0) goto TOO_SMALL;

		offset[iTag] = (ReAlignSettings[iTag].AlignmentRange[0] - prf.ExtraLeftGenomeShift - ReAlignSettings[iTag].ProfileRange[0]);

		continue;

		TOO_SMALL:
		ReAlignSettings[iTag].Cigar[0] = (unsigned char) ReAlignSettings[iTag].TagLength;
		ReAlignSettings[iTag].Cigar[1] = 'S';
		ReAlignSettings[iTag].Cigar[2] = kCIGARTERMINATOR;
		offset[iTag] = /*prf.ExtraLeftGenomeShift*/ 0;
	}

	unsigned char tag[kMaxReadLen * 2]; // hold the tag revised with indels for display on genome
	unsigned char qual[kMaxReadLen * 2]; // hold the tag quality revised with indels for display on genome
#ifdef PRINTOUT_DEBUG
	printf("CIG 1: REV %i : ", ReAlignSettings[0].revNmer); printCigar(&ReAlignSettings[0]);
	printf("CIG 2: REV %i : ", ReAlignSettings[1].revNmer); printCigar(&ReAlignSettings[1]);

// 	printf("\nADNI\n");
// 	print_pair_ADNIview_g(tpp, i, dpp);
#endif
	int orientation = 0;

	if ( (ReAlignSettings[0].Genome + offset[0]) <= (ReAlignSettings[1].Genome + offset[1]) ){
		GenerateTagString(&ReAlignSettings[0], tag, qual, tag1, dpp->qual1, &(dpp->taglen1));
		printf("%.*s\t%s\t%s\t%s\t%s\t%u\t%d\t%d\t%d\t%lu\n",
		       dpp->hhlen1 - 1, dpp->hdr1 + 1, tag, qual, tag1, dpp->qual1, tpp->ppd[0][i].tag1pos-1+offset[0],
		       kIsPairedEnd, dpp->taglen1,orientation,dpp->ordinal);
		GenerateTagString(&ReAlignSettings[1], tag, qual, tag2, dpp->qual2, &(dpp->taglen2));
		printf("%.*s\t%s\t%s\t%s\t%s\t%u\t%d\t%d\t%d\t%lu\n",
		       dpp->hhlen2 - 1, dpp->hdr2 + 1, tag, qual, tag2, dpp->qual2, tpp->ppd[0][i].tag1pos-1+delta+offset[1],
		       kIsPairedEnd, dpp->taglen2,orientation,dpp->ordinal);
	}
	else {
		GenerateTagString(&ReAlignSettings[1], tag, qual, tag2, dpp->qual2, &(dpp->taglen2));
		printf("%.*s\t%s\t%s\t%s\t%s\t%u\t%d\t%d\t%d\t%lu\n",
		       dpp->hhlen2 - 1, dpp->hdr2 + 1, tag, qual, tag2, dpp->qual2, tpp->ppd[0][i].tag1pos-1+delta+offset[1],
		       kIsPairedEnd, dpp->taglen2,orientation,dpp->ordinal);
		GenerateTagString(&ReAlignSettings[0], tag, qual, tag1, dpp->qual1, &(dpp->taglen1));
		printf("%.*s\t%s\t%s\t%s\t%s\t%u\t%d\t%d\t%d\t%lu\n",
		       dpp->hhlen1 - 1, dpp->hdr1 + 1, tag, qual, tag1, dpp->qual1, tpp->ppd[0][i].tag1pos-1+offset[0],
		       kIsPairedEnd, dpp->taglen1,orientation,dpp->ordinal);
	}
}
//---------------------------------------------------------------
#endif

static void
output_p(TLBDATA *tpp, int bamStyle, int lowercase, void (*print_single) (TLBDATA *, unsigned int, decodedSingle_p_t), void (*print_pair) (TLBDATA *, unsigned int, decodedPair_p_t))
{
	if ((tpp->header.flags_7_0 & kTLBHflagPairedReads) == 0)
	{
		decodedSingle_t ds;
		int i;
		// decompress and print
		for (i=0; i<tpp->cnt; i++)
		{
			decodeSingle(tpp,i,&ds,bamStyle);
			print_single(tpp,i,&ds);
		}
	}
	else
	{
		decodedPair_t dp;
		int i;
		// decompress and print
		for (i=0; i<tpp->cnt; i++)
		{
			decodePair(tpp,i,&dp,bamStyle);
			print_pair(tpp,i,&dp);
		}
	}
}
//---------------------------------------------------------------

static void
output_m(TLBDATA *tpp, int bamStyle, int lowercase, void (*print_single) (TLBDATA *, unsigned int, decodedSingle_p_t), void (*print_pair) (TLBDATA *, unsigned int, decodedPair_p_t))
{
	if ((tpp->header.flags_7_0 & kTLBHflagPairedReads) == 0)
	{
		decodedSingle_t ds;
		int i;
		int mmpos[kMAXnbMatches] = { 0 };
		// decompress and print
		for (i=0; i<tpp->cnt; i++)
		{
			decodeSingle(tpp,i,&ds,bamStyle);
			mmpos[0] = applyDiffSingle(&ds,tpp->diff[0],mmpos[0],0);
			// applyDiff updates mmposStart for 0; need to do same for the other matches
			for (unsigned int j = 1; j < tpp->nbMatches; j++)
			{
				ds.mmposStart[j] = mmpos[j];
				while (tpp->diff[j][mmpos[j]] != kMMTERMINATOR)
					mmpos[j] += 1;
				mmpos[j] += 1;
			}
			print_single(tpp,i,&ds);
		}
	}
	else
	{
		decodedPair_t dp;
		int i;
		int mmpos[kMAXnbMatches] = { 0 };
		// decompress and print
		for (i=0; i<tpp->cnt; i++)
		{
			decodePair(tpp,i,&dp,bamStyle);
			mmpos[0] = applyDiffPair(&dp,tpp->diff[0],mmpos[0],0);
			// applyDiff updates mmposStart for 0; need to do same for the other matches
			for (unsigned int j = 1; j < tpp->nbMatches; j++)
			{
				dp.mmposStart[j] = mmpos[j];
				while (tpp->diff[j][mmpos[j]] != kMMTERMINATOR)
					mmpos[j] += 1;
				mmpos[j] += 1;
				while (tpp->diff[j][mmpos[j]] != kMMTERMINATOR)
					mmpos[j] += 1;
				mmpos[j] += 1;
			}
			print_pair(tpp,i,&dp);
		}
	}
}
//---------------------------------------------------------------

static void
output_mN(TLBDATA *tpp, int bamStyle, int lowercase, void (*print_single) (TLBDATA *, unsigned int, decodedSingle_p_t), void (*print_pair) (TLBDATA *, unsigned int, decodedPair_p_t))
{
	if ((tpp->header.flags_7_0 & kTLBHflagPairedReads) == 0)
	{
		decodedSingle_t ds;
		int i;
		int mmpos[kMAXnbMatches] = { 0 };
		// decompress and print
		for (i=0; i<tpp->cnt; i++)
		{
			decodeSingle(tpp,i,&ds,bamStyle);
			mmpos[0] = applyNdiffSingle(&ds,tpp->diff[0],mmpos[0]);
			// applyDiff updates mmposStart for 0; need to do same for the other matches
			for (unsigned int j = 1; j < tpp->nbMatches; j++)
			{
				ds.mmposStart[j] = mmpos[j];
				while (tpp->diff[j][mmpos[j]] != kMMTERMINATOR)
					mmpos[j] += 1;
				mmpos[j] += 1;
			}
			print_single(tpp,i,&ds);
		}
	}
	else
	{
		decodedPair_t dp;
		int i;
		int mmpos[kMAXnbMatches] = { 0 };
		// decompress and print
		for (i=0; i<tpp->cnt; i++)
		{
			decodePair(tpp,i,&dp,bamStyle);
			mmpos[0] = applyNdiffPair(&dp,tpp->diff[0],mmpos[0]);
			// applyDiff updates mmposStart for 0; need to do same for the other matches
			for (unsigned int j = 1; j < tpp->nbMatches; j++)
			{
				dp.mmposStart[j] = mmpos[j];
				while (tpp->diff[j][mmpos[j]] != kMMTERMINATOR)
					mmpos[j] += 1;
				mmpos[j] += 1;
				while (tpp->diff[j][mmpos[j]] != kMMTERMINATOR)
					mmpos[j] += 1;
				mmpos[j] += 1;
			}
			print_pair(tpp,i,&dp);
		}
	}
} /* output_fastq_mN */
//---------------------------------------------------------------

static void
output_g(TLBDATA *tpp, int bamStyle, int lowercase, void (*print_single) (TLBDATA *, unsigned int, decodedSingle_p_t), void (*print_pair) (TLBDATA *, unsigned int, decodedPair_p_t))
{
	if ((tpp->header.flags_7_0 & kTLBHflagPairedReads) == 0)
	{
		decodedSingle_t ds;
		int i;
		int mmpos[kMAXnbMatches] = { 0 };
		int cigarpos[kMAXnbMatches] = { 0 };
		// decompress and print
		for (i=0; i<tpp->cnt; i++)
		{
			decodeSingle(tpp,i,&ds,bamStyle);
			cigarpos[0] = applyCigarSingle(&ds,tpp->cigar[0],cigarpos[0]);
			// applyCigar updates cigarposStart for 0; need to do same for the other matches
			for (unsigned int j = 1; j < tpp->nbMatches; j++)
			{
				ds.cigarposStart[j] = cigarpos[j];
				while (tpp->cigar[j][cigarpos[j]] != kCIGARTERMINATOR)
					cigarpos[j] += 1;
				cigarpos[j] += 1;
			}
			mmpos[0] = applyDiffSingle(&ds,tpp->diff[0],mmpos[0],lowercase);
			// applyDiff updates mmposStart for 0; need to do same for the other matches
			for (unsigned int j = 1; j < tpp->nbMatches; j++)
			{
				ds.mmposStart[j] = mmpos[j];
				while (tpp->diff[j][mmpos[j]] != kMMTERMINATOR)
					mmpos[j] += 1;
				mmpos[j] += 1;
			}
			print_single(tpp,i,&ds);
		}
	}
	else
	{
		decodedPair_t dp;
		int i;
		int mmpos[kMAXnbMatches] = { 0 };
		int cigarpos[kMAXnbMatches] = { 0 };
		// decompress and print
		for (i=0; i<tpp->cnt; i++)
		{
			decodePair(tpp,i,&dp,bamStyle);
			cigarpos[0] = applyCigarPair(&dp,tpp->cigar[0],cigarpos[0]);
			// applyCigar updates cigarposStart for 0; need to do same for the other matches
			for (unsigned int j = 1; j < tpp->nbMatches; j++)
			{
				dp.cigarposStart[j] = cigarpos[j];
				while (tpp->cigar[j][cigarpos[j]] != kCIGARTERMINATOR)
					cigarpos[j] += 1;
				cigarpos[j] += 1;
				while (tpp->cigar[j][cigarpos[j]] != kCIGARTERMINATOR)
					cigarpos[j] += 1;
				cigarpos[j] += 1;
			}
			mmpos[0] = applyDiffPair(&dp,tpp->diff[0],mmpos[0],lowercase);
			// applyDiff updates mmposStart for 0; need to do same for the other matches
			for (unsigned int j = 1; j < tpp->nbMatches; j++)
			{
				dp.mmposStart[j] = mmpos[j];
				while (tpp->diff[j][mmpos[j]] != kMMTERMINATOR)
					mmpos[j] += 1;
				mmpos[j] += 1;
				while (tpp->diff[j][mmpos[j]] != kMMTERMINATOR)
					mmpos[j] += 1;
				mmpos[j] += 1;
			}
			print_pair(tpp,i,&dp);
		}
	}
}
//---------------------------------------------------------------

static void
output_u(TLBDATA *tpp, int bamStyle, int lowercase, void (*print_single) (TLBDATA *, unsigned int, decodedSingle_p_t), void (*print_pair) (TLBDATA *, unsigned int, decodedPair_p_t))
{
	if ((tpp->header.flags_7_0 & kTLBHflagPairedReads) == 0)
	{
		decodedSingle_t ds;
		int i;
		int mmpos = 0;
		// decompress and print
		for (i=0; i<tpp->cnt; i++)
		{
			mmpos = decodeSingleUnmapped(tpp,i,&ds,mmpos,bamStyle);
			print_single(tpp,i,&ds);
		}
	}
	else
	{
		decodedPair_t dp;
		int i;
		int mmpos = 0;
		// decompress and print
		for (i=0; i<tpp->cnt; i++)
		{
			mmpos = decodePairUnmapped(tpp,i,&dp,mmpos,bamStyle);
			print_pair(tpp,i,&dp);
		}
	}
}
//---------------------------------------------------------------

static int
decompress(char *fn, OPTIONS *opt)
{
	int fd = open(fn, O_RDONLY);
	TLBHEADER th;
	if (fd < 0)
	{
		fprintf(stderr,"couldn't open file %s\n",fn);
		perror("open");
		return 1;
	}
	// FIXME - not all that clean...  no better way yet
	gFilterSAM = opt->filterSAM;
	TLBDATA td;
	allocBlock(&td, true);
	while (1) {
		ssize_t res = read(fd,&th,sizeof(TLBHEADER));
		if (res == 0)
			break;
		if (res != sizeof(TLBHEADER))
		{
			fprintf(stderr,"Wrong size %ld while reading %s as %d\n",res,fn,fd);
			perror("read");
			close(fd);
			return 1;
		}
		unsigned char hcs = simple8bitCS((unsigned char *) &th, sizeof(TLBHEADER));
		if (hcs != 0)
		{
			fprintf(stderr,"Checksum mismatch: %02x %02x\n",th.headerCS,hcs);
			close(fd);
			return 1;
		}
		// do we want this block ?
		void (*outFunc) (TLBDATA *, int, int, void (*) (TLBDATA *, unsigned int, decodedSingle_p_t), void (*) (TLBDATA *, unsigned int, decodedPair_p_t)) = NULL;
		void (*printSingle) (TLBDATA *, unsigned int, decodedSingle_p_t);
		void (*printPair) (TLBDATA *, unsigned int, decodedPair_p_t);
		if ((th.flags_15_8 & kTLBHflagUnmappedBlock) != 0
				|| (th.flags_15_8 & kTLBHflagChimeraBlock) != 0
				|| (th.flags_15_8 & kTLBHflagHalfmapBlock) != 0)
		{
			if ((th.flags_15_8 & kTLBHflagUnmappedBlock) != 0 && opt->decompressUnmapped == 1)
			{
				outFunc = output_u;
				printSingle = PrintFuncs[outSelect].print_single_u;
				printPair = PrintFuncs[outSelect].print_pair_u;
			}
			else if ((th.flags_15_8 & kTLBHflagChimeraBlock) != 0 && opt->decompressChimeric == 1)
			{
				outFunc = output_u;
				printSingle = PrintFuncs[outSelect].print_single_u;
				printPair = PrintFuncs[outSelect].print_pair_u;
			}
			else if ((th.flags_15_8 & kTLBHflagHalfmapBlock) != 0 && opt->decompressHalfmapped == 1)
			{
				outFunc = output_u;
				printSingle = PrintFuncs[outSelect].print_single_u;
				printPair = PrintFuncs[outSelect].print_pair_u;
			}
		}
		else if (th.chr >= opt->fromchr && th.chr <= opt->tochr && !(opt->toPos < th.minPos || opt->fromPos > th.maxPos))
		{

			//fprintf(stderr,"outSelect = %d  block: chr %d:%d %d\n",outSelect,th.chr,th.minPos,th.maxPos);
			if (PrintFuncs[outSelect].rangeSelect && th.minPos > opt->toPos)
			{
				// FIXME: only valid for sorted files.
				close(fd);
				return 0;
			}

			if ((th.flags_7_0 & kTLBHflagPairPosBlock) != 0)
			{
				if ((th.flags_7_0 & 0xf0) == 0 && opt->decompressPerfectMatches == 1)
				{
					outFunc = output_p;
					printSingle = PrintFuncs[outSelect].print_single_p;
					printPair = PrintFuncs[outSelect].print_pair_p;
				}
				if ((th.flags_7_0 & 0xf0) == kTLBHflagNMismatchBlock && opt->decompressMismatchesN == 1)
				{
					outFunc = output_mN;
					printSingle = PrintFuncs[outSelect].print_single_mN;
					printPair = PrintFuncs[outSelect].print_pair_mN;
				}
				if ((th.flags_7_0 & 0xf0) == kTLBHflagMismatchBlock && opt->decompressMismatches == 1)
				{
					outFunc = output_m;
					printSingle = PrintFuncs[outSelect].print_single_m;
					printPair = PrintFuncs[outSelect].print_pair_m;
				}
				if ((th.flags_7_0 & 0xf0) == (kTLBHflagMismatchBlock | kTLBHflagCigarBlock) && opt->decompressAligned == 1)
				{
					outFunc = output_g;
					printSingle = PrintFuncs[outSelect].print_single_g;
					printPair = PrintFuncs[outSelect].print_pair_g;
				}
			}
		}
		if (outFunc != NULL)
		{
			// we do...
			if (readDecompressBlock(fd, &th, &td, 0) != 0)
			{
				close(fd);
				return 1;
			}
			outFunc(&td,PrintFuncs[outSelect].bamStyle,PrintFuncs[outSelect].lowercase,printSingle,printPair);
		}
		else
		{
			// nope, read and discard (so that we can stream...)
			size_t size = (size_t) th.blockLength - sizeof(TLBHEADER);
			assert(size <= kMaxBlockLength);
			void *data = td.diskBuffer;
			ssize_t res = read(fd,data,size);
			while (size > 0 && res != size)
			{
				if (res < 0)
				{
					fprintf(stderr,"%lu %lu\n",res,size);
					perror("skip read:");
					exit(1);
				}
				size -= res;
				res = read(fd,data,size);
			}
		}
	}

	// clean up
	close(fd);
	freeBlock(&td);

	return 0;

} // decompress
//---------------------------------------------------------------

static void
print_single_raw(TLBDATA *tpp, unsigned int i, decodedSingle_p_t dpp)
{
	fputs("Caramba olé, ca manque!!\n", stdout);
}
//---------------------------------------------------------------

static void
print_pair_raw(TLBDATA *tpp, unsigned int i, decodedPair_p_t dpp)
{
	const unsigned long ordinal = dpp->ordinal;
	printf("%1$lu00 Ordinal %1$lu", ordinal);
	if (tpp->header.flags_7_0 & kTLBHflagPairedReads)
		printf(" flagPairedReads");
	if (tpp->header.flags_7_0 & kTLBHflagFixedLength)
		printf(" FixedLen %u", tpp->lengthInfo);
	else
		printf(" VariableLen");
	if (tpp->header.flags_15_8 & kTLBHflagUnmappedBlock)
		printf(" flagUnmapped");
	if (tpp->header.flags_15_8 & kTLBHflagChimeraBlock)
		printf(" flagChimera");
	if (tpp->header.flags_15_8 & kTLBHflagHalfmapBlock)
		printf(" flagHalfmap");
	if (tpp->header.flags_7_0 & kTLBHflagPairChrPosBlock)
		printf(" flagPairChrPos");
	if (tpp->header.flags_7_0 & kTLBHflagPairPosBlock)
		printf(" flagPairPos");
	if (tpp->header.flags_7_0 & kTLBHflagNMismatchBlock)
		printf(" flagNMismatch");
	if (tpp->header.flags_7_0 & kTLBHflagMismatchBlock)
		printf(" flagMismatch");
	if (tpp->header.flags_7_0 & kTLBHflagCigarBlock)
		printf(" flagCigar");
	printf(" nbMatches %u headerChr %08X\n", tpp->nbMatches, tpp->header.chr);
	printf("%1$lu01 Head1 %2$u %3$.*2$s (%4$u %5$.*4$s)\n", ordinal, dpp->hlen1, dpp->hdr1, dpp->hhlen1 - 1, dpp->hdr1 + 1);
	printf("%1$lu02 Read1 %2$u %3$.*2$s\n", ordinal, dpp->taglen1, dpp->tag1);
	printf("%1$lu03 Qual1 %2$u %3$.*2$s\n", ordinal, dpp->taglen1, dpp->qual1);
	printf("%1$lu04 Head2 %2$u %3$.*2$s (%4$u %5$.*4$s)\n", ordinal, dpp->hlen2, dpp->hdr2, dpp->hhlen2 - 1, dpp->hdr2 + 1);
	printf("%1$lu05 Read2 %2$u %3$.*2$s\n", ordinal, dpp->taglen2, dpp->tag2);
	printf("%1$lu06 Qual2 %2$u %3$.*2$s\n", ordinal, dpp->taglen2, dpp->qual2);
	unsigned int lno = 7;
	for (unsigned int j = 0; j < tpp->nbMatches; j++)
	{
		printf("%1$lu%2$02u SAM%3$u flag1 %4$u flag2 %5$u ilen %6$d chr %7$u\n", ordinal, lno++, j, dpp->flag1[j], dpp->flag2[j], dpp->ilen[j], dpp->chr[j]);
		if (tpp->header.flags_7_0 & kTLBHflagPairPosBlock)
		{
			printf("%1$lu%2$02u PP%3$u tag1pos %4$u tag2offset %5$08X delta %6$d rev1 %7$u rev2 %8$u chr2 %9$u\n",
				ordinal, lno++, j, tpp->ppd[j][i].tag1pos, tpp->ppd[j][i].tag2offset,
				dpp->delta[j], dpp->reverseTAG1[j], dpp->reverseTAG2[j], (tpp->ppd[j][i].tag2offset & k32b2nd_MATCH_CHR_MASK) >> 19);
		}
		if (tpp->header.flags_7_0 & kTLBHflagPairChrPosBlock)
		{
			printf("%1$lu%2$02u PCP%3$u tag1pos %4$u tag2pos %5$u chr %6$08X chr1 %7$u rev1 %8$u chr2 %9$u rev2 %10$u\n",
				ordinal, lno++, j, tpp->pcpd[j][i].tag1pos, tpp->pcpd[j][i].tag2pos, tpp->pcpd[j][i].chr,
				GET_TAG1_CHR(tpp->pcpd[j][i].chr), dpp->reverseTAG1[j], GET_TAG2_CHR(tpp->pcpd[j][i].chr), dpp->reverseTAG2[j]);
		}
		if (tpp->header.flags_7_0 & kTLBHflagNMismatchBlock)
		{
			printf("%1$lu%2$02u MMN%3$u READ1 ", ordinal, lno++, j);
			const unsigned char *mm = tpp->diff[j];
			mm += dpp->mmposStart[j];
			while (*mm != kMMTERMINATOR) {
				const unsigned int pos  = (unsigned int) mm[0];
				printf("%u ", pos);
				mm += 1;
			}
			printf("END ; READ2 ");
			mm++;
			while (*mm != kMMTERMINATOR) {
				const unsigned int pos  = (unsigned int) mm[0];
				printf("%u ", pos);
				mm += 1;
			}
			printf("END\n");
		}
		if (tpp->header.flags_7_0 & kTLBHflagMismatchBlock)
		{
			printf("%1$lu%2$02u MM%3$u READ1 ", ordinal, lno++, j);
			const unsigned char *mm = tpp->diff[j];
			mm += dpp->mmposStart[j];
			while (*mm != kMMTERMINATOR) {
				if (*mm == kMMSOFTCLIP) {
					mm++;
					printf("SOFTCLIP ");
					while (*mm != kMMTERMINATOR) fputc((int) *mm++, stdout);
					fputc((int) ' ', stdout);
					break;
				}
				const unsigned int pos  = (unsigned int) mm[0];
				const int letter = (int) mm[1];
				printf("%u %c ", pos, letter);
				mm += 2;
			}
			printf("END ; READ2 ");
			mm++;
			while (*mm != kMMTERMINATOR) {
				if (*mm == kMMSOFTCLIP) {
					mm++;
					printf("SOFTCLIP ");
					while (*mm != kMMTERMINATOR) fputc((int) *mm++, stdout);
					fputc((int) ' ', stdout);
					break;
				}
				const unsigned int pos  = (unsigned int) mm[0];
				const int letter = (int) mm[1];
				printf("%u %c ", pos, letter);
				mm += 2;
			}
			printf("END\n");
		}
		if (tpp->header.flags_7_0 & kTLBHflagCigarBlock)
		{
			printf("%1$lu%2$02u CIG%3$u READ1 ", ordinal, lno++, j);
			const unsigned char * cigar = tpp->cigar[j];
			cigar += dpp->cigarposStart[j];
			while (*cigar != kCIGARTERMINATOR) {
				const unsigned int len  = (unsigned int) cigar[0];
				const int type = (int) cigar[1];
				printf("%u%c ", len, type);
				cigar += 2;
			}
			printf("END ; READ2 ");
			cigar++;
			while (*cigar != kCIGARTERMINATOR) {
				const unsigned int len  = (unsigned int) cigar[0];
				const int type = (int) cigar[1];
				printf("%u%c ", len, type);
				cigar += 2;
			}
			printf("END\n");
		}
	}
}
//---------------------------------------------------------------

static void
loadTrim(const char *fn)
{
#ifdef BIN_OUTPUT
	trim = malloc(kMaxTrimPairs * sizeof(trimData_t));
	if (trim == NULL)
	{
		perror("loadTrim : calloc");
		exit(1);
	}
	int fd = open(fn,O_RDONLY);
	if (fd == -1)
	{
		perror("open");
		exit(1);
	}
	void *ptr = trim;
	size_t cnt = kMaxTrimPairs * sizeof(trimData_t);
	ssize_t res = read(fd,ptr,cnt);
	while (res != cnt)
	{
		if (res <= 0)
		{
			perror("read");
			fprintf(stderr,"problem reading %s : %ld\n",fn,res);
			exit(1);
		}
		ptr += res;
		cnt -= res;
		res = read(fd,ptr,cnt);
	}
	close(fd);
#else
	trim = calloc(kMaxTrimPairs,sizeof(trimData_t));
	if (trim == NULL)
	{
		perror("loadTrim : calloc");
		exit(1);
	}
	FILE *f = fopen(fn,"r");
	if (f == NULL)
	{
		fprintf(stderr,"Can't open trim file %s\n",fn);
		exit(1);
	}
	unsigned long ord;
	unsigned int lenL1, utL1, lenR1, utR1, lenL2, utL2, lenR2, utR2;
	while (fscanf(f,"%lu %u %u %u %u %u %u %u %u",&ord,&lenL1,&utL1,&lenR1,&utR1,&lenL2,&utL2,&lenR2,&utR2) == 9)
	{
		assert(ord < kMaxTrimPairs);
		trim[ord].len[0] = lenL1;
		trim[ord].cnt[0] = utL1;
		trim[ord].len[1] = lenR1;
		trim[ord].cnt[1] = utR1;
		trim[ord].len[2] = lenL2;
		trim[ord].cnt[2] = utL2;
		trim[ord].len[3] = lenR2;
		trim[ord].cnt[3] = utR2;
	}
	fclose(f);
#endif
}
//---------------------------------------------------------------

int
main (int argc, char **argv)
{
	int c;
	int inputIsDirectory;
	char rsltfile[512];
	char inputdir[512];
	char genomefile[512];
#ifndef __APPLE__
	const char * restrict AlignerScoreConfFile = NULL;
#endif
	int err;
	int chr = -1;
	OPTIONS options;
	memset(&options,0,sizeof(OPTIONS));

	options.fromPos = 0;
	options.toPos = INT_MAX;  // this one is to large or negative ?  and proidices no output... 2147483647;

	/* --------- process arguments ----------------------*/

	inputIsDirectory = 0;
	rsltfile[0] = 0;
	strcpy(genomefile,"/tmp/nodelete/nguex/data/hg19.bin");
#ifdef __APPLE__
	while ((c = getopt (argc, argv, "r:i:g:C:P:o:v:acdhmnpt:u")) != -1)
#else
	while ((c = getopt (argc, argv, "r:i:g:C:P:o:v:acdhmnpt:uR:D:f")) != -1)
#endif
	switch (c)
	{
		case 'r':
			strcpy(rsltfile,optarg);
			break;

		case 'i':
			strcpy(inputdir,optarg);
			inputIsDirectory = 1;
			break;

		case 'C':
			sscanf(optarg,"%d",&chr);
			break;

		case 'P':
			sscanf(optarg,"%d..%d",&options.fromPos,&options.toPos);
			break;

		case 'g':
			strcpy(genomefile,optarg);
			break;

		case 'p':
			options.decompressPerfectMatches = 1;
			break;

		case 'a':
			options.decompressAligned = 1;
			break;

		case 'm':
			options.decompressMismatches = 1;
			break;

		case 'n':
			options.decompressMismatchesN = 1;
			break;

		case 'u':
			options.decompressUnmapped = 1;
			break;

		case 'h':
			options.decompressHalfmapped = 1;
			break;

		case 'c':
			options.decompressChimeric = 1;
			break;

		case 'd':
			debug = 1;
			break;

		case 't':
			loadTrim(optarg);
			break;

		case 'v':
			sscanf(optarg,"%d",&verbose);
			break;

		case 'o':
			{
				int i;
				for (i = 0; i < kOUTTYPESNB; i++) {
					if (strcasecmp(optarg, PrintFuncs[i].name) == 0)
						break;
				}
				if (i < kOUTTYPESNB)
					outSelect = i;
				else {
					fprintf(stderr, "Output type %s is not known\n", optarg);
					fprintf(stderr, "  known types:");
					for (i = 0; i < kOUTTYPESNB; i++)
						fprintf(stderr, " %s", PrintFuncs[i].name);
					fprintf(stderr, "\n");
					return 1;
				}
			}
			break;

#ifndef __APPLE__
		case 'R':
			AlignerScoreConfFile = optarg;
			if (loadAlignerScores(AlignerScoreConfFile) < 0) {
				fputs("Aligner configuration file failed to load\n", stderr);
				return 1;
			}
			break;

		case 'D':
			dumpAlignerScores(optarg);
			return 0;
			break;
#endif

		case 'f':
			options.filterSAM = 1;
			break;

		case '?':
			break;
	}


	if (inputIsDirectory && PrintFuncs[outSelect].rangeSelect)
	{
		char kind;

		if (options.decompressPerfectMatches)
			kind = 'p';
		else if (options.decompressMismatchesN)
			kind = 'N';
		else if (options.decompressMismatches)
			kind = 'm';
		else if (options.decompressAligned)
			kind = 'g';

		if ((options.decompressPerfectMatches + options.decompressMismatchesN + options.decompressMismatches + options.decompressAligned) != 1)
		{
			printf("only options -p -n -m -a are valid for ADNIview. Choose only one of them at the same time\n");
			return(1);
		}
		sprintf(rsltfile,"%s/chr%d%c_%c.gtl",inputdir,chr,'a' + options.fromPos / kCountsPerChrPart,kind);
	}


	if (rsltfile[0] == 0)
	{
		printf("usage:\n\n");
		printf("decompress [-g genomefile]  [-C chromosome] [-r resultfile | -i resultdir] [-p][-m][-a][-u][-h][-c] [-o outputType] [-v level]\n\n");
		printf("           -g genomefile         : full path of the genome encoded by encode_genome\n");
		printf("                                   default: '%s' \n",genomefile);
		printf("           -C chromosome         : 1..25 to restrict to a specific chromosme, default means all\n");
		printf("           -P fromPos..ToPos     : to restrict to a specific chromosme region, default means all\n");
		printf("           -r resultfile         : name of the result file (produced by match)\n");
		printf("           -i resultdir          : name of the directory containing sorted indexed results\n");
		printf("                                   note: only valid along with -o ADNIview\n");
		printf("           -p                    : decompress tags with perfect match\n");
		printf("           -m                    : decompress tags with mismatch\n");
		printf("           -n                    : decompress tags whose mismatches are only N\n");
		printf("           -a                    : decompress tags with alignments\n");
		printf("           -u                    : decompress unmapped tags\n");
		printf("           -h                    : decompress half mapped tags\n");
		printf("           -c                    : decompress chimeric tags\n");
		printf("           -o outputType         : output type [%s]\n", PrintFuncs[outSelect].name);
		printf("                                   note: for ADNIview, GTL file needs to be sorted !  only one kind of tag per call !!\n");
		printf("           -f                    : filter SAM output - do not output pairs that extend beyond reference template\n");
		printf("           -d                    : debug mode, prints line number of original fastQ input\n");
		printf("                                   tag1 in stdo tag2 in stderr\n");
		printf("                                   ./GTLdecompress -r /tmp/test_all.gtl -p -n -m -a -u -h -c -d >r1.fq 2>r2.fq\n");
#ifndef __APPLE__
		printf("           -R alignconffile      : aligner configuration file\n");
		printf("           -D alignconffile      : dump alignment score JSON file\n");
#endif
		printf("           -v level              : verbose level\n");
		printf("\n");
		printf(GTL_VERSION_FULL_STRING "\n");
		printf("(c) Nicolas Guex & Christian Iseli & Cie 2014-2019\n");
		return(1);
	}


	char *configfn = strdup(genomefile);
	unsigned int len = strlen(configfn);
	memcpy(configfn + len - 3,"cfg",3);

	err = loadGenome(genomefile);
	if (err)
		goto bail;

	initializeVirtualChromosomes(configfn,virtchr);
	initRevCompTable(gIUPACrevcomp);

	if (chr == -1)
	{
		options.fromchr = 1;
		options.tochr = chrMax;
	}
	else
		options.fromchr = options.tochr = chr;

#ifndef __APPLE__
	if (strcmp(PrintFuncs[outSelect].name, "realign") == 0) {
		memset(&ReAlignSettings[0], 0, 2*sizeof(GlobalData_t));
		ReAlignSettings[0].TagLength = kMaxReadLen;
		ReAlignSettings[0].GenomeLength = kGAPPED_ALIGN_GENOME_LENGTH + prf.ExtraLeftGenomeShift;

		if (cpu_std_border.allocStorage(&ReAlignSettings[0]) != 0) {
			fputs("Unable to allocate storage for realigner\n", stderr);
			goto bail;
		}

		ReAlignSettings[1].TagLength = kMaxReadLen;
		ReAlignSettings[1].GenomeLength = kGAPPED_ALIGN_GENOME_LENGTH + prf.ExtraLeftGenomeShift;
		if (cpu_std_border.allocStorage(&ReAlignSettings[1]) != 0) {
			fputs("Unable to allocate storage for realigner\n", stderr);
			goto bail;
		}
	}
#endif

	//fprintf(stderr, "chr %d %d %d\n", options.fromchr,options.fromPos,options.toPos);
	gFromPos = options.fromPos;
	gToPos = options.toPos;
	err = decompress(rsltfile,&options);

	err = 0;

	/* ------------------------ clean up and quit ------------------- */

bail:;
#ifndef __APPLE__
	if (strcmp(PrintFuncs[outSelect].name, "realign") == 0) {
		cpu_std_border.freeStorage(&ReAlignSettings[0]);
		cpu_std_border.freeStorage(&ReAlignSettings[1]);
	}
#endif

	if (genome_tbl != MAP_FAILED)
		munmap(genome_tbl,kGENOME_DATA_SIZE);

	if (verbose > 0)
	{
		fprintf(stderr,"pos\tMMcnt\n");
		for (unsigned int i = 0; i <= maxMMpos; i++)
			fprintf(stderr,"%u\t%u\n",i,statMMpos[i]);
	}
	return(err);

} /* main */
//---------------------------------------------------------------

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
