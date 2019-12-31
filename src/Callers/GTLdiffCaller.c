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

#include <sys/io.h>
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
#include "GTLcaller_common.h"

//=================================================================================================

#define kMAXSAMPLE 2
#define kMAXCALLS 128

//=================================================================================================

typedef struct _delVar {
	unsigned int len;
	unsigned int cnt[kMAXSAMPLE];
	unsigned char isAllele[kMAXSAMPLE];
} delVar_t;

typedef struct _insVar {
	char nt[kMAXINDELLEN];
	unsigned int cnt[kMAXSAMPLE];
	unsigned char isAllele[kMAXSAMPLE];
} insVar_t;

typedef struct _MVar {
	unsigned int cnt[kMAXSAMPLE];
	unsigned char isAllele[kMAXSAMPLE];
} MVar_t;

typedef struct _calledVar {
	unsigned int pos;
	unsigned int nbD;
	unsigned int nbI;
	unsigned int nbDal;
	unsigned int nbIal;
	unsigned int nbM;
	unsigned int maxDelLen;
	char delRef[kMAXINDELLEN];
	delVar_t del[kMAXCALLS];
	insVar_t ins[kMAXCALLS];
	MVar_t M[4];
	unsigned int refCnt[kMAXSAMPLE];
	unsigned char isRefAllele[kMAXSAMPLE];
	// special case for inserts, where ref is not covered by the reads with inserts
	unsigned int refICnt[kMAXSAMPLE];
	unsigned char isRefIAllele[kMAXSAMPLE];
	char refNt;
	char refINt;
} calledVar_t, *calledVar_p_t;

typedef struct _alCnt_t {
	unsigned int al;
	unsigned int cnt;
	unsigned int called;
} alCnt_t, *alCnt_p_t;

//=================================================================================================

//static struct timeval ts;
//static struct timeval te;
static int verbose = 0;
static int debug = 0;
static int minAllCov = 8;
static int minVarCov = 4;

//=================================================================================================

static int
is_significant_p(coverage_p_t cov)
{
	int ret = 0;
	if (cov->var != NULL)
	{
		unsigned int allCov = cov->cntM + cov->cntP;
		unsigned int varCov = 0;
		unsigned int j;
		for (j = 0; j < 4; j++)
		{
			allCov += cov->var->ACGTP[j] + cov->var->ACGTM[j];
			varCov += cov->var->ACGTP[j] + cov->var->ACGTM[j];
		}
		for (j = 0; j < kMAXNBINDELS; j++)
		{
			if (cov->var->del[j].cnt == 0)
				break;
			allCov += cov->var->del[j].cnt;
			varCov += cov->var->del[j].cnt;
		}
		for (j = 0; j < kMAXNBINDELS; j++)
		{
			if (cov->var->ins[j].cnt == 0)
				break;
			// no allCov in this case...
			varCov += cov->var->ins[j].cnt;
		}
		if (allCov >= minAllCov && varCov >= minVarCov)
			ret = 1;
	}
	return ret;
}
//---------------------------------------------------------------

static int
is_diff_p(coverage_p_t cov1, coverage_p_t cov2)
{
	unsigned int j;
	unsigned int tot1 = cov1->cntP + cov1->cntM;
	unsigned int tot2 = cov2->cntP + cov2->cntM;
	for (j = 0; j < 4; j++)
	{
		tot1 += cov1->var->ACGTP[j] + cov1->var->ACGTM[j];
		if (cov2->var != NULL)
		{
			tot2 += cov2->var->ACGTP[j] + cov2->var->ACGTM[j];
		}
	}
	for (j = 0; j < 4; j++)
	{
		if (((cov1->var->ACGTP[j] > 1 && cov1->var->ACGTM[j] > 1) || (tot1 >= minAllCov && ((cov1->var->ACGTP[j] + cov1->var->ACGTM[j]) * 100) / tot1 >= 33))
				&& cov1->var->ACGTP[j] + cov1->var->ACGTM[j] >= minVarCov
				&& ((cov1->var->ACGTP[j] + cov1->var->ACGTM[j]) * 100) / tot1 >= 8
				&& tot2 >= minAllCov
				&& (cov2->var == NULL || cov2->var->ACGTP[j] + cov2->var->ACGTM[j] == 0)
				&& (cov2->var == NULL || ((cov2->var->ACGTP[j] + cov2->var->ACGTM[j]) * 100) / tot2 <= 1))
			return 1;
	}
	return 0;
}
//---------------------------------------------------------------

static void
doCall(gtl_variant_p_t var, unsigned int cntP, unsigned int cntM, calledVar_p_t call, unsigned int sample, unsigned int genomepos, unsigned int i, OPTIONS *opt)
{
	if (sample == 0)
	{
		call->pos = i;
		call->refINt = genome_tbl[genomepos+i-1];
		call->refNt = genome_tbl[genomepos+i];
	}
	if (var == NULL)
	{
		call->refCnt[sample] = cntP + cntM;
		call->refICnt[sample] = call->refCnt[sample];
		if (call->refCnt[sample] >= opt->reportVCF)
			call->isRefAllele[sample] = 1;
		call->isRefIAllele[sample] = call->isRefAllele[sample];
		return;
	}
	unsigned int j;
	unsigned int total = 0;
	var->ACGTP[gNT2bits[genome_tbl[genomepos+i]]] += cntP;
	var->ACGTM[gNT2bits[genome_tbl[genomepos+i]]] += cntM;
	for (j = 0; j < 4; j++)
	{
		var->ACGTP[j] += var->ACGTM[j];
		total += var->ACGTP[j];
	}
	for (j = 0; j < kMAXNBINDELS; j++)
	{
		if (var->del[j].cnt == 0)
			break;
		total += var->del[j].cnt;
	}
	call->refCnt[sample] = var->ACGTP[gNT2bits[genome_tbl[genomepos+i]]];
	if (call->refCnt[sample] >= opt->reportVCF && call->refCnt[sample] >= total * opt->minorAllelePct / 100)
		call->isRefAllele[sample] = 1;
	call->refICnt[sample] = call->refCnt[sample];
	// report indels first, so that vcf output is ordered by pos
	for (j = 0; j < kMAXNBINDELS; j++)
	{
		if (var->ins[j].cnt == 0)
			break;
		unsigned int k = 0;
		// reads part of an insert do not count towards reference allele, even though the reference allele is "covered"...
		if (call->refICnt[sample] >= var->ins[j].cnt)
			call->refICnt[sample] -= var->ins[j].cnt;
		else
			call->refICnt[sample] = 0;
		while (k < call->nbI && strncmp(call->ins[k].nt,var->ins[j].nt,38) != 0)
			k += 1;
		if (k >= call->nbI)
		{
			if (k >= kMAXCALLS)
			{
				fprintf(stderr,"nbI >= kMAXCALLS\n");
				exit(1);
			}
			call->nbI += 1;
			strncpy(call->ins[k].nt,var->ins[j].nt,38);
		}
		call->ins[k].cnt[sample] = var->ins[j].cnt;
		if (var->ins[j].cnt >= opt->reportVCF && var->ins[j].cnt >= total * opt->minorAllelePct / 100)
		{
			call->ins[k].isAllele[sample] = 1;
			call->nbIal += 1;
		}
	}
	if (call->refICnt[sample] >= opt->reportVCF && call->refICnt[sample] >= total * opt->minorAllelePct / 100)
		call->isRefIAllele[sample] = 1;
	for (j = 0; j < kMAXNBINDELS; j++)
	{
		if (var->del[j].cnt == 0)
			break;
		if (var->del[j].len > call->maxDelLen)
		{
			strncpy(call->delRef,(const char *) genome_tbl + genomepos + i - 1,var->del[j].len + 1);
			call->maxDelLen = var->del[j].len;
		}
		unsigned int d = 0;
		while (d < call->nbD && call->del[d].len != var->del[j].len)
			d += 1;
		if (d >= call->nbD)
		{
			if (d >= kMAXCALLS)
			{
				fprintf(stderr,"nbD >= kMAXCALLS\n");
				exit(1);
			}
			call->nbD += 1;
			call->del[d].len = var->del[j].len;
		}
		call->del[d].cnt[sample] = var->del[j].cnt;
		if (var->del[j].cnt >= opt->reportVCF && var->del[j].cnt >= total * opt->minorAllelePct / 100)
		{
			call->del[d].isAllele[sample] = 1;
			call->nbDal += 1;
		}
	}
	for (j = 0; j < 4; j++)
	{
		if (gNT2bits[call->refNt] == j)
			continue;
		call->M[j].cnt[sample] = var->ACGTP[j];
		if (var->ACGTP[j] >= opt->reportVCF && var->ACGTP[j] >= total * opt->minorAllelePct / 100)
		{
			call->M[j].isAllele[sample] = 1;
			// FIXME - only works while we only deal with 2 samples
			if (sample == 0 || call->M[j].isAllele[0] == 0)
				call->nbM += 1;
			assert(call->nbM <= 3);
		}
	}
}
//---------------------------------------------------------------

static int
cmpAl(const void *p1, const void *p2)
{
	alCnt_p_t a = (alCnt_p_t) p1;
	alCnt_p_t b = (alCnt_p_t) p2;
	if (a->cnt < b->cnt)
		return 1;
	if (a->cnt > b->cnt)
		return -1;
	if (a->al < b->al)
		return -1;
	if (a->al > b->al)
		return 1;
	return 0;
}
//---------------------------------------------------------------

static void
printDel(calledVar_p_t call, OPTIONS *opt)
{
	// determine longest called
	call->maxDelLen = 0;
	for (unsigned int i = 0; i < call->nbD; i++)
		if (call->del[i].isAllele[0] > 0 || call->del[i].isAllele[1] > 0)
			if (call->del[i].len > call->maxDelLen)
				call->maxDelLen = call->del[i].len;
	printf("%s\t%u\t.\t%.*s\t",
		virtchr[opt->chr].AC,
		call->pos - 1,
		call->maxDelLen + 1,
		call->delRef);
	unsigned int v = 0;
	unsigned int zeroCalled = call->isRefAllele[0];
	unsigned int total[kMAXSAMPLE] = { call->refCnt[0], call->refCnt[1] };
	for (unsigned int i = 0; i < call->nbD; i++)
		if (call->del[i].isAllele[0] > 0 || call->del[i].isAllele[1] > 0)
		{
			printf("%s%.*s",
				(v > 0) ? "," : "",
				call->maxDelLen + 1 - call->del[i].len,
				call->delRef);
			if (call->del[i].isAllele[0] > 0)
				zeroCalled = 1;
			v += 1;
			total[0] += call->del[i].cnt[0];
			total[1] += call->del[i].cnt[1];
		}
	unsigned int somatic = 0;
	if (zeroCalled)
		for (unsigned int i = 0; i < call->nbD; i++)
			if (call->del[i].isAllele[1] > 0 && call->del[i].cnt[0] <= ((call->del[i].cnt[1] * total[0]) / total[1]) * opt->minorAllelePct / 100)
				somatic = 1;
	printf("\t.\tPASS\t%s\tGT:DP:AD",
		somatic ? "SOMATIC" : ".");
	const unsigned int nbAl = v + 1;
	alCnt_t ac[nbAl];
	for (unsigned int s = 0; s < kMAXSAMPLE; s++)
	{
		ac[0].al = 0;
		ac[0].cnt = call->refCnt[s];
		ac[0].called = call->isRefAllele[s];
		v = 1;
		for (unsigned int i = 0; i < call->nbD; i++)
			if (call->del[i].isAllele[0] > 0 || call->del[i].isAllele[1] > 0)
			{
				ac[v].al = v;
				ac[v].cnt = call->del[i].cnt[s];
				ac[v].called = call->del[i].isAllele[s];
				v += 1;
			}
		qsort(ac,nbAl,sizeof(alCnt_t),cmpAl);
		if (ac[1].called)
			printf("\t%u/%u:%u:%u",
				(ac[0].al < ac[1].al) ? ac[0].al : ac[1].al,
				(ac[0].al < ac[1].al) ? ac[1].al : ac[0].al,
				total[s],
				call->refCnt[s]);
		else
			printf("\t%u/%u:%u:%u",
				ac[0].called ? ac[0].al : 0,
				ac[0].called ? ac[0].al : 0,
				total[s],
				call->refCnt[s]);
		for (unsigned int i = 0; i < call->nbD; i++)
			if (call->del[i].isAllele[0] > 0 || call->del[i].isAllele[1] > 0)
				printf(",%u",
					call->del[i].cnt[s]);
	}
	printf("\n");
}
//---------------------------------------------------------------

static void
printIns(calledVar_p_t call, OPTIONS *opt)
{
	printf("%s\t%u\t.\t%c\t",
		virtchr[opt->chr].AC,
		call->pos - 1,
		call->refINt);
	unsigned int v = 0;
	unsigned int zeroCalled = call->isRefIAllele[0];
	unsigned int total[kMAXSAMPLE] = { call->refICnt[0], call->refICnt[1] };
	for (unsigned int i = 0; i < call->nbI; i++)
		if (call->ins[i].isAllele[0] > 0 || call->ins[i].isAllele[1] > 0)
		{
			printf("%s%c%s",
				(v > 0) ? "," : "",
				call->refINt,
				call->ins[i].nt);
			if (call->ins[i].isAllele[0] > 0)
				zeroCalled = 1;
			v += 1;
			total[0] += call->ins[i].cnt[0];
			total[1] += call->ins[i].cnt[1];
		}
	unsigned int somatic = 0;
	if (zeroCalled)
		for (unsigned int i = 0; i < call->nbI; i++)
			if (call->ins[i].isAllele[1] > 0 && call->ins[i].cnt[0] <= ((call->ins[i].cnt[1] * total[0]) / total[1]) * opt->minorAllelePct / 100)
				somatic = 1;
	printf("\t.\tPASS\t%s\tGT:DP:AD",
		somatic ? "SOMATIC" : ".");
	const unsigned int nbAl = v + 1;
	alCnt_t ac[nbAl];
	for (unsigned int s = 0; s < kMAXSAMPLE; s++)
	{
		ac[0].al = 0;
		ac[0].cnt = call->refICnt[s];
		ac[0].called = call->isRefIAllele[s];
		v = 1;
		for (unsigned int i = 0; i < call->nbI; i++)
			if (call->ins[i].isAllele[0] > 0 || call->ins[i].isAllele[1] > 0)
			{
				ac[v].al = v;
				ac[v].cnt = call->ins[i].cnt[s];
				ac[v].called = call->ins[i].isAllele[s];
				v += 1;
			}
		qsort(ac,nbAl,sizeof(alCnt_t),cmpAl);
		if (ac[1].called)
			printf("\t%u/%u:%u:%u",
				(ac[0].al < ac[1].al) ? ac[0].al : ac[1].al,
				(ac[0].al < ac[1].al) ? ac[1].al : ac[0].al,
				total[s],
				call->refICnt[s]);
		else
			printf("\t%u/%u:%u:%u",
				ac[0].called ? ac[0].al : 0,
				ac[0].called ? ac[0].al : 0,
				total[s],
				call->refICnt[s]);
		for (unsigned int i = 0; i < call->nbI; i++)
			if (call->ins[i].isAllele[0] > 0 || call->ins[i].isAllele[1] > 0)
				printf(",%u",
					call->ins[i].cnt[s]);
	}
	printf("\n");
}
//---------------------------------------------------------------

static void
printMM(calledVar_p_t call, OPTIONS *opt)
{
	printf("%s\t%u\t.\t%c\t",
		virtchr[opt->chr].AC,
		call->pos,
		call->refNt);
	unsigned int v = 0;
	unsigned int zeroCalled = call->isRefAllele[0];
	unsigned int total[kMAXSAMPLE] = { call->refCnt[0], call->refCnt[1] };
	for (unsigned int i = 0; i < 4; i++)
		if (call->M[i].isAllele[0] > 0 || call->M[i].isAllele[1] > 0) // will be wrong for the ref allele...
		{
			printf("%s%c",
				(v > 0) ? "," : "",
				"ACGT"[i]);
			if (call->M[i].isAllele[0] > 0)
				zeroCalled = 1;
			v += 1;
			total[0] += call->M[i].cnt[0];
			total[1] += call->M[i].cnt[1];
		}
	unsigned int somatic = 0;
	if (zeroCalled)
		for (unsigned int i = 0; i < 4; i++)
			if (call->M[i].isAllele[1] > 0 && call->M[i].cnt[0] <= ((call->M[i].cnt[1] * total[0]) / total[1]) * opt->minorAllelePct / 100)
				somatic = 1;
	printf("\t.\tPASS\t%s\tGT:DP:AD",
		somatic ? "SOMATIC" : ".");
	const unsigned int nbAl = v + 1;
	alCnt_t ac[nbAl];
	for (unsigned int s = 0; s < kMAXSAMPLE; s++)
	{
		ac[0].al = 0;
		ac[0].cnt = call->refCnt[s];
		ac[0].called = call->isRefAllele[s];
		v = 1;
		for (unsigned int i = 0; i < 4; i++)
			if (call->M[i].isAllele[0] > 0 || call->M[i].isAllele[1] > 0)
			{
				ac[v].al = v;
				ac[v].cnt = call->M[i].cnt[s];
				ac[v].called = call->M[i].isAllele[s];
				v += 1;
			}
		qsort(ac,nbAl,sizeof(alCnt_t),cmpAl);
		if (ac[1].called)
			printf("\t%u/%u:%u:%u",
				(ac[0].al < ac[1].al) ? ac[0].al : ac[1].al,
				(ac[0].al < ac[1].al) ? ac[1].al : ac[0].al,
				total[s],
				call->refCnt[s]);
		else
			printf("\t%u/%u:%u:%u",
				ac[0].called ? ac[0].al : 0,
				ac[0].called ? ac[0].al : 0,
				total[s],
				call->refCnt[s]);
		for (unsigned int i = 0; i < 4; i++)
			if (call->M[i].isAllele[0] > 0 || call->M[i].isAllele[1] > 0)
				printf(",%u",
					call->M[i].cnt[s]);
	}
	printf("\n");
}
//---------------------------------------------------------------

static void
diff(coverage_p_t cov1, coverage_p_t cov2, reference_p_t ref1, reference_p_t ref2, OPTIONS *opt)
{
	unsigned int i;
	unsigned int notN = 0;
	unsigned int callable = 0;
	unsigned int genomepos = (virtchr[opt->chr].chr << 28) + virtchr[opt->chr].offset;
	if (opt->reportVCF == 0)
	{
		printf("chr\tchrPos\n");
		for (i = 1; i <= virtchr[opt->chr].len; i++)
		{
			if (genome_tbl[genomepos+i] != 'N')
				notN += 1;
			unsigned int tot1 = cov1[i].cntP + cov1[i].cntM;
			unsigned int tot2 = cov2[i].cntP + cov2[i].cntM;
			if (cov1[i].var != NULL)
				for (unsigned int j = 0; j < 4; j++)
					tot1 += cov1[i].var->ACGTP[j] + cov1[i].var->ACGTM[j];
			if (cov2[i].var != NULL)
				for (unsigned int j = 0; j < 4; j++)
					tot2 += cov2[i].var->ACGTP[j] + cov2[i].var->ACGTM[j];
			if (tot1 >= minAllCov && tot2 >= minAllCov)
				callable += 1;
			if (cov1[i].var == NULL && cov2[i].var == NULL)
				continue;
			//int is_var1 = is_significant_p(cov1 + i);
			int is_var2 = is_significant_p(cov2 + i);
			int printit = 0;
			//if (is_var1)
			//	printit = is_diff_p(cov1 + i, cov2 + i);
			if (is_var2)
				printit = is_diff_p(cov2 + i, cov1 + i);
			if (printit)
				printf("%s\t%d\n",virtchr[opt->chr].SAMname,i);
		}
		fprintf(stderr,"chr:%u len:%u notN:%u callable:%u (%.2f %.2f)\n",opt->chr,virtchr[opt->chr].len,notN,callable,(callable * 100.0)/(virtchr[opt->chr].len * 1.0),(callable * 100.0)/(notN * 1.0));
		return;
	}
	// compute median and iqr
	unsigned int covCntNb1, median1, iqr251, iqr751;
	unsigned int covCntNb2, median2, iqr252, iqr752;
	unsigned int covCntNbSel1, medianSel1, iqr25Sel1, iqr75Sel1;
	unsigned int covCntNbSel2, medianSel2, iqr25Sel2, iqr75Sel2;
	computeCov(cov1, opt, NULL, &covCntNb1, &median1, &iqr251, &iqr751);
	computeCov(cov2, opt, NULL, &covCntNb2, &median2, &iqr252, &iqr752);
	if (opt->selectRegions[0] != 0)
		computeCov(cov1, opt, cov1, &covCntNbSel1, &medianSel1, &iqr25Sel1, &iqr75Sel1);
	if (opt->selectRegions[0] != 0)
		computeCov(cov2, opt, cov1, &covCntNbSel2, &medianSel2, &iqr25Sel2, &iqr75Sel2);

	printf("##fileformat=VCFv4.1\n");
	if (opt->selectRegions[0] != 0)
	{
		printf("##SAMPLE=<ID=%s,CONTIG=%s,IQR25=%u,MEDIAN=%u,IQR75=%u,CNT=%u,SEL=%s>\n",opt->sampleName1,virtchr[opt->chr].AC,iqr25Sel1,medianSel1,iqr75Sel1,covCntNbSel1,opt->selectRegions);
		printf("##SAMPLE=<ID=%s,CONTIG=%s,IQR25=%u,MEDIAN=%u,IQR75=%u,CNT=%u,SEL=%s>\n",opt->sampleName2,virtchr[opt->chr].AC,iqr25Sel2,medianSel2,iqr75Sel2,covCntNbSel2,opt->selectRegions);
	}
	printf("##SAMPLE=<ID=%s,CONTIG=%s,IQR25=%u,MEDIAN=%u,IQR75=%u,CNT=%u>\n",opt->sampleName1,virtchr[opt->chr].AC,iqr251,median1,iqr751,covCntNb1);
	printf("##SAMPLE=<ID=%s,CONTIG=%s,IQR25=%u,MEDIAN=%u,IQR75=%u,CNT=%u>\n",opt->sampleName2,virtchr[opt->chr].AC,iqr252,median2,iqr752,covCntNb2);
	printf("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth\">\n");
	printf("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
	printf("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Per-sample read depths for each allele\">\n");
	printf("##INFO=<ID=varType,Number=1,Type=String,Description=\"Type of variant\">\n");
	for (unsigned int c = 1; c <= chrMax; c += 1)
		printf("##contig=<ID=%s,length=%u,assembly=%s>\n",virtchr[c].AC,virtchr[c].len,opt->refName);
	printf("##phasing=none\n");
	printf("##source=SelectVariants\n");
	printf("##variants_justified=left\n");
	printf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\t%s\n",opt->sampleName1,opt->sampleName2);
	for (i = 0; i <= virtchr[opt->chr].len; i++)
	{
		if (gNT2bits[genome_tbl[genomepos+i]] == kNOT_ACGT)
		{
			free(cov1[i].var);
			cov1[i].var = NULL;
			free(cov2[i].var);
			cov2[i].var = NULL;
			continue;
		}
		calledVar_t call;
		memset(&call,0,sizeof(calledVar_t));
		doCall(cov1[i].var,cov1[i].cntP,cov1[i].cntM,&call,0,genomepos,i,opt);
		doCall(cov2[i].var,cov2[i].cntP,cov2[i].cntM,&call,1,genomepos,i,opt);
		if (call.nbDal + call.nbIal + call.nbM == 0)
		{
			free(cov1[i].var);
			cov1[i].var = NULL;
			free(cov2[i].var);
			cov2[i].var = NULL;
			continue;
		}
		if (call.nbDal > 0)
			printDel(&call,opt);
		if (call.nbIal > 0)
			printIns(&call,opt);
		if (call.nbM > 0)
			printMM(&call,opt);
		if (cov1[i].var != NULL)
			reportRef(opt->ref1,cov1[i].var,genome_tbl[genomepos+i],ref1,i,opt);
		if (cov2[i].var != NULL)
			reportRef(opt->ref2,cov2[i].var,genome_tbl[genomepos+i],ref2,i,opt);
	}
}
//---------------------------------------------------------------

int
main (int argc, char **argv)
{
	char rsltfile1[512];
	char rsltfile2[512];
	char genomefile[512];
	OPTIONS options;
	const char * restrict AlignerScoreConfFile = NULL;
	int c;
	int err;
	memset(&options,0,sizeof(OPTIONS));
	options.maxPairLen = kMAXPAIRLEN;

	/* --------- process arguments ----------------------*/

	rsltfile1[0] = 0;
	rsltfile2[0] = 0;
	strcpy(genomefile,"/tmp/nodelete/nguex/data/hg19.bin");
	opterr = 0;
	while ((c = getopt (argc, argv, "r:s:g:C:v:dacf:N:O:b:m:V:R:S:W:F:D:")) != -1)
	switch (c)
	{
	  case 'r':
			strcpy(rsltfile1,optarg);
			break;

	  case 's':
			strcpy(rsltfile2,optarg);
			break;

		case 'R':
			options.ref1 = fopen(optarg,"w");
			if (options.ref1 == NULL)
				perror("fopen");
			break;

		case 'S':
			options.ref2 = fopen(optarg,"w");
			if (options.ref2 == NULL)
				perror("fopen");
			break;

	  case 'C':
			sscanf(optarg,"%d",&options.chr);
			break;

		case 'g':
			strcpy(genomefile,optarg);
			break;

		case 'a':
			options.additionalFilter = 1;
			break;

		case 'c':
			options.noCloneFilter = 1;
			break;

		case 'f':
			options.filterQualValue = optarg[0] - '#';
			break;

	  case 'd':
			debug = 1;
			break;

	  case 'v':
			sscanf(optarg,"%d",&verbose);
			break;

		case 'N':
			strncpy(options.sampleName1,optarg,64);
			break;

		case 'O':
			strncpy(options.sampleName2,optarg,64);
			break;

		case 'b':
			strcpy(options.selectRegions,optarg);
			break;

		case 'm':
			options.minorAllelePct = atoi(optarg);
			break;

		case 'V':
			options.reportVCF = atoi(optarg);
			break;

		case 'F':
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

		case 'W':
			options.maxPairLen = atoi(optarg);
			break;
	}

	if (rsltfile1[0] == 0 || rsltfile2[0] == 0 || options.chr == -1)
	{
		printf("usage:\n\n");
		printf("GTLdiffCaller [-g genomefile] -C chromosome -r resultfile1 -s resultfile2 [-v level]\n\n");
		printf("           -g genomefile         : full path of the genome encoded by encode_genome\n");
		printf("                                   default: '%s' \n",genomefile);
		printf("           -b <file>             : selected regions to compute coverage (bed-style)\n");
		printf("           -C chromosome         : select chromosome to analyze\n");
		printf("           -m <int>              : minor allele minimum percent coverage\n");
		printf("           -a                    : perform additional filtering, discard pairs with more than 2 indel or indel longer than 40 nt (default is off)\n");
		printf("           -c                    : do not filter clones (default is to filter them out)\n");
		printf("           -f <char>             : filter out nucleotides where the quality value is *below* the specified character [%c]\n", options.filterQualValue + '#');
		printf("           -r resultfile1        : name of the first result file (produced by match)\n");
		printf("           -s resultfile2        : name of the second result file (produced by match)\n");
		printf("           -R reffile1           : output reference for variant calls of first sample\n");
		printf("           -S reffile2           : output reference for variant calls of second sample\n");
		printf("           -V <int>              : report as VCF, insertions, deletions and SNPs; value is minimum observation count\n");
		printf("           -N sampleName         : name of first sample to put in VCF file\n");
		printf("           -O sampleName         : name of second sample to put in VCF file\n");
		printf("           -d                    : debug mode...\n");
		printf("           -v level              : verbose level\n");
		printf("           -F alignconffile      : aligner configuration file\n");
		printf("           -D alignconffile      : dump alignment score JSON file\n");
		printf("           -W <int>              : maximum length of 'well-behaved' pairs (distance and relative strand-orientation); 0 means no check [%u]\n", options.maxPairLen);
		printf("\n");
		printf(GTL_VERSION_FULL_STRING "\n");
		printf("(c) Nicolas Guex & Christian Iseli & Cie 2014-2018\n");
		return(1);
	}

	if (options.sampleName1[0] == 0)
		snprintf(options.sampleName1, 64, "%s", rsltfile1);
	if (options.sampleName2[0] == 0)
		snprintf(options.sampleName2, 64, "%s", rsltfile2);

	char *configfn = strdup(genomefile);
	unsigned int len = strlen(configfn);
	memcpy(configfn + len - 3,"cfg",3);

	err = loadGenome(genomefile,&options);
	if (err)
		goto bail;

	initializeVirtualChromosomes(configfn,virtchr);
	initRevCompTable(gIUPACrevcomp);

	memset(gNT2bits,kNOT_ACGT,128);
	gNT2bits['A'] = 0x0;
	gNT2bits['C'] = 0x1;
	gNT2bits['G'] = 0x2;
	gNT2bits['T'] = 0x3;

	coverage_p_t cov1 = calloc(virtchr[options.chr].len + kMaxReadLen + 1,sizeof(coverage_t));
	// FIXME - maybe someday look into circular DNA, but for now just allocate some extra space at the end of the piece
	reference_p_t ref1 = calloc((virtchr[options.chr].len + 1025) >> 8,sizeof(reference_t));
	loadSelected(cov1,&options);

	int (*decompressFct)(const char * const restrict, OPTIONS*, coverage_p_t, reference_p_t) = (!AlignerScoreConfFile) ? decompress : decompress_realign;

	err = decompressFct(rsltfile1,&options,cov1,ref1);
	if (err)
		goto bail;
	coverage_p_t cov2 = calloc(virtchr[options.chr].len + kMaxReadLen + 1,sizeof(coverage_t));
	reference_p_t ref2 = calloc((virtchr[options.chr].len + 1025) >> 8,sizeof(reference_t));
	err = decompressFct(rsltfile2,&options,cov2,ref2);
	if (err)
		goto bail;

	diff(cov1,cov2,ref1,ref2,&options);

	free(cov1);
	free(cov2);
	err = 0;

	/* ------------------------ clean up and quit ------------------- */

bail:;

	if (genome_tbl != MAP_FAILED)
		munmap(genome_tbl,kGENOME_DATA_SIZE);

	return(err);

} /* main */
//---------------------------------------------------------------

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
