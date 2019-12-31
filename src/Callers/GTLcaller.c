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

#define _GNU_SOURCE     /* Expose declaration of tdestroy() */
#include <search.h>
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
#include "fastalign.h"
#include "GTLcaller_common.h"

#include "htslib/vcf.h"
#include "htslib/vcfutils.h"


//=================================================================================================

//static struct timeval ts;
//static struct timeval te;
static int verbose = 0;
static int debug = 0;

//=================================================================================================

static void
outputPWM(OPTIONS *opt, unsigned int pwm, coverage_p_t cov, unsigned int genomepos)
{
	if (opt->chr != opt->PWMpos[pwm].chr)
		return;
	unsigned int min = 1;
	if (opt->PWMpos[pwm].pos > 60)
		min = opt->PWMpos[pwm].pos - 60;
	unsigned int max = opt->PWMpos[pwm].pos + 60;
	if (max > virtchr[opt->chr].len)
		max = virtchr[opt->chr].len;
	char fn[544];
	sprintf(fn,"%schr%u_%u.pwm",opt->outbase,opt->chr,min);
	FILE *f = fopen(fn,"w");
	if (!f)
	{
		fprintf(stderr,"Error: Cannot create %s\n",fn);
		perror("fopen");
		exit(1);
	}
	unsigned int i;
	fprintf(f,"PO\tA\tC\tG\tT\n");
	for (i = min; i <= max; i++)
	{
		fprintf(f,"%02u",i - min + 1);
		unsigned int j;
		for (j = 0; j < 4; j++)
		{
			unsigned int cnt = 0;
			if (genome_tbl[genomepos+i] == "ACGT"[j])
				cnt += cov[i].cntP + cov[i].cntM;
			if (cov[i].var != NULL)
				cnt += cov[i].var->ACGTP[j] + cov[i].var->ACGTM[j];
			fprintf(f,"\t%u",cnt);
		}
		fprintf(f,"\n");
	}
	fclose(f);
}
//---------------------------------------------------------------

static void
printWiggle(coverage_p_t cov, OPTIONS *opt)
{
	FILE *f = fopen (opt->wiggleBed, "w");
	if (!f)
	{
		perror("fopen");
		fprintf(stderr, "Error:Cannot open %s\n",opt->wiggleBed);
		exit(1);
	}
	fprintf(f,"browser position %s:1-%u\n",virtchr[opt->chr].AC,virtchr[opt->chr].len);
	fprintf(f,"browser hide all\n");
	fprintf(f,"browser squish sibGene refGene\n");
	fprintf(f,"browser dense rmsk\n");
	fprintf(f,"track type=wiggle_0 name=\"%s_%s\" description=\"Coverage %s\" visibility=full color=0,200,100 maxHeightPixels=100:50:20 graphType=points priority=10\n",opt->sampleName1,virtchr[opt->chr].AC,opt->sampleName1);
	unsigned int prev = 0, start;
	for (unsigned int i = 0; i <= virtchr[opt->chr].len; i++)
	{
		unsigned int c = cov[i].cov[0] + cov[i].cov[1] + cov[i].cov[2];
		if (c < kMINCOVWIGGLEREPORT)
			c = 0;
		if (c != prev)
		{
			if (prev != 0)
			{
				fprintf(f,"%s\t%u\t%u\t%u\n",virtchr[opt->chr].AC,start-1,i-1,prev);
			}
			prev = c;
			start = i;
		}
	}
	fclose(f);
}
//---------------------------------------------------------------

static int
cmpint(const void *p1, const void *p2)
{
	unsigned int a = * (unsigned int *) p1;
	unsigned int b = * (unsigned int *) p2;
	if (a > b)
		return 1;
	if (a < b)
		return -1;
	return 0;
}
//---------------------------------------------------------------

static void
getCovStat(coverage_p_t cov, unsigned int start, unsigned int end, covStat_p_t csp)
{
	const unsigned int aLen = end - start + 1;
	unsigned int a[aLen];
	unsigned int total;
	for (unsigned int k = 0; k < 3; k++)
	{
		total = 0;
		for (unsigned int i = 0; i < aLen; i++)
		{
			a[i] = cov[i + start].cov[k];
			total += a[i];
		}
		qsort(a, aLen, sizeof(unsigned int), cmpint);
		csp->max[k] = a[aLen - 1];
		csp->mean[k] = total / aLen;
		csp->median[k] = a[aLen / 2];
	}
	total = 0;
	for (unsigned int i = 0; i < aLen; i++)
	{
		a[i] = cov[i + start].cov[0] + cov[i + start].cov[1] + cov[i + start].cov[2];
		total += a[i];
	}
	qsort(a, aLen, sizeof(unsigned int), cmpint);
	csp->max[3] = a[aLen - 1];
	csp->mean[3] = total / aLen;
	csp->median[3] = a[aLen / 2];
}
//---------------------------------------------------------------

typedef struct _exonList {
	unsigned int p1, p2;
} exonList_t, *exonList_p_t;

static void
getElCovStat(coverage_p_t cov, exonList_p_t elp, unsigned int cnt, const unsigned int aLen, covStat_p_t csp)
{
	unsigned int a[aLen];
	unsigned int total;
	for (unsigned int k = 0; k < 3; k++)
	{
		total = 0;
		unsigned int i = 0;
		for (unsigned int e = 0; e < cnt; e++)
		{
			unsigned int start = elp[e].p1;
			unsigned int end = elp[e].p2;
			while (start <= end)
			{
				assert(i < aLen);
				a[i] = cov[start].cov[k];
				total += a[i];
				start += 1;
				i += 1;
			}
		}
		qsort(a, aLen, sizeof(unsigned int), cmpint);
		csp->max[k] = a[aLen - 1];
		csp->mean[k] = total / aLen;
		csp->median[k] = a[aLen / 2];
	}
	total = 0;
	unsigned int i = 0;
	for (unsigned int e = 0; e < cnt; e++)
	{
		unsigned int start = elp[e].p1;
		unsigned int end = elp[e].p2;
		while (start <= end)
		{
			a[i] = cov[start].cov[0] + cov[start].cov[1] + cov[start].cov[2];
			total += a[i];
			start += 1;
			i += 1;
		}
	}
	qsort(a, aLen, sizeof(unsigned int), cmpint);
	csp->max[3] = a[aLen - 1];
	csp->mean[3] = total / aLen;
	csp->median[3] = a[aLen / 2];
}
//---------------------------------------------------------------

static void
printTromer(coverage_p_t cov, reference_p_t ref, OPTIONS *opt)
{
	char tromefn[512];
	strcpy(tromefn,opt->expressionBase);
	strcpy(tromefn + (strlen(tromefn) - 9),"_trome.txt");
	FILE *ft = fopen(tromefn,"w");
	if (!ft)
	{
		perror("fopen");
		fprintf(stderr, "Error: Cannot open %s\n",tromefn);
		exit(1);
	}
	FILE *tr = fopen(opt->tromerTranscripts, "r");
	const unsigned int bufSize = 512*1024;
	char s[bufSize];
	exonList_t el[512];
	while (fgets(s,bufSize,tr) != NULL)
	{
		// >map|NC_000010_0_0|HTR045305|+ E0,E5,E6 NC_000010.11[10748..11242,12798..13000,15263..15996] (LOC102723376) Homo sapiens uncharacterized LOC102723376 (LOC102723376), long non-coding RNA. 3PE_BX089548; LEN=1432
		// >map|NC_000010_1_0|HTR034524|- E4,E2,E0 NC_000010.11[13787..14299,16502..16544,20786..20815] unknown 3P073890; LEN=586
		char ID[64], name[64];
		unsigned int cnt = 0;
		unsigned int len;
		if (strncmp(s,">map|",5))
			continue;
		unsigned int i;
		for (i = 0; i < 63 && s[i+5] != '|'; i++)
			ID[i] = s[i+5];
		ID[i] = 0;
		char *start = strchr(s,'[');
		if (start == NULL)
			continue;
		start += 1;
		char *end;
		do
		{
			el[cnt].p1 = strtoul(start,&end,10);
			if (*end != '.')
			{
				fprintf(stderr,"error parsing %s\n",s);
				exit(1);
			}
			start = end+2;
			el[cnt].p2 = strtoul(start,&end,10);
			if (*end != ',' && *end != ']')
			{
				fprintf(stderr,"error parsing %s\n",s);
				exit(1);
			}
			start = end+1;
			cnt += 1;
			if (cnt >= 512)
			{
				fprintf(stderr,"too many exons in %s\n",tromefn);
				exit(1);
			}
		} while (*end != ']');
		start = strchr(start,'(');
		i = 0;
		if (start != NULL)
		{
			start += 1;
			while (i < 63 && *start != ';' && *start != ')')
				name[i++] = *start++;
			name[i] = 0;
		}
		else
			strcpy(name,"Unknown");
		start = strstr(s," LEN=");
		if (start == NULL)
		{
			fprintf(stderr,"missing LEN= info in %s\n",s);
			exit(1);
		}
		len = strtoul(start+5,NULL,10);
		covStat_t cs;
		getElCovStat(cov,el,cnt,len,&cs);
		if (cs.max[3] == 0)
			continue;
		unsigned int readCount = 0;
		unsigned int p1 = el[0].p1;
		unsigned int p2 = el[cnt-1].p2;
		unsigned int posH = p1 >> 8;
		if (posH > 0)
			posH -= 1;
		while (posH <= (p2 >> 8))
		{
			reference_p_t r = ref + posH;
			unsigned int offset1 = p1 - (posH << 8);
			if (p1 < (posH << 8))
				offset1 = 0;
			unsigned int offset2 = p2 - (posH << 8);
			for (unsigned int i = 0; i < r->cnt; i++)
			{
				if (r->rr[i].offset <= offset2 && r->rr[i].offset + r->rr[i].len > offset1)
					readCount += 1;
			}
			posH += 1;
		}
		fprintf(ft,"%s\t%u\t%u\t%s\t%s\t%u\t%u\t%u",virtchr[opt->chr].AC,el[0].p1,el[cnt-1].p2,ID,name,len,cnt,readCount);
		for (unsigned int k = 0; k < 4; k++)
			fprintf(ft,"\t%u\t%u\t%u",cs.max[k],cs.mean[k],cs.median[k]);
		fputc('\n',ft);
	}
	fclose(tr);
	fclose(ft);
}
//---------------------------------------------------------------

static void
printExpr(coverage_p_t cov, reference_p_t ref, OPTIONS *opt)
{
	char wcovfn[512];
	strcpy(wcovfn,opt->expressionBase);
	strcpy(wcovfn + (strlen(wcovfn) - 9),"_wcov.txt");
	FILE *fw = fopen(wcovfn,"w");
	if (!fw)
	{
		perror("fopen");
		fprintf(stderr, "Error: Cannot open %s\n",wcovfn);
		exit(1);
	}
	unsigned int pos = 1;
	while (pos + 199 <= virtchr[opt->chr].len)
	{
		covStat_t cs;
		getCovStat(cov,pos,pos+199,&cs);
		fprintf(fw,"%s\t%u\t%u",virtchr[opt->chr].AC,pos,pos+199);
		for (unsigned int k = 0; k < 4; k++)
			fprintf(fw,"\t%u\t%u\t%u",cs.max[k],cs.mean[k],cs.median[k]);
		fputc('\n',fw);
		pos += 100;
	}
	fclose(fw);
	FILE *fg = fopen(opt->expressionBase,"w");
	if (!fg)
	{
		perror("fopen");
		fprintf(stderr, "Error: Cannot open %s\n",opt->expressionBase);
		exit(1);
	}
	char exonfn[512];
	strcpy(exonfn,opt->expressionBase);
	strcpy(exonfn + (strlen(exonfn) - 9),"_exon.txt");
	FILE *fe = fopen(exonfn,"w");
	if (!fe)
	{
		perror("fopen");
		fprintf(stderr, "Error: Cannot open %s\n",exonfn);
		exit(1);
	}
	gzFile gff = gzopen(opt->gencodeGFF3, "rb");
	char s[512];
	while (gzgets(gff,s,512) != NULL)
	{
		char chr[32], elt[32];
		unsigned int p1, p2;
		int res = sscanf(s,"%31s %*s %31s %u %u",chr,elt,&p1,&p2);
		if (res != 4)
			continue;
		if (strcmp(chr,virtchr[opt->chr].AC))
			continue;
		if (strcmp(elt,"gene") && strcmp(elt,"exon"))
			continue;
		char ID[64], name[64];
		char *start = strstr(s,"gene_id=");
		if (start == NULL)
			continue;
		start += 8;
		size_t len = strcspn(start,";\t\n");
		strncpy(ID,start,len);
		ID[len] = 0;
		start = strstr(s,"gene_name=");
		if (start != NULL)
		{
			start += 10;
			len = strcspn(start,";\t\n");
			strncpy(name,start,len);
			name[len] = 0;
		}
		else
			name[0] = 0;
		covStat_t cs;
		getCovStat(cov,p1,p2,&cs);
		if (strcmp(elt,"gene") == 0)
		{
			// we assume we see the gene line before the exons line...
			fprintf(fg,"%s\t%u\t%u\t%s\t%s",chr,p1,p2,ID,name);
			for (unsigned int k = 0; k < 4; k++)
				fprintf(fg,"\t%u\t%u\t%u",cs.max[k],cs.mean[k],cs.median[k]);
			fputc('\n',fg);
			continue;
		}
		if (cs.max[3] == 0)
			continue;
		unsigned int exon = 0;
		start = strstr(s,"exon_number=");
		if (start != NULL)
		{
			start += 12;
			exon = strtoul(start,NULL,0);
		}
		fprintf(fe,"%s\t%u\t%u\t%s\t%s\t%u",chr,p1,p2,ID,name,exon);
		for (unsigned int k = 0; k < 4; k++)
			fprintf(fe,"\t%u\t%u\t%u",cs.max[k],cs.mean[k],cs.median[k]);
		unsigned int posH = p1 >> 8;
		if (posH > 0)
			posH -= 1;
		char sep = '\t';
		while (posH <= (p2 >> 8))
		{
			reference_p_t r = ref + posH;
			unsigned int offset1 = p1 - (posH << 8);
			if (p1 < (posH << 8))
				offset1 = 0;
			unsigned int offset2 = p2 - (posH << 8);
			for (unsigned int i = 0; i < r->cnt; i++)
			{
				if (r->rr[i].offset <= offset2 && r->rr[i].offset + r->rr[i].len > offset1)
				{
					char c;
					if (r->rr[i].second == 0)
					{
						if (r->rr[i].reverse == 0)
							c = 'L';
						else
							c = 'l';
					}
					else
					{
						if (r->rr[i].reverse == 0)
							c = 'r';
						else
							c = 'R';
					}
					fprintf(fe,"%c%lu%c",sep,(unsigned long) r->rr[i].ordinal,c);
					sep = ';';
				}
			}
			posH += 1;
		}
		fputc('\n',fe);
	}
	gzclose(gff);
	fclose(fe);
	fclose(fg);
}
//---------------------------------------------------------------

static void
report(OPTIONS *opt, coverage_p_t cov, reference_p_t ref)
{
	unsigned int genomepos = (virtchr[opt->chr].chr << 28) + virtchr[opt->chr].offset;
	unsigned int i;
	if (opt->PWMcnt > 0)
	{
		for (i = 0; i < opt->PWMcnt; i++)
			outputPWM(opt, i, cov, genomepos);
		return;
	}

	if (opt->wiggleBed[0] != 0)
		printWiggle(cov,opt);

	if (opt->gencodeGFF3[0] != 0)
		printExpr(cov,ref,opt);

	if (opt->tromerTranscripts[0] != 0)
		printTromer(cov,ref,opt);

	// compute median and iqr
	unsigned int covCntNb, median, iqr25, iqr75;
	unsigned int covCntNbSel, medianSel, iqr25Sel, iqr75Sel;
	computeCov(cov, opt, NULL, &covCntNb, &median, &iqr25, &iqr75);
	if (opt->selectRegions[0] != 0)
		computeCov(cov, opt, cov, &covCntNbSel, &medianSel, &iqr25Sel, &iqr75Sel);

	if (opt->reportVCF > 0)
	{
		printf("##fileformat=VCFv4.1\n");
		if (opt->selectRegions[0] != 0)
			printf("##SAMPLE=<ID=%s,CONTIG=%s,IQR25=%u,MEDIAN=%u,IQR75=%u,CNT=%u,SEL=%s>\n",opt->sampleName1,virtchr[opt->chr].AC,iqr25Sel,medianSel,iqr75Sel,covCntNbSel,opt->selectRegions);
		printf("##SAMPLE=<ID=%s,CONTIG=%s,IQR25=%u,MEDIAN=%u,IQR75=%u,CNT=%u>\n",opt->sampleName1,virtchr[opt->chr].AC,iqr25,median,iqr75,covCntNb);
		printf("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth\">\n");
		printf("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
		printf("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Per-sample read depths for each allele\">\n");
		printf("##INFO=<ID=varType,Number=1,Type=String,Description=\"Type of variant\">\n");
		for (unsigned int c = 1; c <= chrMax; c += 1)
			printf("##contig=<ID=%s,length=%u,assembly=%s>\n",virtchr[c].AC,virtchr[c].len,opt->refName);
		printf("##phasing=none\n");
		printf("##source=SelectVariants\n");
		printf("##variants_justified=left\n");
		printf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n",opt->sampleName1);
		for (i = 0; i <= virtchr[opt->chr].len; i++)
		{
			if (cov[i].var != NULL && gNT2bits[genome_tbl[genomepos+i]] == kNOT_ACGT)
			{
				free(cov[i].var);
				cov[i].var = NULL;
			}
			if (cov[i].var != NULL)
			{
				int called = 0;
				unsigned int j;
				unsigned int total = 0;
				cov[i].var->ACGTP[gNT2bits[genome_tbl[genomepos+i]]] += cov[i].cntP;
				cov[i].var->ACGTM[gNT2bits[genome_tbl[genomepos+i]]] += cov[i].cntM;
				for (j = 0; j < 4; j++)
					total += cov[i].var->ACGTP[j] + cov[i].var->ACGTM[j];
				int max = -1;
				int min = -1;
				for (j = 0; j < 4; j++)
				{
					cov[i].var->ACGTP[j] += cov[i].var->ACGTM[j];
					if (max == -1 || cov[i].var->ACGTP[j] > cov[i].var->ACGTP[max])
						max = j;
					if (min == -1 || cov[i].var->ACGTP[j] < cov[i].var->ACGTP[min])
						min = j;
				}
				// report indels first, so that vcf output is ordered by pos
				for (j = 0; j < kMAXNBINDELS; j++)
				{
					if (cov[i].var->del[j].cnt == 0)
						break;
					if (cov[i].var->del[j].cnt >= opt->reportVCF && cov[i].var->del[j].cnt >= total * opt->minorAllelePct / 100)
					{
						char *gt = "1/1";
						if (cov[i].var->ACGTP[gNT2bits[genome_tbl[genomepos+i]]] >= opt->reportVCF && cov[i].var->ACGTP[gNT2bits[genome_tbl[genomepos+i]]] >= (total + cov[i].var->del[j].cnt) * opt->minorAllelePct / 100)
							gt = "0/1";
						called += 1;
						printf("%s\t%u\t.\t%.*s\t%c\t255\tPASS\t.\tGT:DP:AD\t%s:%u:%u,%u\n",
							virtchr[opt->chr].AC,
							i - 1,
							cov[i].var->del[j].len + 1,
							genome_tbl + genomepos + i - 1,
							genome_tbl[genomepos + i - 1],
							gt,
							cov[i].var->ACGTP[gNT2bits[genome_tbl[genomepos+i]]]+cov[i].var->del[j].cnt,
							cov[i].var->ACGTP[gNT2bits[genome_tbl[genomepos+i]]],
							cov[i].var->del[j].cnt);
					}
				}
				for (j = 0; j < kMAXNBINDELS; j++)
				{
					if (cov[i].var->ins[j].cnt == 0)
						break;
					if (cov[i].var->ins[j].cnt >= opt->reportVCF && cov[i].var->ins[j].cnt >= total * opt->minorAllelePct / 100)
					{
						char *gt = "1/1";
						// reads part of an insert do not count towards reference allele, even though the reference allele is "covered"...
						unsigned int refCnt = 0;
						if (cov[i].var->ACGTP[gNT2bits[genome_tbl[genomepos+i]]] >= cov[i].var->ins[j].cnt)
							refCnt = cov[i].var->ACGTP[gNT2bits[genome_tbl[genomepos+i]]] - cov[i].var->ins[j].cnt;
						if (refCnt >= opt->reportVCF && refCnt >= total * opt->minorAllelePct / 100)
							gt = "0/1";
						called += 1;
						printf("%s\t%u\t.\t%c\t%c%.*s\t255\tPASS\t.\tGT:DP:AD\t%s:%u:%u,%u\n",
							virtchr[opt->chr].AC,
							i - 1,
							genome_tbl[genomepos + i - 1],
							genome_tbl[genomepos + i - 1],
							38,
							cov[i].var->ins[j].nt,
							gt,
							refCnt + cov[i].var->ins[j].cnt,
							refCnt,
							cov[i].var->ins[j].cnt);
					}
				}
				if (cov[i].var->ACGTP[max] >= opt->reportVCF)
				{
					int a[4] = {max, -1, -1, min};
					for (j = 0; j < 4; j++)
						if (j != a[0] && j != a[3])
						{
							if (a[1] == -1)
								a[1] = j;
							else
								a[2] = j;
						}
					if (cov[i].var->ACGTP[a[1]] < cov[i].var->ACGTP[a[2]])
					{
						int tem = a[1];
						a[1] = a[2];
						a[2] = tem;
					}
					if (gNT2bits[genome_tbl[genomepos+i]] == max)
					{
						// major allele is reference; check if second is high enough and report
						if (cov[i].var->ACGTP[a[1]] >= opt->reportVCF && cov[i].var->ACGTP[a[1]] >= total * opt->minorAllelePct / 100)
						{
							called += 1;
							printf("%s\t%u\t.\t%c\t%c\t255\tPASS\t.\tGT:DP:AD\t0/1:%u:%u,%u\n",
								virtchr[opt->chr].AC,
								i,
								genome_tbl[genomepos+i],
								"ACGT"[a[1]],
								cov[i].var->ACGTP[max] + cov[i].var->ACGTP[a[1]],
								cov[i].var->ACGTP[max],
								cov[i].var->ACGTP[a[1]]);
						}
					}
					else
					{
						if (a[1] == gNT2bits[genome_tbl[genomepos+i]] || !(cov[i].var->ACGTP[a[1]] >= opt->reportVCF && cov[i].var->ACGTP[a[1]] >= total * opt->minorAllelePct / 100))
						{
							// major allele is alternate; report
							char *gt = "1/1";
							if (cov[i].var->ACGTP[a[1]] >= opt->reportVCF && cov[i].var->ACGTP[a[1]] >= total * opt->minorAllelePct / 100)
								gt = "0/1";
							called += 1;
							printf("%s\t%u\t.\t%c\t%c\t255\tPASS\t.\tGT:DP:AD\t%s:%u:%u,%u\n",
								virtchr[opt->chr].AC,
								i,
								genome_tbl[genomepos+i],
								"ACGT"[max],
								gt,
								cov[i].var->ACGTP[gNT2bits[genome_tbl[genomepos+i]]] + cov[i].var->ACGTP[max],
								cov[i].var->ACGTP[gNT2bits[genome_tbl[genomepos+i]]],
								cov[i].var->ACGTP[max]);
						}
						else
						{
							// reference is neither major nor minor
							called += 1;
							printf("%s\t%u\t.\t%c\t%c,%c\t255\tPASS\t.\tGT:DP:AD\t1/2:%u:%u,%u,%u\n",
								virtchr[opt->chr].AC,
								i,
								genome_tbl[genomepos+i],
								"ACGT"[max],
								"ACGT"[a[1]],
								cov[i].var->ACGTP[a[1]] + cov[i].var->ACGTP[max],
								cov[i].var->ACGTP[gNT2bits[genome_tbl[genomepos+i]]],
								cov[i].var->ACGTP[max],
								cov[i].var->ACGTP[a[1]]);
						}
					}
				}
				if (called == 0)
				{
					free(cov[i].var);
					cov[i].var = NULL;
				}
			}
		}
		for (i = 0; i <= virtchr[opt->chr].len; i++)
		{
			if (cov[i].var != NULL)
				reportRef(stderr,cov[i].var,genome_tbl[genomepos+i],ref,i,opt);
		}
		return;
	}

	printf("position\tref.nt\tcoverage\tA.cnt\tC.cnt\tG.cnt\tT.cnt\trank.1.idx\trank.2.idx\trank.3.idx\trank.4.idx\tflag\trank.1.cnt\trank.2.cnt\trank.3.cnt\trank.4.cnt\n");
	for (i = 0; i <= virtchr[opt->chr].len; i++)
	{
		if (cov[i].var != NULL && gNT2bits[genome_tbl[genomepos+i]] != kNOT_ACGT)
		{
			unsigned int j;
			unsigned int total = 0;
			cov[i].var->ACGTP[gNT2bits[genome_tbl[genomepos+i]]] += cov[i].cntP;
			cov[i].var->ACGTM[gNT2bits[genome_tbl[genomepos+i]]] += cov[i].cntM;
			for (j = 0; j < 4; j++)
				total += cov[i].var->ACGTP[j] + cov[i].var->ACGTM[j];
			int max = -1;
			int min = -1;
			for (j = 0; j < 4; j++)
			{
				cov[i].var->ACGTP[j] += cov[i].var->ACGTM[j];
				if (max == -1 || cov[i].var->ACGTP[j] > cov[i].var->ACGTP[max])
					max = j;
				if (min == -1 || cov[i].var->ACGTP[j] < cov[i].var->ACGTP[min])
					min = j;
			}
			// discard if highest count is less than 2
			if (cov[i].var->ACGTP[max] < 2)
				continue;
			int a[4] = {max, -1, -1, min};
			for (j = 0; j < 4; j++)
				if (j != a[0] && j != a[3])
				{
					if (a[1] == -1)
						a[1] = j;
					else
						a[2] = j;
				}
			if (cov[i].var->ACGTP[a[1]] < cov[i].var->ACGTP[a[2]])
			{
				int tem = a[1];
				a[1] = a[2];
				a[2] = tem;
			}
			printf("%u\t%c\t%u",i,genome_tbl[genomepos+i],total);
			for (j = 0; j < 4; j++)
				printf("\t%u",cov[i].var->ACGTP[j]);
			for (j = 0; j < 4; j++)
				printf("\t%u",a[j]);
			printf("\t%u",cov[i].var->ACGTP[a[1]] > 0 && cov[i].var->ACGTP[a[1]] == cov[i].var->ACGTP[a[2]] ? 1 : 0);
			for (j = 0; j < 4; j++)
				printf("\t%u",cov[i].var->ACGTP[a[j]]);
			printf("\n");
		}
	}
} // report
//---------------------------------------------------------------

static void
myReport(OPTIONS *opt, coverage_p_t cov, reference_p_t ref)
{
	unsigned int genomepos = (virtchr[opt->chr].chr << 28) + virtchr[opt->chr].offset;
	unsigned int i;
	if (opt->PWMcnt > 0)
	{
		for (i = 0; i < opt->PWMcnt; i++)
			outputPWM(opt, i, cov, genomepos);
		return;
	}

	if (opt->wiggleBed[0] != 0)
		printWiggle(cov,opt);

	if (opt->gencodeGFF3[0] != 0)
		printExpr(cov,ref,opt);

	if (opt->tromerTranscripts[0] != 0)
		printTromer(cov,ref,opt);

	// compute median and iqr
	unsigned int covCntNb, median, iqr25, iqr75;
	unsigned int covCntNbSel, medianSel, iqr25Sel, iqr75Sel;
	computeCov(cov, opt, NULL, &covCntNb, &median, &iqr25, &iqr75);
	if (opt->selectRegions[0] != 0)
		computeCov(cov, opt, cov, &covCntNbSel, &medianSel, &iqr25Sel, &iqr75Sel);

	if (opt->reportVCF > 0)
	{
		printf("##fileformat=VCFv4.1\n");
		if (opt->selectRegions[0] != 0)
			printf("##SAMPLE=<ID=%s,CONTIG=%s,IQR25=%u,MEDIAN=%u,IQR75=%u,CNT=%u,SEL=%s>\n",opt->sampleName1,virtchr[opt->chr].AC,iqr25Sel,medianSel,iqr75Sel,covCntNbSel,opt->selectRegions);
		printf("##SAMPLE=<ID=%s,CONTIG=%s,IQR25=%u,MEDIAN=%u,IQR75=%u,CNT=%u>\n",opt->sampleName1,virtchr[opt->chr].AC,iqr25,median,iqr75,covCntNb);
		printf("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth\">\n");
		printf("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
		printf("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Per-sample read depths for each allele\">\n");
		printf("##INFO=<ID=varType,Number=1,Type=String,Description=\"Type of variant\">\n");
		for (unsigned int c = 1; c <= chrMax; c += 1)
			printf("##contig=<ID=%s,length=%u,assembly=%s>\n",virtchr[c].AC,virtchr[c].len,opt->refName);
		printf("##phasing=none\n");
		printf("##source=SelectVariants\n");
		printf("##variants_justified=left\n");
		printf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n",opt->sampleName1);
		for (i = 0; i <= virtchr[opt->chr].len; i++)
		{
			if (cov[i].var != NULL && gNT2bits[genome_tbl[genomepos+i]] == kNOT_ACGT)
			{
				free(cov[i].var);
				cov[i].var = NULL;
			}
			if (cov[i].var != NULL)
			{
				int called = 0;
				unsigned int j;

				cov[i].var->ACGTP[gNT2bits[genome_tbl[genomepos+i]]] += cov[i].cntP;
				cov[i].var->ACGTM[gNT2bits[genome_tbl[genomepos+i]]] += cov[i].cntM;
				int max = -1;
				int min = -1;
				for (j = 0; j < 4; j++)
				{
					cov[i].var->ACGTP[j] += cov[i].var->ACGTM[j];
					if (max == -1 || cov[i].var->ACGTP[j] > cov[i].var->ACGTP[max])
						max = j;
					if (min == -1 || cov[i].var->ACGTP[j] < cov[i].var->ACGTP[min])
						min = j;
				}

				// Compute coverage (total)
				unsigned int total = 0U;
				for (j = 0; j < 4; j++) total += cov[i].var->ACGTP[j];
				for (j = 0; j < kMAXNBINDELS; j++) {
					if (cov[i].var->del[j].cnt == 0) break;
					total += cov[i].var->del[j].cnt;
				}
				for (j = 0; j < kMAXNBINDELS; j++)
				{
					if (cov[i].var->ins[j].cnt == 0) break;
					total += cov[i].var->ins[j].cnt;
				}

				const unsigned int MinorAlleleCoverage = total * opt->minorAllelePct / 100;

				// report indels first, so that vcf output is ordered by pos
				for (j = 0; j < kMAXNBINDELS; j++)
				{
					if (cov[i].var->del[j].cnt == 0)
						break;
					if (cov[i].var->del[j].cnt >= opt->reportVCF && cov[i].var->del[j].cnt >= MinorAlleleCoverage)
					{
						char *gt = "1/1";
						const unsigned int refCnt = cov[i].var->ACGTP[gNT2bits[genome_tbl[genomepos+i]]];
						if (refCnt >= opt->reportVCF && refCnt >= MinorAlleleCoverage)
							gt = "0/1";
						called += 1;
						printf("%s\t%u\t.\t%.*s\t%c\t255\tPASS\t.\tGT:DP:AD\t%s:%u:%u,%u\n",
							virtchr[opt->chr].AC,
							i - 1,
							cov[i].var->del[j].len + 1,
							genome_tbl + genomepos + i - 1,
							genome_tbl[genomepos + i - 1],
							gt,
							refCnt+cov[i].var->del[j].cnt,
							refCnt,
							cov[i].var->del[j].cnt);
					}
				}
				for (j = 0; j < kMAXNBINDELS; j++)
				{
					if (cov[i].var->ins[j].cnt == 0)
						break;
					if (cov[i].var->ins[j].cnt >= opt->reportVCF && cov[i].var->ins[j].cnt >= MinorAlleleCoverage)
					{
						char *gt = "1/1";
						// reads part of an insert do not count towards reference allele, even though the reference allele is "covered"...
						const unsigned int refCnt = cov[i].var->ACGTP[gNT2bits[genome_tbl[genomepos+i]]];
						if (refCnt >= opt->reportVCF && refCnt >= MinorAlleleCoverage)
							gt = "0/1";
						called += 1;
						printf("%s\t%u\t.\t%c\t%c%.*s\t255\tPASS\t.\tGT:DP:AD\t%s:%u:%u,%u\n",
							virtchr[opt->chr].AC,
							i - 1,
							genome_tbl[genomepos + i - 1],
							genome_tbl[genomepos + i - 1],
							38,
							cov[i].var->ins[j].nt,
							gt,
							refCnt + cov[i].var->ins[j].cnt,
							refCnt,
							cov[i].var->ins[j].cnt);
					}
				}
				if (cov[i].var->ACGTP[max] >= opt->reportVCF)
				{
					int a[4] = {max, -1, -1, min};
					for (j = 0; j < 4; j++)
						if (j != a[0] && j != a[3])
						{
							if (a[1] == -1)
								a[1] = j;
							else
								a[2] = j;
						}
					if (cov[i].var->ACGTP[a[1]] < cov[i].var->ACGTP[a[2]])
					{
						int tem = a[1];
						a[1] = a[2];
						a[2] = tem;
					}
					if (gNT2bits[genome_tbl[genomepos+i]] == max)
					{
						// major allele is reference; check if second is high enough and report
						if (cov[i].var->ACGTP[a[1]] >= opt->reportVCF && cov[i].var->ACGTP[a[1]] >= MinorAlleleCoverage)
						{
							called += 1;
							printf("%s\t%u\t.\t%c\t%c\t255\tPASS\t.\tGT:DP:AD\t0/1:%u:%u,%u\n",
								virtchr[opt->chr].AC,
								i,
								genome_tbl[genomepos+i],
								"ACGT"[a[1]],
								cov[i].var->ACGTP[max] + cov[i].var->ACGTP[a[1]],
								cov[i].var->ACGTP[max],
								cov[i].var->ACGTP[a[1]]);
						}
					}
					else
					{
						if (a[1] == gNT2bits[genome_tbl[genomepos+i]] || !(cov[i].var->ACGTP[a[1]] >= opt->reportVCF && cov[i].var->ACGTP[a[1]] >= MinorAlleleCoverage))
						{
							// major allele is alternate; report
							char *gt = "1/1";
							if (cov[i].var->ACGTP[a[1]] >= opt->reportVCF && cov[i].var->ACGTP[a[1]] >= MinorAlleleCoverage)
								gt = "0/1";
							called += 1;
							printf("%s\t%u\t.\t%c\t%c\t255\tPASS\t.\tGT:DP:AD\t%s:%u:%u,%u\n",
								virtchr[opt->chr].AC,
								i,
								genome_tbl[genomepos+i],
								"ACGT"[max],
								gt,
								cov[i].var->ACGTP[gNT2bits[genome_tbl[genomepos+i]]] + cov[i].var->ACGTP[max],
								cov[i].var->ACGTP[gNT2bits[genome_tbl[genomepos+i]]],
								cov[i].var->ACGTP[max]);
						}
						else
						{
							// reference is neither major nor minor
							called += 1;
							printf("%s\t%u\t.\t%c\t%c,%c\t255\tPASS\t.\tGT:DP:AD\t1/2:%u:%u,%u,%u\n",
								virtchr[opt->chr].AC,
								i,
								genome_tbl[genomepos+i],
								"ACGT"[max],
								"ACGT"[a[1]],
								cov[i].var->ACGTP[a[1]] + cov[i].var->ACGTP[max],
								cov[i].var->ACGTP[gNT2bits[genome_tbl[genomepos+i]]],
								cov[i].var->ACGTP[max],
								cov[i].var->ACGTP[a[1]]);
						}
					}
				}
				if (called == 0)
				{
					free(cov[i].var);
					cov[i].var = NULL;
				}
			}
		}
		for (i = 0; i <= virtchr[opt->chr].len; i++)
		{
			if (cov[i].var != NULL)
				reportRef(stderr,cov[i].var,genome_tbl[genomepos+i],ref,i,opt);
		}
		return;
	}

	printf("position\tref.nt\tcoverage\tA.cnt\tC.cnt\tG.cnt\tT.cnt\trank.1.idx\trank.2.idx\trank.3.idx\trank.4.idx\tflag\trank.1.cnt\trank.2.cnt\trank.3.cnt\trank.4.cnt\n");
	for (i = 0; i <= virtchr[opt->chr].len; i++)
	{
		if (cov[i].var != NULL && gNT2bits[genome_tbl[genomepos+i]] != kNOT_ACGT)
		{
			unsigned int j;
			unsigned int total = 0;
			cov[i].var->ACGTP[gNT2bits[genome_tbl[genomepos+i]]] += cov[i].cntP;
			cov[i].var->ACGTM[gNT2bits[genome_tbl[genomepos+i]]] += cov[i].cntM;
			for (j = 0; j < 4; j++)
				total += cov[i].var->ACGTP[j] + cov[i].var->ACGTM[j];
			int max = -1;
			int min = -1;
			for (j = 0; j < 4; j++)
			{
				cov[i].var->ACGTP[j] += cov[i].var->ACGTM[j];
				if (max == -1 || cov[i].var->ACGTP[j] > cov[i].var->ACGTP[max])
					max = j;
				if (min == -1 || cov[i].var->ACGTP[j] < cov[i].var->ACGTP[min])
					min = j;
			}
			// discard if highest count is less than 2
			if (cov[i].var->ACGTP[max] < 2)
				continue;
			int a[4] = {max, -1, -1, min};
			for (j = 0; j < 4; j++)
				if (j != a[0] && j != a[3])
				{
					if (a[1] == -1)
						a[1] = j;
					else
						a[2] = j;
				}
			if (cov[i].var->ACGTP[a[1]] < cov[i].var->ACGTP[a[2]])
			{
				int tem = a[1];
				a[1] = a[2];
				a[2] = tem;
			}
			printf("%u\t%c\t%u",i,genome_tbl[genomepos+i],total);
			for (j = 0; j < 4; j++)
				printf("\t%u",cov[i].var->ACGTP[j]);
			for (j = 0; j < 4; j++)
				printf("\t%u",a[j]);
			printf("\t%u",cov[i].var->ACGTP[a[1]] > 0 && cov[i].var->ACGTP[a[1]] == cov[i].var->ACGTP[a[2]] ? 1 : 0);
			for (j = 0; j < 4; j++)
				printf("\t%u",cov[i].var->ACGTP[a[j]]);
			printf("\n");
		}
	}
} // report
//---------------------------------------------------------------

static void
analyzeVCF(OPTIONS *opt, coverage_p_t cov, reference_p_t ref, const char * const restrict VCF)
{

	inline unsigned int getACGTIndex(const char val) {
		unsigned int index;
		switch (val) {
			case 'A': index = 0; break;
			case 'C': index = 1; break;
			case 'G': index = 2; break;
			case 'T': index = 3; break;
#ifndef NDEBUG
			default:
				fprintf(stderr, "Unknown character %c given to getACGTIndex\n", (int) val);
				exit(1);
#endif
		}

		return index;
	}

	const unsigned int genomepos = (virtchr[opt->chr].chr << 28) + virtchr[opt->chr].offset;

	////////////////////////////////////////////////////////////////////////////////////////////
	// Open VCF file
	htsFile * inf = bcf_open(VCF, "r");
	if (inf == NULL) {
		exit(1);
	}

	////////////////////////////////////////////////////////////////////////////////////////////
	// Read header and sequences' name
	bcf_hdr_t *hdr = bcf_hdr_read(inf);
	int nseq;
	const char **seqnames = bcf_hdr_seqnames(hdr, &nseq);
	if (seqnames == NULL) {
		fputs("Unable to get sequence names in VCF file\n", stderr);
		exit(1);
	}
	int WantedSequenceId = -1;
	{
		char tID[3];
		snprintf(tID,3, "%i", opt->chr);
		for (int i = 0; i < nseq; i++) {
			if (strcmp(tID, seqnames[i]) == 0) {
				WantedSequenceId = i;
				break;
			}
		}
	}

	if (WantedSequenceId == -1) {
		fprintf(stderr, "Unable to get chromosome %i sequence in VCF file\n", opt->chr);
		exit(1);
	}

	////////////////////////////////////////////////////////////////////////////////////////////
	// Init record
	bcf1_t *rec = bcf_init();
	if (rec == NULL) {
		fputs("Unable to initialize record\n", stderr);
		exit(1);
	}

	////////////////////////////////////////////////////////////////////////////////////////////
	// Loop on VCF file

	// coverage data for each call
	int ndp_arr = 0;
	int ndp     = 0;
	int *dp     = NULL;
	// genotype data for each call
	// genotype arrays are twice as large as
	// the other arrays as there are two values for each sample
	int ngt_arr = 0;
	int ngt     = 0;
	int *gt     = NULL;

	while (bcf_read(inf, hdr, rec) == 0) {
		if (rec->rid == WantedSequenceId) {
			////////////////////////////////////////////////////////////////////////////////////////
			// Get VCF data
			ndp = bcf_get_format_int32(hdr, rec, "DP", &dp, &ndp_arr);
			ngt = bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt_arr);
			static const char *types[] = {"SNP", "DELETION", "INSERTION" };
			unsigned int CurrentTypeID = 0;

			assert(ngt == 2);

			bcf_get_variant_types(rec);
			if (!bcf_is_snp(rec)) {
				if (rec->d.var[bcf_gt_allele(gt[0])].n < 0 || rec->d.var[bcf_gt_allele(gt[1])].n < 0) {
					CurrentTypeID = 1;
				}
				else {
					CurrentTypeID = 2;
				}
			}

			const int depth = (ndp > 0 ) ? *dp : -1;

			////////////////////////////////////////////////////////////////////////////////////////
			// Get the genome location data
			/* Position is a nightmare !!! */
			unsigned int RealMemoryPosition = 1 + rec->pos;
			if (CurrentTypeID != 0) RealMemoryPosition += 1; /* God knows why, this is necessary !!! */

			coverage_t * restrict Coverage = cov + RealMemoryPosition;
			gtl_variant_t * const restrict var = Coverage->var;

			////////////////////////////////////////////////////////////////////////////////////////
			// Compute coverage
			unsigned int nSeenTags = Coverage->cntP + Coverage->cntM;
			if (var) {
				for (int j=0; j<4; j++) nSeenTags += var->ACGTP[j];
				for (int j=0; j<4; j++) nSeenTags += var->ACGTM[j];
				for (int j = 0; j < kMAXNBINDELS; j++) {
					if (var->del[j].cnt == 0) break;
					nSeenTags += var->del[j].cnt;
				}
				for (int j = 0; j < kMAXNBINDELS; j++)
				{
					if (var->ins[j].cnt == 0) break;
					nSeenTags += var->ins[j].cnt;
				}
			}

			////////////////////////////////////////////////////////////////////////////////////////
			// Output informations
			printf("%u\t%u\t%s\t%s\t%s/%s\t%u/%u\t%i\t%c\t", opt->chr,
						 RealMemoryPosition, types[CurrentTypeID], rec->d.allele[0],
						 rec->d.allele[bcf_gt_allele(gt[0])], rec->d.allele[bcf_gt_allele(gt[1])],
						 bcf_gt_allele(gt[0]), bcf_gt_allele(gt[1]), depth,
						 (CurrentTypeID != 0) ? (int) genome_tbl[genomepos+RealMemoryPosition-1]
						                      : (int) genome_tbl[genomepos+RealMemoryPosition]);

			if (nSeenTags != 0) {
				printf("%u\t{%u,%u}", nSeenTags, Coverage->cntP,Coverage->cntM);
				if (var) {
					const unsigned int MinorAlleleCoverage = nSeenTags * opt->minorAllelePct / 100;
					////////////////////////////////////////////////////////////////////////////////////
					// PRINT WHAT WE HAVE
					var->ACGTP[gNT2bits[genome_tbl[genomepos+RealMemoryPosition - (CurrentTypeID != 0)]]] += Coverage->cntP;
					var->ACGTM[gNT2bits[genome_tbl[genomepos+RealMemoryPosition - (CurrentTypeID != 0)]]] += Coverage->cntM;
					printf("\t{[%u,%u,%u,%u],", var->ACGTP[0], var->ACGTP[1], var->ACGTP[2], var->ACGTP[3]);
					printf("[%u,%u,%u,%u]}", var->ACGTM[0], var->ACGTM[1], var->ACGTM[2], var->ACGTM[3]);
					for (int j=0; j<4; j++) var->ACGTP[j] += var->ACGTM[j];

					fputs("\t{", stdout);
					const unsigned char * const restrict CurrentGenomePos = &genome_tbl[genomepos+RealMemoryPosition];
					for (int j=0; j<kMAXNBINDELS; j++) {
						if (var->del[j].cnt == 0) break;
						printf("(%u,%.*s)", var->del[j].cnt, (int) var->del[j].len, CurrentGenomePos);
					}
					fputs("}\t{", stdout);
					for (int j=0; j<kMAXNBINDELS; j++) {
						if (var->ins[j].cnt == 0) break;
						printf("(%u,%.*s)", var->ins[j].cnt, 38, var->ins[j].nt);
					}
					fputc('}', stdout);

					int FoundError = 0;
					////////////////////////////////////////////////////////////////////////////////////
					// THIS IS SUPPOSED TO BE A SNP
					if (CurrentTypeID == 0) {
						/* Get the most common DNA base */
// 						unsigned int BestBaseId = 0;
// 						for (unsigned int j = 1; j < 4; j++)
// 						{
// 							if (var->ACGTP[j] > var->ACGTP[BestBaseId]) BestBaseId = j;
// 						}
						/* GT = x/y
						 * Testing y first
						 */
						const unsigned int index = getACGTIndex(rec->d.allele[bcf_gt_allele(gt[1])][0]);
// 						if (index != BestBaseId) {
// 							fputs("\tBAD_GT_2", stdout);
// 						}
						if (var->ACGTP[index] < MinorAlleleCoverage) {
							printf("\tMAF_THRESHOLD_GT_2(%f)", (float) var->ACGTP[index]/(float) nSeenTags);
							FoundError = 1;
						}

						if (gt[1] != gt[0]) {
							/* Testing x */
							const unsigned int index2 = getACGTIndex(rec->d.allele[bcf_gt_allele(gt[0])][0]);
							if (var->ACGTP[index2] < MinorAlleleCoverage) {
								printf("\tMAF_THRESHOLD_GT_1(%f)", (float) var->ACGTP[index2]/(float) nSeenTags);
								FoundError = 1;
							}
						}
					}
					////////////////////////////////////////////////////////////////////////////////////
					// THIS IS SUPPOSED TO BE A DELETION
					else if ( CurrentTypeID == 1) {
						if (rec->d.var[bcf_gt_allele(gt[1])].type != VCF_INDEL) {
							const unsigned int index = getACGTIndex(rec->d.allele[bcf_gt_allele(gt[1])][0]);
							if (var->ACGTP[index] < MinorAlleleCoverage) {
								printf("\tMAF_THRESHOLD_GT_2(%f)", (float) var->ACGTP[index]/(float) nSeenTags);
								FoundError = 1;
							}
						}
						else {
							int index = -1;
							const unsigned short DelLength = -rec->d.var[bcf_gt_allele(gt[1])].n;
							for (unsigned int j=0; j<kMAXNBINDELS; j++) {
								if (var->del[j].cnt == 0) break;
								if (var->del[j].len == DelLength) index = j;
							}
							if (index == -1) {
								fputs("\tDELETION_NOT_FOUND_GT_2", stdout);
								FoundError = 1;
							}
							else if (var->del[index].cnt < MinorAlleleCoverage){
								printf("\tMAF_THRESHOLD_GT_2(%f)", (float) var->del[index].cnt/(float) nSeenTags);
								FoundError = 1;
							}
						}

						if (gt[1] != gt[0]) {
							if (rec->d.var[bcf_gt_allele(gt[0])].type != VCF_INDEL) {
								const unsigned int index = getACGTIndex(rec->d.allele[bcf_gt_allele(gt[0])][0]);
								if (var->ACGTP[index] < MinorAlleleCoverage) {
									printf("\tMAF_THRESHOLD_GT_1(%f)", (float) var->ACGTP[index]/(float) nSeenTags);
									FoundError = 1;
								}
							}
							else {
								int index = -1;
								const unsigned short DelLength = -rec->d.var[bcf_gt_allele(gt[0])].n;
								for (unsigned int j=0; j<kMAXNBINDELS; j++) {
									if (var->del[j].cnt == 0) break;
									if (var->del[j].len == DelLength) index = j;
								}
								if (index == -1) {
									fputs("\tDELETION_NOT_FOUND_GT_1", stdout);
									FoundError = 1;
								}
								else if (var->del[index].cnt < MinorAlleleCoverage){
									printf("\tMAF_THRESHOLD_GT_1(%f)", (float) var->del[index].cnt/(float) nSeenTags);
									FoundError = 1;
								}
							}
						}
					}
					////////////////////////////////////////////////////////////////////////////////////
					// THIS IS SUPPOSED TO BE AN INSERTION
					else {
						if (rec->d.var[bcf_gt_allele(gt[1])].type != VCF_INDEL) {
							const unsigned int index = getACGTIndex(rec->d.allele[bcf_gt_allele(gt[1])][0]);
							if (var->ACGTP[index] < MinorAlleleCoverage) {
								printf("\tMAF_THRESHOLD_GT_2(%f)", (float) var->ACGTP[index]/(float) nSeenTags);
								FoundError = 1;
							}
						}
						else {
							int index = -1;
							for (unsigned int j=0; j<kMAXNBINDELS; j++) {
								if (var->ins[j].cnt == 0) break;
								if (strcmp(var->ins[j].nt, &(rec->d.allele[bcf_gt_allele(gt[1])][1])) == 0) index = j;
							}
							if (index == -1) {
								fputs("\tINSERTION_NOT_FOUND_GT_2", stdout);
								FoundError = 1;
							}
							else if (var->ins[index].cnt < MinorAlleleCoverage){
								printf("\tMAF_THRESHOLD_GT_2(%f)", (float) var->ins[index].cnt/(float) nSeenTags);
								FoundError = 1;
							}
						}
						if (gt[1] != gt[0]) {
							if (rec->d.var[bcf_gt_allele(gt[0])].type != VCF_INDEL) {
								const unsigned int index = getACGTIndex(rec->d.allele[bcf_gt_allele(gt[0])][0]);
								if (var->ACGTP[index] < MinorAlleleCoverage) {
									printf("\tMAF_THRESHOLD_GT_1(%f)", (float) var->ACGTP[index]/(float) nSeenTags);
									FoundError = 1;
								}
							}
							else {
								int index = -1;
							for (unsigned int j=0; j<kMAXNBINDELS; j++) {
									if (var->ins[j].cnt == 0) break;
									if (strcmp(var->ins[j].nt, &(rec->d.allele[bcf_gt_allele(gt[0])][1])) == 0) index = j;
								}
								if (index == -1) {
									fputs("\tINSERTION_NOT_FOUND_GT_1", stdout);
									FoundError = 1;
								}
								else if (var->del[index].cnt < MinorAlleleCoverage){
									printf("\tMAF_THRESHOLD_GT_1(%f)", (float) var->del[index].cnt/(float) nSeenTags);
									FoundError = 1;
								}
							}
						}
					}
					/* If no error was found, we probably got an extra variant from the reference */
					if ( FoundError == 0 ) {
						fputs("\tEXTRA_VARIANT", stdout);
					}

					fputc('\n', stdout);
				}
				else {
					printf("\tVARIANT_NOT_FOUND\n");
				}
			}
			else {
				printf("NO_COVERAGE\n");
			}
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////
	// Free memory
	free(gt);
	free(dp);
	free(seqnames);
	bcf_hdr_destroy(hdr);
	bcf_destroy(rec);

	////////////////////////////////////////////////////////////////////////////////////////////
	// Close VCF file
	bcf_close(inf);
}
//---------------------------------------------------------------

static void
loadPos(const char *fname, OPTIONS *opt)
{
	char linbuf[kMaxLineBuf];
	FILE *f = fopen (fname, "r");
	if (!f)
	{
		fprintf(stderr, "Error:Cannot open %s\n",fname);
		exit(1);
	}
	fgets(linbuf,kMaxLineBuf,f); // skip header.
	while(!feof(f))
	{
		char tmpchr[32];
		fgets(linbuf,kMaxLineBuf,f);
		if (feof(f))
			break;
		sscanf(&linbuf[0],"%s\t%u",tmpchr,&opt->PWMpos[opt->PWMcnt].pos);
		if ((tmpchr[0] == 'c') && (tmpchr[1] == 'h') && (tmpchr[2] == 'r'))
		{
			if (tmpchr[3] == 'X') opt->PWMpos[opt->PWMcnt].chr = 23;
			else if (tmpchr[3] == 'Y') opt->PWMpos[opt->PWMcnt].chr = 24;
			else opt->PWMpos[opt->PWMcnt].chr = atoi(&tmpchr[3]);
		}
		else opt->PWMpos[opt->PWMcnt].chr = atoi(tmpchr);
		if (++opt->PWMcnt == kMAXPWMPOS)
			break;
	}
	fclose(f);
}
//---------------------------------------------------------------

int
main (int argc, char **argv)
{
	static char default_genome[] = "/tmp/nodelete/nguex/data/hg19.bin";
	OPTIONS options;
	const char * restrict rsltfile = NULL;
	const char * restrict genomefile = default_genome;
	const char * restrict AlignerScoreConfFile = NULL;
	const char * restrict AnalysisVCFFile = NULL;
	int c, err;

	memset(&options,0,sizeof(OPTIONS));
	options.maxPairLen = kMAXPAIRLEN;
	strcpy(options.outbase,"/tmp/");

	/* --------- process arguments ----------------------*/

	opterr = 0;
	while ((c = getopt (argc, argv, "r:g:C:v:dp:o:acf:V:m:N:s:w:b:G:e:t:F:D:T:A:W:")) != -1)
	switch (c)
	{
		case 'p':
			loadPos(optarg,&options);
			break;

		case 'o':
			strcpy(options.outbase,optarg);
			break;

		case 'r':
			rsltfile = optarg;
			break;

		case 'b':
			strcpy(options.wiggleBed,optarg);
			break;

		case 'G':
			strcpy(options.gencodeGFF3,optarg);
			break;

		case 't':
			strcpy(options.tromerTranscripts,optarg);
			break;

		case 'e':
			strcpy(options.expressionBase,optarg);
			break;

		case 's':
			strcpy(options.selectRegions,optarg);
			break;

		case 'w':
			options.selectAdd = atoi(optarg);
			break;

		case 'C':
			sscanf(optarg,"%d",&options.chr);
			break;

		case 'g':
			genomefile = optarg;
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

		case 'V':
			options.reportVCF = atoi(optarg);
			break;

		case 'm':
			options.minorAllelePct = atoi(optarg);
			break;

		case 'N':
			strncpy(options.sampleName1,optarg,64);
			break;

		case 'v':
			sscanf(optarg,"%d",&verbose);
			break;

		case 'F':
			AlignerScoreConfFile = optarg;
			if (loadAlignerScores(AlignerScoreConfFile) < 0) {
				fputs("Aligner configuration file failed to load\n", stderr);
				return 1;
			}
			options.Realign = 1;
			break;

		case 'D':
			dumpAlignerScores(optarg);
			return 0;
			break;

		case 'T':
		{
			int dummy = atoi(optarg);
			if (dummy > 0)
				nThreads = (unsigned int) dummy;
			else {
				fputs("Thread number must be positive\n", stderr);
				exit(1);
			}
		}
		break;

		case 'A':
			AnalysisVCFFile = optarg;
			break;

		case 'W':
			options.maxPairLen = atoi(optarg);
			break;
	}

	if (rsltfile == 0 || options.chr == -1)
	{
		printf("usage:\n\n");
		printf("GTLcaller [options] -C chromosome -r resultfile\n\n");
		printf("           -g genomefile         : full path of the genome encoded by encode_genome\n");
		printf("                                   default: '%s' \n",genomefile);
		printf("           -C chromosome         : select chromosome to analyze\n");
		printf("           -r resultfile         : name of the result file (produced by match)\n");
		printf("           -a                    : perform additional filtering, discard pairs with more than 2 indel or indel longer than 40 nt (default is off)\n");
		printf("           -c                    : do not filter clones (default is to filter them out)\n");
		printf("           -f <char>             : filter out nucleotides where the quality value is *below* the specified character [%c]\n", options.filterQualValue + '#');
		printf("           -d                    : debug mode...\n");
		printf("           -V <int>              : report as VCF, insertions, deletions and SNPs; value is minimum observation count\n");
		printf("           -m <int>              : minor allele minimum percent coverage\n");
		printf("           -N sampleName         : sample name to put in VCF file\n");
		printf("           -o <path>             : base part of output files for PWM [%s]\n",options.outbase);
		printf("           -p <file>             : generate PWM output (for weblogo) based on list of chr:pos found in supplied file\n");
		printf("           -s <file>             : selected regions to compute coverage (bed-style)\n");
		printf("           -w <int>              : augment selected regions 5' and 3' by this number of nt [%u]\n",options.selectAdd);
		printf("           -b <file>             : generate bed/wiggle file with covered regions\n");
		printf("           -G <file>             : read gene and exon definitions from this file, use with -e option\n");
		printf("           -t <file>             : read tromer transcript definitions from this file, use with -e option\n");
		printf("           -e <file>             : generate gene expression and exon coverage files (expected to end in _gene.txt and _exon.txt\n");
		printf("           -v level              : verbose level\n");
		printf("           -F alignconffile      : aligner configuration file\n");
		printf("           -D alignconffile      : dump alignment score JSON file\n");
		printf("           -T <uint>             : number of thread to dispatch work on\n");
		printf("           -A <vcf file>         : analyze and report position given in vcf file\n");
		printf("           -W <int>              : maximum length of 'well-behaved' pairs (distance and relative strand-orientation); 0 means no check [%u]\n", options.maxPairLen);
		printf("\n");
		printf(GTL_VERSION_FULL_STRING "\n");
		printf("(c) Nicolas Guex, Christian Iseli & Cie 2014-2017\n");
		return(1);
	}

	if (options.sampleName1[0] == 0)
		snprintf(options.sampleName1, 64, "%s", rsltfile);

	char *configfn = strdup(genomefile);
	const size_t len = strlen(configfn);
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

	coverage_p_t cov = calloc(virtchr[options.chr].len + kMaxReadLen + 1,sizeof(coverage_t));
	// FIXME - maybe someday look into circular DNA, but for now just allocate some extra space at the end of the piece
	reference_p_t ref = calloc((virtchr[options.chr].len + 1025) >> 8,sizeof(reference_t));
	loadSelected(cov,&options);

	if (!AlignerScoreConfFile) {
		if (nThreads == 0U) {
			if ((err = decompress(rsltfile,&options,cov,ref)) != 0) goto bail;
		}
		else {
			if ((err = decompress_PairedThreaded(rsltfile,&options,cov,ref)) != 0) goto bail;
		}
	}
	else if (AlignerScoreConfFile && nThreads == 0U) {
		if ((err = decompress_realign(rsltfile,&options,cov,ref)) != 0) goto bail;
	}
	else {
		if ((err = decompress_PairedThreaded(rsltfile,&options,cov,ref)) != 0) goto bail;
	}

	if (AnalysisVCFFile == NULL) {
		if (!AlignerScoreConfFile)
			report(&options,cov,ref);
		else
			myReport(&options,cov,ref);
	}
	else {
		analyzeVCF(&options, cov, ref, AnalysisVCFFile);
	}
	free(cov);
	free(ref);

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
