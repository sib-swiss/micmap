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

	common routines used by various callers

	(c) N.Guex and C.Iseli 2017

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
#include <pthread.h>
#include <signal.h> 
#include <semaphore.h>
#include <sched.h> 
#include <assert.h>
#include <errno.h>

// basic compression
#include <zlib.h>

#include "virt_chr.h"
#include "GTL.h"
#include "GTLcaller_common.h"

//#define BYPASS_COMMON
//#define ORDINAL_TO_PRINT 532719543
//#define POSITION_TO_PRINT 2295852

//=================================================================================================

unsigned char *genome_tbl;	// contain genome in 16 virtual chromosomes.
unsigned char gNT2bits[128];
VIRTUALCHR virtchr[kMAXCHRcnt];
unsigned int nThreads = 0U;

//=================================================================================================
static int
bt_node_compare(const void *a, const void *b)
{
  const bt_node_p_t na = (bt_node_p_t) a, nb = (bt_node_p_t) b;
	if (na->refNo < nb->refNo)
		return -1;
	if (na->refNo > nb->refNo)
		return 1;
	return 0;
}
//---------------------------------------------------------------

int
loadGenome(const char * const restrict genomefile, OPTIONS *opt)
{
	int err = 1;
	int sf = 0;
	struct stat genomeStat;

		/* ----------------------- load genome data ----------------------*/

		sf = open (genomefile, O_RDONLY);
		if (!sf)
		{
			printf("can't open bin file\n");
			goto bail;
		}

		if (fstat (sf,&genomeStat) != 0)
		{
			printf("cannot stat file %s\n", genomefile);
			goto bail;
		}
		if (genomeStat.st_size != kGENOME_DATA_SIZE)
		{
			printf("genomeStat.st_size %lx != %lx\n",genomeStat.st_size,kGENOME_DATA_SIZE);
			goto bail;
		}

		genome_tbl = (unsigned char *)mmap (genome_tbl, kGENOME_DATA_SIZE, PROT_READ, MAP_PRIVATE, sf, 0);
		if (genome_tbl == MAP_FAILED)
		{
			printf("mapping failed.\n");
			goto bail;
		}
		if (opt->refName[0] == 0)
		{
			// get reference assembly name
			char *s = strrchr(genomefile,'/');
			if (s != NULL)
			{
				s += 1;
				unsigned int l = strlen(s);
				if (l > 4)
					l -= 4;
				if (l >= 64)
					l = 63;
				strncpy(opt->refName,s,l);
				opt->refName[63] = 0;
			}
		}

		err = 0;

	bail:;

		return(err);

} // loadGenome
//---------------------------------------------------------------

static void
decodeSingle(TLBDATA *tpp, unsigned int i, decodedSingle_p_t dsp)
{
	if (tpp->header.flags_7_0 & kTLBHflagFixedLength)
	{
		dsp->taglen1 = tpp->lengthInfo;
	}
	else
	{
		dsp->taglen1 = tpp->len[i];
	}
	char *next = strchr(tpp->hdr, '\n');
	tpp->hdr = next + 1;
	dsp->ordinal = 0;
	for (unsigned int j = 0; j < 8; j++)
		dsp->ordinal = (dsp->ordinal << 8) | (unsigned long) *((unsigned char *) tpp->hdr++);

	if (tpp->ppd[0][i].tag2offset & k32bREVCOMP_TAG1_BIT)
		dsp->reverseTAG1 = 1;
	else
		dsp->reverseTAG1 = 0;

	dsp->genomepos = (virtchr[tpp->header.chr].chr << 28) + virtchr[tpp->header.chr].offset + tpp->ppd[0][i].tag1pos;
}
//---------------------------------------------------------------

static void
decodePair(TLBDATA *tpp, unsigned int i, decodedPair_p_t dpp)
{
	if (tpp->header.flags_7_0 & kTLBHflagFixedLength)
	{
		dpp->taglen1 = tpp->lengthInfo;
		dpp->taglen2 = tpp->lengthInfo;
	}
	else
	{
		dpp->taglen1 = tpp->len[2*i];
		dpp->taglen2 = tpp->len[2*i+1];
	}
	memset(dpp->UMI, 0, 8);
	// FIXME - very rigid at this point, fixed to 8 nt
	if (tpp->hdr[9] == ':' && (tpp->hdr[1] == 'A' || tpp->hdr[1] == 'C' || tpp->hdr[1] == 'G' || tpp->hdr[1] == 'T' || tpp->hdr[1] == 'N'))
		memcpy(dpp->UMI, tpp->hdr + 1, 8);
	char *next = strchr(tpp->hdr, '\t');
	tpp->hdr = next + 1;
	next = strchr(tpp->hdr, '\n');
	tpp->hdr = next + 1;
	dpp->ordinal = 0;
	for (unsigned int j = 0; j < 8; j++)
		dpp->ordinal = (dpp->ordinal << 8) | (unsigned long) *((unsigned char *) tpp->hdr++);

	if (tpp->ppd[0][i].tag2offset & k32bREVCOMP_TAG1_BIT)
		dpp->reverseTAG1 = 1;
	else
		dpp->reverseTAG1 = 0;

	if (tpp->ppd[0][i].tag2offset & k32bREVCOMP_TAG2_BIT)
		dpp->reverseTAG2 = 1;
	else
		dpp->reverseTAG2 = 0;

	dpp->delta = (tpp->ppd[0][i].tag2offset & k32bOFFSET_VALUE_MASK);
	if (tpp->ppd[0][i].tag2offset & k32bNEGATIVE_OFFSET_BIT)
		dpp->delta = -dpp->delta;

	dpp->genomepos = (virtchr[tpp->header.chr].chr << 28) + virtchr[tpp->header.chr].offset + tpp->ppd[0][i].tag1pos;
}
//---------------------------------------------------------------

static int
doFilterSingle(unsigned char *diff, unsigned int mmpos, unsigned char *cigar, unsigned int cp)
{
	if (cigar != NULL)
	{
		unsigned int nbID = 0;
		while (cigar[cp] != kCIGARTERMINATOR)
		{
			unsigned int cpos = cigar[cp++];
			char code = cigar[cp++];
			if (code == 'D' || code == 'I')
			{
				if (cpos > 40)
					return 1;
				nbID += 1;
			}
		}
		if (nbID > 2)
			return 1;
	}
#if 0
	if (diff != NULL)
	{
		unsigned int cur = 0;
		unsigned int cnt = 0;
		while (diff[mmpos] != kMMTERMINATOR && diff[mmpos] != kMMSOFTCLIP)
		{
			unsigned int dpos = diff[mmpos++];
			char nt = diff[mmpos++];
			if (dpos != cur + 1)
				cnt += 1;
			cur = dpos;
		}
		if (cnt > 4)
			return 1;
	}
#endif
	return 0;
}
//---------------------------------------------------------------

static int
doFilterPair(unsigned char *diff, unsigned int mmpos, unsigned char *cigar, unsigned int cigarpos)
{
	if (doFilterSingle(diff, mmpos, cigar, cigarpos))
		return 1;
	if (cigar != NULL)
	{
		while (cigar[cigarpos] != kCIGARTERMINATOR)
			cigarpos += 1;
		cigarpos += 1;
	}
#if 0
	if (diff != NULL)
	{
		while (diff[mmpos] != kMMTERMINATOR)
			mmpos += 1;
		mmpos += 1;
	}
#endif
	if (doFilterSingle(diff, mmpos, cigar, cigarpos))
		return 1;
	return 0;
}
//---------------------------------------------------------------

// we keep the min and max pos affected here; when secondRead is true, we avoid counting coverage a second time
static unsigned int minCovPos;
static unsigned int maxCovPos;

static void
cover(coverage_p_t cov, reference_p_t ref,
	const unsigned long ordinal,
	const unsigned int tagPos,
	unsigned int tagLen,
	const unsigned int fragStart,
	const unsigned int fragEnd,
	int reverse,
	unsigned char *diff,
	unsigned int *mmpos,
	unsigned char *cigar,
	unsigned int *cigarpos,
	const unsigned char *qual,
	int ignore,
	int secondRead,
	OPTIONS *opt)
{
	unsigned int j;
	unsigned int pos = tagPos;
	unsigned int cp = *cigarpos;
	unsigned char cc[512]; // the expanded cigar codes
	unsigned char cct[tagLen]; // the expanded cigar codes wrt tag positions
	unsigned int ccp[tagLen]; // genomic positions corresponding to the tag positions
	char ic[tagLen]; // the inserted nucleotides
	memset(cc,'M',tagLen*2);
	memset(cct,'M',tagLen);
	memset(ic,'X',tagLen);
	// need to init the genomic positions, except when there is an explicit cigar... oh well...
	for (j = 0; j < tagLen; j++)
		ccp[j] = tagPos + j;
	if (secondRead == 0)
	{
		minCovPos = 0xFFFFFFFF;
		maxCovPos = 0;
	}
	
	if (cigar != NULL)
	{
		unsigned int a = 0;
		unsigned int r = 0;
		unsigned int gpos = tagPos;

		while (cigar[cp] != kCIGARTERMINATOR)
		{
			unsigned int cpos = cigar[cp++];
			char code = cigar[cp++];
			if (code == 'D')
				tagLen += cpos; // hmm...  I think this is correct
			// while(cpos-- >0)
			// FIXME - something is wrong with some cigar codes, e.g. chr3 in VP15 of MAD data
			while(cpos-- >0 && a < tagLen - 1)
			{
				cc[a++] = code;
				if (code != 'D')
				{
					ccp[r] = gpos;
					cct[r++] = code;
				}
				if (code != 'I')
					gpos += 1;
			}
		}
	}
	if (!ignore)
	{
		if (ref != NULL)
		{
			unsigned int posH = pos >> 8;
			if (ref[posH].cnt == ref[posH].max)
			{
				ref[posH].max += 128;
				ref[posH].rr = realloc(ref[posH].rr,ref[posH].max * sizeof(refRead_t));
			}
			ref[posH].rr[ref[posH].cnt].offset = pos & 0xff;
			ref[posH].rr[ref[posH].cnt].len = tagLen;
			ref[posH].rr[ref[posH].cnt].reverse = reverse;
			ref[posH].rr[ref[posH].cnt].second = secondRead;
			ref[posH].rr[ref[posH].cnt].ordinal = ordinal;
			ref[posH].cnt += 1;
		}
		unsigned int rpos = 0;
#ifdef ORDINAL_TO_PRINT
		if (ordinal == ORDINAL_TO_PRINT) {
			fprintf(stderr,"     "); 
			for (int p=0; p<tagLen; p++) fputc('1'+ (p % 10), stderr);
			fprintf(stderr,"\nSTA: %.*s\nCIG: ", tagLen, cc);
			printCigarStr(&cigar[*cigarpos], stderr);
		}
#endif
		if (reverse)
		{
			for (j = 0; j < tagLen; j++)
			{
				switch (cc[j])
				{
					case 'S': rpos += 1; // rest is same as D: no coverage, fall through
					case 'D': pos += 1; break;
					case 'M':
						if (secondRead == 0)
						{
							if (pos < minCovPos)
								minCovPos = pos;
							if (pos > maxCovPos)
								maxCovPos = pos;
						}
						if (secondRead == 0 || pos < minCovPos || pos > maxCovPos)
						{
							if (qual[tagLen - rpos - 1] - '#' >= opt->filterQualValue && pos >= fragStart && pos <= fragEnd)
								cov[pos].cntM += 1;
#ifdef POSITION_TO_PRINT
							if ( pos == POSITION_TO_PRINT) {
								fprintf(stderr, "Found state M, tag = ?, ordinal=%lu, reverse=%u, qual=%c, P=%u, M=%u",
								        ordinal, reverse, (int) qual[tagLen-rpos-1], cov[pos].cntP,  cov[pos].cntM);
								if (cov[pos].var) {
									for (int l=0; l<4; l++) fprintf(stderr, "\t%u", cov[pos].var->ACGTP[l]);
									for (int l=0; l<4; l++) fprintf(stderr, "\t%u", cov[pos].var->ACGTM[l]);
								}
								fputc('\n', stderr);
							}
#endif
						}
						pos += 1;
						rpos += 1;
						break;
					case 'I': rpos += 1; break;
				}
			}
		}
		else
		{
			for (j = 0; j < tagLen; j++)
			{
				switch (cc[j])
				{
					case 'S': rpos += 1; // rest is same as D: no coverage, fall through
					case 'D': pos += 1; break;
					case 'M':
						if (secondRead == 0)
						{
							if (pos < minCovPos)
								minCovPos = pos;
							if (pos > maxCovPos)
								maxCovPos = pos;
						}
						if (secondRead == 0 || pos < minCovPos || pos > maxCovPos)
						{
							if (qual[rpos] - '#' >= opt->filterQualValue && pos >= fragStart && pos <= fragEnd)
								cov[pos].cntP += 1;
#ifdef POSITION_TO_PRINT
							if ( pos == POSITION_TO_PRINT) {
																fprintf(stderr, "Found state M, tag = ?, ordinal=%lu, reverse=%u, qual=%c, P=%u, M=%u",
								        ordinal, reverse, (int) qual[tagLen-rpos-1], cov[pos].cntP,  cov[pos].cntM);
								if (cov[pos].var) {
									for (int l=0; l<4; l++) fprintf(stderr, "\t%u", cov[pos].var->ACGTP[l]);
									for (int l=0; l<4; l++) fprintf(stderr, "\t%u", cov[pos].var->ACGTM[l]);
								}
								fputc('\n', stderr);
							}
#endif
						}
						pos += 1;
						rpos += 1;
						break;
					case 'I': rpos += 1; break;
				}
			}
		}
	}
	if (diff != NULL)
	{
		while (diff[*mmpos] != kMMTERMINATOR)
		{
			if (diff[*mmpos] == kMMSOFTCLIP)
			{
				// for our purposes, we just need to discard the clipped part
				*mmpos += 1;
				while (diff[*mmpos] != kMMTERMINATOR)
					*mmpos += 1;
				break; // kMMSOFTCLIP is the last thing in MMstring. only one occurence possible. we are done.
			}
			unsigned int dpos = diff[(*mmpos)++];
			char nt = diff[(*mmpos)++];
			if (cct[dpos] == 'S')
				continue;
			if (cct[dpos] == 'I')
			{
				ic[dpos] = nt;
				continue;
			}
			if (cct[dpos] != 'M')
			{
				fprintf(stderr,"HWHAP: cigarcode should be M at %u, but is '%c' at ordinal %lu\n",dpos,cct[dpos],ordinal);
				exit(1);
			}
			if (gNT2bits[(unsigned char) nt] == kNOT_ACGT)
				continue;
			if (!ignore)
			{
				unsigned int gpos = ccp[dpos];
				if (secondRead == 0)
				{
					if (gpos < minCovPos)
						minCovPos = gpos;
					if (gpos > maxCovPos)
						maxCovPos = gpos;
				}
				if (secondRead == 0 || gpos < minCovPos || gpos > maxCovPos)
				{
					unsigned int q = qual[reverse ? tagLen - dpos - 1 : dpos] - '#';
					if (q >= opt->filterQualValue && gpos >= fragStart && gpos <= fragEnd)
					{
						if (cov[gpos].var == NULL)
						{
							cov[gpos].var = malloc(sizeof(gtl_variant_t));
							memset(cov[gpos].var,0,sizeof(gtl_variant_t));
						}
						if (ref != NULL)
						{
							bt_node_p_t n = malloc(sizeof(bt_node_t));
							n->refNo = ordinal;
							n->var = gNT2bits[(unsigned char) nt];
							bt_node_p_t *key = tsearch(n,&cov[gpos].var->ref,bt_node_compare);
							assert(key != NULL);
							if (*key != n)
								free(n);
						}
						if (reverse)
						{
							if (cov[gpos].cntM > 0)
								cov[gpos].cntM -= 1;
							cov[gpos].var->ACGTM[gNT2bits[(unsigned char) nt]] += 1;
						}
						else
						{
							if (cov[gpos].cntP > 0)
								cov[gpos].cntP -= 1;
							cov[gpos].var->ACGTP[gNT2bits[(unsigned char) nt]] += 1;
						}
					}
#ifdef POSITION_TO_PRINT
					if ( gpos == POSITION_TO_PRINT) {
						fprintf(stderr, "Found state m, tag = %c, ordinal=%lu, reverse=%u, qual=%c", (int) nt, ordinal, reverse, (int) qual[reverse ? tagLen - dpos - 1 : dpos]);
						if (cov[gpos].var) {
							for (int l=0; l<4; l++) fprintf(stderr, "\t%u", cov[gpos].var->ACGTP[l]);fprintf(stderr, "\t\t%u\t", cov[gpos].cntP);
							for (int l=0; l<4; l++) fprintf(stderr, "\t%u", cov[gpos].var->ACGTM[l]);fprintf(stderr, "\t\t%u\t", cov[gpos].cntM);
						}
						fputc('\n', stderr);
					}
#endif
				}
			}
		}
		*mmpos += 1;
	}
	if (cigar != NULL)
	{
		unsigned int gpos = tagPos;
		unsigned int rpos = 0;

		while (cigar[*cigarpos] != kCIGARTERMINATOR)
		{
			unsigned int cpos = cigar[(*cigarpos)++];
			char code = cigar[(*cigarpos)++];
			if (!ignore)
			{
				switch(code)
				{
					case 'S':
					case 'M': gpos += cpos; rpos += cpos; break;
					case 'D':
					{
						if (gpos >= fragStart && gpos <= fragEnd)
						{
							if (secondRead == 0)
							{
								if (gpos < minCovPos)
									minCovPos = gpos;
								if (gpos > maxCovPos)
									maxCovPos = gpos;
							}
							if (secondRead == 0 || gpos < minCovPos || gpos > maxCovPos)
							{
								if (cov[gpos].var == NULL)
								{
									cov[gpos].var = malloc(sizeof(gtl_variant_t));
									memset(cov[gpos].var,0,sizeof(gtl_variant_t));
								}
								if (ref != NULL)
								{
									bt_node_p_t n = malloc(sizeof(bt_node_t));
									n->refNo = ordinal;
									n->var = 4; // del
									bt_node_p_t *key = tsearch(n,&cov[gpos].var->ref,bt_node_compare);
									assert(key != NULL);
									if (*key != n)
										free(n);
								}
								if (reverse)
								{
									// was delM here
									for (j = 0; j < kMAXNBINDELS; j++)
										if (cov[gpos].var->del[j].cnt == 0 || cov[gpos].var->del[j].len == cpos)
											break;
									if (j == kMAXNBINDELS)
									{
										fprintf(stderr,"HWHAP: too many dels of different lengths at %u\n",cpos);
										for (j = 0; j < kMAXNBINDELS; j++)
											fprintf(stderr,"DM%u %u %u\n",j,cov[gpos].var->del[j].len,cov[gpos].var->del[j].cnt);
										exit(1);
									}
									cov[gpos].var->del[j].len = cpos;
									cov[gpos].var->del[j].cnt += 1;
								}
								else
								{
									// was delP here
									for (j = 0; j < kMAXNBINDELS; j++)
										if (cov[gpos].var->del[j].cnt == 0 || cov[gpos].var->del[j].len == cpos)
											break;
									if (j == kMAXNBINDELS)
									{
										fprintf(stderr,"HWHAP: too many dels of different lengths at %u\n",cpos);
										for (j = 0; j < kMAXNBINDELS; j++)
											fprintf(stderr,"DP%u %u %u\n",j,cov[gpos].var->del[j].len,cov[gpos].var->del[j].cnt);
										exit(1);
									}
									cov[gpos].var->del[j].len = cpos;
									cov[gpos].var->del[j].cnt += 1;
								}
							}
						}
						gpos += cpos;
						break;
					}
					case 'I':
					{
						if (gpos >= fragStart && gpos <= fragEnd)
						{
							// FIXME - this is wrong... seems we have fixed 38 in the structure definition...
							if (cpos > kMAXNBINDELS)
							{
								fprintf(stderr,"HWHAP: insert is too long: %u\n",cpos);
								exit(1);
							}
							if (secondRead == 0)
							{
								if (gpos < minCovPos)
									minCovPos = gpos;
								if (gpos > maxCovPos)
									maxCovPos = gpos;
							}
							if (secondRead == 0 || gpos < minCovPos || gpos > maxCovPos)
							{
								if (cov[gpos].var == NULL)
								{
									cov[gpos].var = malloc(sizeof(gtl_variant_t));
									memset(cov[gpos].var,0,sizeof(gtl_variant_t));
								}
								if (ref != NULL)
								{
									bt_node_p_t n = malloc(sizeof(bt_node_t));
									n->refNo = ordinal;
									n->var = 5; // ins
									bt_node_p_t *key = tsearch(n,&cov[gpos].var->ref,bt_node_compare);
									assert(key != NULL);
									if (*key != n)
										free(n);
								}
								if (reverse)
								{
									// was insM here
									for (j = 0; j < kMAXNBINDELS; j++)
										if (cov[gpos].var->ins[j].cnt == 0 || (cov[gpos].var->ins[j].nt[cpos] == 0 && strncmp(cov[gpos].var->ins[j].nt,ic+rpos,cpos) == 0))
											break;
									if (j == kMAXNBINDELS)
									{
										fprintf(stderr,"HWHAP: too many inss of different contents\n");
										exit(1);
									}
									memcpy(cov[gpos].var->ins[j].nt,ic+rpos,cpos);
									cov[gpos].var->ins[j].cnt += 1;
								}
								else
								{
									// was insP here
									for (j = 0; j < kMAXNBINDELS; j++)
										if (cov[gpos].var->ins[j].cnt == 0 || (cov[gpos].var->ins[j].nt[cpos] == 0 && strncmp(cov[gpos].var->ins[j].nt,ic+rpos,cpos) == 0))
											break;
									if (j == kMAXNBINDELS)
									{
										fprintf(stderr,"HWHAP: too many inss of different contents\n");
										exit(1);
									}
									memcpy(cov[gpos].var->ins[j].nt,ic+rpos,cpos);
									cov[gpos].var->ins[j].cnt += 1;
								}
							}
						}
						rpos += cpos;
						break;
					}
				}
			}
		}
		*cigarpos += 1;
	}
}
//---------------------------------------------------------------

static void
reportRef1(FILE *f, gtl_variant_p_t var, char gNT, reference_p_t ref, unsigned int offset)
{
	for (unsigned int i = 0; i < ref->cnt; i++)
	{
		if (ref->rr[i].offset <= offset && ref->rr[i].offset + ref->rr[i].len > offset)
		{
			char c;
			if (ref->rr[i].second == 0)
			{
				if (ref->rr[i].reverse == 0)
					c = 'L';
				else
					c = 'l';
			}
			else
			{
				if (ref->rr[i].reverse == 0)
					c = 'r';
				else
					c = 'R';
			}
			bt_node_t n;
			n.refNo = ref->rr[i].ordinal;
			bt_node_p_t *key = tfind(&n, &var->ref, bt_node_compare);
			fprintf(f,"%lu\t%c\t%c\n",(unsigned long) ref->rr[i].ordinal,c,(key == NULL) ? gNT : "ACGTDI"[(*key)->var]);
		}
	}
}
//---------------------------------------------------------------

void
reportRef(FILE *f, gtl_variant_p_t var, char gNT, reference_p_t ref, unsigned int pos, OPTIONS *opt)
{
	if (f == NULL)
		return;
	fprintf(f,"%d\t%u\n",opt->chr,pos);
	unsigned int posH = pos >> 8;
	if (posH > 0)
		reportRef1(f, var, gNT, ref + posH - 1, pos - ((posH - 1) << 8));
	reportRef1(f, var, gNT, ref + posH, pos & 0xff);
	fprintf(f,"END\n");
}
//---------------------------------------------------------------

void
loadSelected(coverage_p_t cov, OPTIONS *opt)
{
	if (opt->selectRegions[0] == 0)
		return;
	char linbuf[kMaxLineBuf];
	FILE *f = fopen (opt->selectRegions, "r");
	if (!f)
	{
		fprintf(stderr, "Error: Cannot open %s\n",opt->selectRegions);
		exit(1);
	}
	while(fgets(linbuf,kMaxLineBuf,f) != NULL)
	{
		char tmpchr[32];
		unsigned int p1, p2;
		sscanf(linbuf,"%s\t%u\t%u",tmpchr,&p1,&p2);
		//if ((tmpchr[0] == 'c') && (tmpchr[1] == 'h') && (tmpchr[2] == 'r'))
		//{
		//}
		if (strcmp(virtchr[opt->chr].AC,tmpchr))
			continue;
		if (p1 > opt->selectAdd)
			p1 -= opt->selectAdd;
		else
			p1 = 0;
		p2 += opt->selectAdd;
		for (unsigned int i = p1; i <= p2; i++)
			cov[i].flags = 1;
	}
	fclose(f);
}
//---------------------------------------------------------------

void
computeCov(coverage_p_t cov, OPTIONS *opt, coverage_p_t refCov, unsigned int *covCntNb, unsigned int *median, unsigned int *iqr25, unsigned int *iqr75)
{
	// compute median and iqr
	unsigned int covCnt[65536];
	*covCntNb = 0;
	memset(covCnt,0,sizeof(covCnt));
	for (unsigned int i = 0; i <= virtchr[opt->chr].len; i++)
	{
		if (refCov != NULL && refCov[i].flags == 0)
			continue;
		unsigned int c = cov[i].cov[0] + cov[i].cov[1] + cov[i].cov[2];
		if (c >= 65536)
			c = 65535;
		if (c > 0)
		{
			covCnt[c] += 1;
			*covCntNb += 1;
		}
	}
	*median = 0xffffffff, *iqr25 = 0xffffffff, *iqr75 = 0xffffffff;
	{
		unsigned int c = 0;
		for (unsigned int i = 0; i < 65536; i++)
		{
			c += covCnt[i];
			if (*iqr25 == 0xffffffff && c >= *covCntNb / 4)
				*iqr25 = i;
			if (*median == 0xffffffff && c >= *covCntNb / 2)
				*median = i;
			if (*iqr75 == 0xffffffff && c >= 3 * *covCntNb / 4)
			{
				*iqr75 = i;
				break;
			}
		}
	}
}
//---------------------------------------------------------------

int
decompress(const char * const restrict fn, OPTIONS *opt, coverage_p_t cov, reference_p_t ref)
{
	int fd = open(fn, O_RDONLY);
	TLBHEADER th;
	int err = 1;
	// discard clones
	unsigned int cpos1 = 0, cpos2 = 0;
	TLBDATA td;
	allocBlock(&td, true);
	unsigned char cUMI[8];
	memset(cUMI,0,8);
	while (1) {
		ssize_t res = read(fd,&th,sizeof(TLBHEADER));
		if (res == 0)
		{
			err = 0;
			break;
		}
		if (res != sizeof(TLBHEADER))
			goto bail;
		unsigned char hcs = simple8bitCS((unsigned char *) &th, sizeof(TLBHEADER));
		if (hcs != 0)
		{
			fprintf(stderr,"Checksum mismatch: %02x %02x\n",th.headerCS,hcs);
			goto bail;
		}
		// do we want this block ?
		if (th.chr == opt->chr)
		{
			if ((th.flags_7_0 & kTLBHflagPairPosBlock) != 0)
			{
				// we count the coverage of Ns in the reads as if they were reference (i.e., exact matches)
				if ((th.flags_7_0 & 0xf0) == 0
						|| (th.flags_7_0 & 0xf0) == kTLBHflagNMismatchBlock
						|| (th.flags_7_0 & 0xf0) == kTLBHflagMismatchBlock
						|| (th.flags_7_0 & 0xf0) == (kTLBHflagMismatchBlock | kTLBHflagCigarBlock))
				{
					// we do...
					err = readDecompressBlock(fd, &th, &td, 0);
					if (err != 0)
						goto bail;
					unsigned int covSelect = 0;
					unsigned char *diff = td.diff[0];
					if ((th.flags_7_0 & 0xf0) == kTLBHflagNMismatchBlock)
						diff = NULL; // no need to take those into account - and they break the current decoder anyway
					if ((th.flags_7_0 & 0xf0) == kTLBHflagMismatchBlock)
						covSelect = 1;
					if ((th.flags_7_0 & 0xf0) == (kTLBHflagMismatchBlock | kTLBHflagCigarBlock))
						covSelect = 2;
					if ((td.header.flags_7_0 & kTLBHflagPairedReads) == 0)
					{
						decodedSingle_t ds;
						unsigned int mmpos = 0;
						unsigned int cigarpos = 0;
						int i;
						// just add the coverage
						for (i=0; i<td.cnt; i++)
						{
							decodeSingle(&td,i,&ds);
							// we have some rare cases where a read starts mapping "before" the beginning of the chromosome (see badPos fastq examples)
							int ignore = 0;
							if (td.ppd[0][i].tag1pos == 0 || td.ppd[0][i].tag1pos > virtchr[opt->chr].len)
								ignore = 1;
							// FIXME - maybe should also check same tag lengths...
							if (!opt->noCloneFilter && td.ppd[0][i].tag1pos == cpos1)
								ignore = 1;
							if (opt->additionalFilter && doFilterSingle(diff,mmpos,td.cigar[0],cigarpos))
								ignore = 1;
							cpos1 = td.ppd[0][i].tag1pos;
							if (!ignore)
							{
								unsigned int start = cpos1;
								unsigned int end = cpos1 + ds.taglen1;
								for (unsigned int k = start; k < end; k++)
									cov[k].cov[covSelect] += 1;
							}
							cover(cov,ref,ds.ordinal,td.ppd[0][i].tag1pos, ds.taglen1, td.ppd[0][i].tag1pos, td.ppd[0][i].tag1pos + ds.taglen1 - 1, ds.reverseTAG1, diff, &mmpos, td.cigar[0], &cigarpos, td.qs, ignore, 0, opt);
							td.qs += ds.taglen1;
						}
					}
					else
					{
						decodedPair_t dp;
						unsigned int mmpos = 0;
						unsigned int cigarpos = 0;
						int i;
						// just add the coverage
						for (i=0; i<td.cnt; i++)
						{
							// NOTE - we have some rare cases where a read starts mapping "before" the beginning of the chromosome (see badPos fastq examples)
							decodePair(&td,i,&dp);
							// only keep well-behaved pairs unless maxPairLen is 0
							int ignore = 0;
							if (opt->maxPairLen != 0 && (abs(dp.delta) >= opt->maxPairLen || dp.reverseTAG1 == dp.reverseTAG2))
								ignore = 1;
							unsigned int tag2pos = td.ppd[0][i].tag1pos + dp.delta;
							if (td.ppd[0][i].tag1pos == 0 || td.ppd[0][i].tag1pos > virtchr[opt->chr].len
									|| tag2pos == 0 || tag2pos > virtchr[opt->chr].len)
								ignore = 1;
							// FIXME - maybe should also check same tag lengths...
							if (!opt->noCloneFilter)
							{
								// check UMI if we think we have one
								unsigned int matchedUMI = 0;
								for (unsigned int k = 0; k < 8; k++)
									if (cUMI[k] == dp.UMI[k])
										matchedUMI += 1;
								if (matchedUMI >= 7 && td.ppd[0][i].tag1pos == cpos1 && tag2pos == cpos2)
									ignore = 1;
								if (matchedUMI >= 7 && td.ppd[0][i].tag1pos == cpos2 && tag2pos == cpos1)
									ignore = 1;
							}
							if (opt->additionalFilter && doFilterPair(diff,mmpos,td.cigar[0],cigarpos))
								ignore = 1;
							cpos1 = td.ppd[0][i].tag1pos;
							cpos2 = tag2pos;
							memcpy(cUMI, dp.UMI, 8);
							unsigned int fragStart = td.ppd[0][i].tag1pos;
							unsigned int fragEnd = tag2pos + dp.taglen2 - 1;
							if (dp.reverseTAG1)
							{
								fragStart = tag2pos;
								fragEnd = td.ppd[0][i].tag1pos + dp.taglen1 - 1;
							}
							if (!ignore)
							{
								for (unsigned int k = fragStart; k <= fragEnd; k++)
									cov[k].cov[covSelect] += 1;
							}
							cover(cov,ref,dp.ordinal,td.ppd[0][i].tag1pos, dp.taglen1, fragStart, fragEnd, dp.reverseTAG1, diff, &mmpos, td.cigar[0], &cigarpos, td.qs, ignore, 0, opt);
							td.qs += dp.taglen1;
#ifdef BYPASS_COMMON
							minCovPos = 0xFFFFFFFFU;
							maxCovPos = 0U;
#endif
							cover(cov,ref,dp.ordinal,tag2pos, dp.taglen2, fragStart, fragEnd, dp.reverseTAG2, diff, &mmpos, td.cigar[0], &cigarpos, td.qs, ignore, 1, opt);
							td.qs += dp.taglen2;
						}
					}
				}
			}
		}
		else
		{
			// nope, read and discard (so that we can stream...)
			size_t size = (size_t) th.blockLength - sizeof(TLBHEADER);
			void *data = malloc(size);
			ssize_t res = read(fd,data,size);
			if (res != size)
			{
				perror("skip read:");
				goto bail;
			}
			free(data);
		}
	}

	err = 0;

bail:;

	// clean up
	close(fd);
	freeBlock(&td);

	return(err);

} // decompress
//---------------------------------------------------------------

/* Extract the Tags and State Sequences, setting score above zero if realignment is useful */  
void 
ExtractStateSequence(decodedPair_p_t dpp, unsigned char * qual,
                     const unsigned char * restrict diff, const unsigned char * restrict cigar,
                     GlobalData_t * const restrict GD)
{
	unsigned char tmpStr[kMaxReadLen];
	
	const unsigned char * restrict copydiff = diff;
	const unsigned char * restrict copycigar = cigar;
	/* NOTE: A loop will do no harm here :) */
	
	{
		unsigned char * const StateSequence1 = GD[0].StateSequence;
		unsigned char * const tag1 = GD[0].Tag;
		const int taglen1 = dpp->taglen1;
		int StateSequenceLength1 = 0;
		int WorthRealign = 0;
		assert(taglen1<kMaxReadLen);
		memcpy(tag1,&genome_tbl[dpp->genomepos],taglen1);
		GD[0].AlignmentRange[0] = 0;
		if (cigar != NULL && *cigar != kCIGARTERMINATOR) {
			int r=0, a=0;
			const unsigned char * restrict genome = &genome_tbl[dpp->genomepos];
			while (*cigar != kCIGARTERMINATOR)
			{
				register int pos = cigar[0];
				const char code  = cigar[1];
				cigar += 2;
				
				{
					assert(StateSequenceLength1+pos < STATE_MEMORY_SIZE);
					unsigned char * const restrict ptr = &StateSequence1[StateSequenceLength1];
					for (int i=0; i<pos; i++) ptr[i] = code;
					StateSequenceLength1 += pos;
				}
				
				switch(code)
				{
					case 'S': // same as 'M'
					case 'M': assert(a+pos<=kMaxReadLen); while(pos-- >0) { tag1[a++] = genome[r]; r++; } break;
					case 'D':                          while(pos-- >0) { r++;                        } break;
					case 'I': assert(a+pos<=kMaxReadLen); while(pos-- >0) { tag1[a++] = 'N';            } break;
				}
			}
			tag1[a] = '\0';
			assert(r >= 1);
			GD[0].AlignmentRange[1] = r - 1;
			WorthRealign = 1;
		}
		else {
			memcpy(tag1,&genome_tbl[dpp->genomepos],taglen1);
			memset(StateSequence1, 'M', taglen1);
			StateSequenceLength1 = taglen1;
			tag1[taglen1] = '\0';
			GD[0].AlignmentRange[1] = taglen1 - 1;
		}
		if (cigar != NULL) 
			if (*cigar == kCIGARTERMINATOR) cigar++;
		assert(StateSequenceLength1 < STATE_MEMORY_SIZE);
		StateSequence1[StateSequenceLength1] = '\0';
			
		if (diff) {
			int count = 0;
			const unsigned char * restrict genome = &genome_tbl[dpp->genomepos];
			unsigned char * restrict StateSequencePtr = StateSequence1;
			while (*diff != kMMTERMINATOR)
			{
				if (*diff == kMMSOFTCLIP)
				{
					// grab the nucleotides clipped in a tmp string
					int i = 0;
					diff++;
					while (*diff != kMMTERMINATOR) { tmpStr[i++] = *diff++ ; }
					// copy the clipped nucleotides at the end of the tag.
					char * const copyto = (char*) ((dpp->reverseTAG1) ? &tag1[0] : &tag1[taglen1-i]);
					memcpy(copyto,tmpStr,i);
					break; // kMMSOFTCLIP is the last thing in MMstring. only one occurence possible. we are done.
				}
				const int pos = diff[0];
				assert(pos<taglen1);
				const unsigned char nt = diff[1];
				tag1[pos] = nt;
				diff += 2;
				{
					while (count < pos) {
						switch(*StateSequencePtr) {
							case 'M': 
							case 'S':
							case 'm':
							case 'n':
							case 'N': genome++;
							case 'I': assert(count <= taglen1); count++; break;
							case 'D': genome++; break;
	#ifndef NDEBUG
							default: 
								fprintf(stderr, "Unknown state %c encountered!\nSTA: %.*s\n",
												(int) StateSequence1[count], StateSequenceLength1, StateSequence1);
								exit(1);
	#endif
						}
						StateSequencePtr++;
					}
					while (*StateSequencePtr == 'D') ++StateSequencePtr;
	// 				printf("count %i, pos %i\n", count , pos);
					if (*StateSequencePtr != 'I') {
						if (nt == 'N') 
							*StateSequencePtr = 'N';
						else if (*genome == 'N') 
							*StateSequencePtr = 'n';
						else
							*StateSequencePtr = 'm';
					}
				}
			}
	// 		printf("\n");
			diff++; // skip the kMMTERMINATOR
			WorthRealign = 1;
		}
#ifdef ORDINAL_TO_PRINT
		if (dpp->ordinal == ORDINAL_TO_PRINT) {
			if (copydiff) {
				fprintf(stderr, "STA: %s\n", StateSequence1);
				fputs("TAG: ", stderr);
				const unsigned char * restrict cptr = tag1; 
				for (int k=0; k<StateSequenceLength1; k++) fputc((StateSequence1[k] == 'D') ? '-' : *cptr++, stderr);
				fputs("\nGEN: ", stderr);
				cptr = &genome_tbl[dpp->genomepos];
#ifdef POSITION_TO_PRINT
				unsigned int location = GD[0].MismatchCount;
				int where = -1;
#endif
				for (int k=0; k<StateSequenceLength1; k++) {
					fputc((StateSequence1[k] == 'I') ? '-' : *cptr++, stderr);
#ifdef POSITION_TO_PRINT
					if (StateSequence1[k] != 'I')
						if (location++ == POSITION_TO_PRINT) where = k;
#endif
				}
#ifdef POSITION_TO_PRINT
				if (where >= 0) {
					fputs("\n     ", stderr);
					for (int k=0; k<where; k++) fputc(' ', stderr);
					fprintf(stderr, "^ %u\n     *=%u == %u\n", POSITION_TO_PRINT, dpp->genomepos, GD[0].MismatchCount);
				}
				else
#endif
				fprintf(stderr, "\n     *=%u == %u\n", dpp->genomepos, GD[0].MismatchCount);
				fprintf(stderr, "CIG: "); 
				if (copycigar) 
					printCigarStr(copycigar, stderr);
				else 
					fputc('\n', stderr);
				fprintf(stderr, "MM : ");
				if (copydiff)
					printMismatchStr(copydiff, stderr);
				else
					fputc('\n', stderr);
			}
		}
#endif
		GD[0].score = WorthRealign;
		GD[0].States = GD[0].StateSequence; 
	}

	GD[0].CigarLength    = (uintptr_t) cigar - (uintptr_t) copycigar; 
	GD[0].MismatchLength = (uintptr_t) diff - (uintptr_t) copydiff;
	GD[0].revNmer        = dpp->reverseTAG1;
	
	copydiff = diff;
	copycigar = cigar;
	 
	{
		unsigned char * const StateSequence2 = GD[1].StateSequence;
		unsigned char * const tag2 = GD[1].Tag;
		const int taglen2 = dpp->taglen2;
		int StateSequenceLength2 = 0;
		int WorthRealign = 0;
		assert(taglen2<kMaxReadLen);
		GD[1].AlignmentRange[0] = 0;
		memcpy(tag2,&genome_tbl[dpp->genomepos + dpp->delta],taglen2);
		if (cigar != NULL && *cigar != kCIGARTERMINATOR) {
			int r=0, a=0;
			const unsigned char * restrict genome = &genome_tbl[dpp->genomepos + dpp->delta];
			while (*cigar != kCIGARTERMINATOR)
			{
				register int pos = cigar[0];
				const char code = cigar[1];
				cigar += 2;
				
				{
					assert(StateSequenceLength2+pos < STATE_MEMORY_SIZE);
					unsigned char * const restrict ptr = &StateSequence2[StateSequenceLength2];
					for (int i=0; i<pos; i++) ptr[i] = code;
					StateSequenceLength2 += pos;
				}
				
				switch(code)
				{
					case 'S': // same as 'M'
					case 'M': assert(a+pos<=kMaxReadLen); while(pos-- >0) { tag2[a++] = genome[r]; r++; } break;
					case 'D':                          while(pos-- >0) { r++;                        } break;
					case 'I': assert(a+pos<=kMaxReadLen); while(pos-- >0) { tag2[a++] = 'N';            } break;
				}
			}
			tag2[a] = '\0';
			assert(r >= 1);
			GD[1].AlignmentRange[1] = r - 1;
			WorthRealign = 1;
		}
		else {
			memcpy(tag2,&genome_tbl[dpp->genomepos + dpp->delta],taglen2);
			memset(StateSequence2, 'M', taglen2);
			StateSequenceLength2 = taglen2;
			tag2[taglen2] = '\0';
			GD[1].AlignmentRange[1] = taglen2 - 1;
		}
		if (cigar != NULL) 
			if (*cigar == kCIGARTERMINATOR) cigar++;
		assert(StateSequenceLength2 < STATE_MEMORY_SIZE);
		StateSequence2[StateSequenceLength2] = '\0';
			
		if (diff) {
			int count = 0;
			const unsigned char * restrict genome = &genome_tbl[dpp->genomepos + dpp->delta];
			unsigned char * restrict StateSequencePtr = StateSequence2;
			while (*diff != kMMTERMINATOR)
			{
				if (*diff == kMMSOFTCLIP)
				{
					// grab the nucleotides clipped in a tmp string
					int i = 0;
					diff++;
					while (*diff != kMMTERMINATOR) { tmpStr[i++] = *diff++ ; }
					// copy the clipped nucleotides at the end of the tag.
					char * const copyto = (char*) ((dpp->reverseTAG2) ? &tag2[0] : &tag2[taglen2-i]);
					memcpy(copyto,tmpStr,i);
					break; // kMMSOFTCLIP is the last thing in MMstring. only one occurence possible. we are done.
				}
				const int pos = diff[0];
				assert(pos<taglen2);
				const unsigned char nt = diff[1];
				tag2[pos] = nt;
				diff += 2;
				{
					while (count < pos) {
						switch(*StateSequencePtr) {
							case 'M': 
							case 'S':
							case 'm':
							case 'n':
							case 'N': genome++;
							case 'I': assert(count <= taglen2); count++; break;
							case 'D': genome++; break;
	#ifndef NDEBUG
							default: 
								fprintf(stderr, "Unknown state %c encountered!\nSTA: %.*s\n",
												(int) StateSequence2[count], StateSequenceLength2, StateSequence2);
								exit(1);
	#endif
						}
						StateSequencePtr++;
					}
					while (*StateSequencePtr == 'D') ++StateSequencePtr;
	// 				printf("count %i, pos %i\n", count , pos);
					if (*StateSequencePtr != 'I') {
						if (nt == 'N') 
							*StateSequencePtr = 'N';
						else if (*genome == 'N') 
							*StateSequencePtr = 'n';
						else
							*StateSequencePtr = 'm';
					}
				}
			}
	// 		printf("\n");
			diff++; // skip the kMMTERMINATOR
			WorthRealign = 1;
		}
#ifdef ORDINAL_TO_PRINT
		if (dpp->ordinal == ORDINAL_TO_PRINT) {
			if (copydiff) {
				fprintf(stderr, "STA: %s\n", StateSequence2);
				fputs("TAG: ", stderr);
				const unsigned char * restrict cptr = tag2; 
				for (int k=0; k<StateSequenceLength2; k++) fputc((StateSequence2[k] == 'D') ? '-' : *cptr++, stderr);
				fputs("\nGEN: ", stderr);
				cptr = &genome_tbl[dpp->genomepos + dpp->delta];
#ifdef POSITION_TO_PRINT
				unsigned int location = GD[1].MismatchCount;
				int where = -1;
#endif
				for (int k=0; k<StateSequenceLength2; k++) {
					fputc((StateSequence2[k] == 'I') ? '-' : *cptr++, stderr);
#ifdef POSITION_TO_PRINT
					if (StateSequence2[k] != 'I')
						if (location++ == POSITION_TO_PRINT) where = k;
#endif
				}
#ifdef POSITION_TO_PRINT
				if (where >= 0) {
					fputs("\n     ", stderr);
					for (int k=0; k<where; k++) fputc(' ', stderr);
					fprintf(stderr, "^ %u\n     *=%u == %u\n", POSITION_TO_PRINT, dpp->genomepos + dpp->delta, GD[1].MismatchCount);
				}
				else
#endif
				fprintf(stderr, "\n     *=%u == %u\n", dpp->genomepos + dpp->delta, GD[1].MismatchCount);
				fprintf(stderr,"CIG: "); 
				if (copycigar) 
					printCigarStr(copycigar, stderr);
				else 
					fputc('\n', stderr);
				fprintf(stderr, "MM : ");
				if (copydiff)
					printMismatchStr(copydiff, stderr);
				else
					fputc('\n', stderr);
			}
		}
#endif
		GD[1].score = WorthRealign;
		GD[1].States = GD[1].StateSequence; 
	}

	GD[1].CigarLength    = (uintptr_t) cigar - (uintptr_t) copycigar; 
	GD[1].MismatchLength = (uintptr_t) diff - (uintptr_t) copydiff; 
	GD[1].revNmer        = dpp->reverseTAG2;
}
//---------------------------------------------------------------

typedef struct flags {
	unsigned char reverse:1;
	unsigned char secondRead:1;
} flags_t;

typedef struct covers {
	unsigned int min;
	unsigned int max;
} Cover_t;

static void
MyCover(coverage_t * const restrict cov, reference_t * const restrict ref, const unsigned long ordinal,
        char * restrict StateSequence, const char * restrict Tag, const unsigned char * restrict Qual, 
        const unsigned int fragStart, const unsigned int fragEnd, Cover_t * const restrict Cover,
        unsigned int TagPos, const flags_t flags, const OPTIONS * const restrict opt
#ifdef POSITION_TO_PRINT
        , const unsigned int ID
#endif
)
{
	// we keep the min and max pos affected here; when secondRead is true, we avoid counting coverage a second time
	const unsigned char QualityThreshold = opt->filterQualValue;
	
	/* WELL WELL, quality is not in the orientation of the genome !!! */
	int QualStep;
	if (!flags.reverse) {
		QualStep = 1;
	}
	else {
		QualStep = -1;
	}
	
#ifdef KEEP_TRACK_OF_ORDINAL
	/* Allocate a key */
	bt_node_p_t n = malloc(sizeof(bt_node_t));
	assert(n);
	n->refNo = ordinal;
#endif
	
	/**************************************** FIRST TAG *************************************************/
	if (flags.secondRead == 0) {
		Cover_t lCover = { .min=0xFFFFFFFFU, .max=0U };
		
		/* StateSequence is null terminated, loop over it */
		while (1) {
			const char CurrentState = *StateSequence;
#ifdef POSITION_TO_PRINT
			unsigned int tpos = 0;
			if ( TagPos == POSITION_TO_PRINT) {
				fprintf(stderr, "Thread %u, Tag 1 Found state %c, tag = %c, ordinal=%lu, reverse=%u, qual=%c, ",
				        ID, (int) CurrentState, (int) *Tag, ordinal, (unsigned int) flags.reverse, (int) *Qual);
				tpos = TagPos;
				if (TagPos >= fragStart && TagPos <= fragEnd) {
					fputs("frag IN, ", stderr);
				}
				else {
					fputs("frag OFF, ", stderr);
				}
				const char quality = (*Qual) - '#';
				if (quality >= QualityThreshold) {
					fputs("quality OK, ", stderr);
				}
				else {
					fputs("quality FAIL, ", stderr);
				}
			}
#endif
			switch (CurrentState) {
				case 'S': {
					/* Check if first match in read */
					if (!flags.reverse) {
						if (TagPos < lCover.min) lCover.min = TagPos;
						if (TagPos > lCover.max) lCover.max = TagPos;
					}
					
					/* move pointers */
					TagPos++;
					Qual += QualStep;
					Tag++;
					StateSequence++;
					break;
				}
				case 'N': {
					/* Check if first match in read */
					if (TagPos < lCover.min) lCover.min = TagPos;
					if (TagPos > lCover.max) lCover.max = TagPos;
					
					/*************************************** TO BE REMOVED ******************************************/
					const char quality = (*Qual) - '#';
					if (quality >= QualityThreshold)
					/* Account the match in the right orientation */
					{
						unsigned short * const Counter = (flags.reverse) ? &cov[TagPos].cntM : &cov[TagPos].cntP;
						*Counter += 1;
					}
					/*************************************** TO BE REMOVED ******************************************/
					
					/* move pointers */
					TagPos++;
					Qual += QualStep;
					Tag++;
					StateSequence++;
					break;
				}
				case 'n':
				case 'm': {
					/* Check if first match in read */
					if (TagPos < lCover.min) lCover.min = TagPos;
					if (TagPos > lCover.max) lCover.max = TagPos;
					
					/* Are we within fragment range */
					if (TagPos >= fragStart && TagPos <= fragEnd) {
						const char quality = (*Qual) - '#';
						/* Are we within quality range */
						if (quality >= QualityThreshold) {
							/* Check that the variant memory is already allocated, if not do it */
							if (!cov[TagPos].var) {
								cov[TagPos].var = calloc(1,sizeof(gtl_variant_t));
								assert(cov[TagPos].var);
							}
							
							/* Get the key to the reference */
							const unsigned char ID = gNT2bits[(unsigned char) *Tag];
#ifdef KEEP_TRACK_OF_ORDINAL
							{
								n->var = ID;
								bt_node_p_t *key = tsearch(n, &cov[TagPos].var->ref, bt_node_compare);
								assert(key != NULL);
								if (*key == n) {
									/* Realloc some memory for the next one */
									n = malloc(sizeof(bt_node_t));
									assert(n);
									n->refNo = ordinal;
								}
							}
#endif
							/* Account for the mismatch base */
							{
								unsigned short * const Counter = (flags.reverse) ? cov[TagPos].var->ACGTM : cov[TagPos].var->ACGTP;
								Counter[ID] += 1;
							}
						}
					}
					
					/* move pointers */
					TagPos++;
					Qual += QualStep;
					Tag++;
					StateSequence++;
					break;
				}
				case 'M': {
					/* Check if first match in read */
					if (TagPos < lCover.min) lCover.min = TagPos;
					if (TagPos > lCover.max) lCover.max = TagPos;
					
					/* Are we within fragment range */
					if (TagPos >= fragStart && TagPos <= fragEnd) {
						const char quality = (*Qual) - '#';
						/* Are we within quality range */
						if (quality >= QualityThreshold) {
							/* Account the match in the right orientation */
							{
								unsigned short * const Counter = (flags.reverse) ? &cov[TagPos].cntM : &cov[TagPos].cntP;
								*Counter += 1;
							}
						}
					}
					
					/* move pointers */
					TagPos++;
					Qual += QualStep;
					Tag++;
					StateSequence++;
					break;
				}
				case 'I': {
					/* Get the length of the insertion, moving StatSequence at the same time */
					const uintptr_t location = (uintptr_t) StateSequence;
					while (*StateSequence == 'I') ++StateSequence;
					const size_t len = (uintptr_t) StateSequence - (uintptr_t) location;
					assert(len > 0);
					
					/* Perform quality check, hence all bases must be above threashold, moving Qual at the same time */ 
					unsigned char ShallWeSkipIt = 0;
					for (size_t k=0; k<len; k++) {
						if (*Qual < QualityThreshold) ShallWeSkipIt = 1;
						Qual += QualStep;
					}
					
					if (!ShallWeSkipIt) {
						/* Are we within range */
						if (TagPos >= fragStart && TagPos <= fragEnd) {
							/* Check that the variant memory is already allocated, if not do it */
							if (!cov[TagPos].var) {
								cov[TagPos].var = calloc(1,sizeof(gtl_variant_t));
								assert(cov[TagPos].var);
							}
							
#ifdef KEEP_TRACK_OF_ORDINAL
							/* Get the key to the reference */
							{
								n->var = 5; // Id for insertion
								bt_node_p_t *key = tsearch(n, &cov[TagPos].var->ref, bt_node_compare);
								assert(key != NULL);
								if (*key == n) {
									/* Realloc some memory for the next one */
									n = malloc(sizeof(bt_node_t));
									assert(n);
									n->refNo = ordinal;
								}
							}
#endif
							
							/* Account the mismatch in the right orientation */
							{
								struct InternalIns * restrict ins = cov[TagPos].var->ins;
								int j;
								for (j=0; j<kMAXNBINDELS; j++, ins++) {
									/* Check that it is not first and not already existing */
									if (ins->cnt == 0 || (ins->nt[len] == 0 && strncmp(ins->nt, Tag, len) == 0) )
										break;
								}
								if (j == kMAXNBINDELS) goto InDelOverflow;
									
								ins->cnt += 1;
								assert(len <= 38);
								for (size_t k=0; k<len; k++) ins->nt[k] = Tag[k];
#ifdef POSITION_TO_PRINT
								if (tpos > 0) fprintf(stderr, " insert=%.*s length=%lu", (int) len, Tag, len);
#endif
							}
						}
					}
					
					/* Moving pointers */
					Tag += len;
					TagPos++;
					break;
				}
				case 'D': {
					/* Get the length of the deletion, moving StatSequence at the same time */
					const uintptr_t location = (uintptr_t) StateSequence;
					while (*StateSequence == 'D') ++StateSequence;
					const size_t len = (uintptr_t) StateSequence - (uintptr_t) location;
					assert(len > 0);
					
					/* Are we within range */
					if (TagPos >= fragStart && TagPos <= fragEnd) {
						/* Check that the variant memory is already allocated, if not do it */
						if (!cov[TagPos].var) {
							cov[TagPos].var = calloc(1,sizeof(gtl_variant_t));
							assert(cov[TagPos].var);
						}
						
#ifdef KEEP_TRACK_OF_ORDINAL
						/* Get the key to the reference */
						{
							n->var = 4; // Id for deletion
							bt_node_p_t *key = tsearch(n, &cov[TagPos].var->ref, bt_node_compare);
							assert(key != NULL);
							if (*key == n) {
								/* Realloc some memory for the next one */
								n = malloc(sizeof(bt_node_t));
								assert(n);
								n->refNo = ordinal;
							}
						}
#endif
						/* Account the mismatch in the right orientation */
						{
							struct InternalDel * restrict del = cov[TagPos].var->del;
							int j;
							for (j=0; j<kMAXNBINDELS; j++, del++) {
								/* Check that it is not first and not already existing */
								if (del->cnt == 0 || del->len == len )
									break;
							}
							if (j == kMAXNBINDELS) goto InDelOverflow;
								
							del->cnt += 1;
							del->len = len;
						}
					}
					/* Moving pointers */
					TagPos += len;
					break;
				}
				case '\0': *Cover = lCover; goto Done;
#ifndef DNDEBUG
				default: 
					fprintf(stderr, "Unknown state encountered '%c'\n", (int) CurrentState);
					exit(1);
#endif
			}
#ifdef POSITION_TO_PRINT
			if ( tpos > 0 ) {
				fprintf(stderr, "P=%u, M=%u", cov[tpos].cntP, cov[tpos].cntM);
				if (cov[tpos].var) {
					for (int l=0; l<4; l++) fprintf(stderr, "\t%u", cov[tpos].var->ACGTP[l]);
					for (int l=0; l<4; l++) fprintf(stderr, "\t%u", cov[tpos].var->ACGTM[l]);
				}
				fputc('\n', stderr);
			}
#endif
		}
	}
	/**************************************** SECOND TAG *************************************************/
	else {
		/* StateSequence is null terminated, loop over it */
		const Cover_t lCover = *Cover;
		while (1) {
			const register char CurrentState = *StateSequence;
#ifdef POSITION_TO_PRINT
			unsigned int tpos = 0;
			if ( TagPos == POSITION_TO_PRINT) {
				fprintf(stderr, "Thread %u, Tag 2 Found state %c, tag = %c, ordinal=%lu, reverse=%u, qual=%c, ",
				        ID, (int) CurrentState, (int) *Tag, ordinal, (unsigned int) flags.reverse, (int) *Qual);
				tpos = TagPos;
				if (TagPos >= fragStart && TagPos <= fragEnd) {
					fputs("frag IN, ", stderr);
				}
				else {
					fputs("frag OFF, ", stderr);
				}
				const char quality = (*Qual) - '#';
				if (quality >= QualityThreshold) {
					fputs("quality OK, ", stderr);
				}
				else {
					fputs("quality FAIL, ", stderr);
				}
				if (TagPos > lCover.max || TagPos < lCover.min) {
					fputs("overlap NO, ", stderr);
				}
				else {
					fprintf(stderr, "overlap YES, range=[%u,%u], ", lCover.min, lCover.max);
				}
			}
#endif
			switch (CurrentState) {
				case 'S': {
					/* move pointers */
					TagPos++;
					Qual += QualStep;
					Tag++;
					StateSequence++;
					break;
				}
				case 'N': {
					/*************************************** TO BE REMOVED ******************************************/
					if (TagPos >= fragStart && TagPos <= fragEnd && (TagPos > lCover.max || TagPos < lCover.min)) {
						const char quality = (*Qual) - '#';
						if (quality >= QualityThreshold)
						/* Account the match in the right orientation */
						{
							unsigned short * const Counter = (flags.reverse) ? &cov[TagPos].cntM : &cov[TagPos].cntP;
							*Counter += 1;
						}
					}
					/*************************************** TO BE REMOVED ******************************************/
					
					/* move pointers */
					TagPos++;
					Qual += QualStep;
					Tag++;
					StateSequence++;
					break;
				}
				case 'n':
				case 'm': {
					/* Are we within fragment range */
					if (TagPos >= fragStart && TagPos <= fragEnd && (TagPos > lCover.max || TagPos < lCover.min)) {
						const char quality = (*Qual) - '#';
						/* Are we within quality range */
						if (quality >= QualityThreshold) {
							/* Check that the variant memory is already allocated, if not do it */
							if (!cov[TagPos].var) {
								cov[TagPos].var = calloc(1,sizeof(gtl_variant_t));
								assert(cov[TagPos].var);
							}
							
							/* Get the key to the reference */
							const unsigned char ID = gNT2bits[(unsigned char) *Tag];
							
#ifdef KEEP_TRACK_OF_ORDINAL
							{
								n->var = ID;
								bt_node_p_t *key = tsearch(n, &cov[TagPos].var->ref, bt_node_compare);
								assert(key != NULL);
								if (*key == n) {
									/* Realloc some memory for the next one */
									n = malloc(sizeof(bt_node_t));
									assert(n);
									n->refNo = ordinal;
								}
							}
#endif
							/* Account for the mismatch base */
							{
								unsigned short * const Counter = (flags.reverse) ? cov[TagPos].var->ACGTM : cov[TagPos].var->ACGTP;
								Counter[ID] += 1;
							}
						}
					}
					
					/* move pointers */
					TagPos++;
					Qual += QualStep;
					Tag++;
					StateSequence++;
					break;
				}
				case 'M': {
					/* Are we within fragment range */
					if (TagPos >= fragStart && TagPos <= fragEnd && (TagPos > lCover.max || TagPos < lCover.min)) {
						const char quality = (*Qual) - '#';
						/* Are we within quality range */
						if (quality >= QualityThreshold) {
							/* Account the match in the right orientation */
							{
								unsigned short * const Counter = (flags.reverse) ? &cov[TagPos].cntM : &cov[TagPos].cntP;
								*Counter += 1;
							}
						}
					}
					
					/* move pointers */
					TagPos++;
					Qual += QualStep;
					Tag++;
					StateSequence++;
					break;
				}
				case 'I': {
					/* Get the length of the insertion, moving StatSequence at the same time */
					const uintptr_t location = (uintptr_t) StateSequence;
					while (*StateSequence == 'I') ++StateSequence;
					const size_t len = (uintptr_t) StateSequence - (uintptr_t) location;
					
					/* Perform quality check, hence all bases must be above threashold, moving Qual at the same time */ 
					unsigned char ShallWeSkipIt = 0;
					for (size_t k=0; k<len; k++) {
						if (*Qual < QualityThreshold) ShallWeSkipIt = 1;
						Qual += QualStep;
					}
					
					if (!ShallWeSkipIt) {
						/* Are we within range */
						if (TagPos >= fragStart && TagPos <= fragEnd && (TagPos > lCover.max || TagPos < lCover.min)) {
							/* Check that the variant memory is already allocated, if not do it */
							if (!cov[TagPos].var) {
								cov[TagPos].var = calloc(1,sizeof(gtl_variant_t));
								assert(cov[TagPos].var);
							}
							
#ifdef KEEP_TRACK_OF_ORDINAL
							/* Get the key to the reference */
							{
								n->var = 5; // Id for insertion
								bt_node_p_t *key = tsearch(n, &cov[TagPos].var->ref, bt_node_compare);
								assert(key != NULL);
								if (*key == n) {
									/* Realloc some memory for the next one */
									n = malloc(sizeof(bt_node_t));
									assert(n);
									n->refNo = ordinal;
								}
							}
#endif
							/* Account the mismatch in the right orientation */
							{
								struct InternalIns * restrict ins = cov[TagPos].var->ins;
								int j;
								for (j=0; j<kMAXNBINDELS; j++, ins++) {
									/* Check that it is not first and not already existing */
									if (ins->cnt == 0 || (ins->nt[len] == 0 && strncmp(ins->nt, Tag, len) == 0) )
										break;
								}
								if (j == kMAXNBINDELS) goto InDelOverflow;
									
								ins->cnt += 1;
								assert(len <= 38);
								for (size_t k=0; k<len; k++) ins->nt[k] = Tag[k];
#ifdef POSITION_TO_PRINT
								if (tpos > 0) fprintf(stderr, " insert=%.*s length=%lu", (int) len, Tag, len);
#endif
							}
						}
					}
					
					/* Moving pointers */
					Tag += len;
					TagPos++;
					break;
				}
				case 'D': {
					/* Get the length of the deletion, moving StatSequence at the same time */
					const uintptr_t location = (uintptr_t) StateSequence;
					while (*StateSequence == 'D') ++StateSequence;
					const size_t len = (uintptr_t) StateSequence - (uintptr_t) location;
					
					/* Are we within range */
					if (TagPos >= fragStart && TagPos <= fragEnd && (TagPos > lCover.max || TagPos < lCover.min)) {
						/* Check that the variant memory is already allocated, if not do it */
						if (!cov[TagPos].var) {
							cov[TagPos].var = calloc(1,sizeof(gtl_variant_t));
							assert(cov[TagPos].var);
						}
						
#ifdef KEEP_TRACK_OF_ORDINAL
						/* Get the key to the reference */
						{
							n->var = 4; // Id for deletion
							bt_node_p_t *key = tsearch(n, &cov[TagPos].var->ref, bt_node_compare);
							assert(key != NULL);
							if (*key == n) {
								/* Realloc some memory for the next one */
								n = malloc(sizeof(bt_node_t));
								assert(n);
								n->refNo = ordinal;
							}
						}
#endif
						/* Account the mismatch in the right orientation */
						{
							struct InternalDel * restrict del = cov[TagPos].var->del;
							int j;
							for (j=0; j<kMAXNBINDELS; j++, del++) {
								/* Check that it is not first and not already existing */
								if (del->cnt == 0 || del->len == len )
									break;
							}
							if (j == kMAXNBINDELS) goto InDelOverflow;
								
							del->cnt += 1;
							del->len = len;
						}
					}
					/* Moving pointers */
					TagPos += len;
					break;
				}
				case '\0': goto Done;
#ifndef DNDEBUG
				default: 
					fprintf(stderr, "Unknown state encountered '%c'\n", (int) CurrentState);
					exit(1);
#endif
			}
#ifdef POSITION_TO_PRINT
			if ( tpos > 0 ) {
				fprintf(stderr, "P=%u, M=%u", cov[tpos].cntP, cov[tpos].cntM);
				if (cov[tpos].var) {
					for (int l=0; l<4; l++) fprintf(stderr, "\t%u", cov[tpos].var->ACGTP[l]);
					for (int l=0; l<4; l++) fprintf(stderr, "\t%u", cov[tpos].var->ACGTM[l]);
				}
				fputc('\n', stderr);
			}
#endif
		}
	}
	
	Done: ;
#ifdef KEEP_TRACK_OF_ORDINAL
	/* Free Key memory */
	free(n);
#endif
	
	return;
	
	InDelOverflow:
		if (StateSequence[-1] == 'I') {
			fprintf(stderr, "HWHAP: too many inss of different contents at genome location %u\n", TagPos);
			struct InternalIns * restrict ins = cov[TagPos].var->ins;
			for (int j=0; j<kMAXNBINDELS; j++, ins++) {
				/* Check that it is not first and not already existing */
				fprintf(stderr, "%i\t%u\t%s\n", j, ins[j].cnt, ins[j].nt);
				if (ins->cnt == 0) break;
			}
		}
		else {
			fprintf(stderr, "HWHAP: too many dels of different sizes at genome location %u\n", TagPos);
			for (int j = 0; j < kMAXNBINDELS; j++)
				fprintf(stderr,"DM%u %u %u\n",j,cov[TagPos].var->del[j].len,cov[TagPos].var->del[j].cnt);
		}
		
		exit(1);
}
//---------------------------------------------------------------

int
decompress_realign(const char * const restrict fn, OPTIONS *opt, coverage_p_t cov, reference_p_t ref)
{
	int fd = open(fn, O_RDONLY);
	TLBHEADER th;
	int err = 1;
	// discard clones
	unsigned int cpos1 = 0, cpos2 = 0;
	TLBDATA td;
	allocBlock(&td, true);
	
	unsigned char tag1[kMaxReadLen];
	unsigned char tag2[kMaxReadLen];
	
	/* Global data for realignment */
	GlobalData_t GD[2];
	GD[0].Tag = tag1;
	GD[0].TagLength = kMaxReadLen;
	GD[0].GenomeLength = kGAPPED_ALIGN_GENOME_LENGTH + prf.ExtraLeftGenomeShift;
	GD[1].Tag = tag2;
	GD[1].TagLength = kMaxReadLen;
	GD[1].GenomeLength = kGAPPED_ALIGN_GENOME_LENGTH + prf.ExtraLeftGenomeShift;
	
	cpu_std_border.allocStorage(&GD[0]);
	cpu_std_border.allocStorage(&GD[1]);
	
	while (1) {
		ssize_t res = read(fd,&th,sizeof(TLBHEADER));
		if (res == 0)
		{
			err = 0;
			break;
		}
		if (res != sizeof(TLBHEADER))
			goto bail;
		unsigned char hcs = simple8bitCS((unsigned char *) &th, sizeof(TLBHEADER));
		if (hcs != 0)
		{
			fprintf(stderr,"Checksum mismatch: %02x %02x\n",th.headerCS,hcs);
			goto bail;
		}
		// do we want this block ?
		if (th.chr == opt->chr)
		{
			if ((th.flags_7_0 & kTLBHflagPairPosBlock) != 0)
			{
				// we count the coverage of Ns in the reads as if they were reference (i.e., exact matches)
				if ((th.flags_7_0 & 0xf0) == 0
						|| (th.flags_7_0 & 0xf0) == kTLBHflagNMismatchBlock
						|| (th.flags_7_0 & 0xf0) == kTLBHflagMismatchBlock
						|| (th.flags_7_0 & 0xf0) == (kTLBHflagMismatchBlock | kTLBHflagCigarBlock))
				{
					// we do...
					err = readDecompressBlock(fd, &th, &td, 0);
					if (err != 0)
						goto bail;
					unsigned int covSelect = 0;
					unsigned char *diff = td.diff[0];
					if ((th.flags_7_0 & 0xf0) == kTLBHflagNMismatchBlock)
						diff = NULL; // no need to take those into account - and they break the current decoder anyway
					if ((th.flags_7_0 & 0xf0) == kTLBHflagMismatchBlock)
						covSelect = 1;
					if ((th.flags_7_0 & 0xf0) == (kTLBHflagMismatchBlock | kTLBHflagCigarBlock))
						covSelect = 2;
					if ((td.header.flags_7_0 & kTLBHflagPairedReads) == 0)
					{
						decodedSingle_t ds;
						unsigned int mmpos = 0;
						unsigned int cigarpos = 0;
						int i;
						// just add the coverage
						for (i=0; i<td.cnt; i++)
						{
							decodeSingle(&td,i,&ds);
							// we have some rare cases where a read starts mapping "before" the beginning of the chromosome (see badPos fastq examples)
							int ignore = 0;
							if (td.ppd[0][i].tag1pos == 0 || td.ppd[0][i].tag1pos > virtchr[opt->chr].len)
								ignore = 1;
							// FIXME - maybe should also check same tag lengths...
							if (!opt->noCloneFilter && td.ppd[0][i].tag1pos == cpos1)
								ignore = 1;
							cpos1 = td.ppd[0][i].tag1pos;
							if (!ignore)
							{
								unsigned int start = cpos1;
								unsigned int end = cpos1 + ds.taglen1;
								for (unsigned int k = start; k < end; k++)
									cov[k].cov[covSelect] += 1;
							}
							cover(cov,ref,ds.ordinal,td.ppd[0][i].tag1pos, ds.taglen1, td.ppd[0][i].tag1pos, td.ppd[0][i].tag1pos + ds.taglen1 - 1, ds.reverseTAG1, diff, &mmpos, td.cigar[0], &cigarpos, td.qs, ignore, 0, opt);
							td.qs += ds.taglen1;
						}
					}
					else
					{
						decodedPair_t dp;
						unsigned int mmpos = 0;
						unsigned int cigarpos = 0;
						int i;
						// just add the coverage
						for (i=0; i<td.cnt; i++)
						{
							// NOTE - we have some rare cases where a read starts mapping "before" the beginning of the chromosome (see badPos fastq examples)
							//decodePair(&td,i,&dp);
							// decodePair(TLBDATA *tpp, unsigned int i, decodedPair_p_t dpp)
							{
								if (td.header.flags_7_0 & kTLBHflagFixedLength)
								{
									dp.taglen1 = td.lengthInfo;
									dp.taglen2 = td.lengthInfo;
								}
								else
								{
									dp.taglen1 = td.len[2*i];
									dp.taglen2 = td.len[2*i+1];
								}
								char *next = strchr(td.hdr, '\t');
								td.hdr = next + 1;
								next = strchr(td.hdr, '\n');
								td.hdr = next + 1;
								dp.ordinal = 0;
								for (unsigned int j = 0; j < 8; j++)
									dp.ordinal = (dp.ordinal << 8) | (unsigned long) *((unsigned char *) td.hdr++);

								if (td.ppd[0][i].tag2offset & k32bREVCOMP_TAG1_BIT)
									dp.reverseTAG1 = 1;
								else
									dp.reverseTAG1 = 0;

								if (td.ppd[0][i].tag2offset & k32bREVCOMP_TAG2_BIT)
									dp.reverseTAG2 = 1;
								else
									dp.reverseTAG2 = 0;

								dp.delta = (td.ppd[0][i].tag2offset & k32bOFFSET_VALUE_MASK);
								if (td.ppd[0][i].tag2offset & k32bNEGATIVE_OFFSET_BIT)
									dp.delta = -dp.delta;

								dp.genomepos = (virtchr[td.header.chr].chr << 28) + virtchr[td.header.chr].offset + td.ppd[0][i].tag1pos;
							}
							
							// What do we have
#ifdef POSITION_TO_PRINT
							GD[0].MismatchCount = td.ppd[0][i].tag1pos;
							GD[1].MismatchCount = td.ppd[0][i].tag1pos + dp.delta;
#endif
							ExtractStateSequence(&dp, td.qs, (diff != NULL) ? &diff[mmpos] : NULL, td.cigar[0] + cigarpos, GD );

							/* Hack to move pointers */
							cigarpos += GD[0].CigarLength + GD[1].CigarLength;
							mmpos    += GD[0].MismatchLength + GD[1].MismatchLength;
							
							// only keep well-behaved pairs for now
							int ignore = (abs(dp.delta) >= kMAXPAIRLEN || dp.reverseTAG1 == dp.reverseTAG2);
							unsigned int tag2pos = td.ppd[0][i].tag1pos + dp.delta;
							ignore |= (td.ppd[0][i].tag1pos == 0 || td.ppd[0][i].tag1pos > virtchr[opt->chr].len || tag2pos == 0 || tag2pos > virtchr[opt->chr].len);

							// FIXME - maybe should also check same tag lengths...
							ignore |= (!opt->noCloneFilter && td.ppd[0][i].tag1pos == cpos1 && tag2pos == cpos2);
							ignore |= (!opt->noCloneFilter && td.ppd[0][i].tag1pos == cpos2 && tag2pos == cpos1);
							
							cpos1 = td.ppd[0][i].tag1pos;
							cpos2 = tag2pos;
							
							if (!ignore) {
								const unsigned char * restrict quality;
								flags_t flags;
								Cover_t Cover;
								unsigned int fragStart, fragEnd;
								const struct GenomePosition { unsigned int start; unsigned int stop; } *LeftTag, *RightTag; 
								
								if (dp.delta >= 0) {
									LeftTag = (const struct GenomePosition*) GD[0].AlignmentRange;
									RightTag = (const struct GenomePosition*) GD[1].AlignmentRange;
								}
								else {
									LeftTag = (const struct GenomePosition*) GD[1].AlignmentRange;
									RightTag = (const struct GenomePosition*) GD[0].AlignmentRange;
								}
								
								/* Compute potential tag overlap for too short DNA tags 
								 * PROBLEMATIC CASE IS <<<<<<XXXXXXXXXX>>>>>>>
								 *              TAG L  <--------------<
								 *              TAG R        >--------------->
								 */
								
								if (LeftTag->stop >= RightTag->start) {
									fragStart = RightTag->start;
									fragEnd   = LeftTag->stop; 
								}
								else {
									fragStart = LeftTag->start;
									fragEnd   = RightTag->stop;
								}
								
								if (dp.reverseTAG1) {
									flags.reverse = 1;
									quality = &td.qs[dp.taglen1 - 1];
									GD[0].revNmer = 1;
								} 
								else {
									flags.reverse = 0;
									quality = td.qs;
									GD[0].revNmer = 0;
								}
								flags.secondRead = 0;
								
								for (unsigned int k = fragStart; k <= fragEnd; k++)
									cov[k].cov[covSelect] += 1;
								
								/* Realign if useful as determined in score attribute in ExtractStateSequence */
								if (GD[0].score > 0) {
#ifdef ORDINAL_TO_PRINT
									if (dp.ordinal == ORDINAL_TO_PRINT) {
										GlobalData_t * ReAlignSettings = GD;
										GD[0].TagLength = dp.taglen1;
										GD[0].Genome = genome_tbl + dp.genomepos;
										if (cpos1 >= prf.ExtraLeftGenomeShift) GD[0].Genome -= prf.ExtraLeftGenomeShift;
										
										cpu_std_border.createMatrix(&GD[0]);
										if (cpu_std_border.getStateSequence(&GD[0]) != 0) {
											fprintf(stderr, "ERROR OOPS CANNOT REALIGN\n");
											exit(1);
										}
										int iTag = 0;
										fprintf(stderr, "TAG %i (len=%zu) @ %i-%i/%i-%i, score=%i, mismatch=%u, sofclip mismatch=%u, RequiredMM=%u, softclip evicted MM=%u,"
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
										fprintf(stderr, "STA: ");
										for ( int space=0; space<ReAlignSettings[iTag].ProfileRange[0]; space++) fputc(' ', stderr);
										fprintf(stderr, "%s\n", ReAlignSettings[iTag].States);
										fprintf(stderr, "TAG: %.*s", (int) ReAlignSettings[iTag].ProfileRange[0], ReAlignSettings[iTag].Tag);
										int l=0, k=0;
										while (ReAlignSettings[iTag].States[l] != 0) {
											if (ReAlignSettings[iTag].States[l] == 'D') 
												fputc(' ', stderr);
											else 
												fputc((int) ReAlignSettings[iTag].Tag[ReAlignSettings[iTag].ProfileRange[0]+(k++)], stderr);
											l++;
										}
										fprintf(stderr, "\nGEN: "); 
										for ( int space=0; space<ReAlignSettings[iTag].ProfileRange[0]; space++) fputc(' ', stderr);
										l=0; k=0;
										while (ReAlignSettings[iTag].States[l] != 0) {
											if (ReAlignSettings[iTag].States[l] == 'I')
												fputc(' ', stderr);
											else
												fputc(ReAlignSettings[iTag].Genome[ReAlignSettings[iTag].AlignmentRange[0]+(k++)], stderr);
											l++;
										}
										fputc('\n', stderr);
									}
#endif
								}
								
								MyCover(cov, ref, dp.ordinal, (char*) GD[0].States, (char*) GD[0].Tag, quality, fragStart, fragEnd, &Cover, cpos1, flags, opt
#ifdef POSITION_TO_PRINT
								        ,0
#endif
								);
								td.qs += dp.taglen1;
								
#ifdef BYPASS_COMMON
								Cover.min = 0xFFFFFFFFU;
								Cover.max = 0U;
#endif
								if (dp.reverseTAG2) {
									flags.reverse = 1;
									quality = &td.qs[dp.taglen2 - 1];
									GD[1].revNmer = 1;
								} 
								else {
									flags.reverse = 0;
									quality = td.qs;
									GD[1].revNmer = 0;
								}
								flags.secondRead = 1;
								
								/* Realign if useful as determined in score attribute in ExtractStateSequence */
								if (GD[1].score > 0) {
#ifdef ORDINAL_TO_PRINT
									if (dp.ordinal == ORDINAL_TO_PRINT) {
										GlobalData_t * ReAlignSettings = GD;
										GD[1].TagLength = dp.taglen2;
										GD[1].Genome = genome_tbl + dp.genomepos + dp.delta;
										if (cpos2 >= prf.ExtraLeftGenomeShift) GD[1].Genome -= prf.ExtraLeftGenomeShift;
										
										cpu_std_border.createMatrix(&GD[1]);
										if (cpu_std_border.getStateSequence(&GD[1]) != 0) {
											fprintf(stderr, "ERROR OOPS CANNOT REALIGN\n");
											exit(1);
										}
										int iTag = 1;
										fprintf(stderr, "TAG %i (len=%zu) @ %i-%i/%i-%i, score=%i, mismatch=%u, sofclip mismatch=%u, RequiredMM=%u, softclip evicted MM=%u,"
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
										fprintf(stderr, "STA: ");
										for ( int space=0; space<ReAlignSettings[iTag].ProfileRange[0]; space++) fputc(' ', stderr);
										fprintf(stderr, "%s\n", ReAlignSettings[iTag].States);
										fprintf(stderr, "TAG: %.*s", (int) ReAlignSettings[iTag].ProfileRange[0], ReAlignSettings[iTag].Tag);
										int l=0, k=0;
										while (ReAlignSettings[iTag].States[l] != 0) {
											if (ReAlignSettings[iTag].States[l] == 'D') 
												fputc(' ', stderr);
											else 
												fputc((int) ReAlignSettings[iTag].Tag[ReAlignSettings[iTag].ProfileRange[0]+(k++)], stderr);
											l++;
										}
										fprintf(stderr, "\nGEN: "); 
										for ( int space=0; space<ReAlignSettings[iTag].ProfileRange[0]; space++) fputc(' ', stderr);
										l=0; k=0;
										while (ReAlignSettings[iTag].States[l] != 0) {
											if (ReAlignSettings[iTag].States[l] == 'I')
												fputc(' ', stderr);
											else
												fputc(ReAlignSettings[iTag].Genome[ReAlignSettings[iTag].AlignmentRange[0]+(k++)], stderr);
											l++;
										}
										fputc('\n', stderr);
									}
#endif
								}
								
								MyCover(cov, ref, dp.ordinal, (char*) GD[1].States,(char*) GD[1].Tag, quality, fragStart, fragEnd, &Cover, cpos2, flags, opt
#ifdef POSITION_TO_PRINT
								        ,0
#endif
								);
								td.qs += dp.taglen2;
								
							}
							else {
								/* move pointers */
								td.qs += dp.taglen1;
								td.qs += dp.taglen2;
							}
						}
					}
				}
			}
		}
		else
		{
			// nope, read and discard (so that we can stream...)
			size_t size = (size_t) th.blockLength - sizeof(TLBHEADER);
			void *data = malloc(size);
			ssize_t res = read(fd,data,size);
			if (res != size)
			{
				perror("skip read:");
				goto bail;
			}
			free(data);
		}
	}

	err = 0;

bail:;

	cpu_std_border.freeStorage(&GD[0]);
	cpu_std_border.freeStorage(&GD[1]);

	// clean up
	close(fd);
	freeBlock(&td);

	return(err);

} // decompress_realign
//---------------------------------------------------------------

typedef struct ThreadArgs {
	coverage_t * cov;
	reference_t * ref;
	OPTIONS *opt;
	unsigned int ID;
} ThreadArgs_t;

static TLBDATA * * BlockData;
static TLBHEADER * BlockHeader;
static sem_t *semThreadStart;
static sem_t *semThreadDone;

static void* CoverBlock(ThreadArgs_t * const args)
{
	const OPTIONS * const restrict opt = args->opt;
	coverage_t * const restrict cov = args->cov;
	reference_t * const restrict ref = args->ref;
	const unsigned int ID = args->ID;
	
	unsigned char tag1[kMaxReadLen];
	unsigned char tag2[kMaxReadLen];
	
	/* Global data for realignment */
	GlobalData_t GD[2];
	GD[0].Tag = tag1;
	GD[0].TagLength = kMaxReadLen;
	GD[0].GenomeLength = kGAPPED_ALIGN_GENOME_LENGTH + prf.ExtraLeftGenomeShift;
	GD[1].Tag = tag2;
	GD[1].TagLength = kMaxReadLen;
	GD[1].GenomeLength = kGAPPED_ALIGN_GENOME_LENGTH + prf.ExtraLeftGenomeShift;
	
	cpu_std_border.allocStorage(&GD[0]);
	cpu_std_border.allocStorage(&GD[1]);
	
	/* Which semaphore to wait for ? */
	sem_t * const restrict semaphore = &semThreadStart[ID];
	
	while (1) {
// 		fprintf(stderr, "Thread %u waiting...\n", ID);
		int ret;
		while ((ret = sem_wait(semaphore)) != 0) {
			if (ret != EINTR) goto bail;
		}
		TLBDATA* const restrict td = BlockData[ID];
		
// 		fprintf(stderr, "Thread %u starts working...\n", ID);
		if (td != NULL) {
			TLBHEADER* const restrict th = &BlockHeader[ID];
			{
				void *dataBlock = td->diskBuffer;
				const unsigned int size = th->blockLength - sizeof(TLBHEADER);
				if (decompressBlock(dataBlock,size,th,td,0)) {
					fputs("Error decompressing block\n", stderr);
					goto bail;
				}
			}
			
			unsigned int covSelect = 0;
			unsigned char * restrict diff = td->diff[0];
			if ((th->flags_7_0 & 0xf0) == kTLBHflagNMismatchBlock) 
				diff = NULL; // no need to take those into account - and they break the current decoder anyway
			if ((th->flags_7_0 & 0xf0) == kTLBHflagMismatchBlock) covSelect = 1;
			if ((th->flags_7_0 & 0xf0) == (kTLBHflagMismatchBlock | kTLBHflagCigarBlock)) covSelect = 2;
			if ((td->header.flags_7_0 & kTLBHflagPairedReads)) {
				decodedPair_t dp;
				unsigned int mmpos = 0U;
				unsigned int cigarpos = 0U;
				// WARNING: discard clones, potentially not correct with threads !!!
				unsigned int cpos1 = 0U, cpos2 = 0U;
				// just add the coverage
				for (int i=0; i<td->cnt; i++)
				{
					/* Basic Decode pair */
					{
						if (td->header.flags_7_0 & kTLBHflagFixedLength)
						{
							dp.taglen1 = td->lengthInfo;
							dp.taglen2 = td->lengthInfo;
						}
						else
						{
							dp.taglen1 = td->len[2*i];
							dp.taglen2 = td->len[2*i+1];
						}
						char *next = strchr(td->hdr, '\t');
						td->hdr = next + 1;
						next = strchr(td->hdr, '\n');
						td->hdr = next + 1;
						dp.ordinal = 0;
						for (unsigned int j = 0; j < 8; j++)
							dp.ordinal = (dp.ordinal << 8) | (unsigned long) *((unsigned char *) td->hdr++);

						if (td->ppd[0][i].tag2offset & k32bREVCOMP_TAG1_BIT)
							dp.reverseTAG1 = 1;
						else
							dp.reverseTAG1 = 0;

						if (td->ppd[0][i].tag2offset & k32bREVCOMP_TAG2_BIT)
							dp.reverseTAG2 = 1;
						else
							dp.reverseTAG2 = 0;

						dp.delta = (td->ppd[0][i].tag2offset & k32bOFFSET_VALUE_MASK);
						if (td->ppd[0][i].tag2offset & k32bNEGATIVE_OFFSET_BIT)
							dp.delta = -dp.delta;

						dp.genomepos = (virtchr[td->header.chr].chr << 28) + virtchr[td->header.chr].offset + td->ppd[0][i].tag1pos;
					}
					
					// What do we have
#ifdef POSITION_TO_PRINT
					GD[0].MismatchCount = td->ppd[0][i].tag1pos;
					GD[1].MismatchCount = td->ppd[0][i].tag1pos + dp.delta;
#endif
					ExtractStateSequence(&dp, td->qs, (diff != NULL) ? &diff[mmpos] : NULL, td->cigar[0] + cigarpos, GD );

					/* Hack to move pointers */
					cigarpos += GD[0].CigarLength + GD[1].CigarLength;
					mmpos    += GD[0].MismatchLength + GD[1].MismatchLength;
					
					// only keep well-behaved pairs for now
					int ignore = (abs(dp.delta) >= kMAXPAIRLEN || dp.reverseTAG1 == dp.reverseTAG2);
					unsigned int tag2pos = td->ppd[0][i].tag1pos + dp.delta;
					ignore |= (td->ppd[0][i].tag1pos == 0 || td->ppd[0][i].tag1pos > virtchr[opt->chr].len || tag2pos == 0 || tag2pos > virtchr[opt->chr].len);

					// FIXME - maybe should also check same tag lengths...
					ignore |= (!opt->noCloneFilter && td->ppd[0][i].tag1pos == cpos1 && tag2pos == cpos2);
					ignore |= (!opt->noCloneFilter && td->ppd[0][i].tag1pos == cpos2 && tag2pos == cpos1);
					
					cpos1 = td->ppd[0][i].tag1pos;
					cpos2 = tag2pos;
					
					if (!ignore) {
						const unsigned char * restrict quality;
						flags_t flags;
						Cover_t Cover;
						unsigned int fragStart, fragEnd;
						unsigned int LeftTagGenLen, RightTagGenLen;
						unsigned int LeftTagGenPos, RightTagGenPos; 
						
						if (dp.delta >= 0) {
							LeftTagGenLen  = GD[0].AlignmentRange[1] - GD[0].AlignmentRange[0];
							LeftTagGenPos  = cpos1;
							RightTagGenLen = GD[1].AlignmentRange[1] - GD[1].AlignmentRange[0];
							RightTagGenPos = cpos2;
						}
						else {
							LeftTagGenLen  = GD[1].AlignmentRange[1] - GD[1].AlignmentRange[0];
							LeftTagGenPos  = cpos2;
							RightTagGenLen = GD[0].AlignmentRange[1] - GD[0].AlignmentRange[0];
							RightTagGenPos = cpos1;
						}
						
						/* Compute potential tag overlap for too short DNA tags 
							* PROBLEMATIC CASE IS <<<<<<XXXXXXXXXX>>>>>>>
							*              TAG L  <--------------<
							*              TAG R        >--------------->
							*/
						
						if (LeftTagGenPos + LeftTagGenLen >= RightTagGenPos) {
							fragStart = RightTagGenPos;
							fragEnd   = LeftTagGenPos + LeftTagGenLen; 
						}
						else {
							fragStart = LeftTagGenPos;
							fragEnd   = RightTagGenPos + RightTagGenLen;
						}
						
						for (unsigned int k = fragStart; k <= fragEnd; k++) cov[k].cov[covSelect] += 1;
						
						if (dp.reverseTAG1) {
							flags.reverse = 1;
							quality = &td->qs[dp.taglen1 - 1];
							GD[0].revNmer = 1;
						} 
						else {
							flags.reverse = 0;
							quality = td->qs;
							GD[0].revNmer = 0;
						}
						flags.secondRead = 0;
#ifndef NDEBUG
						if (cpos1 < cpos2) {
							if (cpos2-cpos1 >= kMAXPAIRLEN) {
								fprintf(stderr, "Potential miss, distance between tags to short, %u < %u\n", cpos2-cpos1, kMAXPAIRLEN);
								exit(1);
							}
						}
						else {
							if (cpos1-cpos2 >= kMAXPAIRLEN) {
								fprintf(stderr, "Potential miss, distance between tags to short, %u < %u\n", cpos1-cpos2, kMAXPAIRLEN);
								exit(1);
							}
						}
#endif
						/* Realign if useful */
						if (GD[0].score > 0) {
#ifdef ORDINAL_TO_PRINT
							if (dp.ordinal == ORDINAL_TO_PRINT) {
								GlobalData_t * ReAlignSettings = GD;
								GD[0].TagLength = dp.taglen1;
								GD[0].Genome = genome_tbl + dp.genomepos;
								if (cpos1 >= prf.ExtraLeftGenomeShift) GD[0].Genome -= prf.ExtraLeftGenomeShift;
								
								cpu_std_border.createMatrix(&GD[0]);
								if (cpu_std_border.getStateSequence(&GD[0]) != 0) {
									fprintf(stderr, "ERROR OOPS CANNOT REALIGN\n");
									exit(1);
								}
								int iTag = 0;
								fprintf(stderr, "TAG %i (len=%zu) @ %i-%i/%i-%i, score=%i, mismatch=%u, sofclip mismatch=%u, RequiredMM=%u, softclip evicted MM=%u,"
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
								fprintf(stderr, "STA: ");
								for ( int space=0; space<ReAlignSettings[iTag].ProfileRange[0]; space++) fputc(' ', stderr);
								fprintf(stderr, "%s\n", ReAlignSettings[iTag].States);
								fprintf(stderr, "TAG: %.*s", (int) ReAlignSettings[iTag].ProfileRange[0], ReAlignSettings[iTag].Tag);
								int l=0, k=0;
								while (ReAlignSettings[iTag].States[l] != 0) {
									if (ReAlignSettings[iTag].States[l] == 'D') 
										fputc(' ', stderr);
									else 
										fputc((int) ReAlignSettings[iTag].Tag[ReAlignSettings[iTag].ProfileRange[0]+(k++)], stderr);
									l++;
								}
								fprintf(stderr, "\nGEN: "); 
								for ( int space=0; space<ReAlignSettings[iTag].ProfileRange[0]; space++) fputc(' ', stderr);
								l=0; k=0;
								while (ReAlignSettings[iTag].States[l] != 0) {
									if (ReAlignSettings[iTag].States[l] == 'I')
										fputc(' ', stderr);
									else
										fputc(ReAlignSettings[iTag].Genome[ReAlignSettings[iTag].AlignmentRange[0]+(k++)], stderr);
									l++;
								}
								fputc('\n', stderr);
							}
#endif
						}
						
						MyCover(cov, ref, dp.ordinal, (char*) GD[0].States, (char*) GD[0].Tag, quality, fragStart, fragEnd, &Cover, cpos1, flags, opt
#ifdef POSITION_TO_PRINT
						        ,ID
#endif
						);
						td->qs += dp.taglen1;
						
#ifdef BYPASS_COMMON
						Cover.min = 0xFFFFFFFFU;
						Cover.max = 0U;
#endif
						if (dp.reverseTAG2) {
							flags.reverse = 1;
							quality = &td->qs[dp.taglen2 - 1];
							GD[1].revNmer = 1;
						} 
						else {
							flags.reverse = 0;
							quality = td->qs;
							GD[1].revNmer = 0;
						}
						flags.secondRead = 1;
						/* Realign if useful */
						if (GD[1].score > 0) {
#ifdef ORDINAL_TO_PRINT
							if (dp.ordinal == ORDINAL_TO_PRINT) {
								GlobalData_t * ReAlignSettings = GD;
								GD[1].TagLength = dp.taglen2;
								GD[1].Genome = genome_tbl + dp.genomepos + dp.delta;
								if (cpos2 >= prf.ExtraLeftGenomeShift) GD[1].Genome -= prf.ExtraLeftGenomeShift;
								
								cpu_std_border.createMatrix(&GD[1]);
								if (cpu_std_border.getStateSequence(&GD[1]) != 0) {
									fprintf(stderr, "ERROR OOPS CANNOT REALIGN\n");
									exit(1);
								}
								int iTag = 1;
								fprintf(stderr, "TAG %i (len=%zu) @ %i-%i/%i-%i, score=%i, mismatch=%u, sofclip mismatch=%u, RequiredMM=%u, softclip evicted MM=%u,"
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
								fprintf(stderr, "STA: ");
								for ( int space=0; space<ReAlignSettings[iTag].ProfileRange[0]; space++) fputc(' ', stderr);
								fprintf(stderr, "%s\n", ReAlignSettings[iTag].States);
								fprintf(stderr, "TAG: %.*s", (int) ReAlignSettings[iTag].ProfileRange[0], ReAlignSettings[iTag].Tag);
								int l=0, k=0;
								while (ReAlignSettings[iTag].States[l] != 0) {
									if (ReAlignSettings[iTag].States[l] == 'D') 
										fputc(' ', stderr);
									else 
										fputc((int) ReAlignSettings[iTag].Tag[ReAlignSettings[iTag].ProfileRange[0]+(k++)], stderr);
									l++;
								}
								fprintf(stderr, "\nGEN: "); 
								for ( int space=0; space<ReAlignSettings[iTag].ProfileRange[0]; space++) fputc(' ', stderr);
								l=0; k=0;
								while (ReAlignSettings[iTag].States[l] != 0) {
									if (ReAlignSettings[iTag].States[l] == 'I')
										fputc(' ', stderr);
									else
										fputc(ReAlignSettings[iTag].Genome[ReAlignSettings[iTag].AlignmentRange[0]+(k++)], stderr);
									l++;
								}
								fputc('\n', stderr);
							}
#endif
						}

						MyCover(cov, ref, dp.ordinal, (char*) GD[1].States,(char*) GD[1].Tag, quality, fragStart, fragEnd, &Cover, cpos2, flags, opt
#ifdef POSITION_TO_PRINT
						        ,ID
#endif
						);
						td->qs += dp.taglen2;
						
					}
					else {
						/* move pointers */
						td->qs += dp.taglen1;
						td->qs += dp.taglen2;
					}
				}
			}
		}
		else {
			break;
		}

// 		fprintf(stderr, "Thread %u done...\n", ID);
		sem_post(&semThreadDone[ID]);
		
	}
// 	fprintf(stderr, "Thread %u, bye bye...\n", ID);
	
	return NULL;
	bail: ;
		exit(1);
}
//---------------------------------------------------------------

static void* RealignCoverBlock(ThreadArgs_t * const args)
{
	const OPTIONS * const restrict opt = args->opt;
	coverage_t * const restrict cov = args->cov;
	reference_t * const restrict ref = args->ref;
	const unsigned int ID = args->ID;
	
	unsigned char tag1[kMaxReadLen];
	unsigned char tag2[kMaxReadLen];
	
	/* Global data for realignment */
	GlobalData_t GD[2];
	GD[0].Tag = tag1;
	GD[0].TagLength = kMaxReadLen;
	GD[0].GenomeLength = kGAPPED_ALIGN_GENOME_LENGTH + prf.ExtraLeftGenomeShift;
	GD[1].Tag = tag2;
	GD[1].TagLength = kMaxReadLen;
	GD[1].GenomeLength = kGAPPED_ALIGN_GENOME_LENGTH + prf.ExtraLeftGenomeShift;
	
	cpu_std_border.allocStorage(&GD[0]);
	cpu_std_border.allocStorage(&GD[1]);
	
	/* Which semaphore to wait for ? */
	sem_t * const restrict semaphore = &semThreadStart[ID];
	
	while (1) {
// 		fprintf(stderr, "Thread %u waiting...\n", ID);
		int ret;
		while ((ret = sem_wait(semaphore)) != 0) {
			if (ret != EINTR) goto bail;
		}
		TLBDATA* const restrict td = BlockData[ID];
		
// 		fprintf(stderr, "Thread %u starts working...\n", ID);
		if (td != NULL) {
			TLBHEADER* const restrict th = &BlockHeader[ID];
			{
				void *dataBlock = td->diskBuffer;
				const unsigned int size = th->blockLength - sizeof(TLBHEADER);
				if (decompressBlock(dataBlock,size,th,td,0)) {
					fputs("Error decompressing block\n", stderr);
					goto bail;
				}
			}
			
			unsigned int covSelect = 0;
			unsigned char * restrict diff = td->diff[0];
			if ((th->flags_7_0 & 0xf0) == kTLBHflagNMismatchBlock) 
				diff = NULL; // no need to take those into account - and they break the current decoder anyway
			if ((th->flags_7_0 & 0xf0) == kTLBHflagMismatchBlock) covSelect = 1;
			if ((th->flags_7_0 & 0xf0) == (kTLBHflagMismatchBlock | kTLBHflagCigarBlock)) covSelect = 2;
			if ((td->header.flags_7_0 & kTLBHflagPairedReads)) {
				decodedPair_t dp;
				unsigned int mmpos = 0;
				unsigned int cigarpos = 0;
				// WARNING: discard clones, potentially not correct with threads !!!
				unsigned int cpos1 = 0, cpos2 = 0;
				// just add the coverage
				for (int i=0; i<td->cnt; i++)
				{
					/* Basic Decode pair */
					{
						if (td->header.flags_7_0 & kTLBHflagFixedLength)
						{
							dp.taglen1 = td->lengthInfo;
							dp.taglen2 = td->lengthInfo;
						}
						else
						{
							dp.taglen1 = td->len[2*i];
							dp.taglen2 = td->len[2*i+1];
						}
						char *next = strchr(td->hdr, '\t');
						td->hdr = next + 1;
						next = strchr(td->hdr, '\n');
						td->hdr = next + 1;
						dp.ordinal = 0;
						for (unsigned int j = 0; j < 8; j++)
							dp.ordinal = (dp.ordinal << 8) | (unsigned long) *((unsigned char *) td->hdr++);

						if (td->ppd[0][i].tag2offset & k32bREVCOMP_TAG1_BIT)
							dp.reverseTAG1 = 1;
						else
							dp.reverseTAG1 = 0;

						if (td->ppd[0][i].tag2offset & k32bREVCOMP_TAG2_BIT)
							dp.reverseTAG2 = 1;
						else
							dp.reverseTAG2 = 0;

						dp.delta = (td->ppd[0][i].tag2offset & k32bOFFSET_VALUE_MASK);
						if (td->ppd[0][i].tag2offset & k32bNEGATIVE_OFFSET_BIT)
							dp.delta = -dp.delta;

						dp.genomepos = (virtchr[td->header.chr].chr << 28) + virtchr[td->header.chr].offset + td->ppd[0][i].tag1pos;
					}
					
					// What do we have
#ifdef POSITION_TO_PRINT
					GD[0].MismatchCount = td->ppd[0][i].tag1pos;
					GD[1].MismatchCount = td->ppd[0][i].tag1pos + dp.delta;
#endif
					ExtractStateSequence(&dp, td->qs, (diff != NULL) ? &diff[mmpos] : NULL, td->cigar[0] + cigarpos, GD );

					/* Hack to move pointers */
					cigarpos += GD[0].CigarLength + GD[1].CigarLength;
					mmpos    += GD[0].MismatchLength + GD[1].MismatchLength;
					
					// only keep well-behaved pairs for now
					int ignore = (abs(dp.delta) >= kMAXPAIRLEN || dp.reverseTAG1 == dp.reverseTAG2);
					unsigned int tag2pos = td->ppd[0][i].tag1pos + dp.delta;
					ignore |= (td->ppd[0][i].tag1pos == 0 || td->ppd[0][i].tag1pos > virtchr[opt->chr].len || tag2pos == 0 || tag2pos > virtchr[opt->chr].len);

					// FIXME - maybe should also check same tag lengths...
					ignore |= (!opt->noCloneFilter && td->ppd[0][i].tag1pos == cpos1 && tag2pos == cpos2);
					ignore |= (!opt->noCloneFilter && td->ppd[0][i].tag1pos == cpos2 && tag2pos == cpos1);
					
					cpos1 = td->ppd[0][i].tag1pos;
					cpos2 = tag2pos;
					
					if (!ignore) {
						
#ifdef ORDINAL_TO_PRINT
						if (dp.ordinal == ORDINAL_TO_PRINT) {
							GlobalData_t * ReAlignSettings = GD;
							for (int iTag = 0; iTag<2; iTag++) {
								fprintf(stderr, "TAG %i (len=%zu) @ %i-%i/%i-%i, realign=%i, mismatch=%u, sofclip mismatch=%u, RequiredMM=%u, softclip evicted MM=%u,"
											"RequireCigar=%u, softclip evicted cigar=%u, softclip boundary @ %i\n",
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
								fprintf(stderr, "STA: ");
								for ( int space=0; space<ReAlignSettings[iTag].ProfileRange[0]; space++) fputc(' ', stderr);
								fprintf(stderr, "%s\n", ReAlignSettings[iTag].States);
								fprintf(stderr, "TAG: %.*s", (int) ReAlignSettings[iTag].ProfileRange[0], ReAlignSettings[iTag].Tag);
								int l=0, k=0;
								while (ReAlignSettings[iTag].States[l] != 0) {
									if (ReAlignSettings[iTag].States[l] == 'D') 
										fputc(' ', stderr);
									else 
										fputc((int) ReAlignSettings[iTag].Tag[ReAlignSettings[iTag].ProfileRange[0]+(k++)], stderr);
									l++;
								}
								fprintf(stderr, "\nGEN: "); 
								for ( int space=0; space<ReAlignSettings[iTag].ProfileRange[0]; space++) fputc(' ', stderr);
								l=0; k=0;
								while (ReAlignSettings[iTag].States[l] != 0) {
									if (ReAlignSettings[iTag].States[l] == 'I')
										fputc(' ', stderr);
									else
										fputc(ReAlignSettings[iTag].Genome[ReAlignSettings[iTag].AlignmentRange[0]+(k++)], stderr);
									l++;
								}
								fputc('\n', stderr);
							}
						}
#endif
						
						/* Realign if useful */
						if (GD[0].score > 0) {
							GD[0].TagLength = dp.taglen1;
							GD[0].Genome = genome_tbl + dp.genomepos;
							if (cpos1 >= prf.ExtraLeftGenomeShift) GD[0].Genome -= prf.ExtraLeftGenomeShift;
							
							cpu_std_border.createMatrix(&GD[0]);
							if (cpu_std_border.getStateSequence(&GD[0]) != 0) {
								fprintf(stderr, "ERROR OOPS CANNOT REALIGN TAG 1 of %lu\n", dp.ordinal);
								exit(1);
							}
							cpos1 += (unsigned int) ( (uintptr_t) (GD[0].Genome + GD[0].AlignmentRange[0]) - (uintptr_t) &genome_tbl[dp.genomepos]);
						}
						
						if (GD[1].score > 0) {
							GD[1].TagLength = dp.taglen2;
							GD[1].Genome = genome_tbl + dp.genomepos + dp.delta;
							if (cpos2 >= prf.ExtraLeftGenomeShift) GD[1].Genome -= prf.ExtraLeftGenomeShift;
							
							cpu_std_border.createMatrix(&GD[1]);
							if (cpu_std_border.getStateSequence(&GD[1]) != 0) {
								fprintf(stderr, "ERROR OOPS CANNOT REALIGN TAG 2 of %lu\n", dp.ordinal);
								exit(1);
							}
							cpos2 += (unsigned int) ( (uintptr_t) (GD[1].Genome + GD[1].AlignmentRange[0]) - (uintptr_t) &genome_tbl[dp.genomepos+dp.delta]);
						}
						
#ifdef ORDINAL_TO_PRINT
							if (dp.ordinal == ORDINAL_TO_PRINT) {
								GlobalData_t * ReAlignSettings = GD;
								for (int iTag = 0; iTag<2; iTag++) {
									fprintf(stderr, "TAG %i (len=%zu) @ %i-%i/%i-%i, score=%i, mismatch=%u, sofclip mismatch=%u, RequiredMM=%u, softclip evicted MM=%u,"
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
									fprintf(stderr, "STA: ");
									for ( int space=0; space<ReAlignSettings[iTag].ProfileRange[0]; space++) fputc(' ', stderr);
									fprintf(stderr, "%s\n", ReAlignSettings[iTag].States);
									fprintf(stderr, "TAG: %.*s", (int) ReAlignSettings[iTag].ProfileRange[0], ReAlignSettings[iTag].Tag);
									int l=0, k=0;
									while (ReAlignSettings[iTag].States[l] != 0) {
										if (ReAlignSettings[iTag].States[l] == 'D') 
											fputc(' ', stderr);
										else 
											fputc((int) ReAlignSettings[iTag].Tag[ReAlignSettings[iTag].ProfileRange[0]+(k++)], stderr);
										l++;
									}
									fprintf(stderr, "\nGEN: "); 
									for ( int space=0; space<ReAlignSettings[iTag].ProfileRange[0]; space++) fputc(' ', stderr);
									l=0; k=0;
									while (ReAlignSettings[iTag].States[l] != 0) {
										if (ReAlignSettings[iTag].States[l] == 'I')
											fputc(' ', stderr);
										else
											fputc(ReAlignSettings[iTag].Genome[ReAlignSettings[iTag].AlignmentRange[0]+(k++)], stderr);
										l++;
									}
									fputc('\n', stderr);
								}
							}
#endif

						const unsigned char * restrict quality;
						flags_t flags;
						Cover_t Cover;
						unsigned int fragStart, fragEnd;
						unsigned int LeftTagGenLen, RightTagGenLen;
						unsigned int LeftTagGenPos, RightTagGenPos; 
						
						if (dp.delta >= 0) {
							LeftTagGenLen  = GD[0].AlignmentRange[1] - GD[0].AlignmentRange[0];
							LeftTagGenPos  = cpos1;
							RightTagGenLen = GD[1].AlignmentRange[1] - GD[1].AlignmentRange[0];
							RightTagGenPos = cpos2;
						}
						else {
							LeftTagGenLen  = GD[1].AlignmentRange[1] - GD[1].AlignmentRange[0];
							LeftTagGenPos  = cpos2;
							RightTagGenLen = GD[0].AlignmentRange[1] - GD[0].AlignmentRange[0];
							RightTagGenPos = cpos1;
						}
						
						/* Compute potential tag overlap for too short DNA tags 
							* PROBLEMATIC CASE IS <<<<<<XXXXXXXXXX>>>>>>>
							*              TAG L  <--------------<
							*              TAG R        >--------------->
							*/
						
						if (LeftTagGenPos + LeftTagGenLen >= RightTagGenPos) {
							fragStart = RightTagGenPos;
							fragEnd   = LeftTagGenPos + LeftTagGenLen; 
						}
						else {
							fragStart = LeftTagGenPos;
							fragEnd   = RightTagGenPos + RightTagGenLen;
						}
						
						for (unsigned int k = fragStart; k <= fragEnd; k++) cov[k].cov[covSelect] += 1;
						
						if (dp.reverseTAG1) {
							flags.reverse = 1;
							quality = &td->qs[dp.taglen1 - 1];
						} 
						else {
							flags.reverse = 0;
							quality = td->qs;
						}
						flags.secondRead = 0;
#ifndef NDEBUG
						/*
						if (cpos1 < cpos2) {
							if (cpos2-cpos1 >= kMAXPAIRLEN) {
								fprintf(stderr, "Potential miss, distance between tags to short, %u > %u\n", cpos2-cpos1, kMAXPAIRLEN);
								exit(1);
							}
						}
						else {
							if (cpos1-cpos2 >= kMAXPAIRLEN) {
								fprintf(stderr, "Potential miss, distance between tags to short, %u > %u\n", cpos1-cpos2, kMAXPAIRLEN);
								exit(1);
							}
						}*/
#endif
						MyCover(cov, ref, dp.ordinal, (char*) GD[0].States, (char*) GD[0].Tag, quality, fragStart, fragEnd, &Cover, cpos1, flags, opt
#ifdef POSITION_TO_PRINT
						        ,ID
#endif
						);
						td->qs += dp.taglen1;
						
#ifdef BYPASS_COMMON
						Cover.min = 0xFFFFFFFFU;
						Cover.max = 0U;
#endif
						if (dp.reverseTAG2) {
							flags.reverse = 1;
							quality = &td->qs[dp.taglen2 - 1];
						} 
						else {
							flags.reverse = 0;
							quality = td->qs;
						}
						flags.secondRead = 1;
						MyCover(cov, ref, dp.ordinal, (char*) GD[1].States,(char*) GD[1].Tag, quality, fragStart, fragEnd, &Cover, cpos2, flags, opt
#ifdef POSITION_TO_PRINT
						        ,ID
#endif
						);
						td->qs += dp.taglen2;
						
					}
					else {
						/* move pointers */
						td->qs += dp.taglen1;
						td->qs += dp.taglen2;
					}
				}
			}
		}
		else {
			break;
		}

// 		fprintf(stderr, "Thread %u done...\n", ID);
		sem_post(&semThreadDone[ID]);
		
	}
// 	fprintf(stderr, "Thread %u, bye bye...\n", ID);
	
	return NULL;
	bail: ;
		exit(1);
}
//---------------------------------------------------------------

int
decompress_PairedThreaded(const char * const restrict fn, OPTIONS *opt, coverage_p_t cov, reference_p_t ref)
{
	if (nThreads & 0x1) {
		fprintf(stderr, "Only even number of threads is accepted, increasing to %u\n", ++nThreads);
	}
	
	int fd = open(fn, O_RDONLY);
	int err = 1;
	ThreadArgs_t args[nThreads];
	pthread_t threads[nThreads];
	pthread_attr_t thread_attr[nThreads];
	cpu_set_t thread_cpu_set[nThreads];

	BlockData = (TLBDATA**) alloca(nThreads*sizeof(TLBDATA*));
	for (unsigned int i=0;i<nThreads; i++) {
		BlockData[i] = (TLBDATA*) malloc(sizeof(TLBDATA));
		if (BlockData[i] == NULL) goto bail;
		allocBlock(BlockData[i], true);
	}
	
	BlockHeader = (TLBHEADER*) malloc(nThreads*sizeof(TLBHEADER));
	if (BlockHeader == NULL) {
		perror("malloc error");
		goto bail;
	}
	
	sem_t * sems = (sem_t*) alloca(2*nThreads*sizeof(sem_t));
	for (unsigned int i=0;i<(2*nThreads); i++) {
		if (sem_init(&sems[i], 0, 0)) {
			perror("sem_init");
			goto bail;
		}
	}
	
	semThreadStart = sems;
	semThreadDone = &sems[nThreads];
	
	/* start threads */
	void* (*WhichCover)(void*) = (void* (*)(void*)) ( opt->Realign  ? RealignCoverBlock : CoverBlock );
	const int nproc = (int) sysconf(_SC_NPROCESSORS_CONF);
	for (unsigned int i=0;i<nThreads; i++) {
		args[i].cov = cov;
		args[i].ref = ref;
		args[i].opt = opt;
		args[i].ID = i;
		pthread_attr_init(&thread_attr[i]);
		CPU_ZERO(&thread_cpu_set[i]);
		CPU_SET((int) (i) % nproc, &thread_cpu_set[i]);
		pthread_attr_setaffinity_np(&thread_attr[i], sizeof(cpu_set_t), &thread_cpu_set[i]);
		
		if (pthread_create(&threads[i], &thread_attr[i], WhichCover, (void*) &args[i])) {
			perror("pthread_create");
			for (unsigned int j=0; j<i; j++) {
				pthread_kill(threads[j], SIGKILL);
			}
			goto bail;
		}
	}
	
	const int MinDistance = prf.ExtraLeftGenomeShift + kGAPPED_ALIGN_GENOME_LENGTH + 2*kMAXPAIRLEN;
	unsigned int BlockCounter = 0U;
	while (1) {
		ssize_t res = read(fd,&BlockHeader[BlockCounter],sizeof(TLBHEADER));
		if (res == 0)
		{
// 			fprintf(stderr, "File fully read\n");
			if (BlockCounter) {
// 				fprintf(stderr, "Sending last even batch\n");
				for (unsigned int i=0; i<BlockCounter; i+=2) {
					if (sem_post(&semThreadStart[i])) {
						perror("sem_post");
						goto bail;
					}
				}
				
// 				fprintf(stderr, "Waiting for even threads to terminate\n");
				for (unsigned int i=0; i<BlockCounter; i+=2) {
					int ret;
					while ( (ret = sem_wait(&semThreadDone[i])) != 0) {
						if (ret != EINTR) goto bail;
					}
				}
				
// 				fprintf(stderr, "Sending last odd batch\n");
				for (unsigned int i=1; i<BlockCounter; i+=2) {
					if (sem_post(&semThreadStart[i])) {
						perror("sem_post");
						goto bail;
					}
				}
				
// 				fprintf(stderr, "Waiting for odd threads to terminate\n");
				for (unsigned int i=1; i<BlockCounter; i+=2) {
					int ret;
					while ( (ret = sem_wait(&semThreadDone[i])) != 0) {
						if (ret != EINTR) goto bail;
					}
				}
			}
			
			/* Ending the remaining threads */
			{
				for (unsigned int i=0; i<nThreads; i++) {
					freeBlock(BlockData[i]);
					BlockData[i] = NULL;
					if (sem_post(&semThreadStart[i])) {
						perror("sem_post");
						goto bail;
					}
				}
				/* wait for completion */
				void * ret;
				for (unsigned int i=0; i<nThreads; i++) pthread_join(threads[i], &ret);
			}
			
			err = 0;
			break;
		}
		if (res != sizeof(TLBHEADER))
			goto bail;
		unsigned char hcs = simple8bitCS((unsigned char *) &BlockHeader[BlockCounter], sizeof(TLBHEADER));
		if (hcs != 0)
		{
			fprintf(stderr,"Checksum mismatch: %02x %02x\n",BlockHeader[BlockCounter].headerCS,hcs);
			goto bail;
		}
		// do we want this block ?
		if (BlockHeader[BlockCounter].chr == opt->chr)
		{
			if ((BlockHeader[BlockCounter].flags_7_0 & kTLBHflagPairPosBlock) != 0)
			{
				// we count the coverage of Ns in the reads as if they were reference (i.e., exact matches)
				if ( (BlockHeader[BlockCounter].flags_7_0 & 0xf0) == 0
					|| (BlockHeader[BlockCounter].flags_7_0 & 0xf0) == kTLBHflagNMismatchBlock
					|| (BlockHeader[BlockCounter].flags_7_0 & 0xf0) == kTLBHflagMismatchBlock
					|| (BlockHeader[BlockCounter].flags_7_0 & 0xf0) == (kTLBHflagMismatchBlock | kTLBHflagCigarBlock))
				{
// 					fprintf(stderr, "Filling in block %u\n", BlockCounter);
					err = readBlock(fd, &BlockHeader[BlockCounter], BlockData[BlockCounter]);
					if (err != 0)
						goto bail;
					
// 					if (BlockData[BlockCounter]->header.flags_7_0 & kTLBHflagPairedReads) {
					if (BlockHeader[BlockCounter].flags_7_0 & kTLBHflagPairedReads) {
						const int CurrentminPos = BlockHeader[BlockCounter].minPos;
						const int CurrentmaxPos = BlockHeader[BlockCounter].maxPos;
						/* Check that we will not have overlap */
						for (unsigned int i=0; i<BlockCounter; i++) {
							const int mindistance = CurrentminPos - BlockHeader[i].maxPos;
							const int maxdistance = BlockHeader[i].minPos - CurrentmaxPos;
							if (abs(mindistance) < MinDistance || abs(maxdistance) < MinDistance) {
								if (!((i+BlockCounter) & 0x1))  {
									fprintf(stderr, "Blocks %u and %u are separated by only %i or %i\n", i, BlockCounter, mindistance, maxdistance);
									goto PerformThese;
								}
							}
						}
						if (++BlockCounter == nThreads) {
						PerformThese: ;
// 						fprintf(stderr, "Starting even threads\n");
							for (unsigned int i=0; i<BlockCounter; i+=2) {
								if (sem_post(&semThreadStart[i])) {
									perror("sem_post");
									goto bail;
								}
							}
							
							/* Wait for the thread to terminate */
// 							fprintf(stderr, "Waiting for even threads to terminate\n");
							for (unsigned int i=0; i<BlockCounter; i+=2) {
								int ret;
								while ( (ret = sem_wait(&semThreadDone[i])) != 0) {
									if (ret != EINTR) goto bail;
								}
							}
							
// 							fprintf(stderr, "Starting odd threads\n");
							for (unsigned int i=1; i<BlockCounter; i+=2) {
								if (sem_post(&semThreadStart[i])) {
									perror("sem_post");
									goto bail;
								}
							}
							
							/* Wait for the thread to terminate */
// 							fprintf(stderr, "Waiting for odd threads to terminate\n");
							for (unsigned int i=1; i<BlockCounter; i+=2) {
								int ret;
								while ( (ret = sem_wait(&semThreadDone[i])) != 0) {
									if (ret != EINTR) goto bail;
								}
							}
							
							/* Move last problematic one to first position, swapping them */
							if (BlockCounter != nThreads) {
								TLBDATA * tdata = BlockData[BlockCounter];
								BlockData[BlockCounter] = BlockData[0];
								BlockData[0] = tdata;
								memcpy(&BlockHeader[0], &BlockHeader[BlockCounter], sizeof(TLBHEADER));
								BlockCounter = 1U;
							}
							else {
								BlockCounter = 0U;
							}
						}
					}
				}
			}
		}
		else
		{
			// nope, read and discard (so that we can stream...)
			size_t size = (size_t) BlockHeader[BlockCounter].blockLength - sizeof(TLBHEADER);
			void *data = malloc(size);
			ssize_t res = read(fd,data,size);
			if (res != size)
			{
				perror("skip read:");
				goto bail;
			}
			free(data);
		}
	}

	err = 0;

bail:;

	// clean up
	close(fd);
	for (unsigned int i=0;i<nThreads; i++) if (BlockData[i]) freeBlock(BlockData[i]);

	return(err);

} // decompress_realign
//---------------------------------------------------------------

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
