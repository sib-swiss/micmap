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

//=================================================================================================

#define kMAXPAIRLEN 512
#define kMaxLineBuf 1024
#define kMAXPWMPOS 1024

//=================================================================================================

//static struct timeval ts;
//static struct timeval te;
static int verbose = 0;
static int debug = 0;
static VIRTUALCHR virtchr[kMAXCHRcnt];

typedef struct POS_struct {
	unsigned int chr;
	unsigned int pos;
} POS;

typedef struct OPTIONS_struct {
	int chr;
	int cloneFilter;
	unsigned int filterQualValue;
	unsigned int PWMcnt;
	unsigned int minCov;
	POS PWMpos[kMAXPWMPOS];
	char outbase[512];
} OPTIONS;

typedef struct _decodedSingle {
	unsigned int genomepos;
	int taglen1;
	int reverseTAG1;
} decodedSingle_t, *decodedSingle_p_t;

typedef struct _decodedPair {
	unsigned int genomepos;
	int taglen1;
	int taglen2;
	int reverseTAG1;
	int reverseTAG2;
	int delta;
	int ilen;
} decodedPair_t, *decodedPair_p_t;

typedef struct _coverage {
	unsigned int nt[4];
} coverage_t, *coverage_p_t;

unsigned char gIUPACrevcomp[128];

//=================================================================================================



static 	unsigned char *genome_tbl;	// contain genome in 16 virtual chromosomes.
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
			printf("genomeStat.st_size %lx != %lx\n",genomeStat.st_size,kGENOME_DATA_SIZE);
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

	if (dpp->delta <= 0) {
		dpp->ilen = dpp->delta + dpp->taglen2;
	} else {
		dpp->ilen = dpp->delta - dpp->taglen1;
	}
}
//---------------------------------------------------------------

#define kNOT_ACGT 0xff
static unsigned char gNT2bits[128];

static void
cover(coverage_p_t cov, unsigned int tagPos, unsigned int tagLen, int reverse, unsigned char *diff, unsigned int *mmpos, unsigned char *cigar, unsigned int *cigarpos, const unsigned char *qual, int ignore, int secondRead, OPTIONS *opt, unsigned int genomepos)
{
	// we keep the min and max pos affected here; when secondRead is true, we avoid counting coverage a second time
	static unsigned int minCovPos;
	static unsigned int maxCovPos;
	unsigned int j;
	unsigned int pos = tagPos;
	unsigned int cp = *cigarpos;
	unsigned char cc[tagLen*2]; // the expanded cigar codes
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
			while(cpos-- >0)
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
		unsigned int rpos = 0;
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
							if (qual[tagLen - rpos - 1] - '#' >= opt->filterQualValue)
								cov[pos].nt[gNT2bits[genome_tbl[genomepos+pos]]] += 1;
							pos += 1;
						}
						else
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
							if (qual[rpos] - '#' >= opt->filterQualValue)
								cov[pos].nt[gNT2bits[genome_tbl[genomepos+pos]]] += 1;
							pos += 1;
						}
						else
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
				fprintf(stderr,"HWHAP: cigarcode should be M at %u, but is %c\n",dpos,cct[dpos]);
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
					if (q >= opt->filterQualValue)
					{
						cov[gpos].nt[gNT2bits[genome_tbl[genomepos+gpos]]] -= 1;
						cov[gpos].nt[gNT2bits[(unsigned char) nt]] += 1;
					}
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
						if (secondRead == 0)
						{
							if (gpos < minCovPos)
								minCovPos = gpos;
							if (gpos > maxCovPos)
								maxCovPos = gpos;
						}
						gpos += cpos;
						break;
					}
					case 'I':
					{
						if (secondRead == 0)
						{
							if (gpos < minCovPos)
								minCovPos = gpos;
							if (gpos > maxCovPos)
								maxCovPos = gpos;
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
outputPWM(OPTIONS *opt, coverage_p_t cov, unsigned int startPos, unsigned int endPos, unsigned int max)
{
	char fn[512];
	sprintf(fn,"%schr%u_%u_%u_%u.pwm",opt->outbase,opt->chr,startPos,endPos,max);
	FILE *f = fopen(fn,"w");
	if (!f)
	{
		fprintf(stderr,"Error: Cannot create %s\n",fn);
		perror("fopen");
		exit(1);
	}
	unsigned int i;
	fprintf(f,"PO\tA\tC\tG\tT\n");
	for (i = startPos; i <= endPos; i++)
		fprintf(f,"%02u\t%u\t%u\t%u\t%u\n",i - startPos + 1,cov[i].nt[0],cov[i].nt[1],cov[i].nt[2],cov[i].nt[3]);
	fclose(f);
}
//---------------------------------------------------------------

static int
decompress(char *fn, OPTIONS *opt)
{
	int fd = open(fn, O_RDONLY);
	TLBHEADER th;
	coverage_p_t cov = calloc(virtchr[opt->chr].len+1,sizeof(coverage_t));
	unsigned int genomepos = (virtchr[opt->chr].chr << 28) + virtchr[opt->chr].offset;
	memset(gNT2bits,kNOT_ACGT,128);
	gNT2bits['A'] = 0x0;
	gNT2bits['C'] = 0x1;
	gNT2bits['G'] = 0x2;
	gNT2bits['T'] = 0x3;
	int err = 1;
	// discard clones
	unsigned int cpos1 = 0, cpos2 = 0;
	TLBDATA td;
	allocBlock(&td, true);
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
					unsigned char *diff = td.diff[0];
					if ((th.flags_7_0 & 0xf0) == kTLBHflagNMismatchBlock)
						diff = NULL; // no need to take those into account - and they break the current decoder anyway
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
							if (opt->cloneFilter && td.ppd[0][i].tag1pos == cpos1)
								ignore = 1;
							cpos1 = td.ppd[0][i].tag1pos;
							cover(cov,td.ppd[0][i].tag1pos, ds.taglen1, ds.reverseTAG1, diff, &mmpos, td.cigar[0], &cigarpos, td.qs, ignore, 0, opt, genomepos);
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
							// only keep well-behaved pairs for now
							int ignore = 0;
							if (abs(dp.delta) >= kMAXPAIRLEN || dp.reverseTAG1 == dp.reverseTAG2)
								ignore = 1;
							unsigned int tag2pos = td.ppd[0][i].tag1pos + dp.delta;
							if (td.ppd[0][i].tag1pos == 0 || td.ppd[0][i].tag1pos > virtchr[opt->chr].len
									|| tag2pos == 0 || tag2pos > virtchr[opt->chr].len)
								ignore = 1;
							// FIXME - maybe should also check same tag lengths...
							if (opt->cloneFilter && td.ppd[0][i].tag1pos == cpos1 && tag2pos == cpos2)
								ignore = 1;
							if (opt->cloneFilter && td.ppd[0][i].tag1pos == cpos2 && tag2pos == cpos1)
								ignore = 1;
							cpos1 = td.ppd[0][i].tag1pos;
							cpos2 = tag2pos;
							cover(cov,td.ppd[0][i].tag1pos, dp.taglen1, dp.reverseTAG1, diff, &mmpos, td.cigar[0], &cigarpos, td.qs, ignore, 0, opt, genomepos);
							td.qs += dp.taglen1;
							cover(cov,tag2pos, dp.taglen2, dp.reverseTAG2, diff, &mmpos, td.cigar[0], &cigarpos, td.qs, ignore, 1, opt, genomepos);
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
	unsigned int i;

	int on = 0;
	unsigned int startPos;
	for (i = 0; i <= virtchr[opt->chr].len; i++)
	{
		unsigned int sum = cov[i].nt[0] + cov[i].nt[1] + cov[i].nt[2] + cov[i].nt[3];
		// look for the amplicon regions
		if (on)
		{
			if (sum == 0)
			{
				unsigned int j;
				unsigned int max = 0;
				unsigned int sum;
				on = 0;
				for (j = startPos; j < i; j++)
				{
					sum = cov[j].nt[0] + cov[j].nt[1] + cov[j].nt[2] + cov[j].nt[3];
					if (sum > max)
						max = sum;
				}
				sum = cov[startPos].nt[0] + cov[startPos].nt[1] + cov[startPos].nt[2] + cov[startPos].nt[3];
				while (sum * 2 < max)
				{
					startPos += 1;
					sum = cov[startPos].nt[0] + cov[startPos].nt[1] + cov[startPos].nt[2] + cov[startPos].nt[3];
				}
				unsigned int endPos = i - 1;
				sum = cov[endPos].nt[0] + cov[endPos].nt[1] + cov[endPos].nt[2] + cov[endPos].nt[3];
				while (sum * 2 < max)
				{
					endPos -= 1;
					sum = cov[endPos].nt[0] + cov[endPos].nt[1] + cov[endPos].nt[2] + cov[endPos].nt[3];
				}
				printf("%u\t%u\t%u\t%u\n",opt->chr,startPos,endPos,max);
				// see if we have holes or SNPs
				unsigned int minVar = (sum * 15) / 100; // 15 % change
				int goingDown = 0;
				unsigned int prevSum = 0;
				int reportMe = 0;
				for (j = startPos; j <= endPos; j++)
				{
					unsigned int k;
					unsigned int cntVar = 0;
					sum = 0;
					for (k = 0; k < 4; k++)
					{
						sum = cov[j].nt[k];
						if (cov[j].nt[k] >= minVar)
							cntVar += 1;
					}
					if (cntVar > 1)
					{
						fprintf(stderr,"%u\t%u\tSNP\n",opt->chr,j);
						reportMe = 1;
					}
					if (sum + minVar <= prevSum)
						goingDown = 1;
					if (sum >= prevSum + minVar && goingDown)
					{
						fprintf(stderr,"%u\t%u\tDEL\n",opt->chr,j);
						reportMe = 1;
					}
				}
				if (reportMe)
					outputPWM(opt, cov, startPos, endPos, max);
			}
		}
		else
		{
			if (sum >= opt->minCov)
			{
				on = 1;
				startPos = i;
			}
		}
	}

bail:;

	// clean up
	close(fd);
	freeBlock(&td);
	free(cov);

	return(err);

} // decompress
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
	register const char* const dummy = fgets(linbuf,kMaxLineBuf,f); // skip header.
	while(!feof(f))
	{
		char tmpchr[32];
		register const char* const dummy2 = fgets(linbuf,kMaxLineBuf,f);
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
	int c;
	char rsltfile[512];
	char genomefile[512];
	int err;
	OPTIONS options;
	memset(&options,0,sizeof(OPTIONS));
	options.minCov = 100;
	strcpy(options.outbase,"/tmp/");

	/* --------- process arguments ----------------------*/

	rsltfile[0] = 0;
	strcpy(genomefile,"/tmp/nodelete/nguex/data/hg19.bin");
	opterr = 0;
	while ((c = getopt (argc, argv, "r:g:C:v:dp:o:cf:m:")) != -1)
	switch (c)
	{
		case 'p':
			loadPos(optarg,&options);
			break;

		case 'o':
			strcpy(options.outbase,optarg);
			break;

		case 'r':
			strcpy(rsltfile,optarg);
			break;

		case 'C':
			sscanf(optarg,"%d",&options.chr);
			break;

		case 'g':
			strcpy(genomefile,optarg);
			break;

		case 'c':
			options.cloneFilter = 1;
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

		case 'm':
			sscanf(optarg,"%d",&options.minCov);
			break;
	}

	if (rsltfile[0] == 0 || options.chr == -1)
	{
		printf("usage:\n\n");
		printf("GTLcaller [options] -C chromosome -r resultfile\n\n");
		printf("           -g genomefile         : full path of the genome encoded by encode_genome\n");
		printf("                                   default: '%s' \n",genomefile);
		printf("           -C chromosome         : select chromosome to analyze\n");
		printf("           -r resultfile         : name of the result file (produced by match)\n");
		printf("           -c                    : do filter clones (default is to not filter them out since we handle amplicons)\n");
		printf("           -f <char>             : filter out nucleotides where the quality value is *below* the specified character [%c]\n", options.filterQualValue + '#');
		printf("           -d                    : debug mode...\n");
		printf("           -m                    : minimum coverage to report [%u]\n",options.minCov);
		printf("           -o <path>             : base part of output files for PWM [%s]\n",options.outbase);
		printf("           -p <file>             : generate PWM output (for weblogo) based on list of chr:pos found in supplied file\n");
		printf("           -v level              : verbose level\n");
		printf("\n");
		printf(GTL_VERSION_FULL_STRING "\n");
		printf("(c) Nicolas Guex & Christian Iseli 2014-2015\n");
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

	err = decompress(rsltfile,&options);

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
