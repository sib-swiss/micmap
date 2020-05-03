/*
 * ------------------------------------------------------------------------------------------------------------------------
 *
 *                                * micmap *
 *      Mapping of short reads in fastq format onto a reference
 *
 *
 *  Copyright (C) SIB  - Swiss Institute of Bioinformatics,  2015-2019 Nicolas Guex, Thierry Schuepbach and Christian Iseli
 *  Copyright (C) UNIL - University of Lausanne, Switzerland 2019-2020 Nicolas Guex and Christian Iseli
 *  Copyright (C) EPFL - Lausanne, Switzerland                    2020 Christian Iseli
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

//=================================================================================================
/*

	This software will analyze a GTL file produced by MicMap using the provided
	genome reference and produce a count of nt per position

*/

/* test:


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
#include <errno.h>

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

#define kMAXPAIRLEN 800
#define kMaxLineBuf 1024

//=================================================================================================

//static struct timeval ts;
//static struct timeval te;
static int verbose = 0;
static int debug = 0;
static VIRTUALCHR virtchr[kMAXCHRcnt];

typedef struct OPTIONS_struct {
	int *modTable;
	char *table;
	unsigned int mod;
	int chr;
	int cloneFilter;
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
	int ilen;
} decodedPair_t, *decodedPair_p_t;

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
decodeOrdinalSingle(TLBDATA *tpp, decodedSingle_p_t dsp)
{
	char *next = strchr(tpp->hdr, '\n');

	tpp->hdr = next + 1;
	dsp->ordinal = 0;
	for (unsigned int i = 0; i < 8; i++)
		dsp->ordinal = (dsp->ordinal << 8) | (unsigned long) *((unsigned char *) tpp->hdr++);
}
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
	decodeOrdinalSingle(tpp,dsp);

	if (tpp->ppd[0][i].tag2offset & k32bREVCOMP_TAG1_BIT)
		dsp->reverseTAG1 = 1;
	else
		dsp->reverseTAG1 = 0;

	dsp->genomepos = (virtchr[tpp->header.chr].chr << 28) + virtchr[tpp->header.chr].offset + tpp->ppd[0][i].tag1pos;
}
//---------------------------------------------------------------

static void
decodeOrdinalPair(TLBDATA *tpp, decodedPair_p_t dpp)
{
	char *next = strchr(tpp->hdr, '\t');

	tpp->hdr = next + 1;
	next = strchr(tpp->hdr, '\n');
	tpp->hdr = next + 1;
	dpp->ordinal = 0;
	for (unsigned int i = 0; i < 8; i++)
		dpp->ordinal = (dpp->ordinal << 8) | (unsigned long) *((unsigned char *) tpp->hdr++);
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
	decodeOrdinalPair(tpp,dpp);

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

static int
isSelectedOrdinal(unsigned long ordinal, OPTIONS *opt)
{
	if (opt->mod < 2)
		return 1;
	unsigned int r = ordinal % opt->mod;
	return opt->modTable[r];
}
//---------------------------------------------------------------

static int
decompress(char *fn, OPTIONS *opt)
{
	int fd = open(fn, O_RDONLY);
	TLBHEADER th;
	unsigned int *cov = calloc(virtchr[opt->chr].len+1,sizeof(unsigned int));
	unsigned int genomepos = (virtchr[opt->chr].chr << 28) + virtchr[opt->chr].offset;
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
					err = readDecompressBlock(fd, &th, &td, 2);
					if (err != 0)
						goto bail;
					unsigned char *diff = td.diff[0];
					if ((th.flags_7_0 & 0xf0) == kTLBHflagNMismatchBlock)
						diff = NULL; // no need to take those into account - and they break the current decoder anyway
					if ((td.header.flags_7_0 & kTLBHflagPairedReads) == 0)
					{
						decodedSingle_t ds;
						int i;
						// just add the coverage
						for (i=0; i<td.cnt; i++)
						{
							decodeSingle(&td,i,&ds);
							if (isSelectedOrdinal(ds.ordinal, opt))
							{
								// we have some rare cases where a read starts mapping "before" the beginning of the chromosome (see badPos fastq examples)
								int ignore = 0;
								if (td.ppd[0][i].tag1pos == 0 || td.ppd[0][i].tag1pos > virtchr[opt->chr].len)
									ignore = 1;
								// FIXME - maybe should also check same tag lengths...
								if (opt->cloneFilter && td.ppd[0][i].tag1pos == cpos1)
									ignore = 1;
								cpos1 = td.ppd[0][i].tag1pos;
								if (!ignore)
								{
									for (unsigned int j = 0; j < ds.taglen1; j++)
										cov[td.ppd[0][i].tag1pos + j] += 1;
								}
							}
							//td.qs += ds.taglen1;
						}
					}
					else
					{
						decodedPair_t dp;
						int i;
						// just add the coverage
						for (i=0; i<td.cnt; i++)
						{
							// NOTE - we have some rare cases where a read starts mapping "before" the beginning of the chromosome (see badPos fastq examples)
							decodePair(&td,i,&dp);
							if (isSelectedOrdinal(dp.ordinal, opt))
							{
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
								if (!ignore)
								{
									unsigned int start = td.ppd[0][i].tag1pos;
									unsigned int end = td.ppd[0][i].tag1pos + dp.taglen1;
									if (tag2pos < start)
										start = tag2pos;
									if (tag2pos + dp.taglen2 > end)
										end = tag2pos + dp.taglen2;
									for (unsigned int j = start; j < end; j++)
										cov[j] += 1;
								}
							}
							//td.qs += dp.taglen1;
							//td.qs += dp.taglen2;
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
			void *rp = data;
			ssize_t rs = size;
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
					fprintf(stderr,"unexpected EOF while reading block in decompress\n");
					goto bail;
				}
				if (res == -1 && (errno == EINTR || errno == EAGAIN))
					continue;
				perror("decompress read:");
				goto bail;
			}
			free(data);
		}
	}

	err = 0;

	if (opt->table != NULL)
	{
		FILE *f = fopen(opt->table, "r");
		char buf[512];
		while (fgets(buf, 512, f))
		{
			unsigned int len = strlen(buf);
			// remove newline
			len -= 1;
			buf[len] = 0;
			unsigned int start, end;
			if (sscanf(buf,"%*s\t%u\t%u\t",&start,&end) == 2)
			{
				unsigned int sum = 0;
				for (unsigned int i = start; i <= end; i++)
					sum += cov[i];
				printf("%s\t%.2f\n",buf,(double) sum / (double) (end - start + 1));
			}
		}
		fclose(f);
	}
	else
	{
		for (unsigned int i = 1; i <= virtchr[opt->chr].len; i++)
			printf("%u\n",cov[i]);
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
parseSelection(OPTIONS *opt, const char *s)
{
	char *e;
	unsigned long int v = strtoul(s, &e, 0);
	if (v < 2 || v > 65536)
	{
		fprintf(stderr,"Bad modulus value %lu\n", v);
		exit(1);
	}
	if (*e != ',')
	{
		fprintf(stderr,"Expected , after modulus is not present : %s\n", e);
		exit(1);
	}
	opt->mod = v;
	opt->modTable = calloc(v, sizeof(int));
	while (*e == ',')
	{
		v = strtoul(e + 1, &e, 0);
		if (v >= opt->mod)
		{
			fprintf(stderr,"Value %lu is out of the %u modulus range\n", v, opt->mod);
			exit(1);
		}
		if (*e != 0 && *e != ',')
		{
			fprintf(stderr,"Expected , or end of string is not present : %s\n", e);
			exit(1);
		}
		opt->modTable[v] = 1;
	}
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

	/* --------- process arguments ----------------------*/

	rsltfile[0] = 0;
	strcpy(genomefile,"/tmp/nodelete/nguex/data/hg19.bin");
	opterr = 0;
	while ((c = getopt (argc, argv, "r:g:C:v:dct:s:")) != -1)
	switch (c)
	{
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

		case 'd':
			debug = 1;
			break;

		case 'v':
			sscanf(optarg,"%d",&verbose);
			break;

		case 't':
			options.table = optarg;
			break;

		case 's':
			parseSelection(&options,optarg);
			break;
	}

	if (rsltfile[0] == 0 || options.chr == -1)
	{
		printf("usage:\n\n");
		printf("GTLsubCoverage [options] -C chromosome -r resultfile\n\n");
		printf("           -g genomefile         : full path of the genome encoded by encode_genome\n");
		printf("                                   default: '%s' \n",genomefile);
		printf("           -C chromosome         : select chromosome to analyze\n");
		printf("           -r resultfile         : name of the result file (produced by match)\n");
		printf("           -c                    : do filter clones (default is to not filter them out since we handle amplicons)\n");
		printf("           -d                    : debug mode...\n");
		printf("           -v level              : verbose level\n");
		printf("           -t tablefile          : path to a table of regions to obtain coverage from (columns 2 and 3 are start/end)\n");
		printf("           -s selection          : select which group of ordinals is included in the counts\n");
		printf("\n");
		printf("selection : <modulus>,<elements of modulus>,...\n");
		printf("e.g., 10,0,1 will do 'ordinal % 10' and keep if the remainder is 0 or 1\n");
		printf("\n");
		printf("(c) Nicolas Guex & Christian Iseli 2014-2020\n");
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
/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
