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

  This software will encode a list of tags that have exact unique matches against a genome reference
  Tags must be 18nt long.
  

  It will produce four tables: start0.A.2nt.bin start0.A.chr.bin start0.A.strand.bin start0.A.valid.bin
  
  1)  a 4Gb table which contains the last 2nt of the tag at the memory address 
      pointed by the 32bits encoding of the 16 first nt.

  2)  a 16Gb table which contains the chr & position of the tag encoded on 16 virtual chromosomes
	  pointed by the 32bits encoding of the 16 first nt.

  3)  a 500Mb table which contains the orientation of each encoded tag with respect to the reference genome.

  4)  a 500Mb table which specifies which tags are valid.


	(c) N.Guex and C.Iseli 2014 & 2015

*/
#define _GNU_SOURCE
#include "config.h"
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <sys/time.h>
#include <zlib.h>
#include "virt_chr.h"
#include "GTLkonst.h"
#include "constants.h"
#include "MATCHkonst.h"
#include "gtlVersion.h"
#define kmaxPAR 16

/*
        cc -m64 -O3 -o encode_UEM encode_UEM.c virt_chr.c -Wall  
*/
//---------------------------------------------------------------

static VIRTUALCHR virtchr[kMAXCHRcnt];
static int verbose = 0;
typedef struct _chrReg_t
{
	unsigned int chr;
	unsigned int start;
	unsigned int end;
} chrReg_t;
typedef struct _PARinfo_t
{
	unsigned int cnt;
	chrReg_t PAR[kmaxPAR];
} PARinfo_t;
typedef struct _params_t
{
	unsigned int *micdata;
	unsigned int *hostdata;
	unsigned char *stranddata;
	unsigned char *validdata;
	PARinfo_t *pi;
	int startingNnt;
	char NT;
} PARAMS;
typedef struct _tagEntry_t
{
	char tag[32];
	char AC[32];
	unsigned int chr;
	unsigned int pos;
	char strand;
} tagEntry_t;
static int whichNntStat[kFollowNtCount];

//---------------------------------------------------------------

void write2NTtable(unsigned int *data,size_t datasize,char *fn)
{
	FILE *f = NULL;

		f = fopen(fn,"wb");
		if (f == NULL)
		{
			printf("Error writing %s\n",fn);
			return;
		}
		printf("writing 2NT table\n");
		fwrite(data,sizeof(unsigned int),datasize,f);
		fclose(f);

	
} /* write2NTtable */
//---------------------------------------------------------------

void writeCHRPOStable(unsigned int *data,size_t datasize,char *fn)
{
	FILE *f = NULL;

	f = fopen(fn,"wb");
	if (f == NULL)
	{
		printf("Error writing %s\n",fn);
		return;
	}
	printf("writing CHRPOS table\n");
	fwrite(data,sizeof(unsigned int),kFollowNtCount*datasize,f);
	fclose(f);

} /* writeCHRPOStable */
//---------------------------------------------------------------

void writeStrandOrValidTable(unsigned char *data,size_t datasize,char *fn)
{
	FILE *f = NULL;

	f = fopen(fn,"wb");
	if (f == NULL)
	{
		printf("Error writing %s\n",fn);
		return;
	}
	printf("writing strand or valid table\n");
	fwrite(data,sizeof(unsigned char),datasize,f);
	fclose(f);

} /* writeStrandOrValidTable */
//---------------------------------------------------------------

static inline void
readTag(gzFile f, tagEntry_t *te)
{
	char buffer[256];
	unsigned int c;
	buffer[0] = 0;
	te->pos = 0;
	char *res = gzgets(f,buffer,256);
	if (res == NULL || buffer[0] == 0) {
		te->tag[0] = 0;
		te->chr = 0;
		return;
	}
	sscanf(buffer,"%s\t%*s\t%[A-Za-z0-9_.][%u..%*u]\t%c",te->tag,te->AC,&te->pos,&te->strand);
	for (c = 1; c <= chrMax; c++)
		if (strcmp(te->AC,virtchr[c].AC) == 0)
			break;
	if (c > chrMax)
	{
		fprintf(stderr,"Unknown AC %s\n",te->AC);
		exit(1);
	}
	te->chr = c;
	if (verbose)
		printf("got %s %s %d %d %c\n",te->tag,te->AC,te->chr,te->pos,te->strand);
}
//---------------------------------------------------------------

static void
encodeEntry(tagEntry_t *te, unsigned int whichNnt, unsigned int *lineok, PARAMS *prm)
{
	unsigned int i;
	unsigned int tagbin = 0;
	unsigned int nextNbin = 0;
	unsigned int *ptr;

	if ((whichNnt < prm->startingNnt) || (whichNnt >= (prm->startingNnt+kFollowNtCount)) || (te->tag[0] != prm->NT)) 
	{ 
		if (verbose > 1)
			printf("whichNnt=%4d ignoring %s\n",whichNnt,te->tag); 
		return;
	}

	for (i = 0; i < 16; i++) 
	{
		tagbin <<= 2;
		switch(te->tag[i]) 
		{
			case 'A':
				break;
			case 'C':
				tagbin |= 1;
				break;
			case 'G':
				tagbin |= 2;
				break;
			case 'T':
				tagbin |= 3;
				break;
		}
	}

	tagbin &= 0x3FFFFFFF;	// make sure the first nucleotide is encoded as 00 (irrespective of actual nucleotide) to be able to use as memory jump table. 
										// we will have one table for each starting nucleotide, so we can infer the actual nucleotide.

	for (i = 16; i < kbaitLen; i++) 
	{
		nextNbin <<= 2;
		switch(te->tag[i]) 
		{
			case 'A':
				break;
			case 'C':
				nextNbin |= 1;
				break;
			case 'G':
				nextNbin |= 2;
				break;
			case 'T':
				nextNbin |= 3;
				break;
		}
	}


	/* encode chr:pos */
	unsigned int pos = te->pos + virtchr[te->chr].offset;
	if (pos & 0xF0000000)
		printf("overflow: chr=%u pos=%u encoded in virt.chr %u\n",te->chr,te->pos,virtchr[te->chr].chr);
	if (verbose)
		printf("chr=%u pos=%u encoded in virt.chr %u virt.pos %u whichnt=%d\n",te->chr,te->pos,virtchr[te->chr].chr,pos,whichNnt);
	unsigned int chr = virtchr[te->chr].chr << 28;
	ptr = prm->hostdata;
	ptr += (whichNnt-prm->startingNnt) * k15NT_ENTRIES_COUNT + tagbin;
	*ptr = chr | pos;

	/* set the reverse flag if necessary */
	if (te->strand == '-')
	{
		unsigned long tmp = (whichNnt-prm->startingNnt) * k15NT_ENTRIES_COUNT + tagbin;  // position as if byte array.
		char setbit = (tmp & 0x07);	// get the position of the bit we need to set to 1.
		prm->stranddata[tmp >> 3] |= (1<<setbit);  // positions divided by 8, as we pack 8 bits in each char.
	}

	/* set the valid flag */
	{
		unsigned long tmp = (whichNnt-prm->startingNnt) * k15NT_ENTRIES_COUNT + tagbin;  // position as if byte array.
		char setbit = (tmp & 0x07);	// get the position of the bit we need to set to 1.
		prm->validdata[tmp >> 3] |= (1<<setbit);  // positions divided by 8, as we pack 8 bits in each char.
	}

	/* encode the current N nt */
	// shift to the appropriate position
	nextNbin <<= (whichNnt-prm->startingNnt) * kFollowNtLen * 2;
	prm->micdata[tagbin] |= nextNbin;
	if (verbose > 1)
		printf(
			"tagbin: %08x |= %08x   now = %08x ; whichNnt: %u ->  ptr=%016lx addr= %016lx, contains %08x\n\n",
			tagbin,
			nextNbin,
			prm->micdata[tagbin],
			whichNnt,
			(long unsigned int)ptr,
			(whichNnt-prm->startingNnt) * k15NT_ENTRIES_COUNT + tagbin,
			prm->hostdata[(whichNnt-prm->startingNnt) * k15NT_ENTRIES_COUNT + tagbin]);
	*lineok += 1;
	whichNntStat[whichNnt-prm->startingNnt] += 1;
}
//---------------------------------------------------------------

static inline int
isPAR(unsigned int chr, unsigned int pos, PARinfo_t *pi)
{
	unsigned int i;
	for (i = 0; i < pi->cnt; i++)
	{
		if (chr == pi->PAR[i].chr && pos >= pi->PAR[i].start && pos <= pi->PAR[i].end)
			return 1;
	}
	return 0;
}
//---------------------------------------------------------------

unsigned int
load(char *fn_base, PARAMS *prm)
{
	char fn1[512];
	char fn2[512];
	gzFile f1 = NULL;
	gzFile f2 = NULL;
	tagEntry_t te1;
	tagEntry_t te2;
	unsigned int linecnt1 = 0;
	unsigned int linecnt2 = 0;
	unsigned int lineok1 = 0;
	unsigned int lineok2 = 0;
	unsigned int PARcnt = 0;

	sprintf(fn1,"%s_all_1_mapped.txt.gz",fn_base);
	sprintf(fn2,"%s_all_2_mapped.txt.gz",fn_base);
	// /data6/hg38/hg38_18nt_all_1_mapped.txt.gz  /data6/hg38/hg38_18nt_all_2_mapped.txt.gz  /data6/hg38/hg38_18nt_all_n_mapped.txt.gz
	f1 = gzopen(fn1,"r");
	if (f1 == NULL)
	{
		fprintf(stderr,"Error opening %s\n",fn1);
		perror("gzopen");
		return(0);
	}
	f2 = gzopen(fn2,"r");
	if (f2 == NULL)
	{
		fprintf(stderr,"Error opening %s\n",fn2);
		perror("gzopen");
		return(0);
	}
	
	readTag(f1,&te1);
	readTag(f2,&te2);
	linecnt2 += 1;
	while (te1.tag[0] != 0 || te2.tag[0] != 0)
	{
		char prevtag[32];
		unsigned int whichNnt = 0;
		prevtag[0] = 0;
		// get tags from uniques
		if (te1.tag[0] != 0 && (te2.tag[0] == 0 || strncmp(te1.tag,te2.tag,16) <= 0))
		{
			linecnt1 += 1;
			strcpy(prevtag,te1.tag);
			encodeEntry(&te1,whichNnt,&lineok1,prm);
			whichNnt += 1;
			readTag(f1,&te1);
			while (te1.tag[0] != 0 && strncmp(te1.tag,prevtag,16) == 0)
			{
				linecnt1 += 1;
				// check if this was already seen
				if (strcmp(te1.tag,prevtag) != 0)
				{
					strcpy(prevtag,te1.tag);
					encodeEntry(&te1,whichNnt,&lineok1,prm);
					whichNnt += 1;
				}
				else
				{
					if (verbose)
						printf("skipping tag already seen %s\n",te1.tag);
				}
				readTag(f1,&te1);
			}
		}
		// add tags from duals unless when they are from the specified pseudo-autosomal region
		// (in principle, skip those from the Y chr)
		if (te2.tag[0] != 0 && (prevtag[0] == 0 || strncmp(prevtag,te2.tag,16) == 0))
		{
			// we want to store only complete pairs; we can have 4^nbNT times 2 (because they are pairs)
			tagEntry_t stack[kFollowNtEnum * 2];
			unsigned int cnt = 0;
			linecnt2 += 1;
			if (isPAR(te2.chr,te2.pos,prm->pi))
			{
				PARcnt += 1;
			}
			else
			{
				strcpy(prevtag,te2.tag);
				stack[cnt++] = te2;
				assert(cnt <= (kFollowNtEnum * 2));
			}
			readTag(f2,&te2);
			while (te2.tag[0] != 0 && strncmp(te2.tag,prevtag,16) == 0)
			{
				linecnt2 += 1;
				// check if this was already seen
				int diffTag = strcmp(te2.tag,prevtag);
				int notSeen = 1;
				// this is to filter out palindromes
				if (diffTag == 0)
				{
					unsigned int i = cnt;
					while (i > 0)
					{
						i -= 1;
						if (strcmp(te2.tag,stack[i].tag) != 0)
							break; // won't find it
						if (te2.chr == stack[i].chr && te2.pos == stack[i].pos)
						{
							// found it
							notSeen = 0;
							break;
						}
					}
				}
				if (notSeen)
				{
					if (isPAR(te2.chr,te2.pos,prm->pi))
					{
						PARcnt += 1;
					}
					else
					{
						if (diffTag != 0 && (cnt & 1) != 0)
						{
							// we start a new tag pair but have an odd number of tags in the stack
							//  => this means the top of stack is a member of the PAR region, and thus unique
							//  => we encode this unique tag and remove it from the stack
							cnt -= 1;
							encodeEntry(stack + cnt,whichNnt,&lineok2,prm);
							whichNnt += 1;
						}
						assert((diffTag == 0 && (cnt & 1) != 0) || ((diffTag != 0) && (cnt & 1) == 0));
						strcpy(prevtag,te2.tag);
						stack[cnt++] = te2;
						assert(cnt <= (kFollowNtEnum * 2));
					}
				}
				else
				{
					if (verbose)
						printf("skipping tag already seen %s\n",te2.tag);
				}
				readTag(f2,&te2);
			}
			if ((cnt & 1) != 0)
			{
				// the top of stack is necessarily a member of the PAR region, and thus unique
				cnt -= 1;
				encodeEntry(stack + cnt,whichNnt,&lineok2,prm);
				whichNnt += 1;
			}
			unsigned int i;
			for (i = 0; i < cnt; i += 2)
			{
				if (whichNnt + 1 == prm->startingNnt + kFollowNtCount)
				{
					// we will not be able to fit the pair, so we quit
					break;
				}
				encodeEntry(stack + i,whichNnt,&lineok2,prm);
				whichNnt += 1;
				encodeEntry(stack + i + 1,whichNnt,&lineok2,prm);
				whichNnt += 1;
			}
		}
	}

	printf("%s : loaded %u lines out of %u\n",fn1,lineok1,linecnt1);
	printf("%s : loaded %u lines out of %u (skipped %u PAR)\n",fn2,lineok2,linecnt2,PARcnt);
	
	gzclose(f1);
	gzclose(f2);
	return(lineok1 + lineok2);
}
//---------------------------------------------------------------

int
main(int argc, char **argv)
{
	unsigned int *following_nt_tbl = NULL;	// contain table of unique exact match 16mers.
	unsigned int *chrpos_tbl = NULL;		// contain table of chr,pos of each 16mer tag (note we have 8 successive such tables, corresponding to each of the 8 2nt followups).
	unsigned char *strand_tbl = NULL;		// contain table of strand for each of the Nmers.
	unsigned char *valid_tbl = NULL;	// tells if a Nmer is valid.
	char inputfile[512];
	char root[512];
	char config[512] = "/tmp/nodelete/nguex/data/hg19.cfg";
	int err;
	int c;
	char NT = '\0';
	int startingNnt = 0;
	PARinfo_t Pi;
	
	/* --------- process arguments */

	opterr = 0;
	strcpy(root,"./");
	inputfile[0] = 0;
	memset(&Pi,0,sizeof(Pi));
	
	while ((c = getopt (argc, argv, "c:i:o:P:n:s:v:")) != -1)
	switch (c)
	{
	  case 'c':
			strcpy(config,optarg);
			break;
	  case 'i':
			strcpy(inputfile,optarg);
			break;

	  case 'n':
			sscanf(optarg,"%c",&NT);
			break;

	  case 'o':
			strcpy(root,optarg);
			break;

	  case 'P':
			if (sscanf(optarg,"%u:%u-%u",&Pi.PAR[Pi.cnt].chr,&Pi.PAR[Pi.cnt].start,&Pi.PAR[Pi.cnt].end) != 3)
			{
				fprintf(stderr,"Failed to grok PAR region %s\n",optarg);
				exit(1);
			}
			Pi.cnt += 1;
			break;

	  case 's':
			sscanf(optarg,"%d",&startingNnt);
			break;

	  case 'v':
			sscanf(optarg,"%d",&verbose);
			break;
	}

	if ( (inputfile[0] == 0) || ((NT != 'A') && (NT != 'C') && (NT != 'G') && (NT != 'T'))  )
	{
		printf("usage:\n\n");
        printf("encode_UEM -c configfile -i inputfile -o output_root_name -n [A|C|G|T|] -s start [-v level]\n\n");
		printf("           -c configfile         : chromosome config file [%s]\n",config);
		printf("           -i inputfile          : gzipped sorted file with unique %dnt for the genome of interest (Chris data)\n",kbaitLen);
		printf("           -o output_root_name   : used for the four output files\n");
		printf("                                   [A|C|G|T].[%dnt|chr|strand|valid].bin will be automatically created.\n",kFollowNtLen);
		printf("           -n nucleotide         : encode only sequences starting with nucleotide\n");
		printf("           -P chr:start-end      : pseudoautosomal region (PAR) to handle as unique instead of dual\n");
		printf("                                   the specified region will be excluded from the table (so in principle choose Y region here)\n");
		printf("                                   (repeat several times to enter several regions)\n");
		printf("           -s start              : encode %d unique %dmers, starting from this value\n",kFollowNtCount,kbaitLen);
		printf("                                   default is zero. A supplementary table using the leftover that could not fit\n");
		printf("                                   can be produced with -s %d  then -s %d  and so on...\n",kFollowNtCount,kFollowNtCount * 2);
		printf("           -v level              : verbose level\n");
		printf("\nExamples:\n\n");
		printf("encode_UEM -i /scratch/fhgfs/chris/hg19/ht_all_mapped.txt.gz   -o /scratch/local/nguex/ht_all_mapped-  -n A\n");
		printf("encode_UEM -i /scratch/fhgfs/chris/hg19/ht_all_1_mapped.txt.gz -o /scratch/local/nguex/   -n A\n\n");
		printf("encode_UEM -i /scratch/fhgfs/chris/hg19/ht_all_1_mapped.txt.gz -o /scratch/local/nguex/start4. -s 4   -n A\n\n");
		printf(GTL_VERSION_FULL_STRING "\n");
		printf("(c) Nicolas Guex & Christian Iseli 2014-2016\n");
		return(1);
	}

	/* ---- do it. */
	
	printf("sizeof(int):%ld bits\n",8*sizeof(unsigned int));
	printf("allocating chrpos_tbl %lu bytes\n",kFollowNtCount*k15NT_ENTRIES_COUNT*sizeof(unsigned int));
	printf("allocating following_nt_tbl %lu bytes\n",k15NT_ENTRIES_COUNT*sizeof(unsigned int));
	printf("allocating strand_tbl %lu bytes\n",kSTRAND_DATA_SIZE*sizeof(char));
	printf("allocating valid_tbl %lu bytes\n",kVALID_DATA_SIZE*sizeof(char));
	err  = posix_memalign((void **)&following_nt_tbl, 64, k15NT_ENTRIES_COUNT*sizeof(unsigned int));
	err += posix_memalign((void **)&chrpos_tbl, 64, kFollowNtCount*k15NT_ENTRIES_COUNT*sizeof(unsigned int)); 
	strand_tbl = calloc(kSTRAND_DATA_SIZE, sizeof(char)); 
	valid_tbl = calloc(kVALID_DATA_SIZE, sizeof(char)); 
	if (err == 0 && (strand_tbl != NULL) && (valid_tbl != NULL))
	{
		PARAMS prm;
		// posix_memalign manpage does not mention that the memory is zeroed
		memset(following_nt_tbl,0,k15NT_ENTRIES_COUNT*sizeof(unsigned int));
		memset(chrpos_tbl,0,kFollowNtCount*k15NT_ENTRIES_COUNT*sizeof(unsigned int));
		prm.pi = &Pi;
		initializeVirtualChromosomes(config,virtchr);
		prm.micdata = following_nt_tbl;
		prm.hostdata = chrpos_tbl;
		prm.stranddata = strand_tbl;
		prm.validdata = valid_tbl;
		prm.startingNnt = startingNnt;
		prm.NT = NT;
		if (load(inputfile,&prm) > 0)
		{
			char fn[512];
			unsigned int i;
			sprintf(fn,"%s%c.%dnt.bin",root,NT,kFollowNtLen);
			write2NTtable (following_nt_tbl,k15NT_ENTRIES_COUNT,fn);
			sprintf(fn,"%s%c.chr.bin",root,NT);
			writeCHRPOStable(chrpos_tbl,k15NT_ENTRIES_COUNT,fn);
			sprintf(fn,"%s%c.strand.bin",root,NT);
			writeStrandOrValidTable(strand_tbl,kSTRAND_DATA_SIZE,fn);
			sprintf(fn,"%s%c.valid.bin",root,NT);
			writeStrandOrValidTable(valid_tbl,kVALID_DATA_SIZE,fn);
			printf("Counts of used whichNnt slots:\n");
			for (i = 0; i < kFollowNtCount; i++)
				printf("%u %u\n",i,whichNntStat[i]);
		}
		free(following_nt_tbl);
		free(chrpos_tbl);
		free(strand_tbl);
		free(valid_tbl);
		return(0);
	} 
	else
	{
		printf("cannot allocate data\n\n");
		return(1);
	}
		
} /* main */
//---------------------------------------------------------------

/* ------------------------------------------------------------------------------------ */
/* vim: tabstop=2 shiftwidth=2
 */
