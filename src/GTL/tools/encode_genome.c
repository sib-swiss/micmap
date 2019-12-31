/*
 * ------------------------------------------------------------------------------------------------------------------------
 *
 *                                * micmap *
 *      Mapping of short reads in fastq format onto a reference
 *
 *
 *  Copyright (C) SIB  - Swiss Institute of Bioinformatics,  2014-2019 Nicolas Guex, Thierry Schuepbach and Christian Iseli
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

  This software will encode the chromosomes in 16 virtual chromosomes.

  It will produce one table:
  

  -  a 4Gb table which contains the nucleotides.

	  /scratch/local/hg19.bin


	(c) N.Guex and C.Iseli 2014

*/
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
#include <ctype.h>
#include <sys/time.h>
#include "virt_chr.h"
#include "gtlVersion.h"

/*
        cc -m64 -O3 -o encode_genome encode_genome.c virt_chr.c  -Wall  
*/
//---------------------------------------------------------------

static VIRTUALCHR virtchr[kMAXCHRcnt];
static int verbose = 0;
unsigned long datasize = 0x100000000;


//---------------------------------------------------------------
void writeGenome(unsigned char *data,size_t datasize,char *fn)
{
	FILE *f = NULL;

		f = fopen(fn,"wb");
		if (f == NULL)
		{
			printf("Error writing %s\n",fn);
			return;
		}
		printf("writing genome in %s\n",fn);
		fwrite(data,sizeof(unsigned char),datasize,f);
		fclose(f);

	
} /* writeGenome */
//---------------------------------------------------------------

unsigned int load(unsigned char *data,VIRTUALCHR *virtchr, char *fn)
{
	FILE *f = NULL;
	unsigned int pos;
	char buffer[256];
	unsigned int linecnt = 0;
	
	f = fopen(fn,"r");
	if (f == NULL)
	{
		printf("Error opening %s\n",fn);
		return(0);
	}

	pos = (virtchr->chr << 28) + virtchr->offset +1;  // to be 1 based.
	printf("Loading virtual chromosome %d from pos %u  addr=%x\n",virtchr->chr,virtchr->offset,pos);
	if (fgets(&buffer[0],256,f) != NULL) { // skip header
	  do
	  {
		  int i;
		  
		  if (fgets(&buffer[0],256,f) != NULL && verbose > 1) {
			  printf("%s",buffer);
		  }
		  if (feof(f))
			  break;
			  
		  i = 0; while(buffer[i] != '\n') { data[pos++] = toupper(buffer[i]); i++; }
		  
		  linecnt++;
	  //	if (linecnt == 1024)
	  //		break;
			  
	  } while(1);
	}

	fclose(f);
	return(linecnt);
	
} /* load */
//---------------------------------------------------------------

int main (int argc, char **argv)
{
	unsigned char *genome_tbl = NULL;		// contain table of nucleotides.
	int err;
	unsigned int chr;
	
	/* --------- process arguments */
	if (argc != 2)
	{
		fprintf(stderr,"Usage: %s <config file>\n%s\n",argv[0], GTL_VERSION_FULL_STRING);
		exit(1);
	}


	printf("sizeof(int):%ld bits\n",8*sizeof(unsigned int));
	printf("allocating genome_tbl %lu bytes\n",1L*datasize*sizeof(char));
	err  = posix_memalign((void **)&genome_tbl, 64, datasize*sizeof(char));
	memset(genome_tbl,'N',datasize);
	if (err == 0)
	{
		initializeVirtualChromosomes(argv[1],virtchr);

#if 0
		for (chr = 1; chr <= 22; chr++)
		{
			char fn[256];
			sprintf(fn,"/db/genome/hg19/chrom%d.seq",chr);
			if (load(genome_tbl,virtchr[chr],fn) == 0)
				printf("Error loading chr %d",chr);
		}
		if (load(genome_tbl,virtchr[23],"/db/genome/hg19/chromX.seq") == 0)
			printf("Error loading chrMt");
		if (load(genome_tbl,virtchr[24],"/db/genome/hg19/chromY.seq") == 0)
			printf("Error loading chrMt");
		if (load(genome_tbl,virtchr[25],"/db/genome/hg19/chromMt.seq") == 0)
			printf("Error loading chrMt");
#endif
		for (chr = 1; chr <= chrMax; chr++)
		{
			char fn[256];
			sprintf(fn,"%s.seq",virtchr[chr].SAMname);
			if (load(genome_tbl,virtchr + chr,fn) == 0)
				printf("Error loading chr %d\n",chr);
		}
			

unsigned int chr = 2;
unsigned int pos = 70014926;
unsigned int genomepos = (virtchr[chr].chr << 28) + virtchr[chr].offset + pos;
printf("chr %d pos %d : %.102s\n",chr,pos,&genome_tbl[genomepos]);

 chr = 2;
pos = 70931321;
genomepos = (virtchr[chr].chr << 28) + virtchr[chr].offset + pos;
printf("chr %d pos %d : %.102s\n",chr,pos,&genome_tbl[genomepos]);

chr = 2;
pos = 70014846;
genomepos = (virtchr[chr].chr << 28) + virtchr[chr].offset + pos;
printf("chr %d pos %d : %.102s\n",chr,pos,&genome_tbl[genomepos]);


chr = 2;
pos = 70014692;
genomepos = (virtchr[chr].chr << 28) + virtchr[chr].offset + pos;
printf("chr %d pos %d : %.102s\n",chr,pos,&genome_tbl[genomepos]);


chr = 2;
pos = 70014926;
genomepos = (virtchr[chr].chr << 28) + virtchr[chr].offset + pos;
printf("chr %d pos %d : %.102s\n",chr,pos,&genome_tbl[genomepos]);


chr = 2;
pos = 70014585;
genomepos = (virtchr[chr].chr << 28) + virtchr[chr].offset + pos;
printf("chr %d pos %d : %.102s\n",chr,pos,&genome_tbl[genomepos]);


chr = 2;
pos = 70014692;
genomepos = (virtchr[chr].chr << 28) + virtchr[chr].offset + pos;
printf("chr %d pos %d : %.102s\n",chr,pos,&genome_tbl[genomepos]);


chr = 2;
pos = 70014971;
genomepos = (virtchr[chr].chr << 28) + virtchr[chr].offset + pos;
printf("chr %d pos %d : %.102s\n",chr,pos,&genome_tbl[genomepos]);

chr = 16;
pos = 70014971;
genomepos = (virtchr[chr].chr << 28) + virtchr[chr].offset + pos;
printf("chr %d pos %d : %.102s\n",chr,pos,&genome_tbl[genomepos]);





		char *of = strdup(argv[1]);
		unsigned int oflen = strlen(of);
		memcpy(of+oflen-3,"bin",3);
		writeGenome (genome_tbl,datasize,of);
		free(genome_tbl);
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
