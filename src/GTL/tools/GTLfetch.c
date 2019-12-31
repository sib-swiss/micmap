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

//=================================================================================================

//static struct timeval ts;
//static struct timeval te;
static int verbose = 0;
static int debug = 0;
static VIRTUALCHR virtchr[kMAXCHRcnt];

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

int
main (int argc, char **argv)
{
	int c;
	char genomefile[512];
	int err;
	int chr = -1;
	unsigned int fromPos = 0;
	unsigned int toPos = INT_MAX;  // this one is to large or negative ?  and proidices no output... 2147483647;

	/* --------- process arguments ----------------------*/

	strcpy(genomefile,"/tmp/nodelete/nguex/data/hg19.bin");
	opterr = 0;
	while ((c = getopt (argc, argv, "g:C:P:v:d")) != -1)
	switch (c)
	{
		case 'C':
			sscanf(optarg,"%d",&chr);
			break;

		case 'P':
			sscanf(optarg,"%d..%d",&fromPos,&toPos);
			break;

		case 'g':
			strcpy(genomefile,optarg);
			break;

		case 'd':
			debug = 1;
			break;

		case 'v':
			sscanf(optarg,"%d",&verbose);
			break;
	}
	if (chr == -1 || toPos < fromPos)
	{
		printf("usage:\n\n");
		printf("GTLfetch [-g genomefile] -C chromosome [-P from..to] [-v level]\n\n");
		printf("           -g genomefile         : full path of the genome encoded by encode_genome\n");
		printf("                                   default: '%s' \n",genomefile);
		printf("           -C chromosome         : select chromosome\n");
		printf("           -P fromPos..ToPos     : to restrict to a specific chromosme region, default means all\n");
		printf("                                   note: for ADNIview, GTL file needs to be sorted !  only one kind of tag per call !!\n");
		printf("           -d                    : debug mode\n");
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

	printf("%.*s\n",toPos - fromPos + 1, genome_tbl + (virtchr[chr].chr << 28) + virtchr[chr].offset + fromPos);

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
