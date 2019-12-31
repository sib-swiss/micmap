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

#include "GTL.h"
#include "gtlVersion.h"

//=================================================================================================

static int
list(const char *fn)
{
	int fd = open(fn, O_RDONLY);
	TLBHEADER th;
	int err = 1;
	TLBDATA td;
	unsigned int blockNo = 0;
	allocBlock(&td, true);
	printf("block\tchr\tminPos\tmaxPos\tcnt\tlen_O\tlen_P\tpair_O\tpair_P\tmm_O\tmm_P\tcigar_O\tcigar_P\tumap_O\tumap_P\thdr_O\thdr_P\tqual_O\tqual_P\tpair2_O\tpair2_P\tmm2_O\tmm2_P\tcigar2_O\tcigar2_P\n");
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
		if (readDecompressBlock(fd,&th,&td,1) != 0)
		{
			err = 1;
			goto bail;
		}
		if ((td.header.chr >> 16) != 0)
			printf("%u\t%u|%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n",
				blockNo++,
				td.header.chr >> 16,td.header.chr & 0xffff, td.header.minPos, td.header.maxPos, td.cnt,
				td.readLen.origLen,td.readLen.packedLen,
				td.pairPos[0].origLen,td.pairPos[0].packedLen,
				td.mismatch[0].origLen,td.mismatch[0].packedLen,
				td.cigarBlock[0].origLen,td.cigarBlock[0].packedLen,
				td.unmapped.origLen,td.unmapped.packedLen,
				td.readHdr.origLen,td.readHdr.packedLen,
				td.readQual.origLen,td.readQual.packedLen,
				td.pairPos[1].origLen,td.pairPos[1].packedLen,
				td.mismatch[1].origLen,td.mismatch[1].packedLen,
				td.cigarBlock[1].origLen,td.cigarBlock[1].packedLen);
		else
			printf("%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n",
				blockNo++,
				td.header.chr, td.header.minPos, td.header.maxPos, td.cnt,
				td.readLen.origLen,td.readLen.packedLen,
				td.pairPos[0].origLen,td.pairPos[0].packedLen,
				td.mismatch[0].origLen,td.mismatch[0].packedLen,
				td.cigarBlock[0].origLen,td.cigarBlock[0].packedLen,
				td.unmapped.origLen,td.unmapped.packedLen,
				td.readHdr.origLen,td.readHdr.packedLen,
				td.readQual.origLen,td.readQual.packedLen,
				td.pairPos[1].origLen,td.pairPos[1].packedLen,
				td.mismatch[1].origLen,td.mismatch[1].packedLen,
				td.cigarBlock[1].origLen,td.cigarBlock[1].packedLen);
	}

	err = 0;

bail:;

	// clean up
	close(fd);

	return(err);

} // sort
//---------------------------------------------------------------

int
main (int argc, char **argv)
{
	int c;
	char rsltfile[512];
	int err;

	/* --------- process arguments ----------------------*/

	rsltfile[0] = 0;
	opterr = 0;
	while ((c = getopt (argc, argv, "r:")) != -1)
	switch (c)
	{
	  case 'r':
			strcpy(rsltfile,optarg);
			break;
	}

	if (rsltfile[0] == 0)
	{
		printf("usage:\n\n");
        printf("list -r resultfile\n\n");
		printf("           -r resultfile         : name of the result file (produced by match)\n");
		printf("\n");
		printf(GTL_VERSION_FULL_STRING "\n");
		printf("(c) Nicolas Guex & Christian Iseli & Thierry Schuepbach 2015-2018\n");
		return(1);
	}

	err = list(rsltfile);

	return(err);

} /* main */
//---------------------------------------------------------------

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
