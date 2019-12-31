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
#define _GNU_SOURCE     /* Expose declaration of tdestroy() */
#include "config.h"
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
#include <semaphore.h>
#include <assert.h>
#include <errno.h>

// basic compression
#include <zlib.h>
#include "Genome.h"
#include "virt_chr.h"
#include "GTL.h"
#include "gtlVersion.h"
#include "GTLdispatch.h"

Genome_t Genome = {
	.FileName = NULL,
	.table = NULL,
	.configBuffer = NULL,
	.tableSize = 0UL,
	.configBufferSize = 0UL
};

static void __attribute__((noreturn)) usage()
{
		printf("usage:\n\n");
		printf("GTLindex -g genomefile  gtl file [-o output path]\n\n");
		printf(GTL_VERSION_FULL_STRING "\n");
		printf("(c) Thierry Schuepbach & Cie 2014-2017\n");
		exit(1);
}
//---------------------------------------------------------------

int main(int argc, char *argv[])
{
	const char * * InputGTLFiles = NULL;
	GTLList_t * restrict GTLInputFileList = NULL;
	const char * restrict outputDirectory = NULL;
	unsigned int nInputGTLFiles = 0U;
	int res = 1, c;
	int nFetchThreads = (int) sysconf(_SC_NPROCESSORS_ONLN);
	_Bool DoNotUseExisting = false;

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Process arguments

	opterr = 0;
	while ((c = getopt (argc, argv, "g:o:f")) != -1)
	switch (c)
	{
		case 'g':
			Genome.FileName = optarg;
			break;
		case 'o':
			outputDirectory = optarg;
			break;
		case 'f':
			DoNotUseExisting = true;
			break;
		case '?':
		case 'h':
		default:
			usage();
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Get input file names
	if (optind < argc) {
		nInputGTLFiles = argc - optind;
		fprintf(stderr,"%i files provided\n", argc - optind);
		InputGTLFiles = alloca(nInputGTLFiles*sizeof(char*));
		for (unsigned int iFile=0; iFile<nInputGTLFiles; iFile++) {
			InputGTLFiles[iFile] = argv[iFile+optind];
			fprintf(stderr,"Input file %2i is %s\n", iFile, InputGTLFiles[iFile]);
		}
	}
	else {
		usage();
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Load genome
	fprintf(stderr, "Loading genome file %s\n", Genome.FileName);
	if (LoadGenome(&Genome) != NULL) {
		fprintf(stderr, "failed to load genome file %s\n", Genome.FileName);
		goto bail;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Prepare dispatcher, indexing file if necessary and sorting
	GTLInputFileList = createGTLList(InputGTLFiles, nInputGTLFiles);
	if (GTLInputFileList == NULL) {
		fputs("Error while preparing GTL input files\n", stderr);
		goto bail;
	}

	if (indexGTLList(GTLInputFileList, Range, DoNotUseExisting, 0)) {
		fputs("Error indexing GTL list\n", stderr);
		goto bail;
	}
	fprintf(stderr, "Found %lu blocks\n", GTLInputFileList->nBlocks);
	if ( GTLInputFileList->PositionType != Range ) {
		fputs("Getting real Block position extrema\n", stderr);
		getRealBlockMinMaxPosition( GTLInputFileList, nFetchThreads, 1 /*DoSort*/);
		fprintf(stderr, "Writing those indices to file...\n");
		if(writeGTLFileIndex(GTLInputFileList, outputDirectory) != 0) {
			fputs("Error in writing indices...\n", stderr);
		}
	}


	res = 0;
bail: ;

	return res;
}
//---------------------------------------------------------------

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
