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
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <stdatomic.h>
#include <pthread.h>
#include <semaphore.h>
#include <assert.h>

#include "config.h"
#include "Genome.h"
#include "virt_chr.h"
#include "GTL.h"
#include "htslib/vcf.h"
#include "htslib/vcfutils.h"

#include "GTLdispatch.h"
#include "pfProfile.h"
#include "pfCompute.h"
#include "AnalyzeTags.h"

#define MINIMUM_APPEARANCE 5
#define MIN_BORDER_DISTANCE 10

struct ExtraData {
	char * ptr;
	unsigned short int * Accounting;
	unsigned short int * PHred;
	unsigned int n;
};

struct Ordering {
	struct ExtraData data;
	unsigned int count;
};

enum StatePriority {
  PRIORITY_MATCH     = 0,
  PRIORITY_INSERTION = 1,
  PRIORITY_DELETION  = 2,
  PRIORITY_EXTRA     = 3
};

#if TAG_ANALYSIS_MODE == NIRVANA
static const char *const types[] = {
	"SNP      ",
	"DELETION ",
	"INSERTION"
};

static const char *const Errors[] = {
	"None                   ",
	"Heterozygote/homozygote",
	"Allele0 not found      ",
	"Allele1 not found      ",
	"Both alleles not found "
};

extern FILE * outfile;
#endif

static const char DNAAlphabet[] = "ACGTN";

extern _Bool OutputVerbose;

static int MyOrdering(const void* A, const void* B)
{
	const struct Ordering * const restrict ta = A, * const restrict tb = B;
	if (ta->count < tb->count)
		return 1;
	else if (ta->count > tb->count) 
		return -1;
	else
	{
		if (ta->data.n > tb->data.n)
			return -1;
		else if (ta->data.n < tb->data.n)
			return 1;
		else 
			return 0;
	}
}
//---------------------------------------------------------------

static int MAFOrdering(const void* A, const void* B)
{
	const struct Ordering * const restrict ta = A, * const restrict tb = B;
	if (ta->count < tb->count)
		return 1;
	else if (ta->count > tb->count) 
		return -1;
	else
		return 0;
}
//---------------------------------------------------------------

static int MyCompare(const void* A, const void* B)
{
	const struct ExtraData *const restrict ta = A, * const restrict tb = B;
	if (ta->n < tb->n)
		return 1;
	else if (ta->n > tb->n) 
		return -1;
	else
		return 0;
}
//---------------------------------------------------------------

void* AnalyzeTags3(threadpool_t * restrict const thpool) 
#include "./AnalyzeTags3.inc"

#define ANALYZE_FULL_OUTPUT
void* AnalyzeTags3_FullOutput(threadpool_t * restrict const thpool)
#include "./AnalyzeTags3.inc"
#undef ANALYZE_FULL_OUTPUT

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
