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
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <sys/time.h>
#include <stdatomic.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>

#include <pthread.h>
#include <semaphore.h>
#include <assert.h>

// basic compression
#include <zlib.h>
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

struct ExtraData {
	char * ptr;
	unsigned int n;
};

struct Ordering {
	struct ExtraData data;
	unsigned int count;
};


enum StatePriority {
  PRIORITY_MATCH     = 1,
  PRIORITY_INSERTION = 2,
  PRIORITY_DELETION  = 0,
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
FILE * outfile = NULL;
#endif


static const char DNAAlphabet[] = "ACGTN";

_Bool OutputVerbose = true;

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

_Bool PerformTwoPass = false;

void* AnalyzeTags(threadpool_t * restrict const thpool) {
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Required Memory
	char * restrict PFResults = NULL;
	char * restrict TransposePFResults = NULL;
	unsigned char * restrict GenomeBuffer = NULL;
	char * restrict Buffer = NULL;
	short int (*CellList)[2] = NULL;
	int * score = NULL;
	union lScores * restrict matrix = NULL;
	int * restrict WORK = NULL;
	FILE * out;

	size_t PFResults_size = 0UL;
	size_t TransposePFResults_size = 0UL;
	size_t score_size = 0UL;
	size_t GenomeBuffer_size = 0UL;
	size_t Buffer_size = 0UL;
	size_t matrix_size = 0UL;
	size_t WORK_size = 0UL;
	
	struct Profile profile;
	char ProfileName[128];
	char Allele0[MAX_ALLELE_SIZE], Allele1[MAX_ALLELE_SIZE];
	unsigned char Indices[kMaxReadLen+1];
	
	unsigned int nRecords=0U, nRecordsOK=0U;

#if TAG_ANALYSIS_MODE == NIRVANA
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
#endif

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Set some constant values
	const size_t tid = (size_t) pthread_self();

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Mark thread as alive (initialized)
	pthread_mutex_lock(&thpool->thcount_lock);
	thpool->num_threads_alive++;
	pthread_mutex_unlock(&thpool->thcount_lock);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Handle jobqueue
	while(thpool->threads_keepalive){
		bsem_wait(thpool->jobqueue_p.has_items);

		if (thpool->threads_keepalive){
			/* Signal I am working */
			pthread_mutex_lock(&thpool->thcount_lock);
			thpool->num_threads_working++;
			pthread_mutex_unlock(&thpool->thcount_lock);

			/* Read job from queue and execute it */
			pthread_mutex_lock(&thpool->jobqueue_p.rwmutex);
			IntervalJob_t * const restrict task = (IntervalJob_t*) jobqueue_pull(&thpool->jobqueue_p);
			pthread_mutex_unlock(&thpool->jobqueue_p.rwmutex);

			/*************************************************************************/
			/*                             START THE TASK                            */
			/*************************************************************************/
			if (task) {
				GTLRawData_t ** ND = (GTLRawData_t **) task->data;
				/*const*/ unsigned int nTags = task->count;
				assert(nTags >= 4);
// 				printf("ANA: GOT %u tags for %u-%u\n", nTags,  task->GenomeStart, task->GenomeEnd);
				unsigned int minPos = -1;
				unsigned int maxPos = 0U;
				for (int iTag=0; iTag<nTags; iTag++) {
					unsigned int uitmp = ND[iTag]->Location;
					if (uitmp < minPos) minPos = uitmp;
					uitmp += ND[iTag]->AlignmentRange[1];
					if (uitmp > maxPos) maxPos = uitmp;
				}
				
				/* Extra margin of both sides of 3 DNA bases */
				maxPos += 2;
				if (task->GenomeStart < minPos) task->GenomeStart = minPos;
				if (task->GenomeEnd > maxPos) task->GenomeEnd = maxPos;
				if (minPos >= 2) minPos -= 2;

#ifdef DEBUG_INFO
				{
					/* Print alignments in global region */
					fprintf(stderr, "\n\n%12lu   ", 0UL);
					unsigned int l;
					for (l=minPos; l<(task->GenomeStart-1); l++) fputc(' ', stderr);
					fputc(']', stderr);

					while (l++ < task->GenomeEnd) fputc('=', stderr);
					fputs("[\n", stderr);

					snprintf(ProfileName, 128, "/tmp/Tags_c%i_%u-%u.fasta", options.chr, minPos, maxPos);
					FILE * tags = fopen(ProfileName, "w");
					for (int iTag=0; iTag<nTags; iTag++) {
						GTLRawData_t * restrict GD = ND[iTag];

						unsigned int TagGenLen   = GD->AlignmentRange[1] - GD->AlignmentRange[0];
						unsigned int TagGenBegin = GD->Location;
						unsigned int TagGenEnd   = TagGenBegin + TagGenLen;

						fprintf(stderr, "%12lu   ", GD->Ordinal);
						unsigned int space;
						for (space= minPos; space<TagGenBegin; space++) fputc('.',stderr);
						fprintf(stderr, "%s", GD->Tag);
						space += GD->TagLength;
						while (space++ <= maxPos) fputc('.',stderr);
						fprintf(stderr, "\t%u-%u\t", GD->AlignmentRange[0], GD->AlignmentRange[1]);
						if (tags) fprintf(tags, ">%lu %u-%u ", GD->Ordinal, GD->AlignmentRange[0], GD->AlignmentRange[1]);
						if (GD->CigarLength > 0) {
							const unsigned char * restrict cigar = GD->Cigar;
							while (cigar[0] != kCIGARTERMINATOR) {
								fprintf(stderr, "%u%c ", (unsigned int) cigar[0], (int) cigar[1]);
								if (tags) fprintf(tags, "%u%c ", (unsigned int) cigar[0], (int) cigar[1]);
								cigar += 2;
							}
						}
						if (tags) {
							fprintf(tags, "\n%.*s\n", (int) GD->TagLength, GD->Tag);
						}
						fputc('\n', stderr);
					}
					if (tags) fclose(tags);
				}
#endif
				unsigned int TotalLength;
				{
					const unsigned int Length = maxPos - minPos + 1;
					unsigned int Deletions[3*Length];
					unsigned int * const Insertions = Deletions + Length;
					unsigned int * const restrict Memory = Deletions + 2*Length;
					memset(Deletions, 0, 3*Length*sizeof(unsigned int));

					for (int iTag=0; iTag<nTags; iTag++) {
						GTLRawData_t * restrict GD = ND[iTag];
						unsigned int * restrict ptr = Deletions;
						assert(GD->Location >= minPos);
						size_t pos = GD->Location - minPos;
						if (GD->CigarLength > 0) {
							const unsigned char * restrict cigar = GD->Cigar;
							unsigned int TagPos = 0U;
							assert(GD->TagLength > 10);
#define MIN_BORDER_DISTANCE 10
							const unsigned int TagPosMax = GD->TagLength - MIN_BORDER_DISTANCE;
							while (cigar[0] != kCIGARTERMINATOR) {
								const unsigned int N = cigar[0];
								const unsigned char type = cigar[1];
								assert(pos <= Length);
								if (type == 'D') {
									if (N > Deletions[pos] && (TagPos >= MIN_BORDER_DISTANCE) && ((TagPos+N) <= TagPosMax)) Deletions[pos] = N;
									pos += N;
								}
								else if (type == 'I') {
									// It is wrong not to take these into account as that would force insertion that we throw away next.
									if (N > Insertions[pos] /*&& (TagPos >= MIN_BORDER_DISTANCE) && ((TagPos+N) <= TagPosMax)*/) Insertions[pos] = N;
									if (N > Memory[pos]) Memory[pos] = N;
									TagPos += N;
								}
								else {
									pos += N;
									TagPos += N;
								}
								cigar += 2;
							}
						}
					}

#ifdef DEBUG_INFO
					fprintf(stderr, "\n\n%12lu   ", 0UL);
					for (unsigned int l=0; l<Length; l++) { 
						int c = ' ';
						if (Deletions[l] > 0 && Insertions[l] > 0) c = (int) '*';
						else if (Deletions[l] > 0) c = (int) '-';
						else if (Insertions[l] > 0) { for (int m=1; m<Insertions[l]; m++) fputc('+', stderr); c = (int) '+'; }
						fputc( c, stderr);
					}
					fputc('\n', stderr);
#endif
					
					TotalLength = Length;
					for (unsigned int k=0; k<Length; k++) TotalLength += Memory[k];

					if (TotalLength > GenomeBuffer_size) {
						GenomeBuffer_size = TotalLength;
						GenomeBuffer = realloc(GenomeBuffer, GenomeBuffer_size*sizeof(char));
						if (GenomeBuffer == NULL) {
							fprintf(stderr, "Unable to allocate memory for the Genome string in AnalyzeTags (%lu)\n", GenomeBuffer_size);
							exit(1);
						}
// 						memset(GenomeBuffer,'?', TotalLength);
					}
					{
						const size_t index = (Genome.virtchr[options.chr].chr << 28) + Genome.virtchr[options.chr].offset + minPos;
						const unsigned char * restrict cptr = Genome.table + index;
						unsigned char * restrict cptro = GenomeBuffer;

						unsigned int k=0U;
						while (k<Length) {
							for (unsigned int j=0; j<Insertions[k]; j++) *cptro++ = ' ';
							*cptro++ = *cptr++;
							k++;
						}
// 						while (++k < TotalLength) { 
// 							assert((uintptr_t) cptro <= (uintptr_t) &GenomeBuffer[TotalLength-1]);
// 							*cptro++ = *cptr++;
// 						}
						while ((uintptr_t) cptro < (uintptr_t) &GenomeBuffer[TotalLength]) *cptro++ = *cptr++;
					}
				}

				/* All in memory */
				{
					/* Check if enough memory, otherwise increase */
					{
						const size_t needed_size = nTags*TotalLength;
						if (needed_size > Buffer_size) {
							Buffer_size = needed_size;
							Buffer = realloc(Buffer, needed_size*sizeof(char));
							if ( Buffer == NULL) {
								fprintf(stderr, "Unable to allocate Buffer memory (%lu)\n", needed_size);
								exit(1);
							}
						}
						memset(Buffer, '.', needed_size*sizeof(char));
					}

					/* Build the buffer space containing all tags in the region */
					_Bool HasGap = false;
					for (int iTag=0; iTag<nTags; iTag++) {
						GTLRawData_t * restrict GD = ND[iTag];
						unsigned int space = 0U;
						{
							unsigned int count = minPos;
							while (count<GD->Location) {
								if (GenomeBuffer[space] == ' ') {
									Buffer[nTags*space+iTag] = ' ';
								}
								else {
									count++;
								}
								space++;
							}
						}
						const unsigned char * restrict cptr = GD->Tag;
						if (GD->CigarLength > 0) {
							const unsigned char * restrict cigar = GD->Cigar;
							while (cigar[0] != kCIGARTERMINATOR) {
								const unsigned int N = cigar[0];
								const unsigned char type = cigar[1];
								unsigned int count = 0U;
								switch (type) {
									case 'S':
									case 'M':
										{
											while (count < N) {
												if (GenomeBuffer[space] == ' ') {
													Buffer[nTags*space+iTag] = ' ';
												}
												else {
													Buffer[nTags*space+iTag] = (*cptr++);
													count++;
												}
												space++;
											}
										}
										break;
									case 'D':
										{
											while (count < N) {
												if (GenomeBuffer[space] == ' ') {
													Buffer[nTags*space+iTag] = ' ';
												}
												else {
													Buffer[nTags*space+iTag] = IF_DEBUG_INFO((N>=longuestAllowedDeletion) ? ' ':) '-';
													count++;
												}
												space++;
											}
											HasGap = true;
										}
										break;
									case 'I':
										{
											while (count<N) {
												Buffer[nTags*space+iTag] = (*cptr++);
												space++; 
												count++;
											}
											HasGap = true;
										}
										break;
								}
								cigar += 2;
							}
						}
						else {
							unsigned int count = 0U;
							while (count <= GD->AlignmentRange[1]) {
								if (GenomeBuffer[space] == ' ') 
									Buffer[nTags*space+iTag] = ' ';
								else {
									Buffer[nTags*space+iTag] = (*cptr++);
									count++;
								}
								space++;
							}
						}
						while (space < TotalLength) {
							if (GenomeBuffer[space] == ' ') Buffer[nTags*space+iTag] = ' ';
							space++;
						}

						if (space > TotalLength) {
							fprintf(stderr, "Space is out of bound %u > %u\n%.*s\n",
							        space, TotalLength, (int) TotalLength, GenomeBuffer);
							for (unsigned int k=0; k<TotalLength; k++) fputc((int) Buffer[k*nTags+iTag], stderr);
							fprintf(stderr, "\nTAG: %.*s\nCIG: ", GD->TagLength, GD->Tag);
							if (GD->CigarLength > 0) {
								const unsigned char * restrict cigar = GD->Cigar;
								while (cigar[0] != kCIGARTERMINATOR) {
									const unsigned int N = cigar[0];
									const unsigned char type = cigar[1];
									fprintf(stderr, "%i%c ", N, (int) type);
									cigar += 2;
								}
							}
							fprintf(stderr, "\n Alignment start @ %u, length = %i, minpos = %u, maxpos = %u\n\n\n",
							        ND[iTag]->Location, ND[iTag]->AlignmentRange[1], minPos, maxPos);

							//assert(false);
							goto TaskDone;
						}
					}

#ifdef DEBUG_INFO
					{
						fprintf(stderr, "\n\n%12lu   ", 0UL);
						for (int j=0; j<TotalLength; j++) {
							fputc((int) GenomeBuffer[j], stderr);
						}
						fputc('\n', stderr);

						for (int iTag=0; iTag<nTags; iTag++) {
							fprintf(stderr, "%12lu %c ", ND[iTag]->Ordinal,(ND[iTag]->revNmer) ? '<' : '>');
							for (int j=0; j<TotalLength; j++) {
								fputc((int) Buffer[j*nTags+iTag], stderr);
							}
							fputc('\n', stderr);
						}
						fputs("\n\n", stderr);
					}
					{
					    fprintf(stderr, "\n\n%12lu   ", 0UL);
					    for (int j=0; j<TotalLength; j++) {
					        fputc((int) GenomeBuffer[j], stderr);
					    }
					    fputc('\n', stderr);

					    for (int iTag=0; iTag<nTags; iTag++) {
					        fprintf(stderr, "%12lu %c ", ND[iTag]->Ordinal,(ND[iTag]->revNmer) ? '<' : '>');
					        _Bool In = false;
									static const char *keys[] = {" ", ".", "_", "\033[32m+\033[0m", "\033[31m-\033[0m" };
					        for (int j=0; j<TotalLength; j++) {
					            const char * restrict ptr;
					            if (Buffer[j*nTags+iTag] == '.' && In)
					                In = false;
					            else if (Buffer[j*nTags+iTag] != '.' && Buffer[j*nTags+iTag] != ' ' && !In)
					                In = true;
					            switch (Buffer[j*nTags+iTag]) {
					            case '.': ptr = keys[0]; break;
					            case ' ': ptr = (GenomeBuffer[j] == ' ') ? ( In ? keys[1] : keys[0]) : keys[0]; break;
					            case '-': ptr = keys[4]; break;
					            default:
					                ptr = (GenomeBuffer[j] == ' ' ) ? keys[3] : keys[2];
					            }
					            fputs(ptr, stderr);
					        }
					        fputc('\n', stderr);
					    }
					    fputs("\n\n", stderr);
					}
#endif

					size_t WindowLength, WindowStart;
					{
						unsigned int k=0;
						const unsigned char * restrict GenomeBufferPtr = &GenomeBuffer[0];
						while (k<(task->GenomeStart-minPos)) {
							if (*GenomeBufferPtr++ != ' ') k++;
							assert((uintptr_t) GenomeBufferPtr < (uintptr_t) &GenomeBuffer[TotalLength]);
						}
						WindowStart = (uintptr_t) GenomeBufferPtr - (uintptr_t) &GenomeBuffer[0];
						const unsigned int Upto = ((maxPos < task->GenomeEnd) ? maxPos : task->GenomeEnd) - minPos;
						while (k<=Upto) {
							if (*GenomeBufferPtr++ != ' ') k++;
							assert((uintptr_t) GenomeBufferPtr < (uintptr_t) &GenomeBuffer[TotalLength]);
						}
						/* WARNING: I BELIEVE THIS IS USELESS !!! */
						while (*GenomeBufferPtr++ == ' ') { // include the last potential insertion
							k++;
							assert((uintptr_t) GenomeBufferPtr < (uintptr_t) &GenomeBuffer[TotalLength]);
						}
						const size_t WindowStop = (uintptr_t) GenomeBufferPtr - (uintptr_t) &GenomeBuffer[0];
						WindowLength = WindowStop - WindowStart - 1; // (WindowStop-1) - WindowStart + 1;
						assert(WindowLength > 0);
					}
					
					///////////////////////////////////////////////////////////////////////////////////////////////////
					// Realign everyone or not
					///////////////////////////////////////////////////////////////////////////////////////////////////
					unsigned int nKeptTags;
					/* Check and if needed alter memory size for alignment */
					const size_t PFResults_ld = (TotalLength + 15UL) & ~(15UL);
					{
						/* WORK and Matrix array */
						{
							size_t dummy = 1UL+TotalLength;
							if ( dummy > WORK_size) {
								WORK_size = dummy;
								if (WORK) _mm_free(WORK);
								size_t dummy2 = 4*dummy*sizeof(int) + 63;
								WORK = _mm_malloc(dummy2, 64);
								if (WORK == NULL) {
									fprintf(stderr, "Unable to allocate the WORK array for alignement, size = %lu\n", dummy);
									exit(1);
								}
								CellList = realloc(CellList, WORK_size*2*sizeof(short int));
								if (CellList == NULL) {
									fprintf(stderr, "Unable to allocate the CellList array for alignement, size = %lu\n", WORK_size*2*sizeof(short int));
									exit(1);
								}
							}

							dummy *= (1UL+kMaxReadLen);
							if (dummy > matrix_size) {
								if (matrix) _mm_free(matrix);
								matrix_size = dummy;
								matrix = _mm_malloc(dummy*sizeof(union lScores), 16);
								if (matrix == NULL) {
									fprintf(stderr, "Unable to allocate the MATRIX array for alignement, size = %lu\n", dummy);
									exit(1);
								}
							}
						}

						/* Resulting strings */
						const size_t lPFResults_size = PFResults_ld*nTags;
						if (lPFResults_size > PFResults_size) {
							PFResults_size = lPFResults_size;
							if (PFResults) _mm_free(PFResults);
							PFResults = (char *) _mm_malloc(lPFResults_size*sizeof(char),16);
							if (PFResults == NULL) {
								fprintf(stderr, "Unable to allocate pfSearch results' memory (%lu b)\n", lPFResults_size);
								exit(1);
							}
						}

						/* Scores */
						if (nTags > score_size) {
							score_size = nTags;
							score = (int*) realloc(score, nTags*sizeof(int));
							if (score == NULL) {
								fputs("Cannot allocate memory for scores\n", stderr);
								exit(1);
							}
						}

						memset(PFResults, ' ', lPFResults_size);
					}

					if (HasGap) {
						/* Generate profile */
						unsigned char FirstPass = 1;
						{
	#ifdef DEBUG_INFO
							snprintf(ProfileName, 128, "/tmp/genome_c%i_%u-%u.prf", options.chr, minPos, maxPos);
	#else
							snprintf(ProfileName, 128, "/tmp/genome_c%i_%u-%u_%lu.prf", options.chr, minPos, maxPos, tid);
	#endif

						NewPass:;
							out = fopen(ProfileName, "w");
							if (out == NULL) {
								fprintf(stderr, "unable to create profile %s\n", ProfileName);
								exit(1);
							}

							/* Header */
							fprintf(out, "ID   Genome_c%i_%u-%u; MATRIX.\n"
													 "AC   G%i_%u;\n"
													 "DT   MICMAP PRF Generator;\n"
													 "DE   Genome - chr %i pos. %u-%u;\n",
													 options.chr, minPos, maxPos, options.chr, minPos, options.chr, minPos, maxPos);
							fprintf(out, "MA   /GENERAL_SPEC: ALPHABET='ACGTN'; LENGTH=%u;\n"
													 "MA   /DISJOINT: DEFINITION=PROTECT; N1=%u; N2=%u;\n"
													 "MA   /NORMALIZATION: MODE=1; FUNCTION=LINEAR; R1=1.0; R2=1.0; TEXT='Dummy';\n"
													 "MA   /CUT_OFF: LEVEL=0; SCORE=%u; N_SCORE=12; MODE=1; TEXT='Dummy';\n",
													 TotalLength, 1, /*WindowStart+WindowLength*/ TotalLength-1 , 100);

#ifndef NDEBUG
							if (WindowStart+WindowLength > TotalLength) {
								fprintf(stderr, "OOPS: total length (%u) < WindowStart (%lu) + WindowLength (%lu) genome start (%u) minpos (%u) genome end (%u) maxpos (%u)\nGEN: %.*s\n",
												TotalLength, WindowStart, WindowLength, task->GenomeStart, minPos,
												task->GenomeEnd, maxPos, (int) TotalLength, GenomeBuffer); 
								assert(0);
							}
#endif

							/* Defaults values, similar to json file except inversion of INSERTION <-> DELETION
							*  "Scores": {
												"Match": 3,
												"MisMatch": -1,
												"MisMatchAtStart": -1,
												"MisMatchAtEnd": -1,
												"Insertion": 0,
												"Deletion": -1
								},
								"Transitions": {
												"MM": 1, "MI": -15, "MD": -3,
												"IM": 0, "II": 0, "ID": -67108864,
												"DM": 0, "DI": -67108864, "DD": 0
								}

							************************************************************************************************
							*                                   TRANSITIONS                                                *
							************************************************************************************************
							*
							* TABLE OF TRANSITIONS
								*
								* CODE    | x < NDIP1  | NDIP1 <= i < NDIP2  | NDIP2 <= i  | part of
								* =============================================================================================
								* _XM     | B1   + BM  |       B1 + BM       |  NLOW + BM  |
								* _XI     | B1   + BI  |       B1 + BI       |  NLOW + BI  | Insertion score EXTRA
								* _XD     | B1   + BD  |       B1 + BD       |  NLOW + BD  |
								*         |            |                     |             |
								* _MX     | NLOW + ME  |       E1 + ME       |   E1  + ME  | Insertion score MATCH
								* _IX     | NLOW + IE  |       E1 + IE       |   E1  + IE  | Insertion score INSERTION
								* _DX     | NLOW + DE  |       E1 + DE       |   E1  + DE  | Insertion score DELETION
								*---------+------------+---------------------+-------------+-----------------------------------
								* _YM     | B0 + BM    |       B0 + BM       |  NLOW + BM  |
								* _YI     | B0 + BI    |       B0 + BI       |  NLOW + BM  | FIRST SEQUENCE PROTEIN
								* _YD     | B0 + BD    |       B0 + BD       |  NLOW + BM  |
								*---------+------------+---------------------+-------------+-----------------------------------
								* _MY     | NLOW + ME  |       E0 + ME       |   E0  + ME  |
								* _IY     | NLOW + IE  |       E0 + IE       |   E0  + IE  | LAST SEQUENCE PROTEIN
								* _DY     | NLOW + DE  |       E0 + DE       |   E0  + DE  |
								*---------+------------+---------------------+-------------+-----------------------------------
								*/

							fprintf(out, "MA   /DEFAULT: I=-1;D=0; B0=0;B1=*; E0=0; E1=*; MM=2;MI=-3;MD=-12;ME=0; IM=0;II=-1;ID=*;IE=0; DM=-3;DI=*;DD=-1;DE=*;\n");

							const char * restrict column = Buffer;
							const char * restrict previousColumn = NULL;
							unsigned char InsertionStarted = 0;
							unsigned char DeletionStarted = 0;
							unsigned char InsertionCnt = 0;
							unsigned char DeletionCnt = 0;
							for (unsigned int iprf=0; iprf<TotalLength; iprf++) {
								int BestBase;
								int Bases[5] = {-1,-1,-1,-1,3};
								if (GenomeBuffer[iprf] != ' ') {
									BestBase = (int) GenomeBuffer[iprf];
									for (int k=0; k<4; k++) Bases[k] = (DNAAlphabet[k] == BestBase) ? 3 : -1;

									/* Looking for genome insertions */
									int HasDeletion = 0;
									int AcceptEnd = 0;
									int AcceptBeginDeletion = 0;
// 									for (int j=0; j<nTags; j++) {
// 										if (column[j] == '-') {
// 											HasDeletion = 1;
// 											if (previousColumn && previousColumn[j] != '-') AcceptBeginDeletion = 1;
// 										}
// 										else if (previousColumn && (previousColumn[j] == '-'))
// 											AcceptEnd = 1;
// 									}

									int doclose = 0;
									if (InsertionStarted) fprintf(out, "CC   ENDING TAGS INSERTIONS %u\n", (unsigned int) InsertionCnt++);
									if (HasDeletion && !DeletionStarted) fprintf(out, "CC   STARTING GENOME INSERTION %u\n", (unsigned int) DeletionCnt);
									if (DeletionStarted && AcceptBeginDeletion) fputs("CC   STARTING ANOTHER GENOME INSERTION\n", out);
									if (InsertionStarted || DeletionStarted || HasDeletion ) {
										fputs("MA   /I:", out);
										doclose=1;
									}
									if (AcceptEnd && DeletionStarted) fputs(" DM=0;", out);
									if (HasDeletion) {
										//fputs(" BD=0; ID=0;", out); 
										fputs(" BD=-1; ID=0; ", out);
// 										if (DeletionStarted)
// 											fputs("MD=0; DD=0; ", out);
// 										else
// 											fputs("MD=2;", out);
										if (AcceptBeginDeletion) fputs("MD=2; ", out);
										if (DeletionStarted) fputs("DD=0; ", out);
										DeletionStarted = 1;
									}
									else {
										DeletionStarted = 0;
									}

									if (InsertionStarted) fputs(" DM=1; DD=0; ", out);

									if (doclose) fputc('\n', out);

									InsertionStarted = 0;
									fprintf(out, "MA   /M: SY='%c'; M=%i,%i,%i,%i,%i;\n", BestBase, Bases[0], Bases[1], Bases[2], Bases[3], Bases[4]);
									if (AcceptEnd && DeletionStarted) fprintf(out, "CC   ENDING GENOME INSERTION %u\n", (unsigned int) DeletionCnt++);
								}
								else {
									int AcceptEnd = 0;
									for (int k=0; k<4; k++) Bases[k] = -100;
									for (int j=0; j<nTags; j++) {
										if (column[j] == ' ')
											AcceptEnd = 1;
										else
											/* It might be worth in the future to handle different mismatch cost within tag insertion based on the coverage */
											for (int k=0; k<4; k++) if (column[j] == DNAAlphabet[k]) Bases[k] = 3;
									}

									BestBase = (int) '?';

									if (!InsertionStarted) fprintf(out, "CC   STARTING TAGS INSERTIONS %u\n", (unsigned int) InsertionCnt);

									//fprintf(out, "MA   /I: BM=1; BD=0; ID=0;");
									//fprintf(out, "MA   /I: BM=-1; BD=0; ID=0;");
									fprintf(out, "MA   /I: BM=0; BD=-10; ID=0;"); 
									fputs(" MM=2;", out);

									if (InsertionStarted || DeletionStarted) 
										fputs(" MD=0; DD=0;", out);
									else 
										fputs(" MD=1;", out);
									if (AcceptEnd) fputs(" DM=0;", out);

									fprintf(out, "\nMA   /M: SY='%c'; M=%i,%i,%i,%i,%i;\n", BestBase, Bases[0], Bases[1], Bases[2], Bases[3], Bases[4]);

									InsertionStarted = 1;
									DeletionStarted = 0;
								}
								previousColumn = column;
								column += nTags;
							}

							fprintf(out, "MA   /I: B0=*; E1=0;\n//\n");
							if (fflush(out) != 0) {
								perror("Flushing profile");
								exit(1);
							}
						}
						
						/* Closing profile file */
						if (fclose(out) != 0) {
							fprintf(stderr,"Closing file %s failed\n", ProfileName);
							perror("Profile closing");
							exit(1);
						}
						
						if (ReadProfile(ProfileName, &profile, true, false) <= 0) {
							fprintf(stderr, "Error reading profile %s generated from interval %u--%u\n", ProfileName, task->GenomeStart, task->GenomeEnd);
							assert(0);
							exit(1);
						}

						/* Perform Alignements */
#ifdef DEBUG_INFO
						{
							fprintf(stderr, "Aligning back and forth %u tags...\n", nTags);
							unsigned int l;
							fputs("               ", stderr);
							for (l=minPos; l<(task->GenomeStart-1); l++) fputc(' ', stderr);
							fputc(']', stderr);
							while (l++ < task->GenomeEnd) fputc('=', stderr);
							fputs("[\n", stderr);
							fprintf(stderr, "               %.*s\n", (int) TotalLength, &GenomeBuffer[0]);
							
							fputs("               ", stderr);
							for (l=0; l<WindowStart-1; l++) fputc(' ', stderr);
							fputc(']', stderr);
							l = 0;
							while (l++ < WindowLength) fputc('=', stderr);
							fputs("[\n", stderr);
						}
#endif
						
						/* Reverse the profile */
						struct Profile * restrict const revPrf = ReverseProfile(&profile);
						if (PrepareExtraTable(revPrf) != 0) {
								fprintf(stderr, "Error preparing extra table in reverse profile %s\n", ProfileName);
								exit(1);
						}
						const unsigned char * const restrict Alphabet_Mapping = profile.Alphabet_Mapping;
						nKeptTags = 0U;
						for (unsigned int iTag=0U; iTag<nTags; iTag++) {
							const size_t SequenceLength = (size_t) ND[iTag]->TagLength;
							assert(SequenceLength < kMaxReadLen);
	//						fprintf(stderr, "Tag %u ", iTag);
							/* Transform Tag into the score indices */
							{
								int counter = 0;
								const unsigned char * const restrict CharPtr = ND[iTag]->Tag;
								for (size_t i=0; i<SequenceLength; ++i) {
									register size_t index = (size_t) ( ( CharPtr[i] >= (unsigned char) 'a' ) ? CharPtr[i] - ((unsigned char) 'a' - (unsigned char) 'A') : CharPtr[i] );
									assert( index >= (size_t) 'A' && index <= (size_t) 'Z');
									Indices[counter++] = Alphabet_Mapping[index - (size_t) 'A'];
								}
							}

							/* OK run the matrix creation */
							Standard_sse41.BuildMatrix(&profile, Indices, matrix, WORK, NULL, 0, SequenceLength);

							/* Get Best Alignment */
							const int nCol = 1 + TotalLength;
							int HasFork = -1;
							int revHasFork = -1;
							int iprf = -1;
							{
								int BestScore = NLOW;
								const union lScores * restrict Line = &matrix[SequenceLength*nCol];
								for (int iCol=1; iCol<nCol; iCol++) {
									const int lscore = GetScore(Line[iCol].Element[EXTRA], &Standard_sse41);
									if (lscore > BestScore) {
										BestScore = lscore;
										HasFork = -1;
										iprf = iCol;
									}
									else if (lscore == BestScore) {
											// Well the not being a deletion is probably not needed again, should check one day!
										if ( ((unsigned int)Line[iCol].Element[EXTRA] & (unsigned int) 0x3) != (unsigned int) PRIORITY_DELETION) HasFork = iCol;
									}
								}
								score[iTag] = BestScore;
	// 							printf("Best alignement ends at %i with score %i\n", iprf, BestScore); fflush(stdout);
							}

							if (iprf == -1) {
								fputs("Alignment found nothing better than NLOW !!!\n", stderr);
								exit(1);
							}

							/* Build the alignment sequence */
							{
								if (iprf < 1) {
	//								fprintf(stderr, " iprf=%i found\n", iprf);
									goto SKIP_THAT_ONE; // could be much more filtered
								}
								const unsigned char * restrict CharPtr = &(ND[iTag]->Tag[SequenceLength - 1]);
								const unsigned char * restrict QualityPtr = &(ND[iTag]->qs[SequenceLength - 1]);
								const char MinQuality = (char) options.filterQualValue;
								assert(nKeptTags < nTags);
								assert((nKeptTags*PFResults_ld+iprf-1) < PFResults_size);
								char * restrict Sequence = &PFResults[nKeptTags*PFResults_ld+iprf-1];
								size_t State = EXTRA;
								int index = SequenceLength;
#ifndef NDEBUG
								int backup = iprf;
#endif
								assert(index < TotalLength);
								short int (* restrict CellListPtr)[2] = CellList;
								while ( iprf >= 0 && index >= 0)
								{
									assert((index*nCol+iprf) >= 0);
									assert((index*nCol+iprf) < (TotalLength+1)*(SequenceLength+1));
									assert((uintptr_t) CellListPtr < (uintptr_t) CellList[WORK_size]);
									const unsigned int Move = matrix[index*nCol+iprf].Element[State] & 0x3;
									CellListPtr[0][0] = (short int) iprf;
									CellListPtr[0][1] = (short int) index;
									CellListPtr++;
									switch (Move)
									{
										case PRIORITY_MATCH:
											--index;
											--iprf;
											State = MATCH;
#ifndef NDEBUG
											if (Sequence < &PFResults[nKeptTags*PFResults_ld] && iprf >= 0){
												printf("%.*s\n%.*s\n", (int) TotalLength, GenomeBuffer, (int) TotalLength, &PFResults[nKeptTags*PFResults_ld]);
												printf("MATCH\n");
												exit(1);
											}
#endif
											assert((uintptr_t) CharPtr >= (uintptr_t) &(ND[iTag]->Tag[0]));
											*Sequence-- = (*QualityPtr-- >= MinQuality) ? *CharPtr : ' ';
											CharPtr--;
											break;
										case PRIORITY_INSERTION:
	//  										fprintf(stderr, "Got an insertion, skipping that tag\n");
											goto SKIP_THAT_ONE;
											break;
										case PRIORITY_DELETION:
											--iprf;
											State = DELETION;
#ifndef NDEBUG
											if (Sequence < &PFResults[nKeptTags*PFResults_ld] && iprf >= 0 ){
												printf("%.*s\n%.*s\n", (int) TotalLength, GenomeBuffer, (int) TotalLength, &PFResults[nKeptTags*PFResults_ld]);
												printf("DELETION %i(%i) %i\n", iprf, backup, index);
												exit(1);
											}
#endif
											*Sequence-- = '-';
											break;
										case PRIORITY_EXTRA:
											goto DONE;
											break;
									}
								}

								fputs("Potential issue as alignment does not start with an ENTRY point in std!\n", stderr);
								goto SKIP_THAT_ONE;
							DONE: ;
								
								/*
								* Attempt to find connection between the different possible alignments
								*/
								if (HasFork >= 0) {
									CellListPtr[0][0] = (short int) -1;
									CellListPtr = CellList;
									iprf = HasFork;
									HasFork = -1;
									index = SequenceLength;
									State = EXTRA;
									while ( iprf >= 0 && index >= 0)
									{
										register const short int liprf = iprf;
										while (CellListPtr[0][0] > liprf) {
											assert((uintptr_t) CellListPtr < (uintptr_t) CellList[WORK_size]);
											++CellListPtr;
										}
										// Insertion are prohibited, hence no need to loop on index.
										if (CellListPtr[0][0] == liprf && CellListPtr[0][1] == (short int) index) {
												HasFork = iprf;
												goto DONE2;
										}
										const unsigned int Move = matrix[index*nCol+iprf].Element[State] & 0x3;
										switch (Move)
										{
											case PRIORITY_MATCH:
												--index;
												--iprf;
												State = MATCH;
												break;
											case PRIORITY_DELETION:
												--iprf;
												State = DELETION;
												break;
											case PRIORITY_EXTRA:
												State = EXTRA;
												goto DONE2;
												break;
											case PRIORITY_INSERTION:
												fprintf(stderr, "A HasFork ended up with insertion, region %u-%u\n", task->GenomeStart, task->GenomeEnd);
												goto SKIP_THAT_ONE;
										}
									}
									DONE2:;
								}
							}

							/* Reverse the Tag indices */
							{
								unsigned char * restrict const CharPtr = Indices;
								unsigned char * restrict BackPtr = &Indices[SequenceLength-1];

								for (size_t i=0; i<SequenceLength/2; ++i) {
										const unsigned char c = CharPtr[i];
										CharPtr[i] = *BackPtr;
										*BackPtr-- = c;
								}
								if (SequenceLength & 0x1) {
										const unsigned char c = CharPtr[SequenceLength/2];
										CharPtr[SequenceLength/2] = *BackPtr;
										*BackPtr = c;
								}
							}

							/* OK run the matrix creation */
							Standard_sse41.BuildMatrix(revPrf, Indices, matrix, WORK, NULL, 0, SequenceLength);

							/* Get Best Alignment */
							iprf = -1;
							{
								int BestScore = NLOW;
								const union lScores * restrict Line = &matrix[SequenceLength*nCol];
								for (int iCol=1; iCol<nCol; iCol++) {
									const int lscore = GetScore(Line[iCol].Element[EXTRA], &Standard_sse41);
									if (lscore > BestScore) {
										BestScore = lscore;
										revHasFork = -1;
										iprf = iCol;
									}
									else if (lscore == BestScore) {
											// Again the deletion check might not be useful, but check !
										if ( ((unsigned int)Line[iCol].Element[EXTRA] & (unsigned int) 0x3) != (unsigned int) PRIORITY_DELETION) {
											revHasFork = iCol;
										}
									}
								}
								
								/*
								* Attempt to find connection between the different possible alignments
								*/
								if (revHasFork >= 0) {
									short int (* restrict CellListPtr)[2] = CellList;
									int index = SequenceLength;
									int State = EXTRA;
									while ( iprf >= 0 && index >= 0)
									{
										assert((uintptr_t) CellListPtr < (uintptr_t) CellList[WORK_size]);
										CellListPtr[0][0] = (short int) iprf;
										CellListPtr[0][1] = (short int) index;
										++CellListPtr;
										const unsigned int Move = matrix[index*nCol+iprf].Element[State] & 0x3;
										switch (Move)
										{
											case PRIORITY_MATCH:
												--index;
												--iprf;
												State = MATCH;
												break;
											case PRIORITY_DELETION:
												--iprf;
												State = DELETION;
												break;
											case PRIORITY_EXTRA:
												State = EXTRA;
												goto DONE3;
												break;
											case PRIORITY_INSERTION:
												fprintf(stderr, "A revHasFork ended up with insertion, region %u-%u\n", task->GenomeStart, task->GenomeEnd);
												goto SKIP_THAT_ONE;
										}
									}
									fputs("Potential issue as alignment does not start with an ENTRY point in std!\n", stderr);
									goto SKIP_THAT_ONE;
									DONE3:;	
									
									CellListPtr[0][0] = (short int) -1;
									CellListPtr = CellList;
									iprf = revHasFork;
									revHasFork = -1;
									index = SequenceLength;
									State = EXTRA;
									while ( iprf >= 0 && index >= 0)
									{
										register const short int liprf = iprf;
										while (CellListPtr[0][0] > liprf) {
											assert((uintptr_t) CellListPtr < (uintptr_t) CellList[WORK_size]);
											++CellListPtr;
										}
										// Insertion are prohibited, hence no need to loop on index.
										if (CellListPtr[0][0] == liprf && CellListPtr[0][1] == (short int) index) {
												revHasFork = iprf;
												goto DONE4;
										}
										const unsigned int Move = matrix[index*nCol+iprf].Element[State] & 0x3;
										switch (Move)
										{
											case PRIORITY_MATCH:
												--index;
												--iprf;
												State = MATCH;
												break;
											case PRIORITY_DELETION:
												--iprf;
												State = DELETION;
												break;
											case PRIORITY_EXTRA:
												State = EXTRA;
												goto DONE4;
												break;
											case PRIORITY_INSERTION:
												fprintf(stderr, "A revHasFork ended up with insertion, region %u-%u\n", task->GenomeStart, 
												task->GenomeEnd);
												goto SKIP_THAT_ONE;
											
										}
									}
									DONE4:;	
								}
							}

							if (iprf == -1) {
								fputs("Alignment found nothing better than NLOW !!!\n", stderr);
								exit(1);
							}

#ifdef DEBUG_INFO
							fprintf(stderr, "%12lu %c ", ND[iTag]->Ordinal,(ND[iTag]->revNmer) ? '<' : '>');
							if (HasFork >= 0) {
								if (revHasFork >= 0) {
										fprintf(stderr,"\033[31m%.*s\033[0m",
													(int) (TotalLength - revHasFork),  &PFResults[(nKeptTags)*PFResults_ld]);
										fprintf(stderr, "\033[32m%.*s\033[0m",
													(int) (HasFork - revHasFork),  &PFResults[(nKeptTags)*PFResults_ld + (TotalLength - revHasFork)]);
										fprintf(stderr, "\033[33m%.*s\033[0m\n",
													(int) (TotalLength - HasFork),  &PFResults[(nKeptTags)*PFResults_ld + HasFork]);
								}
								else {
										fprintf(stderr, "\033[32m%.*s\033[0m",
													(int) HasFork,  &PFResults[(nKeptTags)*PFResults_ld]);
										fprintf(stderr, "\033[33m%.*s\033[0m\n",
													(int) (TotalLength - HasFork),  &PFResults[(nKeptTags)*PFResults_ld + HasFork]);
								}
							}
							else {
								if (revHasFork >= 0) {
									fprintf(stderr,"\033[31m%.*s\033[0m",
													(int) (TotalLength - revHasFork),  &PFResults[(nKeptTags)*PFResults_ld]);
									fprintf(stderr, "\033[32m%.*s\033[0m\n",
													(int) revHasFork,  &PFResults[(nKeptTags)*PFResults_ld + (TotalLength - revHasFork)]);
								}
								else {
									fprintf(stderr, "\033[32m%.*s\033[0m\n",
													(int) TotalLength,  &PFResults[(nKeptTags)*PFResults_ld]);
								}
							}
#endif
							/* Remove the fork at the beginning */
							if (revHasFork >= 0) {
								memset(&PFResults[(nKeptTags)*PFResults_ld],' ', TotalLength - revHasFork);
							}
							
							/* Remove the fork at the end */
							if (HasFork >= 0) {
								memset(&PFResults[(nKeptTags)*PFResults_ld+HasFork],' ', TotalLength - HasFork);
							}

							nKeptTags++;
							continue;
							
							SKIP_THAT_ONE: ;

	// 						printf("ANA: %u-%u Removing tag %3u %s HasFork=%i revHasFork=%i\n", task->GenomeStart,
	// 						       task->GenomeEnd, iTag, ND[iTag]->Tag, HasFork, revHasFork);

							/* Clean before reusing, you idiot !!! */
							memset(&PFResults[nKeptTags*PFResults_ld], ' ', PFResults_ld);
						}
						/* Destroy profile */
						FreeProfile(revPrf, true);
						FreeProfile(&profile, false);

						/* Remove the profile */
#ifndef DEBUG_INFO
						remove(ProfileName);
#endif

						if (FirstPass && PerformTwoPass) {
							FirstPass = 0;
#ifdef DEBUG_INFO
							snprintf(ProfileName, 128, "/tmp/genome_c%i_%u-%u-secondpass.prf", options.chr, minPos, maxPos);
#else
							snprintf(ProfileName, 128, "/tmp/genome_c%i_%u-%u-secondpass-%lu.prf", options.chr, minPos, maxPos, tid);
#endif

// 							fprintf(stderr, "Do it again...\n");
							nTags = nKeptTags;

							/* Reload data into Buffer ( column based ) */
// 							memset(Buffer, '.', nTags*TotalLength*sizeof(char));
							for (int k=0; k<nKeptTags; k++) {
								//fprintf(stderr, "%.*s\n",  (int) PFResults_ld, &PFResults[k*PFResults_ld]);
								for (size_t l=0; l<TotalLength; l++) {
									Buffer[l*nTags+k] = ( PFResults[k*PFResults_ld+l] == ' ') ? '.' : PFResults[k*PFResults_ld+l];
// 									fputc(Buffer[l*nTags+k], stderr);
								}
// 								fputc('\n', stderr);
							}
							goto NewPass;
						}
					}
					else {
						nKeptTags = nTags;
						const char MinQuality = (char) options.filterQualValue;
						for (unsigned int iTag=0U; iTag<nTags; iTag++) {
							const size_t TagLength = ND[iTag]->TagLength;
							const unsigned char * const restrict CharPtr = ND[iTag]->Tag;
							const unsigned char * restrict QualityPtr = ND[iTag]->qs;
							const size_t loffset = ND[iTag]->Location - minPos;
							assert(loffset>=0 && loffset<PFResults_ld);
							char * const restrict copyto = &PFResults[iTag*PFResults_ld+loffset];
							assert(TagLength < kMaxReadLen);
							for (size_t i=0; i<TagLength; i++) {
								copyto[i] = (QualityPtr[i] >= MinQuality) ? CharPtr[i] : ' ';
							}
#ifdef DEBUG_INFO
							fprintf(stderr, "%12lu %c ", ND[iTag]->Ordinal,(ND[iTag]->revNmer) ? '<' : '>');
							fprintf(stderr, "\033[32m%.*s\033[0m\n",
													(int) TotalLength,  &PFResults[iTag*PFResults_ld]);
#endif
						}
					}
					
					///////////////////////////////////////////////////////////////////////////////////////////////////
					// Analysis
					///////////////////////////////////////////////////////////////////////////////////////////////////
					/* Do what you have to do .... */
// 					printf("ANA: Kept %u on %u-%u\n", nKeptTags, task->GenomeStart, task->GenomeEnd);
					if (nKeptTags > 0U) {
						unsigned int WindowStartGenomeAddress;
						const unsigned int InitialWindowStart = WindowStart;
						const unsigned int InitialWindowLength = WindowLength;
						struct ExtraData *nValuableBasesPerTag = malloc(nKeptTags * sizeof(struct ExtraData));
						
						/* Check Window borders, we cannot have neither insertions nor deletions cut in half */
						{
							/* Keeping genome extension ( insertion) */ 
							char * WindowToKeep = Buffer;
							for (unsigned int k=0; k<TotalLength; k++) WindowToKeep[k] = (GenomeBuffer[k] == ' ') ? 1 : 0;
							
							/* Keeping Genome Deletions */
							const char * restrict ptr = PFResults;
							for (unsigned int k=0; k<nKeptTags; k++) {
								for (unsigned int l=0; l<TotalLength; l++) WindowToKeep[l] |= (ptr[l] == '-') ? 1 : 0;
								ptr += PFResults_ld;
							}

							/* Right Border */
							unsigned int WindowEnd;
							{
								ssize_t dummy = WindowStart + WindowLength;
								while (dummy < TotalLength && WindowToKeep[dummy] == 1) { ++dummy;}
								WindowEnd = dummy - 1;
							}

							/* Left Border */
							{
								ssize_t dummy = WindowStart-1;
								while (dummy >= 0L && WindowToKeep[dummy] == 1) { --dummy;}
								dummy++;
// 								const int diff = WindowStart - (int) dummy; 
// 								WindowStartGenomeAddress = task->GenomeStart - diff;
								WindowStartGenomeAddress = task->GenomeStart;
								for (int k=WindowStart-1; k>(int) dummy; k--) {
									if (GenomeBuffer[k] != ' ') WindowStartGenomeAddress--;
								}
								WindowStart = dummy + 1;
							}
							assert(WindowEnd > WindowStart);
							WindowLength = WindowEnd - WindowStart + 1;
						}

						/* Make a transpose copy of the results */
						const size_t TransposePFResults_ld = ((size_t) nKeptTags + 15UL) & ~(15UL);
						{
							const size_t Size = TransposePFResults_ld*WindowLength;
							if (Size > TransposePFResults_size) {
								if (TransposePFResults) _mm_free(TransposePFResults);
								TransposePFResults_size = Size;
								TransposePFResults = (char *) _mm_malloc(Size*sizeof(char),16);
								if (TransposePFResults == NULL) {
									fprintf(stderr, "Unable to allocate memory for the transpose (%lu), %u tags, window of %zu\n", Size, nTags, WindowLength);
									exit(1);
								}
							}

							for (int k=0; k<nKeptTags; k++) {
								for (size_t l=0; l<WindowLength; l++) {
									TransposePFResults[l*TransposePFResults_ld+k] = PFResults[k*PFResults_ld+l+WindowStart];
								}
							}
						}

// 							for (unsigned int k=0; k<nKeptTags; k++) {
// 								for (unsigned int l=0; l<WindowStart; l++) fputc(' ', stdout);
// 								for (unsigned int l=0; l<WindowLength; l++) 
// 									fputc((int) TransposePFResults[l*TransposePFResults_ld+k], stdout);
// 								fputc('\n', stdout);
// 							}

						/* Cleaning the data according to a minimum of occurence (4) */
						unsigned char HasData[WindowLength];
						size_t nValuableBases = 0;
						for (size_t l=0; l<WindowLength; l++) {
							const char GenomeBase = GenomeBuffer[WindowStart+l];
							unsigned int Counters[5] = { 0U,0U,0U,0U,0U };
							char * restrict col = &TransposePFResults[l*TransposePFResults_ld];
							unsigned int Useless = 0U;
							for (size_t k=0; k<nKeptTags; k++) {
								register unsigned char value = col[k];
								if (value == '-') 
									Counters[4]++;
								else if (value == ' ' || value == 'N')
									++Useless;
								else {
									assert(((value>>1) & 0b11) <= 3);
									Counters[(value>>1) & 0b11]++;
								}
							}

							/* Adding a dynamic threshold based upon total number */
							const unsigned int total = Counters[0] + Counters[1] + Counters[2] + Counters[3] + Counters[4];
							const unsigned int MinDyn = (20*total)/100;
							for (int k=0;k<5;k++) if (Counters[k] < MinDyn) {
								Useless += Counters[k];
								Counters[k] = 0U;
							}

							for (int k=0;k<5;k++) if (Counters[k] < MINIMUM_APPEARANCE) { 
								Useless += Counters[k]; 
								Counters[k] = 0U;
							}

							if ((GenomeBase == ' ') && ((Counters[4] + MINIMUM_APPEARANCE) >= (nKeptTags - Useless))) Counters[4] = 0U;

							for (size_t k=0; k<nKeptTags; k++) {
								register unsigned char value = col[k];
								if (value == '-' && Counters[4] == 0U) col[k] = 'x';
								else if (value == GenomeBase && value != ' ') col[k] = '.';
								else if (value >= 'A' && Counters[(value>>1) & 0b11] == 0U) col[k] = '?';
							}

							/* Is it valuable or not as a whole ? */
							if (GenomeBase != ' ') {
								HasData[l] = '0';
								if (GenomeBase != 'N') {
									const unsigned char GenomeBaseIndex = (GenomeBase>>1) & 0b11;
									for (unsigned char k=0; k<5; k++) {
										if ( k != GenomeBaseIndex && Counters[k]) HasData[l] = '1'; 
									}
								}
							}
							else {
								HasData[l] = (Counters[0] || Counters[1] || Counters[2] || Counters[3] || Counters[4]) ? '1' : '0';
							}
							if (HasData[l] == '1') nValuableBases++; 
						}

						
#ifdef DEBUG_INFO
						{
							unsigned int l = 0U;
							fputs("               ", stderr);
							while (l++ < InitialWindowStart-1) fputc(' ', stderr);
							fputc(']', stderr);
							l = 0U;
							while (l++ < InitialWindowLength) fputc('=', stderr);
							fputs("[\n", stderr);
							
							fputs("               ", stderr);
							for (l=0; l<WindowStart-1; l++) fputc(' ', stderr);
							fputc(']', stderr);
							l = 0;
							while (l++ < WindowLength) fputc('-', stderr);
							fputs("[\n", stderr);
							fputs("               ", stderr);
							for(unsigned int k=0; k<WindowStart; k++) fputc(' ', stderr);
							fprintf(stderr, "%.*s\n", (int) WindowLength, HasData);
// 							for(unsigned int k=0; k<WindowStart; k++) fputc(' ', stdout);
// 							printf("%.*s\n", (int) WindowLength, &GenomeBuffer[WindowStart]);
// 							printf("Found %i valuable changes on %u regions\n", nValuableBases, nRegions);
						}
#endif
						/* Avoiding regions not starting in mandated genome range */
						{
							int offset = InitialWindowStart - WindowStart;
							if (HasData[offset] == '1' && HasData[offset-1] == '1') {
								while (offset < WindowLength && HasData[offset] == '1') ++offset;
							}
							else if (HasData[offset] == '1') {
								--offset;
							}
							while (offset >= 0 ) HasData[offset--] = '0';
							
							offset = (InitialWindowStart + InitialWindowLength) - WindowStart;
							while (offset < WindowLength && HasData[offset] == '1') ++offset;
							while (offset < WindowLength) HasData[offset++] = '0';
							
						}
						
						nValuableBases = 0U;
						for (size_t l=0; l<WindowLength; l++) if (HasData[l] == '1') ++nValuableBases;

						
						/* Counting regions with SNPs or InDels, correcting the useless insertion  */
						unsigned int nRegions = (HasData[0] == '1') ? 1U : 0U;
						{
							unsigned char previous = HasData[0];
							unsigned int k = 0U;
							const unsigned char * const restrict LocalGenomeWindow = &GenomeBuffer[WindowStart];
							while (k<WindowLength) {
								const char c = HasData[k];
								if (previous == '0' && c == '1') {
									/* Correcting the useless insertion to take them into account, not splitting deletions */
									int l = (int) k - 1;
									while (l >= 0 && LocalGenomeWindow[l] == ' ' && HasData[l] == '0') l--;
									if (HasData[l] == '1') {
										while(++l < k) { HasData[l] = '1'; nValuableBases++; } 
									}
									else
										nRegions++;
								}
								k++;
								previous = c;
							}
						}

						if (nValuableBases == 0LU) goto TaskDone;

#ifdef DEBUG_INFO
						{
							unsigned int l = 0U;
							fputs("               ", stderr);
							while (l++ < InitialWindowStart-1) fputc(' ', stderr);
							fputc(']', stderr);
							l = 0U;
							while (l++ < InitialWindowLength) fputc('=', stderr);
							fputs("[\n", stderr);
							
							fputs("               ", stderr);
							for (l=0; l<WindowStart-1; l++) fputc(' ', stderr);
							fputc(']', stderr);
							l = 0;
							while (l++ < WindowLength) fputc('-', stderr);
							fputs("[\n", stderr);
							fputs("               ", stderr);
							for(unsigned int k=0; k<WindowStart; k++) fputc(' ', stderr);
							fprintf(stderr, "%.*s\n", (int) WindowLength, HasData);
// 							for(unsigned int k=0; k<WindowStart; k++) fputc(' ', stdout);
// 							printf("%.*s\n", (int) WindowLength, &GenomeBuffer[WindowStart]);
// 							printf("Found %i valuable changes on %u regions\n", nValuableBases, nRegions);
						}
#endif
// 						printf("Before reg : %.*s\n", (int) WindowLength, HasData);
						struct RegionsData {
							unsigned int genomepos;
							int covered;
							unsigned short int length;
							unsigned short int windowpos;
						} RegData[/*nRegions*/kMaxReadLen/3];
						assert(nRegions < kMaxReadLen/3);
						char GenomeWindow[WindowLength];
// 						printf("After reg : %.*s\n", (int) WindowLength, HasData);
						{
							unsigned char previous = '0';
							struct RegionsData * restrict ptr = RegData;
							unsigned short int len = 0;
							unsigned int l = 0U;
							char * restrict GenomeWindowptr = GenomeWindow;
							for (unsigned int k=0; k<WindowLength; k++) {
								const char c = HasData[k];
								assert(GenomeWindowptr <= &GenomeWindow[WindowLength-1]);
								if (c == '1') *GenomeWindowptr++ = GenomeBuffer[WindowStart+k];
								if (c == '1' && previous == '0') {
									ptr->genomepos = WindowStartGenomeAddress + l;
									ptr->windowpos = k;
									ptr++;
									len = 1;
								}
								else if ( previous == '1' && c == '0') {
									ptr[-1].length = len;
								}
								else {
									len++;
								}
								previous = c;
								if (GenomeBuffer[WindowStart + k] != ' ') l++;
								ptr->covered = '0';
							}
							if (previous == '1') ptr[-1].length = len;
							// We will get back to that !!! assert(HasData[WindowLength] == '0');
						}
						
#ifdef DEBUG_INFO
						{
							for (unsigned int k=0; k<nRegions; k++) {
								fprintf(stderr, "%12u   ", RegData[k].genomepos);
								unsigned short int l = 0U;
								while (l++ < (WindowStart + RegData[k].windowpos)) fputc(' ', stderr);
								fputc('|', stderr);
								if (RegData[k].length > 1) {
									l = 2;
									while (l++ < RegData[k].length) fputc('-', stderr);
									fputc('|', stderr);
								}
								fputc('\n', stderr);
							}
						}
#endif

// 						printf("Before concat : %.*s\n", (int) WindowLength, HasData);
						/* Concatenating tags */
						const size_t Informations_ld = (nValuableBases+15UL) & ~(15UL);
						assert(nKeptTags*Informations_ld < Buffer_size);
						char * const restrict Informations = Buffer;
						char * restrict InformationsPtr = Informations;
						unsigned int nValuableTags = 0U;
						for (size_t k=0; k<nKeptTags; k++) {
							unsigned int pos = 0U;
							unsigned int valuable = 0U;
							for (size_t l=0; l<WindowLength; l++) {
								if (HasData[l] == '1') {
									const char currentBase = TransposePFResults[l*TransposePFResults_ld+k];
									InformationsPtr[pos++] = currentBase;
									if ( currentBase != ' ' && currentBase != 'x' && currentBase != '?') valuable++;
								}
							}
// 							printf("'%.3s'\t%.*s\n", InformationsPtr, (int) WindowLength, HasData);
							if (valuable) {
								nValuableBasesPerTag[nValuableTags].ptr = InformationsPtr;
								nValuableBasesPerTag[nValuableTags++].n = valuable;
								InformationsPtr += Informations_ld;
							}
						}

						if (nValuableTags == 0U) goto TaskDone;
						qsort(nValuableBasesPerTag, nValuableTags, sizeof(struct ExtraData), MyCompare);
						
#ifdef DEBUG_INFO
						fputs("Initial:\n", stderr);
						for (int k=0; k<nValuableTags; k++) {
							fprintf(stderr, "%.*s\t%u/%zu\n", (int) nValuableBases, nValuableBasesPerTag[k].ptr, nValuableBasesPerTag[k].n,nValuableBases);
						}
#endif
						unsigned int nSignificative = 1U;
						struct Ordering Significatives[nValuableTags];
						/* Grouping identical ones */
						{
							Significatives[0].data = nValuableBasesPerTag[0];
							Significatives[0].count = 1U;
							for (unsigned int k=1; k<nValuableTags; k++) {
								int found = 0;
								for (unsigned int l=0; l<nSignificative; l++) { 
									if (strncmp(nValuableBasesPerTag[k].ptr, Significatives[l].data.ptr, (size_t) nValuableBases) == 0) {
										Significatives[l].count++;
										found = 1;
										break;
									}
								}
								
								if (!found) {
										Significatives[nSignificative].data = nValuableBasesPerTag[k];
										Significatives[nSignificative++].count = 1U;
								}
							}
							qsort(Significatives, nSignificative, sizeof(struct Ordering), MyOrdering);
						}
#ifdef DEBUG_INFO
						fputs("\nGouping identical\n", stderr);
						for (int k=0; k<nSignificative; k++) {
							fprintf(stderr, "%.*s\t%u/%zu\t%u\n",(int) nValuableBases, Significatives[k].data.ptr, Significatives[k].data.n, nValuableBases, Significatives[k].count);
						}
#endif
						/* Place partial that uniquely match a category with no error */
						{
							for (unsigned int k=1; k<nSignificative; k++) {
								const char * const restrict Test = Significatives[k].data.ptr;
								if (Test == NULL) continue;
								unsigned int Best = 0U;
								int id = -1;
								for (unsigned int l=0; l<k; l++) {
									const char * const restrict Src = Significatives[l].data.ptr;
									if (Src) {
										unsigned int count = 0U;
										
										for (size_t m=0; m<nValuableBases; m++) {
											unsigned int keepIf;
											keepIf  = (Test[m] == '?' && Src[m] != '-');
											keepIf |= (Test[m] == Src[m] && Test[m] != ' ' && Test[m] != 'x');
											if (keepIf) ++count;
										}
										if (Best < count) { Best = count; id = l; }
										else if (Best == count) { id = -1; }
		// 								printf("%.*s %.*s %u %i\n", nValuableBases, Test, nValuableBases, Src, Best, id);
									}
								}
								if (id >= 0 && Best >= Significatives[k].data.n) {
									Significatives[id].count += Significatives[k].count;
									Significatives[k].data.ptr = NULL;
									Significatives[k].data.n = 0U;
									Significatives[k].count = 0U;
								}
							}
							
							qsort(Significatives, nSignificative, sizeof(struct Ordering), MyOrdering);
						}
						
#ifdef DEBUG_INFO
						{
							fputs("\nPlace partial that uniquely match a category with no error\n", stderr);
							int count = 0U;
							const float invnTags = 1.0f/(float) nValuableTags;
							for (int k=0; k<nSignificative; k++) {
								if (Significatives[k].data.ptr) {
									fprintf(stderr, "%.*s\t%u/%zu\t%u\t%f\n", (int) nValuableBases, Significatives[k].data.ptr, Significatives[k].data.n, nValuableBases, Significatives[k].count, (float) Significatives[k].count * invnTags);
									count += Significatives[k].count;
								}
							}
							fprintf(stderr, "Total tags %u/%u\n", count, nValuableTags);
						}
#endif
						
						/* Assign those below MAF (%u) to groups with up to 1 error if sufficient valuable bases */
						if (nValuableBases > 5) {
							const unsigned int MAF = (unsigned int) (0.15f*(float) nValuableTags); 
							for (unsigned int k=1; k<nSignificative; k++) {
								const char * const Test = Significatives[k].data.ptr;
								if (Test == NULL || Significatives[k].count >= MAF) continue;
								unsigned int Best = 0U;
								int id = -1;
								for (unsigned int l=0; l<k; l++) {
									const char * const Src = Significatives[l].data.ptr;
									if (Src) {
										unsigned int count = 0U;
										
										for (size_t m=0; m<nValuableBases; m++) {
											unsigned int keepIf;
											keepIf  = (Test[m] == '?' && Src[m] != '-');
											keepIf |= (Test[m] == Src[m] && Test[m] != ' ' && Test[m] != 'x');
											if (keepIf) ++count;
										}
		// 								 printf("%.23s %.23s %u\n", Test, Src, count);
										if (Best < count) { Best = count; id = l;}
										else if (Best == count) { id = -1; }
									}
								}
								if (id >= 0 && (1+Best) >= Significatives[k].data.n) {
									Significatives[id].count += Significatives[k].count;
									Significatives[k].data.ptr = NULL;
									Significatives[k].data.n = 0U;
									Significatives[k].count = 0U;
								}
							}
							
							qsort(Significatives, nSignificative, sizeof(struct Ordering), MyOrdering);
#ifdef DEBUG_INFO
							{
								fprintf(stderr, "\nAssign those below MAF (%u) to groups with up to 1 error\n", MAF);
								int count = 0U;
								const float invnTags = 1.0f/(float) nValuableTags;
								for (int k=0; k<nSignificative; k++) {
									if (Significatives[k].data.ptr) {
										fprintf(stderr, "%.*s\t%u/%zu\t%u\t%f\n", (int) nValuableBases, Significatives[k].data.ptr, Significatives[k].data.n, nValuableBases, Significatives[k].count, (float) Significatives[k].count * invnTags);
										count += Significatives[k].count;
									}
								}
								fprintf(stderr, "Total tags %u/%u\n", count, nValuableTags);
							}
#endif
						}

						/* Place partial that uniquely match a category with no error */
						{
							for (unsigned int k=1; k<nSignificative; k++) {
								const char * const restrict Test = Significatives[k].data.ptr;
								if (Test == NULL) continue;
								unsigned int Best = 0U;
								int id = -1;
								for (unsigned int l=0; l<k; l++) {
									const char * const restrict Src = Significatives[l].data.ptr;
									if (Src) {
										unsigned int count = 0U;
										
										for (size_t m=0; m<nValuableBases; m++) {
											unsigned int keepIf;
											keepIf  = (Test[m] == '?' && Src[m] != '-');
											keepIf |= (Test[m] == Src[m] && Test[m] != ' ' && Test[m] != 'x');
											if (keepIf) ++count;
										}
										if (Best < count) { Best = count; id = l; }
										else if (Best == count) { id = -1; }
									}
								}
								if (id >= 0 && Best >= Significatives[k].data.n) {
									Significatives[id].count += Significatives[k].count;
									Significatives[k].data.ptr = NULL;
									Significatives[k].data.n = 0U;
									Significatives[k].count = 0U;
								}
							}
							
							qsort(Significatives, nSignificative, sizeof(struct Ordering), MAFOrdering);
						}
						
#ifdef DEBUG_INFO
						{
							fputs("\nPlace partial that uniquely match a category with no error\n", stderr);
							int count = 0U;
							const float invnTags = 1.0f/(float) nValuableTags;
							for (int k=0; k<nSignificative; k++) {
								if (Significatives[k].data.ptr) {
									fprintf(stderr, "%.*s\t%u/%zu\t%u\t%f\n", (int) nValuableBases, Significatives[k].data.ptr, Significatives[k].data.n, nValuableBases, Significatives[k].count, (float) Significatives[k].count * invnTags);
									count += Significatives[k].count;
								}
							}
							fprintf(stderr, "Total tags %u/%u\n", count, nValuableTags);
						}
#endif
						
						/* Keeping... */
						int AllOK = 1;
						int AvailableAlleles = 0;
						char MAFs[64];
						MAFs[0]='(';
						{
							IF_DEBUG_INFO(fprintf(stderr, "\nKeeping MAF >= 0.15 ...\n\n"));
							assert(nValuableTags > 0);
							const float invnTags = 1.0f/(float) nValuableTags;
							int count = 0;
							char * restrict pMafs = &MAFs[1];
							for (int k=0; k<nSignificative; k++) {
								if (Significatives[k].data.ptr ) {
									const float val = (float) Significatives[k].count * invnTags;
									if (val >= 0.15f) {
										IF_DEBUG_INFO(fprintf(stderr, "%.*s\t%u/%zu\t%u\t%f @ %u", (int) nValuableBases, Significatives[k].data.ptr,Significatives[k].data.n,      nValuableBases, Significatives[k].count, val, k));
										if (k < 10) pMafs += sprintf(pMafs, "%4.2f,", val);
										if (++count > 2) {
											int pos=0;
											int AllCovered = 1;
											for (unsigned int iRegion=0; iRegion<nRegions; iRegion++) {
												int covered = 0;
												RegData[iRegion].covered = '0';
												{
													int isSpace = 1;
													const char * restrict const tptr = Significatives[k].data.ptr + pos;
													for (int l=0; l<RegData[iRegion].length; l++) isSpace &= (tptr[l] == ' ');
													if (isSpace) {
														covered = 1;
														RegData[iRegion].covered = '1';
														goto Decide;
													}
												}
												for (int l=0; l<2; l++) {
													if (strncmp(Significatives[k].data.ptr + pos, Significatives[l].data.ptr + pos, RegData[iRegion].length) == 0) {
														RegData[iRegion].covered = '1';
														covered = 1;
													}
												}
											Decide:;
												if (!covered) {
													AllCovered = 0;
													break;
												}
												pos += RegData[iRegion].length;
											}
											if (!AllCovered) {
												IF_DEBUG_INFO(fputs(" Possibilitied NOT covered by kept allele !!!\n", stderr));
												AllOK = 0;
											}
											else {
												IF_DEBUG_INFO(fputs(" All fine :)\n", stderr));
											}
										}
										else { 
											IF_DEBUG_INFO(fputc('\n', stderr));
										}
									}
								}
							}
							AvailableAlleles = count;
							
							assert(pMafs < &MAFs[50]);
							if (pMafs[-1] == ',') {
								pMafs[-1] = ')';
								pMafs[0] = '\0';
							}
							else {
								pMafs[0] = ')';
								pMafs[1] = '\0';
							}
						}
						
#if TAG_ANALYSIS_MODE == NIRVANA
						/* Generating VCF...  */
						bcf1_t * const restrict rec = task->VCFRecord;
						assert(rec);
						ngt = bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt_arr);
						assert(ngt == 2);
						
						unsigned int CurrentTypeID = 0U;
						if (!bcf_is_snp(rec)) {
							if (rec->d.var[bcf_gt_allele(gt[0])].n < 0 || rec->d.var[bcf_gt_allele(gt[1])].n < 0) {
								CurrentTypeID = 1;
							}
							else {
								CurrentTypeID = 2;
							}
						}
						unsigned int RealMemoryPosition = 1 + rec->pos;
						if (CurrentTypeID != 0) RealMemoryPosition += 1; /* God knows why, this is necessary !!! */
						
						const int VCFIsHeterozygote = (bcf_gt_allele(gt[0]) != bcf_gt_allele(gt[1])); 
						int FoundVCFLocation = 0;
#else
						int FoundVCFLocation = 1;
						IF_DEBUG_INFO(fprintf(stderr, "\nVCFs:\n"));
#endif
						if (AvailableAlleles) {
							int nPos = 0;
							for (unsigned int iRegion=0; iRegion<nRegions; iRegion++) {
								/* Correcting borders */
								int len = RegData[iRegion].length;
								int k=nPos; assert(nPos < nValuableBases);
								IF_DEBUG_INFO(fprintf(stdout, "Region %i, alleles=%i, pos=%u, length = %i\n", iRegion, AvailableAlleles, RegData[iRegion].genomepos + k, len));
								if (len > 1) {
									if (AvailableAlleles > 1) {
										while (len > 0 && Significatives[0].data.ptr[k] == '.' && Significatives[1].data.ptr[k] == '.') {
											++k; --len; assert(k<nValuableBases || len <= 0);
										}
										while (len > 0 && Significatives[0].data.ptr[k+len-1] == '.' && Significatives[1].data.ptr[k+len-1] == '.') {
											--len; assert(len >= 0);
										}
									}
									else {
										while (len > 0 && Significatives[0].data.ptr[k] == '.') {
											++k; --len; assert(k < nValuableBases || len <= 0);
										}
										while (len > 0 && Significatives[0].data.ptr[k+len-1] == '.') --len;
									}
									IF_DEBUG_INFO(fprintf(stderr, "Region %i, border correction: pos=%u, length=%u\n", iRegion, RegData[iRegion].genomepos + k, len));
								}
								
								if (AllOK) RegData[iRegion].covered = '1';
								
								assert( (k+len) <= nValuableBases);
								unsigned char DoTest = 0;
#if TAG_ANALYSIS_MODE == NIRVANA
								if ( (RegData[iRegion].genomepos + (k - nPos)) == RealMemoryPosition) DoTest++;
								if ( (RegData[iRegion].genomepos + (k - nPos)) <= RealMemoryPosition
									&& (RegData[iRegion].genomepos + (k - nPos) + len) >= RealMemoryPosition)
									DoTest++;
#else
								DoTest = 1;
#endif
								if (DoTest) {
									FoundVCFLocation = 1;
									/* Remove empty strings */
									int Heterozygote = 0;
									{
										int l=0;
										while(l<len) {
											assert((k+l) < nValuableBases);
											if (Significatives[0].data.ptr[k+l] != ' ') break; else  ++l;
										}
										const int removeFirst = (l>=len); 
										if (AvailableAlleles > 1) {
											l=0;
											while(l<len) {
												assert((k+l) < nValuableBases);
												if (Significatives[1].data.ptr[k+l] != ' ') break; else ++l;
											}
											if (l>=len ) {
												if (removeFirst) {
													IF_DEBUG_INFO(fprintf(stdout, "Region %i: first allele removed, second with space %i>=%i\n", iRegion, l, len));
													goto NextRegion;
												}
											}
										}
										else {
											if (removeFirst) {
												IF_DEBUG_INFO(fprintf(stdout, "Region %i: Single allele with space %i>=%i\n", iRegion, l, len));
												goto NextRegion;
											}
										}
										
										if (removeFirst) {
											l=0;
											while(l<len) { 
												assert((k+l) < nValuableBases);
												Significatives[0].data.ptr[k+l] = Significatives[1].data.ptr[k+l];
												l++;
											}
										}
										else if (AvailableAlleles > 1) {
												if (strncmp(&(Significatives[0].data.ptr[k]), &(Significatives[1].data.ptr[k]), len) != 0)
													Heterozygote = 1;
										}
									}
									
									IF_DEBUG_INFO(printf("Region %i AL 0 : \'%.*s\'\t", iRegion, len, &(Significatives[0].data.ptr[k])));
									IF_DEBUG_INFO(if (AvailableAlleles > 1) printf("AL 1 : \'%.*s\'\t", len, &(Significatives[1].data.ptr[k])));
									IF_DEBUG_INFO(printf("Ext. genome : \'%.*s\'\n", len, &GenomeWindow[k]));
									IF_DEBUG_INFO(printf("Region %i heterozygote = %u\n", iRegion, (unsigned int) Heterozygote));
									
									const char * const restrict Gen = &GenomeWindow[k];
									char GT[4] = { '0', '/', '!', '\0' };
									assert(len < MAX_ALLELE_SIZE);
									char Previously[len+2];
									char * restrict ptrAllele = Allele0;
									int Second = 0;
									const char * restrict Src = &Significatives[0].data.ptr[k];
									int InDel = (Heterozygote && Significatives[1].data.ptr[k] == '-') ? 1 : 0;
									char * restrict PosBackup = Allele0;
								NextAllele: ;
									int l=0;
									if (len > 1 || Gen[0] == ' ' || Src[0] == '-' || InDel) { ptrAllele++; InDel = 1; }
									unsigned int nPoints=0U, nDashesOrX=0U;
									while(l<len) {
										if (Src[l] == '.') {
											assert (ptrAllele < (PosBackup+MAX_ALLELE_SIZE));
											*ptrAllele++ = Gen[l];
											++nPoints;
										}
										else if (Src[l] != '-' && Src[l] != 'x') {
											assert (ptrAllele < (PosBackup+MAX_ALLELE_SIZE));
											*ptrAllele++ = Src[l];
										}
										else if (Gen[l] == ' ')
											++nDashesOrX;
										++l;
									}
									if ((nPoints == len) || (nDashesOrX == len)) {
										//ptrAllele = PosBackup;
										GT[2*Second] = '0';
									}
									else
										GT[2*Second] = GT[0] + 1;
									
									assert (ptrAllele < (PosBackup+MAX_ALLELE_SIZE));
									*ptrAllele = '\0';
									
									IF_DEBUG_INFO(printf("ALLELE %i %s\n", Second, ptrAllele-len));

									if (Heterozygote && !Second) {
										PosBackup = Allele1;
										ptrAllele = Allele1;
										Src = &Significatives[1].data.ptr[k];
										Second = 1;
										goto NextAllele;
									}
									else if (!Heterozygote)
										GT[2] = GT[0];

									/* Build alleles' string */
									{
										int m=0;
										l=-1;
										const unsigned char * const restrict ptr = &GenomeBuffer[WindowStart + RegData[iRegion].windowpos + (k-nPos)];
										if (InDel) {
											while (ptr[l] == ' ') --l;
											Previously[m++] = ptr[l];
											Allele0[0] = ptr[l];
											Allele1[0] = ptr[l];
										}
										l=0;
										while (l<len) if (ptr[l++] != ' ') Previously[m++] = ptr[l-1];
										Previously[m] = '\0';
									}
									
									Heterozygote = (GT[0] != GT[2]);
									unsigned int ok = 1;
									int ErrorId = 0;
#if TAG_ANALYSIS_MODE == NIRVANA
									const unsigned int Position = RegData[iRegion].genomepos + (k - nPos) - InDel;
									const unsigned int VCFPosition = rec->pos + 1U;
									if (Heterozygote != VCFIsHeterozygote) {
										/* check that the VCF info is not in both alleles */
										const char * restrict SrcPtr = rec->d.allele[bcf_gt_allele(gt[1])];
										size_t vcfLen = strlen(SrcPtr);
										int check = 1;
										if (!VCFIsHeterozygote) {
											if (VCFPosition<Position)
												check = 0;
											else {
												size_t l = 0UL;
												do {
													assert(((VCFPosition - Position) + l) < MAX_ALLELE_SIZE); // len does not have the previous character accounted
													if (Allele0[(VCFPosition - Position) + l] != Allele1[(VCFPosition - Position) + l]) check = 0;
												} while (++l<vcfLen); 
											}
										}
										if (check == 0 || VCFIsHeterozygote) {
											ErrorId = 1;
											atomic_fetch_add(&nZygoteError, 1);
											ok = 0;
										}
									}
									
									if (Heterozygote && ok) {
										int found = 0;
										if (Position <= VCFPosition) {
											/* Seeking Allele 0 */
											const char * restrict SrcPtr = rec->d.allele[bcf_gt_allele(gt[0])];
											size_t vcfLen = strlen(SrcPtr);
											int check0 = 1, check1 = 1;
											size_t l = 0UL;
											do {
												assert(((VCFPosition - Position) + l) < MAX_ALLELE_SIZE); // len does not have the previous character accounted
												if (Allele0[(VCFPosition - Position) + l] != SrcPtr[l]) check0 = 0;
												if (Allele1[(VCFPosition - Position) + l] != SrcPtr[l]) check1 = 0;
											} while (++l<vcfLen); 
											found += check0;
											found += check1;
											
											if (found == 0) {
												ok = 0;
												ErrorId = 2;
												atomic_fetch_add(&nAllele0NotFound, 1);
											}
											
											/* Seeking Allele 1 */
											found = 0;
											SrcPtr = rec->d.allele[bcf_gt_allele(gt[1])];
											vcfLen = strlen(SrcPtr);
											
											check0 = 1, check1 = 1;
											l = 0UL;
											do {
												assert(((VCFPosition - Position) + l) < MAX_ALLELE_SIZE); // len does not have the previous character accounted
												if (Allele0[(VCFPosition - Position) + l] != SrcPtr[l]) check0 = 0;
												if (Allele1[(VCFPosition - Position) + l] != SrcPtr[l]) check1 = 0;
											} while (++l<vcfLen); 
											found += check0;
											found += check1;
											
											if (found == 0) {
												ok = 0;
												ErrorId = 3;
												atomic_fetch_add(&nAllele1NotFound, 1);
											}
										}
									}
									
									if (ok) {
										atomic_fetch_add(&nRecordsFound, 1);
										if (outfile) {
											fprintf(outfile, "%u\t%u\t.\t%s\t%s/%s\t50\t%s:%u MATCHES %u\t%s\t%s\t%s/%s\t%u/%u\n",
															options.chr, RegData[iRegion].genomepos + (k - nPos) - InDel, Previously, Allele0, Allele1, GT, nValuableTags,
															RealMemoryPosition, types[CurrentTypeID], rec->d.allele[0], rec->d.allele[bcf_gt_allele(gt[0])], rec->d.allele[bcf_gt_allele(gt[1])], bcf_gt_allele(gt[0]), bcf_gt_allele(gt[1]));
										}
									}
									else {
										if (Allele0[0] != '\0' && !ok) {
											fprintf(stdout, "\x1b[2K\x1b[GError %s : %u\t%u\t.\t%s\t%s/%s\t50\t%s:%u is not %u\t%s\t%s\t%s/%s\t%u/%u\n",
															Errors[ErrorId],
															options.chr, RegData[iRegion].genomepos + (k - nPos) - InDel, Previously, Allele0, Allele1, GT, nValuableTags,
															RealMemoryPosition, types[CurrentTypeID], rec->d.allele[0], rec->d.allele[bcf_gt_allele(gt[0])], rec->d.allele[bcf_gt_allele(gt[1])], bcf_gt_allele(gt[0]), bcf_gt_allele(gt[1]));
											fflush(stdout);
										}
									}
									break;
#else
									char Info[64];
									snprintf(Info, 64, "AlleleOK=%i,nAlleles=%u,MAF=%s", AllOK, AvailableAlleles, MAFs); 
									/* VCF printf("%s\t%u\t.\t%.*s\t%c\t255\tPASS\t.\tGT:DP:AD\t%s:%u:%u,%u\n", */
									IF_DEBUG_INFO(printf("Region start=%u, k=%u, nPos=%u, Indel=%u\n", 
									                     RegData[iRegion].genomepos, k, nPos, InDel));
									const unsigned int Position =  RegData[iRegion].genomepos + (k - nPos) - InDel;
									if ( Position <= task->GenomeEnd && Position >= task->GenomeStart) {
										if (Heterozygote) {
											if (GT[0] != '0' && GT[2] != '0')
												fprintf(stdout, "%s\t%u\t.\t%s\t%s,%s\t50\tPASS\t%s,RegionOK=%c,Window=%u-%u\tGT:DP:AD\t%s:%u:%u,%u\n",
												        Genome.virtchr[options.chr].AC, Position, Previously, Allele0, Allele1, Info, RegData[iRegion].covered, task->GenomeStart, task->GenomeEnd, GT, nValuableTags, Significatives[0].count, Significatives[1].count);
											else {
												const char * const ctmp = (GT[0] != '0') ? Allele0 : Allele1;
												const int itmp = (GT[0] != '0') ? 1 : 0;
												fprintf(stdout, "%s\t%u\t.\t%s\t%s\t50\tPASS\t%s,RegionOK=%c,Window=%u-%u\tGT:DP:AD\t%s:%u:%u,%u\n",
												        Genome.virtchr[options.chr].AC, Position, Previously, ctmp, Info, RegData[iRegion].covered, task->GenomeStart, task->GenomeEnd, GT, nValuableTags,
												        Significatives[1-itmp].count, Significatives[itmp].count);
											}
										}
										else {
											if (GT[0] != '0') {
												fprintf(stdout, "%s\t%u\t.\t%s\t%s\t50\tPASS\t%s,RegionOK=%c,Window=%u-%u\tGT:DP:AD\t%s:%u:0,%u\n",
																Genome.virtchr[options.chr].AC, Position, Previously, Allele0, Info, RegData[iRegion].covered, task->GenomeStart, task->GenomeEnd, GT, nValuableTags,
																Significatives[0].count);
											}
#ifdef DEBUG_INFO
											else {
												fprintf(stdout, "GEE, I will miss region %i !!!\n\tGT: %c/%c\n\tAllele 0: %s\n\tAllele 1: %s\n",
																iRegion, (int) GT[0], (int) GT[2], Allele0, Allele1);
											}
#endif
										}
									}
#ifdef DEBUG_INFO
									else {
											fprintf(stdout, "Region %i out of bound %u > %u > %u\n", iRegion, task->GenomeStart, Position, task->GenomeEnd);
									}
#endif
#endif
								}
#ifdef DEBUG_INFO	
							goto noprint;
							NextRegion:
								fprintf(stdout, "Region %i was skipped\n", iRegion); 
							noprint:
								nPos += RegData[iRegion].length;
#else
							NextRegion:
								nPos += RegData[iRegion].length;
#endif
							}
						}
#if TAG_ANALYSIS_MODE == NIRVANA
						if (FoundVCFLocation == 0) {
							atomic_fetch_add(&nLocationNotFound, 1);
							if (types[CurrentTypeID] != 0) RealMemoryPosition--;
							fprintf(stdout, "\x1b[2K\x1b[GError location not found %u\t%s\t%s\t%s/%s\t%u/%u\n", RealMemoryPosition, types[CurrentTypeID], rec->d.allele[0], rec->d.allele[bcf_gt_allele(gt[0])], rec->d.allele[bcf_gt_allele(gt[1])], bcf_gt_allele(gt[0]), bcf_gt_allele(gt[1]));
						}
#endif
						free(nValuableBasesPerTag);
					}
#if TAG_ANALYSIS_MODE == NIRVANA
					else {
						atomic_fetch_add(&nRemovedHasFork, 1);
						bcf1_t * const restrict rec = task->VCFRecord;
						assert(rec);
						ngt = bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt_arr);
						assert(ngt == 2);
						
						unsigned int CurrentTypeID = 0U;
						if (!bcf_is_snp(rec)) {
							if (rec->d.var[bcf_gt_allele(gt[0])].n < 0 || rec->d.var[bcf_gt_allele(gt[1])].n < 0) {
								CurrentTypeID = 1;
							}
							else {
								CurrentTypeID = 2;
							}
						}
						unsigned int RealMemoryPosition = 1 + rec->pos;
						if (CurrentTypeID != 0) RealMemoryPosition += 1; /* God knows why, this is necessary !!! */
							
						fprintf(stdout, "\x1b[2K\x1b[GError data forked %u\t%s\t%s\t%s/%s\t%u/%u\n", RealMemoryPosition, types[CurrentTypeID], rec->d.allele[0], rec->d.allele[bcf_gt_allele(gt[0])], rec->d.allele[bcf_gt_allele(gt[1])], bcf_gt_allele(gt[0]), bcf_gt_allele(gt[1]));
					}
#elif defined(DEBUG_INFO)
					else {
							fprintf(stderr, "No more tags !!!\n");
					}
#endif
				}
			
			TaskDone:
				/* Free Tags allocated space */
				for (unsigned int k=0; k<nTags; k++) free(task->data[k]);

#if TAG_ANALYSIS_MODE == NIRVANA
				atomic_fetch_add(&nVCFRecord, 1);
				/* clear bcf record */
				bcf_clear(task->VCFRecord);
#endif
				/* put back the Interval memory slot to the available queue */
				pthread_mutex_lock(&thpool->donequeue_p.rwmutex);
				jobqueue_push(&thpool->donequeue_p, (job_t*) task);
				pthread_mutex_unlock(&thpool->donequeue_p.rwmutex);
			}
		
			/* Signal I am no longer working */
			pthread_mutex_lock(&thpool->thcount_lock);
			thpool->num_threads_working--;
			pthread_mutex_unlock(&thpool->thcount_lock);
		}
	}
	
	/* Signal I am done */
	pthread_mutex_lock(&thpool->thcount_lock);
	thpool->num_threads_alive--;
	pthread_mutex_unlock(&thpool->thcount_lock);
	
	if (PFResults) _mm_free(PFResults);
	if (TransposePFResults) _mm_free(TransposePFResults);
	if (score) free(score);
	if (Buffer) free(Buffer);
	if (matrix) _mm_free(matrix);
	if (WORK) _mm_free(WORK);
	if (CellList) free(CellList);
		
	return (void*) 0;
}
//---------------------------------------------------------------

/* vim: tabstop=2 shiftwidth=2
 */
