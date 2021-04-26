/*
 * ------------------------------------------------------------------------------------------------------------------------
 *
 *                                * micmap *
 *      Mapping of short reads in fastq format onto a reference
 *
 *
 *  Copyright (C) SIB  - Swiss Institute of Bioinformatics,                2015-2019 Nicolas Guex, Thierry Schuepbach and Christian Iseli
 *  Copyright (C) UNIL - University of Lausanne, Switzerland                    2019 Nicolas Guex and Christian Iseli
 *  Copyright (C) EPFL - Ecole Polytechnique Fédérale de Lausanne, Switzerland  2021 Christian Iseli
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
 *      Contacts:   Nicolas.Guex@unil.ch and christian.iseli@epfl.ch
 *      Repository: https://github.com/sib-swiss/micmap
 *
 * ------------------------------------------------------------------------------------------------------------------------
 */
#ifndef ANALYZE_FULL_OUTPUT
#define INFO_SIZE 128
#define SUM_UP(dest, src)
#else
#define INFO_SIZE 512
#define SUM_UP(dest, src) {\
	size_t m = 0UL;\
	do {\
		__m128i __sum = _mm_load_si128((__m128i*) ((dest).Accounting + 8UL*m));\
		__sum = _mm_adds_epu16(__sum, *((__m128i*) ((src).Accounting + 8UL*m)));\
		_mm_store_si128((__m128i*) ((dest).Accounting + 8UL*m), __sum);\
	} while (++m < nRegions4);\
	m = 0UL;\
	do {\
		__m128i __sum = _mm_load_si128((__m128i*) ((dest).PHred + 8UL*m));\
		__sum = _mm_adds_epu16(__sum, *((__m128i*) ((src).PHred + 8UL*m)));\
		_mm_store_si128((__m128i*) ((dest).PHred + 8UL*m), __sum);\
	} while (++m < nValuableBases2);\
}
#endif
{
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Required Memory
	char * restrict PFResults = NULL;
	char * restrict TransposePFResults = NULL;
	char * restrict GenomeBuffer = NULL;
	char * restrict Buffer = NULL;
	short int (*CellList)[2] = NULL;
	int * score = NULL;
	union lScores * restrict matrix = NULL;
	int * restrict WORK = NULL;
	FILE * out;
	struct Deletions {
		unsigned int Length;
		unsigned int count;
	} * restrict Deletions = NULL;
#ifdef ANALYZE_FULL_OUTPUT
	unsigned short int * restrict RegionsStat = NULL;
#endif

	size_t PFResults_size = 0UL;
	size_t TransposePFResults_size = 0UL;
	size_t score_size = 0UL;
	size_t GenomeBuffer_size = 0UL;
	size_t Buffer_size = 0UL;
	size_t matrix_size = 0UL;
	size_t WORK_size = 0UL;
	size_t Deletions_size = 0UL;
#ifdef ANALYZE_FULL_OUTPUT
	size_t RegionsStat_Size = 0UL;
#endif

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
	const float MAF = ((float) ((OPTIONS*) (thpool->common))->minorAllelePct) * 0.01f;
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
				const unsigned int nTags = task->count;
				assert(nTags >= 4);
#ifdef DEBUG_INFO
 				printf("ANA: GOT %u tags for %u-%u\n", nTags,  task->GenomeStart, task->GenomeEnd);
#endif
				unsigned int minPos = -1;
				unsigned int maxPos = 0U;
				for (unsigned int iTag=0U; iTag<nTags; iTag++) {
					unsigned int uitmp = ND[iTag]->Location;
					if (uitmp < minPos) minPos = uitmp;
					uitmp += ND[iTag]->AlignmentRange[1];
					if (uitmp > maxPos) maxPos = uitmp;
				}

				/* Extra margin of both sides of 2 DNA bases */
				maxPos += 2;

				assert(maxPos > minPos);
				assert(minPos <= task->GenomeEnd);
				assert(maxPos >= task->GenomeStart);
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
					if ((size_t) Length > Deletions_size) {
						Deletions = (struct Deletions*) realloc(Deletions, Length*sizeof(struct Deletions));
						if (Deletions == NULL) {
							fprintf(stderr, "Unable to allocation Deletions of size %u in AnalyzeTags3\n", Length);
							exit(1);
						}
						Deletions_size = (size_t) Length;
					}
					unsigned int Insertions[2*Length];
					unsigned int * const restrict Memory = Insertions + Length;
					memset(Insertions, 0, 2*Length*sizeof(unsigned int));
					memset(Deletions, 0, Length*sizeof(struct Deletions));

					for (int iTag=0; iTag<nTags; iTag++) {
						GTLRawData_t * restrict GD = ND[iTag];
// 						unsigned int * restrict ptr = Deletions;
						assert(GD->Location >= minPos);
						size_t pos = GD->Location - minPos;
						if (GD->CigarLength > 1) {
							const unsigned char * restrict cigar = GD->Cigar;
							unsigned int TagPos = 0U;
							assert(GD->TagLength > 10);
							const unsigned int TagPosMax = GD->TagLength - MIN_BORDER_DISTANCE;
							while (cigar[0] != kCIGARTERMINATOR) {
								const unsigned int N = cigar[0];
								const unsigned char type = cigar[1];
								assert(pos <= Length);
								if (type == 'D') {
									if ((GD->MismatchLength <= 1) || (N < 20 && (TagPos >= MIN_BORDER_DISTANCE) && ((TagPos+N) <= TagPosMax))) {
										if (N > Deletions[pos].Length) Deletions[pos].Length = N;
										Deletions[pos].count++;
									}
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
						if (Deletions[l].Length > 0 && Insertions[l] > 0) c = (int) '*';
						else if (Deletions[l].Length > 0) c = (int) '-';
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
						const char * restrict cptr = (char *) Genome.table + index;
						char * restrict cptro = GenomeBuffer;

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
							Buffer = (char*) realloc(Buffer, needed_size*sizeof(char));
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
						const char * restrict cptr = (char *) GD->Tag;
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
						fprintf(stderr, "\n\n%12s   ", "GENOME");
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
									static const char *mods[] = {"\033[36mA\033[0m", "\033[31mC\033[0m", "\033[35mG\033[0m", "\033[33mT\033[0m" };
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
											case 'A':
					                ptr = (GenomeBuffer[j] == ' ') ? keys[3] : (GenomeBuffer[j] == Buffer[j*nTags+iTag]) ? keys[2] : mods[0]; break;
											case 'C':
					                ptr = (GenomeBuffer[j] == ' ') ? keys[3] : (GenomeBuffer[j] == Buffer[j*nTags+iTag]) ? keys[2] : mods[1]; break;
											case 'G':
					                ptr = (GenomeBuffer[j] == ' ') ? keys[3] : (GenomeBuffer[j] == Buffer[j*nTags+iTag]) ? keys[2] : mods[2]; break;
											case 'T':
					                ptr = (GenomeBuffer[j] == ' ') ? keys[3] : (GenomeBuffer[j] == Buffer[j*nTags+iTag]) ? keys[2] : mods[3]; break;
					            default:
					                ptr = (GenomeBuffer[j] == ' ') ? keys[3] : keys[2];
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
						unsigned int k = 0U;
						const char * restrict GenomeBufferPtr = &GenomeBuffer[0];
						assert((task->GenomeStart-minPos) <= TotalLength);
						while (k<(task->GenomeStart-minPos)) {
							if (*GenomeBufferPtr++ != ' ') k++;
							assert((uintptr_t) GenomeBufferPtr < (uintptr_t) &GenomeBuffer[TotalLength]);
						}
						WindowStart = (uintptr_t) GenomeBufferPtr - (uintptr_t) &GenomeBuffer[0];
						const unsigned int Upto = ((maxPos < task->GenomeEnd) ? maxPos : task->GenomeEnd) - minPos;
						while (k<=Upto) {
							if (*GenomeBufferPtr++ != ' ') k++;
							assert((uintptr_t) GenomeBufferPtr <= (uintptr_t) &GenomeBuffer[TotalLength]);
						}
						/* WARNING: I BELIEVE THIS IS USELESS !!! */
// 						while (*GenomeBufferPtr++ == ' ') { // include the last potential insertion
// 							assert((uintptr_t) GenomeBufferPtr <= (uintptr_t) &GenomeBuffer[TotalLength]);
// 						}
						const size_t WindowStop = (uintptr_t) GenomeBufferPtr - (uintptr_t) &GenomeBuffer[0];
						WindowLength = WindowStop - WindowStart - 1; // (WindowStop-1) - WindowStart + 1;
						assert(WindowLength > 0);
					}

					///////////////////////////////////////////////////////////////////////////////////////////////////
					// Realign everyone or not
					///////////////////////////////////////////////////////////////////////////////////////////////////
					unsigned int nKeptTags;
					/* Check and if needed alter memory size for alignment */
#ifdef ANALYZE_FULL_OUTPUT
#define EXTRA_COLUMN 1UL
#else
#define EXTRA_COLUMN 0UL
#endif
					const size_t PFResults_ld = (TotalLength + EXTRA_COLUMN +15UL) & ~(15UL);
#undef EXTRA_COLUMN
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
						const size_t lPFResults_size = (PFResults_ld*nTags + 63) & ~(63);
						if (lPFResults_size > PFResults_size) {
							PFResults_size = lPFResults_size;
							if (PFResults) _mm_free(PFResults);
#ifndef ANALYZE_FULL_OUTPUT
							PFResults = (char *) _mm_malloc(lPFResults_size*sizeof(char),64);
#else
							PFResults = (char *) _mm_malloc(2*lPFResults_size*sizeof(char),64);
#endif
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
#ifndef ANALYZE_FULL_OUTPUT
						memset(PFResults, ' ', lPFResults_size);
#else
						memset(PFResults, ' ', 2*lPFResults_size);
#endif
					}

					if (HasGap) {
						/* Generate profile */
// 						unsigned char FirstPass = 1;
						{
#ifdef DEBUG_INFO
							snprintf(ProfileName, 128, "/tmp/genome_ext_c%i_%u-%u.prf", options.chr, minPos, maxPos);
#else
							snprintf(ProfileName, 128, "/tmp/genome_ext_c%i_%u-%u_%lu.prf", options.chr, minPos, maxPos, tid);
#endif

// 						NewPass:;
							out = fopen(ProfileName, "w");
							if (out == NULL) {
								fprintf(stderr, "unable to create profile %s\n", ProfileName);
								exit(1);
							}

							/* Header */
							fprintf(out, "ID   Genome_ext_c%i_%u-%u; MATRIX.\n"
							             "AC   G%i_%u;\n"
							             "DT   MICMAP PRF Generator;\n"
							             "DE   Genome - chr %i pos. %u-%u;\n",
							             options.chr, minPos, maxPos, options.chr, minPos, options.chr, minPos, maxPos);
							fprintf(out, "MA   /GENERAL_SPEC: ALPHABET='%s'; LENGTH=%u;\n"
							             "MA   /DISJOINT: DEFINITION=PROTECT; N1=%u; N2=%u;\n"
							             "MA   /NORMALIZATION: MODE=1; FUNCTION=LINEAR; R1=1.0; R2=1.0; TEXT='Dummy';\n"
							             "MA   /CUT_OFF: LEVEL=0; SCORE=%u; N_SCORE=12; MODE=1; TEXT='Dummy';\n",
							             DNAAlphabet, TotalLength, 1, /*WindowStart+WindowLength*/ TotalLength-1 , 100);

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

							fprintf(out, "MA   /DEFAULT: I=-1;D=0; B0=0;B1=*; E0=0; E1=*; MM=2;MI=-3;MD=-12;ME=0; IM=*;II=%i;ID=*;IE=0; DM=-3;DI=*;DD=-1;DE=*;\n", NLOW);

							const char * restrict column = Buffer;
							const char * restrict previousColumn = NULL;
							unsigned char InsertionStarted = 0;
							unsigned char DeletionStarted = 0;
							int NextDeletionEnd = -1;
							int RealPosition = 0;

							for (unsigned int iprf=0; iprf<TotalLength; iprf++) {
								int Bases[5] = {-1,-1,-1,-1,3};
								/* Matches and Deletions paths */
								if (GenomeBuffer[iprf] != ' ') {
									const char BestBase = GenomeBuffer[iprf];
									for (int k=0; k<4; k++) Bases[k] = (DNAAlphabet[k] == BestBase) ? 3 : -1;

									/* Do we have an open genome deletion to close */
									int DeletionEnd = 0;
									if (RealPosition == NextDeletionEnd) {
										DeletionEnd = 1;
										DeletionStarted = 0;
									}

									/* Looking for genome deletion */
									unsigned char BeginDeletion = 0;
									unsigned char FavorMatch = 0;
									if (Deletions[RealPosition].count > 0) {
										BeginDeletion = 1;
										const int uidummy = RealPosition + (int) Deletions[RealPosition].Length;
										NextDeletionEnd = (uidummy>NextDeletionEnd) ? uidummy : NextDeletionEnd;
 										if ((Deletions[RealPosition].Length == 1) && Deletions[RealPosition].count < 8) FavorMatch = 1;
										fprintf(out, "CC   STARTING GENOME DELETION FOR %u BASES, current %u, real %u, stop @ %u, favorizing match %u\n", Deletions[RealPosition].Length, iprf, RealPosition, NextDeletionEnd, (unsigned int) FavorMatch);
									}
									unsigned char HasDeletion = 0;
									int IntermediateDeletionEnd = 0;
 									for (int j=0; j<nTags; j++) {
										if (column[j] == '-') HasDeletion = 1;
 										else if (previousColumn && (previousColumn[j] == '-')) IntermediateDeletionEnd = 1;
 									}

									fputs("MA   /I: ", out);
									if (DeletionEnd || IntermediateDeletionEnd) fputs("DM=0; ", out);
									if (DeletionStarted) fputs("DD=0; ", out);
									if (BeginDeletion) {
										fputs("MD=2; ", out);
										DeletionStarted = 1;
									}

									//if (InsertionStarted) fputs(" DM=1; DD=0; ", out);
									if (InsertionStarted) {
										fputs("IM=2; DM=0; ", out);
										if (DeletionStarted) FavorMatch = 1;
									}

									fprintf(out,"MI=%i; I=%i,%i,%i,%i,%i;\n", NLOW, NLOW, NLOW, NLOW, NLOW, NLOW);

									InsertionStarted = 0;
									if (FavorMatch) {
										/* The idea is to avoid deletion when the purpose is to have snps */
										fprintf(out, "MA   /M: SY='%c'; M=4,4,4,4,3;\n", BestBase);
									}
									else {
										fprintf(out, "MA   /M: SY='%c'; M=%i,%i,%i,%i,%i;\n", BestBase, Bases[0], Bases[1], Bases[2], Bases[3], Bases[4]);
									}
// 									if (InsertionEnd) fputs("CC   ENDING GENOME INSERTION\n", out);
									if (DeletionEnd) fputs("CC   ENDING GENOME DELETION\n", out);
									RealPosition++;
								}
								/* Insertions paths */
								else {
									int AcceptEnd = 0;
									Bases[0] = -100;
									Bases[1] = -100;
									Bases[2] = -100;
									Bases[3] = -100;

									for (size_t j=0; j<nTags; j++) {
										register const char value = column[j];

										if (value == ' ') {
											if (previousColumn && previousColumn[j] != ' ') AcceptEnd = 1;
										}
										else {
											if (DNAAlphabet[0] == value) Bases[0] = 3;
											if (DNAAlphabet[1] == value) Bases[1] = 3;
											if (DNAAlphabet[2] == value) Bases[2] = 3;
											if (DNAAlphabet[3] == value) Bases[3] = 3;
										}
									}

// 									if (Bases[0] == -100 && Bases[1] == -100 && Bases[2] == -100 && Bases[3] == -100) {
// 										printf("Insertion at location %i has not data in %s: \'", iprf, DNAAlphabet);
// 										for (int j=0; j<nTags; j++) fputc((int) column[j], stdout);
// 										fputs("\'\n", stdout);
// 									}
									if (!InsertionStarted) fputs("CC   STARTING TAGS INSERTIONS\n", out);

									fprintf(out, "MA   /I: BI=0; BM=%i; BD=%i; MM=%i; ", NLOW, NLOW, NLOW);
									fprintf(out, "MI=%i; ", (InsertionStarted) ? NLOW : 2);
									fprintf(out, "II=%i; ", (InsertionStarted) ? 2 : NLOW);
									fprintf(out, "ID=%i; ", (AcceptEnd) ? 0 : NLOW);
									fprintf(out, "DI=%i; ", (DeletionStarted) ? 0 : NLOW);
									if (DeletionStarted)  {
										if (!InsertionStarted) fputs("MD=2; ", out);
										fprintf(out, "DM=%i; DD=0; ", NLOW);
									}
									else {
										if (!InsertionStarted)
											fputs("MD=2; ", out);
										else
											fputs("DD=0; ", out);
									}
									if (InsertionStarted) fprintf(out, "DE=%i; ", NLOW);

									fprintf(out, "I=%i,%i,%i,%i,%i;\n", Bases[0], Bases[1], Bases[2], Bases[3], Bases[4]);

									/*if (InsertionStarted || DeletionStarted)
										fputs(" MD=0; DD=0;", out);
									else
										fputs(" MD=1;", out);
									if (AcceptEnd) fputs(" DM=0;", out);*/

									fprintf(out, "MA   /M: SY='?'; M=%i,%i,%i,%i,%i;\n", NLOW, NLOW, NLOW, NLOW, NLOW);

									InsertionStarted = 1;
									//DeletionStarted = 0;
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
#ifdef DEBUG_INFO
							fprintf(stderr, "Tag %u ", iTag);
#endif
							/* Transform Tag into the score indices */
							{
								int counter = 0;
								const char * const restrict CharPtr = (char *) ND[iTag]->Tag;
								for (size_t i=0; i<SequenceLength; ++i) {
									register size_t index = (size_t) ( ( CharPtr[i] >= (unsigned char) 'a' ) ? CharPtr[i] - ((unsigned char) 'a' - (unsigned char) 'A') : CharPtr[i] );
									assert( index >= (size_t) 'A' && index <= (size_t) 'Z');
									Indices[counter++] = Alphabet_Mapping[index - (size_t) 'A'];
								}
							}

							/* OK run the matrix creation */
							ExtendedProfile_sse41.BuildMatrix(&profile, Indices, matrix, WORK, NULL, 0, SequenceLength);

							/* Get Best Alignment */
							const int nCol = 1 + TotalLength;
							int HasFork = -1;
							int revHasFork = -1;
							int iprf = -1;
							{
								int BestScore = NLOW;
								const union lScores * restrict Line = &matrix[(SequenceLength)*nCol];
								for (int iCol=1; iCol<nCol; iCol++) {
									const int lscore = GetScore(Line[iCol].Element[EXTRA], &ExtendedProfile_sse41);
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
#ifdef DEBUG_INFO
 								printf("Best alignement ends at %i with score %i\n", iprf, BestScore); fflush(stdout);
#endif
							}

							if (iprf == -1) {
								fputs("Alignment found nothing better than NLOW !!!\n", stderr);
#ifdef DEBUG_INFO
								const char From[] = "DMIX";
								for (int iCol=0; iCol<nCol; iCol++) {
									printf("%i\t%8i %c\t%8i %c\t%8i %c\t%8i %c\n", iCol,
									       GetScore(matrix[(0)*nCol+iCol].Element[MATCH], &ExtendedProfile_sse41),
												 (int) From[matrix[(0)*nCol+iCol].Element[MATCH] & 0x3],
									       GetScore(matrix[(0)*nCol+iCol].Element[INSERTION], &ExtendedProfile_sse41),
												 (int) From[matrix[(0)*nCol+iCol].Element[INSERTION] & 0x3],
												 GetScore(matrix[(0)*nCol+iCol].Element[DELETION], &ExtendedProfile_sse41),
												 (int) From[matrix[(0)*nCol+iCol].Element[DELETION] & 0x3],
												 GetScore(matrix[(0)*nCol+iCol].Element[EXTRA], &ExtendedProfile_sse41),
												 (int) From[matrix[(0)*nCol+iCol].Element[EXTRA] & 0x3]);
								}
								exit(1);
#endif
							}

							/* Build the alignment sequence */
							{
								if (iprf < 1) {
#ifdef DEBUG_INFO
									fprintf(stderr, " iprf=%i found\n", iprf);
#endif
									goto SKIP_THAT_ONE; // could be much more filtered
								}
								const char * restrict CharPtr = (char *) &(ND[iTag]->Tag[SequenceLength - 1]);
								const char * restrict QualityPtr = (char *) &(ND[iTag]->qs[SequenceLength - 1]);
								const char MinQuality = (char) options.filterQualValue;
								assert(nKeptTags < nTags);
								assert((nKeptTags*PFResults_ld+iprf-1) < PFResults_size);
								char * restrict Sequence = &PFResults[nKeptTags*PFResults_ld+iprf-1]; // INVESTIGATE MINUS ONE
								size_t State = EXTRA;
								int index = SequenceLength;
#ifndef NDEBUG
								int backup = iprf;
#endif
								assert(index < TotalLength);
								short int (* restrict CellListPtr)[2] = CellList;
								while ( iprf >= 0 && index >= 0) {
									assert((index*nCol+iprf) >= 0);
									assert((index*nCol+iprf) < (TotalLength+1)*(SequenceLength+1));
									assert((uintptr_t) CellListPtr < (uintptr_t) CellList[WORK_size]);
									const unsigned int Move = matrix[index*nCol+iprf].Element[State] & 0x3;
									CellListPtr[0][0] = (short int) iprf;
									CellListPtr[0][1] = (short int) index;
									CellListPtr++;
									switch (Move) {
										case PRIORITY_MATCH:
											--index;
											--iprf;
											State = MATCH;
#ifndef NDEBUG
											if (Sequence < &PFResults[nKeptTags*PFResults_ld] && iprf >= 0) {
												printf("%.*s\n%.*s\n", (int) TotalLength, GenomeBuffer, (int) TotalLength, &PFResults[nKeptTags*PFResults_ld]);
												printf("MATCH\n");
												exit(1);
											}
#endif
											assert((uintptr_t) CharPtr >= (uintptr_t) &(ND[iTag]->Tag[0]));
#ifdef ANALYZE_FULL_OUTPUT
											Sequence[PFResults_size] = *QualityPtr;
#endif
											*Sequence-- = (*QualityPtr-- >= MinQuality) ? *CharPtr : ' ';
											CharPtr--;
											break;
										case PRIORITY_INSERTION:
											--index;
											--iprf;
											State = INSERTION;
#ifndef NDEBUG
											if (Sequence < &PFResults[nKeptTags*PFResults_ld] && iprf >= 0) {
															printf("%.*s\n%.*s\n", (int) TotalLength, GenomeBuffer, (int) TotalLength, &PFResults[nKeptTags*PFResults_ld]);
															printf("MATCH\n");
															exit(1);
											}
#endif
											assert((uintptr_t) CharPtr >= (uintptr_t) &(ND[iTag]->Tag[0]));
#ifdef ANALYZE_FULL_OUTPUT
											Sequence[PFResults_size] = *QualityPtr;
#endif
											*Sequence-- = (*QualityPtr-- >= MinQuality) ? *CharPtr : ' ';
											CharPtr--;
											break;
										case PRIORITY_DELETION:
											--iprf;
											State = DELETION;
#ifndef NDEBUG
											if (Sequence < &PFResults[nKeptTags*PFResults_ld] && iprf >= 0 ) {
												printf("%.*s\n%.*s\n", (int) TotalLength, GenomeBuffer, (int) TotalLength, &PFResults[nKeptTags*PFResults_ld]);
												printf("DELETION %i(%i) %i\n", iprf, backup, index);
												exit(1);
											}
#endif
#ifdef ANALYZE_FULL_OUTPUT
											Sequence[PFResults_size] = '\0';
#endif
											*Sequence-- = '-';
											break;
										case PRIORITY_EXTRA:
// #ifndef NDEBUG
// 											if (Sequence < &PFResults[nKeptTags*PFResults_ld] && iprf >= 0) {
// 												printf("%.*s\n%.*s\n", (int) TotalLength, GenomeBuffer, (int) TotalLength, &PFResults[nKeptTags*PFResults_ld]);
// 												printf("MATCH\n");
// 												exit(1);
// 											}
// #endif
// 											assert((uintptr_t) CharPtr >= (uintptr_t) &(ND[iTag]->Tag[0]));
// 											*Sequence-- = (*QualityPtr-- >= MinQuality) ? *CharPtr : ' ';
// 											CharPtr--;
											goto DONE;
											break;
									}
								}

								fputs("Potential issue as alignment does not start with an ENTRY point in std!\n", stderr);
#ifdef DEBUG_INFO
								const char From[] = "DMIX";
								for (int iCol=0; iCol<nCol; iCol++) {
									printf("%i\t%8i %c\t%8i %c\t%8i %c\t%8i %c\n", iCol,
									       GetScore(matrix[(iCol)*nCol+iCol].Element[MATCH], &ExtendedProfile_sse41),
												 (int) From[matrix[(iCol)*nCol+iCol].Element[MATCH] & 0x3],
									       GetScore(matrix[(iCol)*nCol+iCol].Element[INSERTION], &ExtendedProfile_sse41),
												 (int) From[matrix[(iCol)*nCol+iCol].Element[INSERTION] & 0x3],
												 GetScore(matrix[(iCol)*nCol+iCol].Element[DELETION], &ExtendedProfile_sse41),
												 (int) From[matrix[(iCol)*nCol+iCol].Element[DELETION] & 0x3],
												 GetScore(matrix[(iCol)*nCol+iCol].Element[EXTRA], &ExtendedProfile_sse41),
												 (int) From[matrix[(iCol)*nCol+iCol].Element[EXTRA] & 0x3]);
								}
								exit(1);
#endif
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
											case PRIORITY_INSERTION:
												--index;
												--iprf;
												State = INSERTION;
												break;
											case PRIORITY_DELETION:
												--iprf;
												State = DELETION;
												break;
											case PRIORITY_EXTRA:
												State = EXTRA;
												goto DONE2;
												break;
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
							ExtendedProfile_sse41.BuildMatrix(revPrf, Indices, matrix, WORK, NULL, 0, SequenceLength);

							/* Get Best Alignment */
							iprf = -1;
							{
								int BestScore = NLOW;
								const union lScores * restrict Line = &matrix[(SequenceLength-1)*nCol];
								for (int iCol=1; iCol<nCol; iCol++) {
									const int lscore = GetScore(Line[iCol].Element[EXTRA], &ExtendedProfile_sse41);
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
												--index;
												--iprf;
												State = INSERTION;
												break;
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
												--index;
												--iprf;
												State = INSERTION;
												break;
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
										fprintf(stderr, "\033[33m%.*s\033[0m",
													(int) (TotalLength - HasFork),  &PFResults[(nKeptTags)*PFResults_ld + HasFork]);
								}
								else {
										fprintf(stderr, "\033[32m%.*s\033[0m",
													(int) HasFork,  &PFResults[(nKeptTags)*PFResults_ld]);
										fprintf(stderr, "\033[33m%.*s\033[0m",
													(int) (TotalLength - HasFork),  &PFResults[(nKeptTags)*PFResults_ld + HasFork]);
								}
							}
							else {
								if (revHasFork >= 0) {
									fprintf(stderr,"\033[31m%.*s\033[0m",
													(int) (TotalLength - revHasFork),  &PFResults[(nKeptTags)*PFResults_ld]);
									fprintf(stderr, "\033[32m%.*s\033[0m",
													(int) revHasFork,  &PFResults[(nKeptTags)*PFResults_ld + (TotalLength - revHasFork)]);
								}
								else {
									fprintf(stderr, "\033[32m%.*s\033[0m",
													(int) TotalLength,  &PFResults[(nKeptTags)*PFResults_ld]);
								}
							}
							fprintf(stderr,"\t%i\n", score[iTag]);
#endif
							/* Remove the fork at the beginning */
							if (revHasFork >= 0) {
								memset(&PFResults[(nKeptTags)*PFResults_ld],' ', TotalLength - revHasFork);
							}

							/* Remove the fork at the end */
							if (HasFork >= 0) {
								memset(&PFResults[(nKeptTags)*PFResults_ld+HasFork],' ', TotalLength - HasFork);
							}

#ifdef ANALYZE_FULL_OUTPUT
							PFResults[nKeptTags*PFResults_ld+(PFResults_ld-1UL)] = (ND[iTag]->revNmer) ? '<' : '>';
#endif

							nKeptTags++;
							continue;

							SKIP_THAT_ONE: ;

#ifdef DEBUG_INFO
							printf("ANA: %u-%u Removing tag %3u %s HasFork=%i revHasFork=%i\n", task->GenomeStart,
										 task->GenomeEnd, iTag, ND[iTag]->Tag, HasFork, revHasFork);
#endif

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
					}
					else {
						nKeptTags = nTags;
						const char MinQuality = (char) options.filterQualValue;
						for (unsigned int iTag=0U; iTag<nTags; iTag++) {
							const size_t TagLength = ND[iTag]->TagLength;
							const char * const restrict CharPtr = (char *) ND[iTag]->Tag;
							const char * restrict QualityPtr = (char *) ND[iTag]->qs;
							const size_t loffset = ND[iTag]->Location - minPos;
							assert(loffset>=0 && loffset<PFResults_ld);
							char * const restrict copyto = &PFResults[iTag*PFResults_ld+loffset];
							assert(TagLength < kMaxReadLen);
#ifndef ANALYZE_FULL_OUTPUT
							for (size_t i=0; i<TagLength; i++) {
								copyto[i] = (QualityPtr[i] >= MinQuality) ? CharPtr[i] : ' ';
							}
#else
							for (size_t i=0; i<TagLength; i++) {
								copyto[i+PFResults_size] = QualityPtr[i];
								copyto[i] = (QualityPtr[i] >= MinQuality) ? CharPtr[i] : ' ';
							}
							PFResults[iTag*PFResults_ld+(PFResults_ld-1UL)] = (ND[iTag]->revNmer) ? '<' : '>';
#endif
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
#ifdef DEBUG_INFO
 					printf("ANA: Kept %u on %u-%u\n", nKeptTags, task->GenomeStart, task->GenomeEnd);
#endif
					if (nKeptTags > 0U) {
						unsigned int WindowStartGenomeAddress;
						const unsigned int InitialWindowStart = WindowStart;
						const unsigned int InitialWindowLength = WindowLength;
						struct ExtraData *nValuableBasesPerTag = malloc(nKeptTags * sizeof(struct ExtraData));

						/* Check Window borders, we cannot have neither insertions nor deletions cut in half */
						{
							/* Keeping genome extension ( insertion) */
							unsigned char * WindowToKeep = (unsigned char *) Buffer;
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
								for (int k=WindowStart-1; k>=(int) dummy; k--) {
									if (GenomeBuffer[k] != ' ') WindowStartGenomeAddress--;
								}
								WindowStart = dummy;
							}
							assert(WindowEnd >= WindowStart);
							WindowLength = WindowEnd - WindowStart + 1;
#ifdef DEBUG_INFO
							printf("WindowStart %u WindowEnd %u WindowLength %u\n", WindowStart, WindowEnd, WindowLength);
#endif
						}

						/* Make a transpose copy of the results */
						const size_t TransposePFResults_ld = ((size_t) nKeptTags + 15UL) & ~(15UL);
						{
							const size_t Size = TransposePFResults_ld*WindowLength;
							if (Size > TransposePFResults_size) {
								if (TransposePFResults) _mm_free(TransposePFResults);
								TransposePFResults_size = Size;
#ifndef ANALYZE_FULL_OUTPUT
								TransposePFResults = (char *) _mm_malloc(Size*sizeof(char),64);
#else
								TransposePFResults = (char *) _mm_malloc(2*Size*sizeof(char),64);
#endif
								if (TransposePFResults == NULL) {
									fprintf(stderr, "Unable to allocate memory for the transpose (%lu), %u tags, window of %lu\n", Size, nTags, WindowLength);
									exit(1);
								}
							}

							for (int k=0; k<nKeptTags; k++) {
								for (size_t l=0; l<WindowLength; l++) {
									TransposePFResults[l*TransposePFResults_ld+k] = PFResults[k*PFResults_ld+l+WindowStart];
#ifdef ANALYZE_FULL_OUTPUT
									TransposePFResults[l*TransposePFResults_ld+k+TransposePFResults_size] = PFResults[k*PFResults_ld+l+WindowStart+PFResults_size];
#endif
								}
							}
						}

#ifdef DEBUG_INFO
						printf("Transposed results kept %u WindowStart %u WindowLength %u\n", nKeptTags, WindowStart, WindowLength);
						for (unsigned int k=0; k<nKeptTags; k++) {
							for (unsigned int l=0; l<WindowStart; l++) fputc(' ', stdout);
							for (unsigned int l=0; l<WindowLength; l++)
								fputc((int) TransposePFResults[l*TransposePFResults_ld+k], stdout);
							fputc('\n', stdout);
						}
#endif

						/* Cleaning the data according to a minimum of occurence : MINIMUM_APPEARANCE */
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
							const unsigned int MinDyn = (options.minorAllelePct*total)/100;
#ifdef DEBUG_INFO
							fprintf(stderr,"total = %u MinDyn = %u\n", total, MinDyn);
#endif
							for (int k=0;k<5;k++) if (Counters[k] < MinDyn) {
								Useless += Counters[k];
								Counters[k] = 0U;
							}

							for (int k=0;k<5;k++) if (Counters[k] < MINIMUM_APPEARANCE) {
								Useless += Counters[k];
								Counters[k] = 0U;
							}

							if ((GenomeBase == ' ') && ((Counters[4] + MINIMUM_APPEARANCE) >= (nKeptTags - Useless))) Counters[4] = 0U;

#if 0
// FIXME - ? I do not understand why the content is modified here
							for (size_t k=0; k<nKeptTags; k++) {
								register unsigned char value = col[k];
								if (value == '-' && GenomeBase != ' ' && Counters[4] == 0U) col[k] = 'x';
								else if (value == GenomeBase && value != ' ') col[k] = '.';
								else if ( value >= 'A' && Counters[(value>>1) & 0b11] == 0U) col[k] = '?';
							}
#endif

							/* Is it valuable or not as a whole ? */
							if (GenomeBase != ' ') {
								HasData[l] = '0';
								if (GenomeBase != 'N') {
									const unsigned char GenomeBaseIndex = (GenomeBase>>1) & 0b11;
									//fprintf(stderr,"l %u %c idx %u Counters %u %u %u %u %u\n", l, GenomeBase, GenomeBaseIndex, Counters[0], Counters[1], Counters[2], Counters[3], Counters[4]);
									for (unsigned char k=0; k<5; k++) {
										if ( k != GenomeBaseIndex && Counters[k]) HasData[l] = '1';
									}
									//fprintf(stderr,"l %u %c idx %u Counters %u %u %u %u %u HasData %c\n",
									//				l, GenomeBase, GenomeBaseIndex, Counters[0], Counters[1], Counters[2], Counters[3], Counters[4], HasData[l]);
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
 							for(unsigned int k=0; k<WindowStart; k++) fputc(' ', stdout);
 							printf("%.*s\n", (int) WindowLength, &GenomeBuffer[WindowStart]);
 							printf("Found %i valuable changes\n", nValuableBases);
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
							const char * const restrict LocalGenomeWindow = &GenomeBuffer[WindowStart];
							while (k<WindowLength) {
								const char c = HasData[k];
								if (previous == '0' && c == '1') {
									/* Correcting the useless insertion to take them into account, not splitting deletions */
									/*int l = (int) k - 1;
									while (l >= 0 && LocalGenomeWindow[l] == ' ' && HasData[l] == '0') l--;
									if (HasData[l] == '1') {
										while(++l < k) { HasData[l] = '1'; nValuableBases++; }
									}
									else*/
										nRegions++;
								}
								k++;
								previous = c;
							}
						}

#if TAG_ANALYSIS_MODE == NIRVANA
						/* Generating VCF...  */
						bcf1_t * const restrict rec = task->VCFRecord;
						assert(rec);
						ngt = bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt_arr);
#ifdef DEBUG_INFO
						printf("ngt = %u, expected 2\n", ngt);
						//assert(ngt == 2);
#endif

						unsigned int CurrentTypeID = 0U;
						if (!bcf_is_snp(rec)) {
							CurrentTypeID = (rec->d.var[bcf_gt_allele(gt[0])].n < 0 || rec->d.var[bcf_gt_allele(gt[1])].n < 0) ? 1 : 2;
						}
						unsigned int RealMemoryPosition = 1 + rec->pos;
						if (CurrentTypeID != 0) RealMemoryPosition += 1; /* God knows why, this is necessary !!! */
#endif

						if (nValuableBases == 0LU) {
#if TAG_ANALYSIS_MODE == NIRVANA
							fprintf(stdout,"\x1b[2K\x1b[GNo variant in interval %u-%u for %u\t%s\t%s\t%s/%s\t%u/%u\n", task->GenomeStart, task->GenomeEnd, RealMemoryPosition, types[CurrentTypeID], rec->d.allele[0], rec->d.allele[bcf_gt_allele(gt[0])], rec->d.allele[bcf_gt_allele(gt[1])], bcf_gt_allele(gt[0]), bcf_gt_allele(gt[1]));
							atomic_fetch_add(&nNoVariant, 1U);
#endif
							goto TaskDone;
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
 							for(unsigned int k=0; k<WindowStart; k++) fputc(' ', stdout);
 							printf("%.*s\n", (int) WindowLength, &GenomeBuffer[WindowStart]);
 							printf("Found %i valuable changes on %u regions\n", nValuableBases, nRegions);
						}
 						printf("Before reg : %.*s\n", (int) WindowLength, HasData);
#endif
						struct RegionsData {
							unsigned int genomepos;
							unsigned int uncovered;
							unsigned short int length;
							unsigned short int windowpos;
						} RegData[/*nRegions*/kMaxReadLen/3];
						assert(nRegions < kMaxReadLen/3);
						char GenomeWindow[WindowLength];
#ifdef DEBUG_INFO
 						printf("After reg : %.*s\n", (int) WindowLength, HasData);
#endif
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
								ptr->uncovered = 0;
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

 						printf("Before concat : %.*s\n", (int) WindowLength, HasData);
#endif

						/* Concatenating tags */
						const size_t Informations_ld = (nValuableBases + 15UL) & ~(15UL);
						assert(nKeptTags*Informations_ld < Buffer_size);
						char * const restrict Informations = Buffer;
						char * restrict InformationsPtr = Informations;
						unsigned int nValuableTags = 0U;

#ifdef ANALYZE_FULL_OUTPUT
						const size_t nRegions4 = (nRegions + 3UL) >> 2;
						const size_t nValuableBases2 = (1+nValuableBases)/2;
						{
							const size_t ultmp = 8UL*(nRegions4 + nValuableBases2)*nKeptTags;
							if (ultmp > RegionsStat_Size) {
								if (RegionsStat) _mm_free(RegionsStat);
								RegionsStat = (unsigned short int *) _mm_malloc(ultmp*sizeof(unsigned short int), 16);
								if (RegionsStat == NULL) {
									fprintf(stderr, "Unable to allocate the region stat array for accounting, size = %lu\n", ultmp*sizeof(unsigned short int));
									exit(1);
								}
								RegionsStat_Size = ultmp;
							}
							memset(RegionsStat, 0, 8UL*(nRegions4)*nKeptTags*sizeof(unsigned short int));
						}
						unsigned short int * restrict PHredStatPtr = &RegionsStat[8UL*nRegions4*nKeptTags];
#endif

						for (size_t k=0; k<nKeptTags; k++) {
							unsigned int pos = 0U;
							unsigned int valuable = 0U;
#ifdef ANALYZE_FULL_OUTPUT
							const __m128i __zero = _mm_setzero_si128();
							for (size_t l=0; l<nValuableBases2; l++) _mm_store_si128((__m128i*) &PHredStatPtr[8UL*l], __zero);
#endif
							for (size_t l=0; l<WindowLength; l++) {
								if (HasData[l] == '1') {
									const char currentBase = TransposePFResults[l*TransposePFResults_ld+k];

#ifdef ANALYZE_FULL_OUTPUT
									size_t StatIndex = 0UL;
									static const unsigned char PHredLimits[3] = { 6, 10, 20 };

									char lcurrentBase;
									unsigned char currentQuality;
									if (currentBase != '-') {
										currentQuality = TransposePFResults[l*TransposePFResults_ld+k+TransposePFResults_size];
										lcurrentBase = currentBase;
									}
									else {
										currentQuality = TransposePFResults[(l-1)*TransposePFResults_ld+k+TransposePFResults_size];
										lcurrentBase = TransposePFResults[(l-1)*TransposePFResults_ld+k];;
									}

									if (lcurrentBase != '-' && lcurrentBase != ' ' && currentQuality >= '\0') {
										while (StatIndex < 3UL) {
											if (currentQuality > PHredLimits[StatIndex])
												StatIndex++;
											else
												break;
										}
										PHredStatPtr[4UL*pos+StatIndex] = 1;
									}
#endif
									InformationsPtr[pos++] = currentBase;
									if ( currentBase != ' ' && currentBase != 'x' && currentBase != '?') valuable++;
								}
							}
#ifdef DEBUG_INFO
 							printf("'%.3s'\t%.*s\n", InformationsPtr, (int) WindowLength, HasData);
#endif
							if (valuable) {
								nValuableBasesPerTag[nValuableTags].ptr = InformationsPtr;
#ifdef ANALYZE_FULL_OUTPUT
								nValuableBasesPerTag[nValuableTags].Accounting = &RegionsStat[(8UL*nRegions4)*nValuableTags];
								nValuableBasesPerTag[nValuableTags].PHred = PHredStatPtr;
								PHredStatPtr += 8UL*nValuableBases2;
								const size_t orientOffset = (PFResults[k*PFResults_ld+PFResults_ld-1UL] == '>') ? 0UL : 1UL;
								assert(PFResults[k*PFResults_ld+PFResults_ld-1UL] == '>' ||  PFResults[k*PFResults_ld+PFResults_ld-1UL] == '<');
								unsigned short int * const restrict Countptr = &RegionsStat[(8UL*nRegions4)*nValuableTags + orientOffset];
								const char * restrict SrcPtr = InformationsPtr;
								for (unsigned int iRegion=0U; iRegion<nRegions; iRegion++) {
									const char * const restrict SrcPtrLimit = SrcPtr + RegData[iRegion].length;
									assert(RegData[iRegion].length > 0);
									do {
										if (*SrcPtr != ' ' && *SrcPtr != 'x' && *SrcPtr != '?') {
											Countptr[2*iRegion] = 1;
											break;
										}
									} while ((uintptr_t) ++SrcPtr < (uintptr_t) SrcPtrLimit);
								}
#endif
								nValuableBasesPerTag[nValuableTags++].n = valuable;
								InformationsPtr += Informations_ld;
							}
						}

						if (nValuableTags == 0U) {
#if TAG_ANALYSIS_MODE == NIRVANA
							fprintf(stdout,"\x1b[2K\x1b[GNo variant in interval %u-%u for %u\t%s\t%s\t%s/%s\t%u/%u\n", task->GenomeStart, task->GenomeEnd, RealMemoryPosition, types[CurrentTypeID], rec->d.allele[0], rec->d.allele[bcf_gt_allele(gt[0])], rec->d.allele[bcf_gt_allele(gt[1])], bcf_gt_allele(gt[0]), bcf_gt_allele(gt[1]));
							atomic_fetch_add(&nNoVariant, 1U);
#endif
							goto TaskDone;
						}
						qsort(nValuableBasesPerTag, nValuableTags, sizeof(struct ExtraData), MyCompare);

#ifdef DEBUG_INFO
						fputs("Initial:\n", stderr);
#ifndef ANALYZE_FULL_OUTPUT
						for (int k=0; k<nValuableTags; k++) {
							fprintf(stderr, "%.*s\t%u/%zu\n", (int) nValuableBases, nValuableBasesPerTag[k].ptr, nValuableBasesPerTag[k].n,nValuableBases);
						}
#else
						for (int k=0; k<nValuableTags; k++) {
							fprintf(stderr, "%.*s\t%u/%zu\t", (int) nValuableBases, nValuableBasesPerTag[k].ptr, nValuableBasesPerTag[k].n,nValuableBases);
							const unsigned short int * restrict Counters = nValuableBasesPerTag[k].Accounting;
							for (unsigned int iRegion=0U; iRegion<nRegions; iRegion++) {
								fprintf(stderr, "%hu,%hu ", Counters[2*iRegion], Counters[2*iRegion+1]);
							}
							Counters = nValuableBasesPerTag[k].PHred;
							for (unsigned int iBase=0U; iBase<nValuableBases; iBase++) {
								fputs(" | ", stderr);
								fprintf(stderr, "%4hu %4hu %4hu %4hu", Counters[4*iBase], Counters[4*iBase+1],Counters[4*iBase+2],Counters[4*iBase+3]);
							}
							fputc('\n', stderr);
						}
#endif
#endif
						unsigned int nSignificative = 1U;
						struct Ordering Significatives[nValuableTags];
						/* Grouping identical ones */
						{
							Significatives[0].data = nValuableBasesPerTag[0];
							Significatives[0].count = 1U;
							for (size_t k=1; k<nValuableTags; k++) {
								int found = 0;
								for (size_t l=0; l<nSignificative; l++) {
									if (strncmp(nValuableBasesPerTag[k].ptr, Significatives[l].data.ptr, (size_t) nValuableBases) == 0) {
										Significatives[l].count++;
										SUM_UP(Significatives[l].data, nValuableBasesPerTag[k]);
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
						fputs("\nGrouping identical\n", stderr);
#ifndef ANALYZE_FULL_OUTPUT
						for (int k=0; k<nSignificative; k++) {
							fprintf(stderr, "%.*s\t%u/%zu\t%u\n",(int) nValuableBases, Significatives[k].data.ptr, Significatives[k].data.n, nValuableBases, Significatives[k].count);
						}
#else
						for (int k=0; k<nSignificative; k++) {
							fprintf(stderr, "%.*s\t%u/%zu\t%u\t",(int) nValuableBases, Significatives[k].data.ptr, Significatives[k].data.n, nValuableBases, Significatives[k].count);
							const unsigned short int * const restrict Counters = Significatives[k].data.Accounting;
							for (unsigned int iRegion=0U; iRegion<nRegions4; iRegion++) {
								fprintf(stderr, "%hu,%hu ", Counters[2*iRegion], Counters[2*iRegion+1]);
							}
							fputc('\n', stderr);
						}
#endif
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
									SUM_UP(Significatives[id].data, Significatives[k].data);
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
							const unsigned int uiMAF = (unsigned int) (MAF*(float) nValuableTags);
							for (unsigned int k=1; k<nSignificative; k++) {
								const char * const Test = Significatives[k].data.ptr;
								if (Test == NULL || Significatives[k].count >= uiMAF) continue;
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
									SUM_UP(Significatives[id].data, Significatives[k].data);
									Significatives[k].data.ptr = NULL;
									Significatives[k].data.n = 0U;
									Significatives[k].count = 0U;
								}
							}

							qsort(Significatives, nSignificative, sizeof(struct Ordering), MyOrdering);
#ifdef DEBUG_INFO
							{
								fprintf(stderr, "\nAssign those below MAF (%5.2f) to groups with up to 1 error\n", MAF);
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
									SUM_UP(Significatives[id].data, Significatives[k].data);
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

						char Info[INFO_SIZE];
						Info[0] = 'W'; Info[1] = 'A'; Info[2] = 'F'; Info[3] = '=';
						char * restrict InfoPtr = &Info[4];
						{
							IF_DEBUG_INFO(fprintf(stderr, "\nKeeping MAF >= %5.2f ...\n\n", MAF));
							assert(nValuableTags > 0);
							const float invnTags = 1.0f/(float) nValuableTags;
							int count = 0;

							float TrashPart = 0.0f;
							for (int k=0; k<nSignificative; k++) {
								if (Significatives[k].data.ptr ) {
									const float val = (float) Significatives[k].count * invnTags;
									if (val >= MAF) {
										IF_DEBUG_INFO(fprintf(stderr, "%.*s\t%u/%zu\t%u\t%f @ %u", (int) nValuableBases, Significatives[k].data.ptr, Significatives[k].data.n,      nValuableBases, Significatives[k].count, val, k));
										if (k < 10) InfoPtr += sprintf(InfoPtr, "%4.2f,", val);
										if (++count > 2) {
											int pos=0;
											int AllCovered = 1;
											for (unsigned int iRegion=0; iRegion<nRegions; iRegion++) {
												int covered = 0;
 												RegData[iRegion].uncovered = 1;
												{
													int isSpace = 1;
													const char * restrict const tptr = Significatives[k].data.ptr + pos;
													for (int l=0; l<RegData[iRegion].length; l++) isSpace &= (tptr[l] == ' ');
													if (isSpace) {
														covered = 1;
														RegData[iRegion].uncovered = 0;
														goto Decide;
													}
												}
												for (int l=0; l<2; l++) {
													if (strncmp(Significatives[k].data.ptr + pos, Significatives[l].data.ptr + pos, RegData[iRegion].length) == 0) {
														RegData[iRegion].uncovered = 0;
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
									else {
										TrashPart += val;
									}
								}
							}
							AvailableAlleles = count;

							assert(InfoPtr < &Info[INFO_SIZE]);
							if (InfoPtr[-1] == ',') InfoPtr--;
							InfoPtr += sprintf(InfoPtr, ";WTRASH=%4.2f;WA=%u;WREG=%u", TrashPart ,AvailableAlleles, nRegions);
						}

#if TAG_ANALYSIS_MODE == NIRVANA
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
								IF_DEBUG_INFO(fprintf(stdout, "--------------------------------------------------------------------------------\n"));
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
									IF_DEBUG_INFO(fprintf(stdout, "Region %i, border correction: pos=%u, length=%u\n", iRegion, RegData[iRegion].genomepos + k, len));
								}
								char * restrict RegInfoPtr = InfoPtr;
								if (RegData[iRegion].uncovered == 1) RegInfoPtr += sprintf(RegInfoPtr,";REGWRN");

								assert( (k+len) <= nValuableBases);
								unsigned char DoTest = 0;
#if TAG_ANALYSIS_MODE == NIRVANA
								if ( (RegData[iRegion].genomepos + (k - nPos)) == RealMemoryPosition) DoTest++;
								if ( (RegData[iRegion].genomepos + (k - nPos)) <= RealMemoryPosition
									&& (RegData[iRegion].genomepos + (k - nPos) + len - 1) >= RealMemoryPosition)
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
									IF_DEBUG_INFO(const char * const toPrint = ptrAllele);
									while(l<len) {
										if (Src[l] == '.') {
											assert (ptrAllele < (PosBackup+MAX_ALLELE_SIZE));
											*ptrAllele++ = Gen[l];
											++nPoints;
										}
										else if (Src[l] != '-' && Src[l] != 'x') {
											assert (ptrAllele < (PosBackup+MAX_ALLELE_SIZE));
											*ptrAllele++ = Src[l];
											if (Src[l] == ' ') ++nDashesOrX;
										}
										else if (Gen[l] == ' ')
// 											++nDashesOrX;
											++nPoints;
										++l;
									}
									if (nPoints == len) {
										//ptrAllele = PosBackup;
										GT[2*Second] = '0';
									}
									else if (nDashesOrX == len) {
										GT[2*Second] = '1';
									}
									else
										GT[2*Second] = GT[0] + 1;

									assert (ptrAllele < (PosBackup+MAX_ALLELE_SIZE));
									*ptrAllele = '\0';

									IF_DEBUG_INFO(printf("ALLELE %i \'%s\', %u (.), %u (-/x/space), GT %c\n", Second, toPrint, nPoints, nDashesOrX, (int) GT[2*Second]));

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
										const char * const restrict ptr = &GenomeBuffer[WindowStart + RegData[iRegion].windowpos + (k-nPos)];
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
// 									if (Heterozygote && (GT[0] != '0' && GT[2] != '0') && *Allele0 == *Allele1) fprintf(stdout, "\x1b[2K\x1b[GOOPS !!!!! \'%s\' and \'%s\'\n", Allele0, Allele1);

									const unsigned int Position =  RegData[iRegion].genomepos + (k - nPos) - InDel;
#if TAG_ANALYSIS_MODE == NIRVANA
									unsigned int ok = 1U;
									int ErrorId = 0;
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

									if (ok) {
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
											if (Heterozygote) {
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
										else {
												fprintf(stderr, "Something went wrong, should not be here!\nPosition is %u, real is %u, VCF is %u\n", Position, RealMemoryPosition, VCFPosition);
												fprintf(stderr, "Found: %u\t%u\t.\t%s\t%s/%s\t%s\nTheory: %u\t%s\t%s\t%s/%s\t%u/%u\n",
															options.chr, RegData[iRegion].genomepos + (k - nPos) - InDel, Previously, Allele0, Allele1, GT,
															RealMemoryPosition, types[CurrentTypeID], rec->d.allele[0], rec->d.allele[bcf_gt_allele(gt[0])], rec->d.allele[bcf_gt_allele(gt[1])], bcf_gt_allele(gt[0]), bcf_gt_allele(gt[1]));
												fprintf(stderr, "Interval was %u-%u\n", task->GenomeStart, task->GenomeEnd);
// 												exit(1);
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
										/*if (Allele0[0] != '\0')*/ {
											fprintf(stdout, "\x1b[2K\x1b[GError %s : %u\t%u\t.\t%s\t%s/%s\t50\t%s:%u is not %u\t%s\t%s\t%s/%s\t%u/%u\n",
															Errors[ErrorId],
															options.chr, RegData[iRegion].genomepos + (k - nPos) - InDel, Previously, Allele0, Allele1, GT, nValuableTags,
															RealMemoryPosition, types[CurrentTypeID], rec->d.allele[0], rec->d.allele[bcf_gt_allele(gt[0])], rec->d.allele[bcf_gt_allele(gt[1])], bcf_gt_allele(gt[0]), bcf_gt_allele(gt[1]));
											fflush(stdout);
										}
									}
									break;
#else
									/* VCF printf("%s\t%u\t.\t%.*s\t%c\t255\tPASS\t.\tGT:DP:AD\t%s:%u:%u,%u\n", */
									IF_DEBUG_INFO(printf("Region start=%u, k=%u, nPos=%u, Indel=%u\n",
									                     RegData[iRegion].genomepos, k, nPos, InDel));

									if ( Position <= task->GenomeEnd && Position >= task->GenomeStart) {
#ifdef ANALYZE_FULL_OUTPUT
										const size_t index = (Genome.virtchr[options.chr].chr << 28) + Genome.virtchr[options.chr].offset;
										const char * restrict GenomePtr = (char *) Genome.table + index + Position;
#endif
										const size_t lenOrig = strlen(Previously);
										const size_t len0 = strlen(Allele0);
										/* ----------------------------- HETEROZYGOTE ----------------------------------------- */
										if (Heterozygote) {
											const size_t len1 = strlen(Allele1);

											/* Verify not being the reference */
											if (len0 == lenOrig && strncmp(Previously, Allele0, len0) == 0) {
												GT[0] = '0';
												if (GT[2] == '2') GT[2] = '1';
											}
											if (len1 == lenOrig && strncmp(Previously, Allele1, len1) == 0) {
												GT[2] = '0';
												if (GT[0] == '2') GT[0] = '1';
											}

#ifdef ANALYZE_FULL_OUTPUT
											const unsigned short int * const restrict PHredStat0 = Significatives[0].data.PHred + nPos*4UL;
											const unsigned short int * const restrict PHredStat1 = Significatives[1].data.PHred + nPos*4UL;
											float sumF0[4] __attribute__((aligned(16)));
											float sumF1[4] __attribute__((aligned(16)));
											const size_t llen0 = (len0>1) ? len0 - InDel : len0;
											__m128i __sum0 = _mm_setzero_si128();
											if (llen0 > 0UL) {
												__sum0 = _mm_loadl_epi64((__m128i*) PHredStat0);
												for (size_t l=1UL; l<len0; l++) {
													__sum0 = _mm_adds_epu16(__sum0, _mm_loadl_epi64((__m128i*) &PHredStat0[4UL*l]));
												};
											}
											const size_t llen1 = (len1>1) ? len1 - InDel : len1;
											__m128i __sum1 = _mm_setzero_si128();
											if (llen1 > 0UL) {
												__sum1 = _mm_loadl_epi64((__m128i*) PHredStat1);
												for (size_t l=1UL; l<len1; l++) {
													__sum1 = _mm_adds_epu16(__sum1, _mm_loadl_epi64((__m128i*) &PHredStat1[4UL*l]));
												};
											}
											const __m128 __invscale0 = _mm_rcp_ps(_mm_set1_ps((float) llen0));
											const __m128 __invscale1 = _mm_rcp_ps(_mm_set1_ps((float) llen1));
											__sum0 = _mm_unpacklo_epi16(__sum0, _mm_setzero_si128());
											__sum1 = _mm_unpacklo_epi16(__sum1, _mm_setzero_si128());
											__m128 __sumF0 = _mm_cvtepi32_ps(__sum0);
											__m128 __sumF1 = _mm_cvtepi32_ps(__sum1);
											__sumF0 = _mm_mul_ps(__sumF0, __invscale0);
											__sumF1 = _mm_mul_ps(__sumF1, __invscale1);
											_mm_store_ps( &sumF0[0], __sumF0);
											_mm_store_ps( &sumF1[0], __sumF1);
#endif

											if (GT[0] != '0' && GT[2] != '0') {
												IF_DEBUG_INFO(fprintf(stdout, "0: \'%s\' (%lu) 1:\'%s\' (%lu) G: \'%s\' (%lu)\n",
																							Allele0, len0, Allele1, len1, Previously, lenOrig));

												/* Based on length analysis */
												if (len0 == len1) {
													if (len0 == lenOrig) {
														for (size_t l=1UL;l<len0; l++) {
#ifdef ANALYZE_FULL_OUTPUT
															char * TmpInfoPtr = RegInfoPtr;
															TmpInfoPtr += sprintf(TmpInfoPtr, ";GEN=\"%10.10s|%1.1s|%10.10s\"", GenomePtr-10UL+l, &Previously[l], GenomePtr + l + 1UL);
															assert(TmpInfoPtr < &Info[INFO_SIZE]);
															const unsigned short int * const restrict uiptr0 = PHredStat0 + (l-1UL)*4UL;
															const unsigned short int * const restrict uiptr1 = PHredStat1 + (l-1UL)*4UL;;
#endif
															if ( (Allele0[l] != Previously[l]) && (Allele1[l] != Previously[l]) ) {
																	if (Allele0[l] != Allele1[l]) {
																		if (Allele0[l] != ' ' && Allele1[l] != ' ') {
#ifdef ANALYZE_FULL_OUTPUT
																			TmpInfoPtr += sprintf(TmpInfoPtr, ";ORI=0,0,%hu,%hu,%hu,%hu",
																														Significatives[0].data.Accounting[2*iRegion], Significatives[0].data.Accounting[2*iRegion+1],
																														Significatives[1].data.Accounting[2*iRegion], Significatives[1].data.Accounting[2*iRegion+1]);
																			TmpInfoPtr += sprintf(TmpInfoPtr, ";PHRED=%hu.0,%hu.0,%hu.0,%hu.0,%hu.0,%hu.0,%hu.0,%hu.0", uiptr0[0], uiptr0[1], uiptr0[2], uiptr0[3], uiptr1[0], uiptr1[1], uiptr1[2], uiptr1[3]);
																			TmpInfoPtr += sprintf(TmpInfoPtr,";GT=1/2;DP=%u;AD=0,%u,%u", nValuableTags,
																							Significatives[0].count, Significatives[1].count);
																			assert(TmpInfoPtr < &Info[INFO_SIZE]);
#endif
																			fprintf(stdout, "%s\t%lu\t.\t%1.1s\t%1.1s,%1.1s\t50\tPASS\t%s;WND=%u-%u\tGT:DP:AD\t1/2:%u:0,%u,%u\n",
																							Genome.virtchr[options.chr].AC, Position+l, &Previously[l], &Allele0[l], &Allele1[l], Info,
																							task->GenomeStart, task->GenomeEnd, nValuableTags,
																							Significatives[0].count, Significatives[1].count);
																		}
																		else if (Allele0[l] == ' ') {
#ifdef ANALYZE_FULL_OUTPUT
																			TmpInfoPtr += sprintf(TmpInfoPtr, ";ORI=0,0,%hu,%hu",
																														Significatives[1].data.Accounting[2*iRegion], Significatives[1].data.Accounting[2*iRegion+1]);
																			const unsigned short int * const restrict uiptr1 = Significatives[1].data.PHred +((l-1UL)+nPos)*4UL;
																			TmpInfoPtr += sprintf(TmpInfoPtr, ";PHRED=%hu.0,%hu.0,%hu.0,%hu.0", uiptr1[0], uiptr1[1], uiptr1[2], uiptr1[3]);
																			TmpInfoPtr += sprintf(TmpInfoPtr,";GT=1/1;DP=%u;AD=0,%u", nValuableTags, Significatives[1].count);
																			assert(TmpInfoPtr < &Info[INFO_SIZE]);
#endif
																			fprintf(stdout, "%s\t%lu\t.\t%1.1s\t%1.1s\t50\tPASS\t%s;WND=%u-%u\tGT:DP:AD\t1/1:%u:0,%u\n",
																							Genome.virtchr[options.chr].AC, Position+l, &Previously[l], &Allele1[l], Info,
																							task->GenomeStart, task->GenomeEnd, nValuableTags, Significatives[1].count);
																		}
																		else {
#ifdef ANALYZE_FULL_OUTPUT
																			TmpInfoPtr += sprintf(TmpInfoPtr, ";ORI=0,0,%hu,%hu",
																														Significatives[0].data.Accounting[2*iRegion], Significatives[0].data.Accounting[2*iRegion+1]);
																			const unsigned short int * const restrict uiptr0 = Significatives[1].data.PHred +((l-1UL)+nPos)*4UL;
																			TmpInfoPtr += sprintf(TmpInfoPtr, ";PHRED=%hu.0,%hu.0,%hu.0,%hu.0", uiptr0[0], uiptr0[1], uiptr0[2], uiptr0[3]);
																			TmpInfoPtr += sprintf(TmpInfoPtr,";GT=1/1;DP=%u;AD=0,%u", nValuableTags, Significatives[0].count);
																			assert(TmpInfoPtr < &Info[INFO_SIZE]);
#endif
																			fprintf(stdout, "%s\t%lu\t.\t%1.1s\t%1.1s\t50\tPASS\t%s;WND=%u-%u\tGT:DP:AD\t1/1:%u:0,%u\n",
																							Genome.virtchr[options.chr].AC, Position+l, &Previously[l], &Allele0[l], Info,
																							task->GenomeStart, task->GenomeEnd, nValuableTags, Significatives[0].count);

																		}
																	}
																	else {
#ifdef ANALYZE_FULL_OUTPUT
																		TmpInfoPtr += sprintf(TmpInfoPtr, ";ORI=0,0,%hu,%hu",
																		                     Significatives[0].data.Accounting[2*iRegion] + Significatives[1].data.Accounting[2*iRegion],
																		                     Significatives[0].data.Accounting[2*iRegion+1] + Significatives[1].data.Accounting[2*iRegion+1]);
																		TmpInfoPtr += sprintf(TmpInfoPtr, ";PHRED=0.0,0.0,0.0,0.0,%hu.0,%hu.0,%hu.0,%hu.0",
																		                      uiptr0[0]+uiptr1[0], uiptr0[1]+ uiptr1[1], uiptr0[2]+uiptr1[2], uiptr0[3]+uiptr1[3]);
																		TmpInfoPtr += sprintf(TmpInfoPtr,";GT=1/1;DP=%u;AD=0,%u", nValuableTags,  Significatives[0].count + Significatives[1].count);
																		assert(TmpInfoPtr < &Info[INFO_SIZE]);
#endif
																		fprintf(stdout, "%s\t%lu\t.\t%1.1s\t%1.1s\t50\tPASS\t%s;WND=%u-%u\tGT:DP:AD\t1/1:%u:0,%u\n",
																						Genome.virtchr[options.chr].AC, Position+l, &Previously[l], &Allele0[l], Info,
																		        task->GenomeStart, task->GenomeEnd, nValuableTags,
																		        Significatives[0].count + Significatives[1].count);
																	}
															}
															else if ( (Allele0[l] != Previously[l]) && (Allele1[l] == Previously[l]) ) {
#ifdef ANALYZE_FULL_OUTPUT
																TmpInfoPtr += sprintf(TmpInfoPtr, ";ORI=%hu,%hu,%hu,%hu",
																											Significatives[1].data.Accounting[2*iRegion], Significatives[1].data.Accounting[2*iRegion+1],
																											Significatives[0].data.Accounting[2*iRegion], Significatives[0].data.Accounting[2*iRegion+1]);
																TmpInfoPtr += sprintf(TmpInfoPtr, ";PHRED=%hu.0,%hu.0,%hu.0,%hu.0,%hu.0,%hu.0,%hu.0,%hu.0",
																	                    uiptr1[0], uiptr1[1], uiptr1[2], uiptr1[3], uiptr0[0], uiptr0[1], uiptr0[2], uiptr0[3]);
																TmpInfoPtr += sprintf(TmpInfoPtr,";GT=0/1;DP=%u;AD=%u,%u", nValuableTags,  Significatives[1].count, Significatives[0].count);
																assert(TmpInfoPtr < &Info[INFO_SIZE]);
#endif
																fprintf(stdout, "%s\t%lu\t.\t%1.1s\t%1.1s\t50\tPASS\t%s;WND=%u-%u\tGT:DP:AD\t0/1:%u:%u,%u\n",
																				Genome.virtchr[options.chr].AC, Position+l, &Previously[l], &Allele0[l], Info,
																				task->GenomeStart, task->GenomeEnd, nValuableTags,
																				Significatives[1].count, Significatives[0].count);
															}
															else if ( (Allele0[l] == Previously[l]) && (Allele1[l] != Previously[l]) ) {
#ifdef ANALYZE_FULL_OUTPUT
																TmpInfoPtr += sprintf(TmpInfoPtr, ";ORI=%hu,%hu,%hu,%hu",
																											Significatives[0].data.Accounting[2*iRegion], Significatives[0].data.Accounting[2*iRegion+1],
																											Significatives[1].data.Accounting[2*iRegion], Significatives[1].data.Accounting[2*iRegion+1]);
																TmpInfoPtr += sprintf(TmpInfoPtr, ";PHRED=%hu.0,%hu.0,%hu.0,%hu.0,%hu.0,%hu.0,%hu.0,%hu.0",
																	                    uiptr0[0], uiptr0[1], uiptr0[2], uiptr0[3], uiptr1[0], uiptr1[1], uiptr1[2], uiptr1[3]);
																TmpInfoPtr += sprintf(TmpInfoPtr,";GT=0/1;DP=%u;AD=%u,%u", nValuableTags,  Significatives[0].count, Significatives[1].count);
																assert(TmpInfoPtr < &Info[INFO_SIZE]);
#endif
																fprintf(stdout, "%s\t%lu\t.\t%1.1s\t%1.1s\t50\tPASS\t%s;WND=%u-%u\tGT:DP:AD\t0/1:%u:%u,%u\n",
																				Genome.virtchr[options.chr].AC, Position+l, &Previously[l], &Allele1[l], Info,
																				task->GenomeStart, task->GenomeEnd, nValuableTags,
																				Significatives[0].count, Significatives[1].count);
															}
														}
													}
													else if (len0 < lenOrig) {
														IF_DEBUG_INFO(fprintf(stdout, "Double deletion !!!\n"));
#ifdef ANALYZE_FULL_OUTPUT
														RegInfoPtr += sprintf(RegInfoPtr, ";GEN=\"%10.10s|%s|%10.10s\"", GenomePtr-10UL, Previously, GenomePtr + strlen(Previously));
														RegInfoPtr += sprintf(RegInfoPtr, ";ORI=0,0,%hu,%hu,%hu,%hu",
																					  			Significatives[0].data.Accounting[2*iRegion], Significatives[0].data.Accounting[2*iRegion+1],
																						  		Significatives[1].data.Accounting[2*iRegion], Significatives[1].data.Accounting[2*iRegion+1]);
														RegInfoPtr += sprintf(RegInfoPtr, ";PHRED=0.0,0.0,0.0,0.0,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f",
																	                    sumF0[0], sumF0[1], sumF0[2], sumF0[3], sumF1[0], sumF1[1], sumF1[2], sumF1[3]);
														RegInfoPtr += sprintf(RegInfoPtr,";GT=%s;DP=%u;AD=0,%u,%u", GT, nValuableTags,  Significatives[0].count, Significatives[1].count);
														assert(RegInfoPtr < &Info[INFO_SIZE]);
#endif
														fprintf(stdout, "%s\t%u\t.\t%s\t%s,%s\t50\tPASS\t%s;WND=%u-%u\tGT:DP:AD\t%s:%u:0,%u,%u\n",
													        Genome.virtchr[options.chr].AC, Position, Previously, Allele0, Allele1, Info,
													        task->GenomeStart, task->GenomeEnd, GT, nValuableTags, Significatives[0].count, Significatives[1].count);
													}
													else {
														IF_DEBUG_INFO(fprintf(stdout, "Double Insertion !!!\n"));
#ifdef ANALYZE_FULL_OUTPUT
														RegInfoPtr += sprintf(RegInfoPtr, ";GEN=\"%10.10s|%s|%10.10s\"", GenomePtr-10UL, Previously, GenomePtr + strlen(Previously));
														RegInfoPtr += sprintf(RegInfoPtr, ";ORI=0,0,%hu,%hu,%hu,%hu",
																							  	Significatives[0].data.Accounting[2*iRegion], Significatives[0].data.Accounting[2*iRegion+1],
																						  		Significatives[1].data.Accounting[2*iRegion], Significatives[1].data.Accounting[2*iRegion+1]);
														RegInfoPtr += sprintf(RegInfoPtr, ";PHRED=0.0,0.0,0.0,0.0,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f",
																	                    sumF0[0], sumF0[1], sumF0[2], sumF0[3], sumF1[0], sumF1[1], sumF1[2], sumF1[3]);
														RegInfoPtr += sprintf(RegInfoPtr,";GT=%s;DP=%u;AD=0,%u,%u", GT, nValuableTags,  Significatives[0].count, Significatives[1].count);
														assert(RegInfoPtr < &Info[INFO_SIZE]);
#endif
														fprintf(stdout, "%s\t%u\t.\t%s\t%s,%s\t50\tPASS\t%s;WND=%u-%u\tGT:DP:AD\t%s:%u:0,%u,%u\n",
													        Genome.virtchr[options.chr].AC, Position, Previously, Allele0, Allele1, Info,
													        task->GenomeStart, task->GenomeEnd, GT, nValuableTags, Significatives[0].count, Significatives[1].count);
													}
												}
												else {
													/* Insertion */
													if ((len0 > lenOrig) && (len1 > lenOrig) ) {
#ifdef ANALYZE_FULL_OUTPUT
														RegInfoPtr += sprintf(RegInfoPtr, ";GEN=\"%10.10s|%s|%10.10s\"", GenomePtr-10UL, Previously, GenomePtr + strlen(Previously));
														RegInfoPtr += sprintf(RegInfoPtr, ";ORI=0,0,%hu,%hu,%hu,%hu",
																						  		Significatives[0].data.Accounting[2*iRegion], Significatives[0].data.Accounting[2*iRegion+1],
																							  	Significatives[1].data.Accounting[2*iRegion], Significatives[1].data.Accounting[2*iRegion+1]);
														RegInfoPtr += sprintf(RegInfoPtr, ";PHRED=0.0,0.0,0.0,0.0,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f",
																	                    sumF0[0], sumF0[1], sumF0[2], sumF0[3], sumF1[0], sumF1[1], sumF1[2], sumF1[3]);
														RegInfoPtr += sprintf(RegInfoPtr,";GT=%s;DP=%u;AD=0,%u,%u", GT, nValuableTags,  Significatives[0].count, Significatives[1].count);
														assert(RegInfoPtr < &Info[INFO_SIZE]);
#endif
														fprintf(stdout, "%s\t%u\t.\t%s\t%s,%s\t50\tPASS\t%s;WND=%u-%u\tGT:DP:AD\t%s:%u:0,%u,%u\n",
																	Genome.virtchr[options.chr].AC, Position, Previously, Allele0, Allele1, Info,
																	task->GenomeStart, task->GenomeEnd, GT, nValuableTags, Significatives[0].count, Significatives[1].count);

													}
													/* Deletions */
													else if ((len0 < lenOrig) && (len1 < lenOrig) ) {
#ifdef ANALYZE_FULL_OUTPUT
														RegInfoPtr += sprintf(RegInfoPtr, ";GEN=\"%10.10s|%s|%10.10s\"", GenomePtr-10UL, Previously, GenomePtr + strlen(Previously));
														RegInfoPtr += sprintf(RegInfoPtr, ";ORI=0,0,%hu,%hu,%hu,%hu",
																						  		Significatives[0].data.Accounting[2*iRegion], Significatives[0].data.Accounting[2*iRegion+1],
																							  	Significatives[1].data.Accounting[2*iRegion], Significatives[1].data.Accounting[2*iRegion+1]);
														RegInfoPtr += sprintf(RegInfoPtr, ";PHRED=0.0,0.0,0.0,0.0,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f",
																	                    sumF0[0], sumF0[1], sumF0[2], sumF0[3], sumF1[0], sumF1[1], sumF1[2], sumF1[3]);
														RegInfoPtr += sprintf(RegInfoPtr,";GT=%s;DP=%u;AD=0,%u,%u", GT, nValuableTags,  Significatives[0].count, Significatives[1].count);
														assert(RegInfoPtr < &Info[INFO_SIZE]);

#endif
														fprintf(stdout, "%s\t%u\t.\t%s\t%s,%s\t50\tPASS\t%s;WND=%u-%u\tGT:DP:AD\t%s:%u:0,%u,%u\n",
																		Genome.virtchr[options.chr].AC, Position, Previously, Allele0, Allele1, Info,
																		task->GenomeStart, task->GenomeEnd, GT, nValuableTags, Significatives[0].count, Significatives[1].count);
													}
													else {
														IF_DEBUG_INFO(fprintf(stdout, "Mixture here!!!\n"));
#ifdef ANALYZE_FULL_OUTPUT
														RegInfoPtr += sprintf(RegInfoPtr, ";ORI=%hu,%hu,%hu,%hu",
																								  Significatives[0].data.Accounting[2*iRegion], Significatives[0].data.Accounting[2*iRegion+1],
																							  	Significatives[1].data.Accounting[2*iRegion], Significatives[1].data.Accounting[2*iRegion+1]);

#endif
														/* Full deletion on allele 0 */
														if (Allele0[1] == '\0' ) {
#ifdef ANALYZE_FULL_OUTPUT
															char * const old = RegInfoPtr;
															RegInfoPtr += sprintf(RegInfoPtr, ";PHRED=0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0;GEN=\"%10.10s|%s|%10.10s\"", GenomePtr-10UL, Previously, GenomePtr + strlen(Previously));
															RegInfoPtr += sprintf(RegInfoPtr,";GT=1|0;DP=%u;AD=%u,%u", nValuableTags,  Significatives[0].count, Significatives[1].count);
															assert(RegInfoPtr < &Info[INFO_SIZE]);
															RegInfoPtr = old;
#endif
															fprintf(stdout, "%s\t%u\t.\t%s\t%s\t50\tPASS\t%s;WND=%u-%u\tGT:DP:AD\t1|0:%u:%u,%u\n",
																			Genome.virtchr[options.chr].AC, Position, Previously, Allele0, Info,
																			task->GenomeStart, task->GenomeEnd, nValuableTags, Significatives[0].count, Significatives[1].count);
															/* Trim left */
															size_t offset = 0UL;
															while(Allele1[offset] == Previously[offset]) offset++;

															/* Trim right
															const int rlen = strlen(Allele1);
															assert(rlen == strlen(Previously));
															for (size_t l=rlen-1; l>offset; l--) {
																if (Allele1[l] == Previously[l]) {
																	Allele1[l] = '\0';
																	Previously[l] = '\0';
																}
																else
																	break;
															}*/
															if (Previously[offset] ==  '\0') offset--;
#ifdef ANALYZE_FULL_OUTPUT
															RegInfoPtr += sprintf(RegInfoPtr, ";PHRED=0.0,0.0,0.0,0.0,%.1f,%.1f,%.1f,%.1f",	sumF1[0], sumF1[1], sumF1[2], sumF1[3]);
															RegInfoPtr += sprintf(RegInfoPtr, ";GEN=\"%10.10s|%s|%10.10s\"", GenomePtr-10UL+offset, &Previously[offset], GenomePtr + strlen(&Previously[offset]));
															RegInfoPtr += sprintf(RegInfoPtr,";GT=0|1;DP=%u;AD=%u,%u", nValuableTags,  Significatives[0].count, Significatives[1].count);
															assert(RegInfoPtr < &Info[INFO_SIZE]);
#endif
															fprintf(stdout, "%s\t%lu\t.\t%s\t%s\t50\tPASS\t%s;WND=%u-%u\tGT:DP:AD\t0|1:%u:%u,%u\n",
																			Genome.virtchr[options.chr].AC, Position+offset, &Previously[offset], &Allele1[offset], Info,
																			task->GenomeStart, task->GenomeEnd, nValuableTags, Significatives[0].count,
																			Significatives[1].count);

														}
														/* Full deletion on allele 1 */
														else if (Allele1[1] == '\0') {
#ifdef ANALYZE_FULL_OUTPUT
															char * const old = RegInfoPtr;
															RegInfoPtr += sprintf(RegInfoPtr, ";PHRED=0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0;GEN=\"%10.10s|%s|%10.10s\"", GenomePtr-10UL, Previously, GenomePtr + strlen(Previously));
															RegInfoPtr += sprintf(RegInfoPtr,";GT=1|0;DP=%u;AD=%u,%u", nValuableTags,  Significatives[0].count, Significatives[1].count);
															assert(RegInfoPtr < &Info[INFO_SIZE]);
															RegInfoPtr = old;
#endif
															fprintf(stdout, "%s\t%u\t.\t%s\t%s\t50\tPASS\t%s;WND=%u-%u\tGT:DP:AD\t1|0:%u:%u,%u\n",
																			Genome.virtchr[options.chr].AC, Position, Previously, Allele1, Info,
																			task->GenomeStart, task->GenomeEnd, nValuableTags, Significatives[0].count, Significatives[1].count);
															/* Trim left */
															int offset = 0;
															while(Allele0[offset] == Previously[offset]) offset++;

															/* Trim right
															const int rlen = (int) strlen(Allele0);
															assert(rlen == (int) strlen(Previously));
															for (int l=rlen-1; l>offset; l--) {
																if (Allele0[l] == Previously[l]) {
																	Allele0[l] = '\0';
																	Previously[l] = '\0';
																}
																else
																	break;
															}*/
															if (Previously[offset] ==  '\0') offset--;
#ifdef ANALYZE_FULL_OUTPUT
															RegInfoPtr += sprintf(RegInfoPtr, ";PHRED=0.0,0.0,0.0,0.0,%.1f,%.1f,%.1f,%.1f",	sumF0[0], sumF0[1], sumF0[2], sumF0[3]);
															RegInfoPtr += sprintf(RegInfoPtr, ";GEN=\"%10.10s|%s|%10.10s\"", GenomePtr-10UL+offset, &Previously[offset], GenomePtr + strlen(&Previously[offset]));
															RegInfoPtr += sprintf(RegInfoPtr,";GT=0|1;DP=%u;AD=%u,%u", nValuableTags,  Significatives[0].count, Significatives[1].count);
															assert(RegInfoPtr < &Info[INFO_SIZE]);
#endif
															fprintf(stdout, "%s\t%u\t.\t%s\t%s\t50\tPASS\t%s;WND=%u-%u\tGT:DP:AD\t0|1:%u:%u,%u\n",
																			Genome.virtchr[options.chr].AC, Position+offset, &Previously[offset], &Allele0[offset], Info,
																			task->GenomeStart, task->GenomeEnd, nValuableTags, Significatives[0].count, Significatives[1].count);


														}
														else {

#ifdef ANALYZE_FULL_OUTPUT
															RegInfoPtr += sprintf(RegInfoPtr, ";PHRED=%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f",
																										sumF0[0], sumF0[1], sumF0[2], sumF0[3], sumF1[0], sumF1[1], sumF1[2], sumF1[3]);

															RegInfoPtr += sprintf(RegInfoPtr, ";GEN=\"%10.10s|%s|%10.10s\"", GenomePtr-10UL, Previously, GenomePtr + strlen(Previously));
															RegInfoPtr += sprintf(RegInfoPtr,";GT=%s;DP=%u;AD=%u,%u", GT, nValuableTags,  Significatives[0].count, Significatives[1].count);
															assert(RegInfoPtr < &Info[INFO_SIZE]);
#endif
															fprintf(stdout, "%s\t%u\t.\t%s\t%s,%s\t50\tPASS\t%s;WND=%u-%u\tGT:DP:AD\t%s:%u:%u,%u\n",
																	Genome.virtchr[options.chr].AC, Position, Previously, Allele0, Allele1, Info,
																	task->GenomeStart, task->GenomeEnd, GT, nValuableTags, Significatives[0].count, Significatives[1].count);

														}
													}
												}
											}
											else {
#ifdef ANALYZE_FULL_OUTPUT
												RegInfoPtr += sprintf(RegInfoPtr, ";GEN=\"%10.10s|%s|%10.10s\"", GenomePtr-10UL, Previously, GenomePtr + strlen(Previously));
												RegInfoPtr += sprintf(RegInfoPtr, ";ORI=%hu,%hu,%hu,%hu",
																						  Significatives[0].data.Accounting[2*iRegion], Significatives[0].data.Accounting[2*iRegion+1],
																						  Significatives[1].data.Accounting[2*iRegion], Significatives[1].data.Accounting[2*iRegion+1]);
												RegInfoPtr += sprintf(RegInfoPtr, ";PHRED=%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f",
																	                    sumF0[0], sumF0[1], sumF0[2], sumF0[3], sumF1[0], sumF1[1], sumF1[2], sumF1[3]);
												RegInfoPtr += sprintf(RegInfoPtr,";GT=%s;DP=%u;AD=%u,%u", GT, nValuableTags,  Significatives[0].count, Significatives[1].count);
												assert(RegInfoPtr < &Info[INFO_SIZE]);
#endif
												const char * const ctmp = (GT[0] != '0') ? Allele0 : Allele1;
												fprintf(stdout, "%s\t%u\t.\t%s\t%s\t50\tPASS\t%s;WND=%u-%u\tGT:DP:AD\t%s:%u:%u,%u\n",
												        Genome.virtchr[options.chr].AC, Position, Previously, ctmp, Info,
												        task->GenomeStart, task->GenomeEnd, GT, nValuableTags,
												        Significatives[0].count, Significatives[1].count);
											}
										}
										/* ----------------------------- HOMOZYGOTE ------------------------------------------ */
										else {
											if (GT[0] != '0') {
												if (len0 != lenOrig) {
#ifdef ANALYZE_FULL_OUTPUT
													RegInfoPtr += sprintf(RegInfoPtr, ";GEN=\"%10.10s|%s|%10.10s\"", GenomePtr-10UL, Previously, GenomePtr + strlen(Previously));
													RegInfoPtr += sprintf(RegInfoPtr, ";ORI=0,0,%hu,%hu", Significatives[0].data.Accounting[2*iRegion], Significatives[0].data.Accounting[2*iRegion+1]);
													const unsigned short int * const restrict uiptr = Significatives[0].data.PHred + nPos*4UL;
													__m128i __sum = _mm_loadl_epi64((__m128i*) uiptr);
													for (size_t l=0UL; l<len0; l++) {
														__sum = _mm_adds_epu16(__sum, _mm_loadl_epi64((__m128i*) &uiptr[4UL*l]));
													};
													const __m128 __invscale = _mm_rcp_ps(_mm_set1_ps((float) len0));
													__sum = _mm_unpacklo_epi16(__sum, _mm_setzero_si128());
													__m128 __sumF = _mm_cvtepi32_ps(__sum);
													__sumF = _mm_mul_ps(__sumF, __invscale);
													float sumF[4] __attribute__((aligned(16)));
													_mm_store_ps( &sumF[0], __sumF);
                          RegInfoPtr += sprintf(RegInfoPtr, ";PHRED=0.0,0.0,0.0,0.0,%.1f,%.1f,%.1f,%.1f", sumF[0], sumF[1], sumF[2], sumF[3]);
													RegInfoPtr += sprintf(RegInfoPtr,";GT=%s;DP=%u;AD=0,%u", GT, nValuableTags,  Significatives[0].count);
													assert(RegInfoPtr < &Info[INFO_SIZE]);
#endif
													fprintf(stdout, "%s\t%u\t.\t%s\t%s\t50\tPASS\t%s;WND=%u-%u\tGT:DP:AD\t%s:%u:0,%u\n",
																	Genome.virtchr[options.chr].AC, Position, Previously, Allele0, Info,
																	task->GenomeStart, task->GenomeEnd, GT, nValuableTags,
																	Significatives[0].count);
												}
												else {
													if (len0 == 1) {
#ifdef ANALYZE_FULL_OUTPUT
                            RegInfoPtr += sprintf(RegInfoPtr, ";GEN=\"%10.10s|%1.1s|%10.10s\"", GenomePtr-10UL, Previously, GenomePtr + 1UL);
														RegInfoPtr += sprintf(RegInfoPtr, ";ORI=%hu,%hu",
                                                  Significatives[0].data.Accounting[2*iRegion], Significatives[0].data.Accounting[2*iRegion+1]);
                            const unsigned short int * const restrict uiptr = Significatives[0].data.PHred +(nPos)*4UL;
                            RegInfoPtr += sprintf(RegInfoPtr, ";PHRED=%hu.0,%hu.0,%hu.0,%hu.0", uiptr[0], uiptr[1], uiptr[2], uiptr[3]);
														RegInfoPtr += sprintf(RegInfoPtr,";GT=1/1;DP=%u;AD=0,%u", nValuableTags,  Significatives[0].count);
                            assert(RegInfoPtr < &Info[INFO_SIZE]);
#endif
                            fprintf(stdout, "%s\t%u\t.\t%1.1s\t%1.1s\t50\tPASS\t%s;WND=%u-%u\tGT:DP:AD\t1/1:%u:0,%u\n",
                                    Genome.virtchr[options.chr].AC, Position, Previously, Allele0, Info,
                                    task->GenomeStart, task->GenomeEnd, nValuableTags,
                                    Significatives[0].count);
													}
													else {
														size_t l = 1UL;
														do {
															if (Allele0[l] != Previously[l]) {
#ifdef ANALYZE_FULL_OUTPUT
																char * TmpInfoPtr = RegInfoPtr;
                            		TmpInfoPtr += sprintf(TmpInfoPtr, ";GEN=\"%10.10s|%1.1s|%10.10s\"", GenomePtr-10UL+l, &Previously[l], GenomePtr + l + 1UL);
																TmpInfoPtr += sprintf(TmpInfoPtr, ";ORI=%hu,%hu",
																									    Significatives[0].data.Accounting[2*iRegion], Significatives[0].data.Accounting[2*iRegion+1]);
																const unsigned short int * const restrict uiptr = Significatives[0].data.PHred +((l-1UL)+nPos)*4UL;
																TmpInfoPtr += sprintf(TmpInfoPtr, ";PHRED=%hu.0,%hu.0,%hu.0,%hu.0", uiptr[0], uiptr[1], uiptr[2], uiptr[3]);
																TmpInfoPtr += sprintf(TmpInfoPtr,";GT=1/1;DP=%u;AD=0,%u", nValuableTags,  Significatives[0].count);
																assert(TmpInfoPtr < &Info[INFO_SIZE]);
#endif
																fprintf(stdout, "%s\t%lu\t.\t%1.1s\t%1.1s\t50\tPASS\t%s;WND=%u-%u\tGT:DP:AD\t1/1:%u:0,%u\n",
																				Genome.virtchr[options.chr].AC, Position+l, &Previously[l], &Allele0[l], Info,
																				task->GenomeStart, task->GenomeEnd, nValuableTags,
																				Significatives[0].count);
															}
														} while (++l<len0);
													}
												}
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
#ifdef DEBUG_INFO
						printf("ngt = %u, expected 2\n", ngt);
						//assert(ngt == 2);
#endif

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

				/* Free extra allocated memory in specific cases */
				if (nTags > 4000) {
					free(Buffer);
					Buffer_size = 0UL;
					Buffer = NULL;
					free(score);
					score_size = 0UL;
					score = NULL;
					_mm_free(PFResults);
					PFResults_size = 0UL;
					PFResults = NULL;
					_mm_free(TransposePFResults);
					TransposePFResults = NULL;
					TransposePFResults_size = 0UL;
#ifdef ANALYZE_FULL_OUTPUT
					_mm_free(RegionsStat);
					RegionsStat = NULL;
					RegionsStat_Size = 0UL;
#endif
				}
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

	if (GenomeBuffer) free(GenomeBuffer);
	if (PFResults) _mm_free(PFResults);
	if (TransposePFResults) _mm_free(TransposePFResults);
	if (score) free(score);
	if (Buffer) free(Buffer);
	if (matrix) _mm_free(matrix);
	if (WORK) _mm_free(WORK);
	if (CellList) free(CellList);
	if (Deletions) free(Deletions);
#ifdef ANALYZE_FULL_OUTPUT
	if (RegionsStat) _mm_free(RegionsStat);
#endif
	return (void*) 0;
}
#undef SUM_UP
#undef INFO_SIZE

/* vim: set filetype=c tabstop=2 softtabstop=2 shiftwidth=2 noexpandtab : */
