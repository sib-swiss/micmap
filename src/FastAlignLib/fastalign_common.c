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
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <assert.h> 
#include "json.h"
#include "fastalign.h"


/* WARNING: DUE TO MASK IN MIC, WE HAVE INVERSED MATCH MISMATCH BIT VALUE */
#ifdef __MIC__ 
#define IS_REAL_MATCH 0
#define IS_MISMATCH   1
#else
#define IS_REAL_MATCH 1
#define IS_MISMATCH   0
#endif
/***************************************************************************/

struct Profile prf = {
	/* Enlarge genome range */
	.ExtraLeftGenomeShift = 25U,
	/* Starting range for allowed jump in */
	.AllowedJumpIn = 0U,
	/* Ending range for jump out */
	.AllowedJumpOut = 0U,
	/* Match */
	.M = ( _M << SCORE_SHIFT) + ( PRIORITY_MATCH << STATE_SHIFT)     + (IS_REAL_MATCH << MISMATCH_SHIFT),
	/* Mismatch */
	.m = ( _m << SCORE_SHIFT) + ( PRIORITY_MATCH << STATE_SHIFT)     + (IS_MISMATCH << MISMATCH_SHIFT),
	.mB = ( _m << SCORE_SHIFT) + ( PRIORITY_MATCH << STATE_SHIFT)     + (IS_MISMATCH << MISMATCH_SHIFT),
	.mE = ( _m << SCORE_SHIFT) + ( PRIORITY_MATCH << STATE_SHIFT)     + (IS_MISMATCH << MISMATCH_SHIFT),
	/* Deletion */
	.D = ( _D << SCORE_SHIFT) + ( PRIORITY_DELETION << STATE_SHIFT)  + (IS_MISMATCH << MISMATCH_SHIFT),
	/* Insertion */
	.I = ( _I << SCORE_SHIFT) + ( PRIORITY_INSERTION << STATE_SHIFT) + (IS_MISMATCH << MISMATCH_SHIFT),
	/* Transitions */
	/*                  D                 X                 M                    I (see vector position) */
	.Match.To     = {  _MD << SCORE_SHIFT, 0 << SCORE_SHIFT, _MM << SCORE_SHIFT,  _MI << SCORE_SHIFT},
	.Insertion.To = { NLOW << SCORE_SHIFT, 0 << SCORE_SHIFT, _IM << SCORE_SHIFT,  _II << SCORE_SHIFT},
	.Deletion.To  = {  _DD << SCORE_SHIFT, 0 << SCORE_SHIFT, _DM << SCORE_SHIFT, NLOW << SCORE_SHIFT}
};

void dumpAlignerScores(const char * FileName)
{
	FILE* const restrict out = fopen(FileName, "w");
	if (out) {
		fprintf(out,"{\n" 
		             "/* Genome Extra space in front for search */\n"
		             "\t\"GenomeExtraFrontSpace\": %u,\n"
		             "\t\"AllowedJumpIn\": %u,\n"
		             "\t\"AllowedJumpOut\": %u,\n"
		             "/* --------- Aligner scores ---------- */\n"
		             "\t\"Scores\": {\n"
		             "\t/* State scores */\n"
		             "\t\t\"Match\": %i,\n"
		             "\t\t\"MisMatch\": %i,\n"
		             "\t\t\"MisMatchAtStart\": %i,\n"
		             "\t\t\"MisMatchAtEnd\": %i,\n"
		             "\t\t\"Insertion\": %i,\n"
		             "\t\t\"Deletion\": %i\n"
		             "\t},\n"
		             "/* State transition scores, note that minimal value (NLOW is %i) */\n"
		             "\t\"Transitions\": {\n"
		             "\t\t\"MM\": %i, \"MI\": %i, \"MD\": %i,\n"
		             "\t\t\"IM\": %i, \"II\": %i, \"ID\": %i,\n"
		             "\t\t\"DM\": %i, \"DI\": %i, \"DD\": %i\n"
		             "\t}\n"
		             "}\n",
		             prf.ExtraLeftGenomeShift, prf.AllowedJumpIn, prf.AllowedJumpOut,
		             ToScore(prf.M),
		             ToScore(prf.m), ToScore(prf.mB), ToScore(prf.mE),
		             ToScore(prf.I), ToScore(prf.D),  NLOW, 
		             ToScore(prf.Match.To[MATCH]),  ToScore(prf.Match.To[INSERTION]),  ToScore(prf.Match.To[DELETION]),
		             ToScore(prf.Insertion.To[MATCH]),  ToScore(prf.Insertion.To[INSERTION]),  ToScore(prf.Insertion.To[DELETION]),
		             ToScore(prf.Deletion.To[MATCH]),  ToScore(prf.Deletion.To[INSERTION]),  ToScore(prf.Deletion.To[DELETION]));
		fclose(out);
	}
	else {
		fprintf(stderr, "Unable to generate JSON file %s\n", FileName);
	}
}
//---------------------------------------------------------------

int loadAlignerScores(const char * FileName)
{
	int err = -1;
	
	json_object * const jBase = json_object_from_file(FileName);
	if (! jBase) {
		fprintf(stderr,"Error loading JSON file %s\n", FileName);
		goto bail;
	}
	
	json_object * jElement;

	json_object_object_get_ex(jBase, "GenomeExtraFrontSpace", &jElement);
	if (jElement == NULL) {
		fprintf(stderr, "Unable to find GenomeExtraFrontSpace within JSON file %s\n", FileName);
		goto bail;
	}
	else {
		if (json_object_is_type(jElement, json_type_int)) {
			const int value = json_object_get_int(jElement);
			if (value >= 0) prf.ExtraLeftGenomeShift = (unsigned int) value; 
		}
		else {
			fputs("Unable to read correct format for GenomeExtraFrontSpace\n", stderr);
		}
	}
	
	json_object_object_get_ex(jBase, "AllowedJumpIn", &jElement);
	if (jElement == NULL) {
		fprintf(stderr, "Unable to find AllowedJumpIn within JSON file %s\n", FileName);
		goto bail;
	}
	else {
		if (json_object_is_type(jElement, json_type_int)) {
			const int value = json_object_get_int(jElement);
			if (value >= 0) prf.AllowedJumpIn = (unsigned int) value; 
		}
		else {
			fputs("Unable to read correct format for AllowedJumpIn\n", stderr);
		}
	}
	
	json_object_object_get_ex(jBase, "AllowedJumpOut", &jElement);
	if (jElement == NULL) {
		fprintf(stderr, "Unable to find AllowedJumpOut within JSON file %s\n", FileName);
		goto bail;
	}
	else {
		if (json_object_is_type(jElement, json_type_int)) {
			const int value = json_object_get_int(jElement);
			if (value >= 0) prf.AllowedJumpOut = (unsigned int) value; 
		}
		else {
			fputs("Unable to read correct format for AllowedJumpOut\n", stderr);
		}
	}
	
#define LOAD(string, where) {\
	json_object_object_get_ex(jElement, string, &jSubElement);\
	if (jSubElement) {\
		if (json_object_is_type(jSubElement, json_type_int)) {\
			where = json_object_get_int(jSubElement);\
		}\
		else {\
			fputs("Unable to read correct format for " string "\n", stderr);\
			goto bail;\
		}\
	}\
	else {\
		fprintf(stderr, "Unable to find " string " within JSON file %s\n", FileName);\
		goto bail;\
	}\
}
	
	json_object_object_get_ex(jBase, "Scores", &jElement);
	if (jElement == NULL) {
		fprintf(stderr, "Unable to find Scores within JSON file %s\n", FileName);
		goto bail;
	}
	else {
		json_object * jSubElement;
		LOAD("Match", prf.M);
		LOAD("MisMatch", prf.m);
		LOAD("MisMatchAtStart", prf.mB);
		LOAD("MisMatchAtEnd", prf.mE);
		LOAD("Insertion", prf.I);
		LOAD("Deletion", prf.D);
	}
	
	json_object_object_get_ex(jBase, "Transitions", &jElement);
	if (jElement == NULL) {
		fprintf(stderr, "Unable to find Transitions within JSON file %s\n", FileName);
		goto bail;
	}
	else {
		json_object * jSubElement;
		LOAD("MM", prf.Match.To[MATCH]);LOAD("MI", prf.Match.To[INSERTION]);LOAD("MD", prf.Match.To[DELETION]);
		LOAD("IM", prf.Insertion.To[MATCH]);LOAD("II", prf.Insertion.To[INSERTION]);LOAD("ID", prf.Insertion.To[DELETION]);
		LOAD("DM", prf.Deletion.To[MATCH]);LOAD("DI", prf.Deletion.To[INSERTION]);LOAD("DD", prf.Deletion.To[DELETION]);
	}
#undef LOAD

	prf.M = ( prf.M << SCORE_SHIFT) + ( PRIORITY_MATCH << STATE_SHIFT)     + (IS_REAL_MATCH << MISMATCH_SHIFT);
	/* Mismatch */
	prf.m  = ( prf.m  << SCORE_SHIFT) + ( PRIORITY_MATCH << STATE_SHIFT)   + (IS_MISMATCH << MISMATCH_SHIFT);
	prf.mB = ( prf.mB << SCORE_SHIFT) + ( PRIORITY_MATCH << STATE_SHIFT)   + (IS_MISMATCH << MISMATCH_SHIFT);
	prf.mE = ( prf.mE << SCORE_SHIFT) + ( PRIORITY_MATCH << STATE_SHIFT)   + (IS_MISMATCH << MISMATCH_SHIFT);
	/* Deletion */
	prf.D = ( prf.D << SCORE_SHIFT) + ( PRIORITY_DELETION << STATE_SHIFT)  + (IS_MISMATCH << MISMATCH_SHIFT);
	/* Insertion */
	prf.I = ( prf.I << SCORE_SHIFT) + ( PRIORITY_INSERTION << STATE_SHIFT) + (IS_MISMATCH << MISMATCH_SHIFT);

#ifdef __MIC__
	for (int i=0; i<4;i++) prf.Match.To[i]     <<= SCORE_SHIFT;
	for (int i=0; i<4;i++) prf.Insertion.To[i] <<= SCORE_SHIFT;
	for (int i=0; i<4;i++) prf.Deletion.To[i]  <<= SCORE_SHIFT;
	
#else
	prf.Match.Vector     = _mm_slli_epi32(prf.Match.Vector, SCORE_SHIFT);
	prf.Insertion.Vector = _mm_slli_epi32(prf.Insertion.Vector, SCORE_SHIFT);
	prf.Deletion.Vector  = _mm_slli_epi32(prf.Deletion.Vector, SCORE_SHIFT);
#endif
	err = 0;
	
	bail:;
	return err;
}
//---------------------------------------------------------------

void printCigar(const GlobalData_t * const restrict data) {
	const unsigned char * restrict ptr = data->Cigar;
	const unsigned char * const restrict Limit = data->Cigar + (2*kMaxEncodedCigar);
	while (*ptr != kCIGARTERMINATOR) {
		if (ptr >= Limit) { printf("END ERROR"); break; }
		printf("%i%c ", (int) ptr[0], ptr[1]); ptr += 2;
	}
	fputc('\n', stdout);
}
//---------------------------------------------------------------

void printMismatch(const GlobalData_t * const restrict data) {
	const unsigned char * restrict ptr = data->Mismatch;
	const unsigned char * const restrict Limit = data->Mismatch + (2*kMaxEncodedMismatch);
	while (*ptr != kMMTERMINATOR) {
		if (ptr >= Limit) { printf("END ERROR"); break; }
		printf("%i%c ", (int) ptr[0], ptr[1]); ptr += 2;
	}
	fputc('\n', stdout);
}
//---------------------------------------------------------------

void printCigarStr(const unsigned char  * const restrict data, FILE* out) {
	const unsigned char * restrict ptr = data;
	const unsigned char * const restrict Limit = data + (2*kMaxEncodedCigar);
	while (*ptr != kCIGARTERMINATOR) {
		if (ptr >= Limit) { fprintf(out, "END ERROR"); break; }
		fprintf(out, "%i%c ", (int) ptr[0], ptr[1]); ptr += 2;
	}
	fputc('\n', out);
}
//---------------------------------------------------------------

void printMismatchStr(const unsigned char * const restrict data, FILE* out) {
	const unsigned char * restrict ptr = data;
	const unsigned char * const restrict Limit = data + (2*kMaxEncodedMismatch);
	while (*ptr != kMMSOFTCLIP && *ptr != kMMTERMINATOR) {
		if (ptr >= Limit) { fprintf(out,"END ERROR"); break; }
		fprintf(out,"%i%c ", (int) ptr[0], ptr[1]); ptr += 2;
	}
	fputc('\n', out);
}
//---------------------------------------------------------------

int EncodeInternalAlignment(GlobalData_t * const restrict data, const unsigned int SoftClipMismatchToRemove )
{
	unsigned char * restrict ptrState = data->States;
	unsigned char * restrict Mismatch = data->Mismatch;
	const unsigned char * const restrict MismatchLimit = Mismatch + (2*kMaxEncodedMismatch);
	unsigned char * restrict Cigar = data->Cigar;
	const unsigned char * restrict const CigarLimit = Cigar + (2*kMaxEncodedCigar);
	const unsigned char * const restrict Tag = &data->Genome[data->AlignmentRange[0]];
	unsigned int id = 0U;
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// WITH SOFTCLIP
	int RemovedMismatch = 0;
	int RemovedCigar = 0;
	int RemoveOther = 0;
	int MismatchToRemove = 0;
	int CigarToRemove = 0;
	
	/* Do it right
	 *                                   + LastDumpCigarLimit ptr = CigarLimit - 4
	 *                                   |             + this is CigarLimit ptr, index 2*kMaxEncodedCigar
	 *                                   |             |
	 * ... | 35 M |  3 I | 10 M |  1 D | 23 M | 45 S | TERMINATOR 
	 * 
	 */
	if ( SoftClipMismatchToRemove > 0 ) {
		const unsigned int TotalRequiredCigCount = data->RequiredCigSize;
		const unsigned int TotalRequiredMMCount = data->RequiredMMSize;
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		// TAG is reversed
		if (data->revNmer) {
			unsigned int pos = 0U;
			int chr_pos_shift = 0;
			char  previousState = *ptrState;
			if (previousState != 'I'&& previousState != 'D') previousState = 'M';
			
			while (RemovedMismatch < SoftClipMismatchToRemove) {
				switch(*ptrState) {
					case 'm':
						RemovedMismatch++;
						++pos;
						break;
					case 'I':
						chr_pos_shift--;
					case 'n':
					case 'N':
						RemoveOther++;
					case 'M':
						++pos;
						break;
					case 'D':
						chr_pos_shift++;
						break;
					default:
						fprintf(stderr, "Unknown state '%i' found @ %i softclip reverse (%u/%u) score %i\n%s\n",
										(int) *ptrState, pos, RemovedMismatch, SoftClipMismatchToRemove, data->score,
										data->States);
						fflush(stderr);
						goto bail;
						break;
				}
				char currentState = *ptrState;
				if (currentState != 'I'&& currentState != 'D') currentState = 'M';
				if (currentState != previousState) {
					RemovedCigar++;
					previousState = currentState;
				}
				ptrState++;
			}
			
			{
				MismatchToRemove = TotalRequiredMMCount - RemovedMismatch - kMaxEncodedMismatch;
				while (RemoveOther < MismatchToRemove) {
					switch(*ptrState) {
						case 'm':
							++pos;
							break;
						case 'I':
							chr_pos_shift--;
						case 'n':
						case 'N':
							RemoveOther++;
						case 'M':
							++pos;
							break;
						case 'D':
							chr_pos_shift++;
							break;
						default:
							fprintf(stderr, "Unknown state '%i' found @ %i softclip reverse (%u/%u) score %i\n%s\n",
											(int) *ptrState, pos, RemoveOther, MismatchToRemove, data->score, data->States);
							fflush(stderr);
							goto bail;
							break;
					}
					char currentState = *ptrState;
					if (currentState != 'I'&& currentState != 'D') currentState = 'M';
					if (currentState != previousState) {
						RemovedCigar++;
						previousState = currentState;
					}
					ptrState++;
				}
			}
			
			{
				CigarToRemove = TotalRequiredCigCount - kMaxEncodedCigar;
				while (RemovedCigar <= CigarToRemove) {
					switch(*ptrState) {
						case 'm':
							++pos;
							break;
						case 'I':
							chr_pos_shift--;
						case 'n':
						case 'N':
						case 'M':
							++pos;
							break;
						case 'D':
							chr_pos_shift++;
							break;
						default:
							fprintf(stderr, "Unknown state '%i' found @ %i softclip reverse (%u/%u) score %i\n%s\n",
											(int) *ptrState, pos, RemovedCigar, CigarToRemove, data->score, data->States);
							fflush(stderr);
							goto bail;
							break;
					}
					char currentState = *ptrState;
					if (currentState != 'I' && currentState != 'D') currentState = 'M';
					if (currentState != previousState) {
						RemovedCigar++;
						previousState = currentState;
					}
					ptrState++;
				}
			}
			
			// Adjust chromosome position
			data->AlignmentRange[0] += chr_pos_shift;
			
			// Add softclip region  
#if kMAXcigarSize <= 2
#error "kMAXcigarSize should be at least 3"
#endif
			Cigar[0] = (unsigned char) pos;
			Cigar[1] = 'S';
			Cigar += 2;
			
			previousState = *ptrState;
			if (previousState != 'I'&& previousState != 'D') previousState = 'M';
			unsigned int StateCount = 0U;
			
			while (*ptrState != '\0') {
				char currentState;
				switch(*ptrState) {
					case 'I':
						if (Mismatch >= MismatchLimit) goto bail;
						Mismatch[0] = pos;
						Mismatch[1] = Tag[pos];
						Mismatch += 2;
						currentState = 'I';
						++pos;
						break;
					case 'm':
					case 'n':
					case 'N':
						if (Mismatch >= MismatchLimit) goto bail;
						Mismatch[0] = pos;
						Mismatch[1] = Tag[pos];
						Mismatch += 2;
					case 'M':
						currentState = 'M';
						++pos;
						break;
					case 'D':
						currentState = 'D';
						break;
					default:
						fprintf(stderr, "Unknown state '%i' found @ %i softclip reverse score %i\n%s\nCIG: ", 
										(int) *ptrState, pos, data->score, data->States);
						unsigned char * ptr = data->Cigar;
						while (ptr < Cigar) { fprintf(stderr, "%i%c ", (int) ptr[0], ptr[1]); ptr += 2; }
						fprintf(stderr, "\nMM : ");
						ptr = data->Mismatch;
						while (ptr < Mismatch) { fprintf(stderr, "%i%c ", (int) ptr[0], ptr[1]); ptr += 2; }
						fprintf(stderr, "\n");
						fflush(stderr);
						goto bail;
						break;
				}
				
				if (previousState != currentState) {
					if (Cigar >= CigarLimit) goto bail;
					assert(StateCount < 256);
					Cigar[0] = ((unsigned char) StateCount);
					Cigar[1] = previousState;
					Cigar += 2;
					previousState = currentState;
					StateCount = 1U;
				}
				else {
					StateCount++;
				}
				ptrState++;
			}
			// terminate 
			*Mismatch = kMMTERMINATOR;

			if (StateCount) {
				if (Cigar >= CigarLimit) goto bail;
				assert(StateCount < 256);
				Cigar[0] = ((unsigned char) StateCount);
				Cigar[1] = previousState;
				Cigar += 2;
			}
			
			if (Cigar > CigarLimit) goto bail;
			*Cigar = kCIGARTERMINATOR;
			
		}
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		// TAG is correctly oriented
		else {
			unsigned int pos = 0U;
			unsigned int KeptMismatch = 0U;
			char previousState = *ptrState;
			if (previousState != 'I'&& previousState != 'D') previousState = 'M';
			unsigned int StateCount = 0U;
			
			const unsigned char * const restrict LastDumpCigarLimit = CigarLimit - 2 - 2;
			while (KeptMismatch <= kMaxMismatch) {
				char currentState = *ptrState;
				if (currentState != 'I'&& currentState != 'D') currentState = 'M';
				if (previousState != currentState) {
					if (Cigar >= LastDumpCigarLimit) {
						if (pos >= data->SoftClipBoundary) {
							goto CheckFurtherMM;
						}
						else {
#ifndef NDEBUG
							printf("Within loop exit at cigar location %lu, pos %u, softclip @ %u\n",
										 (uintptr_t) Cigar - (uintptr_t) data->Cigar,
										 pos, data->SoftClipBoundary);
#endif
							goto bail;;
						}
					}
					assert(StateCount < 256);
					Cigar[0] = ((unsigned char) StateCount);
					Cigar[1] = previousState;
					Cigar += 2;
					previousState = currentState;
					StateCount = 1U;
				}
				else {
					StateCount++;
				}
				
				switch(*ptrState) {
					case 'I':
						if (Mismatch >= MismatchLimit) {
							if (pos >= data->SoftClipBoundary) {
								data->MismatchCount = KeptMismatch;
								/* Correct cigar since that base has been accounted */
								StateCount--;
								goto CheckFurtherMM;
							}
							else {
								goto bail;
							}
						}
						Mismatch[0] = pos;
						Mismatch[1] = Tag[pos];
						Mismatch += 2;
						pos++;
						break;
					case 'm':
						if (++KeptMismatch > kMaxMismatch) {
							/* Correct cigar if necessary */
							StateCount--;
							goto CheckFurtherMM;
						}
					case 'n':
					case 'N':
						if (Mismatch >= MismatchLimit) {
							if (pos >= data->SoftClipBoundary) {
								/* Correct cigar since that base has been accounted */
								StateCount--;
								goto CheckFurtherMM;
							}
							else {
								goto bail;
							}
						}
						Mismatch[0] = pos;
						Mismatch[1] = Tag[pos];
						Mismatch += 2;
					case 'M':
						pos++;
						break;
					case 'D':
						break;
					default:
						{
							fprintf(stderr, "Unknown state '%i' found @ %i softclip (%u/%u) score %i\n%s\nCIG: ",
											(int) *ptrState, pos, KeptMismatch, kMaxMismatch, data->score, data->States);
							unsigned char * ptr = data->Cigar;
							while (ptr < Cigar) { fprintf(stderr, "%i%c ", (int) ptr[0], ptr[1]); ptr += 2; }
							fprintf(stderr, "\nMM : ");
							ptr = data->Mismatch;
							while (ptr < Mismatch) { fprintf(stderr, "%i%c ", (int) ptr[0], ptr[1]); ptr += 2; }
							fprintf(stderr, "\n");
							fflush(stderr);
							goto bail;
						}
						break;
				}
				ptrState++;
			}
			// terminate 
		CheckFurtherMM:;
			if (StateCount > 0U) {
				if (previousState != 'D') {
					if (Cigar >= CigarLimit) goto bail;
					assert(StateCount < 256);
					Cigar[0] = ((unsigned char) StateCount);
					Cigar[1] = 'M';
					Cigar += 2;
				}
			}
			*Mismatch = kMMTERMINATOR;
			
			// Add softclip region 
			if (Cigar > CigarLimit) goto bail;
			const int taglen = (int) data->TagLength;
			assert(taglen > pos);
			Cigar[0] = (unsigned char) (taglen - pos); // potentially +1
			Cigar[1] = 'S';
			Cigar += 2;
		  
			if (Cigar > CigarLimit) goto bail;
			*Cigar = kCIGARTERMINATOR;
		}
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// WITHOUT SOFTCLIP
	else {
		char previousState = *ptrState;
		if (previousState != 'I'&& previousState != 'D') previousState = 'M';
		unsigned int StateCount = 0U;
		unsigned char pos = 0;
		while (*ptrState != '\0') {
			char currentState;
			switch(*ptrState) {
				case 'I':
					if (Mismatch >= MismatchLimit) goto bail;
					Mismatch[0] = pos;
					Mismatch[1] = Tag[pos];
					Mismatch += 2;
					currentState = 'I';
					++pos;
					break;
				case 'N':
				case 'n':
				case 'm':
					if (Mismatch >= MismatchLimit) goto bail;
					Mismatch[0] = pos;
					Mismatch[1] = Tag[pos];
					Mismatch += 2;
				case 'M':
					currentState = 'M';
					++pos;
					break;
				case 'D':
					currentState = 'D';
					break;
				default:
					fprintf(stderr, "Unknown state '%c' found @ %i no softclip score %i\nSTA: %.*s\n",
									*ptrState, pos, data->score, (int) pos, data->States);
					fprintf(stderr, "State string starting at 0x%16.16zx\n"
					                "State Memory starting at 0x%16.16zx up to 0x%16.16zx\n"
													"Current State pointer    0x%16.16zx\nCIG: ",
													(uintptr_t) data->States, (uintptr_t) data->StateSequence,
													(uintptr_t) &(data->StateSequence[STATE_MEMORY_SIZE-1]), (uintptr_t) ptrState);
					unsigned char * ptr = data->Cigar;
					while (ptr < Cigar) { fprintf(stderr, "%i%c ", (int) ptr[0], ptr[1]); ptr += 2; }
					fprintf(stderr, "\nMM : ");
					ptr = data->Mismatch;
					while (ptr < Mismatch) { fprintf(stderr, "%i%c ", (int) ptr[0], ptr[1]); ptr += 2; }
					fprintf(stderr, "\n");
					fflush(stderr);
					goto bail;
					break;
			}
			

			if (previousState != currentState) {
				if (Cigar >= CigarLimit) goto bail;
				assert(StateCount < 256);
				Cigar[0] = ((unsigned char) StateCount);
				Cigar[1] = previousState;
				Cigar += 2;
				previousState = currentState;
				StateCount = 1U;
			}
			else {
				StateCount++;
			}
			ptrState++;
		}
		
		*Mismatch = kMMTERMINATOR;

		if (StateCount) {
			if (Cigar >= CigarLimit) goto bail;
			assert(StateCount < 256);
			Cigar[0] = ((unsigned char) StateCount);
			Cigar[1] = previousState;
			Cigar += 2;
		}
		
		if (Cigar > CigarLimit) goto bail;
		*Cigar = kCIGARTERMINATOR;
	}
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// COPY ALIGNMENT INFO TO TRANSFER BUFFER
	
	const unsigned int ciglen = 1 + (unsigned int) ((uintptr_t) Cigar - (uintptr_t) &(data->Cigar[0]));
	if (ciglen == 3) data->Cigar[0]    = kCIGARTERMINATOR;
	
	return 0;
	
	bail:;
		printf("Unable to encode %u mismatch with score %i\n", data->MismatchCount, data->score);
		return -1;
}
//---------------------------------------------------------------

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
