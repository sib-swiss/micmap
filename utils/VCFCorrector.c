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
#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>

char * Buffer = NULL;

static void CorrectAD(char * const Value)
{
		if (Value[0] == '1' && Value[1] == '/' && Value[2] == '2') {
				fprintf(stderr, "FOUND 1/2: %s%s", Value, Buffer);
				char * restrict ptr = Value;
				while(*ptr != ':') ptr++;
				++ptr;
				while(*ptr != ':') ptr++;
				++ptr;
				char * const start = ptr;
				int count = 0;
				while(*ptr != '\n') if (*ptr++ == ',') count++; 
				if (count == 1) {
					ptr[3] = '\0';
					while ((uintptr_t)ptr >= (uintptr_t) start) ptr[2] = *ptr--;
 					start[0] = '0';
 					start[1] = ',';
					fprintf(stderr, "Missing number in AD, adding 0, : %s", Value);
				}
		}
}

int main (int argc, char *argv[])
{
	size_t BufferSize = 0UL;
	int res = 1;
	
	if (argc != 2) {
		fprintf(stderr, "Usage: %s vcf file\n", argv[0]);
		goto bail;
	}
	
	FILE * const restrict in = fopen(argv[1], "r");
	if (in == NULL) {
		fprintf(stderr, "Unable to open vcf file %s\n", argv[1]);
		goto bail;
	}
	unsigned int line = 0U;
	while (!feof(in)) {
		line++;
		const ssize_t nread = getline(&Buffer, &BufferSize, in);
		if (nread == -1) {
			if (!feof(in)) fprintf(stderr, "getline error on line %u\n", line);
			break;
		}
		
		/* loop on header */
		if (Buffer[0] == '#') { 
			printf("%s", Buffer);
			continue;
		}
		const char *chr, *pos, *name, *wt_allele, *alt_allele, *quality, *filter, *info, *format;
 		char* LineRemaining = NULL;
		{
			char * restrict ptr = Buffer;
			chr = ptr;
			while (*ptr != '\t') ++ptr;
			pos = ++ptr;
			while (*ptr != '\t') ++ptr;
			name = ++ptr;
			while (*ptr != '\t') ++ptr;
			wt_allele = ++ptr;
			while (*ptr != '\t') ++ptr;
			alt_allele = ++ptr;
			while (*ptr != '\t') ++ptr;
			quality = ++ptr;
			while (*ptr != '\t') ++ptr;
			filter = ++ptr;
			while (*ptr != '\t') ++ptr;
			info = ++ptr;
			while (*ptr != '\t') ++ptr;
			format = ++ptr;
			while (*ptr != '\t') ++ptr;
			LineRemaining = ++ptr;
		}
				
		/* Homozygote or heterozygote */
		const char * alt_allele2 = NULL;
		{
			const char * restrict ptr = alt_allele;
			if (*ptr == ',') {
					fprintf(stderr, "Found comma instead of first alternate allele on line %u\n\t%s\n", line, Buffer);
					continue;
			}
			while (*ptr != '\t') {
				if ( *ptr == ',' ) { 
					alt_allele2 = ptr + 1;
					break;
				}
				++ptr;
			}
		}
		
		const int chr_len = (uintptr_t) pos - (uintptr_t) chr - 1;
		const int pos_len = (uintptr_t) name - (uintptr_t) pos - 1;
		const int name_len = (uintptr_t) wt_allele - (uintptr_t) name - 1;
		const int wt_allele_len = (uintptr_t) alt_allele - (uintptr_t) wt_allele - 1;
		const int alt_allele_len = (uintptr_t) quality - (uintptr_t) alt_allele - 1;
		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Check some frequent errors
		
		/* 1. missing alt allele */
		if (*alt_allele == '\t') {
			fprintf(stderr, "No allele on line %u\n\t%s\n", line, Buffer);
			continue;
		}
		
		/* 2. ? characters */
		{
			const char * restrict ptr = alt_allele;
			int errFound = 0;
			while (*ptr != '\t') {
				if ( *ptr == '?' ) {
					fprintf(stderr, "'?' symbol found in allele on line %u\n\t%s\n", line, Buffer);
					errFound = 1;
					break;
				}
				++ptr;
			}
			if (errFound) continue;
		}
		
		/* 3. space characters */
		{
			const char * restrict ptr = alt_allele;
			int errFound = 0;
			while (*ptr != '\t') {
				if ( *ptr == ' ' ) {
					fprintf(stderr, "Space found in allele on line %u\n\t%s\n", line, Buffer);
					errFound = 1;
					break;
				}
				++ptr;
			}
			if (errFound) continue;
		}
		
		/* 4. Region= undefined character */
		{
				char * restrict ptr = (char*) info;
				while (ptr < (format-(1+9))) {
						if (strcmp(ptr, "RegionOK=") == 0) {
								if (ptr[9] < '0' || ptr[9] > '1') {
										fprintf(stderr, "REGION OK UNDEFINED VALUE (ascii %i) IN %s", (int) ptr[9], Buffer);
										ptr[9] = '?';
								}
						}
						ptr++;
				}
		}
		
		/* 5. Identical to the reference */
// 		const int NumericalPos = atoi(pos);
		//////////////////////////////////////////////////////////////////////////////////////////////////
		// HETEROZYGOTE
		if (alt_allele2) {
			const int alt_allele1_len = (uintptr_t) alt_allele2 - (uintptr_t) alt_allele - 1; 
			const int alt_allele2_len = (uintptr_t) quality - (uintptr_t) alt_allele2 - 1;
			/*printf("%.*s//%.*s//%.*s//%.*s//%.*s\n",
			 *				       chr_len, chr, 
			 *				       pos_len, pos,
			 *				       wt_allele_len, wt_allele,
			 *				       alt_allele1_len, alt_allele,
			 *				       alt_allele2_len, alt_allele2);*/
			///////////////////////////////////////////////////////////////////////////////////////////////
			// IDENTICAL LENGTH
			if ( alt_allele1_len == alt_allele2_len && alt_allele1_len == wt_allele_len ) {
				if (strncmp(wt_allele, alt_allele, wt_allele_len) == 0) {
					if (strncmp(wt_allele, alt_allele2, wt_allele_len) == 0) {
						fprintf(stderr, "Identical strings: '%.*s' <-> '%.*s' <-> '%.*s'\nRemoving: %s", wt_allele_len, wt_allele, alt_allele_len, alt_allele, alt_allele2_len, alt_allele2, Buffer);
						continue;
					}
				}
// 				if (wt_allele_len > 1) {
// 					char CopyRemain[4096];
// 					LineRemaining[-1] = '\0';
// 					for (int i=0; i<wt_allele_len; i++) {
// 						if (wt_allele[i] != alt_allele[i]) {
// 							if (wt_allele[i] != alt_allele2[i]) {
// 								if (alt_allele[i] == alt_allele2[i]) {
// 									/* Correct GT */
// 									strcpy(CopyRemain, LineRemaining);
// 									char * restrict ptr = CopyRemain;
// 									//fprintf(stdout, "IN: %s", Buffer);
// 									while (ptr[3] != '\n') {
// 										if (ptr[0] == '1' && ptr[1] == '/' && ptr[2] == '2') {
// 											ptr[0] = '0';
// 											ptr[2] = '1';
// 											break;
// 										}
// 										ptr++;
// 									}
// 									fprintf(stdout, "%.*s\t%u\t%.*s\t%c\t%c\t%s\t%s",
// 													chr_len, chr,
// 													NumericalPos+i,
// 													name_len, name,
// 													(int) wt_allele[i],
// 													(int) alt_allele[i],
// 													quality,
// 													CopyRemain);
// 								}
// 								else {
// 									CorrectAD(LineRemaining);
// 									fprintf(stdout, "%.*s\t%u\t%.*s\t%c\t%c,%c\t%s\t%s", 
// 													chr_len, chr,
// 													NumericalPos+i,
// 													name_len, name,
// 													(int) wt_allele[i],
// 													(int) alt_allele[i],
// 													(int) alt_allele2[i],
// 													quality,
// 													LineRemaining);
// 								}
// 							}
// 							else {
// 								if (wt_allele[i] == alt_allele2[i] ) {
// 									/* Correct GT */
// 									strcpy(CopyRemain, LineRemaining);
// 									char * restrict ptr = CopyRemain;
// 									//fprintf(stdout, "IN: %s", Buffer);
// 									while (ptr[3] != '\n') {
// 										if (ptr[0] == '1' && ptr[1] == '/' && ptr[2] == '2') {
// 											ptr[0] = '0';
// 											ptr[2] = '1';
// 											break;
// 										}
// 										ptr++;
// 									}
// 									LineRemaining[-1] = '\0';
// 									fprintf(stdout, "%.*s\t%u\t%.*s\t%c\t%c\t%s\t%s",
// 													chr_len, chr,
// 													NumericalPos+i,
// 													name_len, name,
// 													(int) wt_allele[i],
// 													(int) alt_allele[i], quality, CopyRemain);
// 									
// 								}
// 								else {
// 									CorrectAD(LineRemaining);
// 									fprintf(stdout, "%.*s\t%u\t%.*s\t%c\t%c\t%s\t%s",
// 													chr_len, chr,
// 													NumericalPos+i,
// 													name_len, name,
// 													(int) wt_allele[i],
// 													(int) alt_allele[i],
// 													quality,
// 													LineRemaining);
// 								}
// 							}
// 						}
// 						else {
// 							if (wt_allele[i] != alt_allele2[i]) {
// 								/* Correct GT */
// 								strcpy(CopyRemain, LineRemaining);
// 								char * restrict ptr = CopyRemain;
// 								//fprintf(stdout, "IN: %s", Buffer);
// 								while (ptr[3] != '\n') {
// 									if (ptr[0] == '1' && ptr[1] == '/' && ptr[2] == '2') {
// 										ptr[0] = '0';
// 										ptr[2] = '1';
// 										break;
// 									}
// 									ptr++;
// 								}
// 								fprintf(stdout, "%.*s\t%u\t%.*s\t%c\t%c\t%s\t%s",
// 												chr_len, chr,
// 												NumericalPos+i,
// 												name_len, name,
// 												(int) wt_allele[i],
// 												(int) alt_allele2[i],
// 												quality,
// 												CopyRemain);
// 							}
// 						}
// 					}
// 					continue;
// 					//fprintf(stderr, "Equal length>1 on heterozygote\n%s", Buffer);
// 					//exit(1);
// 				}
			}
// 			
// 			///////////////////////////////////////////////////////////////////////////////////////////////
// 			// DIFFERENT LENGTHS
// 			int minLength = (wt_allele_len < alt_allele1_len) ? wt_allele_len : alt_allele1_len;
// 			minLength = (minLength < alt_allele2_len) ? minLength : alt_allele2_len;
// 			if (minLength > 1) {
// 				int i=1;
// 				while (i<minLength-1) {
// 					if (wt_allele[i] != alt_allele[i]) break;
// 					if (wt_allele[i] != alt_allele2[i]) break;
// 					++i;
// 				}
// 				int needGTchange = 0;
// 				
// 				if (i > 1) {
// 					if ( wt_allele_len == alt_allele1_len && strncmp(wt_allele+i, alt_allele+i, wt_allele_len-i) == 0 ) {
// 						fprintf(stderr, "THIS IS NO LONGER HETEROZYGOUS !!! ALTER A\n%s", Buffer); 
// 						needGTchange = 1;
// 						
// 					}
// 					if ( wt_allele_len == alt_allele2_len && strncmp(wt_allele+i, alt_allele2+i, wt_allele_len-i) == 0 ) {
// 						fprintf(stderr, "THIS IS NO LONGER HETEROZYGOUS !!! ALTER B\n%s", Buffer); 
// 						needGTchange = 2;
// 					}
// 					LineRemaining[-1] = '\0';
// 					if (needGTchange > 0) {
// 						char * restrict ptr = LineRemaining;
// 						while (ptr[3] != '\n') {
// 							if (ptr[0] == '1' && ptr[1] == '/' && ptr[2] == '2') {
// 								if (needGTchange == 1) {
// 									ptr[0] = '0';
// 									ptr[2] = '1';
// 								}
// 								else {
// 									ptr[0] = '1';
// 									ptr[2] = '0';
// 								}
// 								/*ptr += 3;
// 								 *									while (*ptr != '\n') {
// 								 *										
// 							}*/
// 								break;
// 							}
// 							ptr++;
// 						}
// 						if (needGTchange == 1) {
// 							fprintf(stdout, "%.*s\t%u\t%.*s\t%.*s\t%.*s\t%s\t%s", 
// 											chr_len, chr,
// 											NumericalPos+i,
// 											name_len, name,
// 											wt_allele_len-i, wt_allele+i,
// 											alt_allele2_len-i, alt_allele2+i, quality, LineRemaining);
// 							
// 						}
// 						else {
// 							fprintf(stdout, "%.*s\t%u\t%.*s\t%.*s\t%.*s\t%s\t%s",
// 											chr_len, chr,
// 											NumericalPos+i,
// 											name_len, name,
// 											wt_allele_len-i, wt_allele+i,
// 											alt_allele1_len-i, alt_allele+i, quality,
// 											LineRemaining);
// 						}
// 						continue;
// 					}
// 					
// 					/* Correct AD 1/2, adding missing number if needed */
// 					CorrectAD(LineRemaining);
// 					
// 					fprintf(stdout, "%.*s\t%u\t%.*s\t%.*s\t%.*s,%.*s\t%s\t%s",
// 									chr_len, chr,
// 									NumericalPos+i,
// 									name_len, name,
// 									wt_allele_len-i, wt_allele+i,
// 									alt_allele1_len-i, alt_allele+i,
// 									alt_allele2_len-i, alt_allele2+i, quality, LineRemaining);
// 					continue;
// 				}
// 				else {
// 					if ( wt_allele_len == alt_allele1_len && strncmp(wt_allele, alt_allele, wt_allele_len) == 0 ) {
// 						fprintf(stderr, "THIS IS NO LONGER HETEROZYGOUS !!! ALTER C\n %s", Buffer);
// 						needGTchange = 1;
// 						
// 					}
// 					if ( wt_allele_len == alt_allele2_len && strncmp(wt_allele, alt_allele2, wt_allele_len-i) == 0 ) {
// 						fprintf(stderr, "THIS IS NO LONGER HETEROZYGOUS !!! ALTER D\n %s", Buffer);
// 						needGTchange = 2;
// 					}
// 					
// 					if (needGTchange > 0) {
// 						LineRemaining[-1] = '\0';
// 						char * restrict ptr = LineRemaining;
// 						while (ptr[3] != '\n') {
// 							if (ptr[0] == '1' && ptr[1] == '/' && ptr[2] == '2') {
// 								if (needGTchange == 1) {
// 									ptr[0] = '0';
// 									ptr[2] = '1';
// 								}
// 								else {
// 									ptr[0] = '1';
// 									ptr[2] = '0';
// 								}
// 								break;
// 							}
// 							ptr++;
// 						}
// 						if (needGTchange == 1) {
// 							fprintf(stdout, "%.*s\t%u\t%.*s\t%.*s\t%.*s\t%s\t%s", 
// 											chr_len, chr,
// 											NumericalPos,
// 											name_len, name,
// 											wt_allele_len, wt_allele,
// 											alt_allele2_len, alt_allele2, 
// 											quality, LineRemaining);
// 							
// 						}
// 						else {
// 							fprintf(stdout, "%.*s\t%u\t%.*s\t%.*s\t%.*s\t%s\t%s",
// 											chr_len, chr,
// 											NumericalPos,
// 											name_len, name,
// 											wt_allele_len, wt_allele,
// 											alt_allele1_len, alt_allele,
// 											quality, LineRemaining);
// 						}
// 						continue;
// 					}
// 				}
// 			}
// 			
// 			/* Correct AD 1/2, adding missing number if needed */
// 			CorrectAD(LineRemaining);
// 			
		}
		//////////////////////////////////////////////////////////////////////////////////////////////////
		// HOMOZYGOTE
		else {
			if ( (wt_allele_len == alt_allele_len) && wt_allele_len > 1 ) {
				if (strncmp(wt_allele, alt_allele, wt_allele_len) == 0) {
					fprintf(stderr, "Identical strings: '%.*s' <-> '%.*s'\nRemoving: %s", wt_allele_len, wt_allele, alt_allele_len, alt_allele, Buffer);
					continue;
				}
// 				else {
// 					fprintf(stderr, "Identical length: '%.*s' <-> '%.*s'\n", wt_allele_len, wt_allele, alt_allele_len, alt_allele);
// 					LineRemaining[-1] = '\0';
// 					for (int i=0; i<wt_allele_len;i++) {
// 						if (wt_allele[i] != alt_allele[i]) {
// 							printf("%.*s\t%u\t%.*s\t%c\t%c\t%s\t%s",
// 										 chr_len, chr,
// 										 NumericalPos + i,
// 										 name_len, name,
// 							       (int) wt_allele[i],
// 										 (int) alt_allele[i],
// 										 quality,
// 										 LineRemaining);
// 						}
// 					}
// 					continue;
// 				}
			}
// 			
// 			const int minLength = (wt_allele_len < alt_allele_len) ? wt_allele_len : alt_allele_len;
// 			if (minLength > 1 ) {
// 				
// 				int i=1;
// 				while (i<minLength-1) {
// 					if (wt_allele[i] != alt_allele[i]) break;
// 					++i;
// 				}
// 				if (i > 1 ) {
// 					LineRemaining[-1] = '\0';
// 					fprintf(stdout, "%.*s\t%u\t%.*s\t%.*s\t%.*s\t%s\t%s",
// 									chr_len, chr,
// 								  NumericalPos+i,
// 								  name_len, name,
// 								  wt_allele_len-i, wt_allele+i,
// 								  alt_allele_len-i, alt_allele+i, quality, 
// 								  LineRemaining);
// 					continue;
// 				}
// 			}
		}
		
		//////////////////////////////////////////////////////////////////////////////////////////////////
		// NO PROBLEM
		printf("%s", Buffer);
		
	}
	
	fclose(in);
	
	bail: ;
	return res;
}

