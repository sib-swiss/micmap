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
#include <string.h>

int main(int argc, char *argv[])
{
	unsigned int line = 0U;
	char ChrString[16];
	size_t Account[24][3];
		
		if (argc != 3) {
			fprintf(stderr, "Usage: %s [bed file] [split size]\n", argv[0]);
			return 1;
		}		
		int SplitSize = atoi(argv[2]);
		
		memset(&Account[0][0], 0, 3*24*sizeof(size_t));
		
		FILE * bedFile = fopen(argv[1], "r");
		if (bedFile) {
			while (!feof(bedFile)) {
				unsigned int ChrRangeBegin, ChrRangeEnd, Chr;
				const int num = fscanf(bedFile, "%s\t%u\t%u\n", ChrString, &ChrRangeBegin, &ChrRangeEnd);
				if (ChrString[0] == 'X') {
					Chr = 23;
				}
				else if (ChrString[0] == 'Y') {
					Chr = 24;
				}
				else {
					Chr = (unsigned int) atoi(ChrString);
				}
				if (num == 3) {
					Account[Chr-1][0]++;
					Account[Chr-1][1] += ChrRangeEnd - ChrRangeBegin + 1;
					line++;
				}
				else {
						fprintf(stderr, "Error parsing line %u\n", line);
						return 1;
				}
				
			}
			
			fclose(bedFile);
		}
		
		size_t Overall[2] = { 0UL, 0UL};
		int nTot = 0;
		for (int i=1; i<=24; i++) {
				Account[i-1][2] = (Account[i-1][1] + SplitSize - 1) / SplitSize;
				printf("Chr %i\t%12lu\t%12lu\t%lu\n", i, Account[i-1][0], Account[i-1][1], Account[i-1][2]);
				Overall[0] += Account[i-1][0];
				Overall[1] += Account[i-1][1];
				nTot += Account[i-1][2];
		}
		printf("\n      \t%'12lu\t%'12lu\t%i\n", Overall[0], Overall[1], nTot);
		
		bedFile = fopen(argv[1], "r");
		int oldChr = -1;
		FILE * out = NULL;
		int Part = 0;
		int index = 0;
		size_t soFar = 0UL;
		char outFname[64];
		size_t Limit = Account[index][1] / Account[index][2];
		if (bedFile) {
			while (!feof(bedFile)) {
				unsigned int ChrRangeBegin, ChrRangeEnd, Chr;
				const int num = fscanf(bedFile, "%s\t%u\t%u\n", ChrString, &ChrRangeBegin, &ChrRangeEnd);
				if (ChrString[0] == 'X') {
					Chr = 23;
				}
				else if (ChrString[0] == 'Y') {
					Chr = 24;
				}
				else {
					Chr = (unsigned int) atoi(ChrString);
				}
				
				if (num == 3) {
					const int Length = ChrRangeEnd - ChrRangeBegin + 1;
					soFar += Length;
					if ( (oldChr == -1) || Chr != oldChr || soFar >= Limit) {
						Part++;
						if (out != NULL) fclose(out);
						snprintf(outFname, 64, "SplitChr_%i_%i.bed", Chr, Part);
						out = fopen(outFname, "w");
						if (out == NULL) return -1;
						if (soFar >= Limit) soFar = 0UL;
						if ( (oldChr != -1) && (Chr != oldChr) ) {
							index++;
							Limit = Account[index][1] / Account[index][2];
							soFar = 0UL;
						}
					}	
					
					fprintf(out, "%i\t%u\t%u\n", Chr, ChrRangeBegin, ChrRangeEnd);
						
					
				}
				else {
						fprintf(stderr, "Error parsing line %u\n", line);
						return 1;
				}
				
				oldChr = Chr;
			}
			
			fclose(out);
			fclose(bedFile);
		}
		
		
		
    return 0;
}
