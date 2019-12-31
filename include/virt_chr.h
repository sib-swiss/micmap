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
#ifndef _VIRT_CHR_H
#define _VIRT_CHR_H

#define kGENOME_DATA_SIZE 0x100000000
#define kMAXCHRcnt 32  /* unmapped=0, chr=1..22,X,Y,M plus some unmapped and extra stuff */

typedef	struct VIRTUALCHR_struct  VIRTUALCHR;
struct VIRTUALCHR_struct
{
	unsigned int chr;
	unsigned int offset;
	unsigned int len;
	char AC[32];
	char SAMname[32];
};
extern unsigned int chrMax;

typedef struct REVVIRTUALCHR_struct
{
	unsigned int nbChr;
	unsigned int chr[4];
} REVVIRTUALCHR;

extern REVVIRTUALCHR *revVchr;
extern unsigned char gIUPACrevcomp[128] __attribute__((aligned(64)));

void initializeVirtualChromosomes(const char *, VIRTUALCHR *);
void freeVirtualChromosome();
void parseVirtualChromosomes(char *, unsigned int, VIRTUALCHR *);
void initRevCompTable(unsigned char *);
#endif

/* vim: tabstop=2 shiftwidth=2
 */
