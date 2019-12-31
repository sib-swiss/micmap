/*
 * ------------------------------------------------------------------------------------------------------------------------
 *
 *                                * micmap *
 *      Mapping of short reads in fastq format onto a reference
 *
 *
 *  Copyright (C) SIB  - Swiss Institute of Bioinformatics,  2014-2019 Nicolas Guex, Thierry Schuepbach and Christian Iseli
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
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include "config.h"
#include "virt_chr.h"

REVVIRTUALCHR *revVchr = NULL;
unsigned int chrMax;
unsigned char gIUPACrevcomp[128] __attribute__((align(64)));

//---------------------------------------------------------------
void initRevCompTable(unsigned char *gIUPACrevcomp)
{
	memset(gIUPACrevcomp,'N',128);
	gIUPACrevcomp['A'] = 'T';
	gIUPACrevcomp['T'] = 'A';
	gIUPACrevcomp['U'] = 'A';
	gIUPACrevcomp['G'] = 'C';
	gIUPACrevcomp['C'] = 'G';
	gIUPACrevcomp['Y'] = 'R';
	gIUPACrevcomp['R'] = 'Y';
	gIUPACrevcomp['S'] = 'S';
	gIUPACrevcomp['W'] = 'W';
	gIUPACrevcomp['K'] = 'M';
	gIUPACrevcomp['M'] = 'K';
	gIUPACrevcomp['B'] = 'V';
	gIUPACrevcomp['D'] = 'H';
	gIUPACrevcomp['H'] = 'D';
	gIUPACrevcomp['V'] = 'B';

	gIUPACrevcomp['a'] = 't';
	gIUPACrevcomp['t'] = 'a';
	gIUPACrevcomp['u'] = 'a';
	gIUPACrevcomp['g'] = 'c';
	gIUPACrevcomp['c'] = 'g';
	gIUPACrevcomp['y'] = 'r';
	gIUPACrevcomp['r'] = 'y';
	gIUPACrevcomp['s'] = 's';
	gIUPACrevcomp['w'] = 'w';
	gIUPACrevcomp['k'] = 'm';
	gIUPACrevcomp['m'] = 'k';
	gIUPACrevcomp['b'] = 'v';
	gIUPACrevcomp['d'] = 'h';
	gIUPACrevcomp['h'] = 'd';
	gIUPACrevcomp['v'] = 'b';
	gIUPACrevcomp['n'] = 'n';

} // initRevCompTable
//---------------------------------------------------------------

// get data from config file
//chr	len	virtChr	offset	AC	SAMname
//1	249250621	1	0	NC_000001.10	chr1
static void
parseLine(char *buf, unsigned int bufLen, VIRTUALCHR *virtchr)
{
	if (bufLen == 0)
		return;
	if (buf[0] == '#' || buf[0] == 'c')
		return;
	char *sep = strchr(buf, '\t');
	if (sep == NULL)
	{
		fprintf(stderr, "expected chr, got nothing\n");
		exit(1);
	}
	unsigned int chr = atoi(buf);
	if (chr > chrMax)
		chrMax = chr;
	char *cur = sep + 1;
	sep = strchr(sep + 1, '\t');
	if (sep == NULL)
	{
		fprintf(stderr, "expected chr len, got nothing\n");
		exit(1);
	}
	unsigned int len = atoi(cur);
	cur = sep + 1;
	sep = strchr(sep + 1, '\t');
	if (sep == NULL)
	{
		fprintf(stderr, "expected virtChr, got nothing\n");
		exit(1);
	}
	unsigned int virtChr = atoi(cur);
	if (virtChr > 15)
	{
		fprintf(stderr, "virtual chr must be between 0 and 15 inclusive, got %d\n",virtChr);
		exit(1);
	}
	cur = sep + 1;
	sep = strchr(sep + 1, '\t');
	if (sep == NULL)
	{
		fprintf(stderr, "expected offset, got nothing\n");
		exit(1);
	}
	unsigned int offset = atoi(cur);
	cur = sep + 1;
	sep = strchr(sep + 1, '\t');
	if (sep == NULL)
	{
		fprintf(stderr, "expected AC, got nothing\n");
		exit(1);
	}
	virtchr[chr].chr = virtChr;
	virtchr[chr].offset = offset;
	virtchr[chr].len = len;
	if (revVchr[virtChr].nbChr > 3)
	{
		fprintf(stderr,"Too many chromosomes map to virtual chr%d\n",virtChr);
		exit(1);
	}
	revVchr[virtChr].chr[revVchr[virtChr].nbChr] = chr;
	revVchr[virtChr].nbChr += 1;
	unsigned int AClen = sep - cur;
	if (AClen > 31)
	{
		fprintf(stderr,"AC is too long (%d); max is 31\n",AClen);
		exit(1);
	}
	memcpy(virtchr[chr].AC,cur,AClen);
	virtchr[chr].AC[AClen] = 0;
	unsigned int SAMnameLen = bufLen - (sep - buf) - 1;
	if (SAMnameLen > 31)
	{
		fprintf(stderr,"AC is too long (%d); max is 31\n",SAMnameLen);
		exit(1);
	}
	memcpy(virtchr[chr].SAMname,sep+1,SAMnameLen);
	virtchr[chr].SAMname[SAMnameLen] = 0;
}
//---------------------------------------------------------------

void
parseVirtualChromosomes(char *buf, unsigned int bufLen, VIRTUALCHR *virtchr)
{
	unsigned int c = 0;
	chrMax = 0U;
	if (revVchr == NULL)
	{
		revVchr = calloc(16,sizeof(REVVIRTUALCHR));
	}
	while (c < bufLen)
	{
		char *eol = strchr(buf + c, '\n');
		unsigned int ll = bufLen - c;
		if (eol != NULL)
			ll = eol - (buf + c);
		parseLine(buf + c, ll, virtchr);
		c += ll + 1;
	}
}
//---------------------------------------------------------------

void
initializeVirtualChromosomes(const char *fn, VIRTUALCHR *virtchr)
{
	int fd = open(fn, O_RDONLY);
	if (fd == -1)
	{
		perror("open");
		exit(1);
	}
	struct stat sb;
	if (fstat(fd,&sb) != 0)
	{
		perror("fstat");
		exit(1);
	}
	char *readBuf = malloc(sb.st_size);
	if (readBuf == NULL)
	{
		perror("malloc");
		exit(1);
	}
	size_t rCnt = sb.st_size;
	char *rPos = readBuf;
	ssize_t res;
	virtchr[0].chr = 0;
	virtchr[0].offset = 0;
	virtchr[0].len = 0;
	strcpy(virtchr[0].AC,"unknown");
	strcpy(virtchr[0].SAMname,"unknown");
	while ((res = read(fd,rPos,rCnt)) != rCnt)
	{
		if (res == 0)
		{
			fprintf(stderr,"we miss %ld bits (got 0)\n", rCnt);
			exit(1);
		}
		if (res == -1)
		{
			if (errno != EINTR)
			{
				perror("read");
				exit(1);
			}
		}
		else
		{
			rPos += res;
			rCnt -= res;
		}
	}
	close(fd);
	parseVirtualChromosomes(readBuf, sb.st_size, virtchr);
	free(readBuf);
}
//---------------------------------------------------------------

void
freeVirtualChromosome()
{
	if (revVchr) free(revVchr);
}
//---------------------------------------------------------------

/* ------------------------------------------------------------------------------------ */
/* vim: tabstop=2 shiftwidth=2
 */
