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
#include <errno.h>

#include <assert.h>
#include "GTL.h"
#include "gtlVersion.h"

typedef struct _elem_t {
	TLBHEADER th;
	PAIRDATA pd[kBlockMaxCount];
  char *fName;
  int fd;
  unsigned int cur;
  unsigned int cnt;
  int unlink;
} elem_t, *elem_p_t;


//=================================================================================================

static void
writeBuffer(TLBHEADER *thp, unsigned int cnt, const char* fname, TLBDATA *wbp)
{
	unsigned int i;
	if (cnt == 0)
		return;
	zeroBlock(wbp);
	memcpy(&wbp->header,thp,sizeof(TLBHEADER));
	wbp->cnt = cnt;
	PAIRDATA **pdp = wbp->pdp;
	// Get the various info from the block components
	for (i = 0; i < cnt; i++)
	{
		if (i == 0)
		{
			wbp->lengthInfo = pdp[i]->len1;
			wbp->nbMatches = pdp[i]->nbMatches;
			assert(pdp[i]->nbMatches != 0);
		}
		assert(wbp->nbMatches == pdp[i]->nbMatches);
		if (pdp[i]->len1 != wbp->lengthInfo || pdp[i]->len2 != wbp->lengthInfo)
			wbp->lengthInfo = 0;
		if ((thp->flags_7_0 & kTLBHflagPairPosBlock) != 0
				|| (thp->flags_7_0 & kTLBHflagPairChrPosBlock) != 0)
		{
			unsigned int minp;
			if ((thp->flags_7_0 & kTLBHflagPairedReads) == 0)
			{
				minp = pdp[i]->pcpd[0].tag1pos;
			}
			else
			{
				if (pdp[i]->pcpd[0].tag1pos != 0)
				{
					if (pdp[i]->pcpd[0].tag2pos != 0)
					{
						if (pdp[i]->pcpd[0].tag1pos <= pdp[i]->pcpd[0].tag2pos)
							minp = pdp[i]->pcpd[0].tag1pos;
						else
							minp = pdp[i]->pcpd[0].tag2pos;
					}
					else
						minp = pdp[i]->pcpd[0].tag1pos;
				}
				else
					minp = pdp[i]->pcpd[0].tag2pos;
			}
			if (wbp->header.minPos > minp)
				wbp->header.minPos = minp;
			if (wbp->header.maxPos < minp)
				wbp->header.maxPos = minp;
		}
		wbp->hdrSize += pdp[i]->hlen1;
		wbp->hdrSize += pdp[i]->hlen2;
		wbp->qualSize += pdp[i]->len1;
		wbp->qualSize += pdp[i]->len2;
		for (unsigned int j = 0; j < kMAXnbMatches; j++)
		{
			wbp->mmSize[j] += pdp[i]->ml[j];
			wbp->cigarSize[j] += pdp[i]->cl[j];
		}
		wbp->p4Size += pdp[i]->pl;
	}
	unsigned char *dataBlock = packBlock(wbp, wbp);
	int fd = open(fname, O_WRONLY|O_CREAT|O_APPEND, 0664);
	if (fd == -1)
	{
		perror("open:");
		exit(1);
	}
	if (write(fd,dataBlock,wbp->header.blockLength) != wbp->header.blockLength)
	{
		perror("write failed");
		exit(1);
	}
	close(fd);
	fprintf(stderr,".");
}
//---------------------------------------------------------------

static void
init_elem(elem_p_t e, TLBDATA *tdp)
{
	e->fd = open(e->fName, O_RDONLY);
	ssize_t res = read(e->fd,&e->th,sizeof(TLBHEADER));
	if (res != sizeof(TLBHEADER))
		exit(1);
	unsigned char hcs = simple8bitCS((unsigned char *) &e->th, sizeof(TLBHEADER));
	if (hcs != 0)
	{
		fprintf(stderr,"Checksum mismatch: %02x %02x\n",e->th.headerCS,hcs);
		exit(1);
	}
	int cnt = readDecodeBlock(e->fd,&e->th,e->pd,tdp);
	if (cnt < 0)
		exit(1);
	e->cnt = cnt;
	e->cur = 0;
}
//---------------------------------------------------------------

static void
get_next(elem_p_t e, TLBDATA *tdp)
{
	if (e->fd == -1)
		return;
	e->cur += 1;
	if (e->cur < e->cnt)
		return;
	ssize_t res = read(e->fd,&e->th,sizeof(TLBHEADER));
	if (res == 0)
	{
		close(e->fd);
		e->fd = -1;
		e->cur = 0;
		e->cnt = 0;
		return;
	}
	if (res != sizeof(TLBHEADER))
		exit(1);
	unsigned char hcs = simple8bitCS((unsigned char *) &e->th, sizeof(TLBHEADER));
	if (hcs != 0)
	{
		fprintf(stderr,"Checksum mismatch: %02x %02x\n",e->th.headerCS,hcs);
		exit(1);
	}
	int cnt = readDecodeBlock(e->fd,&e->th,e->pd,tdp);
	if (cnt < 0)
		exit(1);
	e->cnt = cnt;
	e->cur = 0;
}
//---------------------------------------------------------------

static void
merge(elem_p_t eTbl, unsigned int eTblLen, const char *oFname, TLBDATA *tdp)
{
	PAIRDATA *pd = tdp->pd;
	PAIRDATA **pdp = tdp->pdp;
	unsigned int i;
	for (i = 0; i < kBlockMaxCount; i++)
		pdp[i] = pd + i;
	unsigned int cnt = 0;
	int more;
	PAIRDATA gMin;
	int gMinSet = 0;
	do {
		PAIRDATA min;
		more = 0;
		for (i = 0; i < eTblLen; i++) {
			elem_p_t e = eTbl + i;
			if (e->fd != -1) {
				if (more) {
					if (comparePAIRDATA(e->pd + e->cur,&min) < 0)
						min = e->pd[e->cur];
				} else {
					more = 1;
					min = e->pd[e->cur];
				}
			}
		}
		if (more) {
			if (!gMinSet) {
				gMinSet = 1;
				gMin = min;
			}
			if (comparePAIRDATA(&min,&gMin) < 0) {
				fprintf(stderr, "merge: bug min < gMin\n");
				//exit(1);
			}
			gMin = min;
			for (i = 0; i < eTblLen; i++) {
				elem_p_t e = eTbl + i;
				while (e->fd != -1 && comparePAIRDATA(&min,e->pd + e->cur) == 0) {
					pd[cnt++] = e->pd[e->cur];
					if (cnt >= kBlockMaxCount)
					{
						writeBuffer(&e->th,cnt,oFname,tdp);
						cnt = 0;
					}
					get_next(e,tdp);
				}
			}
		}
	} while (more);
	writeBuffer(&eTbl->th,cnt,oFname,tdp);
	fprintf(stderr,"\n");
}
//---------------------------------------------------------------

int
main (int argc, char **argv)
{
	int c;
	unsigned int i;
	char oFname[512] = "";
	/* --------- process arguments ----------------------*/

	while ((c = getopt (argc, argv, "o:s:z:")) != -1)
	switch (c)
	{
		case 'o':
			sscanf(optarg,"%511s",oFname);
			break;

		case 's':
			gCompress.scoredef = compress_find(optarg);
			break;

		case 'z':
			gCompress.def = compress_find(optarg);
			break;
	}

	if (oFname[0] == 0 || gCompress.scoredef < 0 || gCompress.def < 0 || argc == optind)
	{
		printf("usage:\n\n");
		printf("GTLmerge -o file ...\n\n");
		printf("           -o file               : output file\n");
		printf("           -z compression        : compression to use for data, default is '%s'\n", compress_name(gCompress.def));
		printf("           -s compression        : compression to use for scores, default is '%s'\n", compress_name(gCompress.scoredef));
		// TODO add an option either generic "level" or per compressor flag (like "hc" for lz4)
		printf("                    available compression :\n");
		for (i = 0; compress_dict[i].name != NULL; i++)
			printf("                     %11.11s : %s\n", compress_dict[i].name, compress_dict[i].description);
		printf("\n");
		printf(GTL_VERSION_FULL_STRING "\n");
		printf("(c) Nicolas Guex & Christian Iseli & Thierry Schuepbach 2015-2018\n");
		return(1);
	}

	elem_p_t eTbl = malloc((argc - optind) * sizeof(elem_t));
	TLBDATA td;
	unsigned int eTblLen = 0;
	allocBlock(&td, true);
	while (optind < argc) {
		eTbl[eTblLen].fName = strdup(argv[optind]);
		if (eTbl[eTblLen].fName == NULL) {
			perror("main: strdup");
			return 1;
		}
		init_elem(eTbl + eTblLen,&td);
		eTblLen += 1;
		optind += 1;
	}
	merge(eTbl,eTblLen,oFname,&td);
	freeBlock(&td);

	return 0;
} /* main */
//---------------------------------------------------------------

/* vim: tabstop=2 shiftwidth=2
 */
