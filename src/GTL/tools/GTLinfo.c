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
#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include "GTL.h"
#include "gtlVersion.h"

typedef struct _decodedPair {
	unsigned long ordinal;
	unsigned int tag1pos;
	unsigned int tag2pos;
	int taglen1;
	int taglen2;
	int reverseTAG1;
	int reverseTAG2;
} decodedPair_t, *decodedPair_p_t;

int main(int argc, char * argv[])
{
	TLBHEADER th;
	TLBDATA td;
	int err=1;

	if (argc!= 2) {
		fprintf(stderr, "Usage: %s [gtl file]\n" GTL_VERSION_FULL_STRING "\n", argv[0]);
		err = 0;
		goto bail;
	}

	const char * const restrict fn = argv[1];

	int fd = open(fn, O_RDONLY);
	allocBlock(&td, true);

	int iBlock = 0U;
	while (1) {
		ssize_t res = read(fd,&th,sizeof(TLBHEADER));
		if (res == 0)
		{
			err = 0;
			break;
		}
		if (res != sizeof(TLBHEADER))
			goto bail;
		unsigned char hcs = simple8bitCS((unsigned char *) &th, sizeof(TLBHEADER));
		if (hcs != 0)
		{
			fprintf(stderr,"Checksum mismatch: %02x %02x\n",th.headerCS,hcs);
			goto bail;
		}

		err = readDecompressBlock(fd, &th, &td, 0);
		if (err != 0)
			goto bail;

		printf("Block %u:\tchr %2u [%u,%u]", iBlock++, th.chr, th.minPos, th.maxPos);

		if ((td.header.flags_7_0 & kTLBHflagPairedReads) == 0) {
			printf(" SINGLE %u tags\n", td.cnt);
		}
		else {
			decodedPair_t dp;
			unsigned int minPos = 0xFFFFFFFF;
			unsigned int maxPos = 0;
			printf(" PAIRED %u tags", td.cnt);
			for (int i=0; i<td.cnt; i++)
			{
				if (td.header.flags_7_0 & kTLBHflagFixedLength)
				{
					dp.taglen1 = td.lengthInfo;
					dp.taglen2 = td.lengthInfo;
				}
				else
				{
					dp.taglen1 = td.len[2*i];
					dp.taglen2 = td.len[2*i+1];
				}
				dp.ordinal = 0;
				for (unsigned int j = 0; j < 8; j++)
					dp.ordinal = (dp.ordinal << 8) | (unsigned long) *((unsigned char *) td.hdr++);

				if (td.ppd[0][i].tag2offset & k32bREVCOMP_TAG1_BIT)
					dp.reverseTAG1 = 1;
				else
					dp.reverseTAG1 = 0;

				if (td.ppd[0][i].tag2offset & k32bREVCOMP_TAG2_BIT)
					dp.reverseTAG2 = 1;
				else
					dp.reverseTAG2 = 0;
				const int delta = td.ppd[0][i].tag2offset & k32bOFFSET_VALUE_MASK;
				if (td.ppd[0][i].tag2offset & k32bNEGATIVE_OFFSET_BIT) {
					dp.tag1pos = td.ppd[0][i].tag1pos - delta;
					dp.tag2pos = td.ppd[0][i].tag1pos;
				}
				else {
					dp.tag1pos = td.ppd[0][i].tag1pos;
					dp.tag2pos = td.ppd[0][i].tag1pos + delta;
				}

				if (dp.tag1pos < minPos) minPos = dp.tag1pos;
				const unsigned int lastBase = dp.tag1pos;// + kGAPPED_ALIGN_GENOME_LENGTH;
				if (lastBase > maxPos) maxPos = lastBase;
			}

			printf("\t[%u:%u]", minPos, maxPos);
			printf("\tHeader (%s)", compress_name(td.readHdr.compressor));
			printf("\tQuality scores (%s)\n", compress_name(td.readQual.compressor));
		}

	}

	err = 0;

	bail: ;
	return err;
}
//---------------------------------------------------------------

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
