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
#ifndef _STREAM_H
#define _STREAM_H
#include "Aligner.h"

typedef struct StreamData {
	const outputPairBlockHeader_t * BlockHeader;
	const AlignmentResults_t * AlnResults[4];
	const char * TagHdr[2];
	const char * TagQS[2];
	const char * Tags[2];
	const char * SoftClips[4];	
	size_t SoftClipLengths[4];
} StreamData_t;

int GetStreamData(const char * BufferPtr, StreamData_t * data);
void PrintAlnResult(const AlignmentResults_t * AlnResult, const Genome_t * const genome, const int taglen);
void PrintStreamData(const StreamData_t * const data, const Genome_t * const genome);

#endif

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
