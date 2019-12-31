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
#include "constants.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "stream.h"

int GetStreamData(const char * BufferPtr, StreamData_t * data)
{
	const uintptr_t Start = (uintptr_t) BufferPtr; 
	BufferPtr += 8;
	const register outputPairBlockHeader_t * const Header = (outputPairBlockHeader_t *) BufferPtr;
	data->BlockHeader = Header;
	BufferPtr += sizeof(outputPairBlockHeader_t);	

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Set the alignment results
	data->AlnResults[0] = NULL;
	data->AlnResults[1] = NULL;
	data->AlnResults[2] = NULL;
	data->AlnResults[3] = NULL;
	{
		register const int nAlignments = (int) Header->nAlignmentResults;
		for (int i=0; i<nAlignments; i++) {
			data->AlnResults[i] = (AlignmentResults_t *) BufferPtr;
			BufferPtr += sizeof(AlignmentResults_t);
		}
	}
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Set the headers
	data->TagHdr[0] = BufferPtr;
	BufferPtr += (size_t) Header->HeaderLength1;
	data->TagHdr[1] = BufferPtr;
	BufferPtr += (size_t) Header->HeaderLength2;
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Set quality scores
	data->TagQS[0] = BufferPtr;
	BufferPtr += (size_t) Header->TagLength1;
	data->TagQS[1] = BufferPtr;
	BufferPtr += (size_t) Header->TagLength2;
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Streamed Next
	data->Tags[0] = NULL;
	data->Tags[1] = NULL;
	if (Header->StreamedNext & 0b00000001) {
		data->Tags[0] = BufferPtr;
		BufferPtr += Header->TagLength1;
	}
	if (Header->StreamedNext & 0b00000010) {
		data->Tags[1] = BufferPtr;
		BufferPtr += Header->TagLength2;
	}
	char maskSoftClip = 0b00010000;
	for (int i=0; i<4; i++) {
		if ((data->AlnResults[i] != NULL) && (Header->StreamedNext & maskSoftClip)) {
			data->SoftClips[i] = BufferPtr;
			size_t SoftClipLength;
			if (!data->AlnResults[i]->revNmer) {
				assert(data->AlnResults[i]->Cigar[(size_t) data->AlnResults[i]->CigarLength-2] == 'S');
			 	SoftClipLength = (size_t) data->AlnResults[i]->Cigar[(size_t) data->AlnResults[i]->CigarLength-3];
			}
			else {
				assert(data->AlnResults[i]->Cigar[1] == 'S');
				SoftClipLength = (size_t) data->AlnResults[i]->Cigar[0];
			}
			data->SoftClipLengths[i] = SoftClipLength;	
			BufferPtr += SoftClipLength;
		}
		else {
			data->SoftClips[i] = NULL;
			data->SoftClipLengths[i] = 0UL;
		}
		maskSoftClip <<= 1;
	}

	//AlignAddressTo(BufferPtr, __alignof__(outputPairBlockHeader_t));
	return (uintptr_t) BufferPtr - Start;
}
//---------------------------------------------------------------

void PrintAlnResult(const AlignmentResults_t * AlnResult, const Genome_t * const genome,
                    const int taglen)
{
	char Buffer[kTagSize+2];
	const unsigned int genomepos = (genome->virtchr[AlnResult->chr].chr << 28) \
	                             + genome->virtchr[AlnResult->chr].offset + AlnResult->chr_pos;
#ifndef NDEBUG
	printf("Genome search range [%u-%u], anchor id %u ", AlnResult->GenomeSeekStart, AlnResult->GenomSeekStop, (unsigned int) AlnResult->AnchorID);								 
#endif
	printf("score=%d rev=%d MisMatch=%d chr=%d genomestart=%u pos=%d history=%d" ,
				 AlnResult->score, AlnResult->revNmer, AlnResult->MismatchCount, AlnResult->chr,
				 genomepos, AlnResult->chr_pos, AlnResult->DecisionTreeHistory);
	fputc('(',stdout);
	for (int m=13; m>=0; m--) {
		const int l = (AlnResult->DecisionTreeHistory & (1<<m) ) ? '1': '0';
		fputc(l, stdout);
		if (m % 4 == 0) fputc(' ',stdout);
	}
	fputs(")\n",stdout);
	
	printf("GEN: %.*s\n", taglen, genome->table + genomepos);
	if (AlnResult->revNmer) {
		printf("REV: ");
		int k=taglen;
		unsigned char * ptr = genome->table + genomepos;
		while (--k>=0) {
			int dna = (int) ptr[k];
			int p = '?';
			switch(dna) {
				case 'A': p = 'T'; break;
				case 'T': p = 'A'; break;
				case 'C': p = 'G'; break;
				case 'G': p = 'C'; break;
				case 'N': p = 'N'; break;
			}
			fputc(p, stdout);
		}
		fputc('\n', stdout);
	}
	int mpos = 0;
#ifdef EXPORT_STATES
	printf("STA: %s\n", AlnResult->States);
#endif
	const unsigned char * const restrict cigar = &(AlnResult->Cigar[0]);;
	printf("CIG: ");
	while (cigar[mpos] != kCIGARTERMINATOR) {
		printf(" %d%1c(%u)", (int) cigar[mpos], (int) cigar[mpos+1], cigar[mpos+1]);
		if (mpos >= kMAXcigarSize) {
			printf(" ERROR");
			break;
		}
		mpos += 2;
	}
	fputs(",\nMM : ", stdout);
	
	const unsigned char * const restrict MM = &(AlnResult->Mismatch[0]);;
	
	mpos = 0;
	while (MM[mpos] != kMMTERMINATOR) { 
		if (MM[mpos] == kMMSOFTCLIP) {
				printf(" SOFTCLIP ");
				mpos++;
				while(MM[mpos] != kMMTERMINATOR) fputc((int) MM[mpos++], stdout);
				break;
		}
		printf(" %d%1c", (int) MM[mpos], (int) MM[mpos+1]);
// 		assert(mpos <= (kMAXmismatchSize-1));
		if (mpos > (kMAXmismatchSize-1)) {
			printf(" ERROR");
			break;
		}
		mpos += 2;
	}
	fputs(",\n",stdout);
}
//---------------------------------------------------------------

void PrintStreamData(const StreamData_t * const data, const Genome_t * const genome) 
{
	const register outputPairBlockHeader_t * const Header = data->BlockHeader;
	printf("Block Header: alignments(%u), streamnext(0x%8.8x), ordinals(%lu,%lu)\n"
	       "       Tag 1: header length(%u), tag length(%hu)\n"
				 "              '%.*s'\n"
				 "              '%.*s'\n"
				 "       Tag 2: header length(%u), tag length(%hu)\n"
				 "              '%.*s'\n"
				 "              '%.*s'\n\n",
				 (unsigned int) Header->nAlignmentResults, Header->StreamedNext, Header->Ordinal1, Header->Ordinal2,
				 (unsigned int) Header->HeaderLength1, Header->TagLength1, 
				 (int) Header->HeaderLength1, data->TagHdr[0], (int) Header->TagLength1, data->TagQS[0],
				 (unsigned int) Header->HeaderLength2, Header->TagLength2, 
				 (int) Header->HeaderLength2, data->TagHdr[1], (int) Header->TagLength2, data->TagQS[1]);
	for (int i=0; i<4; i++) {
		if (data->SoftClips[i]) printf("SoftClip %i: %.*s\n", i, (int) data->SoftClipLengths[i], data->SoftClips[i]);

	}	
	if (data->AlnResults[0]) PrintAlnResult(data->AlnResults[0], genome, (int) Header->TagLength1);
	if (data->AlnResults[1]) PrintAlnResult(data->AlnResults[1], genome, (int) Header->TagLength2);
	if (data->AlnResults[2]) PrintAlnResult(data->AlnResults[2], genome, (int) Header->TagLength1);
	if (data->AlnResults[3]) PrintAlnResult(data->AlnResults[3], genome, (int) Header->TagLength2);
}
//---------------------------------------------------------------

/* vim: tabstop=2 shiftwidth=2
 */
