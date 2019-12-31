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
#define _GNU_SOURCE
#include <sched.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <fcntl.h>
#include <unistd.h>
#include <pthread.h>
#include <semaphore.h>
#include <assert.h>
#include <errno.h>
#include <signal.h>

#include "Aligner.h"
#include "Writer.h"
#include "statistics.h"
#include "GTL.h"
#include "ALIGNkonst.h"

#ifndef NDEBUG
#include <ctype.h>
#include "stream.h"
#endif

//--------------------------------------------------------------- 
// DEFINITIONS
//--------------------------------------------------------------- 
#define AlignAddressTo(ptr, req_aln) {\
	ptr = (char*) (((uintptr_t) (ptr) + (req_aln-1)) & ~(req_aln-1));\
}

//--------------------------------------------------------------- 
// STRUCTURE DEFINITIONS 
//--------------------------------------------------------------- 
typedef struct StringLength_s {
	unsigned char *stream;
	unsigned int length;
} StringLength_t;

typedef struct MatchData_s {
	StringLength_t cigar;
	StringLength_t mismatch;
	StringLength_t packedData;
} MatchData_t;

typedef struct MATCHES_struct {
	PAIRCHRPOSDATA pcpd[kMAXnbMatches];
	const char *hdr;
	const char *qs;
	unsigned long ordinal;
	unsigned int mml[kMAXnbMatches];
	unsigned int cl[kMAXnbMatches];
	unsigned int pl;
	unsigned int nbMatches;
	unsigned char MMstring[kMAXnbMatches][4*kTagSize+2];
	unsigned char cigar[kMAXnbMatches][2*kMAXcigarSize];
	unsigned char p4nt[kTagSize/2];
} MATCHES;

//--------------------------------------------------------------- 
// EXTERNS 
//---------------------------------------------------------------
extern struct timeval ts;
extern int verbose;

//--------------------------------------------------------------- 
// GLOBALS 
//--------------------------------------------------------------- 
int * restrict timesWritten = NULL;

//--------------------------------------------------------------- 
// LOCALS 
//--------------------------------------------------------------- 
static unsigned char gNT2bits[128];
static const VIRTUALCHR * restrict virtchr = NULL;
static unsigned int chrPartOffset[kMAXCHRcnt];
static unsigned int internalbuf;

//--------------------------------------------------------------- 
// LOCAL FCTS 
//---------------------------------------------------------------
static inline void PrintTag(const AlignmentResults_t * const restrict AlnResult, const unsigned char * const restrict tag )
{
	extern char genomefile[512];
	char Buffer[kTagSize+2];
	
	const unsigned int genomepos = (virtchr[AlnResult->chr].chr << 28) \
	                              + virtchr[AlnResult->chr].offset + AlnResult->chr_pos - 1;
	
	FILE * gf = fopen(genomefile, "r");
	if (gf == NULL) {
		fputs("Cannot open genome file\n", stdout);
		return;
	}
	
	const int taglen = (int) tag[kTagSize-1];
	fseek(gf, (long) genomepos, SEEK_SET);
	const size_t nread = fread(Buffer, sizeof(char), taglen+2, gf);
	fclose(gf);
	
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
	
	if (AlnResult->revNmer) {
		printf("TAG: %.*s\n", taglen, tag);
		printf("REV: ");
		int k=taglen;
		while (--k>=0) {
			int dna = (int) tag[k];
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
	else 
		printf("TAG: %.*s\n", taglen, tag);
	
	printf("GEN:%.*s\nCIG:", taglen+2, Buffer);
	int mpos = 0;
	
	const unsigned char * const restrict cigar = &(AlnResult->Cigar[0]);;
	
	while (cigar[mpos] != kCIGARTERMINATOR) {
		printf(" %d%1c(%u)", (int) cigar[mpos], (int) cigar[mpos+1], cigar[mpos+1]);
		if (mpos >= kMAXcigarSize) {
			printf(" ERROR");
			break;
		}
		mpos += 2;
	}
	fputs(",\nMM :", stdout);
	
	const unsigned char * const restrict MM = &(AlnResult->Mismatch[0]);;
	
	mpos = 0;
	while (MM[mpos] != kMMTERMINATOR) { 
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

static unsigned int
StatCount(const AlignmentResults_t * const restrict AlnResult1, const AlignmentResults_t * const restrict AlnResult2,
 				  const int maxTagScore1, const int maxTagScore2,
					const unsigned int DecisionTreeHistory_OR, const unsigned int DecisionTreeHistory_AND,
					COUNTS_TYPE_t * const CountPtr)
{
	unsigned int err = 0U;
	/* No Soft clip */
	if (!(DecisionTreeHistory_OR & kSoftClipRescued)) {
		/* No state change */
		if (AlnResult1->Cigar[0] == kCIGARTERMINATOR && AlnResult2->Cigar[0] == kCIGARTERMINATOR ) {
			const unsigned char MMstring1 = AlnResult1->Mismatch[0];
			const unsigned char MMstring2 = AlnResult2->Mismatch[0];
			/* Do we have mismatches ? */
			if ((AlnResult1->MismatchCount == 0) && (AlnResult2->MismatchCount == 0)) {
				if (AlnResult1->score == maxTagScore1 && AlnResult2->score == maxTagScore2) {
					if (DecisionTreeHistory_AND & kPerfectMatch) {
						CountPtr->pppf++;
						if (MMstring1 != kMMTERMINATOR) {
							printf("pppf error: ");
							err |= 0x1;
						}
						if (MMstring2 != kMMTERMINATOR) {
							printf("pppf error: ");
							err |= 0x2;
						}
					}
					/* At least 1 tag has some N */
					else {
						CountPtr->mNppf++;
						if ((MMstring1 == kMMTERMINATOR) && (MMstring2 == kMMTERMINATOR)) {
							printf("mNppf error should be pppf: ");
							err |= 0x3;
						}
					}
				}
				else {
					CountPtr->genomeN++;
				}
			}
			else {
				if ((MMstring1 == kMMTERMINATOR) && (MMstring2 == kMMTERMINATOR)) {
					printf("mppf error: ");
					err |= 0x3;
				}
				if (DecisionTreeHistory_OR & kContainsN) {
					CountPtr->mppfN++;
				}
				else {
					CountPtr->mppf++;
				}
			}
		}
		/* Not perfect match with state change*/
		else {
			if (DecisionTreeHistory_OR & kContainsN) {
				CountPtr->gppfN++;
			}
			else {
				CountPtr->gppf++;
			}
		}
	}
	/* This was rescued by Soft clip */
	else {
		/* No state change */
		if (AlnResult1->Cigar[0] == kCIGARTERMINATOR && AlnResult2->Cigar[0] == kCIGARTERMINATOR ) {
			if (DecisionTreeHistory_OR & kContainsN) {
				CountPtr->smppfN++;
			}
			else {
				CountPtr->smppf++;
			}
		}
		else {
			if (DecisionTreeHistory_OR & kContainsN) {
				CountPtr->sgppfN++;
			}
			else {
				CountPtr->sgppf++;
			}
		}
	}
	return err;
}
//---------------------------------------------------------------

static void 
ReGenerateTag(const Genome_t * const restrict Genome, const AlignmentResults_t * const restrict AlnResult,
              char * const restrict Tag, const char * restrict * SoftClips, const size_t TagLength)
{
	unsigned char tmpStr[kTagSize] __attribute__((aligned(16)));
	const size_t offset = (Genome->virtchr[AlnResult->chr].chr << 28) + Genome->virtchr[AlnResult->chr].offset + AlnResult->chr_pos;
	const unsigned char * restrict GenomeString = &(Genome->table[offset]);
	
	memcpy(Tag, GenomeString, TagLength);
	Tag[TagLength] = '\0';
	unsigned int cigarlength = 0U;
	const unsigned char * const restrict cigar = AlnResult->Cigar;
	int HasSoftClip = 0;
	if (cigar != NULL && *cigar != kCIGARTERMINATOR) {
		int r=0, a=0;
// 			printf("CIG: ");
		while (cigar[cigarlength] != kCIGARTERMINATOR)
		{
			register int pos = cigar[cigarlength];
			const char code  = cigar[cigarlength+1];
// 				printf("%u%c ", (unsigned int) cigar[cigarlength], (int) cigar[cigarlength+1]); 
			cigarlength += 2;

			switch(code)
			{
				case 'S': HasSoftClip = 1; // same as 'M'
				case 'M': assert(a+pos<=kTagSize); while(pos-- >0) { Tag[a++] = GenomeString[r]; r++; } break;
				case 'D':                          while(pos-- >0) { r++;                             } break;
				case 'I': assert(a+pos<=kTagSize); while(pos-- >0) { Tag[a++] = 'N';                  } break;
			}
		}
// 			printf("\n");
		Tag[a] = '\0';
		assert(r >= 1);
	}

	unsigned int count = 0U;
	const unsigned char * const restrict diff = AlnResult->Mismatch;
	if (diff) {
		while (diff[count] != kMMTERMINATOR) {
			const unsigned int pos = (unsigned int) diff[count];
			assert(pos<TagLength);
			const unsigned char nt = diff[count+1];
			Tag[pos] = nt;
			count += 2;
		}
// 		printf("\n");
	}
	
	if (HasSoftClip) {
		if (AlnResult->revNmer) {
			assert(cigar[1] == 'S');
			int cnt = (int) cigar[0];
			char * restrict dest = Tag;
			const char * restrict src = *SoftClips;
			while(cnt-- > 0) *dest++ = *src++;
			*SoftClips = src;
		}
		else {
			assert(cigar[cigarlength-1] == 'S');
			int cnt = (int) cigar[cigarlength-2];
			char * restrict dest = &Tag[TagLength-cnt];
			const char * restrict src = *SoftClips;
			while(cnt-- > 0) *dest++ = *src++;
			*SoftClips = src;
		}
	}
}
//---------------------------------------------------------------

static void
addMatch(FlushQueue *fq, TLBDATA **bp, MATCHES * const restrict mp, const outputPairBlockHeader_t * const Header)
{
	TLBDATA * const restrict block = *bp;
	assert(mp->nbMatches == 1 || mp->nbMatches == 2);
	// keep track if all reads have the same length
	if (block->cnt == 0) {
		block->lengthInfo = Header->TagLength1;
		block->nbMatches = mp->nbMatches;
	}
	if (block->nbMatches != mp->nbMatches) printf("Block made for %u matches\n", block->nbMatches);
	assert(block->nbMatches == mp->nbMatches);
//	fprintf(stderr, "%d %d %d %ld %ld\n", block->lengthInfo, Header->TagLength1, len2, strlen(mp->qs), strlen(qs2));
	//if (pcpdp != NULL)
	if (mp->pcpd[0].chr != kUNMAPPED) {
		if (mp->nbMatches > 1) {
			// for multiple matches, put the lowest first
			assert((block->header.flags_15_8 & kTLBHflagMultipleMatches) != 0);
			if (comparePAIRCHRPOSDATA(mp->pcpd, mp->pcpd + 1) > 0) {
				// we need to swap...
				PAIRCHRPOSDATA t_pcpd;
				unsigned char t_MMstring[4*kTagSize+2];
				unsigned char t_cigar[2*kMAXcigarSize];
				unsigned int tem;
				t_pcpd = mp->pcpd[0];
				memcpy(t_MMstring,mp->MMstring[0],mp->mml[0]);
				memcpy(t_cigar,mp->cigar[0],mp->cl[0]);
			 	mp->pcpd[0] = mp->pcpd[1];
				memcpy(mp->MMstring[0],mp->MMstring[1],mp->mml[1]);
				memcpy(mp->cigar[0],mp->cigar[1],mp->cl[1]);
				mp->pcpd[1] = t_pcpd;
				memcpy(mp->MMstring[1],t_MMstring,mp->mml[0]);
				memcpy(mp->cigar[1],t_cigar,mp->cl[0]);
				tem = mp->mml[0];
				mp->mml[0] = mp->mml[1];
				mp->mml[1] = tem;
				tem = mp->cl[0];
				mp->cl[0] = mp->cl[1];
				mp->cl[1] = tem;
			}
		}
		unsigned int minp;
		assert(((mp->pcpd[0].chr >> 16) & 0x7FFF) <= kMAXCHRcnt);
#ifndef NDEBUG
		if (mp->nbMatches > 1) {
			unsigned int t = (mp->pcpd[1].chr >> 16) & 0x7FFF;
			if (t >= kMAXCHRcnt) {
				printf("Got an alternate tag 1 chromosome number %u, ordinal %lu\n", t, mp->ordinal);
				assert(0);
			}
			t = mp->pcpd[1].chr & 0x7FFF;
			if (t >= kMAXCHRcnt) {
				printf("Got an alternate tag 2 chromosome number %u, ordinal %lu\n", t, mp->ordinal);
				assert(0);
			}
		}
#endif
		assert(((mp->pcpd[0].chr >> 16) & 0x7FFF) == block->header.chr || block->header.chr >= kMAXCHRcnt);
		if ((block->header.flags_7_0 & kTLBHflagPairedReads) == 0) {
			minp = mp->pcpd[0].tag1pos;
		}
		else {
			if (mp->pcpd[0].tag1pos != 0) {
				if (mp->pcpd[0].tag2pos != 0) {
					if (mp->pcpd[0].tag1pos <= mp->pcpd[0].tag2pos)
						minp = mp->pcpd[0].tag1pos;
					else
						minp = mp->pcpd[0].tag2pos;
				}
				else
					minp = mp->pcpd[0].tag1pos;
			}
			else
				minp = mp->pcpd[0].tag2pos;
		}
		if (block->header.minPos > minp) block->header.minPos = minp;
		if (block->header.maxPos < minp) block->header.maxPos = minp;
		for (unsigned int i = 0; i < mp->nbMatches; i++) block->pd[block->cnt].pcpd[i] = mp->pcpd[i];
	}
	else {
		block->header.minPos = 0;
		block->header.maxPos = 0;
	}
	assert(mp->hdr[Header->HeaderLength1 - 1] != '\n');
	block->pd[block->cnt].ordinal = mp->ordinal;
	block->pd[block->cnt].len1 = Header->TagLength1;
	block->pd[block->cnt].hlen1 = Header->HeaderLength1;
	memcpy(block->pd[block->cnt].hdr1,mp->hdr,Header->HeaderLength1);
	block->hdrSize += Header->HeaderLength1;
	memcpy(block->pd[block->cnt].qs1,mp->qs,Header->TagLength1);
	block->qualSize += Header->TagLength1;
	
	if ((block->header.flags_7_0 & kTLBHflagPairedReads) != 0) {
		const char *hdr2 = mp->hdr + Header->HeaderLength1;
		const char *qs2 = mp->qs + Header->TagLength1;
		if (Header->TagLength1 != block->lengthInfo || Header->TagLength2 != block->lengthInfo) block->lengthInfo = 0;
		assert(Header->HeaderLength2 == 0 || hdr2[Header->HeaderLength2 - 1] != '\n');
		block->pd[block->cnt].len2 = Header->TagLength2;
		block->pd[block->cnt].hlen2 = Header->HeaderLength2;
		memcpy(block->pd[block->cnt].hdr2,hdr2,Header->HeaderLength2);
		block->hdrSize += Header->HeaderLength2;
		memcpy(block->pd[block->cnt].qs2,qs2,Header->TagLength2);
		block->qualSize += Header->TagLength2;
	}
	else {
		if (Header->TagLength1 != block->lengthInfo) block->lengthInfo = 0;
	}
	for (unsigned int i = 0; i < mp->nbMatches; i++) {
		memcpy(block->pd[block->cnt].MMstring[i],mp->MMstring[i],mp->mml[i]);
		block->mmSize[i] += mp->mml[i];
		block->pd[block->cnt].ml[i] = mp->mml[i];
		memcpy(block->pd[block->cnt].cigar[i],mp->cigar[i],mp->cl[i]);
		block->cigarSize[i] += mp->cl[i];
		block->pd[block->cnt].cl[i] = mp->cl[i];
	}
	memcpy(block->pd[block->cnt].packed4nt,mp->p4nt,mp->pl);
	block->pd[block->cnt].pl = mp->pl;
	block->p4Size += mp->pl;
	block->cnt += 1;
	if (block->cnt >= kBlockMaxCount) flushBlock(fq,bp);
}
//---------------------------------------------------------------

static void
outputToPairedBuffer(GTLBUFFERS * const restrict G, char * const restrict Buffer,
                     COUNTS * const restrict cnts, const Genome_t * const restrict Genome, const size_t BufferSize)
{
	char GeneratedTags[2][kTagSize];
	const char * BufferPtr = Buffer;
	const uintptr_t BufferLimit = (uintptr_t) Buffer + BufferSize - 7; // 7 is for alignment padding
#ifndef NDEBUG
	const char * volatile CurrentBuffer;
	const char * volatile PastBuffer = Buffer;
#endif		
	while ((uintptr_t) BufferPtr < BufferLimit) {
#ifndef NDEBUG
		CurrentBuffer = BufferPtr;
#endif
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Define and initialize data to be fille
		MATCHES M;
		// FIXME - maybe can initialize only parts...
		// assume tagpos1 and tagpos2 are set to 0 below
		memset(&M,0,sizeof(MATCHES));
		
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Set the info header after aligning correctly
		//printf("Aligning results to %lu bytes, initially 0x%lx", __alignof__(outputPairBlockHeader_t), BufferPtr);
		//AlignAddressTo(BufferPtr, __alignof__(outputPairBlockHeader_t));
		//printf(", now 0x%lx\n", BufferPtr); fflush(stdout);
		if (memcmp(BufferPtr, "NEW_ALN\0", 8) != 0) {
			printf("Inconsistency @ buffer address 0x%lx\n", (uintptr_t) BufferPtr); fflush(stdout);
#ifndef NDEBUG
			printf("Error came from previous data located @ 0x%lx\n", (uintptr_t) PastBuffer); fflush(stdout);
			//AlignAddressTo(PastBuffer, __alignof__(outputPairBlockHeader_t));
			StreamData_t stream;
			const int PairBufferLength = GetStreamData(PastBuffer, &stream);
			PrintStreamData(&stream, Genome);
			PastBuffer += PairBufferLength;
			printf("\n\nLast pair used %i bytes of the buffer, thus providing this pair starting position @ 0x%lx\n",
			       PairBufferLength, (uintptr_t) PastBuffer);
			printf("data: ");
			for (int i=0; i<16; i++) {
				printf("%i ", (int) PastBuffer[i]);
				if (isprint((int) PastBuffer[i])) printf("(%c) ", (int) PastBuffer[i]);
			}
			fputc('\n', stdout);
			
			AlignAddressTo(PastBuffer, __alignof__(outputPairBlockHeader_t));
			printf("Access to next aligned address would be 0x%lx\ndata: ", (uintptr_t) PastBuffer);
			for (int i=0; i<16; i++) {
				printf("%i ", (int) PastBuffer[i]);
				if (isprint((int) PastBuffer[i])) printf("(%c) ", (int) PastBuffer[i]);
			}
			fputc('\n', stdout);
#endif
			assert(0);
			exit(1);
		}
		//printf("BufferPtr = 0x%lx\n", BufferPtr);
		BufferPtr += 8UL;
		const outputPairBlockHeader_t * const Header = (outputPairBlockHeader_t *) BufferPtr;
		//printf("BufferPtr = 0x%lx\n", BufferPtr);
		BufferPtr += sizeof(outputPairBlockHeader_t);	
	
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Set the alignment results
		const AlignmentResults_t * AlnResults[4] = { NULL, NULL, NULL, NULL };
		{
			register const int nAlignments = (int) Header->nAlignmentResults;
			for (int i=0; i<nAlignments; i++) {
				AlnResults[i] = (AlignmentResults_t *) BufferPtr;
				//printf("BufferPtr = 0x%lx\n", BufferPtr);
				BufferPtr += sizeof(AlignmentResults_t);
			}
		}
		//printf("BufferPtr = 0x%lx\n", BufferPtr);
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Set the headers
		M.hdr = BufferPtr;
		BufferPtr += (size_t) ((size_t) Header->HeaderLength1 + (size_t) Header->HeaderLength2);
		//printf("BufferPtr = 0x%lx\n", BufferPtr);
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Set quality scores
		M.qs = BufferPtr;
		BufferPtr += (size_t) ((size_t) Header->TagLength1 + (size_t) Header->TagLength2);
		//printf("Tag 1: %hu, Tag2: %hu , Header 1: %u, Header 2: %u\n", Header->TagLength1, Header->TagLength2, (unsigned int) Header->HeaderLength1, (unsigned int) Header->HeaderLength2);
		//printf("BufferPtr = 0x%lx\n", BufferPtr);
	
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Set the ordinal
		M.ordinal = Header->Ordinal1;
		
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Analyze the type of alignment here
		const char * Tags[2] = { NULL, NULL };
		unsigned int nNotMapped = 0U;
		if (Header->StreamedNext & 0b00000001) {
				Tags[0] = BufferPtr;
				BufferPtr += Header->TagLength1;
				//printf("BufferPtr = 0x%lx\n", BufferPtr);
				nNotMapped++;
		}
		if (Header->StreamedNext & 0b00000010) {
				Tags[1] = BufferPtr;
				BufferPtr += Header->TagLength2;
				//printf("BufferPtr = 0x%lx\n", BufferPtr);
				nNotMapped++;
		}
		
		const char * restrict SoftClips = (Header->StreamedNext & 0b11110000) ? BufferPtr : NULL;
		
		//printf("Next aln should be @ 0x%lx\n", BufferPtr);	
		{
			switch (nNotMapped) {
				case 0:
					// check if
					// - both reads are mapped on the same chr
					// - they are not too far apart i.e., let's go with 10 bits for now (1023), maybe revisit for RNAseq
					{
						assert(AlnResults[0] != NULL && AlnResults[1] != NULL);
						/*const*/ unsigned int DecisionTreeHistory_OR  = AlnResults[0]->DecisionTreeHistory | AlnResults[1]->DecisionTreeHistory;
						/*const*/ unsigned int DecisionTreeHistory_AND = AlnResults[0]->DecisionTreeHistory & AlnResults[1]->DecisionTreeHistory;
				
						const int maxTagScore1 = (_M + _MM) * (int) Header->TagLength1;
						const int maxTagScore2 = (_M + _MM) * (int) Header->TagLength2;

						assert(DecisionTreeHistory_AND & kSuccessfullyEncoded);
						COUNTS_TYPE_t * restrict CountPtr = NULL;
						if ((AlnResults[0]->chr == AlnResults[1]->chr) && (labs((long) AlnResults[0]->chr_pos - (long) AlnResults[1]->chr_pos) <= 0x3ff)) {
							CountPtr = &(cnts->Pair[AlnResults[0]->chr]);
						}
						// - both reads are mapped on the same chr but too far appart to be connected
						else if (AlnResults[0]->chr == AlnResults[1]->chr) {
							CountPtr = &(cnts->TooFar[AlnResults[0]->chr]);
						}
						// - not mapped on the same chromosome, chimera
						else {
							CountPtr = &(cnts->Chimera);
						}
						StatCount(AlnResults[0], AlnResults[1], maxTagScore1, maxTagScore2,
											DecisionTreeHistory_OR, DecisionTreeHistory_AND, CountPtr);
					}
					break;
				case 1:
					{
						// - half mapped / grosse magouille pour reutiliser le code de StatCount
						assert((AlnResults[0]->DecisionTreeHistory & kSuccessfullyEncoded));
						const int maxTagScore1 = (_M + _MM) * (int) Header->TagLength1;
						StatCount(AlnResults[0], AlnResults[0], maxTagScore1, maxTagScore1,	AlnResults[0]->DecisionTreeHistory,
											AlnResults[0]->DecisionTreeHistory, &(cnts->HalfMap[AlnResults[0]->chr]));
					}
					break;
				case 2:
					cnts->upf++;
					break;
				default: /* should be none here */
					printf("THERE IS AN ERROR!!!\n");
					//PrintTag(&(TagResult1->Best), Tag1);
					//PrintTag(&(TagResult2->Best), Tag2);
					fflush(stdout);
					exit(1);
			}
				
//			if (verbose & 0x4) {
//				printf("\nBest 1: "); PrintTag(&(TagResult1->Best), Tag1);
//				printf("Best 2: "); PrintTag(&(TagResult2->Best), Tag2);
//				if (TagResult1->Alternate.chr != kUNMAPPED) {
//					printf("Alternate 1: "); PrintTag(&(TagResult1->Alternate), Tag1);
//					printf("Alternate 2: "); PrintTag(&(TagResult2->Alternate), Tag2);
//				}
//			}
			
			/*for (unsigned int i=0; i<nNotMapped; i++) {
				if (AlnResults[i] && AlnResults[i]->chr > 0x3FF) {
 					PrintTag(AlnResults[i], Tag1);
				}
			}*/
		}
		
		M.pcpd[0].chr = kUNMAPPED;
		M.nbMatches = 1;
		// M.nbMatches = 0; performed by memset at the beginning

		// check if
		// - both reads are mapped on the same chr
		// - they are not too far apart i.e., let's go with 10 bits for now (1023), maybe revisit for RNAseq
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// MAPPED ON THE SAME CHROMOSOME AT DISTANCE BELOW 1024
		if ((AlnResults[0] && AlnResults[1]) && (AlnResults[0]->chr == AlnResults[1]->chr)
		 && (abs((int)AlnResults[0]->chr_pos - (int)AlnResults[1]->chr_pos) <= 0x3ff))
		{
			unsigned int chr = AlnResults[0]->chr;
			M.pcpd[0].tag1pos = AlnResults[0]->chr_pos;
			M.pcpd[0].tag2pos = AlnResults[1]->chr_pos;
			unsigned int blockId = AlnResults[0]->chr_pos; //M.pcpd[0].tag1pos;
			if (M.pcpd[0].tag2pos < blockId) blockId = M.pcpd[0].tag2pos;
			blockId /= kCountsPerChrPart;
			blockId += chrPartOffset[chr];
			assert(blockId < G->OUTPUTFILEScnt);			
			
			M.pcpd[0].chr = (chr << 16) | chr;
			if (AlnResults[0]->revNmer) M.pcpd[0].chr |= 0x80000000;
			if (AlnResults[1]->revNmer) M.pcpd[0].chr |= 0x00008000;

			if ((AlnResults[2] && AlnResults[3]) && (AlnResults[2]->chr == AlnResults[3]->chr)
			 && (abs((int)AlnResults[2]->chr_pos - (int)AlnResults[3]->chr_pos) <= 0x3ff))
			{
				unsigned int chr2 = AlnResults[2]->chr;
				unsigned int blockId2;
				M.nbMatches = 2;
				M.pcpd[1].tag1pos = AlnResults[2]->chr_pos;
				blockId2 = M.pcpd[1].tag1pos;
				M.pcpd[1].tag2pos = AlnResults[3]->chr_pos;
				if (M.pcpd[1].tag2pos < blockId2)
					blockId2 = M.pcpd[1].tag2pos;
				blockId2 /= kCountsPerChrPart;
				blockId2 += chrPartOffset[chr2];
				M.pcpd[1].chr = chr2 << 16 | chr2;
				if (AlnResults[2]->revNmer) M.pcpd[1].chr |= 0x80000000;
				if (AlnResults[3]->revNmer) M.pcpd[1].chr |= 0x00008000;
				if (blockId2 < blockId) {
					chr = chr2;
					blockId = blockId2;
				}
			}
// 			else
// 				extrslts[extrsltstart+i].Alternate.chr = kUNMAPPED;

			// !!!! should also test that length is exactly the same as 'normal' one as it is not stored.
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// PERFECT MATCHES
			if (   (AlnResults[0]->Mismatch[0] == kMMTERMINATOR) && (AlnResults[1]->Mismatch[0] == kMMTERMINATOR)
			    && (AlnResults[0]->Cigar[0] == kCIGARTERMINATOR) && (AlnResults[1]->Cigar[0] == kCIGARTERMINATOR))
			{
				TLBDATA** ptr;
				if (M.nbMatches == 2) {
					assert((AlnResults[2]->Mismatch[0] == kMMTERMINATOR) && (AlnResults[3]->Mismatch[0] == kMMTERMINATOR)
			        && (AlnResults[2]->Cigar[0] == kCIGARTERMINATOR) && (AlnResults[3]->Cigar[0] == kCIGARTERMINATOR));
					ptr = G->MperfectBlock + blockId;
				}
				else
					ptr = G->perfectBlock + blockId;
				
				addMatch(&G->fq, ptr, &M, Header);
			}
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// IMPERFECT MATCHES (with mismatches, but no alignment (no CIGAR))
			else if ( (AlnResults[0]->Cigar[0] == kCIGARTERMINATOR) && (AlnResults[1]->Cigar[0] == kCIGARTERMINATOR)
			       && (AlnResults[2] == NULL || ( (AlnResults[2]->Cigar[0] == kCIGARTERMINATOR) && (AlnResults[3]->Cigar[0] == kCIGARTERMINATOR) ) ) )
			{
				unsigned char MMstringAlwaysN[kMAXnbMatches][2*((2*kMaxEncodedMismatch)+1)]; // need 2x the number of mismatches (for the pos + nt) 
				unsigned int lN[kMAXnbMatches] = {0};
				int mmpos = 0;
				int MMisAlwaysN = 1;

				// check if 1st tag has mismatches
				while (AlnResults[0]->Mismatch[mmpos] != kMMTERMINATOR)
				{
					MMstringAlwaysN[0][lN[0]++] = AlnResults[0]->Mismatch[mmpos];		// store only position.
					M.MMstring[0][M.mml[0]++] = AlnResults[0]->Mismatch[mmpos++];
					M.MMstring[0][M.mml[0]++] = AlnResults[0]->Mismatch[mmpos];
					if (AlnResults[0]->Mismatch[mmpos++] != 'N')
						MMisAlwaysN = 0;
					if (mmpos == (2*kMaxEncodedMismatch)) {
						assert(AlnResults[0]->Mismatch[mmpos] == kMMTERMINATOR);
						break;
					}
				}
				M.MMstring[0][M.mml[0]++] = kMMTERMINATOR;
				MMstringAlwaysN[0][lN[0]++] = kMMTERMINATOR;

				// check if 2nd tag has mismatches
				mmpos = 0;
				while (AlnResults[1]->Mismatch[mmpos] != kMMTERMINATOR)
				{
					MMstringAlwaysN[0][lN[0]++] = AlnResults[1]->Mismatch[mmpos];		// store only position.
					M.MMstring[0][M.mml[0]++] = AlnResults[1]->Mismatch[mmpos++];
					M.MMstring[0][M.mml[0]++] = AlnResults[1]->Mismatch[mmpos];						
					if (AlnResults[1]->Mismatch[mmpos++] != 'N')
						MMisAlwaysN = 0;
					if (mmpos == (2*kMaxEncodedMismatch)) {
						assert(AlnResults[1]->Mismatch[mmpos] == kMMTERMINATOR);
						break;
					}
				} 
				M.MMstring[0][M.mml[0]++] = kMMTERMINATOR;
				MMstringAlwaysN[0][lN[0]++] = kMMTERMINATOR;

				if (AlnResults[2]) {
					// handle 2nd match
					// check if 1st tag has mismatches
					mmpos = 0;
					while (AlnResults[2]->Mismatch[mmpos] != kMMTERMINATOR) 
					{
						MMstringAlwaysN[1][lN[1]++] = AlnResults[2]->Mismatch[mmpos];		// store only position.
						M.MMstring[1][M.mml[1]++] = AlnResults[2]->Mismatch[mmpos++];
						M.MMstring[1][M.mml[1]++] = AlnResults[2]->Mismatch[mmpos];
						if (AlnResults[2]->Mismatch[mmpos++] != 'N')
							MMisAlwaysN = 0;
						if (mmpos == (2*kMaxEncodedMismatch)) {
							assert(AlnResults[2]->Mismatch[mmpos] == kMMTERMINATOR);
							break;
						}
					}
					M.MMstring[1][M.mml[1]++] = kMMTERMINATOR;
					MMstringAlwaysN[1][lN[1]++] = kMMTERMINATOR;

					// check if 2nd tag has mismatches
					mmpos = 0;
					while (AlnResults[3]->Mismatch[mmpos] != kMMTERMINATOR)
					{
						MMstringAlwaysN[1][lN[1]++] = AlnResults[3]->Mismatch[mmpos];		// store only position.
						M.MMstring[1][M.mml[1]++] = AlnResults[3]->Mismatch[mmpos++];
						M.MMstring[1][M.mml[1]++] = AlnResults[3]->Mismatch[mmpos];						
						if (AlnResults[3]->Mismatch[mmpos++] != 'N')
							MMisAlwaysN = 0;
						if (mmpos == (2*kMaxEncodedMismatch)) {
							assert(AlnResults[3]->Mismatch[mmpos] == kMMTERMINATOR);
							break;
						}
					} 
					M.MMstring[1][M.mml[1]++] = kMMTERMINATOR;
					MMstringAlwaysN[1][lN[1]++] = kMMTERMINATOR;
				}

				if (MMisAlwaysN) {
					memcpy(M.MMstring[0],MMstringAlwaysN[0],lN[0]);
					M.mml[0] = lN[0];
					if (AlnResults[2]) {
						memcpy(M.MMstring[1],MMstringAlwaysN[1],lN[1]);
						M.mml[1] = lN[1];
						addMatch(&G->fq, G->MmismatchNBlock + blockId, &M, Header);
					}
					else
						addMatch(&G->fq, G->mismatchNBlock + blockId, &M, Header);
				}
				else {
					if (AlnResults[2])
						addMatch(&G->fq, G->MmismatchBlock + blockId, &M, Header);
					else
						addMatch(&G->fq, G->mismatchBlock + blockId, &M, Header);
				}
			}
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// CIGAR WITH OR WITHOUT MICMATCHES
			else
			{
				// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
				// check if 1st tag has mismatches
				int mmpos = 0;
				while (AlnResults[0]->Mismatch[mmpos] != kMMTERMINATOR)
				{
					M.MMstring[0][M.mml[0]++] = AlnResults[0]->Mismatch[mmpos++];
					M.MMstring[0][M.mml[0]++] = AlnResults[0]->Mismatch[mmpos++];
					if (mmpos == (2*kMaxEncodedMismatch)) {
						assert(AlnResults[0]->Mismatch[mmpos] == kMMTERMINATOR);
						break;
					}
				}
				M.MMstring[0][M.mml[0]++] = kMMTERMINATOR;

				// check if 1st tag has cigar
				mmpos = 0;
				while (AlnResults[0]->Cigar[mmpos] != kCIGARTERMINATOR)
				{
					if (AlnResults[0]->Cigar[mmpos+1] == 'S') // soft clip detected. must add corresponding nucleotides in the MMstring.
					{
						assert(Header->StreamedNext & 0b00010000);
						M.MMstring[0][M.mml[0]-1] = kMMSOFTCLIP;
						for (int k=0; k<AlnResults[0]->Cigar[mmpos]; k++) { M.MMstring[0][M.mml[0]++] = *BufferPtr; BufferPtr++; } 
						M.MMstring[0][M.mml[0]++] = kMMTERMINATOR;
					}
					M.cigar[0][M.cl[0]++] = AlnResults[0]->Cigar[mmpos++];
					M.cigar[0][M.cl[0]++] = AlnResults[0]->Cigar[mmpos++];
					if (mmpos > 2*kMaxEncodedCigar) {
						printf("********************** OVERFLOW **********************\n"); 
// 						PrintTag(&(extrslts[extrsltstart+i].Best), &tags[tagbufstart+i*kTagSize]);
// 						PrintTag(&(extrslts[extrsltstart+i+1].Best), &tags[tagbufstart+i*kTagSize+kTagSize]);
						fflush(stdout);
						assert(0==1);
						break;
					}
				}
				M.cigar[0][M.cl[0]++] = kCIGARTERMINATOR;

				// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
				// check if 2nd tag has mismatches
				mmpos = 0;
				while (AlnResults[1]->Mismatch[mmpos] != kMMTERMINATOR)
				{
					M.MMstring[0][M.mml[0]++] = AlnResults[1]->Mismatch[mmpos++];
					M.MMstring[0][M.mml[0]++] = AlnResults[1]->Mismatch[mmpos++];						
					if (mmpos == (2*kMaxEncodedMismatch)) {
						assert(AlnResults[1]->Mismatch[mmpos] == kMMTERMINATOR);
						break;
					}
				} 
				M.MMstring[0][M.mml[0]++] = kMMTERMINATOR;

				// check if 2nd tag has cigar
				mmpos = 0;
				while (AlnResults[1]->Cigar[mmpos] != kCIGARTERMINATOR)
				{
					if (AlnResults[1]->Cigar[mmpos+1] == 'S') // soft clip detected. must add corresponding nucleotides in the MMstring.
					{
						//printf("tag2 softclip of len %d\n",extrslts[extrsltstart+i+1].cigar[mmpos]);
						int taglen = Header->TagLength2;
						int k;
						assert(Header->StreamedNext & 0b00100000);
						M.MMstring[0][M.mml[0]-1] = kMMSOFTCLIP;
						for (int k=0; k<AlnResults[1]->Cigar[mmpos]; k++) { M.MMstring[0][M.mml[0]++] = *BufferPtr; BufferPtr++; } 
						M.MMstring[0][M.mml[0]++] = kMMTERMINATOR;
					}
					M.cigar[0][M.cl[0]++] = AlnResults[1]->Cigar[mmpos++];
					M.cigar[0][M.cl[0]++] = AlnResults[1]->Cigar[mmpos++];
					if (mmpos > 2*kMaxEncodedCigar) {
						printf("********************** OVERFLOW **********************\n"); 
// 						PrintTag(&(extrslts[extrsltstart+i].Best), &tags[tagbufstart+i*kTagSize]);
// 						PrintTag(&(extrslts[extrsltstart+i+1].Best), &tags[tagbufstart+i*kTagSize+kTagSize]);
						fflush(stdout);
						assert(0==1);
						break;
					}
				} 
				M.cigar[0][M.cl[0]++] = kCIGARTERMINATOR;

				if (AlnResults[2]) {
					// handle 2nd match
					// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
					// check if 1st tag has mismatches
					mmpos = 0;
					while (AlnResults[2]->Mismatch[mmpos] != kMMTERMINATOR) {
						M.MMstring[1][M.mml[1]++] = AlnResults[2]->Mismatch[mmpos++];
						M.MMstring[1][M.mml[1]++] = AlnResults[2]->Mismatch[mmpos++];
						if (mmpos == (2*kMaxEncodedMismatch)) {
							assert(AlnResults[2]->Mismatch[mmpos] == kMMTERMINATOR);
							break;
						}
					}
					M.MMstring[1][M.mml[1]++] = kMMTERMINATOR;

					// check if 1st tag has cigar
					mmpos = 0;
					while (AlnResults[2]->Cigar[mmpos] != kCIGARTERMINATOR) {
						if (AlnResults[2]->Cigar[mmpos+1] == 'S') // soft clip detected. must add corresponding nucleotides in the MMstring.
						{
							assert(Header->StreamedNext & 0b01000000);
							M.MMstring[1][M.mml[1]-1] = kMMSOFTCLIP;
							for (int k=0; k<AlnResults[2]->Cigar[mmpos]; k++) { M.MMstring[1][M.mml[1]++] = *BufferPtr; BufferPtr++; } 
							M.MMstring[1][M.mml[1]++] = kMMTERMINATOR;
						}
						M.cigar[1][M.cl[1]++] = AlnResults[2]->Cigar[mmpos++];
						M.cigar[1][M.cl[1]++] = AlnResults[2]->Cigar[mmpos++];
						if (mmpos > 2*kMaxEncodedCigar) {
							assert(0==1);
							break;
						}
					}
					M.cigar[1][M.cl[1]++] = kCIGARTERMINATOR;

					// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
					// check if 2nd tag has mismatches
					mmpos = 0;
					while (AlnResults[3]->Mismatch[mmpos] != kMMTERMINATOR) {
						M.MMstring[1][M.mml[1]++] = AlnResults[3]->Mismatch[mmpos++];
						M.MMstring[1][M.mml[1]++] = AlnResults[3]->Mismatch[mmpos++];						
						if (mmpos == (2*kMaxEncodedMismatch)) {
							assert(AlnResults[3]->Mismatch[mmpos] == kMMTERMINATOR);
							break;
						}
					} 
					M.MMstring[1][M.mml[1]++] = kMMTERMINATOR;

					// check if 2nd tag has cigar
					mmpos = 0;
					while (AlnResults[3]->Cigar[mmpos] != kCIGARTERMINATOR) {
						if (AlnResults[3]->Cigar[mmpos+1] == 'S') // soft clip detected. must add corresponding nucleotides in the MMstring.
						{
							//printf("tag2 softclip of len %d\n",extrslts[extrsltstart+i+1].cigar[mmpos]);
							assert(Header->StreamedNext & 0b10000000);
							M.MMstring[1][M.mml[1]-1] = kMMSOFTCLIP;
							for (int k=0; k<AlnResults[3]->Cigar[mmpos]; k++) { M.MMstring[1][M.mml[1]++] = *BufferPtr; BufferPtr++; } 
							M.MMstring[1][M.mml[1]++] = kMMTERMINATOR;
						}
						M.cigar[1][M.cl[1]++] = AlnResults[3]->Cigar[mmpos++];
						M.cigar[1][M.cl[1]++] = AlnResults[3]->Cigar[mmpos++];						
						if (mmpos > 2*kMaxEncodedCigar) {
							assert(0 == 1);
							break;
						}
					} 
					M.cigar[1][M.cl[1]++] = kCIGARTERMINATOR;
				}

#ifdef DEBUG_SOFTCLIP_PROCESS
				printf("CIGAR= ");
				for (mmpos = 0; mmpos < M.cl[0]; mmpos++)
					printf("%d:%c|",(int)M.cigar[0][mmpos],M.cigar[0][mmpos]);

				printf(".  MM= ");
				mmpos = 0;
				for (mmpos = 0; mmpos < M.mml[0]; mmpos++)
					printf("%d:%c|",(int)M.MMstring[0][mmpos],M.MMstring[0][mmpos]);
				printf(".\n");
#endif


// 				cnts->gppf[chr]++;
				if (M.nbMatches == 2)
					addMatch(&G->fq, G->MalignBlock + blockId, &M, Header);
				else
					addMatch(&G->fq, G->alignBlock + blockId, &M, Header);
			}
		}
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// MAPPED TOO FAR AWAY
		/*else if ((AlnResults[0] && AlnResults[1]) && (AlnResults[0]->chr == AlnResults[1]->chr)) {
			//printf("Writer got a paired spaced by %i\n", abs((int)AlnResults[0]->chr_pos - (int)AlnResults[1]->chr_pos));
			for (int i=0; i<(int)Header->nAlignmentResults; i++) {
				if (AlnResults[i]->Cigar[0] != kCIGARTERMINATOR) {
					unsigned int ExtraSpace = 0;
					if (AlnResults[i]->revNmer != 0) {
						if (AlnResults[i]->Cigar[1] == 'S') ExtraSpace = (unsigned int) AlnResults[i]->Cigar[0]; 	
					}
					else {
						if (AlnResults[i]->Cigar[AlnResults[i]->CigarLength-2] == 'S') 
							ExtraSpace = (unsigned int) AlnResults[i]->Cigar[AlnResults[i]->CigarLength-3];
					}
					BufferPtr += ExtraSpace;
				}
			}
			goto TooFarJoin;
		}*/
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// MAP TOO FAR OR NOT MAPPED ON SAME CHROMOSOME OR UNMAPPED.
		else
		{	
			unsigned char *packed4ntPtr = M.p4nt;
			const AlignmentResults_t * restrict AlnTag0, * restrict AlnTag1;
			const unsigned int variation = (unsigned int) (Header->StreamedNext & 0b00000011);
			
			switch (variation) {
				case 0b00:
					AlnTag0 = AlnResults[0];
					AlnTag1 = AlnResults[1];
					break;
				case 0b01:
					AlnTag0 = NULL;
					AlnTag1 = AlnResults[0];
					break;
				case 0b10:
					AlnTag0 = AlnResults[0];
					AlnTag1 = NULL;
					break;
				case 0b11:
					AlnTag0 = NULL;
					AlnTag1 = NULL;
					break;
			}
						
			// handle gapped alignments as unmapped at this stage
			if (AlnTag0 != NULL && AlnTag0->Cigar[0] != kCIGARTERMINATOR) {
				ReGenerateTag(Genome, AlnTag0, GeneratedTags[0], /*&SoftClips*/ &BufferPtr, Header->TagLength1);
				assert(Tags[0] == NULL);
				Tags[0] = GeneratedTags[0];
				AlnTag0 = NULL;
			}
			
			if (AlnTag1 != NULL && AlnTag1->Cigar[0] != kCIGARTERMINATOR) {
				ReGenerateTag(Genome, AlnTag1, GeneratedTags[1], /*&SoftClips*/ &BufferPtr, Header->TagLength2);
				assert(Tags[1] == NULL);
				Tags[1] = GeneratedTags[1];
				AlnTag1 = NULL;
			}
			
			if (AlnTag0 == NULL) {
				assert(Tags[0]);
				unsigned int shift = 8;
				const char * restrict tp = Tags[0];
				unsigned int pos = 0;
				*packed4ntPtr = 0;
				const unsigned int taglen = Header->TagLength1;
				while(pos < taglen) {
					if (shift == 0) {
						shift = 8;
						packed4ntPtr += 1;
						*packed4ntPtr = 0;
						M.pl += 1;
					}
					shift -= 2;
					if (gNT2bits[*tp] == kNOT_ACGT) {
						M.MMstring[0][M.mml[0]++] = pos;
						M.MMstring[0][M.mml[0]++] = *tp;
					}	
					else
						*packed4ntPtr |= (gNT2bits[*tp] << shift);
					tp++;
					pos++;
				}
				M.pl += 1;
				packed4ntPtr += 1;
				shift = 8;
				M.MMstring[0][M.mml[0]++] = kMMTERMINATOR;
			}
			else {
				// check if 1st tag has mismatches
				int mmpos = 0;
				while (AlnTag0->Mismatch[mmpos] != kMMTERMINATOR) {
					M.MMstring[0][M.mml[0]++] = AlnTag0->Mismatch[mmpos++];
					M.MMstring[0][M.mml[0]++] = AlnTag0->Mismatch[mmpos++];
				}
				M.MMstring[0][M.mml[0]++] = kMMTERMINATOR;
			}

			if (AlnTag1 == NULL) {
				assert(Tags[1]);
				unsigned int shift = 8;
				const char * restrict tp = Tags[1];
				const unsigned int taglen = Header->TagLength2;
				unsigned int pos = 0;
				*packed4ntPtr = 0;
				while(pos < taglen)
				{
					if (shift == 0)	{
						shift = 8;
						packed4ntPtr += 1;
						*packed4ntPtr = 0;
						M.pl += 1;
					}
					shift -= 2;
					if (gNT2bits[*tp] == kNOT_ACGT) {
						M.MMstring[0][M.mml[0]++] = pos;
						M.MMstring[0][M.mml[0]++] = *tp;
					}	
					else				
						*packed4ntPtr |= (gNT2bits[*tp] << shift);
					tp++;
					pos++;
				}
				M.pl += 1;
				M.MMstring[0][M.mml[0]++] = kMMTERMINATOR;
			}
			else {
				// check if 2nd tag has mismatches
				int mmpos = 0;
				while (AlnTag1->Mismatch[mmpos] != kMMTERMINATOR) {
					M.MMstring[0][M.mml[0]++] = AlnTag1->Mismatch[mmpos++];
					M.MMstring[0][M.mml[0]++] = AlnTag1->Mismatch[mmpos++];						
				} 
				M.MMstring[0][M.mml[0]++] = kMMTERMINATOR;
			}

			if ((AlnResults[2] != NULL) && (Header->StreamedNext & 0b01000000)) {
				assert(AlnResults[2]->Cigar[0] != kCIGARTERMINATOR);
				unsigned int ExtraSpace = 0U;
				if (AlnResults[2]->revNmer != 0) {
					assert(AlnResults[2]->Cigar[1] == 'S');
					ExtraSpace = (unsigned int) AlnResults[2]->Cigar[0]; 	
				}
				else {
					assert(AlnResults[2]->Cigar[AlnResults[2]->CigarLength-2] == 'S'); 
					ExtraSpace = (unsigned int) AlnResults[2]->Cigar[AlnResults[2]->CigarLength-3];
				}
				BufferPtr += ExtraSpace;
			}

			if ((AlnResults[3] != NULL) && (Header->StreamedNext & 0b10000000)) {
                                assert(AlnResults[3]->Cigar[0] != kCIGARTERMINATOR);
                                unsigned int ExtraSpace = 0U;
                                if (AlnResults[3]->revNmer != 0) {
                                        assert(AlnResults[3]->Cigar[1] == 'S'); 
                                        ExtraSpace = (unsigned int) AlnResults[3]->Cigar[0]; 
                                }       
                                else {  
                                        assert(AlnResults[3]->Cigar[AlnResults[3]->CigarLength-2] == 'S');
                                        ExtraSpace = (unsigned int) AlnResults[3]->Cigar[AlnResults[3]->CigarLength-3];
                                }       
                                BufferPtr += ExtraSpace;
                        }

			//////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// check if CHIMERA
			if (AlnTag0 && AlnTag1) {
				assert(M.pl == 0);
				unsigned int blockId = AlnTag0->chr - 1;
				if (AlnTag0->chr != AlnTag1->chr) {
					unsigned int diff = abs((int) AlnTag0->chr - (int) AlnTag1->chr);
					if (diff > 8) diff = 8;
					blockId = chrMax - 1 + diff;
				}
				M.pcpd[0].chr = (AlnTag0->chr << 16) + AlnTag1->chr;
				if (AlnTag0->revNmer) M.pcpd[0].chr |= 0x80000000;
				if (AlnTag1->revNmer) M.pcpd[0].chr |= 0x00008000;
				// FIXME : a position of 0 can happen, e.g. for mitochondria...  not sure how to properly handle at this point
				//assert(extrslts[extrsltstart+i].pos != 0 && extrslts[extrsltstart+i+1].pos != 0);
				M.pcpd[0].tag1pos = AlnTag0->chr_pos;
				M.pcpd[0].tag2pos = AlnTag1->chr_pos;
				
				addMatch(&G->fq, G->chimeraBlock + blockId, &M, Header);
			}
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// check if completely unmapped
			else if ((AlnTag0 == NULL) && (AlnTag1 == NULL))
			{
				assert(Tags[0]); assert(Tags[1]);
				unsigned int blockId = ((gNT2bits[ Tags[0][0] ] & 0x3) << 2) | (gNT2bits[ Tags[1][0] ] & 0x3);
				addMatch(&G->fq, G->unmappedBlock + blockId, &M, Header);
			}
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// half mapped.
			else
			{
				assert((AlnTag0 == NULL && AlnTag1 != NULL) || (AlnTag0 != NULL && AlnTag1 == NULL));
				unsigned int blockId;
				if (AlnTag0) {
					M.pcpd[0].chr = (AlnTag0->chr << 16);
					if ( AlnTag0->revNmer) M.pcpd[0].chr |= 0x80000000;
					M.pcpd[0].tag1pos =  AlnTag0->chr_pos;
					assert(AlnTag0->chr > 0);
					blockId =  AlnTag0->chr - 1U;
				}
				else {
					M.pcpd[0].chr = AlnTag1->chr;
					if ( AlnTag1->revNmer) M.pcpd[0].chr |= 0x00008000;
					M.pcpd[0].tag2pos =  AlnTag1->chr_pos;
					assert(AlnTag1->chr > 0);
					blockId =  AlnTag1->chr - 1U;
				}
				assert(M.nbMatches == 1);
				assert(blockId < chrMax);
				addMatch(&G->fq, G->halfmapBlock + blockId, &M, Header);
			}
		}
		AlignAddressTo(BufferPtr, __alignof__(outputPairBlockHeader_t));
#ifndef NDEBUG
		PastBuffer = CurrentBuffer;
#endif
	}
}
//---------------------------------------------------------------

static void initPackNT2bitsTable()
{
	memset(gNT2bits,kNOT_ACGT,128);
	gNT2bits['A'] = 0x0;
	gNT2bits['C'] = 0x1;
	gNT2bits['G'] = 0x2;
	gNT2bits['T'] = 0x3;
} 
//---------------------------------------------------------------

//--------------------------------------------------------------- 
// FUNCTIONS 
//---------------------------------------------------------------

//---------------------------------------------------------------
// This thread is doing and infinite loop
// listening for blocks to flush to the disk
//
// it stops when it receives a block with 0 entries
//
void *flushBlockThread(GTLBUFFERS * const restrict G)
{
	TLBDATA wb;
	allocBlock(&wb, true);
	// --------------------  infinite listening loop ---------------------------------------
	while (1)
	{
		// receive next block to flush
		TLBDATA * const restrict block = get_block(&G->fq);
		if (block == NULL) 	return (void *) -1;
		if (block->cnt == 0) return 0;

		// do the job, sort the stuff according to position
		for (unsigned int i = 0; i < block->cnt; i++) block->pdp[i] = block->pd + i;
		
		// only actually sort if we have position info
		if ((block->header.flags_7_0 & kTLBHflagPairPosBlock) != 0)
			qsort(block->pdp,block->cnt,sizeof(PAIRDATA *),comparePAIRDATAptr);
		else if ((block->header.flags_7_0 & kTLBHflagPairChrPosBlock) != 0)
			qsort(block->pdp,block->cnt,sizeof(PAIRDATA *),comparePAIRUDATAptr);
		else if ((block->header.flags_15_8 & kTLBHflagUnmappedBlock) != 0)
			qsort(block->pdp,block->cnt,sizeof(PAIRDATA *),compareUNMAPDATAptr);
		unsigned char *dataBlock = packBlock(block, &wb);

		char fname[1024];
		if ((wb.header.chr & 0xffff) == kHalfMapChrID) {
			unsigned int chr = wb.header.chr >> 16;
			snprintf(fname,1024,"%s_HM_chr%d.gtl",G->seedfn,chr);
		}
		else if ((wb.header.chr & 0xffff) == kUnmappedChrID) {
			if ((wb.header.flags_7_0 & kTLBHflagPairedReads) == 0) {
				unsigned int c = (wb.header.chr >> 16) & 3;
				snprintf(fname,1024,"%s_UM_%c.gtl",G->seedfn,"ACGT"[c]);
			}
			else {
				unsigned int c1 = (wb.header.chr >> 18) & 3;
				unsigned int c2 = (wb.header.chr >> 16) & 3;
				snprintf(fname,1024,"%s_UM_%c%c.gtl",G->seedfn,"ACGT"[c1],"ACGT"[c2]);
			}
		}
		else if ((wb.header.chr & 0xffff) == kChimeraChrID) {
			unsigned int chr = wb.header.chr >> 16;
			if (chr > chrMax)
				snprintf(fname,1024,"%s_C_d%d.gtl",G->seedfn,chr - chrMax);
			else
				snprintf(fname,1024,"%s_C_chr%d.gtl",G->seedfn,chr);
		}
		else {
			if ((block->header.flags_15_8 & kTLBHflagMultipleMatches) != 0)
				snprintf(fname,1024,"%s_M_chr%d%c.gtl",G->seedfn,wb.header.chr,'a' + wb.header.minPos / kCountsPerChrPart);
			else
				snprintf(fname,1024,"%s_chr%d%c.gtl",G->seedfn,wb.header.chr,'a' + wb.header.minPos / kCountsPerChrPart);
		}
		//printf("locking mutex to write and return buffer to free list\n"); fflush(stdout);
		if (pthread_mutex_lock(&G->fq.mutex) != 0) {
			perror("pthread_mutex_lock: ");
			exit(1);
		}

#ifdef LESS_WRITING_NO_UNMAPPED
		if ((wb.header.chr & 0xffff) != kUnmappedChrID)
		{
#endif
			int fd = open(fname, O_WRONLY|O_CREAT|O_APPEND, 0664);
			if (fd == -1) {
				perror("open:");
				exit(1);
			}
			if (write(fd,dataBlock,wb.header.blockLength) != wb.header.blockLength) {
				perror("write failed");
				exit(1);
			}
			//fsync(fd);   // FIXME   was here before the O_SYNC
			close(fd);
#ifdef LESS_WRITING_NO_UNMAPPED
		}
#endif

		atomic_fetch_add(G->WrittenCnt, block->cnt);

		block->cnt = 0;
		block->hdrSize = 0;
		block->qualSize = 0;
		for (unsigned int i = 0; i < kMAXnbMatches; i++) {
			block->mmSize[i] = 0;
			block->cigarSize[i] = 0;
		}
		block->p4Size = 0;
		G->fq.freeBlock[G->fq.nbFreeBlocks++] = block;
		if (pthread_mutex_unlock(&G->fq.mutex)) {
			perror("pthread_mutex_unlock: ");
			exit(1);
		}
	} // infinite loop
	freeBlock(&wb);
	
} // flushBlockThread
//---------------------------------------------------------------

//---------------------------------------------------------------
// This is a single thread in charge of writing mapping results
// It merges the three tables created by the readFASTQ thread
// hdr tags quality
// and adds the genomic position present in the extrslts table
//
// This is an infinite loop which waits for buffers to be treated
// WARNING: buffers will not necessarily be written in the order they were filled.
//
void *Writer(WriterPool_t * const restrict Pool)
{	
	GTLBUFFERS G;
	COUNTS cnts;
	pthread_t bufferthread[Pool->FlushThreadPerWriter];
	const char * const restrict seedfn = Pool->seedfn;
	unsigned int WriterId;
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// Initialization
	{
		const pthread_t Me = pthread_self();
		WriterId = 0U;
		while (WriterId < Pool->num_threads) {
			if (Pool->threads[WriterId] == Me) break;
			WriterId++;
		}
		memset(&G,0,sizeof(GTLBUFFERS));
		pthread_mutex_init(&G.fq.mutex,NULL);
		pthread_cond_init(&G.fq.not_empty_cond,NULL);
		pthread_cond_init(&G.fq.not_full_cond,NULL);
		sprintf(&G.seedfn[0],"%s_%c", seedfn, (char) 'A' + WriterId);
		//cnts = Pool->cnts + WriterId;
		memset(&cnts,0,sizeof(COUNTS));

		initBuffers(&G, internalbuf, chrMax, chrPartOffset, Pool->PairedEnd ? 0 : 1);
		G.WrittenCnt = &(Pool->WrittenAlignments);
	}
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// Launch buffer threads
#ifdef USE_AFFINITY
	pthread_attr_t * restrict ThreadAttr;
	if (Pool->Affinities && Pool->Affinities[WriterId].Slaves) {
		ThreadAttr = (pthread_attr_t*) alloca(Pool->FlushThreadPerWriter*sizeof(pthread_attr_t));
		for (int cpu = 0; cpu < Pool->FlushThreadPerWriter; cpu++) {
			pthread_attr_init(&ThreadAttr[cpu]);
			if (pthread_attr_setaffinity_np(&ThreadAttr[cpu], sizeof(Affinity_Mask_t), (cpu_set_t*) &(Pool->Affinities[WriterId].Slaves[cpu]))) {
				fputs("Unable to set flush buffer thread affinity\n", stderr);
				exit(1);
			}
			if (pthread_create (&bufferthread[cpu], &ThreadAttr[cpu], (void* (*)(void*)) flushBlockThread, (void*) &G)) {
				printf("Error: Failed creating flushBlockThread thread %d\n",cpu);
				exit(1);
			}
		}
	}
	else
#endif
	{
		for (int cpu = 0; cpu < Pool->FlushThreadPerWriter; cpu++) {
			if (pthread_create (&bufferthread[cpu], NULL, (void* (*)(void*)) flushBlockThread, (void*) &G)) {
				printf("Error: Failed creating flushBlockThread thread %d\n",cpu);
				exit(1);
			}
		}
	}
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// Mark thread as alive (initialized)
	atomic_fetch_add(&(Pool->num_threads_alive), 1);
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// Paired end
	if (Pool->PairedEnd) {
		while(Pool->threads_keepalive) {
			bsem_wait(Pool->ToBeProcessed_q->has_items);
			
			if (Pool->threads_keepalive){
				/* Annonce thread working */
				atomic_fetch_add(&(Pool->num_threads_working), 1);
				
				/* Read job from queue and execute it */
				pthread_mutex_lock(&(Pool->ToBeProcessed_q->rwmutex));
				AlignerOutputBlock_t * const restrict task = (AlignerOutputBlock_t *) jobqueue_pull(Pool->ToBeProcessed_q);
				pthread_mutex_unlock(&(Pool->ToBeProcessed_q->rwmutex));
				
				if (task) {
					if (verbose & 0x1)
						printf("Writer %u got a block of %lu bytes to dump\n", WriterId, task->length);
					if (task->length > 0) {
						//////////////////////////////////////////////////////////////////////////////////////////////
						// Empty the aligner output block
						outputToPairedBuffer(&G, task->data, &cnts, Pool->Genome, task->length); // also generate TL blocks
					}
					/*else {
						Pool->threads_keepalive = 0;
						bsem_post_all(Pool->ToBeProcessed_q->has_items);
					}*/
					//////////////////////////////////////////////////////////////////////////////////////////////
					// Put the block back in its memory slot queue
					task->length = 0UL;
					atomic_store(&(task->Stored),0UL);
					if (task->BelongsTo_q) {
						pthread_mutex_lock(&(task->BelongsTo_q->rwmutex));
						assert(task->BelongsTo_q->has_items != NULL);
						jobqueue_push(task->BelongsTo_q, (job_t*) task);
						pthread_mutex_unlock(&(task->BelongsTo_q->rwmutex));
						task->BelongsTo_q = NULL; // ??? just to be sure it is not reused
					}
				}
				
				/* Annonce thread not working animore*/
				atomic_fetch_sub(&(Pool->num_threads_working), 1);
			}
		}
		printf("Writer %u leaving while queue has length %i\n", WriterId, Pool->ToBeProcessed_q->len); fflush(stdout);
		Pool->cnts[WriterId] = cnts;
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// Single
	else {
		
	}
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// Flush all remaining buffers
	flushAllBuffers(&G);
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// Wait until all workers are terminated
	for (int cpu = 0; cpu < Pool->FlushThreadPerWriter; cpu++) {
		if (pthread_join (bufferthread[cpu], NULL)) {
			printf("Error: Failed bufferthread pthread_join %d\n",cpu);
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// Mark thread as done
	atomic_fetch_sub(&(Pool->num_threads_alive), 1);

	printf("Writer %u redeemed its flush workers, closing...\n", WriterId);	
	return (void*) 0;
}
//---------------------------------------------------------------

WriterPool_t* createWriterPool(const Genome_t * const restrict Genome,
                               jobqueue_t * const restrict AlignerProcessedQueue,
                               const char * const restrict seedfn, COUNTS * const restrict cnts,
					                     const cpu_pool_t * const restrict Affinities,
                               const bool PairedEnd, const unsigned int N, const unsigned int nFlushBuffers)
{
	WriterPool_t * const restrict data = (WriterPool_t*) malloc(sizeof(WriterPool_t));
	if (data == NULL) goto fin;
	
	data->threads = (pthread_t*) malloc(N*sizeof(pthread_t));
	if (data->threads == NULL) goto fin1;
	
	data->cnts = cnts;
	data->num_threads = N;
	data->PairedEnd = PairedEnd;
	data->seedfn = seedfn;
	data->ToBeProcessed_q = AlignerProcessedQueue,
	data->threads_keepalive = 1;
	data->FlushThreadPerWriter = nFlushBuffers;
	atomic_store(&(data->num_threads_alive), 0);
	atomic_store(&(data->num_threads_working), 0);
	atomic_store(&(data->WrittenAlignments), 0);
	data->Genome = Genome;
#ifdef USE_AFFINITY
	data->Affinities = Affinities;
#endif

	initPackNT2bitsTable();
	virtchr = Genome->virtchr;

	internalbuf = 0U;
	assert(chrMax < kMAXCHRcnt);
	for (int i = 1; i <= chrMax; i++) {
		chrPartOffset[i] = internalbuf;
		internalbuf += virtchr[i].len / kCountsPerChrPart;
		internalbuf += 1;
	}

	pthread_attr_t * restrict ThreadAttrPtr;
#ifdef USE_AFFINITY
	if (Affinities) {
		pthread_attr_t attr;
		pthread_attr_init(&attr);
		for (unsigned int i=0; i<N; i++) {
			if (pthread_attr_setaffinity_np(&attr, sizeof(Affinity_Mask_t), (cpu_set_t*) &(Affinities->Master))) {
				fputs("Unable to set writer thread affinity\n", stderr);
				goto fin2;
			}
			if (pthread_create(data->threads + i, &attr, (void* (*)(void*)) Writer, (void*) data) != 0) {
				printf("%s: error creating thread\n", __FUNCTION__);
				int d = i;
				while (--d >= 0) pthread_kill(data->threads[d], SIGKILL);
				goto fin2;
			}
		}
	}
	else
#endif
	{
		for (unsigned int i=0; i<N; i++) {
			if (pthread_create(data->threads + i, NULL, (void* (*)(void*)) Writer, (void*) data) != 0) {
				printf("%s: error creating thread\n", __FUNCTION__);
				int d = i;
				while (--d >= 0) pthread_kill(data->threads[d], SIGKILL);
				goto fin2;
			}
		}
	}

	return data;
	fin2:
		free(data->threads);
	fin1:
		free(data);
	fin:
		return NULL;
}
//---------------------------------------------------------------

void freeWriterPool(WriterPool_t * restrict Pool)
{
	free(Pool->threads);
	free(Pool);
	virtchr = NULL;
}
//---------------------------------------------------------------

void terminateWriterPool(WriterPool_t * restrict Pool)
{
	/* End each thread 's infinite loop */
	Pool->threads_keepalive = 0;

	/* Give one second to kill idle threads */
	double TIMEOUT = 1.0;
	time_t start, end;
	double tpassed = 0.0;
	time (&start);
	while (tpassed < TIMEOUT && atomic_load(&(Pool->num_threads_alive))){
		assert(Pool->ToBeProcessed_q->has_items != NULL);
		bsem_post_all(Pool->ToBeProcessed_q->has_items);
		time (&end);
		tpassed = difftime(end,start);
	}

	/* Poll remaining threads */
	struct timespec polling_interval = { .tv_sec = 1, .tv_nsec = 0 };
	while (atomic_load(&(Pool->num_threads_alive)) > 0) {
		assert(Pool->ToBeProcessed_q->has_items != NULL);
		bsem_post_all(Pool->ToBeProcessed_q->has_items);
		nanosleep(&polling_interval, NULL);
	}

}
//---------------------------------------------------------------

/* vim: tabstop=2 shiftwidth=2
 */
