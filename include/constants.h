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
#ifndef _CONSTANTS_H
#define _CONSTANTS_H

//--------------------------------------------------------------- 
// FASTQ READER
//--------------------------------------------------------------- 
#define kFASTQReaderLineBufferSize 1048576 
#define kMaxTagsPerBuffer (2*65536)

//--------------------------------------------------------------- 
// TAGS 
//--------------------------------------------------------------- 
// TODO: this is still in GTLkonst.h
#define kTagSize 192 /* Should be kept a multiple of 16 for alignment and packing */
#define kHdrSize  64 /* Should be kept a multiple of  8 for alignment and packing */

//--------------------------------------------------------------- 
// DECODER 
//--------------------------------------------------------------- 							
/* note: as each nt is encoded on 2 bits, at most eight 16mers can fit in a 256
 * bit register register, but we repeat the operation twice to match the
 * capability of the MIC
 */
#define kMAX16mersToTest 16
/* how many hits (at most) to report per tag */
#define kDesiredHitCnt 16
/* the length of the tags stored in the main table */
#define kbaitLen 18
/* the number of following nucleotides after 16 to compose the full bait length */
#define kFollowNtLen (kbaitLen-16)
/* how many of the following nt elements can we pack in a 32-bit word */
#define kFollowNtCount (16/kFollowNtLen)
/* count of elements in a full table of 16 nucleotides encoded as 2-bit values is full 32-bit */
#define k16NT_ENTRIES_COUNT (1UL << (2 * 16))
/* table of 16nt is 4 giga entries, but since we fix the first NT (usually to A) we can divide
 * by 4 to go down to 1 giga entries 
 */
#define k15NT_ENTRIES_COUNT (k16NT_ENTRIES_COUNT >> 2)
// 1 giga entries times 4 bytes => 4 gigabytes
#define kFOLLOWING_NT_DATA_SIZE (k15NT_ENTRIES_COUNT * 4)
// 1 giga entries times 8 2-mers per 32-bit word divided by 8 since we only need 1 bit per entry
#define kSTRAND_DATA_SIZE (k15NT_ENTRIES_COUNT * kFollowNtCount / 8)
#define kVALID_DATA_SIZE kSTRAND_DATA_SIZE
/* To prevent faslse sharing, we should make sure batch of tags are on separate cache lines.
	 Therefore decoder should work on aligned arries on 64 bytes with batch such that 
	     (kDesiredHitCnt*sizeof(DecoderElement_t)) = (kDesiredHitCnt*8) 
	 is a multiple of the cache line (64b)
	 */
#define kDesiredHitCnt 16        /* how many hits (at most) to report per tag */
#define kDecoderTagBatch 4096

//--------------------------------------------------------------- 
// ALIGNER 
//---------------------------------------------------------------
#define kGAPPED_ALIGN_GENOME_LENGTH  256  /* slow but local gapped alignment */
#define kMaxGenomeChunkSize (1024-kTagSize)
#define kMinTagLenNotSoftClipped 32 /* avoid soft clipping of whole tag and aligning it in stupid places (32 - 8 is leaving only 24 perfect matches...) */

#define kMaxEncodedCigar 15
#define kMAXcigarSize (2*kMaxEncodedCigar+1+1)
#define kMaxMismatch  8
#define kMaxEncodedMismatch 18           /* keep more room to be able to encode N as they do not count as mismatches. */
#define kMAXmismatchSize (2*kMaxEncodedMismatch)+1
#define kMMTERMINATOR 255 /* end of mismatch encoding; WARN: means that max taglength is 253 !!! */ 
#define kMMSOFTCLIP 254 /* end of mismatch encoding; WARN: means that max taglength is 253 !!! */ 
#define kCIGARTERMINATOR 255
#define kUNMAPPED 999

#define kSuccessfullyEncoded 		1
#define kSoftClipRescued 				2
#define kDirectFromShift 				4
#define kMMEncoded 							8
#define kContainsInDel		 			16
#define kHasGoneThruAlignment 	32
#define kPerfectMatch 					64
#define kBestOfTags 						128
#define kLargeAlignmentSearch 	256
#define kAllIsN 								512
#define kContainsN 							1024
#define kTrimmed                2048


#endif

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
