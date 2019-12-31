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
#define MIC_MATCHER_LISTENING_PORT 2048
//(SCIF_PORT_RSVD+2)
//#define CHECK_DOUBLE_MAPPING

#define kALL_DONE 0xFFFFFFFF /* message sent to MIC at end of process */

#define kNOT_ACGT 0xFF

#ifndef kMaxTagsPerBuffer
// FIXME - this is also defined in constants.h
#define kMaxTagsPerBuffer (2*65536)
#endif

#define k_rsltTogenome_PROCESSINGSLAVESCNT 8
#define k_rsltTogenome_TAGcntBySlave (kMaxTagsPerBuffer/k_rsltTogenome_PROCESSINGSLAVESCNT)

#define kbaitLen 18                      /* the length of the tags stored in the main table */
#define kFollowNtLen (kbaitLen-16)       /* the number of following nucleotides after 16 to compose the full bait length */
#define kFollowNtCount (16/kFollowNtLen) /* how many of the following nt elements can we pack in a 32-bit word */
#define kFollowNtEnum ((1 << kFollowNtLen) * (1 << kFollowNtLen)) /* number of combinations is 4^kFollowNtLen */
#define kDesiredHitCnt 16                /* how many hits (at most) to report per tag */
#define kMaxMismatch  8                  /* full processer version uses only this value unlike older version using kMaxEncodedMismatch */
#define kMaxEncodedMismatch 18           /* keep more room to be able to encode N as they do not count as mismatches. */
#define kMAXmismatchSize (2*kMaxEncodedMismatch)+1

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

//#define kFailedGappedAlignment 16
//====================== Size of our tables read from disk ========================================

// count of elements in a full table of 16 nucleotides encoded as 2-bit values is full 32-bit
#define k16NT_ENTRIES_COUNT (1UL << (2 * 16))
// table of 16nt is 4 giga entries, but since we fix the first NT (usually to A) we can divide by 4 to go down to 1 giga entries
#define k15NT_ENTRIES_COUNT (k16NT_ENTRIES_COUNT >> 2)
// 1 giga entries times 4 bytes => 4 gigabytes
#define kFOLLOWING_NT_DATA_SIZE (k15NT_ENTRIES_COUNT * 4)
// 1 giga entries times 8 2-mers per 32-bit word divided by 8 since we only need 1 bit per entry
#define kSTRAND_DATA_SIZE (k15NT_ENTRIES_COUNT * kFollowNtCount / 8)
#define kVALID_DATA_SIZE kSTRAND_DATA_SIZE
#define kUNMAPPED 999

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
