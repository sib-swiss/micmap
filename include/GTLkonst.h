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
#ifndef GTLKONST_H
#define GTLKONST_H

// Tag information held in chr unsigned int
/* 32 bit unsigned integer holding the chromosome number within block is decomposed as follow:
 * 31               15              0
 * |................|................|
 * |-----TAG 1------|-----TAG 2------|
 * |x...............|x...............| <- is tag reversed?  : 1 bit
 * |.x..............|.x..............| <- tag contains N?   : 1 bit
 * |..x.............|..x.............| <- tag mapped        : 1 bit
 * |...x............|...x............| <- tag softclipped?  : 1 bit
 * |....x...........|....x...........| <- tag is encoded?   : 1 bit
 * |.....x..........|.....x..........| <- MM contains N?    : 1 bit
 * |......x.........|......x.........| <- MM is only N?     : 1 bit
 * |.......x........|.......x........| <- Has a cigar ?     : 1 bit
 * |........x.......|........x.......| <- Has a MM string?  : 1 bit
 * |.........x......|.........x......| <- Has Packed data?  : 1 bit
 * |..........x.....|..........x.....| <- still available   : 1 bits
 * |...........xxxxx|...........xxxxx| <- chromosome number : 5 bits
 */

#define TAG_IS_REVERSED 0x00008000
#define TAG_CONTAINS_N  0x00004000
#define TAG_MAPPED      0x00002000
#define TAG_SOFTCLIPPED 0x00001000
#define TAG_IS_ENCODED  0x00000800
#ifndef __FULL_PROCESSER__
#define TAG_CHR_MASK    0x00007FFF
#else
#define TAG_CHR_MASK    0x0000001F
#endif
#define MM_CONTAINS_N   0x00000400
#define MM_IS_ONLY_N    0x00000200

#define TAG_HAS_CIGAR   0x00000100
#define TAG_HAS_MM      0x00000080
#define TAG_HAS_PACKED  0x00000040

#define GET_TAG1_CHR(a)      (((a) >> 16 ) & TAG_CHR_MASK) 
#define GET_TAG2_CHR(a)      (((a)       ) & TAG_CHR_MASK)
#define IS_TAG1_REVERSED(a)  ((a) & (TAG_IS_REVERSED << 16))
#define IS_TAG2_REVERSED(a)  ((a) & (TAG_IS_REVERSED      ))
#define TAG1_CONTAINS_N(a)   ((a) & (TAG_CONTAINS_N << 16))
#define TAG2_CONTAINS_N(a)   ((a) & (TAG_CONTAINS_N      ))
#define TAG1_IS_MAPPED(a)    ((a) & (TAG_MAPPED << 16))
#define TAG2_IS_MAPPED(a)    ((a) & (TAG_MAPPED      ))
#define TAG1_SOFTCLIPPED(a)  ((a) & (TAG_SOFTCLIPPED << 16))
#define TAG2_SOFTCLIPPED(a)  ((a) & (TAG_SOFTCLIPPED      ))
#define TAG1_ENCODED(a)      ((a) & (TAG_IS_ENCODED << 16))
#define TAG2_ENCODED(a)      ((a) & (TAG_IS_ENCODED      ))

#define TAG1_MM_HAS_N(a)     ((a) & (MM_CONTAINS_N << 16))
#define TAG2_MM_HAS_N(a)     ((a) & (MM_CONTAINS_N      ))
#define TAG1_MM_IS_ONLY_N(a) ((a) & (MM_IS_ONLY_N << 16))
#define TAG2_MM_IS_ONLY_N(a) ((a) & (MM_IS_ONLY_N      ))

#define SET_TAG1_CHR(a,b)    a |= (((b) & TAG_CHR_MASK) << 16)
#define SET_TAG2_CHR(a,b)    a |= (((b) & TAG_CHR_MASK)      )
#define SET_TAG1_REVERSED(a) a |= (TAG_IS_REVERSED << 16)
#define SET_TAG2_REVERSED(a) a |= (TAG_IS_REVERSED      )
#define SET_TAG1_FLAGS(a,b)  a |= ((b) << 16)
#define SET_TAG2_FLAGS(a,b)  a |= ((b)      )

// Transport Layer output blocks
#define kBlockMaxCount 65536 // hmm, for starters, why not...  This is roughly 10 kilobase of reference covered at a depth of 1000 using two reads of 100 nt

#define kMaxHdrSize 128
#define kMaxReadLen 256 // use this, for example in GTL structures to make sure we can read any files, also those with long reads (up to 256...)

// bits [7..0]
#define kTLBHflagBlockType       0x01
#define kTLBHflagPairedReads     0x02
#define kTLBHflagFixedLength     0x04
#define kTLBHflagPairPosBlock    0x08
//#define kTLBHflagChimeraPosBlock 0x10
#define kTLBHflagPairChrPosBlock 0x10
#define kTLBHflagNMismatchBlock  0x20
#define kTLBHflagMismatchBlock   0x40
#define kTLBHflagCigarBlock      0x80

// bits [15..8]
#define kTLBHflagUnmappedBlock   0x01
#define kTLBHflagHeaderBlock     0x02
#define kTLBHflagQualityBlock    0x04
#define kTLBHflagMultipleMatches 0x08
#define kTLBHflagSortedByPos     0x10
// FIXME - this is a bit strange as the 2 below are in fact the same block as the Unmapped... ???
#define kTLBHflagChimeraBlock    0x20  // FIXME !!  could be stored on two bits as a choice with    kTLBHflagUnmappedBlock  ?
#define kTLBHflagHalfmapBlock    0x40

#define kUnmappedChrID 0 /* because it makes sense */
#define kHalfMapChrID 50 /* because it is only 50% mappable */
#define kChimeraChrID 57 /* because Chimeras are found in the Fabulae 57 by Hyginus, according to Wikipedia. */

//=================================================================================================

// The maximum length of a read... actually last 2 bytes are needed by length and hight-quality length
// Use this in the matcher, and where size matters for speed
#ifndef kTagSize
#define kTagSize 192  /* max tag size not sure if really needed to be aligned on 64bits */
// FIXME - this constant would probably be better written as kMaxGapLen + taglen + kMaxGapLen and become variable
// need to revisit after FDA challenge, at the moment kMaxGapLen is hidden as 32 around line 1900 of align_MIC.c
// also need to make sure the pftools aligner can deal with a variable instead of an aligned constant
#endif
#if kMaxReadLen < kTagSize
#error kMaxReadLen must be >= kTagSize
#endif

#define kHdrSize 64   /* max tagname size */
#define kMaxEncodedCigar 15
#define kMAXcigarSize (2*kMaxEncodedCigar+1+1)

// max number of matches in case we record multiple matches (kTLBHflagMultipleMatches)
#define kMAXnbMatches 2

#define kMMTERMINATOR    255 /* end of mismatch encoding; WARN: means that max taglength is 253 !!! */ 
#define kMMSOFTCLIP      254 /* end of mismatch encoding; WARN: means that max taglength is 253 !!! */ 
#define kCIGARTERMINATOR 255

#define k32bREVCOMP_TAG1_BIT    0x80000000  
#define k32bREVCOMP_TAG2_BIT    0x40000000 
#define k32bNEGATIVE_OFFSET_BIT 0x20000000 
#define k32b2nd_MATCH_CHR_MASK  0x1FF80000 // seems maximum number of chromosomes is 630 (in a fern)
#define k32bOFFSET_VALUE_MASK   0x0007FFFF

#define kFLUSHQUEUEBUFFERSIZE 8
#define kFLUSHTHREADCOUNT 4

#define kCountsPerChrPart 0x4000000

#endif /* GTLKONST_H */

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
