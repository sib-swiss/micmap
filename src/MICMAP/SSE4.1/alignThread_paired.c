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
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <sys/mman.h>
#include <unistd.h>
#include <pthread.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <smmintrin.h>
#include "ALIGNkonst.h"
#include "virt_chr.h"
#include "../cpuThreadAlignPool.h"

#define MYTEST 0

//========================================================================================================
// GLOBAL DEFINITIONS
//========================================================================================================
#if (kTagSize & 0x3F) != 0
#error "kTAGSize must be a multiple of 64"
#endif

#if kDesiredHitCnt != 16
#error "kDesiredHitCnt has to be 16 in this file"
#endif

#define __WITH_SOFTCLIP_EXTENSION__
// WARNING: This could be not enough as long search may overflow here
#define STATE_MEMORY_SIZE ((3*kTagSize)/2)

#if kMAXmismatchSize > STATE_MEMORY_SIZE
#error "Allowed state memory string size is not sufficient for allowed resulting kMAXmismatchSize"
#endif

#if (_II == 0) && (_DD == 0) && (_MM == 0 )
//#define PFTOOLS_SIMPLIFIED
#else
#warning "Simplified version of the aligner not possible due to non zero extensions II, DD or MM"
#endif

/* how many trials at most we are going to test. */
#define kMaxTrials  (kDesiredHitCnt*kDesiredHitCnt+2*kDesiredHitCnt) 
#define ALIGNER_THREAD_TAGS_BATCH_SIZE 256
#define HUGE_PAGE_SIZE (2*1024*1024)

//========================================================================================================
// BENCHMARKS
//========================================================================================================
// #define BENCHMARKS
#ifdef BENCHMARKS
#include "benchmarks.h"
#else
#define DO_TIME(command,array,index) command;
#endif

//========================================================================================================
// STRUCTURES
//========================================================================================================
typedef struct InternalAlignment_s {
	char * restrict States;
	int score;
	unsigned int DecisionTreeHistory;
	unsigned int chr_pos;
#ifndef NDEBUG
	unsigned int GenomeSeekStart;       // starting position to seek tag
	unsigned int GenomSeekStop;         // ending position to seek tag
#endif
	unsigned short chr;
	unsigned short MismatchCount;
	unsigned char TotalRequiredMMCount;
	unsigned char TotalRequiredCigCount;
	unsigned char revNmer;
	unsigned char SoftClipMismatchToRemove;
	char StateMemory[STATE_MEMORY_SIZE];
#ifndef NDEBUG
	unsigned char SoftClipEvictedMMSize;
	unsigned char SoftClipEvictedCigSize;
	unsigned char AnchorId;
#endif
#ifdef DEBUG_TAG_MAPPING_PROCESS
	unsigned char AnchorPositionInTag;
#endif
} InternalAlignment;

#define __DECLARE_DEFINITIONS__
#ifdef PFTOOLS_SIMPLIFIED
#include "pftools_std_simplified.h"
#else
#include "pftools_std.h"
#endif
#undef __DECLARE_DEFINITIONS__

//========================================================================================================
// EXTERN CONSTANTS
//========================================================================================================
#include "../global_var.h"
extern unsigned int gOutputBlocksCnt;

//========================================================================================================
// LOCAL CONSTANTS
//========================================================================================================

#define __DECLARE_STATIC_VARIABLES__
#ifdef PFTOOLS_SIMPLIFIED
#include "pftools_std_simplified.h"
#else
#include "pftools_std.h"
#endif
#undef __DECLARE_STATIC_VARIABLES__

//========================================================================================================
// LOCAL FUNCTIONS
//========================================================================================================
#define __EXTRA_FUNCTIONS__
#ifdef PFTOOLS_SIMPLIFIED
#include "pftools_std_simplified.h"
#else
#include "pftools_std.h"
#endif
#undef __EXTRA_FUNCTIONS__

/* Include common function to paired and single */
#include "alignThread_common.c"


// REVERSE COMPLEMENT
#define ReverseComplementOutputAligned(in,out,length) Reverse(in,out,gIUPACrevcomp,length)
static inline void __attribute__((always_inline))
Reverse(const unsigned char * const restrict in, unsigned char * const restrict out,
        const unsigned char * const restrict LUT, const size_t length)
{
	int lN = (int) length - 1;
	int li = 0;
	__assume_aligned(out, 16);
	while (lN >= 0) {
		out[li++] = LUT[in[lN--]];
	}
}
//---------------------------------------------------------------

//========================================================================================================
// GLOBAL FUNCTIONS
//========================================================================================================
#if kMaxGenomeChunkSize <= kGAPPED_ALIGN_GENOME_LENGTH
#error "kGenomeChunkSize should be bigger than kGAPPED_ALIGN_GENOME_LENGTH"
#endif

#if (16*(kTagSize>>4)) != kTagSize
#error "kTagSize should be a multiplier of 16"
#endif

#if ((16+16+kTagSize+1)*(kMaxGenomeChunkSize + 1 + kTagSize + 16 + 16)*16) > 2*HUGE_PAGE_SIZE
#error "Allocated buffer for matrix on Huge-TLB is not big enough"
#endif

void* sse41_pair(cpuAlignerPool_t * const restrict data)
{
	/********************************************************************************************************************************/
	/*                                                   INTERNAL MACROS                                                            */
	/********************************************************************************************************************************/

	/********************************************************************************************************************************/
	/*                                                   INTERNAL STRUCTURES                                                        */
	/********************************************************************************************************************************/	
	typedef struct Trial_s {
		InternalAlignment * restrict tag1;
		InternalAlignment * restrict tag2;
	} Trial_t;
	
	/********************************************************************************************************************************/
	/*                                                   STATIC VARIABLES                                                           */
	/********************************************************************************************************************************/
	
	/********************************************************************************************************************************/
	/*                                                STACK ARRAY VARIABLES                                                         */
	/********************************************************************************************************************************/
	// Hit anchor choices data storage
	InternalAlignment tagchoices[4*kDesiredHitCnt];
	
	// Offsets within tagchoices according to the score vector position for the kept anchors
	int vOffsets[2*kDesiredHitCnt] __attribute__((aligned(64)));
	
	// Score vector for ordering of kept anchors
	int vScores[2*kDesiredHitCnt] __attribute__((aligned(64)));
	
	// Reverse Tag storage
	unsigned char RevTagSpace[2][kTagSize] __attribute__((aligned(64)));
	
#ifdef PURIFY_ANCHORS
	// Genome position of each anchor to allow filtering of similar regions
	long tagGenomePosition[2*kMaxTrials]; 
#endif
	
	// Pair trials
	Trial_t trials[kMaxTrials];

	// Alignments 
	AlignmentResults_t localResults[4];
	
	// OutputHeader
	struct outputPairBlockHeader outputHeader;
	
	// Processer statistics 
	union Statistics stats;
	stats.Vector[0] = _mm_setzero_si128();
	stats.Vector[1] = _mm_setzero_si128();
	stats.Vector[2] = _mm_setzero_si128();
	stats.Vector[3] = _mm_setzero_si128();
	
	// Error return code
	int err = 0;
	
	/********************************************************************************************************************************/
	/*                                                    OTHER VARIABLES                                                           */
	/********************************************************************************************************************************/
	const size_t Me = (size_t) pthread_self();
#ifdef BENCHMARKS
	tic_t BenchmarkCounters[MAX_PROCESS_TOPICS];
	memset(BenchmarkCounters, 0, MAX_PROCESS_TOPICS*sizeof(tic_t));
#endif
	int PosixAllocation = 0;
	// Alignment storage array, 3 times because of insertion, deletion, match
	__m128i * restrict IDM_Matrix = (__m128i*) mmap(NULL, 2*HUGE_PAGE_SIZE,
	                                                PROT_READ | PROT_WRITE,
	                                                MAP_PRIVATE | MAP_ANONYMOUS | MAP_POPULATE | MAP_HUGETLB,
	                                                -1, 0);
	if (IDM_Matrix == MAP_FAILED) {
		fprintf(stderr, "Aligner thread %lu is unable to allocate Huge TLB, trying standard!\n", Me);
		PosixAllocation = posix_memalign((void**)&IDM_Matrix, HUGE_PAGE_SIZE, 2*HUGE_PAGE_SIZE);
		if (PosixAllocation == 0) {
			PosixAllocation = 1;
			//if (madvise(IDM_Matrix, 2*HUGE_PAGE_SIZE, MADV_HUGEPAGE)) perror("madvise");
		}
		else {
			perror("posix_memalign");
			fprintf(stderr, "Aligner thread %lu is unable to allocate matrix memory!\n", Me);
			fflush(stderr);
			goto zeend;
		}
	}
	int * const restrict WorkSpace = _mm_malloc(2*(kTagSize+1)*sizeof(union sIOP)+63, 64);
	if (WorkSpace == NULL) {
		printf("Aligner thread %lu is unable to allocate WorkSpace!\n", Me);
		fflush(stdout);
		munmap(IDM_Matrix, 2*HUGE_PAGE_SIZE);
		goto zeend;
	}
#define __DECLARE_VARIABLES__
#ifdef PFTOOLS_SIMPLIFIED
#include "pftools_std_simplified.h"
#else
#include "pftools_std.h"
#endif
#undef __DECLARE_VARIABLES__
	
	/********************************************************************************************************************************/
	/*                                                     THREAD ACTIVE                                                            */
	/********************************************************************************************************************************/
	const unsigned int MyID = atomic_fetch_add(&data->num_threads_alive, 1);
	//printf("Appolo Aligner %u ready for work...\n", MyID); fflush(stdout);

	/********************************************************************************************************************************/
	/*                                                 LOOP UNTIL TOLD TO QUIT                                                      */
	/********************************************************************************************************************************/	
	int CurrentInputIndex;
	atomic_int * JobIndex = &(data->JobIndex);
	atomic_int * OutIndex = &(data->OutIndex);
	cpuAlignerJob_t * const restrict Jobs = data->Jobs;
	const unsigned char * const genome_tbl = data->AlignerBlock->Genome->table;
	const VIRTUALCHR * const restrict virtchr = data->AlignerBlock->Genome->virtchr;
	const unsigned long OutputBufferSize = (unsigned long) data->AlignerBlock->OutputBlockSize;
	cpuSpace_t * const restrict OutputBuffers = data->OutputBuffers;
	jobqueue_t * const restrict ToBeWritten_q = &(data->AlignerBlock->Pool->ToBeWritten_q);
	AlignerPool_t * const restrict AlignPool = data->AlignerBlock->Pool;
	const unsigned int GenomeChunkSize = data->AlignerBlock->GenomeChunkSize;
	
	while(1)
	{
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// Keep working on that batch until done
			CurrentInputIndex = atomic_load(JobIndex);
			if (CurrentInputIndex >= 0) {
			Compute:;
				const unsigned int TagCount = Jobs[CurrentInputIndex].TagCnt;
				do {
					/* Randomly choose a WorkLoad */
					unsigned int WorkLoad = ALIGNER_THREAD_TAGS_BATCH_SIZE;
					unsigned long int TreatedSoFar = atomic_fetch_add(&(Jobs[CurrentInputIndex].Treated), WorkLoad);
						
					///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					// ALIGN
					if ((TreatedSoFar) < TagCount) {
						/* Advise other to fetch new batch */
						const int IamLast = ((TreatedSoFar+WorkLoad) >= TagCount);
						if (IamLast) {
							int itmp = CurrentInputIndex;
							WorkLoad = TagCount - (unsigned int) TreatedSoFar;
							// We only change to Empty if JobIndex is still ours 
							atomic_compare_exchange_strong(JobIndex, &itmp, -1);
						}
											
						/* 
						 * =========================================================================================================================
						 * Perform your partial part 
						 * =========================================================================================================================
						 */
						assert(CurrentInputIndex >= 0);
						//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
						// Process all tags within the batch
						//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
						unsigned int KeptAnchorCount, EvenAnchorCount;
						const ReaderTag_t * const restrict TagsData = Jobs[CurrentInputIndex].TagsData + TreatedSoFar;
						const DecoderElement_t (* const restrict DecoderData)[kDesiredHitCnt] = Jobs[CurrentInputIndex].Anchors + TreatedSoFar;
						stats.Counter[FullAlignment] += WorkLoad ;	
						for (unsigned int i = 0; i < WorkLoad; i++) {
							unsigned char * restrict ptrRevTagSpace;
							unsigned int EvenOrOddOffset;
							if (i & 1) {
								ptrRevTagSpace  = &RevTagSpace[1][0];
								EvenOrOddOffset = kDesiredHitCnt;
							}
							else {
								ptrRevTagSpace  = &RevTagSpace[0][0];
								EvenOrOddOffset = 0U;
								KeptAnchorCount = 0U;
							}
							
							//Indicate that the reverse tag space is not computed yet 
							*ptrRevTagSpace = '\0';
							
							//Clean score vector
							_mm_store_si128((__m128i*) &vScores[EvenOrOddOffset   ], __NLOW);
							_mm_store_si128((__m128i*) &vScores[EvenOrOddOffset+ 4], __NLOW);
							_mm_store_si128((__m128i*) &vScores[EvenOrOddOffset+ 8], __NLOW);
							_mm_store_si128((__m128i*) &vScores[EvenOrOddOffset+12], __NLOW);
							
							{
								const register __m128i __Zero = _mm_setzero_si128();
								_mm_store_si128((__m128i*) &vOffsets[EvenOrOddOffset   ], __Zero);
								_mm_store_si128((__m128i*) &vOffsets[EvenOrOddOffset+ 4], __Zero);
								_mm_store_si128((__m128i*) &vOffsets[EvenOrOddOffset+ 8], __Zero);
								_mm_store_si128((__m128i*) &vOffsets[EvenOrOddOffset+12], __Zero);
							}
							
							// Tag length
							const int taglen = (int) TagsData[i].AlignLen;
						
							//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
							// Loop on the anchor and test them, perfect and 1 mismatch anchor are moved to top
							//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
							for (int hitnum = 0; hitnum < kDesiredHitCnt; hitnum++) {
								InternalAlignment * restrict alntrial = tagchoices + KeptAnchorCount;
								assert(KeptAnchorCount < 2*kDesiredHitCnt);
								
								//------------------------------------------------------------------------------------------------------------------------
								// CHECK IF WE HAVE A HIT OR NOT : GOT ONE :-)
								//------------------------------------------------------------------------------------------------------------------------
								if (WHICH16MER_pos(DecoderData[i][hitnum].which16mer) != 0) {
									unsigned int chr = WHICH16MER_chr(DecoderData[i][hitnum].which16mer);	// virtual chr
									unsigned int pos = WHICH16MER_pos(DecoderData[i][hitnum].which16mer);
									//----------------------------------------------------------------------------------------------------------------------
									// DECODE HIT GENOME POSITION
									//----------------------------------------------------------------------------------------------------------------------
									assert(revVchr[chr].nbChr > 0);
									
									if (revVchr[chr].nbChr == 1) {
										chr = revVchr[chr].chr[0];
									}
									else {
										unsigned int k;
										for (k = 0; k < revVchr[chr].nbChr; k++)
										{
											const unsigned int c = revVchr[chr].chr[k];
											const unsigned int offset = virtchr[c].offset;
											const unsigned int len = virtchr[c].len;
											// a bit of a hack, because some reads map before the actual start of a chromosome...
											if (offset <= pos + taglen && offset + len >= pos) break;
										}
										if (k >= revVchr[chr].nbChr)
										{
#ifdef DEBUG_TAG_MAPPING_PROCESS
											fprintf(stderr,"chr %u pos %u taglen %u\n",chr,pos,taglen);
											for (int kk = 0; kk < revVchr[chr].nbChr; kk++)
											{
												unsigned int cc = revVchr[chr].chr[kk];
												fprintf(stderr,"k %u c %u offset %u len %u\n",kk,cc,virtchr[cc].offset,virtchr[cc].len);
											}
#endif
											k = revVchr[chr].nbChr - 1;
										}
										assert(k < revVchr[chr].nbChr);
										chr = revVchr[chr].chr[k];
										if (pos > virtchr[chr].offset)
											pos -= virtchr[chr].offset;
										else
											pos = 0;
									}
									
									//-----------------------------------------------------------------------------------------------------------------------
									// SET GENOME POSITION OF TAG, EVICT NEARBY ANCHORS i.e. distance <=32, could do much better here TS
									// This could be DANGEROUS AS WE CURRENTLY DO NO EXTEND THE REGION TO CHECK !!!
									//-----------------------------------------------------------------------------------------------------------------------
									const unsigned int BaseChrpos = (virtchr[chr].chr << 28) + virtchr[chr].offset;
									const unsigned char revNmer = (DecoderData[i][hitnum].info & 0x01) ^ ((DecoderData[i][hitnum].info >> 1) & 0x01);
									
#ifdef PURIFY_ANCHORS
									{
										const long GenomePos = (((long) (revNmer)) << 32) | ((long) (BaseChrpos + pos));
										int Eject = 0;
										#pragma noprefetch
										for (unsigned int iKeptAnchor=0; iKeptAnchor<KeptAnchorCount; iKeptAnchor++) {
											Eject |= (labs(GenomePos - tagGenomePosition[iKeptAnchor]) <= 128);
											if (Eject) break;
										}
										if (Eject) {
											stats.Counter[HalfMapped]++;
										}
										tagGenomePosition[KeptAnchorCount] = GenomePos;
									}
#endif
#ifdef DEBUG_TAG_MAPPING_PROCESS
									printf("DBG: i=%d hitnum=%d tc=%d chr%d:%d\n\n",i,hitnum,KeptAnchorCount,chr,pos); fflush(stdout);
#endif
									// Most of the following paragraph is propably useless, come back to it later TS
									alntrial->chr_pos             = pos;
									alntrial->score               = NLOW;
									alntrial->DecisionTreeHistory = 0U;
									alntrial->chr                 = chr;
									alntrial->revNmer             = revNmer;
									alntrial->StateMemory[0]      = '\0';
#ifdef DEBUG_TAG_MAPPING_PROCESS
									alntrial->AnchorPositionInTag = (DecoderData[i][hitnum].info >> 8) & 0xFF;
#endif
#ifndef NDEBUG
									alntrial->AnchorId = (unsigned char) hitnum;
#endif

#ifdef DEBUG_TAG_MAPPING_PROCESS
									if (DecoderData[i][hitnum].info & 0x01) {
										if ((rslt[i][hitnum].info >> 1) & 0x01) // revNmer, but the tag was also reversed !!
											printf("Case revNmer  and revtag  chrpos=%d  matchpos=%d\n",alntrial->chr_pos,alntrial->AnchorPositionInTag);
										else
											printf("Case revNmer  and fwdtag  chrpos=%d  matchpos=%d\n",alntrial->chr_pos,alntrial->AnchorPositionInTag);
									}
									else {
										if ((rslt[i][hitnum].info >> 1) & 0x01)  // fwdNmer, but the tag was reversed !!
											printf("Case fwdNmer  and revtag  chrpos=%d  matchpos=%d\n",alntrial->chr_pos,alntrial->AnchorPositionInTag);
										else
											printf("Case fwdNmer  and fwdtag  chrpos=%d  matchpos=%d\n",alntrial->chr_pos,alntrial->AnchorPositionInTag);
									}
									fflush(stdout);
#endif
									//-----------------------------------------------------------------------------------------------------------------------
									// REVERT TAG AND STORE IT IF NEEDED
									//-----------------------------------------------------------------------------------------------------------------------
									const unsigned char * restrict ptrTag; __assume_aligned(ptrTag, 64);
									if (!(alntrial->revNmer)) {
										ptrTag = TagsData[i].Tag;
									}
									else {
										if (*ptrRevTagSpace == '\0') {
											ReverseComplementOutputAligned(TagsData[i].Tag, ptrRevTagSpace, taglen);
										}
										ptrTag = ptrRevTagSpace;
									} 
									
									//-----------------------------------------------------------------------------------------------------------------------
									// DO SHIFTS ALONG THE GENOME AND SCORE THAT HIT
									//-----------------------------------------------------------------------------------------------------------------------
									unsigned int mismatchCount;
									DO_TIME( mismatchCount = shiftTagAgainstGenome(ptrTag,taglen, &genome_tbl[BaseChrpos+pos], alntrial),\
																															 BenchmarkCounters, PROCESS_TOPIC_SHIFT_TAG);

								//-----------------------------------------------------------------------------------------------------------------------
								// CHECK IF DIRECT HIT, CYCLE NEXT ON ISSUE 
								//-----------------------------------------------------------------------------------------------------------------------
								// if perfect match or 1 mismatch, verify we can encode it (N's)
								if (mismatchCount <= 1U) {
									EncodeMM(ptrTag, alntrial, &genome_tbl[BaseChrpos + alntrial->chr_pos], taglen);
									///////////////////////////////////////////
									// if encoding fails drop it !!!
									///////////////////////////////////////////
#if MYTEST									
									if (TagsData[i].ordinal == 0)  printf("SSE %u anchor %i MisMatch: %u\n", i, hitnum, mismatchCount);
#endif
									if ((alntrial->DecisionTreeHistory & kSuccessfullyEncoded)) {
										stats.Counter[directhit]++;
									}
									else
										continue;
								}
								//-----------------------------------------------------------------------------------------------------------------------
								// CHECK IF AN ALIGNMENT COULD ENHANCE THAT SCORE, CYCLE NEXT ON ISSUE
								//-----------------------------------------------------------------------------------------------------------------------
								else if (taglen < kMinTagLenNotSoftClipped) {
									continue;
								}
								else {
									const unsigned int genomeCheckStart = (pos > 32) ? BaseChrpos + (pos - 32) : (BaseChrpos + 1);
									const char * restrict const genome = (char *) ((uintptr_t) &genome_tbl[genomeCheckStart] & ~(15));
									const size_t GenomeLength = kGAPPED_ALIGN_GENOME_LENGTH + (size_t) ((uintptr_t) &genome_tbl[genomeCheckStart] & (15));
									const size_t MatrixLD     = 1UL + taglen;
									const size_t ReadLength   = (size_t) taglen;
									__m128i __BestScore, __location;
									
									// FILL IN MATRIX
#ifdef BENCHMARKS
									unsigned long s0 = timer();
#endif
									#define __FILL_MATRIX__
									#ifdef PFTOOLS_SIMPLIFIED
									#include "pftools_std_simplified.h"
									#else
									#include "pftools_std.h"
									#endif
									#undef __FILL_MATRIX__
#ifdef BENCHMARKS
									BenchmarkCounters[PROCESS_TOPIC_UNGAP_MATRIX].tics += timer() - s0;
									BenchmarkCounters[PROCESS_TOPIC_UNGAP_MATRIX].counter++;
#endif
									__BestScore = _mm_srai_epi32(__BestScore, SCORE_SHIFT);
									int Scores[4] __attribute__((aligned(16)));
									_mm_store_si128((__m128i*) Scores, __BestScore);
									
// 									const size_t num = (ReadLength-1) & (size_t) 0xF;
									int iseq; //= Id[num];
									unsigned int EndOnMatchOrDeletion; // = MorD[num];
									if (Scores[MATCH] >= Scores[DELETION]) {
										iseq = _mm_extract_epi32(__location, MATCH);
										EndOnMatchOrDeletion = 0U;
									}
									else {
										iseq = _mm_extract_epi32(__location, DELETION);
										EndOnMatchOrDeletion = 1U;
										Scores[MATCH] = Scores[DELETION];
									}
									int Align[2];
									char * restrict ptrState = &(alntrial->StateMemory[STATE_MEMORY_SIZE-1]);
									const char * const restrict StateLimit = &(alntrial->StateMemory[0]);
									*ptrState = '\0';
									
									/* Limit SoftClip Size according to kMinTagLenNotSoftClipped */
									assert(ReadLength >= kMinTagLenNotSoftClipped);
									int SoftClipBoundary = (TagsData[i].ValidQualityLen<kMinTagLenNotSoftClipped) ? kMinTagLenNotSoftClipped : (int) TagsData[i].ValidQualityLen;
									if (alntrial->revNmer) {
										SoftClipBoundary = ReadLength-SoftClipBoundary;
									}
									int iprf = ReadLength;
									unsigned int SoftClipMismatchCount = 0U;
									int ret;
									unsigned int MismatchCount  = 0U;
									unsigned char ContainsInDel = 0;
#ifdef BENCHMARKS
									unsigned long s1 = timer();
#endif
									unsigned int RequiredMMSize = 0U;
									unsigned int SoftClipEvictedMMSize = 0U;
									unsigned int RequiredCigSize = (TagsData[i].TagLen == TagsData[i].AlignLen) ? 1U : 2U;
									unsigned int SoftClipEvictedCigSize = 0U;
									
									#define __GET_ALIGNMENT_STATE__
									#ifdef PFTOOLS_SIMPLIFIED
									#include "pftools_std_simplified.h"
									#else
									#include "pftools_std.h"
									#endif
									#undef __GET_ALIGNMENT_STATE__
									
#ifdef BENCHMARKS
									const unsigned long se = timer();
									BenchmarkCounters[PROCESS_TOPIC_UNGAP_ALIGN].tics += se - s1;
									BenchmarkCounters[PROCESS_TOPIC_UNGAP].tics += se - s0;
									BenchmarkCounters[PROCESS_TOPIC_UNGAP_ALIGN].counter++;
									BenchmarkCounters[PROCESS_TOPIC_UNGAP].counter++;
#endif
#if MYTEST
									if (TagsData[i].ordinal == 0) {
										printf("SSE %u anchor %i OK: %i MisMatch: %u MM: %u CIG: %u SoftClip: %u Boudnary: %u genome: %u state: %s\n",
										       i, hitnum, ret, MismatchCount, RequiredMMSize, RequiredCigSize, SoftClipMismatchCount, SoftClipBoundary,
										       (unsigned int) ((intptr_t) (genome + Align[0]) - (intptr_t) &genome_tbl[BaseChrpos]), ptrState);										
									}
#endif
									if (ret == 0) {
										// POTENTIAL SOFTCLIP RESCUE
										if (MismatchCount > kMaxMismatch) {
											assert(MismatchCount >= SoftClipMismatchCount);
											alntrial->SoftClipMismatchToRemove = MismatchCount - kMaxMismatch;
											if (!(alntrial->revNmer)) {
												// check we will be able to encode, keep in mind that SoftClipEvictedCigSize in this orientation
												// cannot be totally removed, only  SoftClipEvictedCigSize-1. In addition, we need to keep space 
												// for the Softclip. Consequently, the maximum allowed size is kMaxEncodedCigar-1 with 
												//     kMaxEncodedCigar-1 >= RequiredCigSize - (SoftClipEvictedCigSize-1)
												assert(RequiredMMSize >= SoftClipEvictedMMSize);
												assert(RequiredCigSize >= SoftClipEvictedCigSize);
												if (   (RequiredMMSize-SoftClipEvictedMMSize) <= kMaxEncodedMismatch
														&& (RequiredCigSize-SoftClipEvictedCigSize) < (kMaxEncodedCigar-1) )
													MismatchCount -= SoftClipMismatchCount;
												else
													MismatchCount = kMaxMismatch+1;
											}
											else {
												// check we will be able to encode
												if (SoftClipEvictedMMSize <= kMaxEncodedMismatch && SoftClipEvictedCigSize <= kMaxEncodedCigar)
													MismatchCount = SoftClipMismatchCount;
												else 
													MismatchCount = kMaxMismatch+1;
											}
											if (MismatchCount < kMaxMismatch) MismatchCount = kMaxMismatch;
											alntrial->DecisionTreeHistory = kSoftClipRescued;
										}
										else if (RequiredMMSize > kMaxEncodedMismatch || RequiredCigSize > kMaxEncodedCigar) {
											continue;
										}

										if (MismatchCount <= kMaxMismatch) {
											assert(Align[0] >=0 ); 
											assert(((intptr_t) genome) + Align[0] > (intptr_t) BaseChrpos);
											alntrial->States               = ptrState;
											alntrial->DecisionTreeHistory |= kHasGoneThruAlignment | (ContainsInDel ? kContainsInDel : 0U);
											alntrial->score                = Scores[MATCH]; // we choose MATCH here to avoid come copy
											alntrial->chr_pos              = (unsigned int) ((intptr_t) (genome + Align[0]) - (intptr_t) &genome_tbl[BaseChrpos]);
											alntrial->MismatchCount        = MismatchCount ;
											alntrial->TotalRequiredMMCount = (unsigned char) RequiredMMSize;
											alntrial->TotalRequiredCigCount= (unsigned char) RequiredCigSize;
#ifndef NDEBUG
											alntrial->SoftClipEvictedCigSize = (unsigned char) SoftClipEvictedCigSize;
											alntrial->SoftClipEvictedMMSize  = (unsigned char) SoftClipEvictedMMSize;
											alntrial->GenomeSeekStart        = (unsigned int)  genomeCheckStart - ((uintptr_t) &genome_tbl[genomeCheckStart] & (15));
											alntrial->GenomSeekStop          = alntrial->GenomeSeekStart + GenomeLength;
#endif
											stats.Counter[Alignment]++;
											if (alntrial->DecisionTreeHistory & kSoftClipRescued) {
												stats.Counter[AlignmentRescued]++;
											}
#ifdef DEBUG_TAG_MAPPING_PROCESS
											printf("PFTOOLS: buf %u chunk %u trial %d/%d GOT IT on chr %d, offset %u with %u mismatches, score was %i\n",
													((BUFFERADDR*)data)->buf, ((BUFFERADDR*)data)->chunk, hitnum, kDesiredHitCnt,
													alntrial->chr, alntrial->chr_pos, MismatchCount, alntrial->score);
#endif
											goto PlaceThatHit;
										}
									}
#ifdef DEBUG_TAG_MAPPING_PROCESS
									printf("PFTOOLS: buf %u chunk %u trial %d/%d failed to align tag score was %i though!\n",
										((BUFFERADDR*)data)->buf, ((BUFFERADDR*)data)->chunk, hitnum, kDesiredHitCnt,Scores[MATCH]);
									fflush(stdout);
#endif
									stats.Counter[mismatch_overflow]++;
									continue;
								}
								//-----------------------------------------------------------------------------------------------------------------------
								// PLACE ALIGNMENT IN ORDER USING SCORE IF ANY FOUND
								//-----------------------------------------------------------------------------------------------------------------------
						PlaceThatHit:;
#ifdef DEBUG_TAG_MAPPING_PROCESS
								{
									printf("AS_IS: i=%d hitnum=%d shiftTagAgainstGenome Score=%d mismatch=%d chr=%d:%d  revNmer=%d revtag=%d matchpos=%d\n\n",
												i,hitnum,alntrial->score, alntrial->MismatchCount,
												alntrial->chr,
												alntrial->chr_pos,
												alntrial->revNmer,
												(rslt[i][hitnum].info >> 1) & 0x01,(rslt[i][hitnum].info >> 8) & 0xFF);
									fflush(stdout);
								}
#endif
								if ((hitnum != 0)) {
// 									__m512i __CurrentScore  = _mm512_extload_epi32(&(alntrial->score), _MM_UPCONV_EPI32_NONE, _MM_BROADCAST_1X16, 0);
// 									__m512i __CurrentOffset = _mm512_extload_epi32(&(KeptAnchorCount), _MM_UPCONV_EPI32_NONE, _MM_BROADCAST_1X16, 0);
// 									__m512i __vScores       = _mm512_load_epi32(&vScores[EvenOrOddOffset]);
// 									__m512i __vOffsets      = _mm512_load_epi32(&vOffsets[EvenOrOddOffset]);
// 									const __mmask16 Bigger  = _mm512_cmpge_epi32_mask(__CurrentScore, __vScores);
// 									__m512i __newScores     = _mm512_mask_mov_epi32(__vScores, _mm512_knot(Bigger), __CurrentScore);
// 									__m512i __newOffsets    = _mm512_mask_mov_epi32(__vOffsets, _mm512_knot(Bigger), __CurrentOffset);
// 									__newScores             = _mm512_mask_alignr_epi32(__vScores, Bigger, __newScores, __CurrentScore, 15);
// 									__newOffsets            = _mm512_mask_alignr_epi32(__vOffsets, Bigger, __newOffsets, __CurrentOffset, 15);
// 									_mm512_store_epi32(&vScores[EvenOrOddOffset], __newScores);
// 									_mm512_store_epi32(&vOffsets[EvenOrOddOffset], __newOffsets);
									
									int * const restrict ptrScore  = &vScores[EvenOrOddOffset];
									int * const restrict ptrOffset = &vOffsets[EvenOrOddOffset];
									const register int Score = alntrial->score;
									int k = kDesiredHitCnt-1;
									while (ptrScore[k] <= Score) {
										if (--k < 0) break;
										ptrScore[k+1]  = ptrScore[k];
										ptrOffset[k+1] = ptrOffset[k]; 
									}
									if (k < (kDesiredHitCnt-1)) {
										ptrScore[k+1]  = Score;
										ptrOffset[k+1] = KeptAnchorCount;
									}
								}
								else {
									vScores[EvenOrOddOffset]  = alntrial->score;
									vOffsets[EvenOrOddOffset] = KeptAnchorCount;  
								}
								KeptAnchorCount++;
							}
							//-------------------------------------------------------------------------------------------------------------------------
							// CHECK IF WE HAVE A HIT OR NOT : NO HIT
							//-------------------------------------------------------------------------------------------------------------------------
							else {
#ifdef DEBUG_TAG_MAPPING_PROCESS
			// 					printf("NO HIT for i=%d hitnum=%d\n\n",i,hitnum); fflush(stdout);
#endif
							}
						}
						//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
						// Loop on the anchor and test them, perfect and 1 mismatch anchor are moved to top
						//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

						//---------------------------------------------------------------------------------------------------------------------------
						// DECISION TREE
						//---------------------------------------------------------------------------------------------------------------------------
						if (!(i&1)) {
							EvenAnchorCount = KeptAnchorCount;
						}
						else { // ONLY WORK ON IMPAIR INDEX NUMBER
							/* 
							* At this stage we should have for each tag a list of (EvenAnchorCount or KeptAnchorCount - EvenAnchorCountpotential) hits
							* with its ordering (vOffsets [0..kDesiredHitCnt) and [kDesiredHitCnt .. 2*kDirectFromShift)) according to mismatch count
							* inside vScores ([0..kDesiredHitCnt) and [kDesiredHitCnt .. 2*kDirectFromShift)).
							*/
													
							unsigned int OddAnchorCount = KeptAnchorCount - EvenAnchorCount; 
							//==========================================================================================================================
							//
							//                                  ***   PREPARING POTENTIAL HIT PAIRS TO TEST  ***
							//
							//==========================================================================================================================
							unsigned long TotalRequiredSpace;
							char * const TempSoftClipString = (char*) IDM_Matrix;
							size_t SoftClipStringLength = 0UL;
							if (KeptAnchorCount > 0U) {
								int tcnt = 0;
								//--------------------------------------------------------------------------------------------------------------------------
								// Merge Odd and Even anchor list where combinations have to be on same chromosome but and on opposing strands.
								// Futhermore the distance in between is at most 2047
								//--------------------------------------------------------------------------------------------------------------------------
								int PairOrdering[kDesiredHitCnt];
								int PairScores[kDesiredHitCnt] __attribute__((aligned(16)));
								for (int i=0; i<kDesiredHitCnt; i++) PairScores[i] = NLOW;
								for (int EvenAnchor=0; EvenAnchor < EvenAnchorCount; EvenAnchor++) {
									InternalAlignment * const restrict EvenAlignment = &tagchoices[vOffsets[EvenAnchor]];
									for (int OddAnchor=0; OddAnchor < OddAnchorCount; OddAnchor++) {
										InternalAlignment * const restrict OddAlignment = &tagchoices[vOffsets[kDesiredHitCnt + OddAnchor]];
										if ( EvenAlignment->chr == OddAlignment->chr \
											&& EvenAlignment->revNmer != OddAlignment->revNmer \
											&& abs((int) EvenAlignment->chr_pos - (int) OddAlignment->chr_pos) <= 0x3FF ) {
											const int PairScore = EvenAlignment->score + OddAlignment->score;

											trials[tcnt].tag1 = EvenAlignment;
											trials[tcnt].tag2 = OddAlignment;

// 											__m512i __CurrentScore  = _mm512_extload_epi32(&score, _MM_UPCONV_EPI32_NONE, _MM_BROADCAST_1X16, 0);
// 											__m512i __CurrentOffset = _mm512_extload_epi32(&tcnt, _MM_UPCONV_EPI32_NONE, _MM_BROADCAST_1X16, 0);
// 											const __mmask16 Bigger  = _mm512_cmpge_epi32_mask(__CurrentScore, __Scores);
// 											__m512i __newScores     = _mm512_mask_mov_epi32(__Scores, _mm512_knot(Bigger), __CurrentScore);
// 											__m512i __newOffsets    = _mm512_mask_mov_epi32(__Ordering, _mm512_knot(Bigger), __CurrentOffset);
// 											__Scores                = _mm512_mask_alignr_epi32(__Scores, Bigger, __newScores, __CurrentScore, 15);
// 											__Ordering              = _mm512_mask_alignr_epi32(__Ordering, Bigger, __newOffsets, __CurrentOffset, 15);
											int k = kDesiredHitCnt-1;
											while (PairScores[k] <= PairScore) {
												if (--k < 0) break;
												PairScores[k+1]   = PairScores[k];
												PairOrdering[k+1] = PairOrdering[k]; 
											}
											if (k<(kDesiredHitCnt-1)) {
												PairScores[k+1]   = PairScore;
												PairOrdering[k+1] = tcnt;
											}
											stats.Counter[quickalign]++;
											
											tcnt++;
										}
									}
								}
								
								//------------------------------------------------------------------------------------------------------------------------
								// if the pairs on same chr were not good, then make sure we still try de novo nearby for each hit from
								// hit tag (ignoring the other tag). might just be the case for palindromes ?
								//------------------------------------------------------------------------------------------------------------------------
								if (tcnt == 0) {
									const int InitialKeptAnchorCount = KeptAnchorCount;
									const int InitialEvenAnchorCount = EvenAnchorCount;
									int IsItEvenTag = 1;
									//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
									// Loop on kept anchors
									//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
									for (int iAnchor=0; iAnchor<InitialKeptAnchorCount; iAnchor++) {
										/* Switch to second tag results */
										if (iAnchor >= InitialEvenAnchorCount) {
											IsItEvenTag = 0;
										}
										
										//----------------------------------------------------------------------------------------------------------------------
										// Check that anchors do not overlap, distance should be at least 100
										// Again this is DANGEROUS as we do not increase the region size for the remaining
										//----------------------------------------------------------------------------------------------------------------------
#ifdef PURIFY_ANCHORS
										{
											int Eject = 0;
											const long CurrentGenomePos = tagGenomePosition[iAnchor];
											#pragma noprefetch
											for (int OtherAnchor=0; OtherAnchor<iAnchor; OtherAnchor++) {
												Eject |= (labs(CurrentGenomePos - tagGenomePosition[OtherAnchor]) <= 300);
												if (Eject) break;
											}
											if (Eject) {
												stats.Counter[ForCompression]++;
											}
										}
#endif
										//----------------------------------------------------------------------------------------------------------------------
										// Prepare the side to seek, and its length
										//----------------------------------------------------------------------------------------------------------------------
										/* Get the storing location */ 
										InternalAlignment * const restrict alntrial = &tagchoices[KeptAnchorCount];
										assert(KeptAnchorCount < 4*kDesiredHitCnt);

										/* Get the base tag and compute the correct orientation and genome location*/
										const InternalAlignment * const restrict alnbase = &tagchoices[iAnchor];
										alntrial->chr = alnbase->chr;
										alntrial->revNmer = !(alnbase->revNmer);
										alntrial->chr_pos = alnbase->chr_pos;
										alntrial->DecisionTreeHistory = 0U;
#ifndef NDEBUG
										alntrial->AnchorId = alnbase->AnchorId + 128;
#endif
										
										const unsigned char * restrict ptrTag;
										size_t ReadLength;
										int SoftClipBoundary;
										if (alnbase->revNmer == 0) {
											register int itmp;
											if (IsItEvenTag) { // Base even and not reversed => look for reversed odd -> index of tags is i
												ReadLength = (size_t) TagsData[i].AlignLen;
												/* Check that the reverse was computed, if not do it */
												if (RevTagSpace[1][0] == '\0') {
													ReverseComplementOutputAligned(TagsData[i].Tag, &RevTagSpace[1][0], ReadLength);
												}
												ptrTag = &RevTagSpace[1][0];
												itmp = (int) TagsData[i].ValidQualityLen;
											}
											else { // Base Odd and not reversed => look for reversed even -> index is (i-1)
												ReadLength = (size_t) TagsData[i-1].AlignLen;
												/* Check that the reverse was computed, if not do it */
												if (RevTagSpace[0][0] == '\0') {
													ReverseComplementOutputAligned(TagsData[i-1].Tag,  &RevTagSpace[0][0], ReadLength);
												}
												ptrTag = &RevTagSpace[0][0];
												itmp = (int) TagsData[i-1].ValidQualityLen;
											}
											itmp = (itmp<kMinTagLenNotSoftClipped) ? kMinTagLenNotSoftClipped : itmp;
											SoftClipBoundary = ReadLength-itmp;
										}
										else {
											register int itmp;
											if (IsItEvenTag) { // Base even and reversed => look for standard odd -> index is i
												ptrTag = TagsData[i].Tag;
												itmp = (int) TagsData[i].ValidQualityLen;
												ReadLength = (size_t) TagsData[i].AlignLen;
											}
											else { // Base odd and reversed => look for standard even -> index is (i-1)
												ptrTag = TagsData[i-1].Tag;
												itmp = (int) TagsData[i-1].ValidQualityLen;
												ReadLength = (size_t) TagsData[i-1].AlignLen;
											}
											SoftClipBoundary =  (itmp<kMinTagLenNotSoftClipped) ? kMinTagLenNotSoftClipped : itmp;
											if (alntrial->chr_pos > (GenomeChunkSize - kTagSize)) {
												alntrial->chr_pos -= (GenomeChunkSize - kTagSize);
											}
											else {
												alntrial->chr_pos = 1;
											}
										}
										
										//----------------------------------------------------------------------------------------------------------------------
										// Perform the alignment (old DoUngap), place potential new alignment into correct lane
										//----------------------------------------------------------------------------------------------------------------------
										{
											unsigned int genomeCheckStart = (virtchr[alntrial->chr].chr << 28) \
																										+ virtchr[alntrial->chr].offset + alntrial->chr_pos;
											const unsigned char * restrict const genome = (unsigned char*) ((uintptr_t) &genome_tbl[genomeCheckStart] & ~(15));
											const size_t GenomeLength = (size_t) GenomeChunkSize + kTagSize + ((uintptr_t) &genome_tbl[genomeCheckStart] & (15));
											const size_t MatrixLD     = 1 + ReadLength;
											
#ifdef DEBUG_TAG_MAPPING_PROCESS
											printf("DO UNGAP SEARCH tag len: %lu chr:%d pos:%u, search length: %lu genomepos=%u\n", ReadLength,
														alntrial->chr, alntrial->chr_pos, GenomeLength, genomeCheckStart);
#endif
											
#ifdef BENCHMARKS
											unsigned long s0 = timer();
#endif
											__m128i __BestScore, __location;
											// FILL IN MATRIX
											#define __FILL_MATRIX__
											#ifdef PFTOOLS_SIMPLIFIED
											#include "pftools_std_simplified.h"
											#else
											#include "pftools_std.h"
											#endif
											#undef __FILL_MATRIX__
												
#ifdef BENCHMARKS
											BenchmarkCounters[PROCESS_TOPIC_CREATE_MATRIX].tics += timer() - s0;
											BenchmarkCounters[PROCESS_TOPIC_CREATE_MATRIX].counter++;
#endif
											__BestScore = _mm_srai_epi32(__BestScore, SCORE_SHIFT);
											int tmpScores[4] __attribute__((aligned(16)));
											_mm_store_si128((__m128i*) tmpScores, __BestScore);
											
											int iseq;
											unsigned int EndOnMatchOrDeletion; // = MorD[num];
											if (tmpScores[MATCH] >= tmpScores[DELETION]) {
												iseq = _mm_extract_epi32(__location, MATCH);
												EndOnMatchOrDeletion = 0U;
											}
											else {
												iseq = _mm_extract_epi32(__location, DELETION);
												EndOnMatchOrDeletion = 1U;
												tmpScores[MATCH] = tmpScores[DELETION];
											}
	
											int Align[2];
											int ret;
											unsigned int MismatchCount = 0U;
				// 							unsigned int InsertionCount = 0U;
											unsigned int SoftClipMismatchCount = 0U;
																						
											int iprf = ReadLength;
											char * restrict ptrState = &(alntrial->StateMemory[STATE_MEMORY_SIZE-1]);
											const char * const restrict StateLimit = &(alntrial->StateMemory[0]);
											*ptrState = '\0';
#ifdef BENCHMARKS
											unsigned long s1 = timer();
#endif
											unsigned int RequiredMMSize = 0U;
											unsigned int SoftClipEvictedMMSize = 0U;
											unsigned int RequiredCigSize = (TagsData[i].TagLen == TagsData[i].AlignLen) ? 1U : 2U;
											unsigned int SoftClipEvictedCigSize = 0U;
											unsigned char ContainsInDel = 0;
											#define OUT OUT2 
											#define DONE DONE2 
											#define __GET_ALIGNMENT_STATE__
											#ifdef PFTOOLS_SIMPLIFIED
											#include "pftools_std_simplified.h"
											#else
											#include "pftools_std.h"
											#endif
											#undef __GET_ALIGNMENT_STATE__
											#undef OUT
											#undef DONE
											
#ifdef BENCHMARKS
											const unsigned long se = timer();
											BenchmarkCounters[PROCESS_TOPIC_GET_ALIGNMENT].tics += se - s1;
											BenchmarkCounters[PROCESS_TOPIC_GET_ALIGNMENT].counter++;
#endif
#if MYTEST
											if (TagsData[i].ordinal == 0) {
												printf("SSE %u anchor %i OK: %i MisMatch: %u MM: %u CIG: %u SoftClip: %u Boundary: %u %s\n",
												       i, iAnchor, ret, MismatchCount, RequiredMMSize, RequiredCigSize, SoftClipMismatchCount, SoftClipBoundary, ptrState);
                      }
#endif
				// 							printf("%u\t%u\n%s\n", MismatchCount, SoftClipMismatchCount, ptrState);
											if (ret == 0 ) {
												// POTENTIAL SOFTCLIP RESCUE
												if (MismatchCount > kMaxMismatch) {
													assert(MismatchCount >= SoftClipMismatchCount);
													alntrial->SoftClipMismatchToRemove = MismatchCount - kMaxMismatch;
													if (!(alntrial->revNmer)) {
														// check we will be able to encode, keep min mind that SoftClipEvictedCigSize in this orientation
														// cannot be totally removed, only  SoftClipEvictedCigSize-1. In addition, we need to keep space 
														// for the Softclip. Consequently, the maximum allowed size is kMaxEncodedCigar-1 with 
														//     kMaxEncodedCigar-1 >= RequiredCigSize - (SoftClipEvictedCigSize-1)
														assert(RequiredMMSize >= SoftClipEvictedMMSize);
														assert(RequiredCigSize >= SoftClipEvictedCigSize);
														if (   (RequiredMMSize-SoftClipEvictedMMSize) <= kMaxEncodedMismatch
																&& (RequiredCigSize-SoftClipEvictedCigSize) < (kMaxEncodedCigar-1) )
																MismatchCount -= SoftClipMismatchCount;
															else 
																MismatchCount = kMaxMismatch+1;
													}
													else {
														// check we will be able to encode
														if (SoftClipEvictedMMSize <= kMaxMismatch && SoftClipEvictedCigSize <= kMaxEncodedCigar)
															MismatchCount = SoftClipMismatchCount;
														else 
															MismatchCount = kMaxMismatch+1;
													}
													if (MismatchCount < kMaxMismatch) MismatchCount = kMaxMismatch;
													alntrial->DecisionTreeHistory = kSoftClipRescued;
												}
												else if (RequiredMMSize > kMaxEncodedMismatch || RequiredCigSize > kMaxEncodedCigar) {
													continue;
												}
												
												if (MismatchCount <= kMaxMismatch) {
#ifdef DEBUG_TAG_MAPPING_PROCESS
													if (Align[0] < 0 ) 
														printf("PFTOOLS : Alignment start is negative %i !!!\n", Align[0]);
#endif			
													alntrial->DecisionTreeHistory |= kHasGoneThruAlignment | kLargeAlignmentSearch | (ContainsInDel ? kContainsInDel : 0U);
													alntrial->score                = tmpScores[MATCH]; // we choose MATCH here to avoid come copy
													alntrial->chr_pos              = (int) alntrial->chr_pos + (Align[0] - (int) ((uintptr_t) &genome_tbl[genomeCheckStart] & (15)));
													if(((int) alntrial->chr_pos + (Align[0] - (int) ((uintptr_t) &genome_tbl[genomeCheckStart] & (15))))<=0) {
														printf("Wrong alignment : genome address 0x%16.16lx, alignment start at %i\n",
																	(uintptr_t) &genome_tbl[genomeCheckStart], Align[0]);
											
													}
													alntrial->MismatchCount         = MismatchCount;
													alntrial->States                = ptrState;
													alntrial->TotalRequiredMMCount  = (unsigned char) RequiredMMSize;
													alntrial->TotalRequiredCigCount = (unsigned char) RequiredCigSize;
#ifndef NDEBUG
													alntrial->SoftClipEvictedCigSize = (unsigned char) SoftClipEvictedCigSize;
													alntrial->SoftClipEvictedMMSize  = (unsigned char) SoftClipEvictedMMSize;
													alntrial->GenomeSeekStart        = (unsigned int)  genomeCheckStart - ((uintptr_t) &genome_tbl[genomeCheckStart] & (15));
													alntrial->GenomSeekStop          = alntrial->GenomeSeekStart + GenomeLength;
#endif
													
													stats.Counter[LargeAlignment]++;
													if (alntrial->DecisionTreeHistory & kSoftClipRescued) {
														stats.Counter[LargeAlignmentRescued]++;
													}
												
													const int PairScore = alnbase->score + alntrial->score;
													if (IsItEvenTag) {
														trials[tcnt].tag1 = (InternalAlignment*) alnbase;
														trials[tcnt].tag2 = alntrial;
														EvenOrOddOffset = kDesiredHitCnt;
														EvenAnchorCount++;
													}
													else {
														trials[tcnt].tag1 = alntrial;
														trials[tcnt].tag2 = (InternalAlignment*) alnbase;
														EvenOrOddOffset = 0U;
														OddAnchorCount++;
													}
												
													{
														int * const restrict ptrScore  = &vScores[EvenOrOddOffset];
														int * const restrict ptrOffset = &vOffsets[EvenOrOddOffset];
														const register int Score = alntrial->score;
														int k = kDesiredHitCnt-1;
														while (ptrScore[k] <= Score) {
															if (--k < 0) break;
															ptrScore[k+1]  = ptrScore[k];
															ptrOffset[k+1] = ptrOffset[k]; 
														}
														if (k<(kDesiredHitCnt-1)) {
															ptrScore[k+1]  = Score;
															ptrOffset[k+1] = KeptAnchorCount;
														}
													}
													KeptAnchorCount++;
													{
														int k = kDesiredHitCnt-1;
														while (PairScores[k] <= PairScore) {
															if (--k < 0) break;
															PairScores[k+1]   = PairScores[k];
															PairOrdering[k+1] = PairOrdering[k]; 
														}
														if (k<(kDesiredHitCnt-1)) {
															PairScores[k+1]   = PairScore;
															PairOrdering[k+1] = tcnt;
														}
													}
													tcnt++;
												}
											}
										}
									}
									//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
									// Loop on kept anchors
									//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
								}
								
								//--------------------------------------------------------------------------------------------------------------------------
								// See if something is left within the allowed number of mismatch
								//--------------------------------------------------------------------------------------------------------------------------
								Trial_t BestPair = { .tag1=NULL, .tag2=NULL };
								Trial_t AlternatePair = { .tag1=NULL, .tag2=NULL };
								if (tcnt > 0) {
									/* Get the best pair */
									BestPair = trials[PairOrdering[0]];
									if (PairScores[1] == PairScores[0]) {
										/* Check we do not get the same location */
										const _Bool AA = trials[PairOrdering[1]].tag1->chr==BestPair.tag1->chr && trials[PairOrdering[1]].tag1->chr_pos==BestPair.tag1->chr_pos;
										const _Bool BB = trials[PairOrdering[1]].tag2->chr==BestPair.tag2->chr && trials[PairOrdering[1]].tag2->chr_pos==BestPair.tag2->chr_pos;
										const _Bool AB = trials[PairOrdering[1]].tag1->chr==BestPair.tag2->chr && trials[PairOrdering[1]].tag1->chr_pos==BestPair.tag2->chr_pos;
										const _Bool BA = trials[PairOrdering[1]].tag2->chr==BestPair.tag1->chr && trials[PairOrdering[1]].tag2->chr_pos==BestPair.tag1->chr_pos;

										if (AA || AB || BB || BA) {
											const long int space0 = labs((long int) BestPair.tag1->chr_pos - (long int) BestPair.tag2->chr_pos);
											const long int space1 = labs((long int) trials[PairOrdering[1]].tag1->chr_pos - (long int) trials[PairOrdering[1]].tag2->chr_pos);
											if (space1 < space0) {
												BestPair =  trials[PairOrdering[1]];
											}
										}
										else {
											AlternatePair = trials[PairOrdering[1]];
											stats.Counter[Alternate]++;
										}
									}
									
				// 				printf("BestPair found with Even score %i and Odd score %i\n", BestPair.tag1->score, BestPair.tag2->score);
									goto DoOutput;
								}
								
								//--------------------------------------------------------------------------------------------------------------------------
								// HALFMAP & CHIMERA : CONDITION IS 1 MISMATCH AT MOST, NO CIGAR, THUS NO TRIM
								//--------------------------------------------------------------------------------------------------------------------------
								if (TagsData[i].TagLen == TagsData[i].AlignLen) {
									/* Loop on the Even side to see if we pass the 1 mismatch threshold, remember those are ordered */
									if (EvenAnchorCount > kDesiredHitCnt) EvenAnchorCount = kDesiredHitCnt;
									for (int iEvenAnchor=0; iEvenAnchor<EvenAnchorCount; iEvenAnchor++) {
										InternalAlignment * const restrict Aln = &tagchoices[vOffsets[iEvenAnchor]];
										if (Aln->MismatchCount <= 1 && !(Aln->DecisionTreeHistory & kContainsInDel)) {
											Aln->DecisionTreeHistory |= kBestOfTags;
											BestPair.tag1 = Aln;
											stats.Counter[BestOfEach]++;
											break;
										}
									}
									
									/* Loop on the Odd side to see if we pass the 1 mismatch threshold, remember those are ordered */
									if (OddAnchorCount > kDesiredHitCnt) OddAnchorCount = kDesiredHitCnt;
									for (int iOddAnchor=0; iOddAnchor<OddAnchorCount; iOddAnchor++) {
										InternalAlignment * const restrict Aln = &tagchoices[vOffsets[kDesiredHitCnt+iOddAnchor]];
										if (Aln->MismatchCount <= 1 && !(Aln->DecisionTreeHistory & kContainsInDel)) {
											Aln->DecisionTreeHistory |= kBestOfTags;
											BestPair.tag2 = Aln;
											stats.Counter[BestOfEach]++;
											break;
										}
									}
								}
								
								//--------------------------------------------------------------------------------------------------------------------------
								// DUMP WHAT WE FOUND, unmapped data if none
								//--------------------------------------------------------------------------------------------------------------------------
							DoOutput:;
								stats.Counter[HalfMapped]++;
								/* First compute the space required for the output */
								outputHeader.Ordinal1          = TagsData[i-1].ordinal;
								outputHeader.Ordinal2          = TagsData[i].ordinal;
								outputHeader.TagLength1        = TagsData[i-1].TagLen;
								outputHeader.TagLength2        = TagsData[i].TagLen;
								outputHeader.HeaderLength1     = TagsData[i-1].HeaderLen;
								outputHeader.HeaderLength2     = TagsData[i].HeaderLen;
								outputHeader.StreamedNext      = 0;
								outputHeader.nAlignmentResults = 0;
								TotalRequiredSpace = sizeof(outputPairBlockHeader_t) + TagsData[i-1].HeaderLen + TagsData[i-1].TagLen \
								                   + TagsData[i].HeaderLen + TagsData[i].TagLen;
								
#define PopulateResult(PairTagPtr, ResAln, TagsData, RevRealTagStringPtr, SoftClipMask) {\
	const unsigned short HardClipSize = (TagsData.TagLen - TagsData.AlignLen);\
	if (!(PairTagPtr->DecisionTreeHistory & kMMEncoded)) {\
		const unsigned char * const restrict ptrTag = (PairTagPtr->revNmer) ? (const unsigned char*) RevRealTagStringPtr : (const unsigned char*) TagsData.Tag;\
		__assume_aligned(ptrTag, 64);\
		const unsigned int SoftClipSpace = EncodeInternalAlignment(ptrTag, PairTagPtr, &(ResAln), &TagsData);\
		if (SoftClipSpace>0U) {\
			TotalRequiredSpace += SoftClipSpace;\
			memcpy(&TempSoftClipString[SoftClipStringLength], TagsData.Tag + (TagsData.TagLen - SoftClipSpace), SoftClipSpace);\
			SoftClipStringLength += SoftClipSpace;\
			outputHeader.StreamedNext |= SoftClipMask;\
		}\
	}\
	else {\
		ResAln.score         = PairTagPtr->score;\
		ResAln.chr           = PairTagPtr->chr;\
		ResAln.chr_pos       = PairTagPtr->chr_pos;\
		ResAln.MismatchCount = PairTagPtr->MismatchCount;\
		ResAln.revNmer       = PairTagPtr->revNmer;\
		if (HardClipSize == (unsigned short) 0) {\
			ResAln.DecisionTreeHistory = PairTagPtr->DecisionTreeHistory;\
			ResAln.Cigar[0]            = kCIGARTERMINATOR;\
			ResAln.CigarLength         = 1;\
		}\
		else {\
			ResAln.DecisionTreeHistory = PairTagPtr->DecisionTreeHistory | kTrimmed;\
			if (PairTagPtr->revNmer) {\
				ResAln.Cigar[0] = (unsigned char) HardClipSize;\
				ResAln.Cigar[1] = 'S';\
				ResAln.Cigar[2] = (unsigned char) TagsData.AlignLen;\
				ResAln.Cigar[3] = 'M';\
				ResAln.chr_pos -= (unsigned int) HardClipSize;\
			}\
			else {\
				ResAln.Cigar[0] = (unsigned char) TagsData.AlignLen;\
				ResAln.Cigar[1] = 'M';\
				ResAln.Cigar[2] = (unsigned char) HardClipSize;\
				ResAln.Cigar[3] = 'S';\
			}\
			ResAln.Cigar[4]     = kCIGARTERMINATOR;\
			ResAln.CigarLength  = 5;\
			TotalRequiredSpace += (size_t) HardClipSize;\
			memcpy(&TempSoftClipString[SoftClipStringLength], TagsData.Tag + TagsData.AlignLen, (size_t) HardClipSize);\
			SoftClipStringLength += (size_t) HardClipSize;\
			outputHeader.StreamedNext |= SoftClipMask;\
		}\
		/* Copy Mismatch string hidden in StateMemory by EncodeMM */\
		unsigned int k = 0;\
		const unsigned char * const restrict ptrState = (unsigned char*) &(PairTagPtr->StateMemory[0]);\
		do {\
			ResAln.Mismatch[k] = ptrState[k];\
			assert(k<(kMAXmismatchSize));\
		} while (ptrState[k++] != kMMTERMINATOR);\
		ResAln.MismatchLength = (unsigned char) k;\
	}\
	/* Correct mismatches in the case hardClip and reverse */\
	if ( ResAln.Mismatch[0] != kMMTERMINATOR && PairTagPtr->revNmer && (HardClipSize > (unsigned short) 0)) {\
		unsigned char * MismatchPtr = &(ResAln.Mismatch[0]);\
		do {\
			*MismatchPtr += HardClipSize;\
			MismatchPtr += 2;\
		} while (*MismatchPtr != kMMTERMINATOR);\
	}\
}
								outputHeader.nAlignmentResults = 0;
								if (BestPair.tag1 != NULL) {
#ifdef EXPORT_STATES
									if (! (BestPair.tag1->DecisionTreeHistory & kMMEncoded)) {
										int t=0;
										while(t<strlen(BestPair.tag1->States)) {
											localResults[outputHeader.nAlignmentResults].States[t] = BestPair.tag1->States[t++];
										}
										while(t<(3*kTagSize/2)) localResults[outputHeader.nAlignmentResults].States[t++] = '\0';
									}
									else {
										memset(localResults[outputHeader.nAlignmentResults].States, 0, 3*kTagSize/2);
									}
#endif
									PopulateResult(BestPair.tag1, localResults[outputHeader.nAlignmentResults], TagsData[i-1], &RevTagSpace[0][0], 0b00010000);
#ifndef NDEBUG
									localResults[outputHeader.nAlignmentResults].GenomeSeekStart = BestPair.tag1->GenomeSeekStart;
									localResults[outputHeader.nAlignmentResults].GenomSeekStop = BestPair.tag1->GenomSeekStop;
									localResults[outputHeader.nAlignmentResults].AnchorID = BestPair.tag1->AnchorId;
#endif
									outputHeader.nAlignmentResults++;
									TotalRequiredSpace += sizeof(AlignmentResults_t);
								}
								else {
									outputHeader.StreamedNext |= 0x1;
									TotalRequiredSpace += TagsData[i-1].TagLen;
#ifndef NDEBUG
									localResults[outputHeader.nAlignmentResults].GenomeSeekStart = 0U;
									localResults[outputHeader.nAlignmentResults].GenomSeekStop = 0U;
									localResults[outputHeader.nAlignmentResults].AnchorID = -1;
#endif
#ifdef EXPORT_STATES
									memset(localResults[outputHeader.nAlignmentResults].States, 0, 3*kTagSize/2);
#endif
								}
								if (BestPair.tag2 != NULL) {
#ifdef EXPORT_STATES
									if (! (BestPair.tag2->DecisionTreeHistory & kMMEncoded)) {
										int t=0;
										while(t<strlen(BestPair.tag2->States)) {
															localResults[outputHeader.nAlignmentResults].States[t] = BestPair.tag2->States[t++];
										}
										while(t<(3*kTagSize/2)) localResults[outputHeader.nAlignmentResults].States[t++] = '\0';
									}
									else {
										memset(localResults[outputHeader.nAlignmentResults].States, 0, 3*kTagSize/2);
									}
#endif
									PopulateResult(BestPair.tag2, localResults[outputHeader.nAlignmentResults],TagsData[i], &RevTagSpace[1][0], 0b00100000);
#ifndef NDEBUG
									localResults[outputHeader.nAlignmentResults].GenomeSeekStart = BestPair.tag2->GenomeSeekStart;
									localResults[outputHeader.nAlignmentResults].GenomSeekStop = BestPair.tag2->GenomSeekStop;
									localResults[outputHeader.nAlignmentResults].AnchorID = BestPair.tag2->AnchorId;
#endif
									outputHeader.nAlignmentResults++;
									TotalRequiredSpace += sizeof(AlignmentResults_t);
								}
								else {
									outputHeader.StreamedNext |= 0x2;
									TotalRequiredSpace += TagsData[i].TagLen;
#ifndef NDEBUG
									localResults[outputHeader.nAlignmentResults].GenomeSeekStart = 0U;
									localResults[outputHeader.nAlignmentResults].GenomSeekStop = 0U;
									localResults[outputHeader.nAlignmentResults].AnchorID = -1;
#endif
#ifdef EXPORT_STATES
									memset(localResults[outputHeader.nAlignmentResults].States, 0, 3*kTagSize/2);
#endif
								}

								if (AlternatePair.tag1 != NULL) {
#ifdef EXPORT_STATES
									if (! (AlternatePair.tag1->DecisionTreeHistory & kMMEncoded)) {
										int t=0;
										while(t<strlen(AlternatePair.tag1->States)) {
											localResults[outputHeader.nAlignmentResults].States[t] = AlternatePair.tag1->States[t++];
										}
										while(t<(3*kTagSize/2)) localResults[outputHeader.nAlignmentResults].States[t++] = '\0';
									}
									else {
										memset(localResults[outputHeader.nAlignmentResults].States, 0, 3*kTagSize/2);
									}
#endif
									PopulateResult(AlternatePair.tag1, localResults[outputHeader.nAlignmentResults], TagsData[i-1], &RevTagSpace[0][0], 0b01000000);
#ifndef NDEBUG
									localResults[outputHeader.nAlignmentResults].GenomeSeekStart = AlternatePair.tag1->GenomeSeekStart;
									localResults[outputHeader.nAlignmentResults].GenomSeekStop = AlternatePair.tag1->GenomSeekStop;
									localResults[outputHeader.nAlignmentResults].AnchorID = AlternatePair.tag1->AnchorId;
#endif
									outputHeader.nAlignmentResults++;
									TotalRequiredSpace += sizeof(AlignmentResults_t);
								}
								
								if (AlternatePair.tag2 != NULL) {
#ifdef EXPORT_STATES
									if (! (AlternatePair.tag2->DecisionTreeHistory & kMMEncoded)) {
										int t=0;
										while(t<strlen(AlternatePair.tag2->States)) {
															localResults[outputHeader.nAlignmentResults].States[t] = AlternatePair.tag2->States[t++];
										}
										while(t<(3*kTagSize/2)) localResults[outputHeader.nAlignmentResults].States[t++] = '\0';
									}
									else {
										memset(localResults[outputHeader.nAlignmentResults].States, 0, 3*kTagSize/2);
									}
#endif
									PopulateResult(AlternatePair.tag2, localResults[outputHeader.nAlignmentResults], TagsData[i], &RevTagSpace[1][0], 0b10000000);
#ifndef NDEBUG
									localResults[outputHeader.nAlignmentResults].GenomeSeekStart = AlternatePair.tag2->GenomeSeekStart;
									localResults[outputHeader.nAlignmentResults].GenomSeekStop = AlternatePair.tag2->GenomSeekStop;
									localResults[outputHeader.nAlignmentResults].AnchorID = AlternatePair.tag2->AnchorId;
#endif
									outputHeader.nAlignmentResults++;
									TotalRequiredSpace += sizeof(AlignmentResults_t);
								}

#undef PopulateResult
							}
							else {
								//-------------------------------------------------------------------------------------------------------------------------
								// IS IT HOPELESS AS NO HIT WERE FOUND ON EACH TAG, move to next pair if so
								//-------------------------------------------------------------------------------------------------------------------------
								/* First compute the space required for the output */
								outputHeader.Ordinal1 = TagsData[i-1].ordinal;
								outputHeader.Ordinal2 = TagsData[i].ordinal;
								outputHeader.TagLength1 = TagsData[i-1].TagLen;
								outputHeader.TagLength2 = TagsData[i].TagLen;
								outputHeader.HeaderLength1 = TagsData[i-1].HeaderLen;
								outputHeader.HeaderLength2 = TagsData[i].HeaderLen;
								outputHeader.nAlignmentResults = 0;
								outputHeader.StreamedNext = 0x3;
								
								TotalRequiredSpace = sizeof(outputPairBlockHeader_t) + TagsData[i-1].HeaderLen + 2*TagsData[i-1].TagLen \
								                   + TagsData[i].HeaderLen + 2*TagsData[i].TagLen;
								stats.Counter[hopeless]++;
							}
							
							/* Extend the required space to account potential alignment for reading back direclty, thus avoiding extra memcpy */ 
							TotalRequiredSpace += 8;
							TotalRequiredSpace = (TotalRequiredSpace + (__alignof__(outputPairBlockHeader_t) - 1)) & ~(__alignof__(outputPairBlockHeader_t) - 1);
							
							/* Get output buffer */
							unsigned long int SpaceUsedSoFar;
							int CurrentOutputIndex;
							while (1) {
								while ((CurrentOutputIndex = atomic_load(OutIndex)) < 0) {
									CurrentOutputIndex = -1;
									if (atomic_compare_exchange_strong(OutIndex,  &CurrentOutputIndex, -2)) {
//  										printf("Thread %u seeking new space\n", MyID);fflush(stdout);
										while (1) {
											bsem_wait(data->space_q.has_items);
											pthread_mutex_lock(&(data->space_q.rwmutex));
											cpuSpace_t * const restrict newSpace = (cpuSpace_t*) jobqueue_pull(&(data->space_q));
											pthread_mutex_unlock(&(data->space_q.rwmutex));
											if ( newSpace != NULL ) {
												assert(newSpace->index >= 0); assert(newSpace->index < data->AlignerBlock->nOutputBlock);
												assert(atomic_load(&(newSpace->Stored)) == 0UL );
												atomic_store(&(newSpace->Used), 0UL);
												atomic_store(OutIndex, newSpace->index);
												CurrentOutputIndex = newSpace->index;
//  												printf("Thread %u pulled space %i\n", MyID, newSpace->index);fflush(stdout);
												break;
											}
										}
									}
								}
								
								assert(TotalRequiredSpace > 0 && TotalRequiredSpace <= OutputBufferSize);
								assert(CurrentOutputIndex >= 0); assert(CurrentOutputIndex < data->AlignerBlock->nOutputBlock);	
								SpaceUsedSoFar = atomic_fetch_add(&(OutputBuffers[CurrentOutputIndex].Used), TotalRequiredSpace);
								if ((SpaceUsedSoFar + TotalRequiredSpace) >= OutputBufferSize) {
									/* WARNING:
										* If still using the current buffer and we are the one just passing over the limit,
										* advise other to fetch new buffer, then dump data 
										*/
									if (SpaceUsedSoFar < OutputBufferSize) {
										int itmp = CurrentOutputIndex;
										atomic_compare_exchange_strong(OutIndex, &itmp, -1);
										{
//  											printf("Appolo Thread %u not enough space in %i, dumping %lu bytes\n", MyID, CurrentOutputIndex, SpaceUsedSoFar);
// 											fflush(stdout);	
											/* Dump data */
											assert(data->space_q.has_items != NULL);
											SendBufferToHost(ToBeWritten_q, &OutputBuffers[CurrentOutputIndex], &(data->space_q), SpaceUsedSoFar);
										}
									}
								}
								else 
									break;
							}
							
							/* Now write down the solutions */
							char * restrict outPtr = OutputBuffers[CurrentOutputIndex].OutBufferPtr + SpaceUsedSoFar;
							memcpy(outPtr, "NEW_ALN\0", 8); outPtr += 8;
							memcpy(outPtr, &outputHeader, sizeof(struct outputPairBlockHeader));
							outPtr += sizeof(struct outputPairBlockHeader);
							if ( outputHeader.nAlignmentResults ) {
								memcpy(outPtr, &localResults[0], outputHeader.nAlignmentResults*sizeof(AlignmentResults_t));
								outPtr += outputHeader.nAlignmentResults*sizeof(AlignmentResults_t);
							}
							memcpy(outPtr, TagsData[i-1].Header, TagsData[i-1].HeaderLen);
							outPtr += TagsData[i-1].HeaderLen;
							memcpy(outPtr, TagsData[i].Header, TagsData[i].HeaderLen);
							outPtr += TagsData[i].HeaderLen;
							memcpy(outPtr, TagsData[i-1].Quality, TagsData[i-1].TagLen);
							outPtr += TagsData[i-1].TagLen;
							memcpy(outPtr, TagsData[i].Quality, TagsData[i].TagLen);
							outPtr += TagsData[i].TagLen;
							if (outputHeader.StreamedNext & 0x1) {
									memcpy(outPtr, TagsData[i-1].Tag, TagsData[i-1].TagLen);
									outPtr += TagsData[i-1].TagLen;
							}
							if (outputHeader.StreamedNext & 0x2) {
									memcpy(outPtr, TagsData[i].Tag, TagsData[i].TagLen);
									outPtr += TagsData[i].TagLen;
							}
							if (SoftClipStringLength) {
								memcpy(outPtr, TempSoftClipString, SoftClipStringLength);
							}
							
							atomic_fetch_add(&(OutputBuffers[CurrentOutputIndex].Stored), TotalRequiredSpace);
							
						} // ONLY WORK ON IMPAIR INDEX NUMBER
					}
					//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
					// Process all tags within the batch
					//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
					
					/* 
					 * =========================================================================================================================
					 * Done aligning partial part 
					 * =========================================================================================================================
					 */
					/* Inform store is complete */
					atomic_fetch_add(&(Jobs[CurrentInputIndex].Stored), (unsigned long int) WorkLoad);

					/* Ask for new job */
					if (IamLast) {
						/* Put back the space of Reader and Decoder to their respective stack */
						SendSpaceBufferBack(data->AlignerBlock, &Jobs[CurrentInputIndex]);
						
//  						printf("Appolo Thread %u job results fully stored, requesting new job for buffer %u\n", MyID, CurrentInputIndex); fflush(stdout);
						if(recvBufferFromHost(AlignPool, &Jobs[CurrentInputIndex]) != 0) {
// 							printf("Appolo Thread %u got end of task\n", MyID); fflush(stdout);
							data->threads_keepalive = 0;
							/* Place the space back into its queue, setting the count to 0 to avoid flushing later */
							Jobs[CurrentInputIndex].TagCnt = 0U;
							goto LoopDone;
						}
						else {
							/* Replace the jobs in the queue */
							assert(Jobs[CurrentInputIndex].TagCnt > 0U);
							int Test = 0;
							while (!atomic_compare_exchange_strong(&(data->job_q.rwmutex), &Test, -1)) {
								Test = 0;
								__asm__ __volatile__ ("pause" ::: "memory");
							}
							pushToQueue(&(data->job_q), (Task_t*) &Jobs[CurrentInputIndex]);
							atomic_store(&(data->job_q.rwmutex), 0);
// 							printf("New batch in buffer %u with %u tags\n", CurrentInputIndex, Jobs[CurrentInputIndex].TagCnt);
						}
					}
				}
			} while (CurrentInputIndex == atomic_load(JobIndex));
// 			printf("Thread 0x%lx done with job %i\n", Me, CurrentInputIndex);fflush(stdout);			
		}
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//  Load next if index is empty
		int Empty = -1;
		if (atomic_compare_exchange_strong(JobIndex,  &Empty, -2)) {
//  			printf("Appolo Thread %u seeking new job\n", MyID);fflush(stdout);
			do  { //while (data->threads_keepalive) {
				int Test = 0;
				while (!atomic_compare_exchange_strong(&(data->job_q.rwmutex), &Test, -1)) {
					Test = 0;
					__asm__ __volatile__("pause" ::: "memory");
				}
				cpuAlignerJob_t * const restrict newJob = (cpuAlignerJob_t*) pullFromQueue(&(data->job_q));
				atomic_store(&(data->job_q.rwmutex), 0);
				
				if ( newJob != NULL ) {
					assert( newJob->TagsData != NULL); 
					atomic_store(JobIndex, newJob->index);
					CurrentInputIndex = newJob->index;
//  					printf("Appolo Thread %u pulled job %i\n", MyID, newJob->index);fflush(stdout);
// 	 				printf("Appolo Thread %u working on job %i\n", MyID, CurrentInputIndex);fflush(stdout);
					goto Compute;
				}
			} while (data->threads_keepalive);
			atomic_store(JobIndex, -1);
			break;
		}
	}
	
	LoopDone: ;

	/********************************************************************************************************************************/
	/*                                                    THREAD INACTIVE                                                           */
	/********************************************************************************************************************************/
	//atomic_fetch_sub(&data->num_threads_alive, 1);
		
	/********************************************************************************************************************************/
	/*                                                          POST                                                                */
	/********************************************************************************************************************************/
	if (!PosixAllocation) {
		/* Release Huge TLB matrix */
		munmap(IDM_Matrix, 2*HUGE_PAGE_SIZE);
	}
	else {
		free(IDM_Matrix);
	}
	/* Release Work space */
	_mm_free(WorkSpace);
	
	/* Pass process statistics to the receiver */	
	for (int i=0; i<4; i++) data->stats[MyID].Vector[i]  = stats.Vector[i];

	
#ifdef BENCHMARKS
	const int id = buf*gSlaveCnt+chunk;
	for (int i=0; i<MAX_PROCESS_TOPICS; i++) {
		ProcessBenchmarks[id].topics[i].tics    = BenchmarkCounters[i].tics;
		ProcessBenchmarks[id].topics[i].counter = BenchmarkCounters[i].counter;
	}
	/* populate usage */
	if (!getrusage(RUSAGE_THREAD, &(ProcessBenchmarks[id].usage))) {
		fprintf(stderr, "Processer %i of receiver %d failed to collect usage informations...\n", id, buf);
	}
#endif
	/********************************************************************************************************************************/
	/*                                                         ERROR                                                                */
	/********************************************************************************************************************************/
zeend:;
 	//printf("Appolo Thread %u leaving...\n", MyID);
	return (void*) (uintptr_t) err;
}
//---------------------------------------------------------------

/* vim: tabstop=2 shiftwidth=2
 */
