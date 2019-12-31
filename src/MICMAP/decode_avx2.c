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
#define _GNU_SOURCE
#include "constants.h"
#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <sys/mman.h> 
#include <immintrin.h>
#include <pthread.h>
#include <semaphore.h>
#include <assert.h>
#include <string.h>

#include "Decoder.h"

//--------------------------------------------------------------- 
// EXTERNS 
//--------------------------------------------------------------- 
int compare_res_element(const DecoderElement_t * const restrict a, const DecoderElement_t * const restrict b);
extern int verbose;

//--------------------------------------------------------------- 
// LOCALS 
//--------------------------------------------------------------- 

//--------------------------------------------------------------- 
// LOCAL FCTS 
//--------------------------------------------------------------- 

//--------------------------------------------------------------- 
// FUNCTIONS 
//---------------------------------------------------------------
void* DecodeTags_avx2(DecoderPool_t * const restrict thpool)
{
	unsigned int mic512dump[16] __attribute__((aligned(64)));
	unsigned int packed16nt[kMAX16mersToTest] __attribute__((aligned(64)));
	unsigned int tagbin[kTagSize];
	int Aloc[kTagSize];
	int Tloc[kTagSize];
	int nt;
	unsigned int packed16ntCnt;
	unsigned int *packed16ntPtr;
	unsigned int packed4nt;
	unsigned int packed4nt_startpos[kMAX16mersToTest];
	__m256i tag16mers;
	__m256i expected4mers;
	__m256i expected4mers2;
	__m256i zero;
	__m256i candidate;
	
	const VALIDDATA * const restrict ValidDataTable = thpool->ValidDataTable;
	const unsigned int * restrict chrpos_tbl[kFollowNtCount];
	const unsigned char * restrict strand_tbl = thpool->strand_tbl;
	for (int i=0; i<kFollowNtCount; i++) chrpos_tbl[i] = thpool->chrpos_tbl[i];
	
	/* Mark thread as alive (initialized) */
	atomic_fetch_add(&(thpool->num_threads_alive), 1);
	
	while(thpool->threads_keepalive) {
		bsem_wait(thpool->PoolJob_q.has_items);
		
		if (thpool->threads_keepalive){
			/* Annonce thread working */
			atomic_fetch_add(&(thpool->num_threads_working), 1);
			
			/* Read job from queue and execute it */
			pthread_mutex_lock(&thpool->PoolJob_q.rwmutex);
			DecoderJobElement_t * const restrict task = (DecoderJobElement_t *) jobqueue_pull(&thpool->PoolJob_q);
			pthread_mutex_unlock(&thpool->PoolJob_q.rwmutex);
			
			/***************************************************************/
			/*                     START THE TASK                          */
			/***************************************************************/
			if (task) {
				const unsigned int nToBeProcessed = task->nToBeProcessed;

				/* ------------------- DECODE THAT BATCH  -------------------- */
				if (nToBeProcessed) {
					const ReaderTag_t  * const restrict Tags = task->Reader.Tags;
					DecoderElement_t (* const restrict rslt)[kDesiredHitCnt] = task->Decoder.Anchors;
					//const int loopMax = (nToBeProcessed < kDecoderTagBatch) ? nToBeProcessed : kDecoderTagBatch;
					for (unsigned int loop = 0; loop < nToBeProcessed; loop++)
					{
						unsigned int taghits = 0;

						// encode the tag and determine all potential starting points for our search
						const unsigned char * const restrict tag = Tags[loop].Tag;

						//	taglen = (unsigned int)strlen(tag);
						const unsigned short taglen = Tags[loop].AlignLen;
						//	int taglen2 = (unsigned int)tag[kTagSize-1];
						//	if (taglen2 != taglen) { printf("%d != %d\n",taglen,taglen2); }
						// we can have very small tags which we have no hope to find, or we can be
						// in the pad-filled part at the end of a buffer, where all tags printed are
						// NNNNNNNNNNNNNNNNNNNN
						if (taglen <= (kbaitLen-1)) {
							for (unsigned int i =0; i < kDesiredHitCnt; i++) {
								rslt[loop][i].which16mer = 0;
								rslt[loop][i].info = 0;
							}
							continue;
						}

						unsigned int ii;
						unsigned int pos = 0U, binpos = 0U, Acnt = 0U, Tcnt = 0U;

						do
						{
							switch(tag[pos])
							{
								case 'A':
									tagbin[binpos] = 0;
									Aloc[Acnt++] = pos;
									break;
								case 'C':
									tagbin[binpos] = 1;
									break;
								case 'G':
									tagbin[binpos] = 2;
									break;
								case 'T':
									tagbin[binpos] = 3;
									break;
								default:
									tagbin[binpos] = 0;
							}
							pos++;
							binpos++;
						} 	while (pos < (kbaitLen-1)) ;
						do
						{
							switch(tag[pos])
							{
								case 'A':
									tagbin[binpos] = 0;
									Aloc[Acnt++] = pos;
									break;
								case 'C':
									tagbin[binpos] = 1;
									break;
								case 'G':
									tagbin[binpos] = 2;
									break;
								case 'T':
									tagbin[binpos] = 3;
									Tloc[Tcnt++] = pos;
									break;
								default:
									tagbin[binpos] = 0;
							}
							pos++;
							binpos++;
						} 	while (pos < (taglen-(kbaitLen-1))) ;
						do
						{
							switch(tag[pos])
							{
								case 'A':
									tagbin[binpos] = 0;
									break;
								case 'C':
									tagbin[binpos] = 1;
									break;
								case 'G':
									tagbin[binpos] = 2;
									break;
								case 'T':
									tagbin[binpos] = 3;
									Tloc[Tcnt++] = pos;
									break;
								default:
									tagbin[binpos] = 0;
							}
							pos++;
							binpos++;
						} 	while (pos < taglen) ;

						// test forward strand candidates
						// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

						int starting16mer = 0;
			more16mer:;
						packed16ntPtr = &packed16nt[0];
						packed16ntCnt = 0;
						for (unsigned int i = starting16mer; i < Acnt; i++)
						{
							if (packed16ntCnt == kMAX16mersToTest) break;
							*packed16ntPtr = 0;
							nt = Aloc[i];
							if ((nt+16) >= taglen) break;
							do {
								*packed16ntPtr <<= 2;
								*packed16ntPtr |= tagbin[nt];
							} while (++nt < Aloc[i]+16);
							
							packed4nt_startpos[packed16ntCnt] = nt;
							packed16ntPtr++;
							packed16ntCnt++;
						}

						/* fill up */
						while(packed16ntCnt < kMAX16mersToTest)
						{
							*packed16ntPtr = 0;
							packed16ntPtr++;
							packed16ntCnt++;
						}

						tag16mers = _mm256_load_si256((__m256i const *)packed16nt);
						expected4mers = _mm256_i32gather_epi32((void*) ValidDataTable->following_nt_tbl,tag16mers, 4);
						zero = _mm256_xor_si256(zero,zero);
						candidate = _mm256_cmpeq_epi32(expected4mers, zero);	// find out which packed16nt has candidate unique match
						unsigned int first8;
						first8 = _mm256_movemask_epi8(candidate); // which of the 8 packed16nt matches.
						//printf("first8=%08x\n",first8);


						tag16mers = _mm256_load_si256((__m256i const *)&packed16nt[8]);

						expected4mers2 = _mm256_i32gather_epi32((void*)ValidDataTable->following_nt_tbl,tag16mers, 4);
						zero = _mm256_xor_si256(zero,zero);
						candidate = _mm256_cmpeq_epi32(expected4mers2, zero);	// find out which packed16nt has candidate unique match
						unsigned int next8;
						next8 = _mm256_movemask_epi8(candidate); // which of the 8 packed16nt matches.
						//printf("next8=%08x\n",next8);


						unsigned long j = ((unsigned long)next8 << 32) | first8;
						if (j != 0xFFFFFFFFFFFFFFFF) // verify if following N nt are matching. we should be able to verctorize and process more than one at the same time ?
						{
							_mm256_store_si256((__m256i *)mic512dump, expected4mers);
							_mm256_store_si256((__m256i *)&mic512dump[8], expected4mers2);


							pos = 0;
							while((j & 1) != 0) { j >>= 4; pos++; if (pos == kMAX16mersToTest) goto noluck; };

			tryAgain:;
							unsigned int rev = ((packed4nt_startpos[pos] & 0x80000000) != 0);
							if (rev == 0) {
								unsigned int match;
								unsigned int curMask;
								char s[256];
								ii = packed4nt_startpos[pos] &= ~0x80000000;

								packed4nt = tagbin[ii];
								packed4nt <<= 2;
								packed4nt |= tagbin[ii+1];
								packed4nt <<= 2;
								packed4nt |= tagbin[ii];
								packed4nt <<= 2;
								packed4nt |= tagbin[ii+1];

								packed4nt = (packed4nt << 8) | packed4nt;       // repeat the 2nt everywhere...
								packed4nt = (packed4nt << 16) | packed4nt;
								packed4nt ^= mic512dump[pos];

								// encode results
								ii -= 16;

								for (match = 0, curMask = (1 << (kFollowNtLen*2)) - 1; match < kFollowNtCount && taghits < kDesiredHitCnt; match += 1, curMask <<= (kFollowNtLen*2))
								{
									if ((packed4nt & curMask) == 0)
									{
										// we have a match
										// check if valid
										unsigned long tmp = match*k15NT_ENTRIES_COUNT + packed16nt[pos]; // position as if byte array.
										unsigned char getbit = (tmp & 0x07);                 // get the position of the bit we need to get.
										if (((ValidDataTable->valid_tbl[tmp >> 3]) >> getbit) & 0x01)       // positions of byte array divided by 8, as we pack 8 bits in each char.
										{
											rslt[loop][taghits].which16mer = packed16nt[pos];
											rslt[loop][taghits].info  = (match << 16);			// whichNnt
											rslt[loop][taghits].info |= (ii << 8);				// tagpos
											rslt[loop][taghits].info |= (pos << 2);				// which of the kMAX16mersToTest was successful
											//rslt[loop].info[taghits] |= (rev << 1);			// tested tag was reversed ?   always false in this context.
	#ifdef DBGHITS
											printf("FWD hit%2d loop=%2d %.100s matchpos=%3d %.*s which16mer=%d : %d match %d\n",taghits,loop,tag,ii,kbaitLen,&tag[ii],packed16nt[pos],pos,match);
	#endif
											taghits++;
										}
										else
											break; // the rest will be invalid as well
									}
								}
							}

							// Awful: need to refactor code, although it is still easy to follow the logic.
							// for now, if the tested candidate 16mer does not exapnd to a Nmer hit, then we go check the next candidate...
							if (taghits < kDesiredHitCnt) {
								do {j >>= 4; pos++; if (pos == kMAX16mersToTest) goto noluck; }  while((j & 1) != 0);
								//printf("j=%08lx pos=%d  goto tryAgain\n",j,pos); fflush(stdout);
								goto tryAgain;
							}
						} // candidate != 0
			noluck:;

						// if still no hit, then encode the next 8 candidates, and try again.
						starting16mer += 16;
						//printf("starting16mer=%d  Acnt=%d  taghits=%d kDesiredHitCnt=%d\n",starting16mer,Acnt,taghits,kDesiredHitCnt); fflush(stdout);
						if ((taghits < kDesiredHitCnt)  && (Acnt > starting16mer)) goto more16mer;
						
						// end of forward
						// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

						//printf("NEXT\n"); fflush(stdout);

						// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
						// try with reverse.  empirical effective gain is only about 4.6% !    :-(
						if (taghits < kDesiredHitCnt)
						{

							starting16mer = 0;
			moreRev16mer:;
							//printf("reverse\n"); fflush(stdout);

							/* reverse */
							packed16ntPtr = &packed16nt[0];
							packed16ntCnt = 0;
							for (unsigned int i = starting16mer; i < Tcnt; i++)
							{
								if (packed16ntCnt == kMAX16mersToTest) break;
								*packed16ntPtr = 0;
								nt = Tloc[i];
								if (nt > kbaitLen) // make sure we still have N nt following...
								{
									do
									{
										*packed16ntPtr <<= 2;
										*packed16ntPtr |= tagbin[nt];
									} while (--nt >= (Tloc[i] - 15));
									packed4nt_startpos[packed16ntCnt] = nt;
									packed4nt_startpos[packed16ntCnt] |= 0x80000000;  // indicate it is rev.
									*packed16ntPtr ^= 0xFFFFFFFF;
									packed16ntPtr++;
									packed16ntCnt++;
								}
							}
							//printf("fill up\n"); fflush(stdout);

							/* fill up */
							while(packed16ntCnt < kMAX16mersToTest)
							{
								*packed16ntPtr = 0;
								packed16ntPtr++;
								packed16ntCnt++;
							}

							//printf("filled\n"); fflush(stdout);
							// now do the real job, checking if we have hits among any of the 16 x 16mers.

							tag16mers = _mm256_load_si256((__m256i const *)packed16nt);
							expected4mers = _mm256_i32gather_epi32((void*)ValidDataTable->following_nt_tbl,tag16mers, 4);
							zero = _mm256_xor_si256(zero,zero);
							candidate = _mm256_cmpeq_epi32(expected4mers, zero);	// find out which packed16nt has candidate unique match
							unsigned int first8;
							first8 = _mm256_movemask_epi8(candidate); // which of the 8 packed16nt matches.
							//printf("first8=%08x\n",first8); fflush(stdout);


							tag16mers = _mm256_load_si256((__m256i const *)&packed16nt[8]);
							expected4mers2 = _mm256_i32gather_epi32((void*)ValidDataTable->following_nt_tbl,tag16mers, 4);
							zero = _mm256_xor_si256(zero,zero);
							candidate = _mm256_cmpeq_epi32(expected4mers2, zero);	// find out which packed16nt has candidate unique match
							unsigned int next8;
							next8 = _mm256_movemask_epi8(candidate); // which of the 8 packed16nt matches.
							//printf("next8=%08x\n",next8); fflush(stdout);

							unsigned long j = ((unsigned long)next8 << 32) | first8;

							if (j != 0xFFFFFFFFFFFFFFFF) // verify if following N nt are matching. we should be able to verctorize and process more than one at the same time ?
							{
								_mm256_store_si256((__m256i *)mic512dump, expected4mers);
								_mm256_store_si256((__m256i *)&mic512dump[8], expected4mers2);

								pos = 0;
								while((j & 1) != 0) { j >>= 4; pos++; if (pos == kMAX16mersToTest) goto noluckRev; };

			tryAgainRev:;
								unsigned int rev = ((packed4nt_startpos[pos] & 0x80000000) != 0);

								assert(rev != 0);
								
								// reverse tag.
								{
									unsigned int match;
									unsigned int curMask;
									char s[256];
									ii = packed4nt_startpos[pos] & ~0x80000000;

									packed4nt = tagbin[ii];
									packed4nt <<= 2;
									packed4nt |= tagbin[ii-1];
									packed4nt <<= 2;
									packed4nt |= tagbin[ii];
									packed4nt <<= 2;
									packed4nt |= tagbin[ii-1];

									packed4nt = (packed4nt << 8) | packed4nt;       // repeat the 2nt everywhere...
									packed4nt = (packed4nt << 16) | packed4nt;
									packed4nt ^= 0xFFFFFFFF;						// reverse complement.
									packed4nt ^= mic512dump[pos];

									// encode results

									ii -= (kFollowNtLen-1); // in fwd we do a -16 because we test the 16mers and then validate with N extra nt.
													// here we reverse comp the tag not sure...
									for (match = 0, curMask = (1 << (kFollowNtLen*2)) - 1; match < kFollowNtCount && taghits < kDesiredHitCnt; match += 1, curMask <<= (kFollowNtLen*2))
									{
										if ((packed4nt & curMask) == 0)
										{
											// we have a match
											// check if valid
											unsigned long tmp = match*k15NT_ENTRIES_COUNT + packed16nt[pos]; // position as if byte array.
											unsigned char getbit = (tmp & 0x07);            // get the position of the bit we need to get.
											if (((ValidDataTable->valid_tbl[tmp >> 3]) >> getbit) & 0x01)   // positions of byte array divided by 8, as we pack 8 bits in each char.
											{
												rslt[loop][taghits].which16mer = packed16nt[pos];
												rslt[loop][taghits].info  = (match << 16);    // whichNnt
												rslt[loop][taghits].info |= (ii << 8);        // tagpos
												rslt[loop][taghits].info |= (pos << 2);       // which of the kMAX16mersToTest was successful
												rslt[loop][taghits].info |= (rev << 1);       // tested tag was reversed
	#ifdef DBGHITS
												printf("REV hit%2d loop=%2d %.100s matchpos=%3d %.*s which16mer=%d : %d match %d\n",taghits,loop,tag,ii,kbaitLen,&tag[ii],packed16nt[pos],pos,match);
	#endif
												taghits++;
											}
											else
												break; // the rest will be invalid as well
										}
									}
								}

								// Awful: need to refactor code, although it is still easy to follow the logic.
								// for now, if the tested candidate 16mer does not exapnd to a Nmer hit, then we go check the next candidate...
								if (taghits < kDesiredHitCnt)
								{
									do {j >>= 4; pos++; if (pos == kMAX16mersToTest) goto noluckRev; }  while((j & 1) != 0);
									//printf("j=%08lx pos=%d  goto tryAgainRev\n",j,pos); fflush(stdout);
									goto tryAgainRev;
								}
							} // candidate != 0
			noluckRev:;

							// if still no hit, then encode the next 16 candidates, and try again.
							starting16mer += 16;
							//printf("starting16mer=%d  Tcnt=%d  taghits=%d kDesiredHitCnt=%d\n",starting16mer,Tcnt,taghits,kDesiredHitCnt); fflush(stdout);
							if ((taghits < kDesiredHitCnt)  && (Tcnt > starting16mer)) goto moreRev16mer;

						}	// if taghits < kDesiredHitCnt

						// end of reverse
						// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	#ifdef RECORD_STATS
						taghits_stats[(((DECODETAGS*)data)->buf)][chunk][taghits]++;
	#endif
						// make sure to flag tag as not matching
						if (taghits < kDesiredHitCnt)
						{
							for (unsigned int i = taghits; i < kDesiredHitCnt; i++)
							{
								rslt[loop][i].which16mer = 0;
								rslt[loop][i].info = 0;
							}
						}

					} // loop

					for (unsigned int loop=0; loop < nToBeProcessed;  loop++)
					{
						unsigned int whichNnt;
						const unsigned int *tbl;
						int hitA;
						int tryHits = 0;
						const int taglen = Tags[loop].AlignLen;
						for (hitA = 0; hitA < kDesiredHitCnt; hitA++)
						{
							if (rslt[loop][hitA].which16mer == 0)
								break; // we reached the end
							whichNnt = (rslt[loop][hitA].info >> 16);
							/* get the reverse flag */
							unsigned long tmp = whichNnt*k15NT_ENTRIES_COUNT + rslt[loop][hitA].which16mer; // position as if byte array.
							unsigned char getbit = (tmp & 0x07);                             // get the position of the bit we need to get.
							if (((strand_tbl[tmp >> 3]) >> getbit) & 0x01)                   // positions of byte array divided by 8, as we pack 8 bits in each char.
								rslt[loop][hitA].info |= 1;                                       // original Nmer is a reverse complement of the genome.
							tbl = chrpos_tbl[whichNnt];
							// check if the chr:starting position has already been seen
							DecoderElement_t RE;
							RE.which16mer = tbl[rslt[loop][hitA].which16mer];
							if (RE.which16mer != 0)
							{
								int j;
								int cmp = 1;
								RE.info = rslt[loop][hitA].info;
								// Do the pos ajustment here, because we want to remove duplicates at this point
								if (RE.info & 1)
								{
									if ((RE.info >> 1) & 0x01)  // revNmer, but the tag was also reversed !!
									{
										if (WHICH16MER_pos(RE.which16mer) > ((RE.info >> 8) & 0xFF))
											RE.which16mer -= ((RE.info >> 8) & 0xFF);
										else
											RE.which16mer = WHICH16MER_chr_part(RE.which16mer) | 1;
									}
									else
									{
										if (WHICH16MER_pos(RE.which16mer) > taglen - kbaitLen - ((RE.info >> 8) & 0xFF))
											RE.which16mer += kbaitLen - taglen + ((RE.info >> 8) & 0xFF);
										else
											RE.which16mer = WHICH16MER_chr_part(RE.which16mer) | 1;
									}
								}
								else
								{
									if ((RE.info >> 1) & 0x01)  // fwdNmer, but the tag was reversed !!
									{
										if (WHICH16MER_pos(RE.which16mer) > taglen - kbaitLen - ((RE.info >> 8) & 0xFF))
											RE.which16mer += kbaitLen - taglen + ((RE.info >> 8) & 0xFF);
										else
											RE.which16mer = WHICH16MER_chr_part(RE.which16mer) | 1;
									}
									else
									{
										if (WHICH16MER_pos(RE.which16mer) > ((RE.info >> 8) & 0xFF))
											RE.which16mer -= ((RE.info >> 8) & 0xFF);
										else
											RE.which16mer = WHICH16MER_chr_part(RE.which16mer) | 1;
									}
								}
								//fprintf(stderr,"%d %d : %u %u %u %u %u %u\n",i,hitA,RE.which16mer & 0x0fffffff,whichNnt,(RE.info >> 8) & 0xff,(RE.info >> 2) & 0x3f,(RE.info >> 1) & 1,RE.info & 1);
								for (j = 0; j < tryHits; j++)
								{
									cmp = compare_res_element(&RE,rslt[loop] + j);
									if (cmp <= 0)
										break;
								}
								// at this point RE is less than or equal to rslt[loop][j], or we have exhausted j
								// if equal we are done - no need to add this one
								if (cmp < 0)
								{
									// we found the spot - need to insert the new element
									assert(j < tryHits);
									unsigned int k = tryHits;
									while (k > j)
									{
										rslt[loop][k] = rslt[loop][k-1];
										k -= 1;
									}
									rslt[loop][j] = RE;
									tryHits += 1;
								}
								if (cmp > 0)
								{
									// we reached the end of the current sorted array - just append the new one
									rslt[loop][tryHits++] = RE;
								}
							}
						}//hitA
						while (tryHits < kDesiredHitCnt)
							rslt[loop][tryHits++].which16mer = 0;
					}

					/* Update the counter for last processing thread to be aware we are done here */
					atomic_fetch_add(task->count, nToBeProcessed);
					
					/* Push the DecoderJob back to empty slot */
					//printf("Task 0x%lx, count @ 0x%lx = %u, in 0x%lx, out 0x%lx\n", task, task->count, atomic_load(task->count), Tags, rslt); fflush(stdout);
					pthread_mutex_lock(&(thpool->PoolMemory_q.rwmutex));
					jobqueue_push(&(thpool->PoolMemory_q), (job_t*) task);
					pthread_mutex_unlock(&(thpool->PoolMemory_q.rwmutex));
				}
				/* -------------- WE ARE LAST, PERFORM POSTPROCESSING -------- */
				else {
					//printf("Task 0x%lx terminating reader job @ 0x%lx\n", task, task->Reader.ReaderJob); fflush(stdout);
					const ReaderJob_t * const restrict ReaderJob = task->Reader.ReaderJob;

					/* Wait until Decoders are done writing */
					register const unsigned int TagsInReaderJob = ReaderJob->TagCnt;
					while (atomic_load(task->count) < TagsInReaderJob) ;
					
					/* Push the DecoderJob further */
					if (verbose & 0x1)
					{
						printf("Task 0x%lx adding jobs with %u tags to Aligner Queue\n", (uintptr_t) task, atomic_load(task->count)); fflush(stdout);
					}
					pthread_mutex_lock(&(thpool->ToBeAligned_q->rwmutex));
					jobqueue_push(thpool->ToBeAligned_q, (job_t*) task);
					pthread_mutex_unlock(&(thpool->ToBeAligned_q->rwmutex));
				}
			}
			
			/* Annonce thread not working animore*/
			atomic_fetch_sub(&(thpool->num_threads_working), 1);
		}
	}
	
	/* Annonce thread about to end */
	atomic_fetch_sub(&(thpool->num_threads_alive), 1);
	//printf("%s: closing...\n", __FUNCTION__);	
	return (void*) 0;
} 
//---------------------------------------------------------------

/* vim: tabstop=2 shiftwidth=2
 */
