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
extern int verbose;

//========================================================================================================
// shiftTagAgainstGenome was modified to directly return the score and 
// according to the scoring system of pftools aligner for further comparison.
// That is: # matches * _M + (taglen - # mathches) * _m 
static unsigned int shiftTagAgainstGenome(const unsigned char * const restrict tag, const int taglen,
                                          const unsigned char * const restrict genome,
                                          InternalAlignment * const restrict tagchoice)
{
	__assume_aligned(tag, 16);
	int tp;
	int from;
#ifdef DEBUG_GENERIC_shiftTagAgainstGenome
	char DBGg[128];
	char DBGt[128];
#endif

#ifdef DEBUG_GENERIC_shiftTagAgainstGenome
	printf("TESTING %s\nchr=%d pos=%9d matchpos=%9d\n\n",tag,tagchoice->chr,tagchoice->chr_pos,
	       tagchoice->best.matchpos);
	fflush(stdout);
#endif

	assert(((uintptr_t) tag & 0xF) == 0); // tag should be aligned on 64B boundary
	const unsigned char * restrict ptrTag = tag;  __assume_aligned(ptrTag, 16);

#ifdef DEBUG_GENERIC_shiftTagAgainstGenome
	strncpy(DBGg,genome,127);
	strncpy(DBGt,ptrTag,taglen);
	DBGg[127] = 0;
	DBGg[taglen] = 0;
	printf("GEN=%.*s\n", taglen, DBGg);
	printf("TAG=%.*s\n\n",taglen, DBGt);
	fflush(stdout);
#endif

	// test perfect match for quick abort.
	for (tp = 0; tp < taglen; tp++)
	{
		if (!((ptrTag[tp] == genome[tp]) || (ptrTag[tp] == 'N')))
			break;
	}
	if (tp == taglen)
	{
#ifdef DEBUG_GENERIC_shiftTagAgainstGenome
		printf("direct hit\n");
#endif
		tagchoice->score         = taglen*_M;
		tagchoice->MismatchCount = 0U;
		// WARNING: WE COMMENT THAT OUT !!!!!!!!!!!!!!!!!!!!!!!!!!
		//tagchoice->chr_pos           = tagchoice->matchpos;
		// WARNING: WE COMMENT THAT OUT !!!!!!!!!!!!!!!!!!!!!!!!!!
		return 0;
	}
	
#ifdef DEBUG_GENERIC_shiftTagAgainstGenome
	strncpy(DBGg,genome,127);
	strncpy(DBGt,ptrTagPtr,taglen);
	DBGg[127] = 0;
	DBGg[taglen] = 0;
	printf("GEN=%.*s\n", taglen, DBGg);
	printf("TAG=%.*s\n\n", taglen, DBGt);
	fflush(stdout);
#endif

	
	// start checking a bit upstream in case a gap would be welcome.
	if (tagchoice->chr_pos >= 16)
		from = -16;
	else
		from = -tagchoice->chr_pos;
		
	/* check best acceptable local ungapped match nearby  */

	unsigned int smallestMismatchcnt = (unsigned int) taglen;
	int bestshift = 0;

	for (int shift = from; shift < 16; shift++)
	{
		unsigned int mismatchcnt = 0;
		for (tp = 0; tp < taglen; tp++)
		{
			if ((ptrTag[tp] != 'N') && (ptrTag[tp] != genome[shift+tp])) mismatchcnt++;
		}

		if (mismatchcnt < smallestMismatchcnt)
		{
			// store current best result
			smallestMismatchcnt = mismatchcnt;
			bestshift = shift;
		}
		
#ifdef DEBUG_GENERIC_shiftTagAgainstGenome
		strncpy(DBGg,&genome[shift],127);
		strncpy(DBGt,ptrTagPtr,taglen);
		DBGg[127] = 0;
		DBGg[taglen] = 0;
		printf("GEN=%.*s\n", taglen, DBGg);
		printf("TAG=%.*s shift=%d match=%d\n\n", taglen, DBGt,shift,match);
		fflush(stdout);
#endif
	} // shift

	tagchoice->MismatchCount = smallestMismatchcnt;
	tagchoice->chr_pos      += bestshift;

	return smallestMismatchcnt;
} // shiftTagAgainstGenome

//========================================================================================================
static inline void __attribute__((always_inline))
EncodeMM(const unsigned char * const restrict tag, InternalAlignment * const restrict alnrslt,
         const unsigned char * const restrict genome, const unsigned int taglen)
{
	/* NOTE: we go in reverse order as the alignement so that we will not need a copy on the
	 *       host anymore.
	 */	
#ifdef DEBUG_TAG_MAPPING_PROCESS
	printf("ENCODING mismatches : %d:%u\n",alnrslt->chr,genomepos);
#endif
	/* ----- check tag against genome ----- */

	const unsigned char * const restrict tp  = tag;
	unsigned char * restrict MMString = (unsigned char*) &(alnrslt->StateMemory[0]);
	// String Limit has to be defined based upon the result structure i.e. EXTRESULTS
	unsigned char * const restrict StringLimit = (unsigned char*) &(alnrslt->StateMemory[(kMAXmismatchSize-1)]);
	unsigned int mmpos = 0;
	unsigned char IsAllN = 1;
	unsigned char ContainsN = 0;
	int maxScore = (_M + _MM)*taglen;
	
	while (mmpos < taglen) {
		register const unsigned char val = tp[mmpos];
		if (val != genome[mmpos]) {
			MMString[0] = (unsigned char) mmpos;
			MMString[1] = val;
			MMString += 2;
			if ( MMString >= StringLimit) return;
			if (val != 'N') {
				IsAllN = 0;
				maxScore += _m - _M;
			}
			else {
				ContainsN = 1;
			}
		}
		mmpos++;
	}

	*MMString = kMMTERMINATOR;
	
	// Set All is N and success encoding
	unsigned int uitmp;
	
	if ((uintptr_t) MMString == (uintptr_t) &(alnrslt->StateMemory[0])) {
		uitmp = kPerfectMatch;
	} 
	else {
		uitmp = (IsAllN) ? kAllIsN : 0;
	}
	if (ContainsN) uitmp |= kContainsN;
	alnrslt->score = maxScore;
	alnrslt->DecisionTreeHistory |= uitmp | kMMEncoded | kSuccessfullyEncoded;
#ifdef DEBUG_TAG_MAPPING_PROCESS
	printf("TAG: %.*s\n   : %.*s\nGEN: %.*s\n", taglen, tag, taglen, tp, taglen, genome);
#endif
} // encodeMM

//========================================================================================================
// EncodeInternalAlignment transforms a state string into the mismatch and cigar strings while 
// adding kMapped, kSuccessfullyEncoded, kAllIsN and kContainsN flags.
// Softclip is handle when kSoftClipRescued is activated but requires both
//  - SoftClipMismatchToRemove
//
// Cigar is set to its terminator when everything is a match or mismatch
static inline unsigned int __attribute__((always_inline))
EncodeInternalAlignment(const unsigned char * const restrict Tag, InternalAlignment * const restrict Trial,
                        AlignmentResults_t * const restrict extrslt, const ReaderTag_t * const restrict TagData )
{
	char * restrict ptrState = Trial->States;
	unsigned char * restrict Mismatch = extrslt->Mismatch;
	const unsigned char * const restrict MismatchLimit = Mismatch + (2*kMaxEncodedMismatch);
	unsigned char * restrict Cigar = extrslt->Cigar;
	const unsigned char * restrict const CigarLimit = Cigar + (2*kMaxEncodedCigar);
	unsigned char IsAllN = 1;
	unsigned char ContainsN = 0;
	unsigned int id = 0U;
	assert(TagData->TagLen >= TagData->AlignLen);
	unsigned int SoftClipSize = TagData->TagLen - TagData->AlignLen;
	int chr_pos_shift = 0;
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// WITH SOFTCLIP
	int RemovedMismatch = 0;
	int RemovedCigar = 0;
	int RemoveOther = 0;
	int MismatchToRemove = 0;
	int CigarToRemove = 0;

	/* Do it right
	 *                                   + LastDumpCigarLimit ptr = CigarLimit - 4
	 *                                   |             + this is CigarLimit ptr, index 2*kMaxEncodedCigar
	 *                                   |             |
	 * ... | 35 M |  3 I | 10 M |  1 D | 23 M | 45 S | TERMINATOR 
	 * 
	 */
	if ( Trial->DecisionTreeHistory & kSoftClipRescued) {
		const unsigned int SoftClipMismatchToRemove = (unsigned int) Trial->SoftClipMismatchToRemove;
		assert(SoftClipMismatchToRemove > 0);
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		// TAG is reversed
		if (Trial->revNmer) {
			unsigned int pos = 0U;
			//chr_pos_shift -= SoftClipSize;	
			char  previousState = *ptrState;
			if (previousState != 'I'&& previousState != 'D') previousState = 'M';
			
			while (RemovedMismatch < SoftClipMismatchToRemove) {
				switch(*ptrState) {
					case 'm':
						RemovedMismatch++;
						if (Tag[pos] == 'N')
							ContainsN = 1;
						else
							IsAllN = 0;
						++pos;
						break;
					case 'I':
						chr_pos_shift--;
					case 'n':
					case 'N':
						RemoveOther++;
						if (Tag[pos] == 'N')
							ContainsN = 1;
						else
							IsAllN = 0;
					case 'M':
						++pos;
						break;
					case 'D':
						chr_pos_shift++;
						break;
					default:
						fprintf(stderr, "Unknown state '%i' found @ %i softclip reverse (%u/%u) score %i\n%s\n",
										(int) *ptrState, pos, RemovedMismatch, SoftClipMismatchToRemove, Trial->score,
										Trial->States);
						fflush(stderr);
						goto bail;
						break;
				}
				char currentState = *ptrState;
				if (currentState != 'I'&& currentState != 'D') currentState = 'M';
				if (currentState != previousState) {
					RemovedCigar++;
					previousState = currentState;
				}
				ptrState++;
			}
			
			{
				MismatchToRemove = Trial->TotalRequiredMMCount - RemovedMismatch - kMaxEncodedMismatch;
				while (RemoveOther < MismatchToRemove) {
					switch(*ptrState) {
						case 'm':
							if (Tag[pos] == 'N')
								ContainsN = 1;
							else
								IsAllN = 0;
							++pos;
							break;
						case 'I':
							chr_pos_shift--;
						case 'n':
						case 'N':
							RemoveOther++;
							if (Tag[pos] == 'N')
								ContainsN = 1;
							else
								IsAllN = 0;
						case 'M':
							++pos;
							break;
						case 'D':
							chr_pos_shift++;
							break;
						default:
							fprintf(stderr, "Unknown state '%i' found @ %i softclip reverse (%u/%u) score %i\n%s\n",
											(int) *ptrState, pos, RemoveOther, MismatchToRemove, Trial->score, Trial->States);
							fflush(stderr);
							goto bail;
							break;
					}
					char currentState = *ptrState;
					if (currentState != 'I'&& currentState != 'D') currentState = 'M';
					if (currentState != previousState) {
						RemovedCigar++;
						previousState = currentState;
					}
					ptrState++;
				}
			}
			
			{
				CigarToRemove = Trial->TotalRequiredCigCount - kMaxEncodedCigar;
				while (RemovedCigar <= CigarToRemove) {
					switch(*ptrState) {
						case 'm':
							if (Tag[pos] == 'N')
								ContainsN = 1;
							else
								IsAllN = 0;
							++pos;
							break;
						case 'I':
							chr_pos_shift--;
						case 'n':
						case 'N':
							if (Tag[pos] == 'N')
								ContainsN = 1;
							else
								IsAllN = 0;
						case 'M':
							++pos;
							break;
						case 'D':
							chr_pos_shift++;
							break;
						default:
							fprintf(stderr, "Unknown state '%i' found @ %i softclip reverse (%u/%u) score %i\n%s\n",
											(int) *ptrState, pos, RemovedCigar, CigarToRemove, Trial->score, Trial->States);
							fflush(stderr);
							goto bail;
							break;
					}
					char currentState = *ptrState;
					if (currentState != 'I' && currentState != 'D') currentState = 'M';
					if (currentState != previousState) {
						RemovedCigar++;
						previousState = currentState;
					}
					ptrState++;
				}
			}

			// Add softclip region 
#if kMAXcigarSize <= 2
#error "kMAXcigarSize should be at least 3"
#endif
			SoftClipSize += pos;
			assert(SoftClipSize < 255);
			Cigar[0] = (unsigned char) SoftClipSize;
			Cigar[1] = 'S';
			Cigar += 2;
			
			previousState = *ptrState;
			if (previousState != 'I'&& previousState != 'D') previousState = 'M';
			unsigned int StateCount = 0U;
			
			while (*ptrState != '\0') {
				char currentState;
				switch(*ptrState) {
					case 'I':
						if (Mismatch >= MismatchLimit) goto bail;
						Mismatch[0] = (unsigned char) pos;
						Mismatch[1] = Tag[pos];
						Mismatch += 2;
						if (Tag[pos] == 'N')
							ContainsN = 1;
						else
							IsAllN = 0;
						currentState = 'I';
						++pos;
						break;
					case 'm':
					case 'n':
					case 'N':
						if (Mismatch >= MismatchLimit) goto bail;
						Mismatch[0] = (unsigned char) pos;
						Mismatch[1] = Tag[pos];
						Mismatch += 2;
						if (Tag[pos] == 'N')
							ContainsN = 1;
						else
							IsAllN = 0;
					case 'M':
						currentState = 'M';
						++pos;
						break;
					case 'D':
						currentState = 'D';
						break;
					default:
						fprintf(stderr, "Unknown state '%i' found @ %i softclip reverse score %i\n%s\nCIG: ", 
										(int) *ptrState, pos, Trial->score, Trial->States);
						unsigned char * ptr = extrslt->Cigar;
						while (ptr < Cigar) { fprintf(stderr, "%i%c ", (int) ptr[0], ptr[1]); ptr += 2; }
						fprintf(stderr, "\nMM : ");
						ptr = extrslt->Mismatch;
						while (ptr < Mismatch) { fprintf(stderr, "%i%c ", (int) ptr[0], ptr[1]); ptr += 2; }
						fprintf(stderr, "\n");
						fflush(stderr);
						goto bail;
						break;
				}
				
				if (previousState != currentState) {
					if (Cigar >= CigarLimit) goto bail;
					assert(StateCount < 256);
					Cigar[0] = ((unsigned char) StateCount);
					Cigar[1] = previousState;
					Cigar += 2;
					previousState = currentState;
					StateCount = 1U;
				}
				else {
					StateCount++;
				}
				ptrState++;
			}
			// terminate 
			*Mismatch = kMMTERMINATOR;

			if (StateCount) {
				if (Cigar >= CigarLimit) goto bail;
				assert(StateCount < 256);
				Cigar[0] = ((unsigned char) StateCount);
				Cigar[1] = previousState;
				Cigar += 2;
			}
			
			if (Cigar > CigarLimit) goto bail;
			*Cigar = kCIGARTERMINATOR;
			
		}
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		// TAG is correctly oriented
		else {
			unsigned int pos = 0U;
			unsigned int KeptMismatch = 0U;
			char previousState = *ptrState;
			if (previousState != 'I'&& previousState != 'D') previousState = 'M';
			unsigned int StateCount = 0U;
			
			const unsigned char * const restrict LastDumpCigarLimit = CigarLimit - 2 - 2;
			while (KeptMismatch <= kMaxMismatch) {
				char currentState = *ptrState;
				if (currentState != 'I'&& currentState != 'D') currentState = 'M';
				if (previousState != currentState) {
					if (Cigar >= LastDumpCigarLimit) {
						if (pos >= TagData->ValidQualityLen) {
							goto CheckFurtherMM;
						}
						else {
#ifndef NDEBUG
							printf("Within loop exit at cigar location %lu, pos %u, softclip @ %u\n",
										 (uintptr_t) Cigar - (uintptr_t) extrslt->Cigar,
										 pos, TagData->ValidQualityLen);
#endif
							goto bail;;
						}
					}
					assert(StateCount < 256);
					Cigar[0] = ((unsigned char) StateCount);
					Cigar[1] = previousState;
					Cigar += 2;
					previousState = currentState;
					StateCount = 1U;
				}
				else {
					StateCount++;
				}
				
				switch(*ptrState) {
					case 'I':
						if (Mismatch >= MismatchLimit) {
							if (pos >= TagData->ValidQualityLen) {
								Trial->MismatchCount = KeptMismatch;
								/* Correct cigar since that base has been accounted */
								StateCount--;
								goto CheckFurtherMM;
							}
							else {
								goto bail;
							}
						}
						Mismatch[0] = (unsigned char) pos;
						Mismatch[1] = Tag[pos];
						Mismatch += 2;
						if (Tag[pos] == 'N')
							ContainsN = 1;
						else
							IsAllN = 0;
						pos++;
						break;
					case 'm':
						if (++KeptMismatch > kMaxMismatch) {
							/* Correct cigar if necessary */
							StateCount--;
							goto CheckFurtherMM;
						}
					case 'n':
					case 'N':
						if (Mismatch >= MismatchLimit) {
							if (pos >= TagData->ValidQualityLen) {
								/* Correct cigar since that base has been accounted */
								StateCount--;
								goto CheckFurtherMM;
							}
							else {
								goto bail;
							}
						}
						Mismatch[0] = (unsigned char) pos;
						Mismatch[1] = Tag[pos];
						Mismatch += 2;
						if (Tag[pos] == 'N')
							ContainsN = 1;
						else
							IsAllN = 0;
					case 'M':
						pos++;
						break;
					case 'D':
						break;
					default:
						{
							fprintf(stderr, "Unknown state '%i' found @ %i softclip (%u/%u) score %i\n%s\nCIG: ",
											(int) *ptrState, pos, KeptMismatch, kMaxMismatch, Trial->score, Trial->States);
							unsigned char * ptr = extrslt->Cigar;
							while (ptr < Cigar) { fprintf(stderr, "%i%c ", (int) ptr[0], ptr[1]); ptr += 2; }
							fprintf(stderr, "\nMM : ");
							ptr = extrslt->Mismatch;
							while (ptr < Mismatch) { fprintf(stderr, "%i%c ", (int) ptr[0], ptr[1]); ptr += 2; }
							fprintf(stderr, "\n");
							fflush(stderr);
							goto bail;
						}
						break;
				}
				ptrState++;
			}
			// terminate 
		CheckFurtherMM:;
			if (StateCount > 0U) {
				if (previousState != 'D') {
					if (Cigar >= CigarLimit) goto bail;
					assert(StateCount < 256);
					Cigar[0] = ((unsigned char) StateCount);
					Cigar[1] = previousState;
					Cigar += 2;
				}
			}
			*Mismatch = kMMTERMINATOR;
			
			// Add softclip region 
			if (Cigar > CigarLimit) goto bail;
			assert(TagData->AlignLen > pos);
			SoftClipSize += (TagData->AlignLen - pos);
			
			Cigar[0] = (unsigned char) (SoftClipSize); // potentially +1
			Cigar[1] = 'S';
			Cigar += 2;
			
			if (Cigar > CigarLimit) goto bail;
			*Cigar = kCIGARTERMINATOR;
		}
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// WITHOUT SOFTCLIP
	else {
		
		if ((SoftClipSize > 0) && Trial->revNmer) {
			chr_pos_shift -= SoftClipSize;
			assert(SoftClipSize < 256);
			Cigar[0] = ((unsigned char) SoftClipSize);
			Cigar[1] = 'S';
			Cigar += 2;
		}
		
		char previousState = *ptrState;
		if (previousState != 'I' && previousState != 'D') previousState = 'M';
		unsigned int StateCount = 0U;
		unsigned int pos = 0U;
		while (*ptrState != '\0') {
			char currentState;
			switch(*ptrState) {
				case 'I':
					if (Mismatch >= MismatchLimit) goto bail;
					Mismatch[0] = (unsigned char) pos;
					Mismatch[1] = Tag[pos];
					Mismatch += 2;
					if (Tag[pos] == 'N')
						ContainsN = 1;
					else
						IsAllN = 0;
					currentState = 'I';
					++pos;
					break;
				case 'N':
				case 'n':
				case 'm':
					if (Mismatch >= MismatchLimit) goto bail;
					Mismatch[0] = (unsigned char) pos;
					Mismatch[1] = Tag[pos];
					Mismatch += 2;
					if (Tag[pos] == 'N')
						ContainsN = 1;
					else
						IsAllN = 0;
				case 'M':
					currentState = 'M';
					++pos;
					break;
				case 'D':
					currentState = 'D';
					break;
				default:
					fprintf(stderr, "Unknown state '%c' found @ %i no softclip score %i\n", *ptrState, pos, Trial->score);
					unsigned char * ptr = extrslt->Cigar;
						while (ptr < Cigar) { fprintf(stderr, "%i%c ", (int) ptr[0], ptr[1]); ptr += 2; }
					fprintf(stderr, "\nMM : ");
					ptr = extrslt->Mismatch;
					while (ptr < Mismatch) { fprintf(stderr, "%i%c ", (int) ptr[0], ptr[1]); ptr += 2; }
					fprintf(stderr, "\n");
					fflush(stderr);
					goto bail;
					break;
			}
			

			if (previousState != currentState) {
				if (Cigar >= CigarLimit) goto bail;
				assert(StateCount < 256);
				Cigar[0] = ((unsigned char) StateCount);
				Cigar[1] = previousState;
				Cigar += 2;
				previousState = currentState;
				StateCount = 1U;
			}
			else {
				StateCount++;
			}
			ptrState++;
		}
		
		*Mismatch = kMMTERMINATOR;

		if (StateCount) {
			if (Cigar >= CigarLimit) goto bail;
			assert(StateCount < 256);
			Cigar[0] = ((unsigned char) StateCount);
			Cigar[1] = previousState;
			Cigar += 2;
		}
		
		if ((SoftClipSize > 0) && !(Trial->revNmer)) {
			assert(SoftClipSize < 256);
			Cigar[0] = ((unsigned char) SoftClipSize);
			Cigar[1] = 'S';
			Cigar += 2;
		}
		
		if (Cigar > CigarLimit) goto bail;
		*Cigar = kCIGARTERMINATOR;
	}
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// COPY ALIGNMENT INFO TO TRANSFER BUFFER
	unsigned int map = Trial->DecisionTreeHistory | kSuccessfullyEncoded;
	if (IsAllN) map |= kAllIsN;
	if (ContainsN) map |= kContainsN;
	if (TagData->TagLen > TagData->AlignLen) map |= kTrimmed;
	
	extrslt->score               = Trial->score;
	extrslt->chr                 = Trial->chr;
	extrslt->chr_pos             = Trial->chr_pos + chr_pos_shift;
	extrslt->MismatchCount       = Trial->MismatchCount;
	extrslt->revNmer             = Trial->revNmer;
#ifdef DEBUG_TAG_MAPPING_PROCESS
	extrslt->AnchorPositionInTag = Trial->AnchorPositionInTag;
#endif
	extrslt->MismatchLength   = 1 + (unsigned char) ((uintptr_t) Mismatch - (uintptr_t) &(extrslt->Mismatch[0]));
	const unsigned int ciglen = 1 + (unsigned int) ((uintptr_t) Cigar - (uintptr_t) &(extrslt->Cigar[0]));
	assert(ciglen >= 3);
	if (ciglen != 3) 
		extrslt->CigarLength = ciglen;
	else { 
		assert(extrslt->Cigar[1] == 'M');
		extrslt->CigarLength = 1;
		extrslt->Cigar[0]    = kCIGARTERMINATOR;
		if (extrslt->MismatchLength == 1) map |= kPerfectMatch;
	}
	
	extrslt->DecisionTreeHistory = map;

	assert(extrslt->Cigar[0] != 0);
#ifndef NDEBUG
	if (Trial->chr > 0x3FF) {
		printf("We have a chromosome number %u !!!\n", Trial->chr);
	}
#endif
	return SoftClipSize;
	
	bail:;
#ifndef NDEBUG
	  char decisionString[16];
		{
			register char * cptr = &decisionString[0];
			for (int i=10; i>=0; i--) {
				char c;
				if (Trial->DecisionTreeHistory & (1 << i))
					c = '1';
				else 
					c = '0';
				*cptr++ = c;
				if ((i % 4) == 0) *cptr++ = ' ';
			}
			*cptr = '\0';
		}
		
		printf("Unable to encode %u mismatch with score %i, decision %u=(%s)...\n"
		       "Mismatch: total requirement %i, removed %i out of %i, count @ softclip location %u\n"
		       "Cigar: total requirement %i, removed %i out of %i, count @ softclip location %u\n"
		       "TAG: %.*s %c\n"
		       "STA: %s\n"
		       "CIG: ",
		       Trial->MismatchCount, Trial->score, Trial->DecisionTreeHistory, decisionString,
		       Trial->TotalRequiredMMCount, RemovedMismatch, MismatchToRemove, Trial->SoftClipEvictedMMSize,
		       Trial->TotalRequiredCigCount, RemovedCigar, CigarToRemove, Trial->SoftClipEvictedCigSize,
		       (int) TagData->AlignLen, Tag, (Trial->revNmer) ? (int)'<' : (int)'>', Trial->States);
		{
			unsigned char * restrict ptr = extrslt->Cigar;
			int c=0;
			while(*ptr != kCIGARTERMINATOR) {
				printf("%u%c ", ptr[0], ptr[1]);
				ptr += 2;
				c += 1;
				if (c >= kMaxEncodedCigar) break; 
			}
			fputs("\nMM : ", stdout);
			ptr = extrslt->Mismatch;
			c = 0;
			while(*ptr != kMMTERMINATOR) {
				printf("%u%c ", ptr[0], ptr[1]);
				ptr += 2;
				c += 1;
				if (c >= kMaxEncodedMismatch) break; 
			}
			fputs("\n\n", stdout);
		}
#else
		printf("Unable to encode %u mismatch with score %i\n", Trial->MismatchCount, Trial->score);
#endif
		fflush(stdout);
		
		extrslt->chr                 = kUNMAPPED;
		extrslt->Cigar[0]            = kCIGARTERMINATOR;
		extrslt->Mismatch[0]         = kMMTERMINATOR;
		extrslt->DecisionTreeHistory = 0U;
		return 0U;
}
/* vim: tabstop=2 shiftwidth=2
 */
