/*
 * ------------------------------------------------------------------------------------------------------------------------
 *
 *                                * micmap *
 *      Mapping of short reads in fastq format onto a reference
 *
 *
 *  Copyright (C) SIB  - Swiss Institute of Bioinformatics,  2015-2019 Nicolas Guex, Thierry Schuepbach and Christian Iseli
 *  Copyright (C) UNIL - University of Lausanne, Switzerland      2019 Nicolas Guex and Christian Iseli
 *  Copyright (C) EPFL - EPFL, Lausanne, Switzerland              2022 Christian Iseli
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
#include "config.h" 
#include "constants.h"
#include <stdlib.h> 
#include <stdio.h> 
#include <sys/types.h> 
#include <sys/stat.h> 
#include <fcntl.h> 
#include <unistd.h> 
#include <string.h> 
#include <errno.h> 
#include <semaphore.h> 
#include <pthread.h> 
#include <assert.h>
#include <sys/mman.h>

#include "Reader.h" 

#define MAP_HUGE_1GB    (30 << MAP_HUGE_SHIFT)
#define MAP_HUGE_2MB    (21 << MAP_HUGE_SHIFT)

//--------------------------------------------------------------- 
// AUTO TRIMMING setup
//--------------------------------------------------------------- 
#define TRIM_PROBELEN 12
#define TRIM_MINLEN 50
#define TRIM_MAXERR 10
const unsigned short TRIM_OFFSET[3] = { 5, 10, 17 };

//--------------------------------------------------------------- 
// STRUCTURE DEFINITIONS 
//--------------------------------------------------------------- 
typedef struct _read_buf_t {
	unsigned int lc;
	unsigned int ic;
	unsigned int ip;
	char in[kFASTQReaderLineBufferSize];
} read_buf_t;

//--------------------------------------------------------------- 
// LOCAL FCTS 
//--------------------------------------------------------------- 
static void read_line_buf(read_buf_t * const  b, const int fd, char * const line,
                          const unsigned int lmax)
{
	b->lc = 0;
	while (1) {
		ssize_t rc;
		if (b->ic > b->ip) {
			char *ptr;
			unsigned int len;
			/* See what's left.  */
			if ((ptr = memchr(b->in + b->ip, '\n', (size_t) (b->ic - b->ip))) != NULL) {
				/* Ok, we have a full line.  */
				len = (unsigned int) (ptr - (b->in + b->ip) + 1);
				if (lmax < b->lc + len + 1) {
					fprintf(stderr, "buffer too small %d < %d\n", lmax, b->lc + len + 1);
					b->lc = 0;
					line[0] = 0;
					return;
				}
				assert(b->ip+len <= kFASTQReaderLineBufferSize);
				assert(b->lc+len < lmax);
				memcpy(line + b->lc, b->in + b->ip, (size_t) len);
				b->lc += len;
				line[b->lc] = 0;
				b->ip += len;
				return;
			}
			/* Copy what's left, then read more.  */
			len = b->ic - b->ip;
			if (lmax < b->lc + len + 1) {
				fprintf(stderr, "buffer too small %d < %d\n", lmax, b->lc + len + 1);
				b->lc = 0;
				line[0] = 0;
				return;
			}
			assert(b->ip+len <= kFASTQReaderLineBufferSize);
			assert(b->lc+len < lmax);
			memcpy(line + b->lc, b->in + b->ip, (size_t) len);
			b->lc += len;
		}
		b->ip = b->ic = 0;
		while ((rc = read(fd, b->in, kFASTQReaderLineBufferSize)) == -1) {
			if (errno != EINTR && errno != EAGAIN) {
				fprintf(stderr, "Could not read: %s(%d)\n", strerror(errno), errno);
				b->lc = 0;
				line[0] = 0;
				return;
			}
		}
		b->ic += (unsigned int) rc;
		if (b->ic == 0) {
			/* Got to the EOF...  */
			line[b->lc] = 0;
			return;
		}
	}
	
} /* read_line_buf */
//---------------------------------------------------------------

//--------------------------------------------------------------- 
// FUNCTIONS 
//--------------------------------------------------------------- 
void* readFASTQ(FASTQ * const fastq)
{
	read_buf_t rb1; // = { .lc=0, .ic=0, .ip=0 };
	read_buf_t rb2; // = { .lc=0, .ic=0, .ip=0 };
	char trash[kTagSize];
	int fd1, fd2;
	

	rb1.lc = 0; rb1.ic = 0; rb1.ip = 0;
	rb2.lc = 0; rb2.ic = 0; rb2.ip = 0;

	//////////////////////////////////////////////////////////////////////////////////////////////////
	// open tag file
	fd1 = open(fastq->tagfile1, O_RDONLY);
	if (fd1 == -1) {
		printf("Error opening %s : %s(%d)\n", fastq->tagfile1, strerror(errno), errno);
		return (void*) 1;
	}

	if (fastq->PairedEnd) {
		/* open tag file */
		fd2 = open(fastq->tagfile2, O_RDONLY);
		if (fd2 == -1) {
			printf("Error opening %s : %s(%d)\n", fastq->tagfile2, strerror(errno), errno);
			close(fd1);
			return (void*) 2;
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////
	// Loop on tags	
	unsigned int totalTags = 0U;
	unsigned long ordinal = fastq->startOrdinal;
	unsigned int tagCnt;
	ReaderJob_t * restrict work;
	while(1)
	{
		tagCnt = 0U;
		
		//////////////////////////////////////////////////////////////////////////////////////////////
		// Aquire empty memory slot
		do {
			bsem_wait(fastq->ReaderMemorySlot.has_items);
			pthread_mutex_lock(&(fastq->ReaderMemorySlot.rwmutex));
			work = (ReaderJob_t *) jobqueue_pull(&(fastq->ReaderMemorySlot));
			pthread_mutex_unlock(&(fastq->ReaderMemorySlot.rwmutex));
		} while (work == NULL);

		//////////////////////////////////////////////////////////////////////////////////////////////
		// Fill in the Block of tags
		ReaderTag_t * restrict CurrentTag = work->Tags;
		do 
		{
			char tmphdr[kHdrSize * 2];
			//////////////////////////////////////////////////////////////////////////////////////////////
			// FASTQ Header
			read_line_buf(&rb1, fd1, tmphdr, kHdrSize * 2);
			if (rb1.lc == 0) break;
			
			/* see if we have an original ordinal to keep */
			char *sptr = NULL;
			for (unsigned int j = 1; j <= 21; j++)
			{
				if (tmphdr[j] == '@')
				{
					sptr = tmphdr + j;
					break;
				}
			}
			if (tmphdr[0] == '@' && sptr != NULL && (sptr - (tmphdr + 1)) <= 20) {
				char *endptr;
				unsigned long val = strtoul(tmphdr + 1,&endptr,0);
				if (endptr == sptr) {
					ordinal = val;
					unsigned int newlen = rb1.lc - (sptr - tmphdr);
					memmove(tmphdr,sptr,newlen);
					rb1.lc = newlen;
				}
			}
			assert(rb1.lc - 1 < kHdrSize);
			tmphdr[rb1.lc - 1] = '\0';
			memcpy(CurrentTag->Header, tmphdr, rb1.lc);

			CurrentTag->HeaderLen = (unsigned short) (rb1.lc-1);
			CurrentTag->ordinal = ordinal;
			
			//////////////////////////////////////////////////////////////////////////////////////////////
			// FASTQ Tag
			read_line_buf(&rb1, fd1, (char *) CurrentTag->Tag, kTagSize);
			CurrentTag->AlignLen = (unsigned short) (rb1.lc-1);
			CurrentTag->TagLen = (unsigned short) (rb1.lc-1);
			CurrentTag->Tag[rb1.lc -1 ] = '\0';

			//////////////////////////////////////////////////////////////////////////////////////////////
			// FASTQ intermediate dummy line
			read_line_buf(&rb1, fd1, trash, kTagSize);

			//////////////////////////////////////////////////////////////////////////////////////////////
			// FASTQ Quality scores
			read_line_buf(&rb1, fd1, (char *) CurrentTag->Quality, kTagSize);
			CurrentTag->Quality[rb1.lc - 1] = '\0';

			// adjust qs to match our expactations (lowest score == '#')
			rb1.lc -= 1;
			if (fastq->MinimalQualityScore != '#') {
				const char Subs = fastq->MinimalQualityScore - '#';
				for (unsigned int i = 0; i < rb1.lc; i++) CurrentTag->Quality[i] -= Subs;
			}
			rb1.lc -= 1; // point to last quality score of tag.
			while(rb1.lc > 0) {
				if (CurrentTag->Quality[rb1.lc] != '#') break;
				rb1.lc--;
			}
			
			// store the length of the "not low quality" tag in penultimate char of tag.
			CurrentTag->ValidQualityLen = (unsigned short) (rb1.lc+1);

			tagCnt++;
			ReaderTag_t * const restrict PreviousTag = CurrentTag;
			CurrentTag++;

			//////////////////////////////////////////////////////////////////////////////////////////////
			// Paired End ?
			if (fastq->PairedEnd) {			
				int pos = 0;

				//////////////////////////////////////////////////////////////////////////////////////////////
				// FASTQ Header
				read_line_buf(&rb2, fd2, tmphdr, kHdrSize);
				assert(rb2.lc > 0 && rb2.lc - 1 < kHdrSize);
				tmphdr[rb2.lc - 1] = '\0';
				// identifies where second hdr starts to differ from first one.
				while(tmphdr[pos] == PreviousTag->Header[pos]) 
				{
					if (pos == (rb2.lc - 1))
						break;
					pos++; 
				};
				// check that first char to differ (pos) is actually the last caracter of both strings
				
				//printf("%c == %c ?   pos=%d  (rb2.lc - 2)=%d  %d\n",tmphdr[pos],PreviousTag->Header[pos],pos,(rb2.lc - 2),(PreviousTag->HeaderLen -1));
				if ((tmphdr[pos] == '2')
				    && (PreviousTag->Header[pos] == '1')
				    && (pos == (rb2.lc - 2))
				    && (pos == (PreviousTag->HeaderLen - 1)))
				{
					// store an empty string, and set length to zero.
					CurrentTag->Header[0] = 0; /* probably not necessary */
					CurrentTag->HeaderLen = 0;

					// delete the last character ('1') of hdr1.
					PreviousTag->Header[pos] = 0; /* probably not necessary */
					PreviousTag->HeaderLen -= 1;
				}
				else
				{
					memcpy(CurrentTag->Header, tmphdr, rb2.lc);
					CurrentTag->HeaderLen = (unsigned short) (rb2.lc-1);
				}

				// a more generic alternative could be to store only part that differ.
				// but we also need to say from which pos it starts to differ.	
				// or to explicitely parse a formatted known header from a list of known formats...	

				CurrentTag->ordinal = ordinal;
				/* put string length as last character (6 times below).  */

				//////////////////////////////////////////////////////////////////////////////////////////////
				// FASTQ Tag
				read_line_buf(&rb2, fd2, (char *) CurrentTag->Tag, kTagSize);
				CurrentTag->AlignLen = (unsigned short) (rb2.lc-1);
				CurrentTag->TagLen = (unsigned short) (rb2.lc-1);
				CurrentTag->Tag[rb2.lc - 1] = '\0';
				
				//////////////////////////////////////////////////////////////////////////////////////////////
				// FASTQ intermediate dummy line
				read_line_buf(&rb2, fd2, trash, kTagSize);
				
				//////////////////////////////////////////////////////////////////////////////////////////////
				// FASTQ quality scores
				read_line_buf(&rb2, fd2, (char *) CurrentTag->Quality, kTagSize);
				CurrentTag->Quality[rb2.lc - 1] = '\0';

				// adjust qs to match our expactations (lowest score == '#')
				rb2.lc -= 1;
				if (fastq->MinimalQualityScore != '#')	{
					const char Subs = fastq->MinimalQualityScore - '#';
					for (unsigned int i = 0; i < rb2.lc; i++) CurrentTag->Quality[i] -= Subs;
				}
				rb2.lc -= 1; // point to last quality score of tag.
				while(rb2.lc > 0) {
					if (CurrentTag->Quality[rb2.lc] != '#') break;
					rb2.lc--;
				}
				
				// store the length of the "not low quality" tag in penultimate char of tag.
				CurrentTag->ValidQualityLen = (unsigned short) (rb2.lc+1);
				// +++
				// Do automated trimming
				int idx = -1;
				unsigned short maxlen = PreviousTag->AlignLen;
				if (CurrentTag->AlignLen < maxlen)
					maxlen = CurrentTag->AlignLen;
				if (maxlen >= TRIM_MINLEN)
				{
					char RevTag[maxlen];
					for (unsigned int i = 0; i < maxlen; i++)
					{
						char c = CurrentTag->Tag[i];
						char r;
						switch (c)
						{
							case 'A':
								r = 'T';
								break;
							case 'C':
								r = 'G';
								break;
							case 'G':
								r = 'C';
								break;
							case 'T':
								r = 'A';
								break;
							default:
								r = 'N';
						}
						RevTag[maxlen - i - 1] = r;
					}
					for (unsigned int trial = 0; trial < 3; trial++)
					{
						for (unsigned int i = TRIM_OFFSET[trial]; i < maxlen - TRIM_PROBELEN; i++)
						{
							if (PreviousTag->Tag[TRIM_OFFSET[trial]] == RevTag[i])
							{
								if (memcmp(PreviousTag->Tag + TRIM_OFFSET[trial], RevTag + i, TRIM_PROBELEN) == 0)
								{
									int tmpIdx = i - TRIM_OFFSET[trial];
									unsigned int score = 0;
									const unsigned maxscore = maxlen - tmpIdx;
									for (unsigned int j = tmpIdx; j < maxlen; j++)
									{
										if (PreviousTag->Tag[j - tmpIdx] == RevTag[j])
											score += 1;
									}
									if (score + TRIM_MAXERR >= maxscore && maxscore >= TRIM_MINLEN)
									{
										idx = tmpIdx;
										break;
									}
								}
							}
						}
						if (idx >= 0)
							break;
					}
					// FIXME - when idx == 0 we do not need to trim, but there is no need to align the second read either since it will be the same as the first
					// ... conceivably, we could also take advantage of pairs where there is a large overlap to perform only 1 alignment ...
					if (idx > 0)
					{
						// we need to trim
						// FIXME - do we really count on the '\0' at the end of the tag and quality strings ?
						PreviousTag->AlignLen = maxlen - idx;
						CurrentTag->AlignLen = maxlen - idx;
						PreviousTag->Tag[PreviousTag->TagLen] = '\0';
						CurrentTag->Tag[CurrentTag->TagLen] = '\0';
						PreviousTag->Quality[PreviousTag->TagLen] = '\0';
						CurrentTag->Quality[CurrentTag->TagLen] = '\0';
						if (PreviousTag->ValidQualityLen > PreviousTag->AlignLen)
							PreviousTag->ValidQualityLen = PreviousTag->AlignLen;
						if (CurrentTag->ValidQualityLen > CurrentTag->AlignLen)
							CurrentTag->ValidQualityLen = CurrentTag->AlignLen;
					}
				}
				// ---
				CurrentTag++;	
				tagCnt++;
			} // fastq->PairedEnd
			
			//////////////////////////////////////////////////////////////////////////////////////////////
			// Increase ordinal in case there is no original one
			ordinal++;
		} while (tagCnt < kMaxTagsPerBuffer);
		
		/* Set the number of tags in this block */
		work->TagCnt = tagCnt;
		totalTags += tagCnt;

		//////////////////////////////////////////////////////////////////////////////////////////////
		// Send Block to Process queue
		//printf("%s pushing new batch\n", __FUNCTION__);
		pthread_mutex_lock(&(fastq->ReaderJobQueue.rwmutex));
		jobqueue_push(&(fastq->ReaderJobQueue), (job_t*) work);
		pthread_mutex_unlock(&(fastq->ReaderJobQueue.rwmutex));

                /////////////////////////////////////////////////////////////////////////////////////////////
		// Is it done?
		if (tagCnt < kMaxTagsPerBuffer) {
			printf("Fastq reader is done reading, overall %u tags were read\n", totalTags);
			break;
		}
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////
	// Send End of Sequence Information
	if (tagCnt > 0) {
		do {
			bsem_wait(fastq->ReaderMemorySlot.has_items);
			pthread_mutex_lock(&(fastq->ReaderMemorySlot.rwmutex));
			work = (ReaderJob_t *) jobqueue_pull(&(fastq->ReaderMemorySlot));
			pthread_mutex_unlock(&(fastq->ReaderMemorySlot.rwmutex));
		} while (work == NULL);
		work->TagCnt = 0U;

		printf("%s pushing End of Sequence\n", __FUNCTION__);
		pthread_mutex_lock(&(fastq->ReaderJobQueue.rwmutex));
		jobqueue_push(&(fastq->ReaderJobQueue), (job_t*) work);
		pthread_mutex_unlock(&(fastq->ReaderJobQueue.rwmutex));	
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////
	// close tag file
	close(fd1);
	if (fastq->PairedEnd) close(fd2);

	//////////////////////////////////////////////////////////////////////////////////////////////////
	// Summary
	if (fastq->PairedEnd) printf("Total read %u pair of tags\n", totalTags/2);
	
	return (void*) 0;

} /* readFASTQ */
//---------------------------------------------------------------

int allocateReaderMemory(FASTQ * const restrict ReaderInputs, const unsigned int nBlocks)
{
	ReaderMemoryBulk_t * const restrict Memory = &ReaderInputs->Memory;
	
	const size_t TLBSize = 1UL << 21;
	size_t MemorySize = nBlocks*sizeof(ReaderTag_t)*kMaxTagsPerBuffer;
	MemorySize = (MemorySize + (TLBSize-1UL)) & ~(TLBSize-1UL);
	
	Memory->Ptr = (char *) mmap(NULL, MemorySize, PROT_READ | PROT_WRITE,
	                            MAP_PRIVATE | MAP_ANONYMOUS | MAP_POPULATE | MAP_HUGETLB,
	                            -1, 0);
	if (Memory->Ptr == MAP_FAILED) {
		fprintf(stderr, "%s: Unable to allocate %lu Mbytes Huge TLB, going for standard!\n", __FUNCTION__, MemorySize >> 20);

		if ( posix_memalign((void**) &(Memory->Ptr), TLBSize, MemorySize) == 0) {
			//if (madvise(Memory->Ptr, MemorySize, MADV_HUGEPAGE)) perror("madvise");
			Memory->IsHugeTLB = 0;
		}
		else {
			perror("posix_memalign");
			fputs("FASTQ Reader unable to allocate memory!\n", stderr);
			fflush(stderr);
			goto bail;
		}
	}
	else 
		Memory->IsHugeTLB = 1;
	Memory->Size = MemorySize;
	
	ReaderJob_t * Jobs = (ReaderJob_t*) malloc(nBlocks*sizeof(ReaderJob_t));
	if (Jobs == NULL) goto bail;
	
	jobqueue_init(&ReaderInputs->ReaderMemorySlot);
	jobqueue_init(&ReaderInputs->ReaderJobQueue);
	
	pthread_mutex_lock(&ReaderInputs->ReaderMemorySlot.rwmutex);
	for (size_t iBlock=0; iBlock<nBlocks; iBlock++) {
		Jobs[iBlock].TagCnt = 0U;
		Jobs[iBlock].Tags = (ReaderTag_t*) (Memory->Ptr + iBlock*kMaxTagsPerBuffer*sizeof(ReaderTag_t));
		Jobs[iBlock].prev = NULL;
		jobqueue_push(&(ReaderInputs->ReaderMemorySlot), (job_t*) &Jobs[iBlock]);
	}
	pthread_mutex_unlock(&ReaderInputs->ReaderMemorySlot.rwmutex);
	
	Memory->Jobs = Jobs;
	
	return 0;
	
	bail: ;
	freeReaderMemory(ReaderInputs);
	
	return 1;
}
//---------------------------------------------------------------

void freeReaderMemory(FASTQ * const restrict ReaderInputs)
{
	ReaderMemoryBulk_t * const restrict Memory = &ReaderInputs->Memory; 
	if (Memory->Ptr) {
		if (Memory->IsHugeTLB)
			munmap(Memory->Ptr, Memory->Size);
		else 
			free(Memory->Ptr);
	}
	
	if (Memory->Jobs) {
		jobqueue_destroy(&ReaderInputs->ReaderMemorySlot);
		jobqueue_destroy(&ReaderInputs->ReaderJobQueue);
		
		free(Memory->Jobs);
	}
}
//---------------------------------------------------------------

/* vim: tabstop=2 shiftwidth=2
 */
