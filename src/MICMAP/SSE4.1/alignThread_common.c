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
//========================================================================================================
// Common static include function
#include "../alignThread_common.h"

//========================================================================================================
// Queues
static inline void __attribute__((always_inline))
SendBufferToHost(jobqueue_t * const restrict ToBeWritten_q, cpuSpace_t * const restrict space,
								 jobqueue_t * const restrict BelongsTo_q, const unsigned long Size)
{
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Make sure everyone has finished writing
	while (atomic_load(&(space->Stored)) < Size) { __asm__ __volatile__ ("pause" ::: "memory"); }
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Put the job into to be written queue
	space->length = (size_t) Size;
	assert(BelongsTo_q->has_items != NULL);
	space->BelongsTo_q = BelongsTo_q;
	//pthread_yield();
	//printf("space->BelongsTo_q = %p\n",space->BelongsTo_q);
	//printf("  prev : %p; data : %p; BelongsTo_q : %p; length : %zu\n",space->prev,space->OutBufferPtr,space->BelongsTo_q,space->length);
	pthread_mutex_lock(&(ToBeWritten_q->rwmutex));
	assert(space->BelongsTo_q == BelongsTo_q);
	jobqueue_push(ToBeWritten_q, (job_t*) space);
	pthread_mutex_unlock(&(ToBeWritten_q->rwmutex));
	
}
//---------------------------------------------------------------

static inline int __attribute__((always_inline))
recvBufferFromHost(AlignerPool_t * const restrict AlignPool, cpuAlignerJob_t * const restrict task)
{
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Make sure everyone has finished reading
	const unsigned long TagCnt = task->TagCnt;
	while (atomic_load(&(task->Stored)) < TagCnt) { __asm__ __volatile__ ("pause" ::: "memory"); }
	
 	////////////////////////////////////////////////////////////////////////////////////////////////////////
 	// Get a new job to align
	pthread_mutex_lock(&(AlignPool->thcount_lock));
	AlignPool->num_threads_working--;
	pthread_mutex_unlock(&(AlignPool->thcount_lock));
	bsem_wait(AlignPool->ToBeAligned_q->has_items);
	pthread_mutex_lock(&(AlignPool->ToBeAligned_q->rwmutex));
	register DecoderJobElement_t * const restrict job = (DecoderJobElement_t *) jobqueue_pull(AlignPool->ToBeAligned_q);
	pthread_mutex_unlock(&(AlignPool->ToBeAligned_q->rwmutex));
	
	if (job != NULL) {
		pthread_mutex_lock(&(AlignPool->thcount_lock));
		AlignPool->num_threads_working++;
		pthread_mutex_unlock(&(AlignPool->thcount_lock));
		task->Anchors = (const DecoderElement_t (*)[16]) job->Decoder.DecoderJob->Anchors;
		task->TagsData = job->Reader.ReaderJob->Tags;
		task->DecoderJob = job;
		task->TagCnt = job->Reader.ReaderJob->TagCnt;
		if (verbose & 0x1)
		  printf("JOB IN: ReaderJob @ 0x%lx, DecoderJob @ 0x%lx\n", (uintptr_t) job->Decoder.DecoderJob, (uintptr_t) job->Reader.ReaderJob);
		__asm __volatile( "sfence" );
		atomic_store(&(task->Stored), 0UL);
		atomic_store(&(task->Treated), 0UL);
		return 0;
	}
	else {
		return 1;
	}
}
//---------------------------------------------------------------

static inline void __attribute__((always_inline))
SendSpaceBufferBack(AlignerBlockArgs_t * const restrict AlignerBlock, cpuAlignerJob_t * const restrict task)
{
	DecoderJobElement_t * const restrict job = task->DecoderJob;
	assert(job);
 
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Make sure everyone has finished reading
	const unsigned long TagCnt = task->TagCnt;
	while (atomic_load(&(task->Stored)) < TagCnt) { __asm__ __volatile__ ("pause" ::: "memory"); }
	 
	/* Put the memory back to its corresponding queue */
	if (AlignerBlock->ReaderMemoryQueue) {
		pthread_mutex_lock(&(AlignerBlock->ReaderMemoryQueue->rwmutex));
		jobqueue_push(AlignerBlock->ReaderMemoryQueue, (job_t*) job->Reader.ReaderJob);
		pthread_mutex_unlock(&(AlignerBlock->ReaderMemoryQueue->rwmutex));
	}
	
	if (AlignerBlock->DecoderMemoryQueue)	{
		pthread_mutex_lock(&(AlignerBlock->DecoderMemoryQueue->rwmutex));
		jobqueue_push(AlignerBlock->DecoderMemoryQueue, (job_t*) job->Decoder.DecoderJob);
		pthread_mutex_unlock(&(AlignerBlock->DecoderMemoryQueue->rwmutex));
	}

	//printf("%s : returning counter @ 0x%lx, count %u\n", __FUNCTION__, job->count, atomic_load(job->count));	
	pthread_mutex_lock(&(AlignerBlock->DecoderInternalJob_q->rwmutex));
	jobqueue_push(AlignerBlock->DecoderInternalJob_q, (job_t*) job);
	pthread_mutex_unlock(&(AlignerBlock->DecoderInternalJob_q->rwmutex));
}
//---------------------------------------------------------------

/* vim: tabstop=2 shiftwidth=2
 */
