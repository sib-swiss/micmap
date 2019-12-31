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
#include "../constants.h"
#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <sys/mman.h>
#include <mpi.h>
#include <assert.h>
#include "Genome.h"
#include "global_var.h"
#include "mpi_transfer.h"
#include "Decoder.h"

static const char MPIGenomeName[] = "MPI distributed genome";

int SendGenome(Genome_t * const restrict genome)
{
	int err = -1;
		
	MPI_Bcast(&(genome->configBufferSize), 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
	MPI_Bcast(&(genome->configBuffer), genome->configBufferSize, MPI_CHARACTER, 0, MPI_COMM_WORLD);
	MPI_Bcast(genome->table, kGENOME_DATA_SIZE>>2, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	
	err = 0;
bail:
	return err;
}
//---------------------------------------------------------------

int ReceiveGenome(Genome_t * const restrict genome)
{
	int err = -1;
	genome->FileName = MPIGenomeName;
	
	MPI_Bcast(&(genome->configBufferSize), 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
	genome->configBuffer = malloc(genome->configBufferSize);
	if (genome->configBuffer == NULL) { 
		perror("malloc"); goto bail;
	}
	MPI_Bcast(&(genome->configBuffer), genome->configBufferSize, MPI_CHARACTER, 0, MPI_COMM_WORLD);
	
	genome->tableSize = kGENOME_DATA_SIZE;
	genome->table = (unsigned char *) mmap (NULL, kGENOME_DATA_SIZE, PROT_READ, MAP_PRIVATE, -1, 0);
	if (genome->table == MAP_FAILED) {
		printf("%s: mapping failed.\n", __FUNCTION__);
		goto bail;
	}
	MPI_Bcast(genome->table, kGENOME_DATA_SIZE>>2, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	
	// parse the data
	genome->virtchr = calloc(kMAXCHRcnt, sizeof(VIRTUALCHR));
	if (genome->virtchr == NULL) {
		fprintf(stderr,"%s: error allocating virtual chromosomes\n", __FUNCTION__);
		goto bail;
	}

	genome->virtchr[0].chr = 0;
	genome->virtchr[0].offset = 0;
	genome->virtchr[0].len = 0;
	strcpy(genome->virtchr[0].AC,"unknown");
	strcpy(genome->virtchr[0].SAMname,"unknown");

	parseVirtualChromosomes(genome->configBuffer, genome->configBufferSize, genome->virtchr);
	
	err = 0;
bail:
	return err;
}
//---------------------------------------------------------------

void* MPIDispatch(AlignerBlockArgs_t * const restrict data)
{
	ThreadItem_t Me;
	AlignerPool_t * const restrict Pool = data->Pool;
	MPI_Status status;
	MPI_Datatype ReaderType, DecoderType;
	int trigger;
	int err = 1;
	printf("This is MPI master dispatcher\n"); fflush(stdout);	
	//////////////////////////////////////////////////////////////////////////////////////////////////
	// Get my PID and add it to the Aligner pool, increase active thread counter
	{
		Me.thread = pthread_self();
		pthread_mutex_lock(&(Pool->threads_q.rwmutex));
		jobqueue_push(&(Pool->threads_q), (job_t*) &Me);
		pthread_mutex_unlock(&(Pool->threads_q.rwmutex));
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////
	// MPI Data types
	MPI_Type_contiguous((int) sizeof(ReaderTag_t), MPI_CHAR, &ReaderType);
	MPI_Type_contiguous((int) sizeof(DecoderElement_t)*kDesiredHitCnt, MPI_CHAR, &DecoderType);
	MPI_Type_commit(&ReaderType);
	MPI_Type_commit(&DecoderType);
		
	//////////////////////////////////////////////////////////////////////////////////////////////////
	// Mark thread as alive (initialized)
	pthread_mutex_lock(&Pool->thcount_lock);
	Pool->num_threads_alive += 1;
	pthread_mutex_unlock(&Pool->thcount_lock);
	
	//////////////////////////////////////////////////////////////////////////////////////////////////
	// Accepting connection and sending data
	memqueue_t * const restrict ReaderMemoryQueue = data->ReaderMemoryQueue;
	memqueue_t * const restrict DecoderMemoryQueue = data->DecoderMemoryQueue;
	memqueue_t * const restrict DecoderInternalMemoryQueue = data->DecoderInternalJob_q;

	while (1) {
		MPI_Recv(&trigger, 1, MPI_INTEGER, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
				
		int SuccessFullySent = 0;
		if (Pool->threads_keepalive) {
			bsem_wait(Pool->ToBeAligned_q->has_items);
			
			if (Pool->threads_keepalive) {
				
				/* Annonce thread working */
				pthread_mutex_lock(&Pool->thcount_lock);
				Pool->num_threads_working++;
				pthread_mutex_unlock(&Pool->thcount_lock);

				/* Read job from queue and execute it */
				pthread_mutex_lock(&Pool->ToBeAligned_q->rwmutex);
				const DecoderJobElement_t * const restrict task = (DecoderJobElement_t *) jobqueue_pull(Pool->ToBeAligned_q);
				pthread_mutex_unlock(&Pool->ToBeAligned_q->rwmutex);
				printf("MPI sender populating nodes %i with %u tags\n", status.MPI_SOURCE, task->Reader.ReaderJob->TagCnt);
				/* --------------------------------------- Tag count ---------------------------------------- */
				const ReaderJob_t * const restrict ReaderData = task->Reader.ReaderJob;
				MPI_Send(&(ReaderData->TagCnt), 1, MPI_UNSIGNED, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
				
				/* ------------------------------------------ Tags ------------------------------------------ */
				assert(ReaderData->TagCnt);
				MPI_Send(ReaderData->Tags, (int) (ReaderData->TagCnt), ReaderType, status.MPI_SOURCE, 2, MPI_COMM_WORLD);
				
				/* Put the memory back to its corresponding queue */
				pthread_mutex_lock(&(ReaderMemoryQueue->rwmutex));
				jobqueue_push(ReaderMemoryQueue, (job_t*) ReaderData);
				pthread_mutex_unlock(&(ReaderMemoryQueue->rwmutex));

				/* ----------------------------------------- Decoder ---------------------------------------- */
				const DecoderMemoryElement_t * const restrict DecoderData = task->Decoder.DecoderJob;
				MPI_Send(DecoderData->Anchors, (int) (ReaderData->TagCnt), DecoderType, status.MPI_SOURCE, 3, MPI_COMM_WORLD);
								
				/* Put the memory back to its corresponding queue */
				pthread_mutex_lock(&(DecoderMemoryQueue->rwmutex));
				jobqueue_push(DecoderMemoryQueue, (job_t*) DecoderData);
				pthread_mutex_unlock(&(DecoderMemoryQueue->rwmutex));
				
				pthread_mutex_lock(&(DecoderInternalMemoryQueue->rwmutex));
				jobqueue_push(DecoderInternalMemoryQueue, (job_t*) task);
				pthread_mutex_unlock(&(DecoderInternalMemoryQueue->rwmutex));
				
				pthread_mutex_lock(&Pool->thcount_lock);
				Pool->num_threads_working--;
				pthread_mutex_unlock(&Pool->thcount_lock);
				
				SuccessFullySent = 1;
			}
		}
		
		/* Annonce end of data */
		if (!Pool->threads_keepalive && (SuccessFullySent == 0)) {
			unsigned int uitmp = 0U;
			MPI_Send(&uitmp, 1, MPI_UNSIGNED, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
			printf("MPI dispatcher closing node %i..\n", status.MPI_SOURCE);
			break;
		}
	}
	printf("MPI Dispatcher sent everything, now closing...\n"); fflush(stdout);
	//////////////////////////////////////////////////////////////////////////////////////////////////
	// Send End Of Sequence to MPI aligner nodes
	unsigned int OK = 0U;
	int TotalMPISlave;
	if (MPI_Comm_size(MPI_COMM_WORLD, &TotalMPISlave) != 0) {
		fprintf(stderr, "Unable to get MPI comm size\n");
		goto bail;
	}
	TotalMPISlave -= 2;
	for (int i=0; i<TotalMPISlave; i++) {
		MPI_Recv(&trigger, 1, MPI_INTEGER, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
		printf("MPI dispatcher closing node %i..\n", status.MPI_SOURCE);
		MPI_Send(&OK, 1, MPI_UNSIGNED, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////
	// Terminate
	
	/* Annonce thread about to end */
	pthread_mutex_lock(&Pool->thcount_lock);
	Pool->num_threads_alive--;
	pthread_mutex_unlock(&Pool->thcount_lock);
	
	err = 0;
bail: ;
	MPI_Type_free(&ReaderType);
	MPI_Type_free(&DecoderType);
	return (void*) ((uintptr_t) err);
}
//---------------------------------------------------------------

TransferMemoryBulk_t * allocateTransferMemory(const size_t nBlocks)
{
	printf("This is allocate transfer memory for %lu blocks\n", nBlocks); fflush(stdout);
	TransferMemoryBulk_t * const restrict data = calloc(1UL, sizeof(TransferMemoryBulk_t));
	if (data == NULL) {
		fputs("Unable to allocate transfer memory...\n", stderr);
		goto bail;
	}

	const size_t MemorySize = (nBlocks*kMaxTagsPerBuffer*(sizeof(ReaderTag_t) + sizeof(DecoderElement_t)*kDesiredHitCnt) + 63UL + 4095UL) & ~(4095UL);
	data->ReaderDecoderMemory = (char *) mmap(NULL, MemorySize, PROT_READ | PROT_WRITE,
	                                          MAP_PRIVATE | MAP_ANONYMOUS | MAP_POPULATE | MAP_HUGETLB,
	                                          -1, 0);
	if (data->ReaderDecoderMemory == MAP_FAILED) {
		fprintf(stderr, "%s: Unable to allocate %lu Mbytes Huge TLB, going for standard!\n", __FUNCTION__, MemorySize >> 20);
		data->ReaderDecoderMemory = (char *) mmap(NULL, MemorySize, PROT_READ | PROT_WRITE,
	                                            MAP_PRIVATE | MAP_ANONYMOUS | MAP_POPULATE,
	                                            -1, 0);
		if (data->ReaderDecoderMemory == MAP_FAILED) goto fin3;
	}
	data->MemorySize = MemorySize;
	ReaderJob_t * const ReaderJobs = (ReaderJob_t*) calloc(nBlocks, sizeof(ReaderJob_t));
	DecoderMemoryElement_t * const DecoderElements = (DecoderMemoryElement_t*) calloc(nBlocks, sizeof(DecoderMemoryElement_t));
	data->jobs = (DecoderJobElement_t*) calloc(nBlocks, sizeof(DecoderJobElement_t));
	if (data->jobs == NULL || ReaderJobs == NULL || DecoderElements == NULL) goto fin2;
	
	jobqueue_init(&(data->Memory_q));
	jobqueue_init(&(data->ToBeAligned_q));
	
	{
		ReaderTag_t * restrict ReaderMemory = (ReaderTag_t*) data->ReaderDecoderMemory;
		DecoderElement_t (* restrict DecoderMemory)[kDesiredHitCnt] = (DecoderElement_t (*)[kDesiredHitCnt]) ( ( (uintptr_t) data->ReaderDecoderMemory + nBlocks*kMaxTagsPerBuffer*sizeof(ReaderTag_t) + 63UL) & ~(63UL) );
		for (size_t i=0; i<nBlocks; i++) {
			DecoderElements[i].Anchors = DecoderMemory;
			ReaderJobs[i].Tags = ReaderMemory;
			data->jobs[i].Reader.ReaderJob = &ReaderJobs[i];
			data->jobs[i].Decoder.DecoderJob = &DecoderElements[i];
			jobqueue_push(&(data->Memory_q), (job_t*) &(data->jobs[i]));
			ReaderMemory += kMaxTagsPerBuffer;
			DecoderMemory += kMaxTagsPerBuffer;
		}
	}
	
	return data;
fin1:
	free(data->jobs);
fin2:
	munmap(data->ReaderDecoderMemory, MemorySize);
fin3:
	free(data);
bail:
	return NULL;
}
//---------------------------------------------------------------

void* MPIReceiver(TransferMemoryBulk_t * const restrict data)
{
	DecoderJobElement_t * restrict work;
	MPI_Datatype ReaderType, DecoderType;
	MPI_Status Status;
	printf("This is MPIReceiver\n");fflush(stdout);	
	
	//////////////////////////////////////////////////////////////////////////////////////////////////
	// MPI Data types
	MPI_Type_contiguous((int) sizeof(ReaderTag_t), MPI_CHAR, &ReaderType);
	MPI_Type_contiguous((int) sizeof(DecoderElement_t)*kDesiredHitCnt, MPI_CHAR, &DecoderType);
	MPI_Type_commit(&ReaderType);
	MPI_Type_commit(&DecoderType);
	
	//////////////////////////////////////////////////////////////////////////////////////////////////
	// Get jobs data
	while(1) {
	
		//////////////////////////////////////////////////////////////////////////////////////////////
		// Aquire empty memory slot
		do {
			bsem_wait(data->Memory_q.has_items);
			pthread_mutex_lock(&(data->Memory_q.rwmutex));
			work = (DecoderJobElement_t *) jobqueue_pull(&(data->Memory_q));
			pthread_mutex_unlock(&(data->Memory_q.rwmutex));
			assert(work);
		} while (work == NULL);
		
		printf("MPI receiver querying master for task\n"); fflush(stdout);	
		//////////////////////////////////////////////////////////////////////////////////////////////
		// Acquire master attention
		MPI_Send(&(data->MPINodeId), 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD);
		//////////////////////////////////////////////////////////////////////////////////////////////
		// Get master node tag count
		MPI_Recv(&(work->nToBeProcessed), 1, MPI_UNSIGNED, 0, 1, MPI_COMM_WORLD, &Status);
		if (work->nToBeProcessed == 0U) break;
		//////////////////////////////////////////////////////////////////////////////////////////////
		// Fill in the Block of tags
		MPI_Recv((void*) work->Reader.ReaderJob->Tags, work->nToBeProcessed, ReaderType, 0, 2, MPI_COMM_WORLD, &Status);
		work->Reader.ReaderJob->TagCnt = work->nToBeProcessed;
		//////////////////////////////////////////////////////////////////////////////////////////////
		// Fill in the Block of decoding anchors
		MPI_Recv(&(work->Decoder.DecoderJob->Anchors[0][0]), work->nToBeProcessed, DecoderType, 0, 3, MPI_COMM_WORLD, &Status);	
		//////////////////////////////////////////////////////////////////////////////////////////////
		// Push to aligner job queue
		pthread_mutex_lock(&(data->ToBeAligned_q.rwmutex));
		jobqueue_push(&(data->ToBeAligned_q), (job_t*) work);
		pthread_mutex_unlock(&(data->ToBeAligned_q.rwmutex));
		printf("MPI node %i got a batch of %u tags\n", data->MPINodeId, work->nToBeProcessed);
		
	}
	printf("MPI node %i receiver ending...\n", data->MPINodeId);
	return (void*) ((uintptr_t) 0);
bail:
	MPI_Type_free(&ReaderType);
	MPI_Type_free(&DecoderType);
	return (void*) ((uintptr_t) 1);
}
//---------------------------------------------------------------

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
