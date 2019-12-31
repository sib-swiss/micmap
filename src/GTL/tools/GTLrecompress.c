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

//=================================================================================================

#define _GNU_SOURCE     /* Expose declaration of tdestroy() */
#include "config.h"
#include <search.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <sys/time.h>

#include <sys/io.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>

#include <pthread.h>
#include <semaphore.h>
#include <assert.h>

#include <errno.h>
#include <stdatomic.h>

// basic compression
#include <zlib.h>

#include "virt_chr.h"
#include "GTL.h"
#include "gtlVersion.h"
#include "GTLdispatch.h"

//=================================================================================================
// Extra parameter type between dispatcher and threads to link working and storing jobs
struct ExtraParam {

};

//=================================================================================================
static int verbose = 0;
static int debug = 0;

//=================================================================================================

static void __attribute__((noreturn)) usage() 
{
		printf("usage:\n\n");
		printf("GTLrecompress [options] GTL files\n\n");
		printf("           -o path               : output base file name generation\n");
		printf("                                   (default is /tmp/original file name\n");
		printf("           -v level              : verbose level\n");
		printf("           -T <uint>             : number of threads\n");
		printf("\n");
		printf(GTL_VERSION_FULL_STRING "\n");
		printf("(c) Thierry Schuepbach & Cie 2014-2018\n");
		exit(1);
}
//---------------------------------------------------------------

static void* recompress(threadpool_t * const restrict thpool)
{
	char outFileName[PATH_MAX];
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Mark thread as alive (initialized)
	pthread_mutex_lock(&thpool->thcount_lock);
	thpool->num_threads_alive++;
	pthread_mutex_unlock(&thpool->thcount_lock);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Handle jobqueue
	while(thpool->threads_keepalive){
		bsem_wait(thpool->jobqueue_p.has_items);

		if (thpool->threads_keepalive){
			/* Signal I am working */
			pthread_mutex_lock(&thpool->thcount_lock);
			thpool->num_threads_working++;
			pthread_mutex_unlock(&thpool->thcount_lock);
			
			/* Read job from queue and execute it */
			pthread_mutex_lock(&thpool->jobqueue_p.rwmutex);
			GTLBlockDispatchFile_t * const restrict task = (GTLBlockDispatchFile_t*) jobqueue_pull(&thpool->jobqueue_p);
			pthread_mutex_unlock(&thpool->jobqueue_p.rwmutex);
			
			/*************************************************************************/
			/*                             START THE TASK                            */
			/*************************************************************************/
			if (task) {
				GTLBlockInfo_t * const restrict Block = &(task->Block);
		
				/*
				 * TotalJobs triggers either:
				 * - TotalJobs = zero for fetchers that will analyze and export tags
				 * - TotalJobs > 0 for last job to collect and send further to caller
				 */
				
				if (task->TotalJobs == 0 ) {
					/* Read and decompress the block from file */
					TLBDATA * const restrict td = (TLBDATA*) malloc(sizeof(TLBDATA)); 
					if(td == NULL) {
						fputs("Unable to allocate memory for TLBDATA\n", stderr);
						exit(1);
					}
					allocBlock(td, true);
					
					/* Save it for other to use */
					Block->whatever = td;
					
					if (preadDecompressBlock(Block->fd, Block->offset + sizeof(TLBHEADER), &(Block->thd), td, 0)) {
						fputs("Error decompressing block\n", stderr);
						exit(1);
					}
									
					/* Recompress it */
					if (compressBlock(td) == NULL) {
						fprintf(stderr, "Error in compressBlock\n");
						exit(1);
					}
					
					/* I have done my part , signal it */
					sem_post(task->semaphore);
				}
				else {
					
					/* Wait for the working threads to perform their task */
					for (unsigned int k=0; k<task->TotalJobs; k++) {
						int res;
						again:;
							res = sem_wait(task->semaphore);
							if (res == EINTR) goto again;
							else if (res != 0) {
								perror("sem_wait");
								exit(1);
							}
					}
					sem_destroy(task->semaphore);
					free(task->semaphore);
					task->semaphore = NULL;
					
					snprintf(outFileName, PATH_MAX, "%s/%s", (char*) task->params, task->FileName);
					
					for (unsigned int k=0; k<task->TotalJobs; k++) {
						/* Get the compressed block */
						GTLBlockInfo_t * const restrict block = &(task->nPreviousBlocks[k]->Block);
						
						/* Get its real data */
						TLBDATA * const restrict td = (TLBDATA*) block->whatever;
						printf("Appending %u bytes to %s\n", td->header.blockLength, outFileName);					
						int fd = open(outFileName, O_WRONLY|O_CREAT|O_APPEND, 0664);
						if (fd == -1) {
							perror("open:");
							exit(1);
						}
						if (write(fd,td->diskBuffer,td->header.blockLength) != td->header.blockLength) {
							perror("write failed");
							exit(1);
						}
						//fsync(fd);   // FIXME   was here before the O_SYNC
						close(fd);
						
						/* Free memory */
						td->cnt = 0;
						td->hdrSize = 0;
						td->qualSize = 0;
						for (unsigned int i = 0; i < kMAXnbMatches; i++) {
							td->mmSize[i] = 0;
							td->cigarSize[i] = 0;
						}
						td->p4Size = 0;
						freeBlock(td);
						free(td);
						
						block->whatever = NULL;
						
						/* Send back previous jobs to done queue */
						pthread_mutex_lock(&(thpool->donequeue_p.rwmutex));
						jobqueue_push(&(thpool->donequeue_p), (job_t*) task->nPreviousBlocks[k]);
						pthread_mutex_unlock(&(thpool->donequeue_p.rwmutex));
					}
					free((void*) task->nPreviousBlocks);
					task->nPreviousBlocks = NULL;
				
					/* Give back the memory task to the available queue */
					pthread_mutex_lock(&(thpool->donequeue_p.rwmutex));
					jobqueue_push(&(thpool->donequeue_p), (job_t*) task);
					pthread_mutex_unlock(&(thpool->donequeue_p.rwmutex));
				}
			}
			
			/* Signal I am no longer working */
			pthread_mutex_lock(&thpool->thcount_lock);
			thpool->num_threads_working--;
			pthread_mutex_unlock(&thpool->thcount_lock);
		}
	}

	/* Signal I am done */
	pthread_mutex_lock(&thpool->thcount_lock);
	thpool->num_threads_alive--;
	pthread_mutex_unlock(&thpool->thcount_lock);
	 
	return (void*) 0;
}
//---------------------------------------------------------------

int main (int argc, char **argv)
{
	const static char DefaultOut[] = "/tmp";
	const char * OutputBaseFileName = DefaultOut;
	const char * * InputGTLFiles = NULL;
	GTLList_t * restrict GTLInputFileList = NULL;
	threadpool_t * restrict thpool = NULL;
	cpu_set_t * affinities = NULL;
	int c, err=1;
	unsigned int nThreads = 4U;
	unsigned int UseAffinities = 0U;
	unsigned int nInputGTLFiles = 0U;
	
	gCompress.scoredef = CMP_LZ4;
	gCompress.def = CMP_BZ2;
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Process arguments
	opterr = 0;
	while ((c = getopt (argc, argv, "o:v:T:ah")) != -1)
	switch (c)
	{
		case 'o':
			OutputBaseFileName = optarg;
			break;
		case 'v':
			sscanf(optarg,"%d",&verbose);
			break;

		case 'T':
			{
				int dummy = atoi(optarg);
				if (dummy > 0) { 
					nThreads = (unsigned int) dummy;
				}
				else {
					fputs("Thread number must be positive\n", stderr);
					exit(1);
				}
			}
			break;
		case 'a':
			UseAffinities = 1;
			break;
		case 'h':
			usage();
			break;
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Check requirements

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Get input file names
	if (optind < argc) {
		nInputGTLFiles = argc - optind;
		if (verbose) fprintf(stderr,"%i files provided\n", argc - optind);
		InputGTLFiles = alloca(nInputGTLFiles*sizeof(char*));
		for (unsigned int iFile=0; iFile<nInputGTLFiles; iFile++) {
			InputGTLFiles[iFile] = argv[iFile+optind];
			if (verbose) fprintf(stderr,"Input file %2i is %s\n", iFile, InputGTLFiles[iFile]);
		}
	}
	else {
		fputs("Expected input files\n", stderr);
		usage();
	}
	
	if (verbose) fprintf(stderr, "Will use %u threads for fetching data and report types\n", nThreads);
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Prepare affinities
	if (UseAffinities) {
		const unsigned int nTotalThreads = nThreads;
		affinities = (cpu_set_t*) malloc(nTotalThreads*sizeof(cpu_set_t));
		if (affinities != NULL) {
			unsigned int iThread;
			const unsigned int maxN = nTotalThreads < (unsigned int) sysconf(_SC_NPROCESSORS_ONLN) ? nTotalThreads :  (unsigned int) sysconf(_SC_NPROCESSORS_ONLN);
			for (iThread=0; iThread<maxN; iThread++) {
					CPU_ZERO(&affinities[iThread]);
					CPU_SET(iThread, &affinities[iThread]);
			}
			
			while (iThread <nTotalThreads) {
				CPU_ZERO(&affinities[iThread]);
				for (int j = 0; j < CPU_SETSIZE; j++) CPU_SET(j, &affinities[iThread]);
				iThread++;
			}
		}
		else {
			fprintf(stderr, "Affinities could not be allocated, keeping on without...\n");
		}
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Prepare GTL files list
	GTLInputFileList = createGTLList(InputGTLFiles, nInputGTLFiles);
	if (GTLInputFileList == NULL) {
		fputs("Error while preparing GTL input files\n", stderr);
		goto bail;
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Create the thread pool with data of type GTLBlockDispatchFile_t
	thpool = createGTLThreadPool(recompress, affinities, sizeof(GTLBlockDispatchFile_t), nThreads, 0);
	
	if (thpool == NULL) {
		fputs("Thread pool creation error\n", stderr);
		goto bail;
	}
		
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Prepare Extra parameters

	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Dispatch files and its blocks
	err = dispatchGTLBlockFile(GTLInputFileList, thpool, -1, true, nThreads, (void*) OutputBaseFileName);
	
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Wait for jobs to terminate and empty queues
	thpool_wait(thpool);
	
	destroyGTLThreadPool(thpool);
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Clean up and quit

bail:;
	fflush(stdout);
	if (GTLInputFileList) destroyGTLList(GTLInputFileList);

	if (affinities) free(affinities);
	
	return err;

} /* main */
//---------------------------------------------------------------

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */

