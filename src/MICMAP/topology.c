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
#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#ifdef USE_AFFINITY
#include <unistd.h>
#include <sched.h>
#endif
#include "constants.h"
#include "json.h"
#include "topology.h"


Topology_t gTopology = {
#ifdef USE_AFFINITY
	.Reader = NULL,
	.MicAligners = NULL,
	.Aligners = NULL,
	.Writers = NULL,
	.Decoders = NULL,
#endif

	.nDecoders = 4U, 
	.nAligners = 1U,
	.nAlignerThreads = 4U,
	.nWriters = 1U,
	.nWriterThreads = 3U,
	.nReaderBlocks = 3U,
	.nDecoderBlocks = 3U,
	.nAlignerInBlocks = 3U,
	.nAlignerOutBlocks = 3U,
	.nAlignerOutBlockSize = 64UL * 1024UL * 1024UL
};

static void printfCPUThreadID(topology_mask_t Ids[3])
{
	fputs("\tSockets:", stdout);
	for (unsigned int j=0; j<8*sizeof(topology_mask_t); ++j) {
		const topology_mask_t mask = (topology_mask_t) 1 << ((sizeof(topology_mask_t)*8-1) - j);
		if ( mask & Ids[0] )
			fputc('1',stdout);
		else
			fputc('0', stdout);
	}
	fputc('\n',stdout);
	fputs("\tCores  :", stdout);
	for (unsigned int j=0; j<8*sizeof(topology_mask_t); ++j) {
		const topology_mask_t mask = (topology_mask_t) 1 << ((sizeof(topology_mask_t)*8-1) - j);
		if ( mask & Ids[1] )
			fputc('1',stdout);
		else
			fputc('0', stdout);
	}
	fputc('\n',stdout);
	fputs("\tThreads:", stdout);
	for (unsigned int j=0; j<8*sizeof(topology_mask_t); ++j) {
		const topology_mask_t mask = (topology_mask_t) 1 << ((sizeof(topology_mask_t)*8-1) - j);
		if ( mask & Ids[2] )
			fputc('1',stdout);
		else
			fputc('0', stdout);
	}
	fputc('\n',stdout);
}
//---------------------------------------------------------------

void printJSON(const SystemInfo * const restrict info, FILE * const restrict fd)
{
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Temporary allocations
	const unsigned int nThreads = info->nOverallCores/(info->nSockets*info->nCores);
	char * const Sockets = (char *) malloc((1+info->nSockets*2)*sizeof(char));
	char * Cores = (char *) malloc((1+info->nCores*3)*sizeof(char));
	char * Threads = (char *) malloc((1+nThreads*2)*sizeof(char));
	
	const unsigned int TotalMicRecvCnt = 0;
	
	if (Sockets == NULL || Cores == NULL || Threads == NULL) {
		fputs("printJSON unable to allocate memory\n", stderr);
		if (Sockets) free(Sockets);
		if (Cores) free(Cores);
		if (Threads) free(Threads); 
		return;
	}
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Prepare CPU Socket, core and thread numbering
	{
		char * restrict ptr = Sockets;
		for(unsigned int iSocket=0; iSocket< info->nSockets; iSocket++) {
			ptr += sprintf(ptr,"%1.1u,", iSocket);
		}
		ptr[-1] = '\0';
	}
	{
		char * restrict ptr = Cores;
		for(unsigned int iCore=0; iCore< info->nCores; iCore++) {
			if (iCore<10)
				ptr += sprintf(ptr,"%1u,", iCore);
			else
				ptr += sprintf(ptr,"%2u,", iCore);
		}
		ptr[-1] = '\0';
	}
	{
		char * restrict ptr = Threads;
		for(unsigned int iThread=0; iThread<nThreads; iThread++) {
			ptr += sprintf(ptr,"%1u,", iThread);
		}
		ptr[-1] = '\0';
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// HEADER
	fprintf(fd, "/* Automatic template generated for hostname %s */\n", info->Nodename);
	fputs("{\n", fd);
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// HOST
	fputs("\t\"HOST\":{\n", fd);
	{
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Reader
		fputs("\t\t\"Reader\":{\n", fd);
		fprintf(fd,"\t\t\t\"Input Blocks\": %u,\n", gTopology.nReaderBlocks);
		fprintf(fd,"\t\t\t\"Sockets\":[%s], \"Cores\":[%s], \"Threads\":[%s]\n", Sockets, Cores, Threads);
		fputs("\t\t},\n", fd);
		
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Writers
		fputs("\t\t\"Writers\":[\n", fd);
		for (unsigned int i=1U; i<gTopology.nWriters;i++) { 
			fputs("\t\t\t{\n", fd);
			fprintf(fd,"\t\t\t\t\"Master\":{ \"Sockets\":[%s], \"Cores\":[%s], \"Threads\":[%s] },\n", Sockets, Cores, Threads);
			fputs("\t\t\t\t\"Flush Buffer Threads\":[\n", fd);
			for (unsigned int j=1U; j<gTopology.nWriterThreads; j++) {
				fprintf(fd,"\t\t\t\t\t{ \"Sockets\":[%s], \"Cores\":[%s], \"Threads\":[%s] },\n", Sockets, Cores, Threads);
			}
			fprintf(fd,"\t\t\t\t\t{ \"Sockets\":[%s], \"Cores\":[%s], \"Threads\":[%s] }\n", Sockets, Cores, Threads);
			fputs("\t\t\t\t]\n", fd);
			fputs("\t\t\t},\n", fd);
		}
		fputs("\t\t\t{\n", fd);
		fprintf(fd,"\t\t\t\t\"Master\":{ \"Sockets\":[%s], \"Cores\":[%s], \"Threads\":[%s] },\n", Sockets, Cores, Threads);
		fputs("\t\t\t\t\"Flush Buffer Threads\":[\n", fd);
		for (unsigned int j=1U; j<gTopology.nWriterThreads; j++) {
			fprintf(fd,"\t\t\t\t\t{ \"Sockets\":[%s], \"Cores\":[%s], \"Threads\":[%s] },\n", Sockets, Cores, Threads);
		}
		fprintf(fd,"\t\t\t\t\t{ \"Sockets\":[%s], \"Cores\":[%s], \"Threads\":[%s] }\n", Sockets, Cores, Threads);
		fputs("\t\t\t\t]\n", fd);
		fputs("\t\t\t}\n", fd);
		fputs("\t\t],\n", fd);
	
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Decoders
		fputs("\t\t\"Decoders\":{\n", fd);
		fprintf(fd, "\t\t\t\"Decoded Blocks\": %u,\n", gTopology.nDecoderBlocks);
		fprintf(fd, "\t\t\t\"Dispatcher Thread\":{ \"Sockets\":[%s], \"Cores\":[%s], \"Threads\":[%s] },\n", Sockets, Cores, Threads);
		fputs("\t\t\t\"Threads\":[\n", fd);
		for (unsigned int i=1; i<gTopology.nDecoders;i++) {
			fprintf(fd,"\t\t\t\t{ \"Sockets\":[%s], \"Cores\":[%s], \"Threads\":[%s] },\n", Sockets, Cores, Threads);
		}
		fprintf(fd,"\t\t\t\t{ \"Sockets\":[%s], \"Cores\":[%s], \"Threads\":[%s] }\n", Sockets, Cores, Threads);
		fputs("\t\t\t]\n", fd);
		fputs("\t\t},\n", fd);
		
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Aligners
		fputs("\t\t\"Aligners\":{\n", fd);
		fprintf(fd, "\t\t\t\"Input Blocks\": %u,\n", gTopology.nAlignerInBlocks);
		fprintf(fd, "\t\t\t\"Output Blocks\": %u,\n", gTopology.nAlignerOutBlocks);
		fprintf(fd, "\t\t\t\"Output Block Size\": %lu,\n", gTopology.nAlignerOutBlockSize);
		fputs("\t\t\t\"Pool\":[\n", fd);
		for (unsigned int i=0; i<gTopology.nAligners;i++) {
			fprintf(fd, "\t\t\t\t{\n\t\t\t\t\t\"Master\":{ \"Sockets\":[%s], \"Cores\":[%s], \"Threads\":[%s] },\n", Sockets, Cores, Threads);
			fputs("\t\t\t\t\t\"Slaves\":[\n", fd);
			for (unsigned int j=1; j<gTopology.nAlignerThreads;j++) {
				fprintf(fd,"\t\t\t\t\t\t{ \"Sockets\":[%s], \"Cores\":[%s], \"Threads\":[%s] },\n", Sockets, Cores, Threads);
			}
			fprintf(fd,"\t\t\t\t\t\t{ \"Sockets\":[%s], \"Cores\":[%s], \"Threads\":[%s] }\n\t\t\t\t\t]\n", Sockets, Cores, Threads);

			fputs("\t\t\t\t}\n", fd);
		}
		fputs("\t\t\t]\n\t\t}\n", fd);
	}
	fputs("\t}\n", fd);
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// FOOTER
	fputs("}\n", fd);
	
	free(Sockets);
	free(Cores);
	free(Threads);
}
//---------------------------------------------------------------

static inline topology_mask_t getCPUBitMask(json_object * const restrict jobj) 
{
	topology_mask_t res = (topology_mask_t) 0;
	assert(json_object_is_type(jobj, json_type_array));
	
	const int arrayLength = json_object_array_length(jobj);
	for (int i=0; i<arrayLength; i++) {
		json_object * const restrict tmp = json_object_array_get_idx(jobj, i);
		if (json_object_is_type(tmp, json_type_int)) {
			const int val = json_object_get_int(tmp);
			assert(val >= 0);
			res |= ((topology_mask_t) 1) << val;
		}
		else {
			fputs("JSON Non numeric value encountered in topology mask!\n", stderr);
		}
	}
	return res;
}
//---------------------------------------------------------------

static void getCPUThreadID(const json_object * const restrict jobj, topology_mask_t Ids[3]) 
{
	Ids[0] = Ids[1] = Ids[2] = 0UL;
	json_object * jSockets;
	json_object_object_get_ex(jobj, "Sockets", &jSockets);
	if (jSockets == NULL) {
		fputs("Unable to find Sockets JSON object\n", stderr);
		exit(1);
	}
	Ids[0] = getCPUBitMask(jSockets);
	
	json_object * jCores;
	json_object_object_get_ex(jobj, "Cores", &jCores);
	if (jCores == NULL) {
		fputs("Unable to find Cores JSON object\n", stderr);
		exit(1);
	}
	Ids[1] = getCPUBitMask(jCores);
	
	json_object * jThreads;
	json_object_object_get_ex(jobj, "Threads", &jThreads);
	if (jThreads == NULL) {
		fputs("Unable to find Threads JSON object\n", stderr);
		exit(1);
	}
	Ids[2] = getCPUBitMask(jThreads);
	
}
//---------------------------------------------------------------

static int getCPUPool(const SystemInfo * const restrict info, json_object * const restrict jobj,
											Affinity_Mask_t * * Threads)
{
	topology_mask_t Ids[3];
	
						if (json_object_is_type(jobj, json_type_array)) {
							const int nSize = json_object_array_length(jobj);
							if (nSize > 0) {
								Affinity_Mask_t * const restrict set = (Affinity_Mask_t*) malloc(nSize*sizeof(Affinity_Mask_t));
								if (set) {
									for (int i=0; i<nSize;i++) {
										const json_object * const restrict tmp = json_object_array_get_idx(jobj, i);
										getCPUThreadID(tmp, Ids);
										getCumulativeMask(info, Ids[0], Ids[1], Ids[2], &set[i]);
										fprintf(stdout,"\tPool %i :\n", i); printfCPUThreadID(Ids);
									}
									*Threads= set;
									return nSize;
								}
							}
						}
	return 0;
}
//---------------------------------------------------------------

int loadTopology(const SystemInfo * const restrict info, const char * const restrict FileName)
{
	int err = 1;
	
	json_object * const jBase = json_object_from_file(FileName);
	if (! jBase) {
		fprintf(stderr,"Error loading JSON file %s\n", FileName);
		goto bail;
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	// Get the cpu HOST informations
	{
		topology_mask_t Ids[3];
		json_object * jHOST;
		json_object_object_get_ex(jBase, "HOST", &jHOST);
		if (jHOST == NULL) {
			fprintf(stderr, "JSON Unable to find HOST within JSON file %s\n", FileName);
			goto bail;
		}
		////////////////////////////////////////////////////////////////////////////////////////////////////
		// Get the Reader 
		json_object * jReader;
		json_object_object_get_ex(jHOST, "Reader", &jReader);
		if (jReader == NULL) {
			fprintf(stderr, "JSON Unable to find Reader within HOST in JSON file %s\n", FileName);
			goto bail;
		}
		getCPUThreadID(jReader, Ids);
		{
			json_object * jValue;
			json_object_object_get_ex(jReader, "Input Blocks", &jValue);
			if (jValue == NULL) {
				fprintf(stderr, "JSON Unable to find object \"Reader/Input Blocks\" within HOST in JSON file %s\n", FileName);
				goto bail;
			}
			if (json_object_is_type(jValue, json_type_int)) {
				const int itmp = json_object_get_int(jValue);
				gTopology.nReaderBlocks = (unsigned int) itmp;
			}
			else {
				fprintf(stderr, "JSON Error object \"Reader/Input Blocks\" is not an integer\n");
				goto bail;
			}
		}
		gTopology.Reader = (Affinity_Mask_t*) malloc(sizeof(Affinity_Mask_t));
		if (gTopology.Reader) 
			getCumulativeMask(info, Ids[0], Ids[1], Ids[2], gTopology.Reader);
		else 
			goto bail;
		
		fputs("Reader:\n", stdout); printfCPUThreadID(Ids);
		
		////////////////////////////////////////////////////////////////////////////////////////////////////
		// Get the Writers 
		json_object * jWriters;
		json_object_object_get_ex(jHOST, "Writers", &jWriters);
		if (jWriters == NULL) {
			fprintf(stderr, "JSON Unable to find Writers within HOST in JSON file %s\n", FileName);
			goto bail;
		}
		int FlushBufferCount;
		if (json_object_is_type(jWriters, json_type_array)) {
			gTopology.nWriters = json_object_array_length(jWriters);
			gTopology.Writers = (cpu_pool_t*) malloc(gTopology.nWriters*sizeof(cpu_pool_t));
			if (gTopology.Writers == NULL) {
				fprintf(stderr, "JSON Unable to allocate memory for writer pool\n");
				goto bail;
			}
			for (unsigned int iWriter=0; iWriter<gTopology.nWriters; iWriter++) {
				json_object * jCurrentWriter = json_object_array_get_idx(jWriters, iWriter);
				json_object * jWriterMaster;
				json_object_object_get_ex(jCurrentWriter, "Master", &jWriterMaster); 
				if (jWriterMaster == NULL) {
					fprintf(stderr, "JSON Unable to find \"Writers/Master\" object\n");
					goto bail;
				}
				getCPUThreadID(jWriterMaster, Ids);
				getCumulativeMask(info, Ids[0], Ids[1], Ids[2], &(gTopology.Writers[iWriter].Master));
				printf("Writer %1.1u Master:\n ", iWriter); printfCPUThreadID(Ids);
				
				json_object * jFlushBuffers;
				json_object_object_get_ex(jCurrentWriter, "Flush Buffer Threads", &jFlushBuffers); 
				if (jFlushBuffers == NULL) {
					fprintf(stderr, "JSON Unable to find \"Writers/Flush Buffer Threads\" object\n");
					goto bail;
				}
				if (!json_object_is_type(jWriters, json_type_array)) {
					fprintf(stderr,"JSON \"Writers/Flush Buffer Threads\" should be an array!\n");
					goto bail;
				}
				const int itmp = json_object_array_length(jFlushBuffers);
				if (iWriter == 0U) 
					FlushBufferCount = itmp;
				else if (itmp != FlushBufferCount) {
					fprintf(stderr, "JSON not corresponding number of flush buffer threads among writers (%i!=%i)\n", itmp, FlushBufferCount);
					goto bail;
				}
				
				gTopology.Writers[iWriter].Slaves = (Affinity_Mask_t*) malloc(FlushBufferCount*sizeof(Affinity_Mask_t));
				if (gTopology.Writers[iWriter].Slaves == NULL) {
					fputs("JSON Unable to allocate memory for Writer flush buffer threads\n", stderr);
					goto bail;
				}
				
				for (int i=0; i<FlushBufferCount; i++) {
					printf("\tWriter Flush Buffer Threads %i:\n", i);
					json_object * const restrict jobj = json_object_array_get_idx(jFlushBuffers, i);
					getCPUThreadID(jobj, Ids);
					getCumulativeMask(info, Ids[0], Ids[1], Ids[2], &(gTopology.Writers[iWriter].Slaves[i]));
					printfCPUThreadID(Ids);
				}
			}
			gTopology.nWriterThreads = FlushBufferCount;
		}
		else {
			fputs("JSON \"Writers\" object should be an array\n", stderr);
			goto bail;
		}
		
		

		////////////////////////////////////////////////////////////////////////////////////////////////////
		// Get the Decoders
		json_object * jDecoders;
		json_object_object_get_ex(jHOST, "Decoders", &jDecoders);
		if (jDecoders == NULL) {
			fprintf(stderr, "JSON Unable to find object \"Decoders\" within HOST in JSON file %s\n", FileName);
			goto bail;
		}
		gTopology.Decoders = (cpu_pool_t*) malloc(sizeof(cpu_pool_t));
		if (gTopology.Decoders == NULL) {
			fprintf(stderr, "JSON cannot allocate memory for decoders\n");
			goto bail;
		}
		{
			json_object * jInputBlocks;
			json_object_object_get_ex(jDecoders, "Decoded Blocks", &jInputBlocks);
			if (jInputBlocks) {
				if (json_object_is_type(jInputBlocks, json_type_int)) {
					const int itmp = json_object_get_int(jInputBlocks);
					if (itmp > 0)
						gTopology.nDecoderBlocks = (unsigned int) itmp;
					else {
						fprintf(stderr, "JSON Bad value %i for \"Aligners/Input Blocks\" in topology file, default to %u\n", itmp, gTopology.nDecoderBlocks);
					}
				}
			}
			else {
				fprintf(stderr, "JSON Cannot find \"Decoders/Decoded Blocks\" in topology file, default to %u\n", gTopology.nDecoderBlocks);
			}
		}
		{
			json_object * jMaster;
			json_object_object_get_ex(jDecoders, "Dispatcher Thread", &jMaster);
			if (!jMaster) {
				fprintf(stderr, "JSON Decoder does not have object \"Dispatcher Thread\"\n");
				goto bail;
			}
			getCPUThreadID(jMaster, Ids);
			getCumulativeMask(info, Ids[0], Ids[1], Ids[2], &(gTopology.Decoders->Master));
			printf("Decoder dispatcherr :\n"); printfCPUThreadID(Ids);
			json_object * jSlaves;
			json_object_object_get_ex(jDecoders, "Threads", &jSlaves);
			if (jSlaves == NULL) {
				fprintf(stderr, "JSON Unable to find object \"Decoders/Threads\"\n");
				goto bail;
			}
			if (json_object_is_type(jSlaves, json_type_array)) {
				gTopology.nDecoders = json_object_array_length(jSlaves);
				if (gTopology.nDecoders > 0) {
					gTopology.Decoders->Slaves = (Affinity_Mask_t*) malloc(gTopology.nDecoders*sizeof(Affinity_Mask_t));
					if (gTopology.Decoders->Slaves == NULL) {
						fprintf(stderr, "JSON unable to allocate memory for decoder slaves\n");
						goto bail;
					}
					for (unsigned int i=0; i<gTopology.nDecoders; i++) {
						printf("\tDecoder worker thread %i:\n", i);
						json_object * const restrict jobj = json_object_array_get_idx(jSlaves, i);
						getCPUThreadID(jobj, Ids);
						getCumulativeMask(info, Ids[0], Ids[1], Ids[2], &(gTopology.Decoders->Slaves[i]));
						printfCPUThreadID(Ids);
					}
				}
			}
		}
		
		////////////////////////////////////////////////////////////////////////////////////////////////////
		// Get the Aligners 
		json_object * jAligners;
		json_object_object_get_ex(jHOST, "Aligners", &jAligners);
		if (jAligners == NULL) {
			fprintf(stderr, "JSON Unable to find object \"Aligners\" within HOST in JSON file %s\n", FileName);
			goto bail;
		}
		
		{
			json_object * jInputBlocks;
			json_object_object_get_ex(jAligners, "Input Blocks", &jInputBlocks);
			if (jInputBlocks) {
				if (json_object_is_type(jInputBlocks, json_type_int)) {
					const int itmp = json_object_get_int(jInputBlocks);
					if (itmp > 0) 
						gTopology.nAlignerInBlocks = (unsigned int) itmp;
					else {
						fprintf(stderr, "JSON Bad value %i for \"Aligners/Input Blocks\" in topology file, default to %u\n", itmp, gTopology.nAlignerInBlocks);
					}
				}
			}
			else {
					fprintf(stderr, "JSON Cannot find \"Aligners/Input Blocks\" in topology file, default to %u\n", gTopology.nAlignerInBlocks);
			}
		}
		
		{
			json_object * jOutputBlocks;
			json_object_object_get_ex(jAligners, "Output Blocks", &jOutputBlocks);
			if (jOutputBlocks) {
				if (json_object_is_type(jOutputBlocks, json_type_int)) {
					const int itmp = json_object_get_int(jOutputBlocks);
					if (itmp > 0) 
						gTopology.nAlignerOutBlocks = (unsigned int) itmp;
					else {
						fprintf(stderr, "JSON Bad value %i for \"Aligners/Output Blocks\" in topology file, default to %u\n", itmp, gTopology.nAlignerOutBlocks);
					}
				}
			}
			else {
					fprintf(stderr, "JSON Cannot find \"Aligners/Output Blocks\" in topology file, default to %u\n", gTopology.nAlignerOutBlocks);
			}
		}
		
		{
			json_object * jOutputBlockSize;
			json_object_object_get_ex(jAligners, "Output Block Size", &jOutputBlockSize);
			if (jOutputBlockSize) {
				if (json_object_is_type(jOutputBlockSize, json_type_int)) {
					const int itmp = json_object_get_int(jOutputBlockSize);
					if (itmp > 0) 
						gTopology.nAlignerOutBlockSize = (off_t) itmp;
					else {
						fprintf(stderr, "JSON Bad value %i for \"Aligners/Output Block Size\" in topology file, default to %zu\n", itmp, gTopology.nAlignerOutBlockSize);
					}
				}
			}
			else {
					fprintf(stderr, "JSON Cannot find \"Aligners/Output Block Size\" in topology file, default to %zu\n", gTopology.nAlignerOutBlockSize);
			}
		}
		
		{
			json_object * jPool;
			json_object_object_get_ex(jAligners, "Pool", &jPool);
			if (json_object_is_type(jPool, json_type_array)) {
				const int nAligners = json_object_array_length(jPool);
				if (nAligners > 0) {
					gTopology.Aligners = (cpu_pool_t*) malloc(nAligners*sizeof(cpu_pool_t));
					if (gTopology.Aligners == NULL) goto bail;

					int PoolSize;
					for (unsigned int i=0; i<nAligners; i++) {
						json_object * restrict jobj = json_object_array_get_idx(jPool, i);
						
						json_object * jMaster;
						json_object_object_get_ex(jobj, "Master", &jMaster);
						if (!jMaster) {
							fprintf(stderr, "JSON Aligner does not have a master\n");
							goto bail;
						}
						getCPUThreadID(jMaster, Ids);
						getCumulativeMask(info, Ids[0], Ids[1], Ids[2], &(gTopology.Aligners[i].Master));
						printf("Aligner Master %1.1u :\n", i); printfCPUThreadID(Ids);
						
						json_object * jSlaves;
						json_object_object_get_ex(jobj, "Slaves", &jSlaves);
						if (!jSlaves) {
							fprintf(stderr, "JSON Aligner does not have slaves\n");
							goto bail;
						}
						
						const int itmp = getCPUPool(info, jSlaves, &(gTopology.Aligners[i].Slaves));
						if (i == 0) 
							PoolSize = itmp;
						else if (itmp != PoolSize) {
							fprintf(stderr, "JSON Aligners' pool of different sizes %i != %i\n", itmp, PoolSize);
							goto bail;
						}
					}
					gTopology.nAlignerThreads = (unsigned int) PoolSize;
				}
				else {
					fprintf(stderr, "JSON Aligners object not of good size %i \n", nAligners);
					goto bail;
				}
			}
			else {
				fputs( "JSON \"Aligners/Pool\" object is not an array\n", stderr);
				goto bail;
			}
		}
	}
	
	err = 0;
	bail:;

	return err;
}
//---------------------------------------------------------------

static void printAffinity(FILE* const out, Affinity_Mask_t * const restrict Mask, const int N)
{
    for (int i=(N-1); i>=0; --i) {
			const int digit = (int) '0' + CPU_ISSET_S(i, sizeof(Affinity_Mask_t), (cpu_set_t*) Mask);
			fputc(digit, out);
    }	
    fputc('\n', out);
}
//---------------------------------------------------------------

void printParameters(FILE* const out)
{
	fprintf(out, "-------------------------------- TOPOLOGY ------------------------------\n");
	fprintf(out, " - Number of memory blocks for reading         : %u\n", gTopology.nReaderBlocks);
	fprintf(out, "   Tags contained by reading block             : %'u\n", kMaxTagsPerBuffer);
	fprintf(out, " - Number of memory blocks for decoding        : %u\n", gTopology.nDecoderBlocks);
	fprintf(out, " - Number of memory blocks for aligner pool in : %u\n", gTopology.nAlignerInBlocks);
	fprintf(out, " - Number of memory blocks for aligner pool out: %u\n", gTopology.nAlignerOutBlocks); 
	fprintf(out, "   Size in Bytes of each aligner pool out block: %'lu\n", gTopology.nAlignerOutBlockSize);
#ifdef USE_AFFINITY
	const int N = (int) sysconf(_SC_NPROCESSORS_ONLN);
	if (gTopology.Reader) {
		fprintf(out, " - Reader thread affinity mask\n\t");
		printAffinity(out, gTopology.Reader, N);
	}
	if (gTopology.Decoders)
	{
		fprintf(out, " - Master Reader To Decoder thread affinity mask\n\t");
		printAffinity(out, &(gTopology.Decoders->Master), N);
	}
#endif
	fprintf(out, " - Number of decoder threads                   : %u\n", gTopology.nDecoders);
#ifdef USE_AFFINITY
	if (gTopology.Decoders) {
		fprintf(out, " - Decoder threads affinity masks\n");
		for (unsigned int i=0; i<gTopology.nDecoders; i++) {
			fputc('\t', out);
			printAffinity(out, &(gTopology.Decoders->Slaves[i]), N);
		}
	}
#endif
	fprintf(out, " - Number of CPU aligners                      : %u\n", gTopology.nAligners);
	fprintf(out, " - Number of CPU worker threads per aligner    : %u\n", gTopology.nAlignerThreads);
#ifdef USE_AFFINITY
	if (gTopology.Aligners) {
		for (unsigned int i=0; i<gTopology.nAligners; i++) {
			fprintf(out, " - Master %u aligner affinity mask\n\t", i);
			printAffinity(out, &(gTopology.Aligners[i].Master), N);
			fprintf(out, " - Aligner %i worker threads affinity masks\n", i);
			for (unsigned int j=0; j<gTopology.nAlignerThreads; j++) {
				fputc('\t', out);
				printAffinity(out, &(gTopology.Aligners[i].Slaves[j]), N);
			}
		}
	}
#endif
	fprintf(out, " - Number of writers                           : %u\n", gTopology.nWriters);
	fprintf(out, " - Number of flush buffer thread per writer    : %u\n", gTopology.nWriterThreads);
#ifdef USE_AFFINITY
	if (gTopology.Writers) {
		for (unsigned int i=0; i<gTopology.nWriters; i++) {
			fprintf(out, " - Master %u writer affinity mask\n\t", i);
			printAffinity(out, &(gTopology.Writers[i].Master), N);
			fprintf(out, " - Writer %i flush buffer threads affinity masks\n", i);
			for (unsigned int j=0; j<gTopology.nWriterThreads; j++) {
				fputc('\t', out);
				printAffinity(out, &(gTopology.Writers[i].Slaves[j]), N);
			}
		}
	}
#endif
	fprintf(out, "------------------------------------------------------------------------\n");
}
//---------------------------------------------------------------

/*
int main(int argc, char *argv[]) 
{
	SystemInfo info;
	getSystemInfo(&info);
	//printSystemInfo(&info);
	FILE * jsonfile = fopen("test.json", "w");
	if (jsonfile) {
		printJSON(&info, jsonfile);
		fclose(jsonfile);
	}
	
	if (loadTopology(&info, "test.json") != 0) {
		fprintf(stderr, "Error loading JSON file\n");
	}
	
	return 0;
}*/
//---------------------------------------------------------------

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
