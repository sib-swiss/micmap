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
#define _GNU_SOURCE
#include "constants.h"
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <pthread.h>
#include <string.h>
#include <assert.h>
#include <signal.h>
#include "Decoder.h"

//--------------------------------------------------------------- 
// DEFINITIONS 
//--------------------------------------------------------------- 
#define FOUR_G (4UL * 1024UL * 1024UL * 1024UL)
#define CORRECT_SIZE (FOUR_G*kFollowNtCount)
#define MAP_ADDR  NULL
//(void*) FOUR_G
#define MAP_FLAGS (MAP_SHARED)
#undef MAP_HUGETLB
#define MAP_HUGETLB 0
#undef MAP_HUGE_1GB
#undef MAP_HUGE_2MB
#define MAP_HUGE_1GB 0
#define MAP_HUGE_2MB 0

//--------------------------------------------------------------- 
// STRUCTURE DEFINITIONS 
//---------------------------------------------------------------
typedef struct _thr_args {
	char * bufPtr;
	size_t size;
	unsigned int chunk;
	const char * fn;
} thr_arg_t, *thr_arg_p_t;

//--------------------------------------------------------------- 
// GLOBALS 
//--------------------------------------------------------------- 
void* addr = MAP_FAILED;
int fd = -1;
char Strand[512];
char ValidTable[512];
const char * MapName = NULL;

//--------------------------------------------------------------- 
// LOCAL FCTS 
//--------------------------------------------------------------- 
static void * read_buffer(void *arg)
{
	thr_arg_p_t tap = (thr_arg_p_t) arg;
	int fd = open(tap->fn, O_RDONLY);
	size_t total = 0;
	if (lseek(fd, tap->size * tap->chunk, SEEK_SET) == -1) {
		fprintf(stderr, "Could not lseek: %s(%d)\n", strerror(errno), errno);
		return (void *) 1;
	}
	while (total < tap->size) {
		ssize_t rc;
		while ((rc = read(fd, tap->bufPtr, 32768)) == -1) {
			if (errno != EINTR && errno != EAGAIN) {
				fprintf(stderr, "Could not read: %s(%d)\n", strerror(errno), errno);
				return (void *) 1;
			}
		}
		total += rc;
		tap->bufPtr += rc;
	}
	close(fd);
	return 0;
} /* read_4G_buffer */
//---------------------------------------------------------------

void termination_handler (int signum)
{
	if (addr != MAP_FAILED) munmap(addr, kSTRAND_DATA_SIZE+CORRECT_SIZE);	
	if (fd >= 0) close(fd);		
	shm_unlink(MapName);
	exit(0);
}
//---------------------------------------------------------------

int main (int argc, char *argv[])
{
		char Strand[512];
		char ValidTable[512];
		struct sigaction new_action, old_action;
		char NT = 'A';
		
		if (argc != 3) {
			fprintf(stderr, "%s shared memory name genome decoding file to load in RAM\n", argv[0]);
			goto bail;
		}
		
		new_action.sa_handler = termination_handler;
		sigemptyset (&new_action.sa_mask);
		new_action.sa_flags = 0;
		
		sigaction (SIGINT, NULL, &old_action);
		if (old_action.sa_handler != SIG_IGN) sigaction (SIGINT, &new_action, NULL);
		sigaction (SIGHUP, NULL, &old_action);
		if (old_action.sa_handler != SIG_IGN) sigaction (SIGHUP, &new_action, NULL);
		sigaction (SIGTERM, NULL, &old_action);
		if (old_action.sa_handler != SIG_IGN) sigaction (SIGTERM, &new_action, NULL);
		
		MapName = argv[1];
		
		fd = open(argv[1], O_EXCL | O_CREAT | O_RDWR, S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH);
		if (fd < 0) {
			perror("shm_open");
			goto bail;
		}
		
		/*
		fd = mkstemp64(MapName);

		if (fd < 0) {
			fprintf(stderr, "mkstemp() failed: %s\n", strerror(errno));
			return -1;
		}

		unlink(MapName);
		*/
		struct stat validStat;
		snprintf(ValidTable, 512, "%s%c.chr.bin",argv[2],NT);
		if (stat (ValidTable,&validStat) != 0) {
			printf("cannot stat file %s\n", ValidTable);
			goto bail;
		}
		if (validStat.st_size != CORRECT_SIZE) {
			printf("validStat.st_size %lx != %lx\n", validStat.st_size, CORRECT_SIZE);
			goto bail;
		}
		snprintf(Strand,512,"%s%c.strand.bin",argv[2],NT);
		if (stat (Strand,&validStat) != 0) {
			printf("cannot stat file %s\n", Strand);
			goto bail;
		}
		if (validStat.st_size != kSTRAND_DATA_SIZE) {
			printf("validStat.st_size %lx != %lx\n", validStat.st_size, kSTRAND_DATA_SIZE);
			goto bail;
		}	
	
		const size_t TotalSize = kSTRAND_DATA_SIZE + CORRECT_SIZE;
	
		ftruncate(fd, TotalSize);
		addr = mmap(MAP_ADDR, TotalSize, 
		            PROT_READ | PROT_WRITE,
		            MAP_FLAGS | MAP_HUGETLB | MAP_HUGE_1GB, fd, 0);
		if (addr == MAP_FAILED) {
			fprintf(stderr, "Unable to get sufficient memory to map %s in Huge TLB 1GB memory : err=%i\n", argv[2], errno);
			perror("mmap");
			addr = mmap(MAP_ADDR, TotalSize,
			            PROT_READ | PROT_WRITE,
			            MAP_FLAGS | MAP_HUGETLB | MAP_HUGE_2MB, fd, 0);
			if (addr == MAP_FAILED) {
				fprintf(stderr, "Unable to get sufficient memory to map %s in Huge TLB 2MB memory\n", argv[2]);
				perror("mmap");
				goto bail;
			}
		}
		
		pthread_t t[kFollowNtCount+1];
		thr_arg_t ta[kFollowNtCount+1];
		char * restrict ptr = addr;
		for (int i = 0; i < kFollowNtCount; i++) {
			ta[i].bufPtr = ptr;
			ta[i].chunk = i;
			ta[i].size = FOUR_G;
			ta[i].fn = ValidTable;
			if (pthread_create(t + i, NULL, read_buffer, ta + i) != 0) {
				fprintf(stderr, "Could not create thread %d: %s(%d)\n", i, strerror(errno), errno);
				goto bail;
			}
			ptr += FOUR_G;
		}
		ta[kFollowNtCount].bufPtr = ptr;
		ta[kFollowNtCount].chunk = 0;
		ta[kFollowNtCount].size = kSTRAND_DATA_SIZE;
		ta[kFollowNtCount].fn = Strand;
		if (pthread_create(t + kFollowNtCount, NULL, read_buffer, ta + kFollowNtCount) != 0) {
			fprintf(stderr, "Could not create thread %d: %s(%d)\n", kFollowNtCount, strerror(errno), errno);
			goto bail;
		}
		
		for (int i = 0; i < (kFollowNtCount+1); i++) {
			if (pthread_join(t[i], NULL) != 0) {
				fprintf(stderr, "Could not join thread %d: %s(%d)\n", i, strerror(errno), errno);
				return 1;
			}
		}
		fputs("Finish loading decoding table...\n", stderr);
		if (msync(addr, TotalSize, MS_SYNC) != 0) {
			perror("msync");
			goto bail;
		} 		
 		for (int i=0; i<256; i++) {
			fprintf(stdout, "%2x ", (int) ((unsigned char*)addr)[CORRECT_SIZE+i]);
			if (((i+1) & 0xF) == 0) fputc('\n', stdout);
                }
		while(1) { sleep(60); }
				
		return 0;
		
	bail:;
	  if (addr != MAP_FAILED) munmap(addr, TotalSize);	
	  if (fd >= 0) close(fd);
	  shm_unlink(argv[1]); 		
	  return 1;
}
//---------------------------------------------------------------

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
