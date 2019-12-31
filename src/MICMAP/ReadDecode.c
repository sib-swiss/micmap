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
#include <sys/statfs.h>
#include <fcntl.h>
#include <errno.h>
#include <pthread.h>
#include <string.h>
#include <assert.h>
#include <signal.h>
#include <hugetlbfs.h> 
#include <libgen.h>
#include "Decoder.h"

//--------------------------------------------------------------- 
// DEFINITIONS 
//--------------------------------------------------------------- 
#define FOUR_G (4UL * 1024UL * 1024UL * 1024UL)
#define CORRECT_SIZE (FOUR_G*kFollowNtCount)
#define MAP_ADDR NULL
//(void*) FOUR_G
#define MAP_FLAGS (MAP_SHARED)

//--------------------------------------------------------------- 
// GLOBALS 
//--------------------------------------------------------------- 
unsigned char * addr = MAP_FAILED;
int fd = -1;

void termination_handler (int signum)
{
  if (addr != MAP_FAILED) munmap(addr, kSTRAND_DATA_SIZE + CORRECT_SIZE);
        if (fd >= 0) close(fd);
  printf("\rQuit cleanly...\n");
        exit(0);
}
//---------------------------------------------------------------

int main (int argc, char *argv[])
{
		struct sigaction new_action, old_action;
		new_action.sa_handler = termination_handler;
                sigemptyset (&new_action.sa_mask);
                new_action.sa_flags = 0;

                sigaction (SIGINT, NULL, &old_action);
                if (old_action.sa_handler != SIG_IGN) sigaction (SIGINT, &new_action, NULL);
                sigaction (SIGHUP, NULL, &old_action);
                if (old_action.sa_handler != SIG_IGN) sigaction (SIGHUP, &new_action, NULL);
                sigaction (SIGTERM, NULL, &old_action);
                if (old_action.sa_handler != SIG_IGN) sigaction (SIGTERM, &new_action, NULL);

		char temp[256];
		strncpy(temp, argv[1], 256);
		char * path = dirname(temp);
		printf("File in directory %s ", path);
		if (hugetlbfs_test_path(path)) {
			
			struct statfs sb;
			int err;
			err = statfs(path, &sb);
			if (err) perror("statfs64");
			const long pgsize = sb.f_bsize;
			printf("is backed by huge TLB of page size %li\n", pgsize);
		}
		else {
			printf("is not  backed by huge TLB\n");
		}

		fd = open(argv[1], O_RDONLY, S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH);
		if (fd < 0) {
			perror("shm_open");
			return 1;
		}
		
		const size_t TotalSize = kSTRAND_DATA_SIZE + CORRECT_SIZE;
	
		addr = mmap(MAP_ADDR, TotalSize, PROT_READ, MAP_FLAGS|MAP_HUGETLB|MAP_HUGE_1GB, fd, 0);
		if (addr == MAP_FAILED) {
			fprintf(stderr, "Unable to get sufficient memory to map %s in memory 1: %i\n", argv[1], errno);
			perror("mmap");
		}
		else goto Next;
		addr = mmap(MAP_ADDR, TotalSize, PROT_READ, MAP_FLAGS|MAP_HUGETLB|MAP_HUGE_2MB, fd, 0);
                if (addr == MAP_FAILED) {
                        fprintf(stderr, "Unable to get sufficient memory to map %s in memory 2: %i\n", argv[1], errno);
                        perror("mmap");
                }
		else goto Next;
		addr = mmap(MAP_ADDR, TotalSize, PROT_READ, MAP_FLAGS|MAP_HUGETLB, fd, 0);
                if (addr == MAP_FAILED) {
                        fprintf(stderr, "Unable to get sufficient memory to map %s in memory 3: %i\n", argv[1], errno);
                        perror("mmap");
                }
		else goto Next;
		addr = mmap(MAP_ADDR, TotalSize, PROT_READ, MAP_FLAGS, fd, 0);
                if (addr == MAP_FAILED) {
                        fprintf(stderr, "Unable to get sufficient memory to map %s in memory 4: %i\n", argv[1], errno);
                        perror("mmap");
                }
	Next:
		if (addr != MAP_FAILED) {	
			for (int i=0; i<256; i++) { 
				fprintf(stdout, "%2x ", (int) addr[CORRECT_SIZE+i]);
				if (((i+1) & 0xF) == 0) fputc('\n', stdout);
			}
			while (1) { sleep(60); }
		}
		close(fd);
		return 0;
		
}
//---------------------------------------------------------------

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
