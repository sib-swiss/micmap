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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <sys/stat.h>
#include <sys/statfs.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include <sys/mman.h>
#include <limits.h>
#include <libgen.h>
#include "config.h"
#include "Genome.h"

#define GENOME_SIZE (4UL*1024*1024*1024)
#if GENOME_SIZE != kGENOME_DATA_SIZE
#error "Various size for genome"
#endif

void *LoadGenome(Genome_t * const restrict genome)
{
	struct statfs64 stfs;
	struct stat64 st;
	char temp[512];
	int fd = -1;

	genome->configBuffer = NULL;
	genome->configBufferSize = 0UL;
	genome->table = NULL;
	genome->tableSize = 0UL;
	genome->virtchr = NULL;
	genome->isHugeTLBFS = false;

	/* ----------------------- load genome data ----------------------*/
	const int gfd = open64 (genome->FileName, O_RDONLY);
	if (!gfd) {
		fprintf(stderr, "%s: can't open bin file\n", __FUNCTION__);
		goto bail;
	}

	if (fstat64 (gfd,&st) != 0) {
		fprintf(stderr, "%s: cannot stat file %s\n", __FUNCTION__, genome->FileName);
		goto bail;
	}
	if (st.st_size != kGENOME_DATA_SIZE) {
		fprintf(stderr, "%s: genomeStat.st_size %lx != %lx\n", __FUNCTION__, st.st_size,kGENOME_DATA_SIZE);
		goto bail;
	}

	strncpy(temp, genome->FileName, 512);
	const char * const path = dirname(temp);
	if (statfs64(path, &stfs)) {
		perror("statfs in loadDecodingTable");
		goto bail;
	}

	const int flags = (stfs.f_type == HUGETLBFS_MAGIC) ? MAP_SHARED : MAP_PRIVATE;

	genome->table = (unsigned char *)mmap64 (genome->table, kGENOME_DATA_SIZE, PROT_READ, flags, gfd, 0);
	if (genome->table == MAP_FAILED) {
		fprintf(stderr, "%s: mapping failed with mmap64 error: %s.\n", __FUNCTION__, strerror(errno));
		goto bail;
	}

	close(gfd);

	/* ----------------------- load genome config ----------------------*/
	// Use same name as genome file, but ending in .cfg instead of .bin
	if (stfs.f_type != HUGETLBFS_MAGIC) {
		strncpy(temp, genome->FileName, 512);
		unsigned int len = strlen(temp);
		memcpy(temp + len - 3,"cfg",3);

		fd = open(temp, O_RDONLY);
		if (fd == -1) { perror("genome open"); goto bail; }
		if (fstat64(fd, &st) != 0) { perror("genome fstat"); goto bail; }
		char *configBuf = malloc(st.st_size);
		if (configBuf == NULL) { perror("malloc"); goto bail; }
		size_t rCnt = st.st_size;
		char *rPos = configBuf;
		ssize_t res;

		while ((res = read(fd,rPos,rCnt)) != rCnt) {
			if (res == 0) {
				fprintf(stderr,"we miss %ld bits (got 0)\n", rCnt);
				goto bail;
			}
			if (res == -1) {
				if (errno != EINTR) {
					perror("read");
					goto bail;
				}
			}
			else {
				rPos += res;
				rCnt -= res;
			}
		}
		close(fd);

		genome->configBuffer = configBuf;
		genome->configBufferSize = st.st_size;
		genome->tableSize = kGENOME_DATA_SIZE;
	}
	else {
		const size_t offset = *((size_t*) &(genome->table[GENOME_SIZE-sizeof(size_t)]));
		genome->configBufferSize = offset;
		genome->configBuffer = (char *) &(genome->table[GENOME_SIZE-sizeof(size_t)-offset]);
		genome->tableSize = GENOME_SIZE;
		genome->isHugeTLBFS = true;
	}

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

	return (void*) NULL;

bail:;
	FreeGenome(genome);
	if (fd != -1) close(fd);

	return (void*) (uintptr_t) 1;
}
//---------------------------------------------------------------

int GenomeToRAM(const char * const restrict path, const char * const restrict genome)
{
	unsigned char * restrict addr = MAP_FAILED;
	struct statfs64 stfs;
	Genome_t lGenome;
	char fn[512];
	char basefn[512];
	int fd = -1;


	if(statfs64(path, &stfs)) {
		perror("statfs64");
		goto bail;
	}

	if (stfs.f_type == HUGETLBFS_MAGIC) {
		strncpy(basefn, genome, 512);
		snprintf(fn, 512, "%s/%s", path, basename(basefn));
		fd = open(fn, O_EXCL | O_CREAT | O_RDWR, S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH);
		if (fd < 0) {
			perror("open");
			goto bail;
		}
		if (ftruncate64(fd, GENOME_SIZE)) {
			fprintf(stderr, "Unable to extend file %s to %lu\n", fn, GENOME_SIZE);
			goto bail;
		}
		addr = mmap64(NULL, GENOME_SIZE, PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0);
		if (addr == MAP_FAILED) { perror("mmap64"); goto bail; }

		lGenome.FileName = genome;

		void * res = LoadGenome(&lGenome);
		if (res != NULL) goto bail;

		if (lGenome.tableSize > GENOME_SIZE) {
			fprintf(stderr, "Genome file %s is bigger than 4GB\n", genome);
			goto bail;
		}
		memcpy(addr, lGenome.table, lGenome.tableSize);
		*((size_t*) &addr[GENOME_SIZE-sizeof(size_t)]) = lGenome.configBufferSize;
		memcpy(&addr[GENOME_SIZE-sizeof(size_t)-lGenome.configBufferSize], lGenome.configBuffer, lGenome.configBufferSize);

		FreeGenome(&lGenome);

		munmap(addr, GENOME_SIZE);
		close(fd);
	}
	else {
		fprintf(stderr, "%s does not point to a huge TLB file system\n", path);
		goto bail;
	}
	return 0;

	bail:
		FreeGenome(&lGenome);
		if (fd > 0) close(fd);
		if (addr != MAP_FAILED) munmap(addr, GENOME_SIZE);
		if (fd > 0) unlink(fn);
		return 1;
}
//---------------------------------------------------------------

void FreeGenome(Genome_t * const genome)
{
	if (!genome->isHugeTLBFS && genome->configBuffer) free(genome->configBuffer);
	if (genome->table != MAP_FAILED) munmap(genome->table, genome->tableSize);
	if (genome->virtchr) free(genome->virtchr);
	freeVirtualChromosome();
}
//---------------------------------------------------------------
/* vim: tabstop=2 shiftwidth=2
 */
