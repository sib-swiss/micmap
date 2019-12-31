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
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <libgen.h> 
#include <zlib.h> 

int main(int argc, char *argv[])
{
	int res = 1;
	if (argc != 3) {
		printf("Usage: %s executable /full/path/and/header_name\n", argv[0]);
		return 1;
	}
	
	const int fd = open(argv[1], O_RDONLY);
	if (fd < 0 ) {
		fprintf(stderr, "Cannot acces file %s\n", argv[1]);
		return 1;
	}
	
	struct stat st;
	const int err = fstat(fd, &st);
	if (res < 0) {
		fprintf(stderr, "Cannot stat file %s\n", argv[1]);
		goto bail;
	}
	
	char header[255];
	char blobfile[255];
	
	snprintf(header, 255, "%s.h", argv[2]);
	FILE * outHeader = fopen(header, "w");
	if (outHeader == NULL) {
		fprintf(stderr, "Cannot create header file %s\n", header);
		goto bail;
	}
	snprintf(blobfile, 255, "%s.bin", argv[2]);
	FILE * outBlob = fopen(blobfile, "wb");
	if (outBlob == NULL) {
		fprintf(stderr, "Cannot create header file %s\n", blobfile);
		goto clean3;
	}
	
	unsigned char * const restrict Code = (unsigned char*) malloc((1+st.st_size)*sizeof(char));
	if (Code == NULL) {
		fprintf(stderr, "Cannot allocate memory to read %s\n", argv[1]);
		goto clean2;
	}
	
	unsigned char * const restrict compressed = (unsigned char*) malloc((2*st.st_size)*sizeof(char));
	if (compressed == NULL) {
		fputs("Cannot allocate memory for compressed code\n", stderr);
		goto clean1;
	}
	
	const ssize_t bytes = read(fd, (void*) Code, st.st_size);
	if (bytes != (ssize_t) st.st_size) {
		fprintf(stderr, "Read error in %s\n", argv[1]);
		goto clean;
	}
	
	size_t psize = 2*st.st_size;
	const int rs = compress2(compressed, &psize, Code, st.st_size, 9);
	if (rs != Z_OK || psize >= st.st_size) {
		fprintf(stderr, "Error in compressing %s\n", argv[1]);
		goto clean;
	}
	
	const size_t compressedSize = psize;
	
	const size_t wbytes = fwrite(compressed, sizeof(char), (size_t) compressedSize, outBlob);
	if (wbytes != (size_t) compressedSize) {
		fprintf(stderr, "Error writnig compressed file %s\n", blobfile);
		goto clean;
	}
	
	char * buffer = strdup(argv[2]);
	char * const SYMBOLName = basename(buffer);
	{
		unsigned char * ptr = (unsigned char*) SYMBOLName;
		while (*ptr != '\0') {
			if (*ptr >= 'a' && *ptr <= 'z') *ptr = *ptr - 'a' + 'A';
			ptr++;
		}
	}
	
	fprintf(outHeader, "#ifndef %s_H_\n#define %s_H_\n\n", SYMBOLName, SYMBOLName);
	fprintf(outHeader, "#define %s_ORIGINAL_SIZE %lu\n#define %s_COMPRESSED_SIZE %lu\n\n",
					SYMBOLName, st.st_size, SYMBOLName, compressedSize);
	fputs("#endif\n\n", outHeader);
	
	free(buffer);
	
	{
		memset(Code, 0, st.st_size*sizeof(char));
		psize = st.st_size;
		const int rs2 = uncompress(Code, &psize, compressed, compressedSize);
		if (rs2 != Z_OK) {
			fprintf(stderr, "Error in redecompressing data...\n");
		}
		snprintf(header, 255, "%s.orig", argv[2]);
		FILE* tst = fopen(header, "wb");
		fwrite(Code, sizeof(char),psize, tst);
		fclose(tst);
	}
	
	res = 0;
	
	clean:
		free(compressed);
	
	clean1:
		free(Code);
	
	clean2:
		fclose(outBlob);
		
	clean3:
		fclose(outHeader);
		
	bail:
		close(fd);
		return res;
}
