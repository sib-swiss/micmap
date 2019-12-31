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
#include <stdint.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>
#include <string.h>
#include <smmintrin.h>
#include <sys/time.h>
#include <alloca.h>
#include <getopt.h>
#include <stdbool.h>
#include <stddef.h>
#include <sched.h>
#include <assert.h>
#include "fastalign.h"

#ifdef __MIC__
#define CLOCK_UNIT 10
extern inline unsigned long __attribute__((always_inline)) timer(void)
{
	
	unsigned int hi, lo;
	
	__asm __volatile__ (
		"xorl %%eax, %%eax \n"
		"cpuid             \n"
		"rdtsc             \n"
		:"=a"(lo), "=d"(hi)
		:
		      :"%ebx", "%ecx"
	);
	return ((unsigned long)hi << (32-CLOCK_UNIT)) | ((unsigned long) lo >> CLOCK_UNIT) ;
}
//---------------------------------------------------------------
#endif

#ifdef USE_PDF
#include "functions.h"
struct IO_PDF Options = { 
	.WithPaths = true,
	.WithProfileScoreGraph = false,
	.WithSequenceScoreGraph = false,
	.WithReverse = false,
	.GivenRegion = false,
	.x = { 0,0 },
	.y = { 0,0 }
};
int PlotAllSequence = 0;
#endif

void PrintParameters()
{
	printf("State priority alignment: Match (%i) Insertion (%i) Deletion (%i) Extra (%i)\n",
				 PRIORITY_MATCH, PRIORITY_DELETION, PRIORITY_INSERTION, PRIORITY_EXTRA);
	printf("Costs: Match (%i) Mismatch (%i) Insertion (%i) Deletion (%i)\n",
				 ToScore(prf.M),ToScore(prf.m), ToScore(prf.D), ToScore(prf.I));
	fputs("Transition costs:\n"
	"                /     M     I     D     X\n", stdout);
	fputs("                M", stdout);
	int value = ToScore(prf.Match.To[MATCH]);
	if (value <= NLOW)
		fputs(" NLOW", stdout);
	else
		printf("%6i", value);
	value = ToScore(prf.Match.To[DELETION]);
	if (value <= NLOW)
		fputs("  NLOW", stdout);
	else
		printf("%6i", value);
	value = ToScore(prf.Match.To[INSERTION]);
	if (value <= NLOW)
		fputs("  NLOW", stdout);
	else
		printf("%6i", value);
	value = ToScore(prf.Match.To[EXTRA]);
	if (value <= NLOW)
		fputs("  NLOW\n", stdout);
	else
		printf("%6i\n", value);
	
	fputs("                I", stdout);
	value = ToScore(prf.Deletion.To[MATCH]);
	if (value <= NLOW)
		fputs("  NLOW", stdout);
	else
		printf("%6i", value);
	value = ToScore(prf.Deletion.To[DELETION]);
	if (value <= NLOW)
		fputs("  NLOW", stdout);
	else
		printf("%6i", value);
	value = ToScore(prf.Deletion.To[INSERTION]);
	if (value <= NLOW)
		fputs("  NLOW", stdout);
	else
		printf("%6i", value);
	value = ToScore(prf.Deletion.To[EXTRA]);
	if (value <= NLOW)
		fputs("  NLOW\n", stdout);
	else
		printf("%6i\n", value);
	
	fputs("                D", stdout);
	value = ToScore(prf.Insertion.To[MATCH]);
	if (value <= NLOW)
		fputs("  NLOW", stdout);
	else
		printf("%6i", value);
	value = ToScore(prf.Insertion.To[DELETION]);
	if (value <= NLOW)
		fputs("  NLOW", stdout);
	else
		printf("%6i", value);
	value = ToScore(prf.Insertion.To[INSERTION]);
	if (value <= NLOW)
		fputs("  NLOW", stdout);
	else
		printf("%6i", value);
	value = ToScore(prf.Insertion.To[EXTRA]);
	if (value <= NLOW)
		fputs("  NLOW\n", stdout);
	else
		printf("%6i\n", value);
	
	fputs("                X", stdout);
	value = ToScore(RowInit[MATCH]);
	if (value <= NLOW)
		fputs("  NLOW", stdout);
	else
		printf("%6i", value);
	value = ToScore(RowInit[DELETION]);
	if (value <= NLOW)
		fputs("  NLOW", stdout);
	else
		printf("%6i", value);
	value = ToScore(RowInit[INSERTION]);
	if (value <= NLOW)
		fputs("  NLOW", stdout);
	else
		printf("%6i", value);
	value = ToScore(RowInit[EXTRA]);
	if (value <= NLOW)
		fputs("  NLOW\n", stdout);
	else
		printf("%6i\n", value);
	
}
//---------------------------------------------------------------

static const char opt_to_test[]= "hHp:r:S:RL:D:"
#ifdef __MIC__
	"st"
#endif
;
static const struct option long_options[] =
{
	{"reverse", no_argument, 				0, 'R'},
	{"softclip", required_argument,	0, 'S'},
	{"prf",			required_argument,	0, 'L'},
	{"dump-prf", required_argument,	0, 'D'},
	#ifdef USE_PDF
	{"pdf", 		required_argument, 	0, 'p'},
	{"region",	required_argument,	0, 'r'},
	#endif
	{"header", 	no_argument, 				0, 'H'},
	{"help", 		no_argument, 				0, 'h'},
#ifdef __MIC__
	{"systolic", no_argument,				0, 's'},
	{"systolic_simplified", no_argument, 0, 't'},
#endif
	{0, 0, 0, 0}
};

static void __attribute__((noreturn)) Usage(FILE * stream)
{
	fputs(
		"FastAlignTest:\n"
		" FastAlignTest [options] tag genome\n\n"
		" Options:\n"
		"   --prf filename [-L]   : load JSON profile parameters\n"
		"   --dump-prf filename   : dump JSON profile parameters\n"
		"   --reverse      [-R]   : reverse the tag\n"
		"   --softclip n   [-S]   : allow a softclip starting from position n\n"
#ifdef __MIC__
		"   --systolic            : compute matrix with systolic version\n"
		"   --systolic-simplified : compute matrix with simplified systolic\n"
#endif
#ifdef USE_PDF
		"   --pdf filename        : generate pdf\n"
		"     --region [a:b][c:d] : specify region to output\n"
#endif
		"   --header              : print parameters value\n"
		"   --help         [-h]   : output command help\n",
		stream);
	exit(1);
}
//---------------------------------------------------------------
	
int main(int argc, char * argv[])
{
	struct timeval t0, t1, t2, t3, t4, t5;
	GlobalData_t data;
#ifndef __MIC__
	ComputeFunctions_t * Fcts = &cpu_std;
#else
	ComputeFunctions_t * Fcts = &mic_systolic;
#endif
	const char * restrict tag = NULL;
	int err = 0;
	char AlgoStateLetter[4];
	
#ifdef USE_PDF
	const char * pdfName = NULL;
#endif
	
	memset(&data, 0, sizeof(GlobalData_t));
	
	/* Setting letter according to tag, not genome I <-> D*/
	AlgoStateLetter[MATCH]     = 'M';
	AlgoStateLetter[INSERTION] = 'D';
	AlgoStateLetter[DELETION]  = 'I';
	AlgoStateLetter[EXTRA]     = 'X';
	
	while (1) {
		/* getopt_long stores the option index here. */
		int option_index = 0;
		
		const int c = getopt_long (argc, argv, opt_to_test, long_options, &option_index);
		
		/* Detect the end of the options. */
		if (c == -1) break;
		switch (c) {
			case 'h': Usage(stdout); break;
			case 'H': PrintParameters(); exit(0); break;
#ifdef USE_PDF
			case 'p': pdfName = optarg; break;
			case 'r':
			{
				int lx[2], ly[2];
				if (sscanf(optarg, "[%i:%i][%i:%i]", &lx[0], &lx[1], &ly[0], &ly[1]) != 4) {
					fprintf(stderr,"Unable to read region interval from %s\n", optarg);
					exit(1);
				}
				Options.GivenRegion = true;
				Options.x[0] = lx[0];
				Options.x[1] = lx[1];
				Options.y[0] = ly[0];
				Options.y[1] = ly[1];
				break;
			}
#endif
#ifdef __MIC__
			case 's':
				Fcts = &mic_systolic;
				break;
			case 't':
				Fcts = &mic_systolic_simplified;
				break;
#endif
			case 'S':
				data.SoftClipBoundary = atoi(optarg);
				break;
			case 'R':
				data.revNmer = 1;
				break;
			case 'L':
				if (loadAlignerScores(optarg) < 0) {
					fputs("Aligner configuration file failed to load\n", stderr);
					return 1;
				}
				break;
			case 'D':
				dumpAlignerScores(optarg);
				return 0;
			default: Usage(stdout);
		}
	}
		
	if (optind == argc) {
		fputs("Error in given options\n", stderr);
		Usage(stderr);
	} else {
		tag = argv[optind];
		data.Genome = (unsigned char *) argv[optind+1];
	}
	
#ifndef __MIC__
	/* Check for Border version depending on the prf values */
	if (prf.AllowedJumpIn != 0U || prf.AllowedJumpOut != 0U) {
		fputs("Detected Border allowed entrance or exit, adjusting to border functions\n", stderr);
		Fcts = &cpu_std_border;
	}
#endif
	
  /* Set thread affinity not to move */
	cpu_set_t affinity_mask;
	CPU_ZERO_S(sizeof(cpu_set_t), &affinity_mask);
	CPU_SET_S(4, sizeof(cpu_set_t), &affinity_mask);
	err = sched_setaffinity(0, sizeof(cpu_set_t), &affinity_mask);
	if (err == -1) {
		perror("Affinity");
	}

  /* Get Reads' length */
	data.TagLength = strlen(tag);
	
	/* Copy tag to correctly aligned memory */
	{
		const size_t kMaxTagSize = (kTagSize + 15UL) & ~(15UL);
		assert(data.TagLength <= kTagSize);
		data.Tag = _mm_malloc((kMaxTagSize)*sizeof(char), ALIGNMENT);
		if (data.Tag == NULL) {
			fputs("Unable to allocate memory for tag in align spot.\n", stderr);
			exit(1);
		}

		memcpy(data.Tag, tag, (data.TagLength+1)*sizeof(char));
	}
	
	/* Get sequence length */
	data.GenomeLength = strlen((char *) data.Genome);
	fprintf(stderr,"Reads of length %zu detected at beginning\n"
	               "Sequence of length %zu detected at beginning\n", data.TagLength, data.GenomeLength); 
#ifdef __MIC__
	/* Copy genome to correctly aligned memory */
	{
		unsigned char * const cptr = data.Genome;
		const size_t kMaxGenomeSize = kGAPPED_ALIGN_GENOME_LENGTH + 16 + 16;
		assert(data.GenomeLength < kGAPPED_ALIGN_GENOME_LENGTH);
		data.Genome = _mm_malloc(kMaxGenomeSize*sizeof(char), ALIGNMENT);
		if (data.Genome == NULL) {
			fputs("Unable to allocate memory for genome in align spot.\n", stderr);
			exit(1);
		}
		memcpy(data.Genome, cptr, (data.GenomeLength+1)*sizeof(char));
	}
#endif
	
	/* Allocating memory */
	if (Fcts->allocStorage(&data) != 0) {
		fprintf(stderr, "Error in storage allocation\n");
	}
	
  gettimeofday(&t0,0);
	const int retval = Fcts->createMatrix(&data);
	gettimeofday(&t1,0);
	
	gettimeofday(&t2,0);
	const int StateVal = Fcts->getStateSequence(&data);
	gettimeofday(&t3,0);
	
	if (StateVal == 0) {
		printf("Alignment found from %i to %i , profile [%i:%i], having %u mismatches, scoring %i.\n"
		       "Cigar size of %u required, mismatch size of %u required.\n"
					 "SoftClip data: boundary=%u, evicted_CIG_size=%u, evicted_MM_size=%u\n",
	       data.AlignmentRange[0], data.AlignmentRange[1], data.ProfileRange[0], data.ProfileRange[1], 
	       data.MismatchCount, data.score,
	       data.RequiredCigSize, data.RequiredMMSize, data.SoftClipBoundary, data.SoftClipEvictedCigSize,
				 data.SoftClipEvictedMMSize);
		fputs("STA: ", stdout);
		for ( int i=0; i<data.ProfileRange[0]; i++) fputc(' ', stdout);
		printf("%s\n", data.States);
	}
	else {
		fprintf(stderr, "Error in getStateSequence\n");
		exit(1);
	}
	
	printf("GEN: ");
	for ( int i=0; i<data.ProfileRange[0]; i++) fputc('-', stdout);
	
	const size_t States = strlen((char *) data.States);
	const char * restrict genptr = (char *) &data.Genome[data.AlignmentRange[0]];
	for (size_t i=0; i<States; i++) {
		if (data.States[i] == 'I')
			fputc('-',stdout);
		else
			fputc(*genptr++, stdout);
	}
	fputc('\n', stdout);
	
	printf("TAG: ");
	const char * restrict tagptr = (char *) &data.Tag[0];
	for ( int i=0; i<data.ProfileRange[0]; i++) fputc((*tagptr++ - 'A' + 'a'),stdout);
	for (size_t i=0; i<States; i++) {
		if (data.States[i] == 'D')
			fputc('-',stdout);
		else if (data.States[i] == 'I')
			fputc((*tagptr++ - 'A' + 'a'),stdout);
		else
			fputc(*tagptr++, stdout);
	}
	const uintptr_t limit = (uintptr_t) &data.Tag[data.TagLength];
	while((uintptr_t) tagptr < limit) fputc((*tagptr++ - 'A' + 'a'),stdout);
	fputc('\n', stdout);
	
	gettimeofday(&t4,0);
	const int retCompute = Fcts->computeCigarAndMismatch(&data);
	gettimeofday(&t5,0);
	
	if ( retCompute == 0) {
		printf("CIG: "); printCigar(&data);
		printf("MM : "); printMismatch(&data);
	}
	else {
		fprintf(stderr, "Error in computeCigarAndMismatch\n");
	}
	
#ifdef USE_PDF
	if (pdfName) {
		PDFOutput(&data, pdfName, (void*) &Options);
	}
#endif
	
  const double mat = 1.e6 * (double) (t1.tv_sec - t0.tv_sec) + (double) (t1.tv_usec - t0.tv_usec);
	const double Ali = 1.e6 * (double) (t3.tv_sec - t2.tv_sec) + (double) (t3.tv_usec - t2.tv_usec);
	const double Str = 1.e6 * (double) (t5.tv_sec - t4.tv_sec) + (double) (t5.tv_usec - t4.tv_usec);
	
	printf("Time (us): Matrix (%lf) Best alignment (%lf) Compute Cigar and Mismatch strings (%lf)\nThroughput: %lf alignment/s\n",
				 mat, Ali,Str, 1.e6/(mat+Ali+Str));
	_mm_free(data.Tag);
#ifdef __MIC__
	_mm_free(data.Genome);
#endif
	Fcts->freeStorage(&data);
	err = 0;
  return err;
}
//---------------------------------------------------------------

/* vim: tabstop=2 shiftwidth=2
 */
