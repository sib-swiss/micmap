/*
 * ------------------------------------------------------------------------------------------------------------------------
 *
 *                                * micmap *
 *      Mapping of short reads in fastq format onto a reference
 *
 *
 *  Copyright (C) SIB  - Swiss Institute of Bioinformatics,  2015-2019 Nicolas Guex, Thierry Schuepbach and Christian Iseli
 *  Copyright (C) UNIL - University of Lausanne, Switzerland 2019-2020 Nicolas Guex and Christian Iseli
 *  Copyright (C) EPFL - Lausanne, Switzerland               2020-2021 Christian Iseli
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
 *      Contacts:   Nicolas.Guex@unil.ch and christian.iseli@epfl.ch
 *      Repository: https://github.com/sib-swiss/micmap
 *
 * ------------------------------------------------------------------------------------------------------------------------
 */
#include "config.h"
#include "constants.h"
#define _GNU_SOURCE
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <sys/time.h>
#include <pthread.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <pthread.h>
#include <string.h>
#include <dlfcn.h>
#include <sys/ioctl.h>
#include <execinfo.h>
#include <signal.h>

#ifdef USE_MPI
#include <mpi.h>
#include "MPI/mpi_transfer.h"
#endif

#include "GTL.h"
#include "gtlVersion.h"
#include "Reader.h"
#include "Writer.h"
#include "Decoder.h"
#include "Aligner.h"
#include "system.h"
#include "cpuThreadAlignPool.h"
#include "topology.h"

//---------------------------------------------------------------------------------------------------------
// DEFINITIONS
//---------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------
// GLOBAL VARIABLES
//---------------------------------------------------------------------------------------------------------
struct timeval ts;
int verbose = 0;
volatile unsigned int gTagCnt = 0U;
cpuAlignerPool_t * CPU_Pools = NULL;


//---------------------------------------------------------------------------------------------------------
// STATIC LOCAL VARIABLES
//---------------------------------------------------------------------------------------------------------
static FILE* logfile = NULL;
static int gWriteQS;
static const char DefaultPrefix[] = "/data6/b37/tbl_18nt/start0.";
static const char DefaultGenomeFile[] = "/data6/b37/b37.bin";
static unsigned char gMinScore = '#';
static volatile int StopTerminal = 0;
static SystemInfo sysInfo;
static FASTQ ReaderInputs;
static DecoderPool_t DecoderPool;
static DecoderMemoryBulk_t * restrict DecoderInOut = NULL;
static AlignerPool_t * restrict AlignerPool = NULL;
static WriterPool_t * restrict WriterPool = NULL;
#ifdef USE_MPI
static TransferMemoryBulk_t * restrict TransferMemory = NULL;
#endif

//---------------------------------------------------------------------------------------------------------
// LOCAL FUNCTION
//---------------------------------------------------------------------------------------------------------
void __attribute__((noreturn)) usage()
{
	printf("usage:\n\n");
	printf("MicMap -1 FASTQ -2 FASTQ -r resultfile\n\n");
	printf("           -g genomefile         : full path of the genome encoded by encode_genome\n");
	printf("                                   default: '%s' \n",DefaultGenomeFile);
	printf("           -1 FASTQ              : pair 1\n");
	printf("           -2 FASTQ              : pair 2\n");
	printf("           -r resultfile         : name of the result file (produced by match)\n");
	printf("           -p prefix             : path and prefix of the tables generated by encode_UEM\n");
	printf("                                   default: '%s' \n",DefaultPrefix);
	printf("           -G chunk size         : genome chunk size, default (576)\n");
	printf("           -n nucleotide         : nucleotide to treat\n");
	printf("           -D <uint>             : number of decoder threads, default (4)\n");
	printf("           -d <uint>             : number of input reader blocks, default (3)\n");
	printf("           -A <uint>             : number of threads per alignment blocks, default (4)\n");
	printf("           -a <uint>             : number of input aligner blocks, default (3)\n");
	printf("           -W <uint>             : number of writer threads, default (2)\n");
	printf("           -w <uint>             : number of output aligner blocks, default (3)\n");
	printf("           -q                    : omit writing of quality score\n");
	printf("           -m <char>             : minimum score character in fastq files [%c]\n",gMinScore);
	printf("           -v <level>            : verbose level\n");
	printf("           -O <unsigned long>    : start ordinal position of the fastq reads [0]\n\n");
	printf("           -x instruction set    : enforces use of instruction set\n");
	printf("                                  "
#ifdef HANDLE_SSE41
	" SSE4.1"
#endif
#ifdef HANDLE_AVX
	" AVX"
#endif
#ifdef HANDLE_AVX2
	" AVX2"
#endif
#ifdef HANDLE_AVX512F
	" AVX512F"
#endif
	"\n");
	printf("           -t                    : dump machine topology in json format\n");
	printf("           -T json file          : read json topology file\n");
#ifndef USE_MPI
	printf("           -s shared location    : stores genome at that location\n");
#endif
	printf("\n");
	printf(GTL_VERSION_FULL_STRING "\n");
	printf("(c) Nicolas Guex, Christian Iseli & Thierry Schuepbach 2014-2020\n");
	exit(1);
}
//---------------------------------------------------------------

/*
- Position the Cursor:
  \033[<L>;<C>H Or \033[<L>;<C>f
  puts the cursor at line L and column C.
- Move the cursor up N lines:
  \033[<N>A
- Move the cursor down N lines:
  \033[<N>B
- Move the cursor forward N columns:
  \033[<N>C
- Move the cursor backward N columns:
  \033[<N>D

- Clear the screen, move to (0,0):
  \033[2J
- Erase to end of line:
  \033[K

- Save cursor position:
  \033[s
- Restore cursor position:
  \033[u
  */
void* TerminalViewer(void* data)
{
	char format[128];
	static const char Line[256] = { [0 ... 255] = '-'};
	static const char * what[] = { "Written pairs of tags", "Written tags"};
	struct timeval StartingTime, Now;
	const char Title[] = "MICMAP - " GTL_VERSION_FULL_STRING;
	struct winsize w;

	const char * const restrict whatPtr = (ReaderInputs.PairedEnd) ? what[0] : what[1];
	gettimeofday(&StartingTime, NULL);
	sleep(20);

	do {
		ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
		const int center = (((int) w.ws_col) - sizeof(Title))/2;
		snprintf(format, 128, "\033[2;0H\033[K\033[%uC%s\033[K", center, Title);

		fprintf(stdout, "\033[s\033[1;0H+%.*s+", (int) (w.ws_col-2 > 256 ? 256 : w.ws_col-2), Line);
		fputs(format, stdout);
		fprintf(stdout, "\033[3;0H\033[K+%.*s+", (int) (w.ws_col-2 > 256 ? 256 : w.ws_col-2), Line);
		//fputs("\033[3;0H\033[K", stdout);
		fprintf(stdout, "\033[4;0H\033[K\033[1B\033[K\033[3CReader blocks                : %u waiting for decoding, %u available out of %u\033[K",
		        ((volatile int) (ReaderInputs.ReaderJobQueue.len)), ((volatile int) (ReaderInputs.ReaderMemorySlot.len)), gTopology.nReaderBlocks);
		fprintf(stdout, "\033[6;0H\033[K\033[3CInternal decoder pool blocks : %u used, %u available\033[K",
		        ((volatile int) (DecoderPool.PoolJob_q.len)), ((volatile int) (DecoderPool.PoolMemory_q.len)));
		fprintf(stdout, "\033[7;0H\033[K\033[3CDecoder blocks               : %u waiting for alignment, %u available out of %u\033[K",
                        ((volatile int) (DecoderInOut->ToBeAligned_q.len)), ((volatile int) (DecoderInOut->DecoderMemory_q->len)), gTopology.nDecoderBlocks);
		fprintf(stdout, "\033[8;0H\033[K\033[3CAligner blocks               : %u waiting for storage\033[K",
                        ((volatile int) (AlignerPool->ToBeWritten_q.len)));
		if (CPU_Pools) {
			fprintf(stdout, "\033[9;0H\033[K\033[3CCPU Aligners pool blocks     : %u waiting input, %u available output\033[K",
			        ((volatile int) (CPU_Pools->job_q.len)), ((volatile int) (CPU_Pools->space_q.len)));
		}
		else {
			fputs("\033[9;0H\033[K", stdout);
		}
		fprintf(stdout, "\033[10;0H\033[K\033[3CWriter blocks                : %u waiting for storage\033[K",
                        ((volatile int) (WriterPool->ToBeProcessed_q->len)));
		fputs("\033[11;0H\033[K", stdout);

		const unsigned long Count = atomic_load(&(WriterPool->WrittenAlignments));
		gettimeofday(&Now, NULL);
		const double elapsetime = (double) ((1.0/60.0)*(Now.tv_sec - StartingTime.tv_sec));
		const double Throughput = (double) Count/elapsetime;
		fprintf(stdout, "\033[12;0H\033[K%s / minute : %lf (%lu/%lf)", whatPtr, Throughput, Count, elapsetime);
		fputs("\033[13;0H\033[K", stdout);
		fprintf(stdout, "\033[14;0H+%.*s+\033[13;0H\033[K\033[u", (int) (w.ws_col-2 > 256 ? 256 : w.ws_col-2), Line);

		fflush(stdout);
		//usleep((useconds_t)10000);
		sleep(1);

	} while (StopTerminal == 0);
	return NULL;
}
//---------------------------------------------------------------

void segfault_handler(int sig) {
  void *array[32];

  // get void*'s for all entries on the stack
  const size_t size = backtrace(array, 32);

  // print out all the frames to stderr
	if (logfile == NULL)
		logfile = stderr;
  fprintf(logfile, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, fileno(logfile));
	time_t t = time(NULL);
	struct tm tm = *localtime(&t);
	fprintf(logfile, "Terminated with error on %d-%d-%d %d:%2.2d:%2.2d.\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
	fclose(logfile);
  exit(1);
}
//---------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------
// MAIN FUNCTION
//---------------------------------------------------------------------------------------------------------

int main (int argc, char *argv[])
{
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Variables
	const char * restrict prefix = DefaultPrefix;
	DecodeTagsPtr DecodeFunction = NULL;
	void * libMicMap = NULL;
	int c;
	int f = 0;
	struct stat s;
	struct stat g;

	char ctmp[512];
	char * restrict rsltfile = NULL;
	char NT = 'A';
	int err = 1;
	int ShowTopology = 0;

#ifndef USE_MPI
	const char * SharedLocation = NULL;
#endif
	Genome_t Genome = {
		.FileName = &DefaultGenomeFile[0],
		.table = NULL,
		.configBuffer = NULL,
		.tableSize = 0UL,
		.configBufferSize = 0UL
	};

	unsigned int GenomeChunkSize = 576U;
	int SingleEnd = 0;  // 0 by default, otherwise treat data as single end, and bypass all paired end treatment.

	char * TopologyConfigurationFile = NULL;
	int PrintOutTopologyTemplate = 0;

	/* Terminal viewer */
	pthread_t * restrict TerminalViewerThread = NULL;

	/* Readers */
	pthread_t ReaderThread, GenomeLoader;
	pthread_attr_t ReaderThreadAttr;
	/* Decoders */
	Affinity_Mask_t * restrict DecoderPoolAffinities;
	/* Readers to decoders */
	ReaderToDecoderArgs_t DecoderCommArgs;
	pthread_t ReaderToDecoder;
#ifdef USE_AFFINITY
	pthread_attr_t ReaderToDecoderAttr;
#endif
	/* Transfer */
#ifdef USE_MPI
	pthread_t Receiver;
#endif
	/* Aligners */
	pthread_t * restrict AlignerThread;
	AlignerBlockArgs_t * restrict AlignerThreadArgs;
#ifdef USE_AFFINITY
	pthread_attr_t * restrict AlignerThreadAttr;
#endif

	enum InstructionSet_t { GENERIC, SSE41, AVX, AVX2, AVX512F, AVX512BW } InstructionSet = GENERIC;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Initialization
	struct timeval StartingTime, EndingTime;
	gettimeofday(&StartingTime, NULL);
#ifdef USE_MPI
	{
		//volatile int i = 0;
		char hostname[256];
		gethostname(hostname, sizeof(hostname));
		printf("PID %d on %s ready for attach\n", getpid(), hostname);
		fflush(stdout);
		//while (0 == i) sleep(5);
	}
	MPI_Init(&argc, &argv);

	int world_rank, world_size;
	if (MPI_Comm_rank(MPI_COMM_WORLD, &world_rank)!=MPI_SUCCESS) {
		fputs("MPI Node failed to get its rank!\n", stderr);
		goto bail;
	}
	if (MPI_Comm_size(MPI_COMM_WORLD, &world_size) != MPI_SUCCESS) {
		fputs("MPI Node failed to get communication size!\n", stderr);
	}
#endif
	gCompress.def = CMP_LZ4;
	gCompress.scoredef = CMP_LZ4;

	gWriteQS = 1;
	memset(&ReaderInputs, 0, sizeof(FASTQ));
	ReaderInputs.MinimalQualityScore = gMinScore;

	getSystemInfo(&sysInfo);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Process command line arguments
	opterr = 0;
	while ((c = getopt (argc, argv, "1:2:v:r:n:p:g:qA:a:D:d:W:w:m:O:G:ItT:x:X:s:")) != -1)
	switch (c)
	{
		case '1':
			ReaderInputs.tagfile1 = optarg;
			break;
		case '2':
			ReaderInputs.tagfile2 = optarg;
			break;
		case 'n':
			sscanf(optarg,"%c",&NT);
			break;
		case 'r':
			rsltfile = optarg;
			break;
		case 'p':
			prefix = optarg;
			break;
		case 'g':
			Genome.FileName = optarg;
			break;
#ifndef USE_MPI
		case 's':
			SharedLocation = optarg;
			break;
#endif
		case 'G':
			{
				int dummy;
				if (sscanf(optarg,"%i",&dummy) != 1) {
					fprintf(stderr, "Cannot read gGenomeChunkSize in %s\n", optarg);
					return 1;
				}
				if (dummy <= 0) {
					fprintf(stderr, "gGenomeChunkSize must be positive (%i)\n", dummy);
					return 1;
				}
				GenomeChunkSize = (unsigned int) dummy;
			}
			break;
		case 'q':
			gWriteQS = 0;
			break;
		case 'A':
			{
				int dummy;
				sscanf(optarg,"%i",&dummy);
				if (dummy > 0) {
					gTopology.nAlignerThreads = (unsigned int) dummy;
				}
				else {
					fputs("Please provide a valid number of threads per aligner\n", stderr);
					return 1;
				}
			}
			break;
		case 'a':
			{
				int dummy;
				sscanf(optarg,"%i",&dummy);
				if (dummy > 0) {
					gTopology.nAlignerInBlocks = (unsigned int) dummy;
				}
				else {
					fputs("Please provide a valid number of input aligner blocks\n", stderr);
					return 1;
				}
			}
			break;
		case 'D':
			{
				int dummy;
				sscanf(optarg,"%i",&dummy);
				if (dummy > 0) {
					gTopology.nDecoders = (unsigned int) dummy;
				}
				else {
					fputs("Please provide a valid number of decoders\n", stderr);
					return 1;
				}
			}
			break;
		case 'd':
			{
				int dummy;
				sscanf(optarg,"%i",&dummy);
				if (dummy > 0) {
					gTopology.nReaderBlocks = (unsigned int) dummy;
				}
				else {
					fputs("Please provide a valid number of input blocks\n", stderr);
					return 1;
				}
			}
			break;
		case 'W':
			{
				int dummy;
				sscanf(optarg,"%i",&dummy);
				if (dummy > 0) {
					gTopology.nWriters = (unsigned int) dummy;
				}
				else {
					fputs("Please provide a valid number of decoders\n", stderr);
					return 1;
				}
			}
			break;
		case 'w':
			{
				int dummy;
				sscanf(optarg,"%i",&dummy);
				if (dummy > 0) {
					gTopology.nAlignerOutBlocks = (unsigned int) dummy;
				}
				else {
					fputs("Please provide a valid number of output aligner blocks\n", stderr);
					return 1;
				}
			}
			break;
		case 'X':
			{
				ssize_t dummy;
				sscanf(optarg,"%li",&dummy);
				if (dummy > 0) {
					gTopology.nAlignerOutBlockSize = (off_t) dummy;
				}
				else {
					fputs("Please provide a valid size for aligner output block\n", stderr);
					return 1;
				}
			}
			break;
		case 'v':
			verbose = atoi(optarg);
			break;
		case 'm':
			ReaderInputs.MinimalQualityScore = (char) ( '#' + (unsigned char) atoi(optarg));
			break;
		case 'O':
			ReaderInputs.startOrdinal = strtoul(optarg, NULL, 0);
			break;
		case 'I':
			ShowTopology = 1;
			break;
		case 't':
			PrintOutTopologyTemplate = 1;
			break;
		case 'T':
			TopologyConfigurationFile = optarg;
			break;
		case 'x':
			{
				if (strncmp(optarg, "SSE4.1", 6) == 0) {
#ifdef HANDLE_SSE41
					if (sysInfo.Extensions & MM_SSE41) {
						InstructionSet = SSE41; goto FOUND;
					}
					else {
						fputs("Current architecture does not support SSE 4.1\n", stderr);
						exit(1);
					}
#else
					fprintf(stderr, "Instruction set %s not built, sorry\n", optarg);
					exit(1);
#endif
				}

				if (strncmp(optarg, "AVX2", 4) == 0) {
#ifdef HANDLE_AVX2
					if (sysInfo.Extensions & MM_AVX2) {
						InstructionSet = AVX2; goto FOUND;
					}
					else {
						fputs("Current architecture does not support AVX2\n", stderr);
						exit(1);
					}
#else
					fprintf(stderr, "Instruction set %s not built, sorry\n", optarg);
					exit(1);
#endif
				}

				if (strncmp(optarg, "AVX512F", 7) == 0) {
#ifdef HANDLE_AVX512F
					if (sysInfo.Extensions & MM_AVX512F) {
						InstructionSet = AVX512F; goto FOUND;
					}
					else {
						fputs("Current architecture does not support AVX 512F\n", stderr);
						exit(1);
					}
#else
					fprintf(stderr, "Instruction set %s not built, sorry\n", optarg);
					exit(1);
#endif
				}

				if (strncmp(optarg, "AVX\0", 4) == 0) {
#ifdef HANDLE_AVX
					if (sysInfo.Extensions & MM_AVX) {
						InstructionSet = AVX; goto FOUND;
					}
					else {
						fputs("Current architecture does not support AVX\n", stderr);
						exit(1);
					}
#else
					fprintf(stderr, "Instruction set %s not built, sorry\n", optarg);
					exit(1);
#endif
				}

				fprintf(stderr, "Unrecognized instruction set %s\n", optarg);
				exit(1);

				FOUND:;
			}
			break;
		default:
			fprintf(stderr,"Argument %s is unknown", optarg);
			exit(1);
	}
	if (ReaderInputs.tagfile2 == NULL)
		SingleEnd = 1;
	else
		ReaderInputs.PairedEnd = true;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Read Topology
	if (TopologyConfigurationFile) {
		if (loadTopology(&sysInfo, TopologyConfigurationFile)) {
				fprintf(stderr,"Error in topology file %s\n", TopologyConfigurationFile);
				exit(1);
		}
		if (verbose) printParameters(stdout);
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Dump Topology
	if (PrintOutTopologyTemplate) {
		printJSON(&sysInfo, stdout);
		exit(0);
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Store genome in memory
#ifndef USE_MPI
	if (SharedLocation) {
		if(DecodingTableToRAM(SharedLocation, prefix, NT)) {
			fputs("Error while loading decoding table to shared memory RAM\n", stderr);
			exit(1);
		}
		if(GenomeToRAM(SharedLocation, Genome.FileName)) {
			fputs("Error while loading genome to shared memory RAM\n", stderr);
			exit(1);
		}
		exit(0);
	}
#endif
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Open Log file, print header
#ifdef USE_MPI
	if (world_rank == 0)
#endif
	{
		if (rsltfile != NULL) {
			snprintf(ctmp, sizeof(ctmp), "%s_log.txt", rsltfile);
			logfile = fopen(ctmp, "w");
			if (logfile == NULL) {
				fprintf(stderr, "Unable to create log file %s\n", ctmp);
				exit(1);
			}
		} else
			logfile = stderr;

		signal(SIGSEGV, segfault_handler);
		signal(SIGINT, segfault_handler);
		signal(SIGABRT, segfault_handler);

		time_t t = time(NULL);
		struct tm tm = *localtime(&t);

		fprintf(logfile,
						"*------------------------------------------------------------------------------*\n"
					  "|                               MICMAP v %i.%i                                   |\n"
					  "|    (c) 2014-2020  Nicolas Guex, Christian Iseli & Thierry Schuepbach         |\n"
					  "*------------------------------------------------------------------------------*\n\n"
					   GTL_VERSION_FULL_STRING "\n\n", GTL_MAJOR_VERSION, GTL_MINOR_VERSION);
		printSystemInfo(logfile, &sysInfo);
		fputc('\n', logfile);
		printParameters(logfile);

		fprintf(logfile, "\n-------------------------- GENOME ---------------------------\n"
		                 "File       : %s\n"
										 "Prefix     : %s\n"
										 "Chunk size : %u\n",
						Genome.FileName, prefix, GenomeChunkSize);

		fprintf(logfile, "\nStarting analysis on %d-%d-%d %d:%2.2d:%2.2d.\n", tm.tm_year + 1900, tm.tm_mon + 1,
						tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
		fflush(logfile);
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Informations
#ifdef USE_MPI
	if (world_rank == 0)
#endif
	{
		/* Clear screen */
		if ((verbose & 0x2) && isatty(fileno(stdout))) fputs("\033[2J", stdout);
		/* Print system information */
		printSystemInfo(stdout, &sysInfo);
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Verifications
	if (GenomeChunkSize > kMaxGenomeChunkSize) {
		fprintf(stderr, "WARN: using maximum genome size of %d instead of %d\n",kMaxGenomeChunkSize,GenomeChunkSize);
		fprintf(logfile, "WARN: using maximum genome size of %d instead of %d\n",kMaxGenomeChunkSize,GenomeChunkSize);
		GenomeChunkSize = kMaxGenomeChunkSize;
	}

#ifdef USE_AFFINITY
	if (sysInfo.nOverallCores > 8UL*sizeof(Affinity_Mask_t)) {
		fprintf(stderr, "%s has more available cores (%u) tham affinity can handle (%lu), rebuild!\n", sysInfo.Nodename, sysInfo.nOverallCores, sizeof(Affinity_Mask_t));
		fprintf(logfile, "%s has more available cores (%u) tham affinity can handle (%lu), rebuild!\n", sysInfo.Nodename, sysInfo.nOverallCores, sizeof(Affinity_Mask_t));
		goto bail;
	}
#endif

	if ( (ReaderInputs.tagfile1 == NULL) || (rsltfile == NULL) ) usage();

	fflush(stdout);

#ifdef USE_MPI
	if (world_rank == 0)
#endif
	{
		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Prepare Reader memory slot queue
#ifdef USE_MPI
		printf("%s is the master node\n", sysInfo.Nodename);
#endif
		if ( allocateReaderMemory(&ReaderInputs, gTopology.nReaderBlocks) != 0) {
			fputs("Error preparing the reader!\n", stderr);
			fputs("Error preparing the reader!\n", logfile);
			goto bail;
		}

		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Start Reader

		pthread_attr_init(&ReaderThreadAttr);
		if (pthread_attr_setstacksize(&ReaderThreadAttr, (size_t) (3*1024*1024)) != 0) {
			fputs("Unable to assign 3MB stack size to Reader thread\n", stderr);
			fputs("Unable to assign 3MB stack size to Reader thread\n", logfile);
			goto bail;
		}
#ifdef USE_AFFINITY
		if (gTopology.Reader) {
			if (pthread_attr_setaffinity_np(&ReaderThreadAttr, sizeof(Affinity_Mask_t), (cpu_set_t*) gTopology.Reader) ) {
				fputs("Unable to set affinity to Reader thread\n", stderr);
				fputs("Unable to set affinity to Reader thread\n", logfile);
				goto bail;
			}
		}
#endif
		if (pthread_create(&ReaderThread, &ReaderThreadAttr, (void* (*)(void*)) readFASTQ, (void*) &ReaderInputs) != 0) {
			fputs("Error starting FASTQ reader thread\n", stderr);
			fputs("Error starting FASTQ reader thread\n", logfile);
			goto bail;
		}

		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Prepare decoder slot queue and global tables
		DecoderPool.ValidDataTable = loadValidTable(prefix, NT);
		if (DecoderPool.ValidDataTable == NULL) {
			fputs("Error loading Valid Table\n", stderr);
			fputs("Error loading Valid Table\n", logfile);
			goto bail;
		}
		if (loadDecodingTable(&DecoderPool, prefix, NT, true) != 0) {
			fputs("Error loading Decoding Table\n", stderr);
			fputs("Error loading Decoding Table\n", logfile);
			goto bail;
		}

		DecoderInOut = allocateDecoderMemory(gTopology.nDecoderBlocks);
		switch (InstructionSet) {
#ifdef HANDLE_AVX2
			case AVX2:
				DecodeFunction = DecodeTags_avx2;
				break;
#endif
			case SSE41:
			default:
				DecodeFunction = DecodeTags;
		}

#ifdef USE_AFFINITY
		if (gTopology.Decoders) {
			DecoderPoolAffinities = gTopology.Decoders->Slaves;
		}
		else {
			DecoderPoolAffinities = NULL;
		}
#else
		DecoderPoolAffinities = NULL;
#endif
		const int res = createDecoderPool(&DecoderPool, DecodeFunction, DecoderPoolAffinities, DecoderInOut, gTopology.nDecoders);
		if (DecoderInOut == NULL || res != 0) {
			fputs("Error preparing decoders\n", stderr);
			fputs("Error preparing decoders\n", logfile);
			goto bail;
		}

		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Start Decoders
		DecoderCommArgs.ReaderData = &ReaderInputs;
		DecoderCommArgs.DecoderPool = &DecoderPool;
		DecoderCommArgs.DecoderMemory_q = DecoderInOut->DecoderMemory_q;
#ifdef USE_AFFINITY
		pthread_attr_t * restrict ReaderToDecoderAttrPtr;
		if (gTopology.Decoders) {
			pthread_attr_init(&ReaderToDecoderAttr);
			if (pthread_attr_setaffinity_np(&ReaderToDecoderAttr, sizeof(Affinity_Mask_t), (cpu_set_t*) &(gTopology.Decoders->Master))) {
				fputs("Unable to set affinity to ReaderToDecoder thread\n", stderr);
				fputs("Unable to set affinity to ReaderToDecoder thread\n", logfile);
				goto bail;
			}
			ReaderToDecoderAttrPtr = &ReaderToDecoderAttr;
		}
		else
			ReaderToDecoderAttrPtr = NULL;
#else
		pthread_attr_t * const ReaderToDecoderAttrPtr = NULL;
#endif
		if (pthread_create(&ReaderToDecoder, ReaderToDecoderAttrPtr, (void* (*)(void*)) ReaderToDecoderThread, (void*) &DecoderCommArgs) != 0) {
			fputs("Error starting reader to decoder thread\n", stderr);
			fputs("Error starting reader to decoder thread\n", logfile);
			goto bail;
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Loading genome, virtual chromosomes
#ifdef USE_MPI
	const _Bool SlaveLoadGenome = true;
	if (world_rank == 0 || SlaveLoadGenome )
#endif
	{
		if (pthread_create(&GenomeLoader, NULL, (void* (*)(void*)) LoadGenome, &Genome) != 0) {
			fputs("Error starting genome loader thread\n", stderr);
			fputs("Error starting genome loader thread\n", logfile);
			goto bail;
		}

		void* tres;
		pthread_join(GenomeLoader, &tres);
		if ( ((uintptr_t) tres) != (uintptr_t) 0) {
			fputs("Error reading genome\n", stderr);
			fputs("Error reading genome\n", logfile);
			goto bail;
		}
	}
#ifdef USE_MPI
	if (! SlaveLoadGenome) {
		if (world_rank == 0) {
			if(SendGenome(&Genome)) {
				fputs("Master node error sending genome\n", stderr);
				fputs("Master node error sending genome\n", logfile);
				goto bail;
			}
		}
		else {
			if(ReceiveGenome(&Genome)) {
				fprintf(stderr, "Node %i failed receiving genome...\n", world_rank);
				goto bail;
			}
		}
	}
#endif
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Prepare Aligner data: creating pool
#ifdef USE_MPI
	if (world_rank == 0)
#endif
	{
		AlignerPool = createAlignerPool(&(DecoderInOut->ToBeAligned_q));
		if (AlignerPool == NULL) {
			fputs("Error creating aligner pool\n", stderr);
			fputs("Error creating aligner pool\n", logfile);
			goto bail;
		}
	}
#ifdef USE_MPI
	else {
		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		printf("%s is a slave node\n", sysInfo.Nodename);
		// Allocate memory to receive data
		TransferMemory = allocateTransferMemory(gTopology.nReaderBlocks);
		if (TransferMemory == NULL) {
			fprintf(stderr, "Error allocating transfer memory on node %i\n", world_rank);
			goto bail;
		}
		AlignerPool = createAlignerPool(&(TransferMemory->ToBeAligned_q));
		if (AlignerPool == NULL) {
			fputs("Error creating aligner pool\n", stderr);
			goto bail;
		}

		TransferMemory->MPINodeId = world_rank;
		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Start MPI Transfer questioning thread
		if (pthread_create(&Receiver, NULL, (void* (*)(void*)) MPIReceiver, (void*) TransferMemory) != 0) {
			fputs("Error creating MPI Receiver thread\n", stderr);
			goto bail;
		}
	}
#endif

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Start Aligners
	AlignerFct alignFct;
	{
		AlignerThreadArgs = (AlignerBlockArgs_t*) alloca(sizeof(AlignerBlockArgs_t));
		AlignerThread = (pthread_t*) alloca(sizeof(pthread_t));
#ifdef USE_AFFINITY
		if (gTopology.Aligners) {
			AlignerThreadAttr = (pthread_attr_t*) alloca(sizeof(pthread_attr_t));
			pthread_attr_init(AlignerThreadAttr);
			if (pthread_attr_setaffinity_np(AlignerThreadAttr, sizeof(Affinity_Mask_t), (cpu_set_t*) &(gTopology.Aligners->Master))) {
				fputs("Unable to set master aligner thread affinity\n", stderr);
				fputs("Unable to set master aligner thread affinity\n", logfile);
				goto bail;
			}
			AlignerThreadArgs[0].Affinities = gTopology.Aligners->Slaves;
		}
		else {
			AlignerThreadAttr = NULL;
			AlignerThreadArgs[0].Affinities = NULL;
		}
#else
		AlignerThreadAttr = NULL;
#endif
		alignFct = cpuBlockAligners;
		switch (InstructionSet) {
#ifdef HANDLE_AVX512F
			case AVX512F:
				AlignerThreadArgs[0].ExtraPtr = (void*) avx512f_pair;
				break;
#endif
			case SSE41:
			default:
				AlignerThreadArgs[0].ExtraPtr = (void*) sse41_pair;
 		}
	}
	AlignerThreadArgs[0].SingleEnd = 0;
	AlignerThreadArgs[0].GenomeChunkSize = GenomeChunkSize;
	AlignerThreadArgs[0].nAligners = gTopology.nAlignerThreads;
	AlignerThreadArgs[0].nInputBlock = gTopology.nAlignerInBlocks;
	AlignerThreadArgs[0].nOutputBlock = gTopology.nAlignerOutBlocks;
	AlignerThreadArgs[0].OutputBlockSize = gTopology.nAlignerOutBlockSize;
#ifdef USE_MPI
	if (world_rank == 0) {
		AlignerThreadArgs[0].Pool = AlignerPool;
		AlignerThreadArgs[0].ReaderMemoryQueue = &(ReaderInputs.ReaderMemorySlot);
		AlignerThreadArgs[0].DecoderMemoryQueue = DecoderInOut->DecoderMemory_q;
		AlignerThreadArgs[0].DecoderInternalJob_q = &(DecoderPool.PoolMemory_q);
		AlignerThreadArgs[0].Genome = &Genome;
		AlignerThreadArgs[0].ID = world_rank;

		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Start MPI Transfer Listening thread
		printf("Master starting MPI dispatcher\n"); fflush(stdout);
		if (pthread_create(&AlignerThread[0], NULL, (void* (*)(void*)) MPIDispatch, (void*) &AlignerThreadArgs[0]) != 0) {
			fputs("Error creating MPI Dsipatcher thread\n", stderr);
			fputs("Error creating MPI Dsipatcher thread\n", logfile);
			goto bail;
		}
	}
	else {
		AlignerThreadArgs[0].Pool = AlignerPool;
		AlignerThreadArgs[0].ReaderMemoryQueue = NULL;
		AlignerThreadArgs[0].DecoderMemoryQueue = NULL;
		AlignerThreadArgs[0].DecoderInternalJob_q = &(TransferMemory->Memory_q);
		AlignerThreadArgs[0].Genome = &Genome;
		AlignerThreadArgs[0].ID = world_rank;
		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Start MPI slaves aligner pool
		printf("Starting %u aligners on node %i\n", AlignerThreadArgs[0].nAligners, world_rank);
		if (pthread_create(&AlignerThread[0], NULL, (void* (*)(void*)) alignFct, (void*) &AlignerThreadArgs[0]) != 0) {
			fputs("Error creating aligner thread\n", stderr);
			goto bail;
		}
	}
#else
	{
		AlignerThreadArgs[0].Pool = AlignerPool,
		AlignerThreadArgs[0].ReaderMemoryQueue = &(ReaderInputs.ReaderMemorySlot);
		AlignerThreadArgs[0].DecoderMemoryQueue = DecoderInOut->DecoderMemory_q;
		AlignerThreadArgs[0].DecoderInternalJob_q = &(DecoderPool.PoolMemory_q);
		AlignerThreadArgs[0].Genome = &Genome;
		AlignerThreadArgs[0].ID = 0;

		if (pthread_create(&AlignerThread[0], AlignerThreadAttr, (void* (*)(void*)) alignFct, (void*) &AlignerThreadArgs[0]) != 0) {
				fputs("Error creating aligner thread\n", stderr);
				fputs("Error creating aligner thread\n", logfile);
				goto bail;
		}
	}
#endif

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Prepare writers slot queue and global tables
#ifndef USE_MPI
	COUNTS * const restrict cnts = (COUNTS*) alloca(gTopology.nWriters*sizeof(COUNTS));
#else
	const register int ln = (world_rank == 0) ? (world_size-1) : gTopology.nWriters;
	COUNTS * const restrict cnts = (COUNTS*) alloca(gTopology.nWriters*sizeof(COUNTS));
#endif
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Start Writers
#ifdef USE_MPI
	char LocalNodeRsltBaseFileName[256];
	if (world_rank != 0)
#endif
	{
#ifdef USE_MPI
		snprintf(&LocalNodeRsltBaseFileName[0], 256, "%s_%i", rsltfile, world_rank);
		rsltfile = LocalNodeRsltBaseFileName;
#endif
		const cpu_pool_t * restrict WriterThreadAttrPtr;
#ifdef USE_AFFINITY
		WriterThreadAttrPtr = (gTopology.Writers) ? gTopology.Writers : NULL;
#else
		WriterThreadAttrPtr = NULL;
#endif
		WriterPool = createWriterPool(&Genome, &(AlignerPool->ToBeWritten_q), rsltfile,
		                              cnts, WriterThreadAttrPtr, ReaderInputs.PairedEnd,
		                              gTopology.nWriters, gTopology.nWriterThreads);
	}
#ifdef USE_MPI
	else
#endif
	{
		////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Start screen data thread
		if ((verbose & 0x2) && isatty(fileno(stdout))) {
			TerminalViewerThread = (pthread_t*) alloca(sizeof(pthread_t));
			if (pthread_create(TerminalViewerThread, NULL, TerminalViewer, NULL)) {
					fprintf(stderr, "Unable to start the terminal viewer\n");
					fprintf(logfile, "Unable to start the terminal viewer\n");
			}
		}

		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Wait for readers to terminate
		void* result;
		pthread_join(ReaderThread, &result);
		if ( (intptr_t) result != (intptr_t) 0 ) {
			fputs("Error in reader thread\n", stderr);
			fputs("Error in reader thread\n", logfile);
			goto bail;
		}

		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Wait for decoders to terminate
		pthread_join(ReaderToDecoder, NULL);
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Wait for aligners to terminate
	waitForAlignerPool(AlignerPool);
	terminateAlignerPool(AlignerPool);
	{
		void* result;
		register unsigned int i=0;
		do {
			pthread_join(AlignerThread[i], &result);
			if ( (intptr_t) result != (intptr_t) 0 ) {
				fputs("Error in main block aligner thread\n", stderr);
				fputs("Error in main block aligner thread\n", logfile);
				goto bail;
			}
		} while (++i < (gTopology.nAligners));
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Wait for writers to terminate
#ifdef USE_MPI
	if (world_rank != 0)
#endif
	{
#ifdef USE_MPI
		/* Wait for Receiver thread to terminate */
		pthread_join(Receiver, NULL);
#endif
		/*while (atomic_load(&(WriterPool->num_threads_alive)) > 0)  {
			sleep(1U);
		}*/
		terminateWriterPool(WriterPool);
		StopTerminal = 1;
		for (unsigned int i=0; i<gTopology.nWriters; i++) pthread_join(WriterPool->threads[i], NULL);
		cumulateStatistics(cnts, gTopology.nWriters);
#ifndef USE_MPI
		if (TerminalViewerThread) pthread_join(*TerminalViewerThread, NULL);

		/* Report */
		gettimeofday(&EndingTime, NULL);
		reportStatistics(cnts, logfile, EndingTime.tv_sec - StartingTime.tv_sec);
		report_compress_stat(logfile);
#else
		/* Send statistics to master node */
		MPI_Send(cnts, sizeof(COUNTS), MPI_CHAR, 0, 0, MPI_COMM_WORLD);
#endif
	}
#ifdef USE_MPI
	else {
		/* Receive all nodes'statistics */
		MPI_Status status;
		for (int iSlave=1; iSlave<world_size; iSlave++) {
			MPI_Recv(&cnts[iSlave-1], sizeof(COUNTS), MPI_CHAR, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
		}

		/* Report */
		cumulateStatistics(cnts, world_size-1);
		reportStatistics(cnts, logfile, EndingTime.tv_sec - StartingTime.tv_sec);
	}
#endif

	err = 0;
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Terminate
bail:;
	if (logfile) {
		time_t t = time(NULL);
		struct tm tm = *localtime(&t);
		if (err == 0) {
			fprintf(logfile, "Succesfully terminated on %d-%d-%d %d:%2.2d:%2.2d.\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
		}
		else {
			fprintf(logfile, "Terminated with error on %d-%d-%d %d:%2.2d:%2.2d.\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
		}
		fclose(logfile);
	}
#ifdef USE_MPI
	if (world_rank == 0)
#endif
	{
		freeReaderMemory(&ReaderInputs);
		freeValidTable(DecoderPool.ValidDataTable);
		freeDecodingTable(&DecoderPool);
	}
#ifdef USE_MPI
	else
#endif
	if (WriterPool) freeWriterPool(WriterPool);

#ifdef USE_MPI
	if (err) MPI_Abort(MPI_COMM_WORLD,1);
	MPI_Finalize();
#endif

	return err;
}
//---------------------------------------------------------------

/* vim: tabstop=2 shiftwidth=2
 */
