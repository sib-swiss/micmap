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
#ifndef _MPI_TRANSFER_H
#define _MPI_TRANSFER_H
#include "Reader.h"
#include "Decoder.h"
#include "Aligner.h"

typedef struct TransferMemoryBulk {
	jobqueue_t ToBeAligned_q;              /* Pool of aligner masters                      */
	memqueue_t Memory_q;                   /* empty memory slot queue                      */
	DecoderJobElement_t * jobs;            /* pointer to all jobs memory, only for freeing */
	char * ReaderDecoderMemory;            /* array of decoder data memory                 */
	size_t MemorySize;                     /* allocated space with mmap                    */
	int MPINodeId;                         /* MPI node number                              */
} TransferMemoryBulk_t;

int ReceiveGenome(Genome_t * const restrict genome);
int SendGenome(Genome_t * const restrict genome);
void* MPIDispatch(AlignerBlockArgs_t * const restrict data);
TransferMemoryBulk_t * allocateTransferMemory(const size_t nBlocks);
void* MPIReceiver(TransferMemoryBulk_t * const restrict data);

#endif /* _MPI_TRANSFER_H */

/* vim: tabstop=2 shiftwidth=2
 */
/* ------------------------------------------------------------------------------------ */
