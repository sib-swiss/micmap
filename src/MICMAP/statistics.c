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
#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include "statistics.h"

void reportStatistics(const COUNTS * const restrict cnts, FILE * const out, unsigned long elapsed)
{
	COUNTS_TYPE_t PairSum;
	COUNTS_TYPE_t TooFarSum;
	COUNTS_TYPE_t HalfMapSum;
	memset(&PairSum, 0, sizeof(COUNTS_TYPE_t));
	memset(&TooFarSum, 0, sizeof(COUNTS_TYPE_t));
	memset(&HalfMapSum, 0, sizeof(COUNTS_TYPE_t));
		
	fprintf(out,
		"\nDistribution on chromosomes"
		"\n   |------------------------------------------- ALIGNMENT STD ------------------------------------------------------------|\n");
	fprintf(out,
		"chr       pppf      mNppf       mppf      mppfN       gppf      gppfN      smppf     smppfN      sgppf     sgppfN   genomeN\n"
	);

	for (int i = 1; i< kMAXCHRcnt;i++)
	{
		const COUNTS_TYPE_t * const ptr = &(cnts->Pair[i]);
		fprintf(out, "%3d  %9u  %9u  %9u  %9u  %9u  %9u  %9u  %9u  %9u  %9u %9u\n", i,
					 ptr->pppf, ptr->mNppf,
					 ptr->mppf, ptr->mppfN,
					 ptr->gppf, ptr->gppfN,
					 ptr->smppf, ptr->smppfN,
					 ptr->sgppf, ptr->sgppfN, ptr->genomeN);
		
		PairSum.pppf    += ptr->pppf;
		PairSum.mNppf   += ptr->mNppf;
		PairSum.mppf    += ptr->mppf;
		PairSum.mppfN   += ptr->mppfN;
		PairSum.gppf    += ptr->gppf;
		PairSum.gppfN   += ptr->gppfN;
		PairSum.smppf   += ptr->smppf;
		PairSum.smppfN  += ptr->smppfN;
		PairSum.sgppf   += ptr->sgppf;
		PairSum.sgppfN  += ptr->sgppfN;
		PairSum.genomeN += ptr->genomeN;
	}
	
	fprintf(out, "\nSUM  %9u  %9u  %9u  %9u  %9u  %9u  %9u  %9u  %9u  %9u %9u\n",
					 PairSum.pppf, PairSum.mNppf,
					 PairSum.mppf, PairSum.mppfN,
					 PairSum.gppf, PairSum.gppfN,
					 PairSum.smppf, PairSum.smppfN,
					 PairSum.sgppf, PairSum.sgppfN, PairSum.genomeN);
	
	const unsigned int MappedOK = PairSum.pppf + PairSum.mNppf + PairSum.mppf + PairSum.mppfN\
	                            + PairSum.gppf + PairSum.gppfN + PairSum.smppf + PairSum.smppfN\
	                            + PairSum.sgppf + PairSum.sgppfN;
	
	fprintf(out,
		"\n   |------------------------------------------- TOO FAR AWAY -------------------------------------------------------------|\n"
		"chr       pppf      mNppf       mppf      mppfN       gppf      gppfN      smppf     smppfN      sgppf     sgppfN   genomeN\n"
	);
	for (int i = 1; i< kMAXCHRcnt;i++)
	{
		const COUNTS_TYPE_t * const ptr = &(cnts->TooFar[i]);
		fprintf(out, "%3d  %9u  %9u  %9u  %9u  %9u  %9u  %9u  %9u  %9u  %9u %9u\n", i,
					 ptr->pppf, ptr->mNppf,
					 ptr->mppf, ptr->mppfN,
					 ptr->gppf, ptr->gppfN,
					 ptr->smppf, ptr->smppfN,
					 ptr->sgppf, ptr->sgppfN, ptr->genomeN);
		
		TooFarSum.pppf    += ptr->pppf;
		TooFarSum.mNppf   += ptr->mNppf;
		TooFarSum.mppf    += ptr->mppf;
		TooFarSum.mppfN   += ptr->mppfN;
		TooFarSum.gppf    += ptr->gppf;
		TooFarSum.gppfN   += ptr->gppfN;
		TooFarSum.smppf   += ptr->smppf;
		TooFarSum.smppfN  += ptr->smppfN;
		TooFarSum.sgppf   += ptr->sgppf;
		TooFarSum.sgppfN  += ptr->sgppfN;
		TooFarSum.genomeN += ptr->genomeN;
	}
	
	fprintf(out, "\nSUM  %9u  %9u  %9u  %9u  %9u  %9u  %9u  %9u  %9u  %9u %9u\n",
					 TooFarSum.pppf, TooFarSum.mNppf,
					 TooFarSum.mppf, TooFarSum.mppfN,
					 TooFarSum.gppf, TooFarSum.gppfN,
					 TooFarSum.smppf, TooFarSum.smppfN,
					 TooFarSum.sgppf, TooFarSum.sgppfN, TooFarSum.genomeN);
	
	
	const unsigned int MappedTooFar = TooFarSum.pppf + TooFarSum.mNppf + TooFarSum.mppf + TooFarSum.mppfN\
	                                + TooFarSum.gppf + TooFarSum.gppfN + TooFarSum.smppf + TooFarSum.smppfN\
	                                + TooFarSum.sgppf + TooFarSum.sgppfN;
	
	
	fprintf(out, 
		"\n   |-------------------------------------------- HALF MAPPED -------------------------------------------------------------|\n"
		"chr       pppf      mNppf       mppf      mppfN       gppf      gppfN      smppf     smppfN      sgppf     sgppfN   genomeN\n"
	);
	for (int i = 1; i< kMAXCHRcnt;i++)
	{
		const COUNTS_TYPE_t * const ptr = &(cnts->HalfMap[i]);
		fprintf(out, "%3d  %9u  %9u  %9u  %9u  %9u  %9u  %9u  %9u  %9u  %9u %9u\n", i,
					 ptr->pppf, ptr->mNppf,
					 ptr->mppf, ptr->mppfN,
					 ptr->gppf, ptr->gppfN,
					 ptr->smppf, ptr->smppfN,
					 ptr->sgppf, ptr->sgppfN, ptr->genomeN);
		
		HalfMapSum.pppf    += ptr->pppf;
		HalfMapSum.mNppf   += ptr->mNppf;
		HalfMapSum.mppf    += ptr->mppf;
		HalfMapSum.mppfN   += ptr->mppfN;
		HalfMapSum.gppf    += ptr->gppf;
		HalfMapSum.gppfN   += ptr->gppfN;
		HalfMapSum.smppf   += ptr->smppf;
		HalfMapSum.smppfN  += ptr->smppfN;
		HalfMapSum.sgppf   += ptr->sgppf;
		HalfMapSum.sgppfN  += ptr->sgppfN;
		HalfMapSum.genomeN += ptr->genomeN;
	}
	
	fprintf(out,"\nSUM  %9u  %9u  %9u  %9u  %9u  %9u  %9u  %9u  %9u  %9u %9u\n",
					 HalfMapSum.pppf, HalfMapSum.mNppf,
					 HalfMapSum.mppf, HalfMapSum.mppfN,
					 HalfMapSum.gppf, HalfMapSum.gppfN,
					 HalfMapSum.smppf, HalfMapSum.smppfN,
					 HalfMapSum.sgppf, HalfMapSum.sgppfN, HalfMapSum.genomeN);
	
	const unsigned int MappedHalf = HalfMapSum.pppf + HalfMapSum.mNppf + HalfMapSum.mppf + HalfMapSum.mppfN\
	                              + HalfMapSum.gppf + HalfMapSum.gppfN + HalfMapSum.smppf + HalfMapSum.smppfN\
	                              + HalfMapSum.sgppf + HalfMapSum.sgppfN;
	
	fprintf(out,
		"\n|--------------------------- CHIMERA ----------------------------|\n"
	  "       pppf      mNppf       mppf       gppf      smppf      sgppf\n"
	);
	fprintf(out, "  %9u  %9u  %9u  %9u  %9u  %9u\n", cnts->Chimera.pppf, cnts->Chimera.mNppf, cnts->Chimera.mppf, cnts->Chimera.gppf, cnts->Chimera.smppf, cnts->Chimera.sgppf);
	const unsigned int MappedChimera = cnts->Chimera.pppf + cnts->Chimera.mNppf + cnts->Chimera.mppf + cnts->Chimera.mppfN\
	                                 + cnts->Chimera.gppf + cnts->Chimera.gppfN + cnts->Chimera.smppf + cnts->Chimera.smppfN\
	                                 + cnts->Chimera.sgppf + cnts->Chimera.sgppfN;
	
	fprintf(out, "\nUnmapped\n  %9u and %9u encoded\n", cnts->upf, cnts->upf_encoded);
	const unsigned int NotKept = HalfMapSum.genomeN;	
	const float invtot = 100.0f/(float) (MappedOK+MappedTooFar+MappedHalf+MappedChimera+cnts->upf+cnts->upf_encoded\
	                        +PairSum.genomeN+TooFarSum.genomeN+HalfMapSum.genomeN);
	fprintf(out, "\nTOTAL:\n"
				 "  Mapped OK     : %9u   %5.2f%%\n"
				 "  Mapped too far: %9u   %5.2f%%\n"
				 "  Mapped half   : %9u   %5.2f%%\n"
				 "  Chimera       : %9u   %5.2f%%\n"
				 "  Unmapped      : %9u   %5.2f%%\n"
				 "  GenomeN       : %9u   %5.2f%%\n"
				 "--------------------------------------\n",
				MappedOK, (float)MappedOK *invtot,
				MappedTooFar, (float)MappedTooFar*invtot,
				MappedHalf, (float) MappedHalf*invtot,
				MappedChimera, (float) MappedChimera*invtot,
				cnts->upf+cnts->upf_encoded+NotKept, (float) (cnts->upf+cnts->upf_encoded+NotKept)*invtot,
				PairSum.genomeN+TooFarSum.genomeN,  (float)(PairSum.genomeN+TooFarSum.genomeN)*invtot
				);
	const unsigned long Count = MappedOK+MappedTooFar+MappedHalf+MappedChimera+cnts->upf+cnts->upf_encoded+NotKept+PairSum.genomeN+TooFarSum.genomeN;
	fprintf(out, "               %12lu pairs\n\n", Count);
	const unsigned int seconds = elapsed % 60;
	const unsigned int minutes = (elapsed / 60) % 60;
	const unsigned int hours = elapsed / 3600;
	const double elapsetime = (double) ((1.0/60.0)*elapsed);
	const double Throughput = (double) Count/elapsetime;
	if (hours == 0)
	{
		if (minutes == 0)
			fprintf(out, "Total Elapsed time : %02u seconds\nThroughput : %.02lf pairs per minute\n\n", seconds, Throughput);
		else
			fprintf(out, "Total Elapsed time : %02u:%02u\nThroughput : %.02lf pairs per minute\n\n", minutes, seconds, Throughput);
	}
	else
		fprintf(out, "Total Elapsed time : %02u:%02u:%02u\nThroughput : %.02lf pairs per minute\n\n", hours, minutes, seconds, Throughput);
}
//---------------------------------------------------------------

/* vim: tabstop=2 shiftwidth=2
 */
