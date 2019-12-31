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
#include <stdio.h>
#include "vcf.h"
#include "vcfutils.h"

void usage() {
        puts(
        "NAME\n"
        "    03_vcf - High quality calls for single sample\n"
        "SYNOPSIS\n"
        "    03_vcf vcf_file sample\n"
        "DESCRIPTION\n"
        "    Given a <vcf file>, extract all calls for <sample> and filter for\n"
        "    high quality, homozygous SNPs. This will omit any positions\n"
        "    that are homozygous ref (0/0) or heterozygous.  The exact filter\n"
        "    used is\n"
        "        FI == 1 & GQ > 20 & GT != '0/0'\n"
        "        [NOTE: FI == 1 implies homozygous call]\n"
        "        \n"
        "    The returned format is\n"
        "        chrom pos[0-based]  REF  ALT GQ|DP\n"
        "    and can be used for Marei's personalizer.py\n"
                );
}



int main(int argc, char **argv) {
        if (argc < 2 || argc > 3 ) {
                usage();
                return 1;
        }
        // counters
        int n    = 0;  // total number of records in file
        int nsnp = 0;  // number of SNP records in file
        int nhq  = 0;  // number of SNPs for the single sample passing filters
        int nseq = 0;  // number of sequences
        // filter data for each call
        int nfi_arr = 0;
        int nfi     = 0;
        int *fi     = NULL;
        // quality data for each call
        int ngq_arr = 0;
        int ngq     = 0;
        int *gq     = NULL;
        // coverage data for each call
        int ndp_arr = 0;
        int ndp     = 0;
        int *dp     = NULL;
        // genotype data for each call
        // genotype arrays are twice as large as
        // the other arrays as there are two values for each sample
        int ngt_arr = 0;
        int ngt     = 0;
        int *gt     = NULL;
        
        // open VCF/BCF file
        //    * use '-' for stdin
        //    * bcf_open will open bcf and vcf files
        //    * bcf_open is a macro that expands to hts_open
        //    * returns NULL when file could not be opened
        //    * by default also writes message to stderr if file could not be found
        htsFile * inf = bcf_open(argv[1], "r");
        if (inf == NULL) {
                return EXIT_FAILURE;
        }
        
        // read header
        bcf_hdr_t *hdr = bcf_hdr_read(inf);
        fprintf(stderr, "File %s contains %i samples\n", argv[1], bcf_hdr_nsamples(hdr));
        // report names of all the sequences in the VCF file
        const char **seqnames = NULL;
        // bcf_hdr_seqnames returns a newly allocated array of pointers to the seq names
        // caller has to deallocate the array, but not the seqnames themselves; the number
        // of sequences is stored in the int pointer passed in as the second argument.
        // The id in each record can be used to index into the array to obtain the sequence
        // name
        seqnames = bcf_hdr_seqnames(hdr, &nseq);
        if (seqnames == NULL) {
                goto error1;
        }
        fprintf(stderr, "Sequence names:\n");
        for (int i = 0; i < nseq; i++) {
                // bcf_hdr_id2name is another way to get the name of a sequence
                fprintf(stderr, "  [%2i] %s (bcf_hdr_id2name -> %s)\n", i, seqnames[i],
                       bcf_hdr_id2name(hdr, i));
        }

        // limit the VCF data to the sample name passed in
//         char * const ctmp = (argc == 3) ? argv[2] : NULL;
// 				if (bcf_hdr_set_samples(hdr, ctmp, 0) < 0) {
// 					perror("bcf_hdr_set_samples");
// 				}
				printf("There are %i samples\n", bcf_hdr_nsamples(hdr));

        // struc for storing each record
        bcf1_t *rec = bcf_init();
        if (rec == NULL) {
                goto error2;
        }
        
        while (bcf_read(inf, hdr, rec) == 0) {
								if (rec->rid != 0) continue;
                n++;
                if (bcf_is_snp(rec)) {
                        nsnp++;
                }
//                 else {
//                         continue;
//                 }
                // the bcf_get_format_int32 function does not appear to reallocate
                // the array it returns for each of the samples on each call. Just
                // needs to be freed in the end. First call to bcf_get_format_*
                // takes care of calling bcf_unpack, which fills the `d` member
                // of bcf1_t
//                 nfi = bcf_get_format_int32(hdr, rec, "FI", &fi, &nfi_arr);
                // GQ can be missing (".") in this VCF file; The htslib version
                // used right now does not return a negative value in that case,
                // so we can't check for it.  As it turns out, all homozygous
                // good quality calls for alt allele have GQ values, so it doesn't matter
                // here, but it's important to keep in mind.
//                 ngq = bcf_get_format_int32(hdr, rec, "GQ", &gq, &ngq_arr);
                 ndp = bcf_get_format_int32(hdr, rec, "DP", &dp, &ndp_arr);
                 ngt = bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt_arr);

//                 if (fi[0] == 1 && gq[0] > 20 && gt[0] != 0 && gt[1] != 0) {
                        nhq++;
//                         printf("chr%s\t%i\t%s\t%s\t%i|%i\n", seqnames[rec->rid],
//                                rec->pos,
//                                rec->d.allele[0],
//                                rec->d.allele[bcf_gt_allele(gt[0])],
//                                gq[0], dp[0]);
												printf("chr%s\t%i\t%s\t%s", seqnames[rec->rid],
                               rec->pos,
                               rec->d.allele[0], rec->d.allele[1]);
												

												if (ngt > 0) {
													printf("\t%s/%s, ngt=%i, gt=", rec->d.allele[bcf_gt_allele(gt[0])],
													rec->d.allele[bcf_gt_allele(gt[1])], ngt);
													for (int p=0; p<ngt; p++) printf("%i,", gt[p]);
													printf("->%u/%u",bcf_gt_allele(gt[0]),bcf_gt_allele(gt[1]));
												}
												
												if (ndp > 0)  printf("\tdepth=%i", *dp);
												else printf("\tndp=%i, ndparr=%i", ndp, ndp_arr);
												
												if (bcf_is_snp(rec)) printf(", SNP");
												bcf_get_variant_types(rec);
												printf(", n_var=%i", rec->d.n_var);
												for(int i=0;i<rec->d.n_var;i++) printf(" %i", rec->d.var[i].n);
												
												fputc('\n', stdout);
//                 }
								 
        }
        fprintf(stderr, "Read %i records %i of which were SNPs\n", n, nsnp);
        fprintf(stderr, "%i records for the selected sample were high quality homozygous ALT SNPs in sample %s\n", nhq, argv[2]);
        free(fi);
        free(gq);
        free(gt);
        free(dp);
        free(seqnames);
        bcf_hdr_destroy(hdr);
        bcf_close(inf);
        bcf_destroy(rec);
        return EXIT_SUCCESS;
error2:
        free(seqnames);
error1:
        bcf_close(inf);
        bcf_hdr_destroy(hdr);
        return EXIT_FAILURE;
}

