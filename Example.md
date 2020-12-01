# Example

[Project on GitHub](https://github.com/sib-swiss/micmap)

[This file](https://github.com/sib-swiss/micmap/blob/master/Example.md)

<a href="http://www.youtube.com/watch?feature=player_embedded&v=mmqRlyJq0_c"
target="_blank"><img src="http://img.youtube.com/vi/mmqRlyJq0_c/0.jpg"
alt="MicMap Setup" width="240" height="180" border="10" /></a>

https://www.youtube.com/watch?v=mmqRlyJq0_c

#  Background info - singularity container

## Initial empty container

```bash
singularity build --sandbox /tmp/MicMap docker://fedora:latest
```

## Add tools from outside the container

```bash
cd src/github/sib/micmap
rm -rf build
mkdir build
cd build
source /opt/intel/bin/compilervars.sh intel64
CC=icc cmake ../ -DGTL_TESTS_PATH=/tmp/GTL -DUSE_RPATH=0 -DCMAKE_INSTALL_RPATH=\$ORIGIN/../lib
ccmake ..
cmake ..
make
make DESTDIR=/tmp/MicMap install
```

## Add tools or tweak the content from within the container

```bash
singularity exec --writable /tmp/MicMap /bin/bash
```

## Generate `.sif` file

```bash
sudo singularity build /tmp/MicMap.sif /tmp/MicMap
```

## Vim tools for singularity config files

```bash
git clone https://github.com/singularityhub/singularity.lang.git
cd singularity.lang/vim
make install
```

#  Getting started - getting the parts

## Retrieve container

```bash
rm -rf /tmp/demo
mkdir /tmp/demo
cd /tmp/demo
wget https://bix.unil.ch/ref/MicMap-2.2.sif
```

## Retrieve reference data

```bash
cd /tmp/demo
wget https://bix.unil.ch/ref/b38.tar.xz
tar xvfJ b38.tar.xz
rm b38.tar.xz
```

## hugepage setup

This needs to be specified in your GRUB config as kernel parameters

```bash
transparent_hugepage=never default_hugepagesz=2M hugepagesz=2M hugepages=1024 hugepagesz=1G hugepages=48
```

This needs to be run as root to mount the hugepage filesystem

```bash
hugeadm --create-global-mounts
chmod 1777 /var/lib/hugetlbfs/global/pagesize-*
```

## Start session in container

```bash
tmux attach
screen -UAa
singularity exec -B "/var/lib/hugetlbfs/global" /tmp/demo/MicMap-2.2.sif /bin/bash
```

# Very simple example

<a href="http://www.youtube.com/watch?feature=player_embedded&v=qMcAJelsK1s"
target="_blank"><img src="http://img.youtube.com/vi/qMcAJelsK1s/0.jpg"
alt="MicMap Simple Example" width="240" height="180" border="10" /></a>

This is fabricated data using `grep NC_000021_593_9 /tmp/demo/b38/mm_RNAseq_clean_b38.txt >/tmp/demo/NC_000021_593_9.txt`

The headers mimic real Illumina data but this is not mandatory.

```bash
singularity exec -B "/var/lib/hugetlbfs/global" /tmp/demo/MicMap-2.2.sif /bin/bash

cat > /tmp/demo/tst_R1.fq <<EOF
@ST-E00137:143:H2JFTALXX:8:1101:1609:2205p
CCCAAGGCCATGCCAGCTCAAAGGTGGTGGGATTCCTTGGCGATGTCCCTGAAGCTCTTATTTCCCACTTCGAAAGAACGAAAAGGGTCTGAATTTTATTTACACATGGATTTAGTTTGCT
+
JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
@ST-E00137:143:H2JFTALXX:8:1101:1649:2205N
GAAGCTCTTATNTCCCACTTCGAAAGAACGAAAAGGGTCTGAATTTTATTTACACATGGATTTAGTTTGCTAATTCTTAGAGATTCCTTTTTTACTGCTAGATCTTTTTTAGTAAACACAA
+
JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
@ST-E00137:143:H2JFTALXX:8:1101:1771:2205m
TGTGTCTGCCAGGAACTGCACTGCATGCGGCCCACAGCCCTGCCTCAGAAGCAGCTTTGGTCTCGTTTGCTGAGACCGCTTGCTGTGCCTCATCAATTTCAGGAGGAAACAGTGAGATGCC
+
JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
@ST-E00137:143:H2JFTALXX:8:1101:1791:2205g
GTGGCCACTTGGCTTAGCCGTGCTCATGGTGCAGACGTGGACGCGGAGGCAGGGAGAGGCTCCGTGACGTCCCCAGGGCCCCAGAACGAGGAAGGAGCGGAGTTGGGATTCCAGCCCAGTT
+
JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
@ST-E00137:143:H2JFTALXX:8:1101:1832:2205h
GTTTTATGCAGTCAGTGCCAGAATTGGATTACATCACCTGGGGGAGGCGCCTTCAGGACCAGGAAATGGTACATGGAGCTTTGAGGTGAAAACTCATCTTAGAAACCCACAAAAAGTCATT
+
JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
@ST-E00137:143:H2JFTALXX:8:1101:1873:2205c
TCCTTGGCGATGTCCCTGAAGCTCTTATTTCCCACTTCGAAAGAACGAAAAGGGTCTGAATTTTATTTACACATGGATTTAGTTTGCTAATTCTTAGAGATTCCTTTTTTACTGCTAGATC
+
JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
@ST-E00137:143:H2JFTALXX:8:1101:1893:2205u
GTGTGTGTGTGTGTGTGTGTATACACACACATACACACACATACACACACATACACACACACATACAGTGTGTACACATACACACACATACACACACATAACACATACACACACATACACA
+
JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
EOF

cat > /tmp/demo/tst_R2.fq <<EOF
@ST-E00137:143:H2JFTALXX:8:1101:1609:2205p
ATTAAATGAGGACTTTTTATAAAGGCTAAGTTTGACTCTGGTTAAAAACTTAATCTGGAAGGCAATGCCATATTTGTATAGAACAAGCACACTTAAAATTACACTTTAATATCTATCTCAA
+
JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
@ST-E00137:143:H2JFTALXX:8:1101:1649:2205N
TATGTGCAATGTCTTACAAGANGATCTATTTTGAAATTAAATGAGGACTTTTTATAAAGGCTAAGTTTGACTCTGGTTAAAAACTTAATCTGGAAGGCAATGCCATATTTGTATAGAACAA
+
JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
@ST-E00137:143:H2JFTALXX:8:1101:1771:2205m
CAAAGCCATCCTGGGGCGTCAACGGCAGAATCCTCACCTACTCCAGTTCACCCTTCAAAAAATAAACCCACATGGTCTCACATATCCATTCCCTTCTTGCAGCCTAACAGCCCTCATGGTT
+
JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
@ST-E00137:143:H2JFTALXX:8:1101:1791:2205g
GACAGCCCCTCTCCACGGAGAGGTACCTGTCACGAGTACTCTGATGCAAGGCCACTGAAGGGCGTTAGGCTTGGGTTGCCCTCTCATCTCAGCCACCTGTTTGGTAACAGAAGGAATCACT
+
JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
@ST-E00137:143:H2JFTALXX:8:1101:1832:2205h
TGTGTATGTGTGTGTATGTGTTATGTGTGTGTATGTGTGTGTATGTGTACACACTGTATGTGTGTGTGTATGTGTGTGTATGTGTGTGTATGTGTGTGTATACACACACACACACACACAC
+
JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
@ST-E00137:143:H2JFTALXX:8:1101:1873:2205c
GATCAGCCAGATTCTCTCGTATGTTTGTTCAGTGCAGACACACAGACGATGGTGGCGCGCTATTGAGAGAAAACAGCTGAAAGTAAAGTAATAATTTTGAGTTTTGCTACTTTTTGGTTAT
+
JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
@ST-E00137:143:H2JFTALXX:8:1101:1893:2205u
TGTGTATGTGTGTGTATGTGTTATGTGTGTGTATGTGTGTGTATGTGTACACACTGTATGTGTGTGTGTATGTGTGTGTATGTGTGTGTATGTGTGTGTATACACACACACACACACACAC
+
JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
EOF

rm /tmp/demo/tst_b38_*
rm -f /var/lib/hugetlbfs/global/pagesize-1GB/*
MicMap -g /tmp/demo/b38/b38.bin -p /tmp/demo/b38/tbl_18nt/start0. -s /var/lib/hugetlbfs/global/pagesize-1GB
MicMap -g /var/lib/hugetlbfs/global/pagesize-1GB/b38.bin -p /var/lib/hugetlbfs/global/pagesize-1GB/start0. -d 1 -a 1 -W 1 -w 1 -1 /tmp/demo/tst_R1.fq -2 /tmp/demo/tst_R2.fq -r /tmp/demo/tst_b38
ls -alth /tmp/demo/|head -20
cat /tmp/demo/tst_b38_log.txt

GTLdecompress -g /tmp/demo/b38/b38.bin -r /tmp/demo/tst_b38_A_chr21a.gtl -p -o SAM
GTLdecompress -g /tmp/demo/b38/b38.bin -r /tmp/demo/tst_b38_A_chr21a.gtl -n -o SAM
GTLdecompress -g /tmp/demo/b38/b38.bin -r /tmp/demo/tst_b38_A_chr21a.gtl -m -o SAM
GTLdecompress -g /tmp/demo/b38/b38.bin -r /tmp/demo/tst_b38_A_chr21a.gtl -a -o SAM
GTLdecompress -g /tmp/demo/b38/b38.bin -r /tmp/demo/tst_b38_A_C_chr21.gtl -c -o SAM
GTLdecompress -g /tmp/demo/b38/b38.bin -r /tmp/demo/tst_b38_A_HM_chr21.gtl -h -o SAM
GTLdecompress -g /tmp/demo/b38/b38.bin -r /tmp/demo/tst_b38_A_UM_GT.gtl -u -o SAM

sortCompressCall -s /tmp/demo -d /tmp/demo -t /tmp/demo -C -D tst_b38
ls -alth /tmp/demo/|head -20
ls /tmp/demo/tst_b38

samtools ADVIEW -a /tmp/demo/b38/GB -G /tmp/demo/b38/b38.bin -1 /tmp/demo/tst_b38 -c 21 -p 42298101

GTLdecompress -g /tmp/demo/b38/b38.bin -i /tmp/demo/tst_b38 -C 21 -P 42287101..42308125 -p -o ADNIview
GTLdecompress -g /tmp/demo/b38/b38.bin -i /tmp/demo/tst_b38 -C 21 -P 42287101..42308125 -n -o ADNIview
GTLdecompress -g /tmp/demo/b38/b38.bin -i /tmp/demo/tst_b38 -C 21 -P 42287101..42308125 -m -o ADNIview
GTLdecompress -g /tmp/demo/b38/b38.bin -i /tmp/demo/tst_b38 -C 21 -P 42287101..42308125 -a -o ADNIview
```

# Grab example from SRA

[sra](https://www.ncbi.nlm.nih.gov/sra)

Search details : `(illumina[All Fields] AND wxs[All Fields]) AND "Homo sapiens"[orgn] AND (cluster_public[prop] AND "biomol dna"[Properties] AND "strategy wxs"[Properties] AND "library layout paired"[Properties] AND "platform illumina"[Properties] AND "filetype fastq"[Properties])`


For example, from biosample [SAMN14857968](https://www.ncbi.nlm.nih.gov/biosample/SAMN14857968) I select 3 samples :

* [SRX8296330](https://www.ncbi.nlm.nih.gov/sra/SRX8296330[accn])
* [SRX8296329](https://www.ncbi.nlm.nih.gov/sra/SRX8296329[accn])
* [SRX8296328](https://www.ncbi.nlm.nih.gov/sra/SRX8296328[accn])

## Retrieve using sra tools

```bash
vdb-config --interactive
mkdir /tmp/Ex
cd /tmp/Ex/
fastq-dump --split-files -F -I -L 5 SRR11742849 SRR11742850 SRR11742851
pigz --best *.fastq&
mkdir FQ
mv SRR117428* FQ
```

# Run MicMap to perform the mapping

```bash
MicMapExome2bam -d /tmp/Ex -r b38 -G P13normal -S P13normal -L SRR11742849 /tmp/Ex/FQ/SRR11742849_?.fastq.gz P13_normal
MicMapExome2bam -d /tmp/Ex -r b38 -G P13lung -S P13lung -L SRR11742850 /tmp/Ex/FQ/SRR11742850_?.fastq.gz P13_lung
MicMapExome2bam -d /tmp/Ex -r b38 -G P13cfDNA -S P13cfDNA -L SRR11742851 /tmp/Ex/FQ/SRR11742851_?.fastq.gz P13_cfDNA
```

# Perform somatic variant detection

From GATK best practice

https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38/

```bash
git clone git@github.com:AstraZeneca-NGS/reference_data.git
grep -v '^chr[0-9XYUn]*_' ~/src/github/reference_data/hg38/bed/Exome-NGv3.bed >~/src/github/reference_data/hg38/bed/Exome-NGv3_fix.bed

wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi

wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi

wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz.tbi

ExomeDiffCall -d /tmp/Ex P13_normal P13_lung ~/src/github/reference_data/hg38/bed/Exome-NGv3_fix.bed

ExomeDiffCall -d /tmp/Ex P13_normal P13_cfDNA ~/src/github/reference_data/hg38/bed/Exome-NGv3_fix.bed
```


# Visualization

```bash
cd /tmp/Ex

samtools ADVIEW -a /data6/b38/GB -G /data6/b38/b38.bin \
  -1 P13_cfDNA_b38 -2 P13_cfDNA_b38.bam -3 P13_normal_b38 -4 P13_normal_b38.bam \
  -c 1 -p 952430 \
  -s <(genRGBsnp <(zgrep -w PASS P13_cfDNA_vs_P13_normal_M2.vcf.gz) \
                 <(grep -w PASS P13_cfDNA_vs_P13_normal.vcf) \
                 <(grep -w PASS P13_cfDNA_vs_P13_normal_old.vcf))

samtools ADVIEW -a /data6/b38/GB -G /data6/b38/b38.bin \
  -1 P13_lung_b38 -2 P13_lung_b38.bam -3 P13_normal_b38 -4 P13_normal_b38.bam \
  -c 1 -p 952430 \
  -s <(genRGBsnp <(zgrep -w PASS P13_lung_vs_P13_normal_M2.vcf.gz) \
                 <(grep -w PASS P13_lung_vs_P13_normal.vcf) \
                 <(grep -w PASS P13_lung_vs_P13_normal_old.vcf))

samtools ADVIEW -a /data6/b38/GB -G /data6/b38/b38.bin \
  -1 P13_lung_b38 -2 P13_lung_b38.bam -3 P13_normal_b38 -4 P13_normal_b38.bam \
  -c 1 -p 952430 \
  -s <(genRGBsnp <(zgrep -w PASS P13_lung_vs_P13_normal_M2.vcf.gz) \
                 <(grep -w PASS P13_lung_vs_P13_normal.vcf) \
                 <(grep -w PASS P13_lung_vs_P13_normal_old.vcf) \
                 | grep -w '10[236]$')

samtools ADVIEW -a /data6/b38/GB -G /data6/b38/b38.bin \
  -1 P13_cfDNA_b38 -2 P13_lung_b38 -3 P13_normal_b38 -4 P13_cfDNA_b38.bam -5 P13_lung_b38.bam -6 P13_normal_b38.bam \
  -c 1 -p 952430 \
  -s <(genRGBsnp <(zgrep -w PASS P13_cfDNA_vs_P13_normal_M2.vcf.gz) \
                 <(grep -w PASS P13_cfDNA_vs_P13_normal.vcf) \
                 <(grep -w PASS P13_cfDNA_vs_P13_normal_old.vcf) \
                 | grep -w '10[26]$')
```

```
         1         2         3         4         5         6         7         8         9        10        11        12
123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
```
