# Singularity container

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

## Retrieve container from GitHub

```bash
cd /tmp
wget https://github.com/sib-swiss/micmap/releases/download/MicMap_2_20200606/MicMap.sif
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
singularity exec -B "/var/lib/hugetlbfs/global,/data6,/scratch/chris" /tmp/MicMap.sif /bin/bash
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
