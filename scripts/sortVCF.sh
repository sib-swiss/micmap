#!/bin/bash

grep "^#" ${1} > ${1/.vcf/}_sorted.vcf
grep -v "^#" ${1} | sort -n -k 2 >> ${1/.vcf/}_sorted.vcf

 
