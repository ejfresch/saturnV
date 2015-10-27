#!/bin/bash


awk 'BEGIN{FS=OFS="\t"}{$1="";sub("\t","")}1' binary_matrix.tsv > binary_formated.tsv

transpose.r
sed -i -e "1d" matrix_transposed.tsv

sed 's/ //g;s/^./>/;s/"/\n/g' matrix_transposed.tsv > binary_matrix.fasta

raxmlHPC-PTHREADS-AVX -f a -p 12345 -s binary_matrix.fasta -x 12345 -# "$2" -m BINGAMMA -T "$1" -n phylogeny

