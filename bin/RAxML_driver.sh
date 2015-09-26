#!/bin/bash


awk 'BEGIN{FS=OFS="\t"}{$1="";sub("\t","")}1' binary_matrix.tsv > binary_formated.tsv

Rscript /home/avincent/Desktop/script_luca_propre/bin/transpose.r
sed -i -e "1d" matrix_transposed.tsv

sed 's/ //g;s/^./>/;s/"/\n/g' matrix_transposed.tsv > binary_matrix.fasta

raxmlHPC-PTHREADS-AVX -m BINGAMMA -s binary_matrix.fasta -p 12345 -n phylogeny -T $1

raxmlHPC-PTHREADS-AVX -m BINGAMMA -s binary_matrix.fasta -p 12345 -b 12345 -# $2 -n bootstrap -T $1

raxmlHPC-PTHREADS-AVX -f b -m BINGAMMA -t RAxML_bestTree.phylogeny -z RAxML_bootstrap.bootstrap -n tree -T $1


