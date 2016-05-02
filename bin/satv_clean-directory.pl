#!/usr/bin/perl

use Getopt::Long;


$method="";

GetOptions("m=s" => \$method,) or die("::usage: $0 -m <method>\n");



print "::removing Phylogeny\n";
`rm -rf Phylogeny`;
print "::removing prot_files\n";
`rm -rf prot_files`;
print "::removing matrix\n";
`rm -rf matrix`;
print "::removing blasted files\n";
`rm -rf *_blasted*`;

print "::removing situation (iter1, iter2 and all)\n";
`rm -rf situation_iter1.txt`;
`rm -rf situation_iter2.txt`;
`rm -rf situation_all.txt`;

print "::removing new_sequences_iter1.faa\n";
`rm -rf new_sequences_iter1.faa`;

print "::removing genomes_to_analyze.txt\n";
`rm -rf genomes_to_analyze.txt`;


if($method eq "strict"){
 print "::removing fasta_subject files\n";
`rm -rf fasta_subject*`;   
 print "::removing screen_paralogs files\n";
`rm -rf screen_paralogs.txt`;   
 print "::removing commands-* files\n";
`rm -rf commands-*`;
 print "::removing usearch-resolve-line-* files\n";
`rm -rf usearch_resolve_line-*`;
 print "::removing pre-table_linked5_strict.tsv\n";
`rm -rf pre-table_linked5_strict.tsv`;   
 print "::removing table_linked5_strict.tsv\n";
`rm -rf table_linked5_strict.tsv`;   

  

}
