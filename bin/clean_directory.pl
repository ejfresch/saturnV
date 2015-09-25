#!/usr/bin/perl

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


