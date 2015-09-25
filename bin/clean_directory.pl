#!/usr/bin/perl

print "::removing Phylogeny\n";
`rm -rf Phylogeny`;
print "::removing prot_files\n";
`rm -rf prot_files`;
print "::removing matrix\n";
`rm -rf matrix`;
print "::removing blasted files\n";
`rm -rf *_blasted*`;
