#!/usr/bin/perl


print "::Creating links for *.faa files\n";
`ln */*.faa .`;

print "::Creating genomes_to_analyze.txt\n";
`ls *.faa > genomes_to_analyze.txt`;

print "::Bye\n";

