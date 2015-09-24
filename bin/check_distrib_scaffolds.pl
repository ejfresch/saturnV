#!/usr/bin/perl

#convert the file to a numerical matrix

use Getopt::Long;

#::usage $0 -d <dir_genomes> -c <n_cpu>

$dir_genomes="";


GetOptions ("d=s" => \$dir_genomes) or die("::usage: $0 -d <dir_genomes>\n");



if(($dir_genomes eq "")){
    print "::usage $0 -d <dir_genomes>\n";
    exit();
}




