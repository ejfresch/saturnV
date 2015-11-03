#!/usr/bin/perl

use Getopt::Long;

$file_in="";
$dir_out="";
$file_ref="";


GetOptions ("in=s" => \$file_in,"ref=s" => \$file_ref,"out=s" => \$dir_out) or die("::usage: $0 -in <fasta_file> -ref <fasta_ref> -out <out_dir>\n");



if(($file_in eq "") or ($file_ref eq "")or ($dir_out eq "")){
    print "::usage: $0 -in <fasta_file> -ref <fasta_ref> -out <out_dir>\n";
    exit();
}

if(!(-e $file_in)){
    print "::input file \"${file_in}\" does not exists!\n";
    exit();

}

if(!(-e $file_ref)){
    print "::reference genome file \"${file_ref}\" does not exists!\n";
    exit();

}

mkdir($dir_out);

#copy the files 
`cp $file_ref $dir_out`;
`cp $file_in $dir_out`;

$name_file_contigs=`basename $file_in`;
chomp($name_file_contigs);

$name_file_ref=`basename $file_ref`;
chomp($name_file_ref);




chdir($dir_out);
`mv $name_file_ref reference.fasta`;
`mv $name_file_contigs contigs.fasta`;


$cmd="python ~/sw/contiguator/CONTIGuator.py  -r reference.fasta -c contigs.fasta";

system($cmd);

$name_strain=$name_file_contigs;
$name_strain=~s/.fasta//;

`mv PseudoContig.fsa ${name_strain}_closed.fasta`;

