#!/usr/bin/perl 

use Getopt::Long;

$input_dir="";
$ref_fasta="";
$out_dir="";


GetOptions ("d=s" => \$input_dir,"r=s" => \$ref_fasta,"out=s" => \$out_dir) or die("::usage: $0 -d <input_dir> -r <ref_fasta> -out <out_dir>\n");

if(($input_dir eq "") or ($ref_fasta eq "") or ($out_dir eq "")){
    print "::usage: $0 -d <input_dir> -r <ref_fasta> -out <out_dir>\n";
    exit();
}



#I exec pangenome_mosaic
print "::Comparing genomes to the reference\n";
$cmd="satv_pangenome-mosaic.pl -d $input_dir -r $ref_fasta -out $out_dir";
system($cmd);

#I reconstruct the mosaic
print "::Generating binary fasta files\n";


$new_ref=`basename $ref_fasta`;
chomp($new_ref);

@files=`ls ${out_dir} |grep "\.mum\$"`;
chomp(@files);

foreach $file (@files){
    print "--generating binary fasta file for $file\n";
    $cmd="satv_reconstruct-mosaic.pl -mum ${out_dir}/$file -r ${out_dir}/${new_ref} -out $out_dir";
    system($cmd);
}


