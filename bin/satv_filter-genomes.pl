#!/usr/bin/perl

use Getopt::Long;


$min=1;
$max="";

GetOptions ("in=s" => \$dir_genomes,"out=s"   => \$out_dir, "min=s"   => \$min, "max=s"   => \$max) or die("::usage: $0 -in <dir_genomes> -out <dir_genomes_filtered> -min <min_num_contigs/scaffolds> -max <max_num_contigs/scaffolds>\n");



if(-e $out_dir){
    print "\n::Be careful! $out_dir already exists. I will add the filtered genomes to this directory, but you have to deal with the genomes that are already present in this directory.\n\n";
}

mkdir($out_dir);

if(($dir_genomes eq "") or ($out_dir eq "")){
    print "::usage: $0 -in <dir_genomes> -out <dir_genomes_filtered> -min <min_num_contigs/scaffolds> -max <max_num_contigs/scaffolds>\n";
    exit();
}

open(IN,"<${dir_genomes}/distribution_fragments_genomes.csv") or die "::I cannot find the distribution_fragments_genomes.csv file in the $dir_genomes folder\n";
$comment=<IN>;
while($line=<IN>){
    chomp($line);
    @data=split(/,/,$line);
    $fragments=$data[1];
    $file_name=$data[0];
    
    if($fragments>=$min){
        if($max eq "")
            {
            print "::${file_name} is ok ($fragments contigs/scaffolds)\n";
                $cmd="cp ${dir_genomes}/${file_name} ${out_dir}/${file_name}";            
                system($cmd);
            }
        if($fragments<=$max)
            {
            print "::${file_name} is ok ($fragments contigs/scaffolds)\n";
            $cmd="cp ${dir_genomes}/${file_name} ${out_dir}/${file_name}";            
            system($cmd);

            }


    
    }


}
  
close(IN);
