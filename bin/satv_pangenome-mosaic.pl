#!/usr/bin/perl

use Parallel::ForkManager;
use Getopt::Long;

$n_cpu=1;
$dir_genomes="";
$out_dir="";

GetOptions ("d=s" => \$dir_genomes,"r=s" => \$ref,"out=s" => \$out_dir,"c=s"   => \$n_cpu) or die("::usage: $0 -d <dir_genomes> -r <ref_genome> -out <out_dir> -c <n_cpu>\n");

if(($dir_genomes eq "") or ($out_dir eq "")){
    print "::usage: $0 -d <dir_genomes> -r <ref_genome> -out <out_dir> -c <n_cpu>\n";
    exit();

}


mkdir($out_dir);


$manager = new Parallel::ForkManager($n_cpu);

$date=`date "+%Y-%m-%d %H:%M:%S"`;
print "::Starting analysis at $date";


chomp($dir_genomes);
print "::DIR_GENOMES=".$dir_genomes."\n";


if((!(-e "${dir_genomes}/$ref")) and (!(-e $ref))){
    print "::Sorry! I cannot find $ref\n";
    exit();
}


@genomes=`ls -1 $dir_genomes`;
chomp(@genomes);

$ref_abs=$ref;
$ref=`basename $ref_abs`;
chomp($ref);
$root_ref=$ref;
$root_ref=~s/.fasta//;
$root_ref=~s/_/-/g;

if(-e "${dir_genomes}/$ref"){
    print "::trimming the name of scaffolds -- $ref\n";
    `satv_trim-name-scaffolds.pl ${dir_genomes}/$ref ${out_dir}/$ref $root_ref`;
    
}else{
    print "::trimming the name of scaffolds -- $ref\n";
    `satv_trim-name-scaffolds.pl ${ref_abs} ${out_dir}/$ref $root_ref`;


}



foreach $genome (@genomes){

    if($genome eq $ref){next;}
    #print $genome."\n";
    $root=$genome;
    $root=~s/.fasta//;
    $root=~s/_/-/g;

    
    if(-e "${out_dir}/${root_ref}_vs_${root}.mum"){
    	print "--I already have already made the comparison -- $ref vs $genome\n";	
    	next;

    }

    

 $manager->start and next;
    
       
    
    
    print "::trimming the name of scaffolds -- $genome\n";
    `satv_trim-name-scaffolds.pl ${dir_genomes}/$genome ${out_dir}/${genome} $root`;
    
    print "::launching mummer -- $ref vs $genome\n";
    
    
    $cmd="mummer -b -c -L ${out_dir}/$ref ${out_dir}/$genome 1>${out_dir}/${root_ref}_vs_${root}.mum 2>>${out_dir}/logs_mummer.txt";
    system($cmd);

  $manager->finish;
}


$manager->wait_all_children;


$date=`date "+%Y-%m-%d %H:%M:%S"`;
print "::Finished! -- $date\n";

