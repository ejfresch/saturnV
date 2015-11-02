#!/usr/bin/perl

use Parallel::ForkManager;
use Getopt::Long;

$dir_genomes="";
$n_cpu=1;
$knife="";

GetOptions ("d=s" => \$dir_genomes,"c=s"   => \$n_cpu ,"k=s"   => \$knife) or die("::usage: $0 -d <dir_genomes> -c <n_cpu> -k <expr>\n");

if(($dir_genomes eq "") or ((!(-e $dir_genomes)) or (!(-d $dir_genomes)))){
    
    if(!(-d $dir_genomes)){
    print "::$dir_genomes is not a directory!\n";
    } 

    print "::usage: $0 -d <dir_genomes> -c <n_cpu> -k <expression>\n";
    exit();
}



$manager = new Parallel::ForkManager($n_cpu);

$date=`date "+%Y-%m-%d %H:%M:%S"`;
print "::Starting analysis at $date";


chomp($dir_genomes);
print "::DIR_GENOMES=".$dir_genomes."\n";

@genomes=`ls -1 $dir_genomes`;
chomp(@genomes);

foreach $genome (@genomes){

    $root=$genome;
    
    if($knife ne ""){
    
    @explode_root=split(/$knife/,$root);
    $root=$explode_root[0];
    }


    $root=~s/.fasta//;
    $root=~s/_/-/g;

    $length_root=length($root);
    if($length_root>15){

    $root=substr($root,0,15);
    print "::[WARGNING] long file name (${genome}). \n--Prokka does not like long names, so I have no choice but cutting your file name. New name: (${root}). You may consider using the option -k to control where the name is cut.\n";

    }
    
    if(-e "${root}/${root}.gff"){
    	print "--I already have the .gff for genome ${root}\n";	
    	next;

    }

 $manager->start and next;
    
       
    
    
    print "::trimming the name of scaffolds -- ${root}\n";
    `satv_trim-name-scaffolds.pl ${dir_genomes}/$genome ${root}.fasta $root`;
    
    print "::launching prokka for genome ${root}\n";
    
    
    `prokka --force --cpus 1 --outdir $root --prefix $root --locustag $root --quiet ${root}.fasta`;
    

  $manager->finish;
}


$manager->wait_all_children;


$date=`date "+%Y-%m-%d %H:%M:%S"`;
print "::Finished! -- $date\n";
print "::Bye\n";

