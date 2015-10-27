#!/usr/bin/perl

use Parallel::ForkManager;
use Getopt::Long;

$dir_genomes="";
$n_cpu=1;

GetOptions ("d=s" => \$dir_genomes,"c=s"   => \$n_cpu) or die("::usage: $0 -d <dir_genomes> -c <n_cpu>\n");

if(($dir_genomes eq "") or ((!(-e $dir_genomes)) or (!(-d $dir_genomes)))){
    
    if(!(-d $dir_genomes)){
    print "::$dir_genomes is not a directory!\n";
    exit();    
    } 

    print "::usage: $0 -d <dir_genomes> -c <n_cpu>\n";
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
    $root=~s/.fasta//;
    $root=~s/_/-/g;
    if(-e "${root}/${root}.gff"){
    	print "--I already have the .gff for genome $genome\n";	
    	next;

    }

 $manager->start and next;
    
       
    
    
    print "::trimming the name of scaffolds -- $genome\n";
    `satv_trim-name-scaffolds.pl ${dir_genomes}/$genome ${genome} $root`;
    
    print "::launching prokka for genome $genome\n";
    
    
    `prokka --force --cpus 1 --outdir $root --prefix $root --locustag $root --quiet $genome`;
    

  $manager->finish;
}


$manager->wait_all_children;


$date=`date "+%Y-%m-%d %H:%M:%S"`;
print "::Finished! -- $date\n";
print "::Bye\n";

