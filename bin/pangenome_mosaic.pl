#!/usr/bin/perl -I /home/lfreschi/tasks/pangenome/script_new/bin

use Parallel::ForkManager;
use Getopt::Long;


GetOptions ("d=s" => \$dir_genomes,"r=s" => \$ref,"out=s" => \$out_dir,"c=s"   => \$n_cpu) or die("::usage: $0 -d <dir_genomes> -r <ref_genome> -out <out_dir> -c <n_cpu>\n");

mkdir($out_dir);


$manager = new Parallel::ForkManager($n_cpu);

$date=`date "+%Y-%m-%d %H:%M:%S"`;
print "::Starting analysis at $date";


chomp($dir_genomes);
print "::DIR_GENOMES=".$dir_genomes."\n";

@genomes=`ls -1 $dir_genomes`;
chomp(@genomes);


$root_ref=$ref;
$root_ref=~s/.fasta//;
$root_ref=~s/_/-/g;

print "::trimming the name of scaffolds -- $ref\n";
`trim_name_scaffolds.pl ${dir_genomes}/$ref $ref $root_ref`;



foreach $genome (@genomes){

    if($genome eq $ref){next;}

    $root=$genome;
    $root=~s/.fasta//;
    $root=~s/_/-/g;

    
    if(-e "${out_dir}/${root_ref}_vs_${root}.mum"){
    	print "--I already have already made the comparison -- $ref vs $genome\n";	
    	next;

    }

    

 $manager->start and next;
    
       
    
    
    print "::trimming the name of scaffolds -- $genome\n";
    `trim_name_scaffolds.pl ${dir_genomes}/$genome ${genome} $root`;
    
    print "::launching mummer -- $ref vs $genome\n";
    
    
    $cmd="mummer -b -c -L $ref $genome 1>${out_dir}/${root_ref}_vs_${root}.mum 2>>logs_mummer.txt";
    system($cmd);

  $manager->finish;
}


$manager->wait_all_children;


$date=`date "+%Y-%m-%d %H:%M:%S"`;
print "::Finished! -- $date\n";
print "::Bye\n";









