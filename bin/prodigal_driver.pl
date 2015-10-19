#!/usr/bin/perl

use Parallel::ForkManager;
use Getopt::Long;


GetOptions ("d=s" => \$dir_genomes,"c=s"   => \$n_cpu) or die("::usage: $0 -d <dir_genomes> -c <n_cpu>\n");

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
	mkdir "$root";
    if(-e "${root}/${root}.faa"){
    	print "--I already have the .gff for genome $genome\n";	
    	next;

    }

 $manager->start and next;
    
       
    
    
    print "::trimming the name of scaffolds -- $genome\n";
    `trim_name_scaffolds.pl ${dir_genomes}/$genome ${genome} $root`;
    
    print "::launching prodigal for genome $genome\n";    
    `prodigal -i $genome -a ${root}/${root}.faa -f gff -o ${root}/${root}.gff >> log_file 2>&1`;
    



	open($faa_clean, '>', "${root}/${root}_clean.faa");
	open(IN1,"${root}/${root}.faa")||die "I cannot open ${root}/${root}.faa";
	$i=1;
   	 while($line1=<IN1>){
     	   chomp($line1);

		if ($line1 =~ /^\>/)
				{
				print $faa_clean ">${root}_seq$i\n";
				$i++;
				}
		else
		{
		print $faa_clean "$line1\n";
		}

	}
	close(IN1);
	close("$faa_clean");
	system("mv ${root}/${root}_clean.faa ${root}/${root}.faa");

  $manager->finish;
}


$manager->wait_all_children;


$date=`date "+%Y-%m-%d %H:%M:%S"`;
print "::Finished! -- $date\n";
print "::Bye\n";

