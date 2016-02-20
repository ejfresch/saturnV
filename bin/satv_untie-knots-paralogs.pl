#!/usr/bin/perl 

use Getopt::Long;
use strict;
use Parallel::ForkManager;


my $input_file="";
my $out_file="";
my $n_cpu=1;
my $identity_orthologs=90;
my $identity_paralogs=90;
my $genome_list="";


GetOptions ("g=s" => \$genome_list,"in=s" => \$input_file,"out=s" => \$out_file,"c=s"   => \$n_cpu,"i=s"   => \$identity_orthologs,"ip=s"   => \$identity_paralogs) or die("::usage: $0 -g <genomes_list> -in <input_file> -c <n_cpu> -i <identity_usearch> -ip <identity_paralogs>\n[ERROR] launch failed! Please check the parameters!\n");

if(($input_file eq "") or ($out_file eq "") or ($genome_list eq "")){

    print "::usage: $0 -g <genomes_list> -in <input_file> -c <n_cpu> -i <identity_usearch> -ip <identity_paralogs> -out <out_file>\n";
    exit();
}

my $manager = new Parallel::ForkManager($n_cpu);



open(IN_FILE,"<$input_file");
open(OUT_FILE,">$out_file");
my $line=<IN_FILE>;
chomp($line);
print OUT_FILE $line."\n";
$line=~s/#//;
my @genomes=split(/\t/,$line);

my $count_overall=1;

while($line=<IN_FILE>){

chomp($line);

    $count_overall++;


    if($line=~/,/){
    

        $manager->start and next;

	    $line=~s/\t/&/g;	
    	print $line."\n";
        my $cmd="satv_resolve-line.pl -g $genome_list -out $out_file -id $count_overall -i $identity_orthologs -ip $identity_paralogs -line \"$line\"";
        system($cmd);


        $manager->finish;

     
    }



    else{print OUT_FILE $line."\n";}



}

$manager->wait_all_children;


close(IN_FILE);
close(OUT_FILE);

print "::merging results\n";

my $cmd="cat usearch_untie_knots_paralogs_table_* >> $out_file";
system($cmd);

