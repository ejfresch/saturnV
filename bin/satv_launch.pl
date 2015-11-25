#!/usr/bin/perl

use Getopt::Long;
use strict;
#::usage $0 -d <dir_genomes> -c <n_cpu>

my $dir_genomes="";
my $n_cpu=1;
my $fast=0;
my $prokka=1;
my $raxml=0;
my $start_run;
my $end_run;
my $run_time;
my $version = '1.0';
my $date;
my $cmd;
my $identity;
my $identity_usearch=0.6;
my $bootstrap_num=100;
my $bootstrap=1;
my $knife="";

GetOptions ("d=s" => \$dir_genomes,"c=s"   => \$n_cpu, "f=s"   => \$fast,"k=s"   => \$knife, "p=s"   => \$prokka, "r=s"   => \$raxml, "i=s"   => \$identity_usearch, "B=s"   => \$bootstrap_num) or die("::usage: $0 -d <dir_genomes> -c <n_cpu> -f <[0|1]> -k <expression>\n[ERROR] launch failed! Please check the parameters!\n");

if($dir_genomes eq ""){
    print "::usage: $0 -d <dir_genomes> -c <n_cpu> -f <[0|1]> -k <expression>\n[ERROR] launch failed! Please check the parameters!\n";
    exit();
}



$date=`date "+%Y-%m-%d %H:%M:%S"`;
# To know the start time of the script
BEGIN { our $start_run = time(); }

print "::3..2..1..and...lift off -- $date";

#prokka annotation -- first stage
    print "::First stage\n";
    

    if($prokka eq "0"){

   	 $cmd="satv_prodigal-driver.pl -d ${dir_genomes} -c $n_cpu";
   	 system($cmd);
    }




    if($prokka eq "1"){

        if($knife ne ""){
    	 $cmd="satv_prokka-driver.pl -d ${dir_genomes} -c $n_cpu -k $knife";
     	 system($cmd);
        
        }
        else{
    	 $cmd="satv_prokka-driver.pl -d ${dir_genomes} -c $n_cpu";
     	 system($cmd);
  

        }

    }

 
    
    if(!(-e "genomes_to_analyze.txt")){

           print "::Interstage -- generating the genomes_to_analyze.txt file\n";
    #I copy all the .faa files and generate the file genomes_to_analyze.txt

        $cmd="cp */*.faa .";
        system($cmd);
        $cmd="ls *.faa > genomes_to_analyze.txt";
        system($cmd);
    }
    else{
        print "::Interstage -- generating the genomes_to_analyze.txt file (skipped. File exists!)\n";

    }
    

    print "::Second stage \n";



    #determining the pangenome --second stage
    if($fast eq "0"){

        $cmd="satv_pangenome4.pl -g genomes_to_analyze.txt -c $n_cpu";
        system($cmd);
    }

     if($fast eq "1"){

        $cmd="satv_pangenome4-fast.pl -g genomes_to_analyze.txt -c $n_cpu -i $identity_usearch";
        system($cmd);
    }



    print "::Interstage -- merging blast results\n";

    #determining the pangenome --second stage
    $cmd='cp situation_iter1.txt situation_all.txt';
    system($cmd);


    $cmd='cat situation_iter2.txt | grep -v "#"  >> situation_all.txt';
    system($cmd);




#mergins the results -- third stage


    print "::Third stage\n";
    
    #determining the pangenome --second stage
    $cmd="satv_merge-data3.pl -i situation_all.txt";
    system($cmd);




	# if the user want a phylogeny generated from the binary matrix
    if($raxml eq "1"){
        system("satv_generate-binary-matrix.pl table_linked3.tsv");
        $cmd="satv_raxml-driver.sh $n_cpu $bootstrap_num";
        system($cmd);

	mkdir 'Phylogeny';
	system("mv RAxML* Phylogeny");
	system("mv binary_matrix.tsv matrix");
	system("mv binary_matrix.fasta matrix");
    }





mkdir 'prot_files';


system("mv *.faa prot_files");
system("mv *.gff prot_files");
system("rm *.udb");






$date=`date "+%Y-%m-%d %H:%M:%S"`;
# To know how many time was done since the script started
$end_run = time();
$run_time = $end_run - our $start_run;
print "::Mission accomplished by SaturnV version $version in $run_time seconds -- $date\n";






