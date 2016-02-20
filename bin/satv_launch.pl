#!/usr/bin/perl

use Getopt::Long;
use strict;


my $dir_genomes="";
my $n_cpu=1;
my $alg="usearch";
my $method="strict";
my $ann="prodigal";
my $raxml=0;
my $start_run;
my $end_run;
my $run_time;
my $date;
my $cmd;
my $identity_orthologs=90;
my $identity_paralogs=90;
my $bootstrap_num=100;
my $bootstrap=1;
my $knife="";
my $force=0;
my $sim=0;

#here are the available algorithms for the search
my %avail_algs=(
    "usearch" => 1,
    "blast" => 1
);

my %avail_methods=(
    "lazy" => 1,
    "strict" => 1,
    "strictest" => 1,
    "centroids" =>1

);

my %avail_ann=(
    "prokka" => 1,
    "prodigal" => 1

);



GetOptions ("d=s" => \$dir_genomes,"c=s"   => \$n_cpu, "ann=s"   => \$ann,"m=s"   => \$method,"a=s"   => \$alg,"k=s"   => \$knife, "r=s"   => \$raxml, "i=s"   => \$identity_orthologs,"ip=s"   => \$identity_paralogs, "B=s"   => \$bootstrap_num,"f=s" => $force, "sim=s"=> \$sim) or die("::usage: $0 -d <dir_genomes> -c <n_cpu> -ann <annotation_software> -m <comparison_method> -a <algorithm> -f <force:[0|1]> -k <expression> -r <raxml:[0|1]> -i <perc_identity> -ip <perc_identity_paralogs> -B <bootstraps>\n[ERROR] launch failed! Please check the parameters!\n");

if($dir_genomes eq ""){
    print "::usage: $0 -d <dir_genomes> -c <n_cpu> -ann <annotation_software> -m <comparison_method> -a <algorithm> -f <force:[0|1]> -k <expression> -r <raxml:[0|1]> -i <perc_identity> -ip <perc_identity_paralogs> -B <bootstraps>\n[ERROR] launch failed! Please check the parameters!\n";
    exit();
}


if(!(exists($avail_algs{$alg}))){
    print "::[ERROR] I do not know the search algorithm you specified\n";
    exit();    
}

elsif(!(exists($avail_ann{$ann}))){
    print "::[ERROR] I do not know the annotation software you specified\n";
    exit();    
}

elsif(!(exists($avail_methods{$method}))){
    print "::[ERROR] I do not know the method you specified\n";
    exit();    
}


$date=`date "+%Y-%m-%d %H:%M:%S"`;
# To know the start time of the script
BEGIN { our $start_run = time(); }

print "::3..2..1..and...lift off -- $date";

#prokka annotation -- first stage
    print "::First stage\n";
    

    if($ann eq "prodigal"){

   	 $cmd="satv_prodigal-driver.pl -d ${dir_genomes} -c $n_cpu";
   	 system($cmd);
    }
    elsif($ann eq "prokka"){

        if($knife ne ""){
    	 $cmd="satv_prokka-driver.pl -d ${dir_genomes} -c $n_cpu -k $knife";
     	 system($cmd);
        
        }
        else{
    	 $cmd="satv_prokka-driver.pl -d ${dir_genomes} -c $n_cpu";
     	 system($cmd);
  

        }

    }

 
    
    if((!(-e "genomes_to_analyze.txt"))or($force eq "1")){

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
    if($method eq "lazy"){

        $cmd="satv_search-pangenome-lazy.pl -g genomes_to_analyze.txt -c $n_cpu -i $identity_orthologs -f $force -a $alg -sim $sim";
        system($cmd);
    }
    elsif($method eq "strict"){

        $cmd="satv_search-pangenome-strict.pl -g genomes_to_analyze.txt -c $n_cpu -i $identity_orthologs -ip $identity_paralogs -f $force -a $alg -sim $sim";
        system($cmd);
    }
    elsif($method eq "strictest"){

        $cmd="satv_search-pangenome-strictest.pl -g genomes_to_analyze.txt -c $n_cpu -i $identity_orthologs -ip $identity_paralogs -f $force -a $alg -sim $sim";
        system($cmd);
    }

    elsif($method eq "centroids"){

        $cmd="satv_search-pangenome-centroids.pl -g genomes_to_analyze.txt -c $n_cpu -i $identity_orthologs -f $force";
        system($cmd);
    }





	# if the user want a phylogeny generated from the binary matrix
    if($raxml eq "1"){
	if (-e 'table_linked5_lazy.tsv')
	{
        system("satv_generate-binary-matrix.pl -tab table_linked5_lazy.tsv -out binary_matrix.tsv");
	}

	if (-e 'table_linked5_strict.tsv')
	{
        system("satv_generate-binary-matrix.pl -tab table_linked5_strict.tsv -out binary_matrix.tsv");
	}

	if (-e 'table_linked5_strictest.tsv')
	{
        system("satv_generate-binary-matrix.pl -tab table_linked5_strictest.tsv -out binary_matrix.tsv");
	}

        $cmd="satv_raxml-driver.sh $n_cpu $bootstrap_num";
        system($cmd);

	mkdir 'Phylogeny';
	system("mv RAxML* Phylogeny");
	system("mv binary_matrix.tsv matrix");
	system("mv binary_matrix.fasta matrix");
    }





mkdir 'prot_files';


system("ln *.faa prot_files");
system("ln */*.gff prot_files");
system("rm *.udb");






$date=`date "+%Y-%m-%d %H:%M:%S"`;
# To know how many time was done since the script started
$end_run = time();
$run_time = $end_run - our $start_run;
print "::Mission accomplished by SaturnV in $run_time seconds -- $date\n";






