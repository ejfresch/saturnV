#!/usr/bin/perl

use Getopt::Long;


$genome_list="";
$identity_orthologs="90";
$n_cpu=1;


GetOptions ("g=s" => \$genome_list,"c=s"   => \$n_cpu,"i=s"   => \$identity_orthologs,"ip=s"   => \$identity_paralogs,"f=s"   => \$force,"a=s"=> \$alg) or die("::usage: $0 -g <genomes_list> -c <n_cpu> -i <perc_identity_orthlogs> -ip <perc_identity_paralogs> -f <force:[0|1]> -a <algorithm>\n[ERROR] launch failed! Please check the parameters!\n");

if($genome_list eq ""){
    print "::usage: $0 -g <genomes_list> -c <n_cpu> -i <perc_identity_orthlogs> -ip <perc_identity_paralogs> -f <force:[0|1]> -a <algorithm>\n";
    exit();
}



$identity_orthologs_usearch=$identity_orthologs/100;

#I merge all .faa files

print "::merging .faa files\n";
@genomes=`cat $genome_list`;
chomp(@genomes);
$cmd="cat ".join(" ",@genomes)."> CONCAT.faa";
print $cmd."\n";
system($cmd);


print "::generating clusters\n";

$cmd="usearch8 -threads $n_cpu -cluster_fast CONCAT.faa -id $identity_orthologs_usearch -uc table_results_centroids.txt";
system($cmd);


open(OUT,">table_links_centroids.txt");

print OUT "#".join("\t",@genomes)."\n";
open(IN,"<table_results_centroids.txt");

%results=();


while($line=<IN>){
    chomp($line);
    $count++;
    @explode=split(/\t/,$line);

    $flag=$explode[0];
    $group=$explode[1];

    if($flag eq "C"){next;}
        
    $query=$explode[8];    
    $hit=$explode[9];
    
    #print $query."\n";
    #print $hit."\n";

    @explode_query=split(/_/,$query);
    @explode_hit=split(/_/,$hit);


    $genome_q=$explode_query[0];
    $genome_h=$explode_hit[0];
 
    $results{$group}{$genome_q}{$query}=1;   

    if($hit ne "*"){
    $results{$group}{$genome_h}{$hit}=1;   
    }

  }


close(IN);

foreach $group(keys(%results)){


my @final_arr=();    
    
    foreach my $genome (@genomes){

    $genome=~s/.faa//g;    
    #print $genome."\n";
    if(exists($results{$group}{$genome})){
       
        my $ref=$results{$group}{$genome};

        my %arr=%$ref;
        push(@final_arr,join(",",sort(keys(%arr))));
    }
    else{push(@final_arr,"-");}

    }

    print OUT join("\t",@final_arr)."\n";



}

close(OUT);




