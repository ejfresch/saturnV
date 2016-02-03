#!/usr/bin/perl -I /home/avincent/Desktop/saturnV/bin

use strict;


use Getopt::Long;
use Graph;


my $genome_list="";
my $n_cpu=1;
my $force=0;
my $identity_orthologs=90;
my $identity_paralogs=90;
my $alg="usearch";

#here are the available algorithms for the search
my %avail_algs=(
    "usearch" => 1,
    "blast" => 1
);


GetOptions ("g=s" => \$genome_list,"c=s"   => \$n_cpu,"i=s"   => \$identity_orthologs,"ip=s"   => \$identity_paralogs,"f=s"   => \$force,"a=s"=> \$alg) or die("::usage: $0 -g <genomes_list> -c <n_cpu> -i <perc_identity_orthlogs> -ip <perc_identity_paralogs> -f <force:[0|1]> -a <algorithm>\n[ERROR] launch failed! Please check the parameters!\n");

if($genome_list eq ""){
    print "::usage: $0 -g <genomes_list> -c <n_cpu> -i <perc_identity_orthlogs> -ip <perc_identity_paralogs> -f <force:[0|1]> -a <algorithm>\n";
    exit();
}

if(!(exists($avail_algs{$alg}))){
    print "::usage: $0 -g <genomes_list> -c <n_cpu> -i <perc_identity_orthlogs> -ip <perc_identity_paralogs> -f <force:[0|1]> -a <algorithm>\n";
    exit();    
}



my $identity_orthologs_usearch=$identity_orthologs/100;




my $date=`date "+%Y-%m-%d %H:%M:%S"`;
print "::Starting analysis at $date";

# concatenation of all faa files 
`cat */*.faa > CONCAT.faa`;

my @genomes=`cat $genome_list`;
chomp(@genomes);


#I generate the db
    if(((!(-e "CONCAT.udb")) or ($force eq "1")) and ($alg eq "usearch")){
    print "::indexing sequences\n";
    `usearch8 -makeudb_usearch CONCAT.faa -output CONCAT.udb >> log_file 2>&1`;
	#print "here\n";
    }
    elsif(((!(-e "CONCAT.pin")) or ($force eq "1")) and ($alg eq "blast")){
    print "::indexing sequences\n";
    `makeblastdb -in CONCAT.faa -dbtype prot > /dev/null`;
    }





#I perform the usearch searches


if($alg eq "usearch"){
        print "::performing $alg searches\n";
`usearch8 -usearch_local CONCAT.faa -threads $n_cpu -db CONCAT.udb -id $identity_orthologs_usearch -userfields query+target+qcov+tcov+id -userout result.aln >> log_file 2>&1`;
        }
        elsif($alg eq "blast"){
              print "::performing $alg searches\n";
             `blastp -query CONCAT.faa -db CONCAT -num_threads 1 -out result.aln -seg no -outfmt 6`;

        }






my %results=();
open(IN,"<result.aln");
while(my $line=<IN>){
     chomp($line);
        #print $line."\n";
        my @data=split(/\t/,$line);

        my $query=$data[0];

        
        my $hit=$data[1];

	my $qcov=$data[2];

	my $tcov=$data[3];

	my $id=$data[4];

	my @explode_genome_h=split(/_/,$hit);
	my $genome_h=$explode_genome_h[0];

	my @explode_genome_q=split(/_/,$query);
	my $genome_q=$explode_genome_q[0];
	
        if($query eq $hit){next;}
    
	if(($qcov < 85) or ($tcov < 85)){next;}

      

           


            if($genome_q eq $genome_h){
                if($id < $identity_paralogs){next;}
                
            }   

            
            if(exists($results{$query}{$genome_h}{"hit"})){

                
                my $prev_hit=$results{$query}{$genome_h}{"hit"};
                my $prev_hit_id=$results{$query}{$genome_h}{"id"};
                
                if($id > $prev_hit_id){


                     $results{$query}{$genome_h}{"hit"}=$hit;
                     $results{$query}{$genome_h}{"id"}=$id;


                }            
                


            }else{            

                $results{$query}{$genome_h}{"hit"}=$hit;
                $results{$query}{$genome_h}{"id"}=$id;

            }
    


    }
   
    




close(IN);





#I write all data
open(OUT,">db_complete_search_filtered.txt");
foreach my $q (keys(%results)){

    my $ref=$results{$q};

    my %hash=%$ref;

    foreach my $k (keys(%hash)){
        print OUT "$q\t".$results{$q}{$k}{"hit"}."\n";

    }



}

close(OUT);


#I create the graph
my $g = Graph->new;


open(IN,"<db_complete_search_filtered.txt");

while(my $line=<IN>){
    chomp($line);
    if($line eq ""){next;}
    my @el_fragmented=split(/\t/,$line);
    #print join("^",@el_fragmented)."\n";        
    $g->add_path(@el_fragmented);

}


close(IN);


print "::perl is determining the strongly connected components of the graph\n";
my @cc=$g->strongly_connected_components();





open(OUT,">table_linked3.tsv");
print OUT "#".join("\t",@genomes)."\n";

my $count=0;
foreach my $c (@cc){
    $count++;
    print "::writing down the strongly connected components of the graph (${count})\n";


    my %hash_results=();
    my @current_cc=@$c;
    foreach my $element (@current_cc){

        if($element eq "-"){next;}
        #print $element."\n";
        my @explode=split(/_/,$element);
        #print $explode[0]."\n";
        $hash_results{$explode[0]}{$element}=1;        


     
    }

    my @final_arr=();    
    
    foreach my $genome (@genomes){

    $genome=~s/.faa//g;    
    #print $genome."\n";
    if(exists($hash_results{$genome})){
       
        my $ref=$hash_results{$genome};

        my %arr=%$ref;
        push(@final_arr,join(",",sort(keys(%arr))));
    }
    else{push(@final_arr,"-");}

    }

    print OUT join("\t",@final_arr)."\n";

}

close(OUT);









print "::Analysis completed at $date";
print "Bye!\n";


