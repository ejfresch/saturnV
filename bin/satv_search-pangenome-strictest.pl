#!/usr/bin/perl

use strict;

use Parallel::ForkManager;
use Getopt::Long;
use Graph;

use Memory::Usage;
my $mu = Memory::Usage->new();

my $genome_list="";
my $n_cpu=1;
my $force=0;
my $identity_orthologs=90;
my $identity_paralogs=90;
my $alg="usearch";
my $file_out="table_linked5_strictest.tsv";
my $sim=0;

#here are the available algorithms for the search
my %avail_algs=(
    "usearch" => 1,
    "blast" => 1
);


GetOptions ("g=s" => \$genome_list,"out=s"   => \$file_out,"c=s"   => \$n_cpu,"i=s"   => \$identity_orthologs,"ip=s"   => \$identity_paralogs,"f=s"   => \$force,"a=s"=> \$alg ,"sim=s"=> \$sim) or die("::usage: $0 -g <genomes_list> -c <n_cpu> -i <perc_identity_orthlogs> -ip <perc_identity_paralogs> -f <force:[0|1]> -a <algorithm>\n[ERROR] launch failed! Please check the parameters!\n");

if($genome_list eq ""){
    print "::usage: $0 -g <genomes_list> -out <file_out> -c <n_cpu> -i <perc_identity_orthlogs> -ip <perc_identity_paralogs> -f <force:[0|1]> -a <algorithm -sim <similarity:[0|1]>>\n";
    exit();
}

if(!(exists($avail_algs{$alg}))){
    print "::usage: $0 -g <genomes_list> -out <file_out> -c <n_cpu> -i <perc_identity_orthlogs> -ip <perc_identity_paralogs> -f <force:[0|1]> -a <algorithm> -sim <similarity:[0|1]\n";
    exit();    
}






my $identity_orthologs_usearch=$identity_orthologs/100;


my $manager = new Parallel::ForkManager($n_cpu);

my $date=`date "+%Y-%m-%d %H:%M:%S"`;
print "::Starting analysis at $date";


my @genomes=`cat $genome_list`;
chomp(@genomes);


#I generate the db

my %db_elements=();
for my $genome (@genomes){


    
    if(((!(-e "${genome}.udb")) or ($force eq "1")) and ($alg eq "usearch")){
    print "::indexing sequences -- ${genome}\n";
    `usearch8 -makeudb_usearch $genome -output ${genome}.udb >> log_file 2>&1`;
	#print "here\n";
    }
    elsif(((!(-e "${genome}.pin")) or ($force eq "1")) and ($alg eq "blast")){
    print "::indexing sequences -- ${genome}\n";
    `makeblastdb -in $genome -dbtype prot > /dev/null`;
    }



    print "determining the length of the sequences contained in genome $genome\n";

    my $id;
    open(IN,"<$genome") or die("::I cannot open the file ${genome}\n");
    while(my $line=<IN>){
        chomp($line);
        
        
        if($line eq ""){next;}
        elsif($line=~/^>/){
            $id=$line;
            $id=~s/>//g;
            $db_elements{$id}{"len"}=0;
            $db_elements{$id}{"gen"}=$genome;


        }
        else{
            $db_elements{$id}{"len"}+=length($line);
            
        }

    }    
    

    close(IN); 



}




#I perform the usearch searches





for(my $i=0;$i<=$#genomes;$i++){

    for(my $j=0;$j<=$#genomes;$j++){
        
        if((!(-e "$genomes[$i]_vs_$genomes[$j]_complete-search.txt"))or ($force eq "1")){

        $manager->start and next;



        if($alg eq "usearch"){
        print "::performing $alg searches -- $genomes[$i] vs $genomes[$j]\n";
`usearch8 -usearch_local $genomes[$i] -threads 1 -db $genomes[$j].udb -id $identity_orthologs_usearch -blast6out $genomes[$i]_vs_$genomes[$j]_complete-search.txt >> log_file 2>&1`;
        }
        elsif($alg eq "blast"){
              print "::performing $alg searches -- $genomes[$i] vs $genomes[$j]\n";
             `blastp -query $genomes[$i] -db $genomes[$j] -num_threads 1 -out $genomes[$i]_vs_$genomes[$j]_complete-search.txt -seg no -outfmt '6 qseqid sseqid pident length ppos'`;

        }




        $manager->finish;
        }

    }



}


$manager->wait_all_children;

`cat *_complete-search.txt > db_complete_search.txt`;



open(O2,">dump_data-searches.txt");

print O2 "#query\thit\tidentity\tlength_ali\tthreshold_query\tthreshold_hit\tdecision\n";

my %results=();


open(IN,"<db_complete_search.txt");
while(my $line=<IN>){
     chomp($line);
        #print $line."\n";
        my @data=split(/\t/,$line);

        my $query=$data[0];

        
        my $hit=$data[1];

        if($query eq $hit){next;}
    
	my $length_query;
	my $genome_q;
	my $length_hit;
	my $genome_h;
	my $perc_identity;
	my $length_alignment;


        if(exists($db_elements{$query})){
            $length_query=$db_elements{$query}{"len"};
            $genome_q=$db_elements{$query}{"gen"};
        }
        else{
            print "::[ERROR] sequence $query not present in db of the hash db_elements\n";
        } 


    
        if(exists($db_elements{$hit})){
            $length_hit=$db_elements{$hit}{"len"};
            $genome_h=$db_elements{$hit}{"gen"};
        }
        else{
            print "::[ERROR] sequence $hit not present in db of the hash db_elements\n";
        } 


	if($sim eq "0"){
        $perc_identity=$data[2]; # identity
	}
	elsif($sim eq "1"){
	$perc_identity=$data[4]; # similarity
	}

        #useful for blast
        if($perc_identity < $identity_orthologs){next;}
        
        $length_alignment=$data[3]; # alignment


        print O2 $query."\t".$hit."\t".$perc_identity."\t".$length_alignment."\t".($length_query*0.85)."\t".($length_hit*0.85)."\t";


        
       # print $query."-".$hit."-"."$identity"."-".$length_alignment."-".($length_query*0.85)."-".($length_hit*0.85)."\n";

        if(($length_alignment >= ($length_query*0.85)) and ($length_alignment >= ($length_hit*0.85)) ){
           
            print O2 "Yes\n"; 

            if($genome_q eq $genome_h){
                if($perc_identity < $identity_paralogs){next;}
                
            }   

            
            if(exists($results{$query}{$genome_h}{"hit"})){

                
                my $prev_hit=$results{$query}{$genome_h}{"hit"};
                my $prev_hit_id=$results{$query}{$genome_h}{"id"};
                
                if($perc_identity > $prev_hit_id){


                     $results{$query}{$genome_h}{"hit"}=$hit;
                     $results{$query}{$genome_h}{"id"}=$perc_identity;


                }            
                


            }else{            
            
                $results{$query}{$genome_h}{"hit"}=$hit;
                $results{$query}{$genome_h}{"id"}=$perc_identity;
                         
            }
        
        }else{print O2 "No\n";}
    


    }
   
    




close(IN);

close(O2);



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


foreach my $element (keys(%db_elements)){
    if(!($g->has_vertex($element))){
        $g->add_vertex($element);        
    }    
}


my $date=`date "+%Y-%m-%d %H:%M:%S"`;
print "::perl is determining the strongly connected components of the graph -- $date";
my @cc=$g->strongly_connected_components();

$mu->record('mem_conn_comp');

my $ref_mem=$mu->state();
my @data_mem=@$ref_mem;
my $ref_data_graph=$data_mem[0];
my @data_graph=@$ref_data_graph;

print "::I used ".sprintf("%.2f",($data_graph[2]/1024))." Mb (Virtual Memory)\n";





my $date=`date "+%Y-%m-%d %H:%M:%S"`;
print "::writing down the strongly connected components of the graph (".($#cc+1).") -- $date";


open(OUT,">${file_out}");
print OUT "#".join("\t",@genomes)."\n";

foreach my $c (@cc){
    
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








my $date=`date "+%Y-%m-%d %H:%M:%S"`;
print "::Analysis completed at $date";
print "Bye!\n";


