#!/usr/bin/perl -I /home/avincent/Desktop/saturnV/bin


use Parallel::ForkManager;
use Getopt::Long;
use Graph;

$genome_list="";
$n_cpu=1;
$force=0;
$identity_orthologs=90;
$identity_paralogs=90;
$alg="usearch";

#here are the available algorithms for the search
%avail_algs=(
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



$identity_orthologs_usearch=$identity_orthologs/100;


$manager = new Parallel::ForkManager($n_cpu);

$date=`date "+%Y-%m-%d %H:%M:%S"`;
print "::Starting analysis at $date";


@genomes=`cat $genome_list`;
chomp(@genomes);


#I generate the db

%db_elements=();
for $genome (@genomes){


    
    if(((!(-e "${genome}.udb")) or ($force eq "1")) and ($alg eq "usearch")){
    print "::indexing sequences -- ${genome}\n";
    `usearch8 -makeudb_usearch $genome -output ${genome}.udb >> log_file 2>&1`;
	print "here\n";
    }
    elsif(((!(-e "${genome}.pin")) or ($force eq "1")) and ($alg eq "blast")){
    print "::indexing sequences -- ${genome}\n";
    `makeblastdb -in $genome -dbtype prot > /dev/null`;
    }



    print "determining the length of the sequences contained in genome $genome\n";
    open(IN,"<$genome") or die("::I cannot open the file ${genome}\n");
    while($line=<IN>){
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





for($i=0;$i<=$#genomes;$i++){

    for($j=0;$j<=$#genomes;$j++){
        
        if((!(-e "$genomes[$i]_vs_$genomes[$j]_complete-search.txt"))or ($force eq "1")){

        $manager->start and next;
        print "::performing usearch searches -- $genomes[$i] vs $genomes[$j] -- algorithm: $alg\n";


        if($alg eq "usearch"){`usearch8 -usearch_local $genomes[$i] -threads 1 -db $genomes[$j].udb -id $identity_orthologs_usearch -blast6out $genomes[$i]_vs_$genomes[$j]_complete-search.txt >> log_file 2>&1`;
        }
        elsif($alg eq "blast"){
             `blastp -query $first_genome -db $genome -num_threads 1 -out $genome[$i]_vs_$genome[$j]_complete-search.txt -seg no -outfmt 6`;

        }




        $manager->finish;
        }

    }



}


$manager->wait_all_children;

`cat *_complete-search.txt > db_complete_search.txt`;


%results=();
open(IN,"<db_complete_search.txt");
while($line=<IN>){
     chomp($line);
        #print $line."\n";
        @data=split(/\t/,$line);

        $query=$data[0];

        
        $hit=$data[1];

        if($query eq $hit){next;}
    
        if(exists($db_elements{$query})){
            $length_hit=$db_elements{$query}{"len"};
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



        $perc_identity=$data[2];

        #useful for blast
        if($perc_identity < $identity_orthologs){next;}
        
        $length_alignment=$data[3];
        
       # print $query."-".$hit."-"."$identity"."-".$length_alignment."-".($length_query*0.85)."-".($length_hit*0.85)."\n";

        if(($length_alignment >= ($length_query*0.85)) and ($length_alignment >= ($length_hit*0.85)) ){
            

            if($genome_q eq $genome_h){
                if($perc_identity < $identity_paralogs){next;}
                
            }   

            
            if(exists($results{$query}{$genome_id}{"hit"})){
                
                $prev_hit=$results{$query}{$genome_id}{"hit"};
                $prev_hit_id=$results{$query}{$genome_id}{"id"};
                
                if($perc_identity > $prev_hit_id){


                     $results{$query}{$genome_id}{"hit"}=$hit;
                     $results{$query}{$genome_id}{"id"}=$perc_identity;


                }            
                


            }else{            

                $results{$query}{$genome_id}{"hit"}=$hit;
                $results{$query}{$genome_id}{"id"}=$perc_identity;

            }
        
        }
    


    }
   
    




close(IN);


#I write all data
open(OUT,">db_complete_search_filtered.txt");
foreach $q (keys(%results)){

    $ref=$results{$q};

    %hash=%$ref;

    foreach $k (keys(%hash)){
        print OUT "$q\t".$results{$q}{$k}{"hit"}."\n";

    }



}

close(OUT);


#I create the graph
my $g = Graph->new;


open(IN,"<db_complete_search_filtered.txt");

while($line=<IN>){
    chomp($line);
    if($line eq ""){next;}
    @el_fragmented=split(/\t/,$line);
    #print join("^",@el_fragmented)."\n";        
    $g->add_path(@el_fragmented);

}


close(IN);


print "::perl is determining the strongly connected components of the graph\n";
@cc=$g->strongly_connected_components();





open(OUT,">table_linked3.tsv");
print OUT join("\t",@genomes)."\n";

$count=0;
foreach $c (@cc){
    $count++;
    print "::writing down the strongly connected components of the graph (${count})\n";


    %hash_results=();
    @current_cc=@$c;
    foreach $element (@current_cc){

        if($element eq "-"){next;}
        print $element."\n";
        @explode=split(/_/,$element);
        print $explode[0]."\n";
        $hash_results{$explode[0]}{$element}=1;        


     
    }

    @final_arr=();    
    
    foreach $genome (@genomes){

    $genome=~s/.faa//g;    
    print $genome."\n";
    if(exists($hash_results{$genome})){
       
        $ref=$hash_results{$genome};

        %arr=%$ref;
        push(@final_arr,join(",",sort(keys(%arr))));
    }
    else{push(@final_arr,"-");}

    }

    print OUT join("\t",@final_arr)."\n";

}

close(OUT);









print "::Analysis completed at $date";
print "Bye!\n";


