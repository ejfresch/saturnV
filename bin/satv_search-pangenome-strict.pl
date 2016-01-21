#!/usr/bin/perl -I /project/rclevesq/users/lfreschi/tasks/pangenome/saturnv/bin

use LibFASTA;
use Parallel::ForkManager;
use Getopt::Long;


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


GetOptions ("g=s" => \$genome_list,"c=s"   => \$n_cpu,"i=s"=> \$identity_orthologs,"ip=s"   => \$identity_paralogs) or die("::usage: $0 -g <genomes_list> -c <n_cpu> -i <perc_identity_orthologs> -ip <perc_identity_paralogs> -f <force:[0|1]> -a <algorithm>\n[ERROR] launch failed! Please check the parameters!\n");

if($genome_list eq ""){

    print "::usage: $0 -g <genomes_list> -c <n_cpu> -i <perc_identity_orthologs> -ip <perc_identity_paralogs> -f <force:[0|1]> -a <algorithm>\n";
    exit();
}


$identity=$identity_orthologs/100;


$manager = new Parallel::ForkManager($n_cpu);

$date=`date "+%Y-%m-%d %H:%M:%S"`;
print "::Starting analysis at $date";


@genomes=`cat $genome_list`;
chomp(@genomes);


#I generate the db
%db_elements=();
@ids_first_genome;

$first_genome=$genomes[0];

for $genome (@genomes){


    
    if(((!(-e "${genome}.udb"))or ($force eq "1")) and ($alg eq "usearch")){
    print "::indexing sequences -- ${genome}\n";
    `usearch8 -makeudb_usearch $genome -output ${genome}.udb >> log_file 2>&1`;

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
            $db_elements{$id}{"seq"}="";
            
            if($genome eq $first_genome){
                push(@ids_first_genome,$id);
            }

        }
        else{
            $db_elements{$id}{"len"}+=length($line);
            $db_elements{$id}{"seq"}.=$line;
        }

    }    
    

    close(IN); 



}








#I perform the usearch searches


$new_blasts=0;



for $genome (@genomes[0..$#genomes]){
    
  if((!(-e "${genome}_blasted.txt"))or ($force eq "1")){  
    
   $new_blasts++;  
  
  $manager->start and next;

    print "::performing usearch searches -- ${first_genome} vs ${genome}\n";

    if($alg eq "usearch"){
    `usearch8 -usearch_local $first_genome -threads 1 -db ${genome}.udb -id $identity -blast6out ${genome}_blasted.txt >> log_file 2>&1`;
    }
    elsif{
    `blastp -query $first_genome -db $genome -num_threads 1 -out ${genome}_blasted.txt -seg no -outfmt 6`;    
    }


  $manager->finish;

    

  }  









}


$manager->wait_all_children;


print "::analyzing the usearch output\n";

#now I can read the usearch output and draw some conclusions


if($new_blasts > 0){


    if(-e "all_blasted.txt"){`rm -rf all_blasted.txt`;}
    `cat *_blasted.txt > all_blasted.txt`;



    #I initialize the hash of the results
    %results=();

    #sequences for which I already found a match and therefore can be excluded from further analysis
    %exclusion_zone=();

    open(IN,"all_blasted.txt")||die "I cannot open all_blasted.txt";
    while($line=<IN>){
        chomp($line);
        #print $line."\n";
        @data=split(/\t/,$line);

        $query=$data[0];
        $exclusion_zone{$query}=1;

        $genome_q=$db_elements{$query}{"gen"};
        #print $genome_q."\n";
        $results{$query}{$genome_q}{$query}=1;


        if(exists($db_elements{$query})){
            $length_query=$db_elements{$query}{"len"};
        }
        else{
            print "::[ERROR] sequence $query not present in db of the hash %db_all_seq\n";
        }
        
        $hit=$data[1];

        if($query eq $hit){next;}
        
        
        if(exists($db_elements{$hit})){
            $length_hit=$db_elements{$hit}{"len"};
            $genome_h=$db_elements{$hit}{"gen"};
            
        }
        else{
            print "::[ERROR] sequence $hit not present in db of the hash %db_all_seq\n";
        } 



        $perc_identity=$data[2];

        if($perc_identity<$identity_orthologs){next;}

        $length_alignment=$data[3];
        
       # print $query."-".$hit."-"."$identity"."-".$length_alignment."-".($length_query*0.85)."-".($length_hit*0.85)."\n";

        if(($length_alignment >= ($length_query*0.85)) and ($length_alignment >= ($length_hit*0.85)) ){
            

            if($genome_q eq $genome_h){
                if($perc_identity < $identity_paralogs){next;}
                
            }        





            
            if(exists($results{$query}{$genome_h}{"hit"})){
                
                $prev_hit=$results{$query}{$genome_h}{"hit"};
                $prev_hit_id=$results{$query}{$genome_h}{"id"};
                
                if($perc_identity > $prev_hit_id){


                     $results{$query}{$genome_h}{"hit"}=$hit;
                     $results{$query}{$genome_h}{"id"}=$perc_identity;
                     $exclusion_zone{$hit}=1;
                     delete($exclusion_zone{$prev_hit});
                }            
                


            }else{            

                $results{$query}{$genome_h}{"hit"}=$hit;
                $results{$query}{$genome_h}{"id"}=$perc_identity;
                $exclusion_zone{$hit}=1;
            }
        
        }
        


    }
    close(IN);
    


    



    #foreach $test (keys(%results)){

    #    print $test." >> ";
    #    $rr=$results{$test};
    #    %hash_test=%$rr;
    #    @keys_test=keys(%hash_test);
    #    print join(",",@keys_test)."\n";  

    #}



    open(OUT,">situation_iter1.txt");
    print OUT "#".join("\t",@genomes)."\n";
    foreach $entry (@ids_first_genome){


         @arr_to_print=();
         push(@arr_to_print,$entry);

        foreach $gen (@genomes[1..$#genomes]){
            $gen_id=$gen;




            if(exists($results{$entry}{$gen_id}{"hit"})){
            
   
            push(@arr_to_print,$results{$entry}{$gen_id}{"hit"});

            }
            else{
                push(@arr_to_print,"-");

            }

        }


        print OUT join("\t",@arr_to_print)."\n";

    }
    close(OUT);



%new_sequences=();








#I write the remaining sequences
open(OUT,">new_sequences_iter1.faa");
foreach $seq (keys(%db_elements)){
    
    if(exists($exclusion_zone{$seq})){
        next;
    }else{
        print OUT ">$seq\n";
        print OUT $db_elements{$seq}{"seq"}."\n";
        $new_sequences{$seq}=1;
    }

}

close(OUT);



}
else{
    print "::I have already seen all these genomes. No need to analyze the same data twice! :)\n";
}
print "::First iteration finished!\n";



#I reinitialize the array of the results
undef %results;
%results=();



#I build the db

if(((!(-e "new_sequences_iter1.faa.udb")) or ($force eq "1")) and ($alg eq "usearch")){
    `usearch8 -makeudb_usearch new_sequences_iter1.faa -output new_sequences_iter1.faa.udb >> log_file 2>&1`;

}
elsif(((!(-e "new_sequences_iter1.faa.pin")) or ($force eq "1")) and ($alg eq "blast")){

     `makeblastdb -in new_sequences_iter1.faa -dbtype prot > /dev/null`;

}



# I run the comparisons

$new_blasts_iter2=0;

for $genome (@genomes){

  if((!(-e "${genome}_blasted_iter2.txt")) or ($force eq "1")){  

  $new_blasts_iter2++;
  $manager->start and next;
    print "::performing usearch searches -- new_sequences_iter1 vs ${genome}\n";



    

    if($alg eq "usearch"){
`usearch8 -usearch_local new_sequences_iter1.faa -threads 1 -db ${genome}.udb -id $convert_id_usearch -blast6out ${genome}_blasted_iter2.txt >> log_file 2>&1`;
    }
    elsif{
    `blastp -query new_sequences_iter1.faa -db $genome -num_threads 1 -out ${genome}_blasted_iter2.txt -seg no -outfmt 6`;    
    }





  $manager->finish;
    


    }

}
$manager->wait_all_children;


print "::analyzing the usearch output\n";

#now I can read the blast output and draw some conclusions

print "new_blasts_iter2:=".$new_blasts_iter2."\n";

if($new_blasts_iter2 > 0){

    if(-e "all_blasted_iter2.txt"){`rm -rf all_blasted_iter2.txt`;}
    `cat *_blasted_iter2.txt > all_blasted_iter2.txt`;


    open(IN,"all_blasted_iter2.txt")||die "I cannot open all_blasted_iter2.txt";
    while($line=<IN>){
        chomp($line);
        #print $line."\n";
        @data=split(/\t/,$line);

        $query=$data[0];

        $genome_q=$db_elements{$query}{"gen"};
        #print $genome_q."\n";
        $results{$query}{$genome_q}{$query}=1;


        if(exists($db_elements{$query})){
            $length_query=$db_elements{$query}{"len"};
        }
        else{
            print "::[ERROR] sequence $query not present in db of the hash %db_elements\n";
        }
    
        $hit=$data[1];
    
        if($query eq $hit){next;}
    
        if(exists($db_elements{$hit})){
            $length_hit=$db_elements{$hit}{"len"};
            $genome_h=$db_elements{$hit}{"gen"};
        }
        else{
            print "::[ERROR] sequence $hit not present in db of the hash %db_elements\n";
        } 



        $perc_identity=$data[2];

        if($perc_identity<$identity_orthologs){next;}
        $length_alignment=$data[3];
        
       # print $query."-".$hit."-"."$identity"."-".$length_alignment."-".($length_query*0.85)."-".($length_hit*0.85)."\n";

        if(($length_alignment >= ($length_query*0.85)) and ($length_alignment >= ($length_hit*0.85)) ){
            

            $genome_h=$db_elements{$hit}{"gen"};
            if($genome_q eq $genome_h){
                if($perc_identity < $identity_paralogs){next;}
                
            }   

            
            if(exists($results{$query}{$genome_h}{"hit"})){
                
                $prev_hit=$results{$query}{$genome_h}{"hit"};
                $prev_hit_id=$results{$query}{$genome_h}{"id"};
                
                if($perc_identity > $prev_hit_id){


                     $results{$query}{$genome_h}{"hit"}=$hit;
                     $results{$query}{$genome_h}{"id"}=$perc_identity;


                }            
                


            }else{            

                $results{$query}{$genome_h}{"hit"}=$hit;
                $results{$query}{$genome_h}{"id"}=$perc_identity;

            }
        
        }
    


    }
    close(IN);


    open(OUT,">situation_iter2.txt");
    print OUT "#".join("\t",@genomes)."\n";

    $ref=&LibFASTA::read_FASTA("new_sequences_iter1.faa"," ","1");
    %new_sequences=%$ref;
    


    foreach $entry (sort(keys(%new_sequences))){


        #print OUT $entry."\n";
        @arr_to_print=();


        foreach $gen (@genomes){
            $gen_id=$gen;


            if(exists($results{$entry}{$gen_id}{"hit"})){

                $target=$results{$entry}{$gen_id}{"hit"};
    
                if(exists($exclusion_zone{$target})){                        
                    push(@arr_to_print,"-");
                }
                else{        
                push(@arr_to_print,$target);
                }
            }
            else{
                push(@arr_to_print,"-");

            }

        }


        print OUT join("\t",@arr_to_print)."\n";

    }
    close(OUT);



}
else{
    print "::I have already made all the blasts of the second iteration. Everything is up-to-date! :)\n";
}




print "::Second iteration finished!\n";
$date=`date "+%Y-%m-%d %H:%M:%S"`;



$cmd="satv_merge-data3.pl -in situation_all.txt -out pre_table_linked.tsv";
system($cmd);


$cmd="satv_untie-knots-paralogs.pl -in pre_table_linked.tsv -out table_linked4.tsv";
system($cmd);




print "::Analysis completed at $date";
print "Bye!\n";


