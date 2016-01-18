#!/usr/bin/perl -I /project/rclevesq/users/lfreschi/tasks/pangenome/saturnv/bin

use LibFASTA;
use Parallel::ForkManager;
use Getopt::Long;


$genome_list="";
$n_cpu=1;
$identity_usearch=90;



GetOptions ("g=s" => \$genome_list,"c=s"   => \$n_cpu,"i=s"   => \$identity_usearch) or die("::usage: $0 -d <genomes_list> -c <n_cpu> -i <identity_usearch>\n");


if($genomes_list eq ""){

    print "::usage: $0 -g <genomes_list> -c <n_cpu> -i <identity_usearch>\n";
    exit();
}



$manager = new Parallel::ForkManager($n_cpu);

$date=`date "+%Y-%m-%d %H:%M:%S"`;
print "::Starting analysis at $date";


@genomes=`cat $genome_list`;
chomp(@genomes);


#I generate the db
for $genome (@genomes){


    
    if(!(-e "${genome}.udb")){
    print "::indexing sequences -- ${genome}\n";
    `usearch8 -makeudb_usearch $genome -output ${genome}.udb >> log_file 2>&1`;
    }
}




$first_genome=$genomes[0];

#I perform the usearch searches


$new_blasts=0;

for $genome (@genomes[0..$#genomes]){
    
  if(!(-e "${genome}_blasted.txt")){  
    
   $new_blasts++;  
  
  $manager->start and next;
    print "::performing usearch searches -- ${first_genome} vs ${genome}\n";
    `usearch8 -usearch_local $first_genome -threads 1 -db ${genome}.udb -id $identity_usearch -blast6out ${genome}_blasted.txt >> log_file 2>&1`;

  $manager->finish;

    

  }  

}


$manager->wait_all_children;




%db_all_seq=();

#I build the hash with the sequences
print "::building the hash with the sequences\n";

@ids_first_genome=();

$ref=&LibFASTA::read_FASTA("$genomes[0]"," ","1");
%current_hash=%$ref;

    foreach $elem (keys(%current_hash)){
        
        if(exists($db_all_seq{$elem})){print "::[ERROR] Hey! There are two or more sequences with the same ID -- ID: $elem !\n";}

        $db_all_seq{$elem}{"seq"}=$current_hash{$elem};
        $db_all_seq{$elem}{"gen"}=$genomes[0];
        push(@ids_first_genome,$elem);

    }



for $genome (@genomes[1..$#genomes]){

$ref=&LibFASTA::read_FASTA("$genome"," ","1");
%current_hash=%$ref;

    foreach $elem (keys(%current_hash)){

        $db_all_seq{$elem}{"seq"}=$current_hash{$elem};
        $db_all_seq{$elem}{"gen"}=$genome;


    }



}

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

        $genome_q=$db_all_seq{$query}{"gen"};
        #print $genome_q."\n";
        $results{$query}{$genome_q}{$query}=1;


        if(exists($db_all_seq{$query})){
            $length_query=length($db_all_seq{$query}{"seq"});
        }
        else{
            print "::[ERROR] sequence $query not present in db of the hash %db_all_seq\n";
        }
        
        $hit=$data[1];
        
        
        if(exists($db_all_seq{$hit})){
            $length_hit=length($db_all_seq{$hit}{"seq"});
            $genome_id=$db_all_seq{$hit}{"gen"};
            
        }
        else{
            print "::[ERROR] sequence $hit not present in db of the hash %db_all_seq\n";
        } 



        $length_alignment=$data[3];
        
       # print $query."-".$hit."-"."$identity"."-".$length_alignment."-".($length_query*0.85)."-".($length_hit*0.85)."\n";

        if(($length_alignment >= ($length_query*0.85)) and ($length_alignment >= ($length_hit*0.85)) ){
            
            $results{$query}{$genome_id}{$hit}=1; 
            $exclusion_zone{$hit}=1;   
        
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




            if(exists($results{$entry}{$gen_id})){
            
            $ref=$results{$entry}{$gen_id};
            %hash_ref=%$ref;


            @hits=sort(keys(%hash_ref));
            push(@arr_to_print,join(",",@hits));

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
foreach $seq (keys(%db_all_seq)){
    
    if(exists($exclusion_zone{$seq})){
        next;
    }else{
        print OUT ">$seq\n";
        print OUT $db_all_seq{$seq}{"seq"}."\n";
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

if(!(-e "new_sequences_iter1.faa.udb")){
    `usearch8 -makeudb_usearch new_sequences_iter1.faa -output new_sequences_iter1.faa.udb >> log_file 2>&1`;


}
# I run the comparisons

$new_blasts_iter2=0;

for $genome (@genomes){

  if(!(-e "${genome}_blasted_iter2.txt")){  

  $new_blasts_iter2++;
  $manager->start and next;
    print "::performing usearch searches -- new_sequences_iter1 vs ${genome}\n";

    `usearch8 -usearch_local new_sequences_iter1.faa -threads 1 -db ${genome}.udb -id $identity_usearch -blast6out ${genome}_blasted_iter2.txt >> log_file 2>&1`;

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

        $genome_q=$db_all_seq{$query}{"gen"};
        #print $genome_q."\n";
        $results{$query}{$genome_q}{$query}=1;


        if(exists($db_all_seq{$query})){
            $length_query=length($db_all_seq{$query}{"seq"});
        }
        else{
            print "::[ERROR] sequence $query not present in db of the hash %db_all_seq\n";
        }
    
        $hit=$data[1];
    
      
    
        if(exists($db_all_seq{$hit})){
            $length_hit=length($db_all_seq{$hit}{"seq"});
            $genome_id=$db_all_seq{$hit}{"gen"};
        }
        else{
            print "::[ERROR] sequence $hit not present in db of the hash %db_all_seq\n";
        } 



        $length_alignment=$data[3];


        if(($length_alignment >= ($length_query*0.85)) and ($length_alignment >= ($length_hit*0.85)) ){
        
            $results{$query}{$genome_id}{$hit}=1; 
            $exclusion_zone{$hit}=1;   
    
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


            if(exists($results{$entry}{$gen_id})){
            $ref_r=$results{$entry}{$gen_id};
            %hash_ref=%$ref_r;


            @hits=sort(keys(%hash_ref));
            push(@arr_to_print,join(",",@hits));

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



print "::Analysis completed at $date";
print "Bye!\n";

