#!/usr/bin/perl -I /home/avincent/Desktop/saturnV/bin

use strict;
use LibFASTA;
use Parallel::ForkManager;
use Getopt::Long;
use Graph;


my $genome_list="";
my $n_cpu=1;
my $force=0;
my $identity_orthologs=90;
my $identity_paralogs=90;
my $alg="usearch";
my $file_out="table_linked5_strict.tsv";
my $sim=0;
my $perc_identity;

#here are the available algorithms for the search
my %avail_algs=(
"usearch" => 1,
"blast" => 1
);


GetOptions ("g=s" => \$genome_list,"out=s"   => \$file_out,"c=s"   => \$n_cpu,"i=s"   => \$identity_orthologs,"ip=s"   => \$identity_paralogs,"f=s"   => \$force,"a=s"=> \$alg, "sim=s"=> \$sim) or die("::usage: $0 -g <genomes_list> -c <n_cpu> -i <perc_identity_orthlogs> -ip <perc_identity_paralogs> -f <force:[0|1]> -a <algorithm> -sim <similarity:[0|1]\n[ERROR] launch failed! Please check the parameters!\n");
if($genome_list eq ""){

    print "::usage: $0 -g <genomes_list> -out <out_file> -c <n_cpu> -i <perc_identity_orthologs> -ip <perc_identity_paralogs> -f <force:[0|1]> -a <algorithm> -sim <similarity:[0|1]\n";
    exit();
}





if(($alg ne "blast")and ($sim eq "1")){
  print "--you cannot use the sim option with the usearch algorithm!\n";
    exit();

}

my $identity=$identity_orthologs/100;


my $manager = new Parallel::ForkManager($n_cpu);

my $date=`date "+%Y-%m-%d %H:%M:%S"`;
print "::Starting analysis at $date";


my @genomes=`cat $genome_list`;
chomp(@genomes);


#I generate the db
my %db_elements=();
my @ids_first_genome;

my $first_genome=$genomes[0];

for my $genome (@genomes){


    
    if((((!(-e "${genome}.udb"))or ($force eq "1"))) and ($alg eq "usearch")){
    print "::indexing sequences -- ${genome}\n";
    `usearch8 -makeudb_usearch $genome -output ${genome}.udb >> log_file 2>&1`;

    }
    elsif((((!(-e "${genome}.pin")) or ($force eq "1"))) and ($alg eq "blast")){
    print "::indexing sequences -- ${genome}\n";
    `makeblastdb -in $genome -dbtype prot > /dev/null`;
    }




    print "--determining the length of the sequences contained in genome $genome\n";

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
            $db_elements{$id}{"seq"}="";
            $db_elements{$id}{"seen"}=0;
            
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


my $new_blasts=0;



for my $genome (@genomes[0..$#genomes]){
    
  if((!(-e "${genome}_blasted.txt"))or ($force eq "1")){  
    
   $new_blasts++;  
  
  $manager->start and next;

    print "::performing $alg searches -- ${first_genome} vs ${genome}\n";

    if($alg eq "usearch"){
    `usearch8 -usearch_local $first_genome -threads 1 -db ${genome}.udb -id $identity -blast6out ${genome}_blasted.txt >> log_file 2>&1`;
    }
    elsif($alg eq "blast"){
    `blastp -query $first_genome -db $genome -num_threads 1 -out ${genome}_blasted.txt -seg no -outfmt '6 qseqid sseqid pident length ppos'`;    
    }


  $manager->finish;

    

  }  









}


$manager->wait_all_children;


print "::analyzing the usearch output\n";

#now I can read the usearch output and draw some conclusions


#I initialize the hash of the results
my %results=();

#sequences for which I already found a match and therefore can be excluded from further analysis
my %exclusion_zone=();

print "--new_blasts:=".$new_blasts."\n";

if($new_blasts > 0){


    if(-e "all_blasted.txt"){`rm -rf all_blasted.txt`;}
    `cat *_blasted.txt > all_blasted.txt`;


    open(IN,"<all_blasted.txt")||die "I cannot open all_blasted.txt";
    while(my $line=<IN>){
        chomp($line);
        #print $line."\n";
        my @data=split(/\t/,$line);

        my $query=$data[0];
        $exclusion_zone{$query}=1;

        my $genome_q=$db_elements{$query}{"gen"};
        #print $genome_q."\n";
        $results{$query}{$genome_q}{$query}=1;


	my $length_query;

        if(exists($db_elements{$query})){
            $length_query=$db_elements{$query}{"len"};
        }
        else{
            print "::[ERROR] sequence $query not present in db of the hash %db_all_seq\n";
        }
        
        my $hit=$data[1];

        if($query eq $hit){next;}
        

	my $length_hit;
	my $genome_h;
        
        if(exists($db_elements{$hit})){
            $length_hit=$db_elements{$hit}{"len"};
            $genome_h=$db_elements{$hit}{"gen"};
            
        }
        else{
            print "::[ERROR] sequence $hit not present in db of the hash %db_all_seq\n";
        } 



               	if($sim eq "0"){
        $perc_identity=$data[2]; # identity
	}
	elsif($sim eq "1"){
	$perc_identity=$data[4]; # similarity
	}

        if($perc_identity<$identity_orthologs){next;}

        my $length_alignment=$data[3];
        
       # print $query."-".$hit."-"."$identity"."-".$length_alignment."-".($length_query*0.85)."-".($length_hit*0.85)."\n";

        if(($length_alignment >= ($length_query*0.85)) and ($length_alignment >= ($length_hit*0.85)) ){
            

            if($genome_q eq $genome_h){
                if($perc_identity < $identity_paralogs){next;}
                
            }        




	    my $prev_hit;
	    my $prev_hit_id;	
            
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
    


    



    #foreach $test (keys(%results)){$cmd="usearch8 -usearch_local new_sequences_iter1.faa -threads 1 -db ${genome}.udb -id $convert_id_usearch -blast6out ${genome}_blasted_iter2.txt >> log_file 2>&1";

    #    print $test." >> ";
    #    $rr=$results{$test};
    #    %hash_test=%$rr;
    #    @keys_test=keys(%hash_test);
    #    print join(",",@keys_test)."\n";  

    #}



    open(OUT,">situation_iter1.txt");
    print OUT "#".join("\t",@genomes)."\n";
    foreach my $entry (@ids_first_genome){


        $db_elements{$entry}{"seen"}=1;

         my @arr_to_print=();
         push(@arr_to_print,$entry);

        foreach my $gen (@genomes[1..$#genomes]){
            my $gen_id=$gen;




            if(exists($results{$entry}{$gen_id}{"hit"})){
            
            $db_elements{$results{$entry}{$gen_id}{"hit"}}{"seen"}=1;
            push(@arr_to_print,$results{$entry}{$gen_id}{"hit"});

            }
            else{
                push(@arr_to_print,"-");

            }

        }


        print OUT join("\t",@arr_to_print)."\n";

    }
    close(OUT);



my %new_sequences=();








#I write the remaining sequences
open(OUT,">new_sequences_iter1.faa");
foreach my $seq (keys(%db_elements)){
    
    if(exists($exclusion_zone{$seq})){
        next;
    }else{
	#print $seq."\n";
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
my %results=();



#I build the db

if(((!(-e "new_sequences_iter1.faa.udb")) or ($force eq "1")) and ($alg eq "usearch")){
    `usearch8 -makeudb_usearch new_sequences_iter1.faa -output new_sequences_iter1.faa.udb >> log_file 2>&1`;

}
elsif(((!(-e "new_sequences_iter1.faa.pin")) or ($force eq "1")) and ($alg eq "blast")){

     `makeblastdb -in new_sequences_iter1.faa -dbtype prot > /dev/null`;

}



# I run the comparisons

my $new_blasts_iter2=0;

for my $genome (@genomes){

  if((!(-e "${genome}_blasted_iter2.txt")) or ($force eq "1")){  

  $new_blasts_iter2++;
  $manager->start and next;
    print "::performing $alg searches -- new_sequences_iter1 vs ${genome}\n";


    if($alg eq "usearch"){
my $cmd="usearch8 -usearch_local new_sequences_iter1.faa -threads 1 -db ${genome}.udb -id $identity -blast6out ${genome}_blasted_iter2.txt >> log_file 2>&1";
	system($cmd);
    }

    elsif($alg eq "blast"){
   my  $cmd="blastp -query new_sequences_iter1.faa -db $genome -num_threads 1 -out ${genome}_blasted_iter2.txt -seg no -outfmt '6 qseqid sseqid pident length ppos'";
	system($cmd);    
    }





  $manager->finish;
    


    }

}
$manager->wait_all_children;


print "::analyzing the usearch output\n";

#now I can read the blast output and draw some conclusions

print "--new_blasts_iter2:=".$new_blasts_iter2."\n";

if($new_blasts_iter2 > 0){

    if(-e "all_blasted_iter2.txt"){`rm -rf all_blasted_iter2.txt`;}
    `cat *_blasted_iter2.txt > all_blasted_iter2.txt`;


    open(IN,"all_blasted_iter2.txt")||die "I cannot open all_blasted_iter2.txt";
    while(my $line=<IN>){
        chomp($line);
        #print $line."\n";
        my @data=split(/\t/,$line);

        my $query=$data[0];

        my $genome_q=$db_elements{$query}{"gen"};
        #print $genome_q."\n";
        $results{$query}{$genome_q}{$query}=1;

	my $length_query;

        if(exists($db_elements{$query})){
            $length_query=$db_elements{$query}{"len"};
        }
        else{
            print "::[ERROR] sequence $query not present in db of the hash %db_elements\n";
        }
    
        my $hit=$data[1];
    
        if($query eq $hit){next;}
    
	my $length_hit;
	my $genome_h;


        if(exists($db_elements{$hit})){
            $length_hit=$db_elements{$hit}{"len"};
            $genome_h=$db_elements{$hit}{"gen"};
        }
        else{
            print "::[ERROR] sequence $hit not present in db of the hash %db_elements\n";
        } 



        	if($sim eq "0"){
        $perc_identity=$data[2]; # identity
	}
	elsif($sim eq "1"){
	$perc_identity=$data[4]; # similarity
	}

        if($perc_identity<$identity_orthologs){next;}
        my $length_alignment=$data[3];
        
       # print $query."-".$hit."-"."$identity"."-".$length_alignment."-".($length_query*0.85)."-".($length_hit*0.85)."\n";

        if(($length_alignment >= ($length_query*0.85)) and ($length_alignment >= ($length_hit*0.85)) ){
            

            my $genome_h=$db_elements{$hit}{"gen"};
            if($genome_q eq $genome_h){
                if($perc_identity < $identity_paralogs){next;}
                
            }   
		
	    my $prev_hit;
	    my $prev_hit_id;

            
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

    my $ref=&LibFASTA::read_FASTA("new_sequences_iter1.faa"," ","1");
    my %new_sequences=%$ref;
    


    foreach my $entry (sort(keys(%new_sequences))){
        
        $db_elements{$entry}{"seen"}=1;

        #print OUT $entry."\n";
        my @arr_to_print=();


        foreach my $gen (@genomes){
            my $gen_id=$gen;


            if(exists($results{$entry}{$gen_id}{"hit"})){

                my $target=$results{$entry}{$gen_id}{"hit"};

                $db_elements{$target}{"seen"}=1;
    
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




`cat situation_iter1.txt > situation_all.txt`;
`cat situation_iter2.txt|tail -n+2 >> situation_all.txt`; 





my $graph = Graph::Undirected->new;


print "::Constructing the graph\n";

open(IN,"<situation_all.txt");
my $comment=<IN>;
chomp($comment);
my $header=$comment;
while(my $line=<IN>){
    chomp($line);
    if($line eq ""){next;}
    my @data=split(/\t/,$line); 
    my %hash_line=();
    foreach my $el (@data){
        my @el_fragmented=split(/,/,$el);
        foreach my $fragment(@el_fragmented){
            if($fragment ne "-"){        
            $hash_line{$fragment}=1;
            } 
        }   
    }
    my @res_path=keys(%hash_line);

    if($#res_path>0){
        $graph->add_path(keys(%hash_line));
    }
    if($#res_path eq 0){
        $graph->add_vertex($res_path[0]);
    }


}



close(IN);

#now I add all other nodes

foreach my $element (keys(%db_elements)){
    if(!($graph->has_vertex($element))){
        $graph->add_vertex($element);        
    }    
}


#print "The graph is $graph\n";
my $date=`date "+%Y-%m-%d %H:%M:%S"`;
chomp($date);
print "::perl is determining the connected components of the graph -- $date\n";
my @cc=$graph->connected_components();

$comment=~s/#//;
$comment=~s/.faa//g;
@genomes=split(/\t/,$comment);


my $date=`date "+%Y-%m-%d %H:%M:%S"`;
chomp($date);
print "::writing down the connected components of the graph (".($#cc+1).") -- ${date}\n";

open(OUT,">pre-table_linked5_strict.tsv");
print OUT $header."\n";


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



my $cmd="satv_untie-knots-paralogs.pl -in pre-table_linked5_strict.tsv -out ${file_out}";
system($cmd);



$date=`date "+%Y-%m-%d %H:%M:%S"`;
print "::Analysis completed at $date";
print "Bye!\n";

