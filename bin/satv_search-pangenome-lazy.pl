#!/usr/bin/perl -I /home/avincent/Desktop/saturnV/bin


use strict;
use LibFASTA;
use Parallel::ForkManager;
use Getopt::Long;


my $genome_list="";
my $n_cpu=1;
my $identity=90;
my $force=0;
my $alg="usearch";
my $file_out="table_linked5_lazy.tsv";
my $sim=0;
my $identity_ali;

#here are the available algorithms for the search
my %avail_algs=(
    "usearch" => 1,
    "blast" => 1
);


GetOptions ("g=s" => \$genome_list,"out=s"   => \$file_out,"c=s"   => \$n_cpu,"i=s"   => \$identity,"f=s"   => \$force,"a=s"=> \$alg, "sim=s"=> \$sim) or die("::usage: $0 -d <genomes_list> -c <n_cpu> -i <perc_identity> -f <force:[0|1]> -a <algorithm>\n");


if($genome_list eq ""){

    print "::usage: $0 -g <genomes_list> -out <out_file> -c <n_cpu> -i <perc_identity> -f <force:[0|1]> -a <algorithm> -sim <similarity:[0|1]\n";
    exit();
}



if(!(exists($avail_algs{$alg}))){
    print "::usage: $0 -g <genomes_list> -out <out_file> -c <n_cpu> -i <perc_identity> -f <force:[0|1]> -a <algorithm> -sim <similarity:[0|1]\n";
    exit();    
}



my $identity_usearch=$identity/100;


my $manager = new Parallel::ForkManager($n_cpu);

my $date=`date "+%Y-%m-%d %H:%M:%S"`;
print "::Starting analysis at $date";


my @genomes=`cat $genome_list`;
chomp(@genomes);


#I generate the db
for my $genome (@genomes){


    
    if(((!(-e "${genome}.udb")) or ($force eq "1"))and ($alg eq "usearch")){
    print "::indexing sequences -- ${genome}\n";
    `usearch8 -makeudb_usearch $genome -output ${genome}.udb >> log_file 2>&1`;
    }

    elsif(((!(-e "${genome}.pin")) or ($force eq "1")) and ($alg eq "blast")){
    print "::indexing sequences -- ${genome}\n";
    `makeblastdb -in $genome -dbtype prot > /dev/null`;
    }



}




my $first_genome=$genomes[0];

#I perform the usearch searches


my $new_blasts=0;

for my $genome (@genomes[0..$#genomes]){
    
  if(((!(-e "${genome}_blasted.txt")) or ($force eq "1")) and ($alg eq "usearch")){  
    
   $new_blasts++;  
  
  $manager->start and next;
    print "::performing $alg searches -- ${first_genome} vs ${genome}\n";
    `usearch8 -usearch_local $first_genome -threads 1 -db ${genome}.udb -id $identity_usearch -blast6out ${genome}_blasted.txt >> log_file 2>&1`;

  $manager->finish;

    

  }  
  elsif(((!(-e "${genome}_blasted.txt")) or ($force eq "1")) and ($alg eq "blast")){  
    
   $new_blasts++;  
  
  $manager->start and next;
    print "::performing $alg searches -- ${first_genome} vs ${genome}\n";
    `blastp -query $first_genome -db $genome -num_threads 1 -out ${genome}_blasted.txt -seg no -outfmt '6 qseqid sseqid pident length ppos'`;

  $manager->finish;

    

  }  



}


$manager->wait_all_children;




my %db_all_seq=();

#I build the hash with the sequences
print "::building the hash with the sequences\n";

my @ids_first_genome=();

my $ref=&LibFASTA::read_FASTA("$genomes[0]"," ","1");
my %current_hash=%$ref;

    foreach my $elem (keys(%current_hash)){
        
        if(exists($db_all_seq{$elem})){print "::[ERROR] Hey! There are two or more sequences with the same ID -- ID: $elem !\n";}

        $db_all_seq{$elem}{"seq"}=$current_hash{$elem};
        $db_all_seq{$elem}{"gen"}=$genomes[0];
        push(@ids_first_genome,$elem);

    }



for my $genome (@genomes[1..$#genomes]){

my $ref=&LibFASTA::read_FASTA("$genome"," ","1");
my %current_hash=%$ref;

    foreach my $elem (keys(%current_hash)){

        $db_all_seq{$elem}{"seq"}=$current_hash{$elem};
        $db_all_seq{$elem}{"gen"}=$genome;


    }



}

print "::analyzing the usearch output\n";

#now I can read the usearch output and draw some conclusions
print "--new_blasts_iter1:=".$new_blasts."\n";

if($new_blasts > 0){


    if(-e "all_blasted.txt"){`rm -rf all_blasted.txt`;}
    `cat *_blasted.txt > all_blasted.txt`;



    #I initialize the hash of the results
    my %results=();

    #sequences for which I already found a match and therefore can be excluded from further analysis
    my %exclusion_zone=();

    open(IN,"all_blasted.txt")||die "I cannot open all_blasted.txt";
    while(my $line=<IN>){
        chomp($line);
        #print $line."\n";
        my @data=split(/\t/,$line);

        my $query=$data[0];
        $exclusion_zone{$query}=1;

       my $genome_q=$db_all_seq{$query}{"gen"};
        #print $genome_q."\n";
        $results{$query}{$genome_q}{$query}=1;

	my $length_query;

        if(exists($db_all_seq{$query})){
            $length_query=length($db_all_seq{$query}{"seq"});
        }
        else{
            print "::[ERROR] sequence $query not present in db of the hash %db_all_seq\n";
        }
        
        my $hit=$data[1];
        
	my $length_hit;
	my $genome_id; #it is the genome of the hit

        
        if(exists($db_all_seq{$hit})){
            $length_hit=length($db_all_seq{$hit}{"seq"});
            $genome_id=$db_all_seq{$hit}{"gen"};
            
        }
        else{
            print "::[ERROR] sequence $hit not present in db of the hash %db_all_seq\n";
        } 


                             	if($sim eq "0"){
        $identity_ali=$data[2]; # identity
	}
	elsif($sim eq "1"){
	$identity_ali=$data[4]; # similarity
	}

        if($identity_ali < $identity){next;}

        my $length_alignment=$data[3];
        
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
    foreach my $entry (@ids_first_genome){


         my@arr_to_print=();
         push(@arr_to_print,$entry);

        foreach my $gen (@genomes[1..$#genomes]){
            my $gen_id=$gen;




            if(exists($results{$entry}{$gen_id})){
            
            my $ref=$results{$entry}{$gen_id};
            my %hash_ref=%$ref;


            my@hits=sort(keys(%hash_ref));
            push(@arr_to_print,join(",",@hits));

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
foreach my $seq (keys(%db_all_seq)){
    
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
my %results=();



#I build the db

if((!(-e "new_sequences_iter1.faa.udb")) or ($force eq "1")){
    `usearch8 -makeudb_usearch new_sequences_iter1.faa -output new_sequences_iter1.faa.udb >> log_file 2>&1`;


}

elsif(((!(-e "new_sequences_iter1.faa.pin")) or ($force eq "1")) and ($alg eq "blast")){

    `makeblastdb -in new_sequences_iter1.faa -dbtype prot > /dev/null`;
}




# I run the comparisons

my $new_blasts_iter2=0;

for my $genome (@genomes){

  if(((!(-e "${genome}_blasted_iter2.txt")) or ($force eq "1"))and ($alg eq "usearch")){  

  $new_blasts_iter2++;
  $manager->start and next;
    print "::performing $alg searches -- new_sequences_iter1 vs ${genome}\n";

    `usearch8 -usearch_local new_sequences_iter1.faa -threads 1 -db ${genome}.udb -id $identity_usearch -blast6out ${genome}_blasted_iter2.txt >> log_file 2>&1`;

  $manager->finish;
    


    }


  elsif(((!(-e "${genome}_blasted_iter2.txt")) or ($force eq "1"))and ($alg eq "blast")){  

  $new_blasts_iter2++;
  $manager->start and next;
    print "::performing $alg searches -- new_sequences_iter1 vs ${genome}\n";

    `blastp -query new_sequences_iter1.faa -db $genome -num_threads 1 -out ${genome}_blasted_iter2.txt -seg no -outfmt '6 qseqid sseqid pident length ppos'`;

 

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

        my $genome_q=$db_all_seq{$query}{"gen"};
        #print $genome_q."\n";
        $results{$query}{$genome_q}{$query}=1;


	my $length_query;

        if(exists($db_all_seq{$query})){
            $length_query=length($db_all_seq{$query}{"seq"});
        }
        else{
            print "::[ERROR] sequence $query not present in db of the hash %db_all_seq\n";
        }
    
        my $hit=$data[1];
    
      
  	my $length_hit; 
  	my $genome_id; 

 
        if(exists($db_all_seq{$hit})){
            $length_hit=length($db_all_seq{$hit}{"seq"});
            $genome_id=$db_all_seq{$hit}{"gen"};
        }
        else{
            print "::[ERROR] sequence $hit not present in db of the hash %db_all_seq\n";
        } 

                       	if($sim eq "0"){
        $identity_ali=$data[2]; # identity
	}
	elsif($sim eq "1"){
	$identity_ali=$data[4]; # similarity
	}

        if($identity_ali < $identity){next;}



        my $length_alignment=$data[3];


        if(($length_alignment >= ($length_query*0.85)) and ($length_alignment >= ($length_hit*0.85)) ){
        
            $results{$query}{$genome_id}{$hit}=1; 

    
        }
    


    }
    close(IN);


    open(OUT,">situation_iter2.txt");
    print OUT "#".join("\t",@genomes)."\n";

    my $ref=&LibFASTA::read_FASTA("new_sequences_iter1.faa"," ","1");
    my %new_sequences=%$ref;
    


    foreach my $entry (sort(keys(%new_sequences))){


        #print OUT $entry."\n";
        my @arr_to_print=();


        foreach my $gen (@genomes){
            my $gen_id=$gen;


            if(exists($results{$entry}{$gen_id})){
            my $ref_r=$results{$entry}{$gen_id};
            my %hash_ref=%$ref_r;


            my @hits=sort(keys(%hash_ref));
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
    print "::I have already made all the searches of the second iteration. Everything is up-to-date! :)\n";
}

print "::Second iteration finished!\n";
$date=`date "+%Y-%m-%d %H:%M:%S"`;


 my $cmd='cp situation_iter1.txt situation_all.txt';
 system($cmd);

 my $cmd='cat situation_iter2.txt |grep -v  "#" >> situation_all.txt';
 system($cmd);



my $cmd="satv_merge-data3.pl -in situation_all.txt -out ${file_out}";
system($cmd);





print "::Analysis completed at $date";
print "Bye!\n";
