#!/usr/bin/perl -I /project/rclevesq/users/lfreschi/tasks/achromo/installation/bin

use Getopt::Long;
use LibFASTA;
use Graph;

$input_file="";
$out_file="";
$n_cpu=1;
$identity_usearch=90;
$identity_paralogs=90;



GetOptions ("in=s" => \$input_file,"out=s" => \$out_file,"c=s"   => \$n_cpu,"i=s"   => \$identity_usearch,"ip=s"   => \$identity_paralogs) or die("::usage: $0 -in <input_file> -c <n_cpu> -i <identity_usearch> -ip <identity_paralogs>\n[ERROR] launch failed! Please check the parameters!\n");

if(($input_file eq "") or ($out_file eq "")){

    print "::usage: $0 -in <input_file> -c <n_cpu> -i <identity_usearch> -ip <identity_paralogs> -out <out_file>\n";
    exit();
}


$identity_usearch=$identity_usearch/100;


open(IN_FILE,"<$input_file");
open(OUT_FILE,">$out_file");

$line=<IN_FILE>;
chomp($line);
print OUT_FILE $line."\n";
$line=~s/#//;
@genomes=split(/\t/,$line);

$count_overall=0;

while($line=<IN_FILE>){
chomp($line);
    
    $count_overall++;


    if($line=~/,/){
    
        print "::analyzing line $count_overall\n";
        @elements=split(/\t/,$line);

        %ids=();

        %genomes_to_scan=();

        for($i=0;$i<=$#elements;$i++){
    
            if($elements[$i] eq "-"){next;}
            $el=$elements[$i];
            @all=split(/,/,$el); 
            foreach $id (@all){
                $ids{$id}=1;

                @explode_id=split(/_/,$id);
                $current_genome=$explode_id[0];
                $genomes_to_scan{$current_genome}{$id}=1;
            
            }

        }
    
        #I need to retrieve the sequences
        %db_elements=();
        foreach $key (keys(%genomes_to_scan)){
            
            $ref=&LibFASTA::read_FASTA("${key}.faa"," ","1");
            %current_hash=%$ref;

            
            $r2=$genomes_to_scan{$key};
            %seq_ids=%$r2;
        
            foreach $s (keys(%seq_ids)){

                if(exists($current_hash{$s})){
                    $current_seq=$current_hash{$s};
                    $current_seq=~s/\*//;

                    $db_elements{$s}{"seq"}=$current_seq;
                    $db_elements{$s}{"gen"}=$key;
                    $db_elements{$s}{"len"}=length($current_seq);


                }
                else{print ":: [WARNING] I did not find a sequence with the following ID ($s) when I read the file ${key}.faa\n";}

            }    

      }

      #Now I create a fasta file with all sequences. Then I create a copy of it and I blast or userach the two files one against the other
  
      open(SEQ,">fasta_query.faa");
      foreach $q (keys(%db_elements)){
            print SEQ ">".$q."\n";
            print SEQ $db_elements{$q}{"seq"}."\n";

      }
      close(SEQ);

      $cmd="cp fasta_query.faa fasta_subject.faa";
      system($cmd);


      print "::performing sequence comparisons with usearch\n";    
      #exit();

    
      $cmd="usearch8 -makeudb_usearch fasta_subject.faa -output fasta_subject.faa.udb >> log_file 2>&1";
      system($cmd);

      $cmd="usearch8 -usearch_local fasta_query.faa -threads 1 -db fasta_subject.faa.udb -id $identity_usearch -blast6out usearch_untie_knots_paralogs_partial_results.txt >> log_file 2>&1";
      system($cmd);
    

    %results=();

    open(SCH,"<usearch_untie_knots_paralogs_partial_results.txt");
    while($line_res=<SCH>){
        chomp($line_res);
        
        @data=split(/\t/,$line_res);
    
        $query=$data[0];

        
        $hit=$data[1];

        if($query eq $hit){next;}
    
        if(exists($db_elements{$query})){
            $length_query=$db_elements{$query}{"len"};
            $genome_q=$db_elements{$query}{"gen"};
        }
        else{print "::[ERROR] sequence $query not present in db of the hash db_elements\n";} 

        if(exists($db_elements{$hit})){
            $length_hit=$db_elements{$hit}{"len"};
            $genome_h=$db_elements{$hit}{"gen"};
        }
        else{print "::[ERROR] sequence $hit not present in db of the hash db_elements\n";} 

        $perc_identity=$data[2];
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


                }            
                


            }else{            

                $results{$query}{$genome_h}{"hit"}=$hit;
                $results{$query}{$genome_h}{"id"}=$perc_identity;

            }
        
        }
        

    }


    close(SCH);


   


    #I write all data
    open(FTR,">usearch_untie_knots_paralogs_db.txt");
    foreach $q (keys(%results)){

        $ref=$results{$q};

        %hash=%$ref;

        foreach $k (keys(%hash)){
            print FTR "$q\t".$results{$q}{$k}{"hit"}."\n";

        }



    }

    close(FTR);

    #print "::Exiting!\n";
    #exit();


    #I create the graph
    my $g = Graph->new;


    open(FTR,"<usearch_untie_knots_paralogs_db.txt");

    while($line_g=<FTR>){
        chomp($line_g);
        if($line_g eq ""){next;}
        @el_fragmented=split(/\t/,$line_g);
        #print join("^",@el_fragmented)."\n";        
        $g->add_path(@el_fragmented);

    }


    close(FTR);


    print "::perl is determining the strongly connected components of the graph\n";
    @cc=$g->strongly_connected_components();

    
    print "::writing down the strongly connected components of the graph (".($#cc+1).")\n";

    #open(D,">>table_linked3_mod.tsv");
    #print D join("\t",@genomes)."\n";

    $count=0;
    foreach $c (@cc){
        $count++;



        %hash_results=();
        @current_cc=@$c;
        foreach $element (@current_cc){

            if($element eq "-"){next;}
            #print $element."\n";
            @explode=split(/_/,$element);
            #print $explode[0]."\n";
            $hash_results{$explode[0]}{$element}=1;        


         
        }

        @final_arr=();    
        
        foreach $genome (@genomes){

        $genome=~s/.faa//g;    
        #print $genome."\n";
        if(exists($hash_results{$genome})){
           
            $ref=$hash_results{$genome};

            %arr=%$ref;
            push(@final_arr,join(",",sort(keys(%arr))));
        }
        else{push(@final_arr,"-");}

        }

        print OUT_FILE join("\t",@final_arr)."\n";

    }

    #close(D);


    




    }



    else{print OUT_FILE $line."\n";}



}

close(IN_FILE);
close(OUT_FILE);

$cmd="rm -rf usearch_untie_knots_paralogs_partial_*";
system($cmd);
