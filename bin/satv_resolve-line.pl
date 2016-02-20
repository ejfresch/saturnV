#!/usr/bin/perl -I /home/lfreschi/tasks/achromo/installation/bin
  
use Getopt::Long;
use LibFASTA;
use Graph;
use strict;


my $out_file="";
my $identity_orthologs=90;
my $identity_paralogs=90;
my $line="";
my $id_line="";
my $genome_list="";


GetOptions ("g=s" => \$genome_list,"out=s" => \$out_file,"i=s"   => \$identity_orthologs,"ip=s"   => \$identity_paralogs,"line=s" => \$line,"id=s" => \$id_line) or die("::usage: $0 -g <genomes_list> -out <out_file> -i <identity_usearch> -ip <identity_paralogs -line <line>\n[ERROR] launch failed! Please check the parameters!\n");

if( ($out_file eq "") or ($line eq "") or ($id_line eq "")){

    print "::usage: $0 -g <genomes_list> -out <out_file> -i <identity_usearch> -ip <identity_paralogs> -line <line>\n";
    exit();
}



my @genomes=`cat $genome_list`;
chomp(@genomes);

$identity_orthologs=$identity_orthologs/100;


$line=~s/&/\t/g;

print "::analyzing line $id_line\n";
my @elements=split(/\t/,$line);

my %ids=();

my %genomes_to_scan=();

for(my $i=0;$i<=$#elements;$i++){
    
	if($elements[$i] eq "-"){next;}
	my $el=$elements[$i];
    my @all=split(/,/,$el); 
    foreach my $id (@all){
         $ids{$id}=1;
        my @explode_id=split(/_/,$id);
        my $current_genome=$explode_id[0];
        $genomes_to_scan{$current_genome}{$id}=1;
            
    }

}
    
#I need to retrieve the sequences
my %db_elements=();
foreach my $key (keys(%genomes_to_scan)){
            
	my $ref=&LibFASTA::read_FASTA("${key}.faa"," ","1");
    my %current_hash=%$ref;
     
	my $r2=$genomes_to_scan{$key};
	my %seq_ids=%$r2;
        
    foreach my $s (keys(%seq_ids)){

    	if(exists($current_hash{$s})){
            my $current_seq=$current_hash{$s};
        	$current_seq=~s/\*//;

            $db_elements{$s}{"seq"}=$current_seq;
            $db_elements{$s}{"gen"}=$key;
            $db_elements{$s}{"len"}=length($current_seq);


        }
        else{print ":: [WARNING] I did not find a sequence with the following ID ($s) when I read the file ${key}.faa\n";}

    }    

}



#Now I create a fasta file with all sequences. Then I create a copy of it and I blast or userach the two files one against the other
  
open(SEQ,">fasta_query_${id_line}.faa");
      foreach my $q (keys(%db_elements)){
            print SEQ ">".$q."\n";
            print SEQ $db_elements{$q}{"seq"}."\n";
      }
close(SEQ);

my $cmd="cp fasta_query_${id_line}.faa fasta_subject_${id_line}.faa";
system($cmd);

print "::performing sequence comparisons with usearch\n";    
    
my $cmd="usearch8 -makeudb_usearch fasta_subject_${id_line}.faa -output fasta_subject_${id_line}.faa.udb >> log_file 2>&1";
system($cmd);

my $cmd="usearch8 -usearch_local fasta_query_${id_line}.faa -threads 1 -db fasta_subject_${id_line}.faa.udb -id $identity_orthologs -maxaccepts 0 -maxrejects 0 -blast6out usearch_resolve_line-${id_line}.txt >> log_file 2>&1";
#print $cmd."\n";
system($cmd);


my %results=();

open(SCH,"<usearch_resolve_line-${id_line}.txt");
while(my $line_res=<SCH>){
	chomp($line_res);

	my @data=split(/\t/,$line_res);
    my $query=$data[0];
	my $hit=$data[1];
	#print "".$line_res."\n";
   
    if($query eq $hit){#print "IF--".$line_res."\n";
		next;
	}



	my $length_query;
	my $genome_q;
	my $length_hit;
	my $genome_h;
	

    
    #print "PASS--".$line_res."\n";
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


    my $perc_identity=$data[2];
    my $length_alignment=$data[3];
    
    #print $query."\t".$hit."\t".$perc_identity."\t".$length_alignment."\t".($length_query*0.85)."\t".($length_hit*0.85)."\n";   
	#print $query."-".$hit."-"."$identity"."-".$length_alignment."-".($length_query*0.85)."-".($length_hit*0.85)."\n";

    if(($length_alignment >= ($length_query*0.85)) and ($length_alignment >= ($length_hit*0.85)) ){
      
		#print $query."\t".$hit."\t".$perc_identity."\t".$length_alignment."\t".($length_query*0.85)."\t".($length_hit*0.85)."\n";
 
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
        }
        else{            
                $results{$query}{$genome_h}{"hit"}=$hit;
                $results{$query}{$genome_h}{"id"}=$perc_identity;
        }
        
    }
        

}


close(SCH);


   


#I write all data
open(FTR,">usearch_untie_knots_paralogs_db_${id_line}.txt");
foreach my $q (keys(%results)){
	if(!(exists($results{$q}))){next;}
    my $ref=$results{$q};
	my %hash=%$ref;
    foreach my $k (keys(%hash)){
    	print FTR "$q\t".$results{$q}{$k}{"hit"}."\n";
	}
}
close(FTR);

#print "::Exiting!\n";
#exit();


#I create the graph
my $g = Graph->new;
open(FTR,"<usearch_untie_knots_paralogs_db_${id_line}.txt");
while(my $line_g=<FTR>){
	chomp($line_g);
    if($line_g eq ""){next;}
    my @el_fragmented=split(/\t/,$line_g);
    #print join("^",@el_fragmented)."\n";        
    $g->add_path(@el_fragmented);
}
close(FTR);

foreach my $element (keys(%db_elements)){
    if(!($g->has_vertex($element))){
        $g->add_vertex($element);        
    }    
}


print "::perl is determining the strongly connected components of the graph\n";
my @cc=$g->strongly_connected_components();

print "::writing down the strongly connected components of the graph (".($#cc+1).")\n";

#open(D,">>table_linked3_mod.tsv");
#print D join("\t",@genomes)."\n";

open(OUT_FILE,">usearch_untie_knots_paralogs_table_${id_line}.txt");

my $count=0;
foreach my $c (@cc){
	my $count++;
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
	    #print ${id_line}.">>".join("\t",@final_arr)."\n";
        print OUT_FILE join("\t",@final_arr)."\n";

}



close(OUT_FILE);



