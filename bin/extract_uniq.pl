#!/usr/bin/perl

#convert the file to a numerical matrix

use Getopt::Long;

#::usage $0 -d <dir_genomes> -c <n_cpu>

$file_in="";
$ref_strain="";

GetOptions ("in=s" => \$file_in) or die("::usage: $0 -in <table_genes>\n");



if(($file_in eq "")){
    print "::usage: $0 -in <table_genes>\n";
    exit();
}



print "::searching for uniq genes\n";

open(IN,"<$file_in");
$comment=<IN>;
$comment=~s/^#//;

@explode=split(/\t/,$comment);

open(OUT,">uniq_genes.txt");


%uniq_genes=();
$uniq_groups=0;

while($line=<IN>){
    @row=();
    
    chomp($line);
    @data=split(/\t/,$line);

    $count_ones=0;
    
    for($i=0;$i<=$#data;$i++){
        if($data[$i] eq "-"){
        
         }
        else{
            $count_ones++;
            push(@row,$explode[$i]);
            
            push(@row,$data[$i]);
            
            @explode_id=split(/,/,$data[$i]);
  
            
        }
    }
    
    if($count_ones eq 1){
        print OUT join(":",@row)."\n";
        
        foreach $gene (@explode_id){
            $uniq_genes{$gene}=1;
        }

    }
    


}

close(IN);
close(OUT);

print "::There are ${uniq_groups} uniq gene groups\n";
@uniq_genes=keys(%uniq_genes);



print "::There are ".($#core_genes+1)." unique genes\n";

