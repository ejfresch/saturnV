#!/usr/bin/perl

use Getopt::Long;

$file_gff="";
$out_dir=".";

GetOptions ("gff=s" => \$file_gff,"out=s" => \$out_dir) or die("::usage: $0 -gff <gff_file> -out <out_dir>\n");



if(($file_gff eq "")){
    print "::usage: $0 -gff <gff_file> -out <out_dir>\n";
    exit();
}


$strain=`basename ${file_gff}`;
chomp($strain);
$strain=~s/.gff//;


open(IN,"<${file_gff}") or die "::I cannot open the ${file_gff} file\n";
$comment=<IN>;


%fragments=();

while($line=<IN>){
    chomp($line);

    $gene="";
   
    if($line=~/##FASTA/){last;}
    else{
        @data=split(/\t/,$line);
        $fragment_name=$data[0];
        #print "FRAG_NAME=".$fragment_name."\n";
 

        $strand=$data[6];
        $details=$data[8];
        @details_exploded=split(/;/,$details);
        

        for($j=0;$j<=$#details_exploded;$j++){
        
            if($details_exploded[$j]=~/^gene=/){
                $gene=$details_exploded[$j];
                $gene=~s/gene=//g;
            }
            

        
        }

        if($gene eq ""){next;}

        if(!(exists($fragments{$fragment_name}))){
            $fragments{$fragment_name}=$strand.$gene;
        }else{
            $fragments{$fragment_name}.="\t".$strand.$gene;
        }       
             
    }

}



close(IN);

open(OUT,">${out_dir}/${strain}.ordered");

foreach $el (keys(%fragments)){
    print OUT $fragments{$el}."\n";

}
close(OUT);

