#!/usr/bin/perl

use Getopt::Long;

$file_gff="";
$coord_start="";
$coord_end="";

GetOptions ("gff=s" => \$file_gff,"start=s" => \$coord_start,"end=s" => \$coord_end) or die("::usage: $0 -gff <gff_file> -start <coord_start> -end <coord_end>\n");


if(($file_gff eq "") or ($coord_start eq "") or ($coord_end eq "")){
    print "::usage: $0 -gff <gff_file> -start <coord_start> -end <coord_end>\n";
    exit();
}

        print "#your coordinates -- start: ${coord_start}; end: ${coord_end}\n";
        %results=();

       $contig_name=""; 


        open(IN,"<${file_gff}") or die "::I cannot open the ${gff} file\n";
        $comment=<IN>;

        while($line=<IN>){
            chomp($line);

            
           
            if($line=~/ID=([0-9A-Za-z\_\-]+);/){
                  $current_id=$1;
                    
                  @explode_line=split(/\t/,$line);
                  $contig=$explode_line[0];
                  $start=$explode_line[3];
                  $end=$explode_line[4];
                  

                  if(($coord_start<=$start) and (($coord_end>=$end))){$str=$current_id."\t$contig\t$start\t$end\n";$results{$str}=1;}
                  if(($coord_start<=$start) and (($coord_end<=$end) and ($coord_end>=$start))){$str=$current_id."\t$contig\t$start\t$end\n";$results{$str}=1;}
                  if((($coord_start>=$start) and ($coord_start<=$end)) and ($coord_end>=$end)){$str=$current_id."\t$contig\t$start\t$end\n";$results{$str}=1;}
  
            }

            




            
        }


close(IN);
       

foreach $key (keys(%results)){
    print $key;

}
    
    
















