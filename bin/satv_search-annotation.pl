#!/usr/bin/perl

use Getopt::Long;

$file_gff="";
$query="";

GetOptions ("gff=s" => \$file_gff,"query=s" => \$query) or die("::usage: $0 -gff <gff_file> -query <query>\n");


if(($file_gff eq "") or ($query eq "")){
    print "::usage: $0 -gff <gff_file> -query <query>\n";
    exit();
}




open(IN,"<${file_gff}") or die "::I cannot open the ${file_gff} file\n";
$comment=<IN>;

print "#ID\tprod\tcontig\tstart\tend\n";
while($line=<IN>){
chomp($line);

            if($line=~/^##/){next;}
            elsif($line=~/^##FASTA/){$switch_fasta=1;next;}
            elsif($switch_fasta eq "1"){last;}
           
            elsif($line=~/$query/i){

                   @explode_line=split(/\t/,$line);
                  $contig=$explode_line[0];
                  $start=$explode_line[3];
                  $end=$explode_line[4];

                   $line=~/ID=([0-9A-Za-z\_\- ]+);/;
                   $ID=$1;

                   $line=~/product=([0-9A-Za-z\_\- ]+)/;
                   $prod=$1;


                 
                print "$ID\t$prod\t$contig\t$start\t$end\n";







            }





            
        }

close(IN);











