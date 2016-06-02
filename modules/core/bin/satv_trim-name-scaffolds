#!/usr/bin/perl


if($#ARGV!=2){
    print "::usage: $0 <file_in> <file_out> <tag>\n";
}


$file_in=$ARGV[0];
$file_out=$ARGV[1];
$tag=$ARGV[2];



open(IN,"<$file_in");
open(OUT,">$file_out");

$count=1;

while($line=<IN>){
    chomp($line);
    if($line=~/^>/){
      print OUT  ">${tag}_$count\n";   
      $count++;                   
    
    }
    else{print OUT $line."\n";}



}



close(IN);
close(OUT);
