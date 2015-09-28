#!/usr/bin/perl

$col=2;

`mkdir out/`;

open(IN,"$ARGV[0]") or die("::I cannot open the file $ARGV[0]\n");
$comment=<IN>;
while($line=<IN>){
    chomp($line);
    @data=split(/\t/,$line);
    $source=$data[0];
    $dest=$data[1].".fasta";    
    `cp ${source} out/${dest}`;    

}

close(IN);

