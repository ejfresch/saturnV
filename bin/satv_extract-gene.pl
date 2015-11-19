#!/usr/bin/perl 


use Getopt::Long;

GetOptions ("tab=s" => \$table, "dir=s" => \$directory) or die("::usage: $0 -tab the table -dir the directory containing the .fna files\n");

mkdir 'extracted_sequences';
open(IN,"$table")||die "I cannot open table_linked3.tsv";
   	 while($line=<IN>){
     	   chomp($line);
	$i++;
	print "line = $i\n";
	if ($line =~ /^\#/){$i--; next}
	

	@line_information=split(/\t/,$line);
	

	open(OUT,">extracted_sequences/seq_$i.fasta");
	foreach(@line_information)
		{
	if ($_ eq '-'){next}
	@file_name=split(/_/,$_);
	$file = "$directory/$file_name[0].fna";
	$pattern = "$_";
	$data = `cat $file`;
	($fa) = $data =~ /(>$pattern[^>]+)/s;
	($fa) =~s/_$file_name[1]//;
	print OUT $fa;

		}
	close (OUT);



	}
close (IN);
