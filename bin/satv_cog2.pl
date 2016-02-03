#!/usr/bin/perl

use Getopt::Long;



$base_path="/project/rclevesq/users/lfreschi/tasks/achromo/installation/bin/database";
my $sequences;
my $line;
my $line2;
my $line3;
my @data;
my @hit_result;
my @data_COG;
my @data_CAT;
my %hash_COG = ();
my %hash_CAT = ();
my %hash_PRODUCT = ();
my $COG;
my $CAT;
my $PRODUCT;
my $report;

GetOptions ("s=s" => \$sequences) or die("::usage: $0 -s the file containing the translated sequences\n");

open($report, '>', 'COG_report.tsv');

if (-e "$base_path/prot2003-2014.fa.udb")
	{ 

	print "## Database found \n";
	}
else
	{
	print "## Database not found, one will be created \n";
	system ("usearch8 -makeudb_usearch $base_path/prot2003-2014.fa -output $base_path/prot2003-2014.fa.udb");
	}

	print "## Performing comparison ...\n";
	system ("usearch8 -usearch_local $sequences -db $base_path/prot2003-2014.fa.udb -id 0.5 -blast6out result_blasted.txt > /dev/null 2>&1");
	print "## Done\n";

	print "## Will put in memory information regarding COG ... \n";
 open(IN,"$base_path/cog2003-2014.csv")||die "I cannot open cog2003-2014.csv";
    while($line=<IN>){
        chomp($line);

        @data_COG=split(/\,/,$line);
  
	
	$hash_COG{ $data_COG[0] } = $data_COG[6];

	}
close(IN);

 open(IN2,"$base_path/cognames2003-2014.tab")||die "I cannot open cognames2003-2014.tab";
    while($line2=<IN2>){
        chomp($line2);

        @data_CAT=split(/\t/,$line2);
  
	
	$hash_CAT{ $data_CAT[0] } = $data_CAT[1];
	$hash_PRODUCT{ $data_CAT[0] } = $data_CAT[2];

	}
close(IN2);
	print "## Done \n";
 open(IN3,"result_blasted.txt")||die "I cannot open result_blasted.txt";
    while($line3=<IN3>){
        chomp($line3);

        @data=split(/\t/,$line3);
   
	@hit_result=split(/\|/,$data[1]);
	

	$COG = "$hash_COG{$hit_result[1]}";

	$CAT = "$hash_CAT{$COG}";
	
	$PRODUCT = "$hash_PRODUCT{$COG}";

	print $report "$data[0]\t$COG\t$CAT\t$PRODUCT\n";
	}
close(IN3);
close($report);