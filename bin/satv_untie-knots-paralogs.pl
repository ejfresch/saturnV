#!/usr/bin/perl -I /project/rclevesq/users/lfreschi/tasks/pangenome/saturnv/bin

use Getopt::Long;
use LibFASTA;

$input_file="";
$out_file="";
$n_cpu=1;
$identity_usearch=90;
$identity_paralogs=90;

GetOptions ("in=s" => \$input_file,"out=s" => \$out_file,"c=s"   => \$n_cpu,"i=s"   => \$identity_usearch,"ip=s"   => \$identity_paralogs) or die("::usage: $0 -in <input_file> -c <n_cpu> -i <identity_usearch> -ip <identity_paralogs>\n[ERROR] launch failed! Please check the parameters!\n");

if($input_file eq ""){

    print "::usage: $0 -in <input_file> -c <n_cpu> -i <identity_usearch> -ip <identity_paralogs> -out <out_file>\n";
    exit();
}



open(IN,"<$input_file");
open(OUT,">$out_file");

$line=<IN>;
chomp($line);
print OUT $line."\n";
$line=~s/#//;
@genomes=split(/\t/,$line);

while($line=<IN>){
chomp($line);
    
@elements=split(/\t/,$line);

$counter=0;
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

    $count=$#all+1;
    $counter+=$count;
}

$expected_sum=$#genomes+1;

if($counter eq $expected_sum){
    print OUT $line;

}
else{
    
    #I need to retrieve the sequences
    %seq=();
    foreach $key (keys(%genomes_to_scan)){
        
        $ref=&LibFASTA::read_FASTA("${key}.faa"," ","1");
        %current_hash=%$ref;

        
        $r2=$genomes_to_scan{$key};
        %seq_ids=%$r2;
    
        foreach $s (keys(%seq_ids)){

            if(exists($current_hash{$s})){
                $current_seq=$current_hash{$s};
                $seq{$s}=$current_seq;
            }else{
                print ":: [WARNING] I did not find a sequence with the following ID ($s) when I read the file ${key}.faa\n";
            }


        }    



    }

    #Now I create a fasta file with all sequences. Then I create a copy of it and I blast or userach the two files one against the other
  
    open(OUT,">fasta_query.faa");
    foreach $q (keys(%seq)){
        print $q."\n";
        print $seq{$q}."\n";

    }
    close(OUT);

    $cmd="cp fasta_query.faa fasta_subject.faa";
    system($cmd);

    
    $cmd="usearch8 -makeudb_usearch fasta_subject.faa -output fasta_subject.faa.udb >> log_file 2>&1";
    system($cmd);

    $cmd="usearch8 -usearch_local fasta_query.faa -threads 1 -db fasta_subject.faa.udb -id $identity_usearch -blast6out .txt >> log_file 2>&1";
    system($cmd);
    
    





}




}

close(IN);
close(OUT);
