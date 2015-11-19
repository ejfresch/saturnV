#!/usr/bin/perl -I /home/lfreschi/tasks/pangenome/scripts/

#I should read the fasta file and for each sequence make an hash of hashes $hash{sequence}{id}

use LibFASTA;
use Getopt::Long;

GetOptions ("l=s" => \$file_in) or die("::usage: $0 -l <list_faa_files>\n");

%hash_seq=();

open(IN,"<$file_in");
while($line=<IN>){
    
    chomp($line);
    print "::reading proteins of $line\n";
    $ref=&LibFASTA::read_FASTA("$line"," ","1");
    %current_hash=%$ref;

    foreach $el (keys(%current_hash)){
        $seq=$current_hash{$el};
        $hash_seq{$seq}{$el}=1;
        #print $seq."=>".$el."\n";

    }




}

close(IN);

print "::printing the elements > 1 copy\n";

foreach $elem (keys(%hash_seq)){
   # print $key."\n";
    $ref=$hash_seq{$elem};
    #print $ref."\n";    
    %hash=%$ref;
    
    @arr=keys(%hash);
    if($#arr>0){
        print join("\t",@arr)."\n";
    }

}









