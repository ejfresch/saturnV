#!/usr/bin/perl

use Getopt::Long;



$file_in="";
$file_out="";

GetOptions ("tab=s" => \$file_in, "out=s" => \$file_out) or die("::usage: $0 -tab <table_genes> -out <file_out>\n");



if(($file_in eq "") or ($file_out eq "")){
    print "::usage: $0 -tab <table_genes> -out <file_out>\n";
    exit();
}



open(IN,"<${file_in}");
@data=<IN>;
chomp(@data);
close(IN);


open(OUT,">${file_out}");
$first_line=shift(@data);
$first_line=~s/^#//;
print OUT "ID\t".$first_line."\n";

foreach $line (@data){
    @explode=split(/\s+/,$line);
    
    $flag=0;
    @results=();    
    $id="";    
    
    foreach $entry (@explode){
        
        if(($flag eq "0") and ($entry ne "-")){
            
            @explode_entry=split(/,/,$entry);
            $id=$explode_entry[0];
            $flag=1;            

        }
        
        if($entry eq "-"){
            push(@results,"0");
        }else{
            push(@results,"1");
        }



    }
    print OUT $id."\t";
    print OUT join("\t",@results)."\n";

}


close(OUT);
