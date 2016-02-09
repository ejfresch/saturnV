#!/usr/bin/perl
use Getopt::Long;

$file_in="";

GetOptions ("tab=s" => \$file_in,) or die("::usage: $0 -tab <file_in> \n[ERROR] launch failed! Please check the parameters!\n");


if($file_in eq ""){
    print "::usage: $0 -tab <file_in> \n[ERROR] launch failed! Please check the parameters!\n";
    exit();
}


@faa=`ls *.faa`;
chomp(@faa);


%list_ids=();

foreach $i (@faa){


    @ids=`cat $i |grep ">"`;
    chomp(@ids);

    foreach $id(@ids){
        $id=~s/>//;
        @explode=split(/ /,$id);
        $id_final=$explode[0];    
        $list_ids{$id_final}=0;    
    }   


}


open(IN,"<${file_in}");
while($line=<IN>){
    if($line=~/^#/){next};

    chomp($line);
    @explode=split(/[\t|,]/,$line);
    foreach $el (@explode){
        if($el eq "-"){next;}
        elsif(exists($list_ids{$el})){$list_ids{$el}=1;}
    
        else{
        print "::I did not find $el in the .faa files\n";
        }
    }

}

close(IN);



foreach $key(keys(%list_ids)){
    if($list_ids{$key} eq "0"){
        print "::I did not find $key in the table you provided\n";
    }

}

