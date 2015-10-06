#!/usr/bin/perl

use Getopt::Long;

$file_in="";
$include="";
$exclude="";
$file_out="";

GetOptions ("tab=s" => \$file_in,"out=s" => \$file_out,"include=s" => \$include, "exclude=s" => \$exclude) or die("::usage: $0 -tab <table_genes> -out <file_out> -include <file_include> -exclude <file_exclude>\n");


if(($file_in eq "") or ($file_out eq "") or (($include eq "") and ($exclude eq ""))){
    print "::usage: $0 -tab <table_genes> -out <file_out> -include <file_include> -exclude <file_exclude>\n";
    exit();
}


if(($include ne "") and ($exclude ne "")){
    print "::usage: $0 -tab <table_genes> -include <file_include> -exclude <file_exclude>\n";
    print "::you cannot use the -include and -exclude option at the same time. It does not make sense.\n";
    exit();
}





$flag_exclude=0;
$flag_include=0;

@include_gen=();
@exclude_gen=();


if($include ne ""){
    
    $flag_include=1;

    open(IN,"<$include")or die "::I cannot open ${include} \n";
    while($line=<IN>){
        chomp($line);
        push(@include_gen,$line);
    }

    close(IN);


}

if($exclude ne ""){
    
    $flag_exclude=1;

    open(IN,"<$exclude")or die "::I cannot open ${exclude} \n";
    while($line=<IN>){
        chomp($line);
        push(@exclude_gen,$line);
    }

    close(IN);


}


@position_genomes=();

open(IN,"<$file_in") or die "::I cannot open the $file_in file\n";

$header=<IN>;
chomp($header);
$header=~s/#//;

@data=split(/\t/,$header);

%hash_genomes=();

for($i=0;$i<$#data;$i++){
    
    $gen=$data[$i];    
    $hash_genomes{$gen}=$i;    

}



if($flag_include eq "1"){

    foreach $el (@include_gen){
        if(exists($hash_genomes{$el})){
            push(@position_genomes,$hash_genomes{$el});
        }        
        else{
        print "::genome $el not found. Skipping this genome\n";        
        }    
    
    }


}



if($flag_exclude eq "1"){
    foreach $el (@exclude_gen){
        delete($hash_genomes{$el});
    }

    @position_genomes=sort(values(%hash_genomes));

}

@position_genomes=sort(@position_genomes);

open(OUT,">${file_out}");

@explode=split(/\t/,$header);
print OUT join("\t",@explode[@position_genomes])."\n";


while($l=<IN>){
    @explode=split(/\t/,$l);
    print OUT join("\t",@explode[@position_genomes])."\n";


}


close(IN);
close(OUT);

