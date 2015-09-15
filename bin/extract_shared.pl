#!/usr/bin/perl

#convert the file to a numerical matrix

use Getopt::Long;

#::usage $0 -d <dir_genomes> -c <n_cpu>

$file_in="";
$ref_strain="";
$inf_lim="";
$sup_lim="";

GetOptions ("in=s" => \$file_in, "ref=s" => \$ref_strain) or die("::usage: $0 -in <table_genes> -inf <num> -sup <num> -ref <strain>\n");



if(($file_in eq "") and ($ref_strain eq "")){
    print "::usage: $0 -in <table_genes> -inf <num> -sup <num> -ref <strain>\n";
    exit();
}



print "::searching for shared genes\n";

open(IN,"<$file_in");



open(OUT,">shared_genes.txt");


$comment=<IN>;
$comment=~s/^#//;
print OUT "#gene\t".$comment;

@explode=split(/\t/,$comment);

if(($inf_lim eq "") or ($sup_lim eq "")){
    $inf_lim=1;
    $sup_lim=$#comment;
}

$index=-1;
for($i=0;$i<=$#explode;$i++){
    if($explode[$i] eq $ref_strain){
        $index=$i;
        break;    
    }
}

if($index eq "-1"){
    print "::Reference strain not found\n";
    exit();
}


$shared_groups=0;
%shared_genes=();

while($line=<IN>){
    @row=();
    
    chomp($line);
    @data=split(/\t/,$line);

    $id=$data[$index];  
    @explode_id=split(/,/,$data[$index]);
  
    
    push(@row,$id);


    $count_ones=0;
    
    foreach $el (@data){
        if($el ne "-"){
            push(@row,"1");
            $count_ones++;        
        }
        else{
            push(@row,"0");
        }
    }
    
    if(($count_ones < ($#data+1)) and ($count_ones > 1)){
        print OUT join("\t",@row)."\n";
        $shared_groups++;
        foreach $gene (@explode_id){
            $shared_genes{$gene}=1;
        }

    }
    


}

close(IN);
close(OUT);

print "::There are ${shared_groups} core gene groups\n";
@core_genes=keys(%shared_genes);


print "::There are ".($#core_genes+1)." core genes\n";

