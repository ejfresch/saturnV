#!/usr/bin/perl

#convert the file to a numerical matrix

use Getopt::Long;

#::usage $0 -d <dir_genomes> -c <n_cpu>

$file_in="";
$ref_strain="";
$inf_lim="";
$sup_lim="";
$file_out="";

GetOptions ("tab=s" => \$file_in, "ref=s" => \$ref_strain,"inf=s" => \$inf_lim, "sup=s" => \$sup_lim, "out=s" => \$file_out) or die("::usage: $0 -tab <table_genes> -inf <num> -sup <num> -ref <strain> -out <file_out>\n");



if(($file_in eq "") or ($ref_strain eq "") or ($file_out eq "")){
    print "::usage: $0 -tab <table_genes> -inf <num> -sup <num> -ref <strain> -out <file_out>\n";
    exit();
}



print "::searching for shared genes\n";

open(IN,"<${file_in}");



open(OUT,">${file_out}");


$comment=<IN>;
chomp($comment);
$comment=~s/^#//;
print OUT "#".$comment."\n";

@explode=split(/\t/,$comment);

if($inf_lim eq ""){
    $inf_lim=2;
}

if($sup_lim eq ""){
    $sup_lim=$#explode;
}


print $inf_lim."-".$sup_lim."\n";


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

    
    chomp($line);
    @data=split(/\t/,$line);

    $id=$data[$index];  
    @explode_id=split(/,/,$data[$index]);
  
    
    $count_ones=0;
    
    foreach $el (@data){
        if($el ne "-"){

            $count_ones++;        
        }
        
    }
    
    if(($count_ones <= ($sup_lim)) and ($count_ones >= ($inf_lim))){
        print OUT $line."\n";
        $shared_groups++;
        foreach $gene (@explode_id){
            $shared_genes{$gene}=1;
        }

    }
    


}

close(IN);
close(OUT);

print "::There are ${shared_groups} shared gene groups overeall (see ${file_out})\n";
@shared_genes=keys(%shared_genes);


print "::There are ".($#shared_genes+1)." shared genes in strain ${ref_strain}\n";

