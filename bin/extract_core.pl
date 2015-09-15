#!/usr/bin/perl

#convert the file to a numerical matrix

use Getopt::Long;

#::usage $0 -d <dir_genomes> -c <n_cpu>

$file_in="";
$ref_strain="";

GetOptions ("in=s" => \$file_in,"ref=s"   => \$ref_strain) or die("::usage: $0 -in <table_genes> -ref <ref_strain>\n");



if(($file_in eq "") or ($ref_strain eq '')){
    print "::usage $0 -in <input_file> -ref <ref_strain>\n";
    exit();
}



print "::searching for core genes\n";

open(IN,"<$file_in");



open(OUT,">core_genes.txt");


$comment=<IN>;
$comment=~s/^#//;
print OUT "#gene\t".$comment;

@explode=split(/\t/,$comment);

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


$core_groups=0;
%core_genes=();

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
    
    if($count_ones eq ($#data+1)){
        print OUT join("\t",@row)."\n";
        $core_groups++;
        foreach $gene (@explode_id){
            $core_genes{$gene}=1;
        }

    }
    


}

close(IN);
close(OUT);

print "::There are ${core_groups} core gene groups\n";
@core_genes=keys(%core_genes);


print "::There are ".($#core_genes+1)." core genes\n";

