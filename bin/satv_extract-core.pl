#!/usr/bin/perl


use Getopt::Long;


$file_in="";
$ref_strain="";
$file_out="";

GetOptions ("tab=s" => \$file_in, "out=s" => \$file_out,"ref=s"   => \$ref_strain) or die("::usage: $0 -tab <table_genes> -ref <ref_strain> -out <out_file>\n");



if(($file_in eq "") or ($ref_strain eq '') or ($file_out eq '')){
    print "::usage $0 -tab <input_file> -ref <ref_strain> -out <out_file>\n";
    exit();
}



print "::searching for core genes\n";

open(IN,"<$file_in");



open(OUT,">${file_out}");


$comment=<IN>;
chomp($comment);
$comment=~s/^#//;
print OUT "#".$comment."\n";

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
    
    if($count_ones eq ($#data+1)){
        print OUT $line."\n";
        $core_groups++;
        foreach $gene (@explode_id){
            $core_genes{$gene}=1;
        }

    }
    


}

close(IN);
close(OUT);

print "::There are ${core_groups} core gene groups overall (see ${file_out})\n";
@core_genes=keys(%core_genes);


print "::There are ".($#core_genes+1)." core genes for strain $ref_strain\n";

