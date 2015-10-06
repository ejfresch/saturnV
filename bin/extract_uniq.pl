#!/usr/bin/perl

#convert the file to a numerical matrix

use Getopt::Long;

#::usage $0 -d <dir_genomes> -c <n_cpu>

$file_in="";
$ref_strain="";

GetOptions ("tab=s" => \$file_in,"ref=s" => \$ref_strain, "out=s" => \$file_out) or die("::usage: $0 -tab <table_genes> -out <file_out> -ref <ref_strain> [optional]\n");



if(($file_in eq "") or ($file_out eq "")){
    print "::usage: $0 -tab <table_genes> -out <file_out> -ref <ref_strain> [optional]\n";
    exit();
}



print "::searching for uniq genes\n";

open(IN,"<$file_in");
$comment=<IN>;
chomp($comment);
$comment=~s/^#//;

@explode=split(/\t/,$comment);



@explode=split(/\t/,$comment);

$index=-1;

if($ref_strain ne ""){

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

}

open(OUT,">${file_out}");
print OUT "#".$comment."\n";

%uniq_genes=();
$uniq_groups=0;

while($line=<IN>){

    $current_pos=-2;
    
    chomp($line);
    @data=split(/\t/,$line);

    $count_ones=0;
    
    for($i=0;$i<=$#data;$i++){
        if($data[$i] eq "-"){
        
         }
        else{
            $current_pos=$i;

            $count_ones++;
          
            @explode_id=split(/,/,$data[$i]);
  
            
        }
    }
    
    if($count_ones eq 1){

        
        if($ref_strain eq ""){
           print OUT $line."\n";
        
            foreach $gene (@explode_id){
                $uniq_genes{$gene}=1;
            }    
        }
        else{
            if($current_pos eq $index){
                           print OUT $line."\n";
        
                foreach $gene (@explode_id){
                    $uniq_genes{$gene}=1;
                }    


            }

        }

    }
    


}

close(IN);
close(OUT);


@uniq_genes=keys(%uniq_genes);


if($ref_strain eq ""){

    print "::There are ".($#uniq_genes+1)." unique genes overall (see ${file_out})\n";
}
else{
    print "::There are ".($#uniq_genes+1)." unique genes in strain $ref_strain (see ${file_out})\n";


}
