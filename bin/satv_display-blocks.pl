#!/usr/bin/perl

use Getopt::Long;


$file_in="";
$ref_strain="";

$data_dir="";

GetOptions ("tab=s" => \$file_in,"ann=s" => \$gff_ann, "out=s" => \$file_out,"include=s" => \$include) or die("::usage: $0 -tab <table_genes> -ann <gff_ref_strain> -out <out_dir> -include <include_file>\n");



if(($file_in eq "") or ($file_out eq "") or ($gff_ann eq "") or ($include eq "")){
    print "::usage: $0 -tab <table_genes> -ann <gff_ref_strain> -out <out_dir> -include <include_file>\n";
    exit();
}



$ref_strain=`basename $gff_ann`;
chomp($ref_strain);

$ref_strain=~s/.gff/.faa/;


open(IN,"<${gff_ann}") or die "::I cannot open the ${gff_ann} file\n";
$comment=<IN>;


%fragments=();

while($line=<IN>){
    chomp($line);
    if($line=~/##sequence-region/){

        
        @data=split(/ /,$line);
        $fragment_name=$data[1];
        print $fragment_name."\n";        
        $last_pos=$data[3];
        $fragments{$fragment_name}{"len"}=$last_pos;
    }
    elsif($line=~/##FASTA/){last;}
    else{
        @data=split(/\t/,$line);
        $fragment_name=$data[0];
        #print "FRAG_NAME=".$fragment_name."\n";
        $start=$data[3];
        $end=$data[4];

        $details=$data[8];
        @details_exploded=split(/;/,$details);
        

        for($j=0;$j<=$#details_exploded;$j++){
        
            if($details_exploded[$j]=~/^ID=/){
                $id=$details_exploded[$j];
                $id=~s/ID=//g;
            }
            

            if($details_exploded[$j]=~/^product=/){

                $product=$details_exploded[$j];            
                $product=~s/product=//g;
            }
        
        }

        $fragments{$fragment_name}{"fs"}{$id}{"start"}=$start;       
        $fragments{$fragment_name}{"fs"}{$id}{"end"}=$end;       
        $fragments{$fragment_name}{"fs"}{$id}{"prod"}=$product;       
    }

}



close(IN);



#I generate the file file with all the features for the other strains


$cmd="satv_extract-subset-matrix.pl -tab $file_in -out slide_matrix.tsv -include $include";
system($cmd);



#I read the file I generated


%db_el_ref=();
%db_el_not_in_ref=();


open(IN,"<slide_matrix.tsv");
$comment=<IN>;
chomp($comment);
$comment=~s/^#//;
@genomes=split(/\t/,$comment);


$index=-1;
for($i=0;$i<=$#genomes;$i++){
    if($genomes[$i] eq $ref_strain){
        $index=$i;
        break;    
    }
}

if($index eq "-1"){
    print "::please insert the strain --${ref_strain}-- in the $include file. Otherwise I do not have the data to reconstruct the situation\n";
    exit();
}


while($line=<IN>){

    chomp($line);
    @explode_line=split(/\t/,$line);

    $test=$explode_line[$index];
    if($test eq "-"){
        for($k=0;$k<=$#explode_line;$k++){
            if($explode_line[$k] eq "-"){next;}
            else{
                $strain=$genomes[$k];
                $gen_el=$explode_line[$k];
                $db_el_not_in_ref{$strain}{$gen_el}=1;       
            }    

        }        
        

    }
    else{

        $current_el=$explode_line[$index];
        @arr_el=split(/,/,$current_el);

        for($k=0;$k<=$#explode_line;$k++){
            
            $curr_strain=$genomes[$k];
            $gen_el=$explode_line[$k];
            if($gen_el ne "-"){
            
                foreach $ref_el (@arr_el){
                    $db_el_ref{$ref_el}{$curr_strain}=$gen_el;    
                }
            }                                    
        
    
        }
    }
   

}


close(IN);



#I write the svg files

#I get the fragments

mkdir($file_out);

foreach $frag (keys(%fragments)){
    open(OUT,">${file_out}/${frag}.svg");
    
    #shape of the canvas        
    print OUT "<svg width=\"10000\" height=\"100\">\n";
    
    #I draw the trail
    #print OUT '<line x1="1" y1="14" x2="10000" y2="14" stroke="black" stroke-width="5"/>'."\n";   
    

    $ref=$fragments{$frag}{"fs"};
    $svg_length_frag=$fragments{$frag}{"len"};
    %hash_frag=%$ref;

    foreach $key (keys(%hash_frag)){

        $svg_start=$hash_frag{$key}{"start"}*10000/$svg_length_frag;
        $svg_end=$hash_frag{$key}{"end"}*10000/$svg_length_frag;
        $svg_len=$svg_end-$svg_start;    



        print OUT '<rect strain="'.$ref_strain.'" id="'.$key.'" prod="'.$hash_frag{$key}{"prod"}.'" x="'.$svg_start.'" y="10" rx="1" ry="1" width="'.$svg_len.'" height="7" fill="green" stroke="black" stroke-width="0.1"/>'."\n";

        $y_coord=10;

        foreach $genome (@genomes){

                        
            if($genome eq $ref_strain){next;}
            $y_coord+=10;

            if(exists($db_el_ref{$key}{$genome})){

                 print OUT '<rect strain="'.$genome.'" id="'.${genome}."-".${key}.'" prod="'.$hash_frag{$key}{"prod"}.'" x="'.$svg_start.'" y="'.$y_coord.'" rx="1" ry="1" width="'.$svg_len.'" height="7" fill="green" stroke="black" stroke-width="0.1"/>'."\n";
            }else{
                print OUT '<rect strain="'.$genome.'" id="'.${genome}."-".${key}.'" prod="'.$hash_frag{$key}{"prod"}.'" x="'.$svg_start.'" y="'.$y_coord.'" rx="1" ry="1" width="'.$svg_len.'" height="7" fill="red" stroke="black" stroke-width="0.1"/>'."\n";
            }


        }



    }

    




    print OUT "</svg>\n";
    

    close(OUT);
=head
    #I generate the html file
    open(OUT,">${file_out}/index.html");

    print OUT "<!DOCTYPE html>\n<html>\n";

    print OUT <<END;
<script src="http://code.jquery.com/jquery-1.9.1.js"></script> 
<script> 
$(function() {
$("#includedContent").load(""); 
}); 
</script>
END

    print OUT "<body>\n"; 
    
       
    print OUT "<body>\n</html>\n";
    

    close(OUT);

=cut



}








