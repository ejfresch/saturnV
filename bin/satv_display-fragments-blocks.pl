#!/usr/bin/perl -I /project/rclevesq/users/lfreschi/tasks/pangenome/saturnv/bin



$dir_web_data=/project/rclevesq/users/lfreschi/tasks/pangenome/saturnv/web


use Getopt::Long;
use LibFASTA;

$file_in="";
$ref_strain="";

$data_dir="";

GetOptions ("tab=s" => \$file_in,"ann=s" => \$gff_ann, "out=s" => \$file_out,"include_mums=s" => \$include_mums,"include=s" => \$include) or die("::usage: $0 -tab <table_genes> -ann <gff_ref_strain> -out <out_dir> -include_mums <include_file> -include <include_file>\n");



if(($file_out eq "") or ($gff_ann eq "") or ($include eq "") or ($file_in eq "") or ($include_mums eq "")){
    print "::usage: $0 -tab <table_genes> -ann <gff_ref_strain> -out <out_dir> -include_mums <include_file> -include <include_file>\n";
    exit();
}



print "::reading the .gff file of the references (${gff_ann})\n";
#I load the data about the reference from the gff file
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
        print "--scaffold: ".$fragment_name."\n";        
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


=head

*the hash %fragments contains the data from the gff file.
*the hash data_strains contains the data about the fragments -- comparison to the ref strain

=cut

print "::reading the fasta binary files\n";
%data_strains=();

@genomes_mums=`cat $include_mums`;
chomp(@genomes_mums);

foreach $genome (@genomes_mums){

    $name_genome=`basename $genome`;
    chomp($name_genome);
    $name_genome=~/_vs_([A-Za-z0-9\-]+)\.cov/;
    $name_genome_ok=$1;
    $name_genome_ok.=".faa";
   
    $ref=&LibFASTA::read_FASTA("$genome"," ","1");
    $data_strains{$name_genome_ok}=$ref;
     print "--file read for $name_genome_ok\n";
}



#I check the array

#print $data_strains{"PAO1_vs_ID84-S8.final.scaffolds.cov"}{"PAO1_1"}."\n";

print "::reading the $file_in matrix and extracting data\n";

$cmd="satv_extract-subset-matrix.pl -tab $file_in -out slide_matrix.tsv -include $include";
system($cmd);


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
mkdir("${file_out}/svg");

foreach $frag (keys(%fragments)){


    $y_coord=10;
    open(OUT,">${file_out}/svg/${frag}.svg");
    
    #shape of the canvas        
    print OUT '<svg xmlns="http://www.w3.org/2000/svg" version="1.1"'." width=\"10000\" height=\"".(50+(20*($#genomes+1)))."\">\n";
    print OUT '<g id="group">'."\n";

    $svg_length_frag=$fragments{$frag}{"len"};
    #I draw the trail
    print OUT '<rect title="strain:'.$ref_strain.'" x="1" y="10" rx="1" ry="1" width="10000" height="7" fill="purple" stroke="black" stroke-width="0.1"/>'."\n";

     $ref=$fragments{$frag}{"fs"};
    %hash_frag=%$ref;

    foreach $key (keys(%hash_frag)){
        
    
        print OUT '<line id="'.$key.'" title="id_ref:'.$key."; prod:".$hash_frag{$key}{"prod"}.'; coord_start:'.$hash_frag{$key}{"start"}.'; coord_end:'.$hash_frag{$key}{"end"}.'" x1="'.($hash_frag{$key}{"start"}*10000/$svg_length_frag).'" y1="19" x2="'.($hash_frag{$key}{"end"}*10000/$svg_length_frag).'" y2="19" stroke="green" stroke-width="3"/>'."\n";   
        
        $y_coord_gen=19;

    foreach $genome (@genomes){

                        
            if($genome eq $ref_strain){next;}
            $y_coord_gen+=15;

            
            

            if(exists($db_el_ref{$key}{$genome})){

                 print OUT '<line title="id_ref:'.$key.'; prod:'.$hash_frag{$key}{"prod"}.'; coord_start:'.$hash_frag{$key}{"start"}.'; coord_end:'.$hash_frag{$key}{"end"}.'" x1="'.($hash_frag{$key}{"start"}*10000/$svg_length_frag).'" y1="'.$y_coord_gen.'" x2="'.($hash_frag{$key}{"end"}*10000/$svg_length_frag).'" y2="'.$y_coord_gen.'" stroke="green" stroke-width="3"/>'."\n";
            }else{
                 print OUT '<line title="id_ref:'.$key.'; prod:'.$hash_frag{$key}{"prod"}.'; coord_start:'.$hash_frag{$key}{"start"}.'; coord_end:'.$hash_frag{$key}{"end"}.'" x1="'.($hash_frag{$key}{"start"}*10000/$svg_length_frag).'" y1="'.$y_coord_gen.'" x2="'.($hash_frag{$key}{"end"}*10000/$svg_length_frag).'" y2="'.$y_coord_gen.'" stroke="green" stroke-width="3"/>'."\n";
            }


        }


    }



    


    foreach $genome (@genomes){
            
   

                        
            if($genome eq $ref_strain){next;}
            print "--drawing genes/fragments of ".$genome."\n";


            $y_coord+=15;
            if(exists($data_strains{$genome}{$frag})){
                #print "GENOME: $genome\n";
                #print "FRAG: $frag\n";
    

                

                #print $genome."\n";
                $string=$data_strains{$genome}{$frag};
                #print "NO. ELEM.:".$#arr."\n";
                @arr_to_print=split(//,$string);
                #print join("+",@arr_to_print)."\n";
                
                #print OUT '<rect '.'" x="'.$hash_frag{$key}{"start"}.'" y="6" rx="1" ry="1" width="'.$len_frag.'" height="7" style="fill:black;stroke:black;stroke-width:2;" />'."\n";
                $count=0;                
                $flag_match=0;
                for($e=0;$e<=$#arr_to_print;$e++){
                
                    if(($arr_to_print[$e] eq "1") and ($flag_match eq "0")){
                         $flag_match=1;
                         $record_start=$e+1;                         
                    }
                    
                    if(($arr_to_print[$e] eq "0") and ($flag_match eq "1")){
                         $flag_match=0;
                         $record_end=$e; 

                        $svg_record_start=$record_start*10000/$svg_length_frag;
                        $svg_record_end=$record_end*10000/$svg_length_frag;
                        $svg_record_len=$svg_record_end-$svg_record_start;


                        print OUT '<rect title="strain:'.$genome.'; coord_start:"'.$record_start.'; coord_end:'.$record_end.'" x="'.$svg_record_start.'" y="'.$y_coord.'" rx="1" ry="1" width="'.($svg_record_len).'" height="7" fill="purple" stroke="black" stroke-width="0.1"/>'."\n";

                        
                    }



                }


            }else{
                print "ERROR\n";
            }
            
           

        }




    print OUT '</g>'."\n";
    print OUT "</svg>\n";
    

    close(OUT);




}


#I copy the things that are needed to run the html stuff

print "::I copy the web data\n";

$cmd="cp -r ${dir_web_data}/css ${file_out}";
system($cmd);

$cmd="cp -r ${dir_web_data}/fonts ${file_out}";
system($cmd);

$cmd="cp -r ${dir_web_data}/js ${file_out}";
system($cmd);

$cmd="cp -r ${dir_web_data}/img ${file_out}";
system($cmd);


print "::I generate the html page\n";

$cmd="cat ${dir_web_data}/segment1.txt > ${file_out}/index.html";
system($cmd);


open(OUT,">>$file_out/index.html");
foreach $frag (keys(%fragments)){

print OUT '<option value="svg/'.$frag.'.svg">'.$frag.'.svg</option>'."\n";

}

$cmd="cat ${dir_web_data}/segment2.txt >> ${file_out}/index.html";
system($cmd);



