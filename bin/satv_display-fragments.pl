#!/usr/bin/perl -I /home/lfreschi/tasks/pangenome/saturnv/bin


use Getopt::Long;
use LibFASTA;

$dir_web_data="/home/lfreschi/tasks/pangenome/saturnv/web";

$file_in="";
$ref_strain="";

$data_dir="";

GetOptions ("ann=s" => \$gff_ann, "out=s" => \$file_out,"include_mums=s" => \$include) or die("::usage: $0 -ann <gff_ref_strain> -out <out_dir> -include_mums <include_file>\n");



if(($file_out eq "") or ($gff_ann eq "") or ($include eq "")){
    print "::usage: $0 -ann <gff_ref_strain> -out <out_dir> -include_mums <include_file>\n";
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

%data_strains=();

@genomes=`cat $include`;

chomp(@genomes);

foreach $genome (@genomes){
    $ref=&LibFASTA::read_FASTA("$genome"," ","1");
    $data_strains{$genome}=$ref;
}



#I check the array

#print $data_strains{"PAO1_vs_ID84-S8.final.scaffolds.cov"}{"PAO1_1"}."\n";




#I write the svg files

#I get the fragments

mkdir($file_out);

foreach $frag (keys(%fragments)){


    $y_coord=10;
    open(OUT,">${file_out}/${frag}.svg");
    
    #shape of the canvas        
    print OUT "<svg width=\"10000\" height=\"".(50+(20*($#genomes+1)))."\">\n";
    $svg_length_frag=$fragments{$frag}{"len"};
    #I draw the trail
    print OUT '<rect strain="'.$ref_strain.'" x="1" y="10" rx="1" ry="1" width="10000" height="7" fill="purple" stroke="black" stroke-width="0.1"/>'."\n";

     $ref=$fragments{$frag}{"fs"};
    %hash_frag=%$ref;

    foreach $key (keys(%hash_frag)){
        
    
        print OUT '<line id="'.$key.'" prod="'.$hash_frag{$key}{"prod"}.'" coord_start="'.$hash_frag{$key}{"start"}.'" coord_end="'.$hash_frag{$key}{"end"}.'" x1="'.($hash_frag{$key}{"start"}*10000/$svg_length_frag).'" y1="19" x2="'.($hash_frag{$key}{"end"}*10000/$svg_length_frag).'" y2="19" stroke="green" stroke-width="3"/>'."\n";   
    }
     $y_coord=15;

    foreach $genome (@genomes){

            $strain_name=$genome;
            $strain_name=~/_vs_([A-Za-z0-9\-\.]+)\.cov/;
            $sname=$1;
            #print "--".$sname."--\n";
                        
            if($genome eq $ref_strain){next;}
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


                        print OUT '<rect strain="'.$sname.'" coord_start="'.$record_start.'" coord_end="'.$record_end.'" x="'.$svg_record_start.'" y="'.$y_coord.'" rx="1" ry="1" width="'.($svg_record_len).'" height="7" fill="purple" stroke="black" stroke-width="0.1"/>'."\n";

                        
                    }



                }


            }else{
                print "ERROR\n";
            }
            

            
        

        }





    print OUT "</svg>\n";
    

    close(OUT);

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




}









