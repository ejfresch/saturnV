#!/usr/bin/perl

use Getopt::Long;


$file_gff="";
$out_dir=".";
$align=0;

GetOptions ("gffs=s" => \$files_gff,"gs=s" => \$genes) or die("::usage: $0 -gffs <gff_file1>,<gff_file2> -genes <gene1>,<gene2> -ali <[0|1]> (default=0)\n");

if(($files_gff eq "") or ($genes eq "")){
    print "::usage: $0 -gffs <gff_file1>,<gff_file2> -genes <gene1>,<gene2> -ali <[0|1]> (default=0)\n";
    exit();
}


@gffs=split(/,/,$files_gff);
@genes=split(/,/,$genes);

%db_genes=();
    
foreach $gene (@genes){



    foreach $gff (@gffs){


       $strain=`basename $gff`;
       chomp($strain); 
       $strain=~s/.gff//; 


       %data_seq=();
       $contig_name=""; 

       $switch_fasta=0;

        open(IN,"<${gff}") or die "::I cannot open the ${gff} file\n";
        $comment=<IN>;

        while($line=<IN>){
            chomp($line);

            
           
            if($line=~/gene=([0-9A-Za-z\_\-]+);/){
                  $current_gene=$1;
                  #print $1;
                  if($current_gene ne $gene){#print ":0\n";
                    next;}
                  #print ":1\n";   
                  @explode_line=split(/\t/,$line);
                  $contig=$explode_line[0];
                  $start=$explode_line[3];
                  $end=$explode_line[4];
  
            }

            if($line=~/^##FASTA/){$switch_fasta=1;next;}
            if($switch_fasta eq "1"){
            
                if($line=~/>/){
                    $contig_name=$line;
                    $contig_name=~s/>//;
                    $data_seq{$contig_name}="";
                    
                }else{

                    $data_seq{$contig_name}.=$line;
                
                }




            }





            
        }



        if(exists($data_seq{$contig})){


            @seq=split(//,$data_seq{$contig});
            $sequence=join("",@seq[($start-1)..($end-1)]);
            
            $name_db=$gene."\$".$strain;
            $db_genes{$gene}{$name_db}=$sequence;
            
        }
        else{
            print "::I could not find the gene \"${gene}\" in the genome \"${strain}\"\n";
        }



        close(IN);

    }
    
    



}



foreach $key(keys(%db_genes)){

    #print ">>".$key."\n";

    open(OUT,">${key}.fasta");
        
        $ref=$db_genes{$key};
      
        %hash=%$ref;
        foreach $el (keys(%hash)){
            print OUT ">".$el."\n";
            print OUT $hash{$el}."\n";
        }
    
    close(OUT);

    if($align eq "1"){

    `muscle -in ${key}.fasta -out ${key}_aligned.fasta`;
    }


}










