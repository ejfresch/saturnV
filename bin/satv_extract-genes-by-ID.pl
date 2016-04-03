#!/usr/bin/perl

use Getopt::Long;

$file_gff="";
$ids="";
$out_file="";
$append=0;

GetOptions ("gff=s" => \$file_gff,"ids=s" => \$ids,"out=s" => \$out_file,"append=s" => \$append) or die("::usage: $0 -gff <gff_file> -ids <id1>,<id2> -out <out_file> -append <[0|1]> (default=0)\n");


if(($file_gff eq "") or ($ids eq "")){
    print "::usage: $0 -gff <gff_file> -ids <id1>,<id2> -out <out_file> -append <[0|1]> (default=0)\n";
    exit();
}

$gff=$file_gff;
@ids=split(/,/,$ids);

%db_ids=();
    
foreach $id (@ids){



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

            
           
            if($line=~/ID=([0-9A-Za-z\_\-]+);/){
                  $current_id=$1;
                  #print $1;
                  if($current_id ne $id){#print ":0\n";
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
            
            $name_db=$id."\$".$strain;
            $db_ids{$id}{$name_db}=$sequence;
            
        }
        else{
            print "::I could not find the id \"${id}\" in the genome \"${strain}\"\n";
        }



        close(IN);


    
    



}



foreach $key(keys(%db_ids)){

    #print ">>".$key."\n";

    if($append eq "0"){
    open(OUT,">${key}.fasta");
    }else{
    open(OUT,">>${key}.fasta");
    }        
        $ref=$db_ids{$key};
      
        %hash=%$ref;
        foreach $el (keys(%hash)){
            print OUT ">".$el."\n";
            print OUT $hash{$el}."\n";
        }
    
    close(OUT);



}










