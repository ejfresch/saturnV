#!/usr/bin/perl

use Getopt::Long;


$res_mum="";
$ref_file="";
$dir_out="";

GetOptions ("mum=s" => \$res_mum,"ref=s" => \$ref_file, "out=s" => \$dir_out) or die("::usage: $0 -mum <results_mummer> -ref <ref_strain> -out <out_dir>\n");


if(($res_mum eq "") or ($ref_file eq "")){
    print "::usage: $0 -mum <results_mummer> -ref <ref_strain> -out <out_dir>\n";
    exit();
}


#I check if there reference is there and I build an hash with all data


if((-e "${dir_genomes}/$ref_file")){
    $ref_file="${dir_genomes}/$ref_file";
   
}

if(!(-e $ref_file)){

   print "::Sorry! I cannot find $ref\n";
   exit();
}


#here is the routine to read a fasta file and return the length of the sequences contained
sub read_FASTA{

	my $file_name=$_[0];
	my $sep=$_[1];
	
	#pos: 1,2,3,..;
	my $pos=$_[2];


	my %sequences=();
	my $id="";
	
	open(IN,"<$file_name")|| die "I cannot open the file $file_name\n";


		while(my $line=<IN>){
		
			chomp($line);
		
			if($line=~/^>/){
			
				if(!($sep eq "")){
					#print "here $sep\n";
					my @arr=split(/[$sep]+/,$line);
				
					#print join("#",@arr)."\n";
					$id=$arr[$pos-1];
					
					}
				else{
					$id=$line;
					
				}	
				
				#print $id."\n";
					
					
				$id=~s/>//;	
				
					
				$sequences{$id}=0;
		
			}
			else{
			
				$sequences{$id}+=length($line);
				
			}
		
		
	
		}

	close(IN);

	my $ref_sequences=\%sequences;
	return $ref_sequences;

}




$ref=&read_FASTA($ref_file," ",1);


if($dir_out ne ""){
    mkdir($dir_out);
}

@genomes=split(/,/,$res_mum);

foreach $genome (@genomes){

    %hash=%$ref;


    #I build the hash of arrays with the sequences

    foreach $key (keys(%hash)){
        #print $key."\n";
    
        $l=$hash{$key};
        @arr = ("0") x $l;
        
        $hash{$key}=\@arr;



    }

    @all_keys=keys(%hash);    

    open(IN,"<$genome");
    while($line=<IN>){
        chomp($line);
        if($line=~/^>/){next;}
        else{
            @current_arr=split(/\s+/,$line);
            #print join("\$",@current_arr)."\n";

            if($#all_keys eq "0"){
                $scaffold=$all_keys[0];
                $pos=$current_arr[1];
                #print $pos."\n";            
                $l_match=$current_arr[3];

            }else{
                $scaffold=$current_arr[1];

                $pos=$current_arr[2];
                #print $pos."\n";            
                $l_match=$current_arr[4];

            }
        
                
            $start=$pos-1;
            $end=($start+$l_match);

            #print $scaffold."&".start."&".$end."\n";   
            if(exists($hash{$scaffold})){
                for($p=$start;$p<=$end;$p++){
                    $hash{$scaffold}[$p]=1;
                }                
                                
            }else{
                print "::$scaffold is not a contig/scaffold present in the file you provided me ($ref_file)\n";
            }            


        }


    }

    close(IN);


    

    if($dir_out ne ""){

        $file_out=$genome;
        $file_out=~s/\.mum/\.cov/;
        $file_out=~s/.+\///;
        
        open(OUT,">${dir_out}/${file_out}");
    }

    $overall_length=0;
    $coverage=0;

    foreach $el (keys(%hash)){
        #print $el."\n";
        $r_arr=$hash{$el};
        @arr=@$r_arr;
    
        if($dir_out ne ""){
            print OUT ">".$el."\n".join("",@arr)."\n";
        }
        
        foreach $elem (@arr){
            if($elem eq "1"){
                $overall_length++;
                $coverage++;   
            }
            elsif($elem eq "0"){$overall_length++;}
            #else{print "::Strange values in the array I use to calculate the coverage of the reference sequence --$elem--.\n";}

        }
           
    
    }

    print "::overall_coverage (${genome}): $coverage / $overall_length (".($coverage/$overall_length*100)."%)\n";

    if($dir_out ne ""){

    close(OUT);

    }

}





