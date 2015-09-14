package LibFASTA;
use strict;
use warnings;


use Exporter qw(import);
 
our @EXPORT_OK = qw(read_FASTA);


#usage: &LibFASTA::read_FASTA("Homo_sapiens.GRCh38.pep.all.fa"," ","1");
# <file>, <separator>, <position>

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
				
					
				$sequences{$id}="";
		
			}
			else{
			
				$sequences{$id}.=$line;
				
			}
		
		
	
		}

	close(IN);

	my $ref_sequences=\%sequences;
	return $ref_sequences;

}