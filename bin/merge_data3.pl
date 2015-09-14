#!/usr/bin/perl


use Graph;
use Getopt::Long;


GetOptions ("i=s" => \$file_in) or die("::usage: $0 -i <input_file> \n");



$graph = Graph::Undirected->new;


print "::constructing the graph\n";
$lines_read=0;

open(IN,"<$file_in")or die "::[ERROR] I cannot open the file $file_in\n";
$comment=<IN>;
chomp($comment);
$header=$comment;
while($line=<IN>){
    chomp($line);
    if($line eq ""){next;}
    @data=split(/\t/,$line); 
    %hash_line=();
    foreach $el (@data){
        if($el ne "-"){        
            $hash_line{$el}=1;
        }    
    }
    $graph->add_path(keys(%hash_line));
    $lines_read++;    
    print "-- $lines_read lines read\n";
    



}



close(IN);

print "\n";
#print "The graph is $graph\n";
print "::perl is determining the connected components of the graph\n";
@cc=$graph->connected_components();

$comment=~s/#//;
$comment=~s/.faa//g;
@genomes=split(/\t/,$comment);



open(OUT,">table_linked3.tsv");
print OUT $header."\n";

$count=0;
foreach $c (@cc){
    $count++;
    print "::writing down the connected components of the graph (${count}/$#cc)\n";


    %hash_results=();
    @current_cc=@$c;
    foreach $element (@current_cc){

        if($element eq "-"){next;}
        #print $element."\n";
        @explode=split(/_/,$element);
        #print $explode[0]."\n";
        $hash_results{$explode[0]}{$element}=1;        


     
    }

    @final_arr=();    
    
    foreach $genome (@genomes){
    #print $genome."\n";
    if(exists($hash_results{$genome})){
       
        $ref=$hash_results{$genome};

        %arr=%$ref;
        push(@final_arr,join(",",sort(keys(%arr))));
    }
    else{push(@final_arr,"-");}

    }

    print OUT join("\t",@final_arr)."\n";

}

close(OUT);








