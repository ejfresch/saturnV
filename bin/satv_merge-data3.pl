#!/usr/bin/perl

use Graph;
use Getopt::Long;

$file_in="";
$file_out="";


GetOptions ("in=s" => \$file_in,"out=s" => \$file_out) or die("::usage: $0 -in <input_file> -out <output_file>\n");


if(($file_in eq "") or ($file_out eq "")){

    print "::usage: $0 -in <input_file> -out <output_file>\n";
    exit();
}


$graph = Graph::Undirected->new;


print "::Constructing the graph\n";


open(IN,"<${file_in}")or die "::[ERROR] I cannot open the file $file_in\n";
$comment=<IN>;
chomp($comment);
$header=$comment;
while($line=<IN>){
    chomp($line);
    if($line eq ""){next;}
    @data=split(/\t/,$line); 
    %hash_line=();
    foreach $el (@data){
        @el_fragmented=split(/,/,$el);
        foreach $fragment(@el_fragmented){
            if($fragment ne "-"){        
            $hash_line{$fragment}=1;
            } 
        }   
    }
    @res_path=keys(%hash_line);

    if($#res_path>0){
        $graph->add_path(keys(%hash_line));
    }
    if($#res_path eq 0){
        $graph->add_vertex($res_path[0]);
    }

    



}



close(IN);


#print "The graph is $graph\n";
my $date=`date "+%Y-%m-%d %H:%M:%S"`;
chomp($date);
print "::perl is determining the connected components of the graph -- $date\n";
@cc=$graph->connected_components();

$comment=~s/#//;
$comment=~s/.faa//g;
@genomes=split(/\t/,$comment);


my $date=`date "+%Y-%m-%d %H:%M:%S"`;
chomp($date);
print "::writing down the connected components of the graph (".($#cc+1).") -- ${date}\n";

open(OUT,">${file_out}");
print OUT $header."\n";


foreach $c (@cc){



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








