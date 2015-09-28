#!/usr/bin/perl

$current_dir=`pwd`;

@files_to_modify=("pangenome4_fast.pl","pangenome4.pl","pangenome_mosaic.pl");

foreach $file (@files_to_modify){
    open(IN,"<$file");  
    @data=<IN>;
    chomp(@data);
    close(IN);


    for($i=0;$i<$#data;$i++){
        
        $line=$data[$i];

        if($line=~/#\!\/usr\/bin\/perl/){
            $data[$i]="#!/usr/bin/perl -I $current_dir";
        }    
    
    }


    open(OUT,">$file");  
    print OUT join("\n",@data);   
    close(OUT);


}






