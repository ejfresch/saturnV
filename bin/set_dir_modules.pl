#!/usr/bin/perl

$current_dir=`pwd`;
$web_dir=$current_dir;
$web_dir=~s/bin/web/;

@files_to_modify=("satv_pangenome4-fast.pl","satv_pangenome4.pl","satv_pangenome-mosaic.pl","satv_display-fragments-blocks.pl","satv_display-fragments.pl","satv_display-blocks.pl");

foreach $file (@files_to_modify){
    open(IN,"<$file");  
    @data=<IN>;
    chomp(@data);
    close(IN);


    for($i=0;$i<$#data;$i++){
        
        $line=$data[$i];

        if($line=~/#\!\/usr\/bin\/perl -I/){
            $data[$i]="#!/usr/bin/perl -I $current_dir";
        }    

        if($line=~/^\$dir_web_data=/){
            
            $data[$i]='$dir_web_data='.$web_dir;
        }    

        
    
    }


    open(OUT,">$file");  
    print OUT join("\n",@data);   
    close(OUT);


}






