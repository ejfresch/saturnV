#!/usr/bin/perl

$current_dir=`pwd`;
chomp($current_dir);
$web_dir=$current_dir;
$web_dir=~s/bin/web/;

@files_to_modify=("satv_pangenome4-fast.pl","satv_pangenome4.pl","satv_pangenome-mosaic.pl","satv_display-fragments-blocks.pl","satv_display-fragments.pl","satv_display-blocks.pl", "satv_cog2.pl");

foreach $file (@files_to_modify){
    open(IN,"<$file");  
    @data=<IN>;
    chomp(@data);
    close(IN);


    for($i=0;$i<$#data;$i++){
        
        $line=$data[$i];

        if($line=~/^\$base_path=/){
             $data[$i]='$base_path="'.$current_dir.'/database";';
        }  
  
	if($line=~/#\!\/usr\/bin\/perl -I/){
            $data[$i]="#!/usr/bin/perl -I $current_dir";
        }    

        if($line=~/^\$dir_web_data=/){
            
            $data[$i]='$dir_web_data="'.$web_dir.'";';
        }    

        
    
    }


    open(OUT,">$file");  
    print OUT join("\n",@data);   
    close(OUT);


}






