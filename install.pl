#!/usr/bin/perl

#::usage $0 <install_dir>
use Getopt::Long;

$install_dir="";

GetOptions ("d=s" => \$install_dir) or die("::usage: $0 -d <install_dir>\n");

if($install_dir eq ""){
    print "::usage: $0 -d <install_dir>\n";
    exit();
}


chdir($install_dir);


#I install prokka
print "::I download prokka\n";
mkdir("prokka");
chdir("prokka");
`wget http://www.vicbioinformatics.com/prokka-1.11.tar.gz`;

`tar -xvzf prokka-1.11.tar.gz`;






#I add the path on 
`echo PATH=${install_dir}:\$PATH >> ~/.bashrc`;
















