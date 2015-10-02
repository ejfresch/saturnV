#!/usr/bin/perl

#::usage $0 <install_dir>
use Getopt::Long;

$install_dir="";

GetOptions ("d=s" => \$install_dir) or die("::usage: $0 -d <install_dir>\n--example: $0 -d /home/luca/saturnv\n\nPlease remember to use an absolute path!\n");

if($install_dir eq ""){
    print "::usage: $0 -d <install_dir>\n";
    print "--example: $0 -d /home/luca/saturnv\n\nPlease remember to use an absolute path!\n";

    exit();
}



#testing all dependencies

print "::I check if the dependences are installed\n";

#R
$ret=system("which Rscript > /dev/null");
if($ret ne 0){
	print "--please install R\n";
    	print "hint: \"sudo apt-get install r-base\" or ask your system administrator\n";	
	exit();
}
else{
	print "--R is installed\n";
}



if(! (-e $install_dir)){
	print "::I create the folder $install_dir\n";

	mkdir($install_dir);
}


print "::I check if usearch is installed\n";

$ret=system("which usearch8 > /dev/null");

if(($ret ne 0) and (!(-e "${install_dir}/usearch/usearch8"))){
	print "--please install usearch\n";
	mkdir("${install_dir}/usearch");

    	print "hint: open your browser, go to http://www.drive5.com/usearch/download.html and follow the instructions.\n*You should get a file named usearch8.1.1756_i86linux32 (the version may vary)\n*Copy it inside the directory $install_dir/usearch/\n*Rename the file usearch8 (\"mv <yourfile> usearch8\")\n*Rerun install2.pl\n";	


	exit();
}
else{
	print "--usearch is installed\n";
}



print "::Copying the bin/ folder to $install_dir\n";
$cmd="cp -r bin ${install_dir}";
system($cmd);

print "::Copying the licence and README files\n";
$cmd="cp -r COPYING ${install_dir}";
system($cmd);
$cmd="cp -r README ${install_dir}";
system($cmd);



print "::Changing directory -- $install_dir\n";
chdir($install_dir);


print "::I check if the remaining dependences are installed\n";

$ret=system("which prokka > /dev/null");
if($ret ne 0){

	print "--I download and install prokka\n";

	mkdir("prokka");
	chdir("prokka");
	`wget http://www.vicbioinformatics.com/prokka-1.11.tar.gz`;
	`tar -xvzf prokka-1.11.tar.gz`;

	chdir("../");


}
else{
	print "--prokka is installed\n";
}





$ret=system("which blastp > /dev/null");
if($ret ne 0){

	#I install blast
	print "--I download and install blast\n";
	mkdir("blast");
	chdir("blast");
	`wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.2.31+-x64-linux.tar.gz`;
	`tar -xvzf ncbi-blast-2.2.31+-x64-linux.tar.gz`;

	chdir("../");


}
else{
	print "--prokka is installed\n";
}



$ret=system("which mummer > /dev/null");
if($ret ne 0){

	$ret_make=system("which make > /dev/null");
	if($ret_make ne 0){
		print "--please install make. I need it to compile mummer\n";
    		print "hint: \"sudo apt-get install make\" or ask your system administrator\n";	
		exit();
	}
	else{
		print "--make is installed\n";
	}

	$ret_csh=system("which csh > /dev/null");
	if($ret_csh ne 0){
		print "--please install csh. I need it to compile mummer\n";
    		print "hint: \"sudo apt-get install csh\" or ask your system administrator\n";	
		exit();
	}
	else{
		print "--csh is installed\n";
	}



	$ret_gcc=system("which g++ > /dev/null");
	if($ret_gcc ne 0){
		print "--please install g++. I need it to compile mummer\n";
    		print "hint: \"sudo apt-get install g++\" or ask your system administrator\n";	
		exit();
	}
	else{
		print "--g++ is installed\n";
	}



	#I install mummer
	print "::I download and install mummer\n";
	mkdir("mummer");
	chdir("mummer");

	`wget http://downloads.sourceforge.net/project/mummer/mummer/3.23/MUMmer3.23.tar.gz`;
	`tar -xvzf MUMmer3.23.tar.gz`;
	chdir("MUMmer3.23");
	`make install`;


	chdir("../../");

}
else{
	print "--prokka is installed\n";
}



#setting up the paths for the modules

print "::I set up the path for LibFASTA\n";
chdir("bin/");
`set_dir_modules.pl`;


print "::I set up the path for the SaturnV binaries\n";

#setting up the paths for saturnV on .bashrc
$cmd="printf \"#SaturnV -- paths to binaries\n\" >> ~/.bashrc";
system($cmd);


$cmd="printf \"export PATH=\$PATH:${install_dir}/bin:${install_dir}/blast/ncbi-blast-2.2.31+/bin:${install_dir}/prokka/prokka-1.11/binaries/linux::${install_dir}/mummer/MUMmer3.23/:${install_dir}/usearch\n\" >> ~/.bashrc";
system($cmd);




print "::Everything is ok and SaturnV is ready on the launchpad\n";












