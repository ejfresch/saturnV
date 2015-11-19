#!/usr/bin/perl

#::usage $0 <install_dir>
use Getopt::Long;

$install_dir="";
$prokka_new_inst=0;


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
	print "Then rerun install2.pl with the same options you provided above ($0 -d ${install_dir})\n";			
	exit();
}
else{
	print "--R is installed\n";
    system("./bin/install-r-packages.R");
}






#perl -- graph
$ret_graph=system("perl -e 'use Graph;' 2>install_err.log > /dev/null");
if($ret_graph ne 0){
	print "--please install the Graph Perl module\n";
    	print "hint: \"sudo apt-get install libgraph-perl\" or ask your system administrator\n\n";
	print "Then rerun install2.pl with the same options you provided above ($0 -d ${install_dir})\n";		
	exit();
}
else{
	print "--the Graph Perl module is installed\n";
}


#perl -- fork manager
$ret_fm=system("perl -e 'use Parallel::ForkManager;' 2>install_err.log > /dev/null");
if($ret_fm ne 0){
	print "--please install the ForkManager Perl module\n";
    	print "hint: \"sudo apt-get install libparallel-forkmanager-perl\" or ask your system administrator\n\n";	
	print "Then rerun install2.pl with the same options you provided above ($0 -d ${install_dir})\n";		
	exit();
}
else{
	print "--the ForkManager Perl module is installed\n";
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

    	print "hint: open your browser, go to http://www.drive5.com/usearch/download.html and follow the instructions.\n*You should get a file named usearch8.1.1756_i86linux32 (NOTE: the version of usearch -- 8.1.1756 in this case -- may vary)\n*Copy it inside the directory $install_dir/usearch/\n*Rename the file usearch8 (\"mv ${install_dir}/usearch/usearch8.1.1756_i86linux32  usearch8\")\n*Be sure that the permissions are correctly set (\"chmod 775 ${install_dir}/usearch/usearch8\")\n*Rerun install2.pl with the same options you provided above ($0 -d ${install_dir})\n";	


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

	#perl -- XML::Simple
	$ret_xml=system("perl -e 'use XML::Simple;' 2>install_err.log > /dev/null");
	if($ret_xml ne 0){
		print "--please install the XML::Simple Perl module\n";
    		print "hint: \"sudo apt-get install libxml-simple-perl\" or ask your system administrator\n\n";
		print "Then rerun install2.pl with the same options you provided above ($0 -d ${install_dir})\n";		
		exit();
	}
	else{
		print "--the XML::Simple Perl module is installed\n";
	}


	#perl -- Bio
	$ret_bp=system("perl -e 'use Bio::Root::Version;' 2>install_err.log > /dev/null");
	if($ret_bp ne 0){
		print "--please install the BioPerl Perl module\n";
    		print "hint: \"sudo apt-get install bioperl\" or ask your system administrator\n";
    		print "NOTE: bioperl consist of a lot of packages. They are not all needed by SaturnV (we only need the modules required to run Prokka). However, we advice you to install all packages since they may be useful for you in the future.\n\n";

		print "Then rerun install2.pl with the same options you provided above ($0 -d ${install_dir})\n";		
		exit();
	}
	else{
		print "--the BioPerl Perl module is installed\n";
	}




	print "--I download and install prokka\n";
	$prokka_new_inst=1;
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
	print "--blastp is installed\n";
}



$ret=system("which mummer > /dev/null");
if($ret ne 0){

	$ret_make=system("which make > /dev/null");
	if($ret_make ne 0){
		print "--please install make. I need it to compile mummer\n";
    		print "hint: \"sudo apt-get install make\" or ask your system administrator\n";	
		print "Then rerun install2.pl with the same options you provided above ($0 -d ${install_dir})\n";		

		exit();
	}
	else{
		print "--make is installed\n";
	}

	$ret_csh=system("which csh > /dev/null");
	if($ret_csh ne 0){
		print "--please install csh. I need it to compile mummer\n";
    		print "hint: \"sudo apt-get install csh\" or ask your system administrator\n";	
		print "Then rerun install2.pl with the same options you provided above ($0 -d ${install_dir})\n";		

		exit();
	}
	else{
		print "--csh is installed\n";
	}



	$ret_gcc=system("which g++ > /dev/null");
	if($ret_gcc ne 0){
		print "--please install g++. I need it to compile mummer\n";
    		print "hint: \"sudo apt-get install g++\" or ask your system administrator\n";	
		print "Then rerun install2.pl with the same options you provided above ($0 -d ${install_dir})\n";		

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
	print "--mummer is installed\n";
}




$ret=system("which CONTIGuator.py 2>install_err.log > /dev/null");
if($ret ne 0){

	#I install blast
	print "--I download and install Contiguator\n";
	mkdir("contiguator");
	chdir("contiguator");
	`wget http://downloads.sourceforge.net/project/contiguator/CONTIGuator_v2.7.tar.gz`;
	`tar -xvzf CONTIGuator_v2.7.tar.gz`;
	`chmod 775 CONTIGuator_v2.7/CONTIGuator.py`;

	chdir("../");


}
else{
	print "--Contiguator is installed\n";

}








#setting up the paths for the modules

print "::I set up the path for LibFASTA\n";
chdir("bin/");
`set_dir_modules.pl`;


print "::I set up the path for the SaturnV binaries\n";

#setting up the paths for saturnV on .bashrc
$cmd="printf \"\n#SaturnV -- paths to binaries\n\" >> ~/.bashrc";
system($cmd);


$cmd="printf \"export PATH=\$PATH:${install_dir}/bin:${install_dir}/blast/ncbi-blast-2.2.31+/bin:${install_dir}/prokka/prokka-1.11/binaries/linux:${install_dir}/prokka/prokka-1.11/bin:${install_dir}/mummer/MUMmer3.23/:${install_dir}/usearch:${install_dir}/contiguator/CONTIGuator_v2.7/\n\n\" >> ~/.bashrc";
system($cmd);


print "::No installation errors reported. I remove install_err.log from ${install_dir}\n";
$cmd="rm -rf ${install_dir}/install_err.log";
system($cmd);


print "::Everything is ok and SaturnV almost ready on the launchpad\n";
print "*Please reload the .bashrc settings (\". ~/.bashrc\")\n";

if($prokka_new_inst eq "1"){
	print "*Then set up the prokka database (prokka --setupdb)\n";
}













