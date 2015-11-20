#!/usr/bin/perl


#testing all dependencies

print "::I check if the dependences are installed\n";

#R
$ret=system("which Rscript > /dev/null");
if($ret ne 0){
	print "--please install R\n";
    	
}
else{
	print "--R (OK)\n";
    system("./bin/test-r-packages.R");


}







#perl -- graph
$ret_graph=system("perl -e 'use Graph;' 2>install_err.log > /dev/null");
if($ret_graph ne 0){
	print "--please install the Graph Perl module\n";
   
}
else{
	print "--Graph (OK)\n";
}


#perl -- fork manager
$ret_fm=system("perl -e 'use Parallel::ForkManager;' 2>install_err.log > /dev/null");
if($ret_fm ne 0){
	print "--please install the ForkManager Perl module\n";
    	
}
else{
	print "--ForkManager (OK)\n";
}


$ret=system("which usearch8 > /dev/null");

if(($ret ne 0)){
	print "--please install usearch\n";

}
else{
	print "--usearch (OK)\n";
}



$ret=system("which prokka > /dev/null");
if($ret ne 0){

	#perl -- XML::Simple
	$ret_xml=system("perl -e 'use XML::Simple;' 2>install_err.log > /dev/null");
	if($ret_xml ne 0){
		print "--please install the XML::Simple Perl module\n";
    	
	}
	else{
		print "--XML::Simple (OK)\n";
	}


	#perl -- Bio
	$ret_bp=system("perl -e 'use Bio::Root::Version;' 2>install_err.log > /dev/null");
	if($ret_bp ne 0){
		print "--please install the BioPerl Perl module\n";
    		
	}
	else{
		print "--BioPerl (OK)\n";
	}




	print "--prokka is not installed, but the installation script can take care of that!\n";


}
else{
	print "--prokka (OK)\n";
}





$ret=system("which blastp > /dev/null");
if($ret ne 0){

	print "--blastp is not installed, but the installation script can take care of that!\n";


}
else{
	print "--blastp (OK)\n";
}



$ret=system("which mummer > /dev/null");
if($ret ne 0){

	$ret_make=system("which make > /dev/null");
	if($ret_make ne 0){
		print "--please install make. I need it to compile mummer\n";
    		
	}
	else{
		print "--make (OK)\n";
	}

	$ret_csh=system("which csh > /dev/null");
	if($ret_csh ne 0){
		print "--please install csh. I need it to compile mummer\n";
    	
	}
	else{
		print "--csh (OK)\n";
	}



	$ret_gcc=system("which g++ > /dev/null");
	if($ret_gcc ne 0){
		print "--please install g++. I need it to compile mummer\n";
    	
	}
	else{
		print "--g++ (OK)\n";
	}



	#I install mummer
	print "--mummer is not installed, but the installation script can take care of that!\n";

}
else{
	print "--mummer (OK)\n";
}



$ret=system("which CONTIGuator.py 2>install_err.log > /dev/null");
if($ret ne 0){

	print "--Contiguator is not installed, but the installation script can take care of that!\n";


}
else{
	print "--Contiguator (OK)\n";
}



print "\n::If you know how to install the missing softwares/packages, go for it! Otherwise, do not worry! The installation script will guide you through the installation process.\n";



print "::Bye!\n";
