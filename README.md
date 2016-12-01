Copyright 2015-2017 -- Luca Freschi, Antony T. Vincent, Julie Jeukens, Jean-Guillaume Emond-Rheault, Irena Kukavica-Ibrulj, Marie-Josée Dupont, Steve J. Charette and Roger C. Levesque

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


TABLE OF CONTENTS

0. Introduction
1. Installation
2. How to determine the core and accessory genomes
3. How to cite SaturnV


0 - INTRODUCTION
------------------
SaturnV is a software that allows to compare two or more genomes in order to determine and study the core and accessory genes. The core genes are the set of genes that are present in all genomes of a given dataset. Unique genes (those present in one genome only) and flexible genes (those present in more than one genome) constitute the accessory genes.

SaturnV takes .fasta files as input (assembled genomes) and produces a matrix of presence/abscence for all the genes that constitute the pangenome of the species considered.

SaturnV is named after the rocket that brought humankind to the Moon in 1969. A bit like its space couterpart, SaturnV relies on several modules that have to function one after the other in order to achieve the objectives of the mission: determining, comparing and displaying core and accessory genes.


SaturnV was mainly written in Perl. Some parts are written in R and python.


1 - INSTALLATION
-----------------
Here is the complete list of dependencies that we are asking you to install in order to run the core module of SaturnV:
* perl + 3 perl libraries (Getopt::Long, Parallel::ForkManager, Graph)
* prokka
* usearch (v8.x)
* GNU parallel

Other softwares you might want to insall:
* last
* blastp
* diamond

NOTE: SaturnV will use usearch as its default method to compare sequences. However, the same task can be performed by other softwares (usearch, blastp, diamond). SaturnV is ready to allow you to choose the software you are more confortable with.

Here is an example of a typical installation process:

--step1: copy the saturnV repository to your computer. Open a terminal emulator, and type:
```
git clone https://github.com/ejfresch/saturnV.git
```

NOTE: you should have git installed! If it is not the case, please install it (hint for Debian/Ubuntu systems: sudo apt-get install git-core).

go into the directory saturnV with the command:
```
cd saturnV
```

--step 2: run the installer and check if some dependencines are missing. If so, please install them and then rerun the installer
```
./install -d /home/lfreschi/sw/saturnv/ -m core
```

--step3: check that everything is OK.
to check that everything went well, type satv_ and then hit the tab key twice. You should be able to see the whole list of SaturnV commands:

```
lfreschi@katak:~/tasks/pangenome/test_saturn/saturnv> satv_
satv_generate_list_genomes              satv_search-pangenome-centroids         satv_search-pangenome-lazy-best-hit
satv_launch                             satv_search-pangenome-laziest           satv_search-pangenome-strictest
satv_prodigal-driver                    satv_search-pangenome-laziest-best-hit  satv_trim-name-scaffolds
satv_prokka-driver                      satv_search-pangenome-lazy 
[...the list continues if you installed two or more modules...]
```

If you are in trouble or you think you have followed all these steps but you just got some complicated error message, just write a message to l.freschi@gmail.com



2 - HOW TO DETERMINE THE CORE AND ACCESSORY GENOMES
----------------------------------------------------

SaturnV comes with a small dataset of Achromobacter genomes (n = 3). You can find them into the directory examples/achromo/. We will use this dataset to go through the SaturnV commands and show you how you can use SaturnV to determine core and accessory genomes.

--step1: create a working directory for the new analysis and create the subdirectory genomes/. Then copy the genomes you want to analyze into genomes/:

I create the working directory
```
mkdir satv_results
```

I move into it
```
cd satv_results
```

I create the genomes/ directory
```
mkdir genomes
```

I copy the genoms inside the genome directory
general synthax: cp <directory_where_you_installed_saturnV>/examples/achromo/* genomes/
in my case:
```
cp ~/sw/saturnv/examples/achromo/* genomes/
```


--step2: launch the analysis

general synthax: satv_launch.pl -d <directory_genomes_to_analyze> -c <cpus_available_for_multithreading> -ann <annotation_software[prokka|prodigal]> -m <clustering_method[lazy|strict|strictest|centroids]> -a <search_algorithm[usearch|blast|last|diamond]> -i <min_perc_identity_orthologs> -ip <min_perc_identity_paralogs>
```
satv_launch.pl -d genomes/ -c 2 -ann prodigal -m lazy -a usearch -i 50 -ip 100
```


--step3: look at the results

the main output is the table_linked5_<method>.tsv. When we launched saturnV, we specified to use the strict method for the clustering step, so the file name will be table_linked5_strict.tsv.

This file is a tab separated value (.tsv) file. In each row there is a gene and in each column there is its ortholog in another genome.

Another output is the annotation. By typing the command "ls" you will see that there are 3 folders, each one with the name of one of the genomes we analyzed. These folders contain the annotation.

NOTE: all files MUST have the extension .fasta

3 - HOW TO CITE SATURNV
-----------------------
If you use SaturnV for your publications, please 
