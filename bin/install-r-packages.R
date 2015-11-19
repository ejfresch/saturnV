#!/usr/bin/Rscript

options(warn=-1)

packages_to_check=c("ggplot2","gplots")


for(package in packages_to_check){


     if(!suppressMessages(require(package,character.only = TRUE))){
        install.packages(package,repos='http://cran.stat.unipd.it/')

    }
    



}
