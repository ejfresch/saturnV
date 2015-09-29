#!/usr/bin/Rscript

library(ggplot2)



args=commandArgs(trailingOnly = TRUE)

dir_genomes=args[1]

#cat("--",dir_genomes,"--","\n",sep="")

if(is.na(dir_genomes)){
    stop("::check_distrib_scaffolds.R <dir_genomes>\n")
  
}

genomes=list.files(dir_genomes,pattern=".fasta")

fragments=c()

for(i in 1:length(genomes)){

    file=paste(dir_genomes,genomes[i],sep="/")
    cmd=paste("cat ",file,' |grep ">"|wc -l',sep="") 
    #cat(cmd)   
    no_of_frag=system(cmd,intern=T)
    cat("--",genomes[i],":",no_of_frag," fragment(s)\n",sep="")
    fragments=c(fragments,no_of_frag)

}


cat=rep("A",length(fragments))

df1=data.frame(cat, fragments)


#print statistics
cat("::quantiles\n")
quantile(as.numeric(fragments))
cat("\n")
cat("::min:\n")
min(as.numeric(fragments))
cat("\n")
cat("::max:\n")
max(as.numeric(fragments))
cat("\n")



#produce a graph with the distribution of the fragments
pdf("distribution_fragments_genomes.pdf")
ggplot(df1, aes(x=fragments)) + geom_histogram(colour="black")

dev.off()


#produce a file with the results

df2=data.frame(genomes,fragments)

write.csv(df2, file="distribution_fragments_genomes.csv", quote=F, row.names=F)


