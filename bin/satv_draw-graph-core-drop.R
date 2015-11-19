#!/usr/bin/Rscript

library(ggplot2)


args=commandArgs(trailingOnly = TRUE)


bin_matrix=args[1]

if(is.na(bin_matrix)){
    stop("::satv_draw-graph-core-drop.R <bin_matrix>\n")
  
}


cat("::reading input file\n",sep="")
df=read.csv(bin_matrix,sep="\t",stringsAsFactors=F)

vect_core=c()
vect_shared=c()
vect_uniq=c()



for(i in 3:ncol(df)){
    
    cat("::iteration: ",i,"/",(ncol(df)-1),"\n",sep="")

    current_df=df[,2:i]
    eval=apply(current_df,1,sum)
    
    core=eval[eval==(i-1)]
    uniq=eval[eval==1]
    shared=eval[eval>1 & eval<(i-1)]

    vect_core=c(vect_core,length(core))
    vect_shared=c(vect_shared,length(shared))
    vect_uniq=c(vect_uniq,length(uniq))
    


}



pdf("graph_core_accessory_genome.pdf")
plot(vect_shared,type="l",col="blue",lwd=3)
lines(vect_core,col="red",lwd=3)
lines(vect_uniq,col="black",lwd=3)

dev.off()


pdf("graph_core_genome.pdf")
plot(vect_core,type="l",col="red",lwd=3)


dev.off()



##I make the graphs with ggplot

vx=seq(1,length(vect_core))

vect_y=c(vect_core,vect_shared,vect_uniq)
vect_x=rep(vx,3)

vect_f=c(rep("core",length(vect_core)),rep("shared",length(vect_shared)),rep("unique",length(vect_uniq)))

df_g=data.frame(vect_x,vect_y,vect_f)

pdf("graph_core_accessory_genome_ggplot.pdf")
ggplot(data=df_g, aes(x=vect_x, y=vect_y, group=vect_f, colour=vect_f)) + 
    geom_line(size=1.5) + labs(x = "No. of genomes", y="No. of genes")+scale_color_manual(values=c("red", "orange", "black"))
dev.off()


df_gc=data.frame(vx,vect_core)

pdf("graph_core_ggplot.pdf")
ggplot(data=df_gc, aes(x=vx, y=vect_core)) + 
    geom_line(size=1.5) + labs(x = "No. of genomes", y="No. of genes")
dev.off()






