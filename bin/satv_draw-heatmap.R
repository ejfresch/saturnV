#!/usr/bin/Rscript



args=commandArgs(trailingOnly = TRUE)


bin_matrix=args[1]

if(is.na(bin_matrix)){
    stop("::satv_draw-heatmap.R <bin_matrix>\n")
  
}



library(gplots)
a=read.csv(bin_matrix,sep="\t",stringsAsFactors=F)

df=a[,2:ncol(a)]

eval=apply(df,1,sum)

n_gen=ncol(df)

df$eval=eval

df2=df[ order(-df[,c("eval")]), ]



df3=data.matrix(df2[,1:n_gen])

mypalette=c("black","red")


png("heatmap_overview_pangenome.png")
heatmap.2(df3,Rowv=FALSE, Colv=TRUE, dendrogram="column", col=mypalette,breaks=c(-0.5,0.5,1.5),density.info="none",trace="none",cexRow=0.9,cexCol=0.9)
		
           
dev.off()



