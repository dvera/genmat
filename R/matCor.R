matCor <-
function( mats, numgroups=3,hcluster=FALSE,plotcolors=c("red","green","black") , legendnames=basename(removeext(mats))){
	library(parallel)
	library(gplots)
	colramp <- colorRampPalette(plotcolors, space = "rgb")
	brks=seq(-1,1,by=2/100)
	hcolors<-colramp(100)
	matnames<-basename(removeext(mats))
	nummats<-length(mats)

	#make list of matrices
	cat("reading in matrices\n")
	matlist<-mclapply(1:nummats,function(x){  read.mat(mats[x]) },mc.cores=detectCores() )
	numcols<-ncol(matlist[[1]])
	#remove rows with no variance or no scores
	cat("finding genes with no variance\n")
	badvarrows<-unique(unlist(lapply(1:nummats, function(x){ which(apply(matlist[[x]],1,sd, na.rm=TRUE)==0) } )))

	cat("finding genes with no scores\n")
	badrows<-unique(unlist(lapply(1:nummats, function(x){ which(rowSums(is.na(matlist[[x]]))==numcols) })))

	badrows<-c(badvarrows,badrows)

	if(length(badrows)>0){
	    cat("removing bad rows\n")
	    matlist<-lapply(1:nummats, function(x){ matlist[[x]][-badrows,] })
	}

	numcols<-ncol(matlist[[1]])
	numgenes<-nrow(matlist[[1]])
	genenames<-row.names(matlist[[1]])

	#make pairwise matrix correlations
	cat("calculating global correlations\n")
	pairs<-expand.grid(1:nummats,1:nummats)
	matcors<-as.numeric(unlist(lapply(1:nrow(pairs),function(x) cor(as.vector(matlist[[pairs[x,1]]]) , as.vector(matlist[[pairs[x,2]]]) , use="complete.obs", method="spearman" ))))
	cormat<-matrix(nrow=nummats,ncol=nummats)
	row.names(cormat)<-legendnames
	colnames(cormat)<-legendnames
	for(i in 1:nrow(pairs)){
		cormat[pairs[i,1],pairs[i,2]]<-matcors[i]
	}
	globaltablename<-paste(matnames[1],"_globalcors.tsv",sep="")
	cat("saving global correlation table to",globaltablename,"\n")
	write.table(cormat,file=globaltablename,sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE)

	#plot pairwise global correlations as a heatmap
	globalhmname<-paste(matnames[1],"_matcors.pdf",sep="")
	cat("saving matrix cross-correlation heatmap to",globalhmname,"\n")
	pdf(file=globalhmname)
	heatmap.2(cormat,Rowv=NA,Colv=NA,dendrogram="none",trace="none",col=hcolors,breaks=brks,main="Global cross-correlations")

	#calculate pairwise gene correlations
	cat("calculating pairwise gene correlations\n")
	genecormat<-as.data.frame(
		lapply(1:nummats,function(y)
			unlist(
				lapply(
					1:numgenes, function(x){
						cor( matlist[[y]][x,],matlist[[1]][x,], use="complete.obs" )
					}
				)
			)
		)

	)
	row.names(genecormat)<-genenames
	colnames(genecormat)<-legendnames
	genecormat<-data.matrix(genecormat)

	genemeta<-data.frame(
	     "symbol"=unlist(lapply(strsplit(row.names(genecormat),"_"), "[", 3 ) ),
	     "id"=unlist(lapply(strsplit(row.names(genecormat),"_"), "[", 1 ) ),
	     "chrom"=unlist(lapply(strsplit(row.names(genecormat),"_"), "[", 2 ) ),
	     stringsAsFactors=FALSE
	)

	goodrows<-complete.cases(genecormat)
	genecormat<-genecormat[goodrows,]
	genemeta<-genemeta[goodrows,]

	#kmeans cluster and plot gene correlations
	cat("kmeans clustering data and plotting heatmap\n")
	k<-kmeans(genecormat,numgroups)
	groupcolors<-unlist(lapply(k$cluster,function(x) rainbow(numgroups)[x]))[order(k$cluster)]
	heatmap.2(genecormat[order(k$cluster),],trace="none",Colv=NA,Rowv=NA,dendrogram="none",RowSideColors=groupcolors,col=hcolors,breaks=brks,labRow=NA,margins=c(10,10),cexCol=1,main=paste("clustered gene-by-gene correlations with",legendnames[1]))
	cat("plotting sorted heatmap\n")
	for(i in 2:nummats){
		heatmap.2(genecormat[order(genecormat[,i]),],trace="none",Colv=NA,Rowv=NA,dendrogram="none",col=hcolors,breaks=brks,labRow=NA,main=paste("gene-by-gene correlations with ",legendnames[1],"sorted on",legendnames[i]))
	}
	if(hcluster==TRUE){
		cat("heirarchical clustering data and plotting heatmap\n")
		heatmap.2(genecormat,Colv="none",dendrogram="row",trace="none",col=hcolors,breaks=brks,labRow=NA,main=paste("hclustered gene-by-gene correlations with ",legendnames[1]))
	}

	#save file
	rbc<-rainbow(nummats-1)
	plot(density(genecormat[,2],from=-1,to=1),col=rbc[1],xlim=c(-1,1),ylim=c(0,10),xlab=paste("gene-by-gene correlation with",legendnames[1]))
	if(nummats>2){
		for(i in 3:(nummats)){
			lines(density(genecormat[,i],from=-1,to=1),col=rbc[i-1])
		}
	}
	legend("topleft",legend=legendnames[2:nummats],col=rbc,lwd=3)

	par(mar=c(7,5,1,1))
	boxplot(genecormat[,2:nummats],ylab=paste("gene-by-gene correlation with",legendnames[1]),las=2)

	dev.off()
	genecormat<-cbind(genecormat,genemeta)
	write.table(genecormat,file=paste(matnames[1],"_gene_correlations.tsv",sep=""),sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE)
}
