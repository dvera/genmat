mat.plotaverages2 <-
function( matfilelist,suffixes,legendnames,plotcolors=rainbow(length(matfilelist)), windowsize=10, matcols=100, xname="Distance from Feature", yname="Average Score"){

	
	# coplot as many data sets as possible

	numlists<-length(matfilelist)
	if(length(suffixes==1)){suffixes<-rep(suffixes,numlists)}
	prefixes<-lapply(1:numlists, function(x){
		remove.suffix(basename(matfilelist[[x]]),suffixes[x])
	})
	feats<-unique(unlist(prefixes))
	matches<-lapply(1:numlists,function(x) match(feats,prefixes[[x]]) )
	pdfname<-paste("mataverages_",length(feats),"centers_",numlists,"scores.pdf",sep="")
	pdf(file=pdfname)
	xaxis<-windowsize*(1:matcols) - (windowsize*matcols)/2
	for(i in 1:length(feats)){
		cat("plotting feature #",i,"/",length(feats),":",feats[i],"\n")
		refs<-unlist(lapply(matches,"[",i))
		nomats<-which(is.na(refs))
		goodmats<-which(1:numlists %ni% nomats)
		curmats<-lapply(goodmats,function(x) read.mat(matfilelist[[x]][matches[[x]][i]]))
		counts<-unlist(lapply(curmats,nrow))
		colmeans<-lapply(curmats,colMeans,na.rm=TRUE)
		all<-unlist(colmeans)
		ylims<-c(min(all),max(all))
		plot(0,type="n",xlim=c(min(xaxis),max(xaxis)),ylim=ylims,main=feats[i],xlab=xname,ylab=yname)
		n<-rep(0,numlists)
		n[goodmats]<-counts
		for(j in 1:numlists){
			if(j %in% goodmats == TRUE){
				lines(xaxis,colmeans[[which(goodmats==j)]],col=plotcolors[j],lwd=3)
			}

		}
		legend("topleft",legend=paste(legendnames,n),col=plotcolors,lwd=3)
	}
	dev.off()
	cat("pdf saved to",pdfname,"\n")
}
