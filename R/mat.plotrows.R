mat.plotrows <-
function( mats,genes,geneid="symbols",sep="_",prunerows=FALSE,lspan=0,normalize=F,legendnames=basename(removeext(mats)),ylims=c(0,"auto"),plotcolors=rainbow(length(mats)),filename="orphan",start="beginning",stop="end",center="TSS",allforward=TRUE,linewidth=3,blockheight=20,grid=TRUE,X=F){
	library(tools)
	library(parallel)
	nummats<-length(mats)

	#read in matrices and get info
	cat("loading matrices...\n")
	matlist<-mclapply(mats,read.mat,mc.cores=detectCores())
	matcols<-unlist(lapply(matlist,ncol))
	numgenes<-length(genes)
	if(numgenes>16 & X==T){X=F;cat("more than 16 genes, forcing X=FALSE\n")}
	matnames<-basename(removeext(mats))
	genenames<-unlist(lapply(strsplit(row.names(matlist[[1]]),sep),"[",4))
	windowsizes=as.numeric( gsub("mat","",file_ext(mats) ) )
	xs<-lapply(1:length(mats),function(x) ((1:matcols[x])-(matcols[x]/2))*windowsizes[x])
	outname<-uniquefilename("mat_plotrows.pdf")

	if(prunerows==TRUE){
		cat("finding empty rows\n")
		badrows<-unique(unlist(lapply(1:nummats, function(x){ which(rowSums(is.na(matlist[[x]]))==matcols[x]) })))
		if(length(badrows)>0){
			cat("removing empty rows\n")
			matlist<-lapply(1:nummats,function(x){
				matlist[[x]][-badrows,]
			})
		}
	}


	#normalize data
	if(normalize==TRUE){
		cat("normalizing data\n")
		normat<-function(curmat){t(simplify2array(lapply(1:nrow(curmat),function(x) curmat[x,]/mean(curmat[x,],na.rm=TRUE))))}
		matlist<-lapply(matlist,normat)
	}

	#find ylims
	cm=1
	if(ylims[1]=="auto" | ylims[2]=="auto"){
		cat("gathering data to calculate ylims\n")
		cm<-do.call(cbind,matlist)
	}

	if(ylims[1]=="auto"){
		cat("calculating ymins\n")
		ymins<-apply(cm,1,min,na.rm=TRUE)
		#print(ymins)
	}
	else{
		ymins<-as.numeric(rep(ylims[1],length(genenames)))
	}
	if(ylims[2]=="auto"){
		cat("calculating ymaxs\n")
		ymaxs<-1.1*apply(cm,1,max,na.rm=TRUE)
		#ok, shprint(ymaxs)
	}
	else{
		ymaxs<-as.numeric(rep(ylims[2],length(genenames)))
	}
	rm(cm)

	#identify genes to plot
	cat("Finding genes\n")
	if(geneid=="symbols"){
		matrows<-unlist(lapply(1:numgenes,function(x) which(genenames==genes[x])[1]))
	}
	if(geneid=="ucsc"){
		matrows<-unlist(lapply(1:numgenes,function(x) which(genenames==genes[x])[1]))
	}
	if(geneid=="rows"){
		matrows<-genes
	}

	if(X==F){
		pdf(file=outname)
	}
	else{
		parnum<-ceiling(sqrt(numgenes))
		par(mfrow=c(parnum,parnum))
	}

	#plot averages
	plot(
			0,
			type="n",
			xlim=c(min(xs[[1]]),max(xs[[1]])),
			ylim=c(mean(ymins),mean(ymaxs)),
			xlab="Distance from TSS (bp)",
			ylab="Average FPM",
			main=paste("Average score for",numgenes,"genes"),
			cex.lab=1.5,
			cex.axis=1.5,
			cex.main=1.5,
			cex.sub=1.5
	)

	for(m in 1:length(mats)){
			x<-1:matcols[m]
			y<-colMeans(matlist[[m]][matrows,],na.rm=TRUE)
			lines(
				xs[[m]],
				y,
				col=plotcolors[m],
				lwd=3
			)
	}

	legend("topright",legend=legendnames,col=plotcolors,lwd=3)



	for(i in 1:numgenes){
		plot(
			0,
			type="n",
			xlim=c(min(xs[[1]]),max(xs[[1]])),
			ylim=c(ymins[matrows[i]],ymaxs[matrows[i]]),
			xlab="Distance from TSS (bp)",
			ylab="Normalized FPM",
			main=genenames[matrows[i]],
			cex.lab=1.5,
			cex.axis=1.5,
			cex.main=1.5,
			cex.sub=1.5
		)

		for(m in 1:length(mats)){
			if(lspan!=0){
				x<-1:matcols[m]
				y<-matlist[[m]][matrows[i],]
				y[-which(is.na(y))]<-loess(y~x,span=lspan)$fitted
			}
			else{
				y<-matlist[[m]][matrows[i],]
			}
			lines(
				xs[[m]],
				y,
				col=plotcolors[m],
				lwd=3
			)
		}

		legend("topright",legend=legendnames,col=plotcolors,lwd=3)
	}
	cat("pdf saved to",outname,"\n")

	if(X==F){dev.off()}
}
