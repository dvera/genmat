mat.plotmats <-
function( matlist , tracknames , colorlist="auto", plottypes = "l" , matnamelist="auto", ymaxs = "auto" , ymins = "auto" , cores = "max", centername="TSS" ){
	library(tools)
	library(parallel)
	if(cores=="max"){cores<-detectCores()-1}
	
	numsets<-length(matlist)
	nummats<-unlist(lapply(matlist,length))
	
	cat("loading matrices\n")
	mats<-lapply(1:numsets, function(s) lapply(1:nummats[s], function(m) read.mat( matlist[[s]][m] ) ) )
	
	if(matnamelist=="auto"){ matnamelist<-lapply(1:numsets, function(s) basename(removeext(matlist[[s]] ) ) ) }
	if(colorlist=="auto"){ colorlist<-lapply(1:numsets, function(s) rainbow(nummats[s]) ) }
	if(length(plottypes)==1){plottypes=rep(plottypes,numsets)}
	
	windowsize=as.numeric( gsub("mat","",file_ext(matlist[[1]][1]) ) )
	matcols<-ncol(mats[[1]][[1]])
	matrows<-nrow(mats[[1]][[1]])
	x<-((1:matcols)-(matcols/2))*windowsize
	xw<-c(min(x),((max(x)-min(x))*0.25+max(x)))
	
	
	cat("calculating y maxs\n")
	matmax<-lapply(1:numsets, function(s) lapply(1:nummats[s], function(m) apply(mats[[s]][[m]],1,max,na.rm=TRUE)  )  )
	yautomax<-lapply(1:numsets, function(s) unlist(mclapply(1:matrows, function(r) max(unlist(lapply(matmax[[s]],"[",r)), na.rm=TRUE) ,mc.cores=cores ) ) )
	
	cat("calculating y mins\n")
	matmin<-lapply(1:numsets, function(s) lapply(1:nummats[s], function(m) apply(mats[[s]][[m]],1,min,na.rm=TRUE)  )  )
	yautomin<-lapply(1:numsets, function(s) unlist(mclapply(1:matrows, function(r) min(unlist(lapply(matmin[[s]],"[",r)), na.rm=TRUE) ,mc.cores=cores ) ) )
	
	pdf(file="plotmatstest.pdf")
	par(mfrow=c(numsets,1),oma=c(4,0,4,0),mar=c(0,4,0,4))
	
	cat("drawing plots\n")
	for(r in 1:matrows){
		for(s in 1:numsets){
			plot(
				0,
				type="n",
				xaxt=if(s==numsets){'s'} else{'n'},
				xlim=xw,
				ylim=c(yautomin[[s]][r],yautomax[[s]][r]),
				ylab=tracknames[s],
				xlab=paste("Distance from",centername,"(bp)")
				#main=row.names(mats[[s]][[m]])[r]
			)
			grid(lty="solid",ny=NA)
			abline(v=max(x),lwd=4,col="black")
			for(m in 1:nummats[s]){
				lines(x,mats[[s]][[m]][r,],col=colorlist[[s]][m],type=plottypes[s])
			}
			legend("topright",legend=matnamelist[[s]],col=colorlist[[s]],lwd=3, cex=0.5, bty="n")
		}
	}
	dev.off()
}
