matHist <-
function( mats, xlims=c("auto","auto"), plotcolors=rainbow(length(mats)), legendnames=basename(removeext(mats)) ){
	library(parallel)
	nummats<-length(mats)
	# todo: indicate dimensions and # of NAs in matrix
	cat("reading in matrices\n")
	matlist<-mclapply( mats, read.mat, mc.cores=detectCores() )

	cm=1
	if(xlims[1]=="auto" | xlims[2]=="auto"){
		cat("calculating xlims\n")
		cm<-do.call(cbind,matlist)
		q<-quantile(cm,probs=c(0.05,0.95),na.rm=T)
	}
	if(xlims[1]=="auto"){xlims[1]=q[1]}
	if(xlims[2]=="auto"){xlims[2]=q[2]}
	xlims<-as.numeric(xlims)
	print(xlims)

	cat("calculating score distributions\n")
	densitylist<-lapply( 1:nummats, function(x) density( matlist[[x]], from=xlims[1], to=xlims[2], na.rm=T))
	ymax<-max(unlist(lapply( 1:nummats, function(x) density( matlist[[x]], from=xlims[1], to=xlims[2],na.rm=T)$y  )),na.rm=TRUE)
	print(ymax)
	cat("plotting score distributions\n")
	plot(0,type="n",ylim=c(0,ymax), xlim=xlims, xlab="Fragment Size", ylab="Density")
	for(i in 1:nummats){
		lines(densitylist[[i]],col=plotcolors[i],lwd=2)
	}
	legend("topright",legend=legendnames, col=plotcolors, lwd=3)
}
