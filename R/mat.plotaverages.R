mat.plotaverages <-
function(
	matrixlist,
	xlabel="Distance from feature (bp)",
	ylabel="Mean score",
	ylims=c("auto","auto"),
	plotcolors=rainbow(length(matrixlist)),
	legendnames=basename(removeext(matrixlist)),
	threads=getOption("threads",1L),
	plottype="l",
	linewidth=2,
	...
){

	# count matrices
	nummats<-length(matrixlist)

	# read in matrices
	matlist<-mclapply(matrixlist,read.mat,mc.cores=threads)
	matcols<-unlist(lapply(matlist,ncol))

	# get column averages for each matrix
	mataverages<-lapply(matlist,colMeans,na.rm=TRUE)

	# find min and max across all matrices
	ymax<-max(unlist(mataverages),na.rm=TRUE)
	ymin<-min(unlist(mataverages),na.rm=TRUE)

	# set default y lims
	if(ylims[1]=="auto"){ylims[1]=ymin}
	if(ylims[2]=="auto"){ylims[2]=ymax}
	ylims<-as.numeric(ylims)

	# find bp size of matrix
	windowsizes<-as.numeric(gsub("mat","",file_ext(matrixlist)))
	xs<-lapply(1:nummats,function(x) ((1:matcols[x])-(matcols[x]/2))*windowsizes[x])

	# plot window
	plot(
    0,
		type="n",
    xlim=c(min(xs[[1]]),max(xs[[1]])),
		ylim=ylims,
    xlab=xlabel,
    ylab=ylabel,
    ...
  )

	#plot averages
	for(k in 1:nummats){
		lines(
		  xs[[k]],
		  mataverages[[k]],
      col=plotcolors[k],
      lwd=linewidth,
      type=plottype
		)
  }

	legend("topleft",legend=legendnames, col=plotcolors, lwd=linewidth, cex=0.75)

}
