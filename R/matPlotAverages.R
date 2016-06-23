matPlotAverages <-
function(
	matrixlist,
	xlabel="Distance from feature (bp)",
	ylabel="Mean score",
	ylims=c("auto","auto"),
	plotcolors=NULL,
	legendnames=NULL,
	threads=getOption("threads",1L),
	plottype="l",
	linewidth=2,
	drawlegend=TRUE,
	mainlabels=NULL,
	legendpos="topleft",
	...
){

	if(!is.list(matrixlist)){
		matrixlist <- list(matrixlist) 
	        if(!is.null(mainlabels)){names(matrixlist) = mainlabels}
	}

	if(is.null(mainlabels)){
		if(!is.null(names(matrixlist))){
			mainlabels=names(matrixlist)
		} else{
			mainlabels=rep("",length(matrixlist))
		}
	}

	if(is.null(plotcolors)){
		plotcolors=lapply(1:length(matrixlist),function(m){ rainbow(length(matrixlist[[m]]))})
	} else if(!is.list(plotcolors)){
		plotcolors=list(plotcolors)
	}

	if(is.null(legendnames)){
		legendnames=lapply(1:length(matrixlist),function(m){ basename(removeext(matrixlist[[m]])) })
	} else if(!is.list(legendnames)){
		legendnames=list(legendnames)
	}

	for(m in 1:length(matrixlist)){

		# count matrices
		nummats<-length(matrixlist[[m]])

		# read in matrices
		matlist<-matRead(matrixlist[[m]],unlistSingle=FALSE)
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
		windowsizes<-as.numeric(gsub("mat","",file_ext(matrixlist[[m]])))
		xs<-lapply(1:nummats,function(x) ((1:matcols[x])-(matcols[x]/2))*windowsizes[x])

		# plot window
		plot(
	    0,
			type="n",
	    xlim=c(min(xs[[1]]),max(xs[[1]])),
			ylim=ylims,
	    xlab=xlabel,
	    ylab=ylabel,
			main=names(mainlabels[[m]]),
	    ...
	  )

		#plot averages
		for(k in 1:nummats){
			lines(
			  xs[[k]],
			  mataverages[[k]],
	      col=plotcolors[[m]][k],
	      lwd=linewidth,
	      type=plottype
			)
	  }

		if(drawlegend){legend(legendpos,legend=legendnames[[m]], col=plotcolors[[m]], lwd=linewidth, cex=0.75)}

	}
}
