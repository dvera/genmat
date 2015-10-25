matQuantileNorm <-
function( mats, normalizeto=NULL, cores ="max" ){
	library(tools)
	library(parallel)
	nummats<-length(mats)
	if(cores=="max"){cores=detectCores()-1}
	if(cores > nummats){cores=nummats}
	matnames<-basename(removeext(mats))
	exts<-file_ext(mats)
	matlist<-mclapply(mats,read.mat,mc.cores=detectCores())

	if(is.null(normalizeto)){
		scores <- as.vector(unlist(matlist))
		scores <- na.omit(scores)
	} else{
		scores <- rep(as.vector(unlist(matlist[normalizeto])),10)
		scores <- na.omit(scores)
	}

	matlist<-mclapply( 1:nummats,function(x) {
			numscores<-length(which(is.na(as.vector(matlist[[x]]))==FALSE))
			samples <- sort(rowMeans(as.data.frame(lapply(1:10, function(x) sort(sample(scores,numscores))))))
			goodcells<-which(is.na(matlist[[x]])==FALSE)
			matlist[[x]][goodcells[order(matlist[[x]][goodcells])]]<-samples
			matlist[[x]]
	}, mc.cores=cores )

	outnames<-paste(matnames,"_qnorm.",exts, sep="")
	mclapply(1:nummats, function(x) write.mat(matlist[[x]],file=outnames[x]) ,mc.cores=detectCores() )
	return(outnames)
}
