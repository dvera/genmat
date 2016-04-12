matIqrNorm <-
function( mats, normalizeto=NULL, cores ="max" ){

  library(tools)
  library(parallel)

  nummats<-length(mats)
  if(cores=="max"){cores=detectCores()-1}
  if(cores > nummats){cores=nummats}

  # get mat info
  matnames<-basename(removeext(mats))
  exts<-file_ext(mats)
  outnames <- paste0( matnames, "_iqrnorm.", exts )

  # read matrices
  matlist<-mclapply(mats,matRead,mc.cores=cores)

  # calculate IQRs
  iqrs <- unlist(mclapply(matlist,IQR,na.rm=T, mc.cores=cores))

  # set target IQR
  if(is.null(normalizeto)){
    iqr <- mean(iqrs)
  } else{
    iqr <- IQR( matlist[normalizeto], na.rm=T )
  }

  # calculate scalars for normalizaion
  scalars <- iqr/iqrs

  # normalize matrices
  matlist<-mclapply( 1:nummats,function(x) {
      nmat <- matlist[[x]] * scalars[x]
      matWrite(nmat,file=outnames[x])
  }, mc.cores=cores )

  # retrn output file names
  return(outnames)

}
