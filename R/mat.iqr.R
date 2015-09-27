mat.iqr <-
function( mats, normalizeto=NULL, cores ="max" ){

  library(tools)
  library(parallel)

  nummats<-length(mats)
  if(cores=="max"){cores=detectCores()-1}
  if(cores > nummats){cores=nummats}
  
  # calculate IQRs
  iqrs <- unlist(mclapply(1:nummats, function(x) IQR( read.mat(mats[x]), na.rm=T), mc.cores=cores ))

  # retrn output file names
  return(iqrs)

}
