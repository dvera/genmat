mat2bg <- function( matfiles , binsize , maxbins , prefix , chromnames = remove.prefix( matfiles , "\\.") , cores="max"){

  options(scipen=9999999)
  if(cores=="max"){ cores <- detectCores()-1 }

  bins <- mclapply(1:length(matfiles), function(i){

    mat <- data.matrix(read.tsv(matfiles[i]))

    if(ncol(mat) != nrow(mat)){
      cat("WARNING: MATRIX",matfiles[i],"DOES NOT HAVE EQUAL DIMENSIONS, TRIMMING\n")
      mindim <- min(dim(mat))
      mat <- mat[ 1:mindim , 1:mindim ]
    }

    numbins <- ncol(mat)

    if( maxbins > (numbins-1)){ maxbins <- numbins-1 }

    flankleft  = rep( c(0:maxbins), each=2 )
    flankright = rotvec( flankleft )
    bincount <- length(diag(mat))

    y <- data.matrix(as.data.frame(lapply(1:maxbins, function(w){
      c( rep(NA, flankleft[w] ) , diag(mat[,w:numbins]) , rep( NA, flankright[w] ) )
    })))

    y[y==0] <- NA


    bgs <- lapply(1:maxbins, function(w){

      bg <- data.frame(
        chromnames[i] ,
        ((0:bincount)*binsize+1)[1:bincount] ,
        (1:bincount)*binsize ,
        y[,w]
      )
      bg[seq(2,bincount,2),2:3] <- bg[seq(2,bincount,2),2:3] - (binsize/2)
      bg <- bg[complete.cases(bg),]
      colnames(bg) <- paste0("V",1:4)
    })

    return(bgs)

  }, mc.cores=cores, mc.preschedule=FALSE)

  distwidth <- nchar(maxbins * binsize)
  dump<-mclapply(1:maxbins, function(w){
    bg<-do.call(rbind,lapply(bins,"[[",w))
    bg<-bg[order(bg[,1],bg[,2]),]
    outname<-paste0(prefix,"_d",sprintf(paste0("%0",distwidth,"d"), binsize*(w-1)),".bg")
    write.tsv(bg,file=outname)
  }, mc.cores=cores, mc.preschedule=FALSE)

}
