matRead <- function( mat, fmat=FALSE , threads=getOption("threads",1L), unlistSingle=TRUE , ... ){
	#if(!fmat){
    m <- mclapply(mat,read.tsv, row.names=1 , mc.cores=threads, mc.preschedule=FALSE,  ... )
		m <- lapply(m, as.matrix)
		if(length(m)==1 & unlistSingle){m<-m[[1]]}
    return(m)
  # } else{
  #   as.matrix(read.table(mat,stringsAsFactors=FALSE,sep="\t",row.names=1, ... ))
  # }

}
