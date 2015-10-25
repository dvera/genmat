matRead <-
function( mat, fmat=FALSE , ... ){
	if(!fmat){
    m <- read_tsv(mat,col_names=FALSE, ... )
    rownames(m) <- m[,1]
    m <- data.matrix(m[,-1])
  } else{
    as.matrix(read.table(mat,stringsAsFactors=FALSE,sep="\t",row.names=1, ... ))
  }

}
