matWrite <-
function( mat, ... ){
	write.table(mat,sep="\t",quote=FALSE,col.names=FALSE, ... )
}
