mat.window <-
function( matfile, operation="mean", windowsize=5){
    library(tools)
    if(ceiling(windowsize/2)==floor(windowsize/2)){stop("windowsize must be an odd number")}
    cat("reading in matrix\n")
    mat<-read.mat(matfile)
    
    cat("gathering info\n")
    flank<-(windowsize-1)/2
    numcols<-ncol(mat)
    matrownames<-row.names(mat)
    
    print(numcols)
    
    cat("windowing matrix\n")
    mat<-simplify2array(lapply(1:numcols,function(x){
	if(x-flank<1){l=1}
	else{l=x-flank}
	if(x+flank>numcols){r=numcols}
	else{r=x+flank}
	rowMeans(mat[,l:r],na.rm=T)
    }))
    
    cat("wrapping up\n")
    row.names(mat)<-matrownames
    ext<-file_ext(matfile)
    outname<-paste(basename(removeext(matfile)),"_",operation,"win",windowsize,".",ext,sep="")
    write.mat(mat,file=outname)
    return(outname)
}
