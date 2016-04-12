matOps <-
function( matlist1 , matlist2=NULL , operation="difference" , outnames = NULL , pattern=NULL , replacement=NULL ){

	library(tools)


	#if(pattern=="" | replacement==""){stop("YOU MUST SPECIFY 'pattern' AND 'replacement' FOR OUTPUT FILE NAMES")}

	# argument check
	if( is.null(outnames)==FALSE & is.null(pattern)==FALSE & is.null(replacement)==FALSE ) { stop("must use pattern/replacement OR outnames, not both")}
	if( length ( operation ) != 1){stop("can only perform 1 operation")}

	if( is.null(outnames)){
		matlist1ext<-file_ext(matlist1)
		matlist2ext<-file_ext(matlist2)
		outnames<-gsub(pattern,replacement,basename(matlist1))
		nametab<-data.frame("MATRIX1"=basename(matlist1),"MATRIX2"=basename(matlist2),"OUTPUT"=outnames,stringsAsFactors=FALSE)
		print(nametab)
		if(length(which(basename(matlist1) %in% outnames)==TRUE)>0){stop("AT LEAST ONE OUTPUT FILE WILL REPLACE AN INPUT FILE")}
	}


	if( length(matlist1)==1 & is.null(matlist2) ){
		if(operation %ni% c("log2","log10","antilog2","antilog10","inverse")) { stop("operations that can be performed on a single matrix include log2, log10, antilog2, antilog10, and inverse")}
		if(is.null(outnames)){outnames <- paste0 (basename(removeext(matlist1)) , "_",operation,".",matlist1ext)}
		mat<-read.mat(matlist1)
		if(sum(grepl(0,mat))>0){cat("WARNING: MATRIX CONTAINS ZEROES THAT WILL BE INFINITE IF CALCULATING A LOGARITHMIC\n")}
		if(operation=="log2"){ mat<-log2(mat) }
		if(operation=="log10"){ mat<-log10(mat) }
		if(operation=="antilog2"){ mat<-2^(mat) }
		if(operation=="antilog10"){ mat<-2^(mat) }
		if(operation=="inverse"){ mat<-1/(mat) }
		matWrite(mat,file=outnames)
	}


	if( length(matlist1)>1 & is.null(matlist2) ){
		if( is.null(outnames) & is.null(pattern) & is.null(replacement) ) { stop("must use pattern/replacement OR outnames")}
		if(operation != "mean") { stop("you can only calculate the mean on a single set of matrices")}
		if(is.null(outnames)){stop("must specify output matrix filename. be sure to match the extension of the inputs")}
		if(length(outnames)!=1){stop("must have only 1 string in 'outname'")}
		mats<-lapply(matlist1,read.mat)
		if(length(unique(lapply(lapply(lapply(mats,dim),as.character),paste,collapse="-"))) != 1 ){stop("dimensions of matrices are not equal")}
		outmat<-Reduce('+',mats)/length(mats)
		matWrite(outmat,file=outnames)
	}

	if( is.null(matlist1)==FALSE & length(matlist2) == length(matlist1) ){
		if(operation %ni% c("log2ratio","ratio","difference","mean")) { stop("operations that can be performed on matrix pairs include log2ratio, ratio, difference, and mean")}
		if( is.null(outnames) & is.null(pattern) & is.null(replacement) ) { stop("must use pattern/replacement OR outnames")}

		nummats<-length(matlist1)

		#check if mats are same dimensions


		if( is.null(outnames)==FALSE){ if(length(outnames) != length(matlist1)){ stop("outnames must be of same length as matlist1")} }


		mats1<-lapply(matlist1,read.mat)
		mats2<-lapply(matlist2,read.mat)

		# check matrix have identical dimensions
		ml1<-unlist(lapply(lapply(lapply(mats1,dim),as.character),paste,collapse="-"))
		ml2<-unlist(lapply(lapply(lapply(mats2,dim),as.character),paste,collapse="-"))
		if(sum(ml1==ml2)!=nummats){stop("matrix dimensions are not identical")}

		dump<-lapply(1:nummats,function(x){
			if(operation=="log2ratio"){mat3<-(log2(mats1[[x]]/mats2[[x]]))}
			if(operation=="ratio"){mat3<-(mats1[[x]]/mats2[[x]])}
			if(operation=="difference"){mat3<-(mats1[[x]]-mats2[[x]])}
			if(operation=="mean"){mat3<-((mats1[[x]]+mats2[[x]])/2)}
			matWrite(mat3,file=outnames[x])
		})

	}
	return(outnames)
}
