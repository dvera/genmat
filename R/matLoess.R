matLoess <-
function( mats, lspan=0.05 ){
	library(parallel)
	matnames<-basename(removeext(mats))
	print(matnames)
	for(i in 1:length(mats)){
		cat("reading mat\n")
		m<-read.mat(mats[i])
		y<-1:ncol(m)
		ml<-data.matrix(t(as.data.frame(mclapply(1:nrow(m),function(x){  m[x,-which(is.na(m[x,]))]<-loess(m[x,]~y,span=lspan)$fitted ; m[x,]   },mc.cores=detectCores()))))
		cat("smoothing\n")
		row.names(ml)<-row.names(m)
		cat("saving mat\n")
		write.mat(ml,file=paste(matnames[i],"_loess",gsub("\\.","-",as.character(lspan)),".mat10",sep=""))
	}
}
