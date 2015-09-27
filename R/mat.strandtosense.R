mat.strandtosense <-
function( topstrandmats,botstrandmats,bed){
	curbed<-read.tsv(bed)
	negrows<-which(curbed[,6]=="-")
	for(i in 1:length(topstrandmats)){
		cat("processing",topstrandmats[i],"\n")
		curplus<-read.mat(topstrandmats[i])
		curminus<-read.mat(botstrandmats[i])
		newplus<-curplus
		newminus<-curminus
		newplus[negrows,]<-curminus[negrows,]
		newminus[negrows,]<-curplus[negrows,]
		write.mat(newplus,file=gsub("\\.mat","_sense.mat",topstrandmats[i]))
		write.mat(newminus,file=gsub("\\.mat","_antisense.mat",botstrandmats[i]))
	}
}
