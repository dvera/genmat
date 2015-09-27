mat.tilegrid <-
function (mat1 , mat2 , scoremats , mat1name = "x" , mat2name = "y" , LOG10=TRUE , abx = 0 , aby = 0 , legendnames = basename(removeext(scoremats)) , plotcolors=rainbow(length(scoremats)) , region = c(-1000,0) , tiles = 3 , ylims = c(-1,1) , breaks=NULL ){


	# load matrices
	s<-read.mat(mat1)
	r<-read.mat(mat2)

	# set quantile breakpoints given # of tiles
	p=c(0,(1:(tiles-1))*1/tiles,1)

	# find scores at break points
	sq<-quantile(s[,1],probs=p,na.rm=T)
	rq<-quantile(r[,1],probs=p,na.rm=T)

	# group scores using break points
	cs<-cut(s[,1],sq,labels=F,include.lowest=T)
	cr<-cut(r[,1],rq,labels=F,include.lowest=T)

	# create table defining group of each row in both matrices
	etab<-data.frame(row.names(s),cs,cr,stringsAsFactors=FALSE)
	colnames(etab)<-c("gene",mat1name,mat1name)

	# save table of groupings
	write.tsv(etab,file=paste0(mat1name,"Vs",mat2name,".tsv"),colnames=T)

	pdf(file=paste0(mat1name,"Vs",mat2name,".pdf"))

	# read score matrices
	sm2<-lapply(scoremats,read.mat)

	# count number of score matrices
	numscoremats<-length(scoremats)

	if(LOG10){

		# take log10 of mat1 and mat2
		ls=log10(s[,1]+1)
		lr=log10(r[,1]+1)

		# set axis labels
		xlabel=paste("log10 (",mat1name,")")
		ylabel=paste("log10 (",mat2name,")")

		# find new breakpoints of mat1 and mat2 after log10
		lsq<-quantile(ls,probs=p,na.rm=T)
		lrq<-quantile(lr,probs=p,na.rm=T)

	} else{

		# use score in first column
		ls=s[,1]
		lr=r[,1]

		# set axis labels
		xlabel=mat1name
		ylabel=mat2name

		# use old break points
		lsq<-sq
		lrq<-rq

	}

	# plot scatterplot of mat1 vs mat2
	plot(ls,lr,pch=16,xlab=xlabel,ylab=ylabel)
	abline(v=lsq[2:(tiles)],h=lrq[2:(tiles)],col="grey50")

	# count number of columns in scoremats
	matcols<-ncol(sm2[[1]])

	# set x axis tick locations
	xa<- ((0:4*(matcols/4)))

	# set y axis tick locations
	ya<-c(ylims[1],(ylims[2]+ylims[1])/2,ylims[2])

	# make pairwise tile indices
	eg<-expand.grid(1:tiles,1:tiles)

	# create layout matrix
	ml<-cbind(rbind(t(matrix(1:(tiles^2) , nrow=tiles)),(tiles^2+1):(tiles^2+tiles)),(tiles^2+tiles+1):((tiles+1)^2))
	ml<-ml[(tiles+1):1,]

	# find indices of bottom and left tiles
	bottoms<-ml[tiles+1,]
	lefts<-ml[,1]

	# set graphical parameters
	layout(ml)
	par(oma=c(5,5,5,5),mar=c(0.2,0.2,0.2,0.2))

	for(i in 1:(tiles^2)){
		cat(i,"\n")
		plot(0,type="n",xlim=c(0,matcols),ylim=ylims,yaxt='n',xaxt='n')
		if(i %in% lefts ){
			axis(2,at=ya,labels=ya)
			mtext(which(i==rev(lefts)),side=2,line=2,cex=2)
		}
		if(i %in% bottoms){
			axis(1,at=xa,labels=xa)
			mtext(i,side=1,line=4,cex=2)
		}
		abline(v=abx,h=aby,col="grey70")
		for(j in 1:numscoremats){
			cat("j=",j,"\n")
			lines(1:matcols,colMeans(sm2[[j]][which(cs==eg[i,1] & cr==eg[i,2]),],na.rm=T),col=plotcolors[j])
		}
		text(xa[1],ya[3],paste0("n=",length(which(cs==eg[i,1] & cr==eg[i,2]))),pos=4)
	}

	for(i in 1:tiles){
		plot(0,type="n",xlim=c(0,matcols),ylim=ylims,yaxt='n',xaxt='n')
		rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "grey90")
		abline(v=abx,h=aby,col="grey70")
		if(i == 1){ axis(2,at=ya,labels=ya) }
		for(j in 1:numscoremats){
			lines(1:matcols,colMeans(sm2[[j]][which(cr==i),],na.rm=T),col=plotcolors[j])
		}
		text(xa[1],ya[3],paste0("n=",length(which(cr==i))),pos=4)
	}
	for(i in 1:tiles){
		plot(0,type="n",xlim=c(0,matcols),ylim=ylims,yaxt='n',xaxt='n')
		rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "grey90")
		abline(v=abx,h=aby,col="grey70")
		if(i == 1){ axis(1,at=xa,labels=xa) }
		for(j in 1:numscoremats){
			lines(1:matcols,colMeans(sm2[[j]][which(cs==i),],na.rm=T),col=plotcolors[j])
		}
		text(xa[1],ya[3],paste0("n=",length(which(cs==i))),pos=4)
	}
	plot(0,type="n",xlim=c(0,matcols),ylim=ylims,yaxt='n',xaxt='n')
	rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "grey90")
	abline(v=abx,h=aby,col="grey70")
	for(j in 1:numscoremats){
			lines(1:matcols,colMeans(sm2[[j]],na.rm=T),col=plotcolors[j])
	}
	text(xa[1],ya[3],paste0("n=",nrow(sm2[[1]])),pos=4)
	legend("topright",legend=legendnames , col=plotcolors , lwd=3 , bty="n")


	if(length(scoremats)==2){
		layout(ml)
		par(oma=c(5,5,5,5),mar=c(0.2,0.2,0.2,0.2))

		all<-as.vector(unlist(sm2))
		sylims<-quantile(all,probs=c(0.0001,0.9999),na.rm=T)
		for(i in 1:(tiles^2)){
			xscores <- as.vector(sm2[[1]][which(cs==eg[i,1] & cr==eg[i,2]),])
			yscores <- as.vector(sm2[[2]][which(cs==eg[i,1] & cr==eg[i,2]),])
			plot(xscores, yscores, col=rgb(0,0,0,20,maxColorValue=255),xlim=sylims,ylim=sylims,pch=20,yaxt='n',xaxt='n')
			if(i %in% lefts ){
				axis(2)
				mtext(which(i==rev(lefts)),side=2,line=2,cex=2)
			}
			if(i %in% bottoms){
				axis(1)
				mtext(i,side=1,line=4,cex=2)
			}
			#abline(v=abx,h=aby,col="grey70")

			text(sylims[2],sylims[1],paste0("n=",length(which(cs==eg[i,1] & cr==eg[i,2]))),pos=4)
			text(sylims[1],sylims[2],cor(xscores,yscores,use="complete.obs"),pos=4)
		}

		for(i in 1:tiles){
			plot(0,type="n",xlim=c(0,matcols),ylim=ylims,yaxt='n',xaxt='n')
			rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "grey90")
			abline(v=abx,h=aby,col="grey70")
			if(i == 1){ axis(2,at=ya,labels=ya) }
			for(j in 1:numscoremats){
				lines(1:matcols,colMeans(sm2[[j]][which(cr==i),],na.rm=T),col=plotcolors[j])
			}
			text(xa[1],ya[3],paste0("n=",length(which(cr==i))),pos=4)
		}
		for(i in 1:tiles){
			plot(0,type="n",xlim=c(0,matcols),ylim=ylims,yaxt='n',xaxt='n')
			rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "grey90")
			abline(v=abx,h=aby,col="grey70")
			if(i == 1){ axis(1,at=xa,labels=xa) }
			for(j in 1:numscoremats){
				lines(1:matcols,colMeans(sm2[[j]][which(cs==i),],na.rm=T),col=plotcolors[j])
			}
			text(xa[1],ya[3],paste0("n=",length(which(cs==i))),pos=4)
		}
		plot(0,type="n",xlim=c(0,matcols),ylim=ylims,yaxt='n',xaxt='n')
		rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "grey90")
		abline(v=abx,h=aby,col="grey70")
		for(j in 1:numscoremats){
				lines(1:matcols,colMeans(sm2[[j]],na.rm=T),col=plotcolors[j])
		}
		text(xa[1],ya[3],paste0("n=",nrow(sm2[[1]])),pos=4)
		legend("topright",legend=legendnames , col=plotcolors , lwd=3 , bty="n")

	}




	dev.off()

}
