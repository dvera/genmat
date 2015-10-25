matHeatmap <-
function( mats=NULL , matnames=NULL , sorton=1 , sort.methods="none" , gene.list = NULL , flip = FALSE , flip.range = c("-50,0","0,50") , flip.method = "max" , align.range = c(-200,0) , align.method = "max" , align = FALSE , sort.ranges = NA , sort.breaks = NA , numgroups=3 , crossmat=NULL , matcolors=NULL , groupcolors=NULL , heatmap.lims="5%,95%" , heatmap.colors="black,white" , agg.lims= NA , fragmats=NULL , fragmatnames=NULL , fragrange=c(0,200) , vplot.colors="black,blue,yellow,red" , vplot.lims="0,95%" , forcescore=TRUE , cores="max" ){


	# TO DO
	# ##############
	# change cores="max" with requestCores()

	# check sorting for expected sorting methods or existing file
	# check that mats is a character vector and all files exist
	# check that all mats have identical row names and nrow
	# check that all colors are R colors
	# make sure align range is within matrix boundaries


	library(gtools)
	library(tools)
	library(parallel)

	pwd<-getwd()

	if(cores=="max"){cores=detectCores()-1}


	if(is.null(mats) & is.null(fragmats)){stop("nothing to plot")}

	# check if provided sort methods are valid
	valid.sort.methods=c("file","none","kmeans","chrom","max","min","mean","median","sd","maxloc","minloc","pregrouped")
	if(sum(sort.methods %in% valid.sort.methods) != length(sort.methods)){stop("unrecognized sort.methods. valid methods include none,kmeans,chrom,max,min,mean,median,sd,maxloc,minloc")}

	if(is.null(mats) == FALSE){

		nummats<-length(mats)

		# get window sizes from file extensions and check if they are good
		winsizes=as.numeric(gsub("mat","",file_ext(mats)))
		if(sum(is.na(winsizes))>0){stop("cannot find window sizes for one or more matrices. window size (in bp) should be in the matrix file extension after 'mat' (for example, matrixname.mat10 for a window size of 10 bp")}

		# #### check mat names for uniqueness and length

		# set the names of matrices using file names if they aren't explicitly defined
		if(is.null(matnames)){ matnames<-basename(removeext(mats)) }

		# set unspecified settings to defaults
		# !!!!!!!!!!!!!!! check to see if length of these arguments are equal to nummats or 1
		if(length(heatmap.lims)==1){heatmap.lims<-rep(heatmap.lims,nummats)}
		if(length(agg.lims)==1){agg.lims<-rep(agg.lims,nummats)}
		if(length(heatmap.colors)==1){heatmap.colors<-rep(heatmap.colors,nummats)}

		# set color gradients for heatmaps
		scolramps<-lapply(heatmap.colors, function(x) colorRampPalette(unlist(strsplit(x,","))))

	}



# ############ !!!!!!!!!!!!!!!!! FIX THIS LATER
	#CHECK FRAGMAT PARAMETERS AND DEFINE FRAGMAT SETTINGS

	if(is.null(fragmats) == FALSE){
		numfragmats<-length(fragmats)
		fragwinsizes=as.numeric(gsub("fmat","",file_ext(fragmats)))
		if(sum(is.na(fragwinsizes))>0){stop("cannot find window sizes for one or more matrices")}
		if(is.null(fragmatnames)){ fragmatnames<-basename(removeext(fragmats)) }
		if(length(vplot.lims)==1){vplot.lims<-rep(vplot.lims,numfragmats)}


		if(length(vplot.colors)==1){vplot.colors<-rep(vplot.colors,numfragmats)}
		vcolramps<-lapply(vplot.colors, function(x) colorRampPalette(unlist(strsplit(x,","))) )

		vplotnames1<-vector(mode="character")
		vplotnames2<-vplotnames1


	}
# ############ !!!!!!!!!!!!!!!!! FIX THIS LATER




# ############ !!!!!!!!!!!!!!!!! FIX THIS LATER

	#load sort table if defined, otherwise sort based on mats defined in sorton
#	if(mode(sort.methods)=="file"){
		# check if file exists and if sort.methods is length 1
#		cat("loading sort table\n")
#		sort.table<-read.table(sorting,quote=F,row.names=1,header=T,sep="\t",stringsAsFactors=F)
#		sort.settings<-data.frame("V1"=colnames(sort.table)[1:(ncol(sort.table)-2)],"V2"="presort","V3"="presort")
#		numgroups<-unlist(lapply(sort.table,length)) #THIS PROBABLY WRONG
#		numsorts<-nrow(sort.settings)
#
	#	#make directory for saving files
	#	if(file.exists("heatmaps")==FALSE){system(command="mkdir heatmaps")}
	#	dname<-paste("heatmaps/",basename(removeext(r)),sep="")
	#	#dname<-uniquefilename(dname)
	#	system(paste("mkdir",dname))
	#	write.table(sort.table,file=paste(dname,"/",basename(dname),".sort",sep=""),sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE)
#	} else{
# ############ !!!!!!!!!!!!!!!!! FIX THIS LATER




		# read matrix for sorting and get dimensions
		sortmat<-read.mat(mats[sorton])
		matcols<-ncol(sortmat)
		matrows<-nrow(sortmat)

		# determine how many sort criteria are defined, and match length of vector defining how many groups to create for each sort criterion
		numsorts<-length(sort.methods)
		length(numgroups)<-numsorts
		numgroups[is.na(numgroups)]<-1

		# if sort range not defined, specify to use entire matrix
		sortranges<-sort.ranges
		sortranges[is.na(sortranges)] <- paste( as.numeric((-1*winsizes[sorton])*(matcols-matcols/2)) , as.numeric(winsizes[sorton]*(matcols-matcols/2)) , sep=",")
		if(length(sortranges)==1){ sortranges <- rep(sortranges , numsorts )}

		# tabulate sort criteria
		sort.settings<-data.frame(
			"METHOD"=sort.methods,
			"START"=as.numeric(unlist(lapply(strsplit(sortranges,","),"[",1))),
			"END"=as.numeric(unlist(lapply(strsplit(sortranges,","),"[",2))),
			stringsAsFactors=FALSE
		)


		# if only 1 sort method, repeat method so sorting script runs OK
		#if(nrow(sort.settings)==1 & sort.settings[1,1] != "kmeans" & sort.settings[1,1] != "chrom"){
		#	sort.settings[2,1]<-sort.settings[1,1]
		#}



		#print sort settings
		cat("\n####################\n")
		cat("SORT CRITERIA\n")
		print(sort.settings)
		cat("####################\n\n")

		#convert bp to column #s in sort settings
		sort.settings[,2]<-sort.settings[,2]/winsizes[sorton] + matcols/2 +1
		sort.settings[,3]<-sort.settings[,3]/winsizes[sorton] + matcols/2
		sort.settings[,2][is.nan(sort.settings[,2])]<-1
		sort.settings[,3][is.nan(sort.settings[,3])]<-1

		#if(normalize==TRUE){
		#	cat("normalizing sort matrix\n")
		#	sortmat<-t(simplify2array(mclapply(1:matrows[r],function(x) sortmat[x,]/mean(sortmat[x,],na.rm=TRUE),mc.cores=cores)))*mean(sortmat,na.rm=TRUE)
		#	row.names(sortmat)<-matlist[[r]]
		#}


		# make sort and group tables
		sort.table<-data.frame(matrix(NA,nrow=matrows,ncol=numsorts),stringsAsFactors=F)
		group.table<-sort.table



		if(align){
			align.start<-align.range[1]/winsizes[sorton] + matcols/2 +1
			align.stop<-align.range[2]/winsizes[sorton] + matcols/2 +1
			allsame <- which(lapply(apply(sortmat[,align.start:align.stop],1,unique),length)==1)

			if(align.method=="max"){
				maxcol <- (apply(sortmat[,align.start:align.stop],1,which.max))
				maxcol[which(unlist(lapply(maxcol,length))==0)]<-round( (align.stop-align.start)/2 )
				maxcol<-unlist(maxcol)
				maxcol[which(is.na(maxcol))] <- round( (align.stop-align.start)/2 )
				maxcol[allsame] <- round( (align.stop-align.start)/2 )
				alignShifts <- winsizes[sorton] * (align.start + maxcol - matcols/2 )
			}
			if(align.method=="min"){
				maxcol <- (apply(sortmat[,align.start:align.stop],1,which.min))
				maxcol[which(unlist(lapply(maxcol,length))==0)]<-round( (align.stop-align.start)/2 )
				maxcol<-unlist(maxcol)
				maxcol[which(is.na(maxcol))] <- round( (align.stop-align.start)/2 )
				maxcol[allsame] <- round( (align.stop-align.start)/2 )
				alignShifts <- winsizes[sorton] * (align.start + maxcol - matcols/2 )
			}
		}

		if(align){
			curalignShifts <- alignShifts / winsizes [ sorton ]

			alignedsortmat <- matrix(NA,nrow=matrows,ncol=3*matcols)
			for(i in 1:matrows){
				alignedstart <- matcols - curalignShifts[i]
				alignedsortmat[i,(alignedstart+1):(alignedstart+matcols)] <- sortmat[i,]
			}
			sortmat.rownames<-row.names(sortmat)
			sortmat<-alignedsortmat[,(matcols+1):(2*matcols)]
			row.names(sortmat)<-sortmat.rownames
			rm(alignedsortmat)
		}



		for(s in 1:numsorts){

			cat("sorting genes\n")
			#sortcol<-vector(mode="numeric",length=matrows)

			#if(s==1){numprevgroups=1} else{ numprevgroups<-length(unique(group.table[,s-1])) }
			if(s==1){numprevgroups=1} else{ numprevgroups<-numgroups[s-1] }

			for(p in 1:numprevgroups){

				if(s==1){grouprows=1:matrows} else{grouprows<-which(group.table[,s-1]==p)}

				#carve out rows and regions to sort
				tmpsortmat<-as.matrix(sortmat[grouprows,sort.settings[s,2]:sort.settings[s,3]])
				submatrows<-nrow(tmpsortmat)

				#convert infinite values to NA
				tmpsortmat[is.infinite(tmpsortmat)]<-NA

				#perform calculation for sorting
				if(sort.settings[s,1] == "sd"){
					#cat("calculating standard deviation\n")
					sortvals<-as.numeric(as.character(apply(tmpsortmat,1,sd, na.rm=TRUE)))
				}
				if(sort.settings[s,1] == "median"){
					#cat("calculating medians\n")
					sortvals<-as.numeric(as.character(apply(tmpsortmat,1,median, na.rm=TRUE)))
				}
				if(sort.settings[s,1] == "max"){
					#cat("calculating maxs\n")
					sortvals<-as.numeric(as.character(apply(tmpsortmat,1,max, na.rm=TRUE)))
				}
				if(sort.settings[s,1] == "min"){
					#cat("calculating mins\n")
					sortvals<-as.numeric(as.character(apply(tmpsortmat,1,min, na.rm=TRUE)))
				}
				if(sort.settings[s,1] == "mean"){
					#cat("calculating means\n")
					sortvals<-as.numeric(as.character(apply(tmpsortmat,1,mean, na.rm=TRUE)))
				}
				if(sort.settings[s,1] == "maxloc"){
					#cat("calculating max locations\n")
					sortvals<-as.numeric(as.character(apply(tmpsortmat,1,which.max)))
				}
				if(sort.settings[s,1] == "minloc"){
					#cat("calculating min locations\n")
					sortvals<-as.numeric(as.character(apply(tmpsortmat,1,which.min)))
				}
				if(sort.settings[s,1] == "chrom"){
					#cat("sorting by chromosome\n")
					chromlist<-unlist(lapply(strsplit(row.names(sortmat),";"),"[",2))
					chroms<-unique(chromlist)
					chroms<-mixedsort(chroms)
					chromnums<-1:length(chroms)
					numgroups[s]<-length(chroms)
					sortvals<-match(chromlist,chroms)
				}
				if(sort.settings[s,1] == "kmeans"){
					#cat("kmeans-clustering data (k=",numgroups[s],")\n",sep="")
					tmpsortmat2<-tmpsortmat
					tmpsortmat2[is.na(tmpsortmat2)]<-0
					sortvals<-kmeans(tmpsortmat2,numgroups[s])$cluster
					rm(tmpsortmat2)
				}
				if(sort.settings[s,1] == "none"){
					#cat("not sorting\n",sep="")
					sortvals<-1:nrow(tmpsortmat)
				}
				if(sort.settings[s,1] == "pregrouped"){
					#cat("not sorting\n",sep="")
					sortvals<-tmpsortmat[,1]
					numgroups[s]<-length(unique(sortvals))
				}


				#assign sortvals to sort.table
				sort.table[grouprows,s] <- sortvals

				# assigned pregrouped data to group table
				if(sort.settings[s,1] == "kmeans" | sort.settings[s,1] == "chrom" | sort.settings[s,1] == "pregrouped"){
					group.table[grouprows,s] <- sortvals
				}
				# if not already pregrouped, group sort values into specified # of groups
				if(sort.settings[s,1] != "kmeans" & sort.settings[s,1] != "chrom" & sort.settings[s,1] != "pregrouped"){

					# #### what if numgroups[s]==1 ???????

					# if breaks not defined for a given sorting criteria
					if(is.na(sort.breaks[s])){


						# find number of genes in an equally-split group
						rowspergroup<-ceiling(submatrows/numgroups[s])

						# define breaks based on number of genes per group
						brks<-c(seq(1,submatrows,rowspergroup),submatrows)

						#assign rows into equally-sized groups
						#group.table[grouprows[order(sortvals)],s]<-cut((1:submatrows)[order(sortvals)],brks,labels=F,include.lowest=T)
						group.table[grouprows,s]<-cut((1:submatrows)[order(order(sortvals,na.last=T))],brks,labels=F,include.lowest=T)

					} else{

						cat("grouping genes into",numgroups[s],"groups based on breaks\n")

						# #set numgroups to length(brks)+1
						numgroups[s]<-length(brks)+1

						brks<-unlist(strsplit(sort.breaks[s],","))

						# if percentiles defined for breaks, calculate break values
						if(grepl("%",brks)){

							# calculate quantile from %s in sort.breaks
							brks[grep("%",brks)]<-quantile(sortvals,probs=as.numeric(gsub("%","",brks[grep("%",brks)])))
						}
						brks<-c(-Inf,brks,Inf)

						# group genes by defined breaks
						group.table[grouprows,s] <- cut(sortvals,brks,include.lowest=T,labels=F)
					}
				}
			}
		}

		# create table of colors
		if(is.null(groupcolors)){ groupcolors<-rainbow(max(numgroups)) }
		getcolor <- function(c){groupcolors[c]}
		color.table <- as.data.frame ( lapply ( group.table , getcolor ) )

		# sort table
		heatmap.order <- order(do.call(order,as.data.frame(cbind(group.table,sort.table[,numsorts]))))

		# name columns of tables
		#colnames(group.table)<-paste("GROUPING",1:numsorts,"-",sort.settings[,1],sort.settings[,3],sort.settings[,3],sep=";")
		#colnames(sort.table)<-paste("SORTVALS",1:numsorts,"-",sort.settings[,1],sort.settings[,3],sort.settings[,3],sep=";")
		#colnames(color.table)<-paste("COLORS",1:numsorts,"-",sort.settings[,1],sort.settings[,3],sort.settings[,3],sep=";")

		#make a table to store order, group, and other info
		heatmap.table <- data.frame(cbind(
			row.names(sortmat),
			heatmap.order,
			color.table,
			sort.table,
			group.table

		))
		colnames(heatmap.table)<-c("gene","order","color",paste0("sort",1:numsorts),paste0("group",1:numsorts))


		# ###### !! NEED FIXING AFTER MAT.MAKE IS CHANGED ####################
		#sort.table$symbol<-unlist(lapply(strsplit(row.names(sort.table),";"), "[", 3 ) )
		#sort.table$id<-unlist(lapply(strsplit(row.names(sort.table),";"), "[", 1 ) )
		#sort.table$chromosome<-unlist(lapply(strsplit(row.names(sort.table),";"), "[", 2 ) )


		#make directory for saving files
		dir.create("heatmaps", showWarnings = FALSE)
		dname<-paste("heatmaps/",basename(matnames[sorton]),"_",numgroups[1],"g_",sort.settings[1,1],"_",sort.settings[1,2],"-",sort.settings[1,3],sep="",collapse="")
		#dname<-uniquefilename(dname)
		dir.create(dname, showWarnings = FALSE)
		write.table(heatmap.table,file=paste0(dname,"/sorting.tsv"),sep="\t",row.names=FALSE,quote=FALSE,col.names=TRUE)
	#}


	if(is.null(mats)==FALSE){
		cat("setting groups and colors\n")
		# ##### get order of preloaded sort table

		# get number of rows of sort table if preloaded
		matrows<-nrow(sort.table)

		# order the table of group info


		# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FIX AND ENABLE THIS SORT SORTFILE CAN BE LOADED. COLOR TABLES GETTING SCRAMBLED
		#heatmap.order<-heatmap.table$heatmap.order
		#table.starts<-1:3*(ncol(heatmap.table)-1)/3
		#group.table<-as.matrix(heatmap.table[,table.starts[2]:(table.starts[2]+1)])
		#color.table<-as.matrix(heatmap.table[,table.starts[3]:(table.starts[3]+1)])
		# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		group.table<-as.matrix(group.table[order(heatmap.order),])
		color.table<-as.matrix(color.table[order(heatmap.order),])


		#group.table<-as.matrix(heatmap.table[,table.starts[2]:(table.starts[2]+1)])[order(heatmap.order,na.last=F),]
		#color.table<-as.matrix(heatmap.table[,table.starts[3]:(table.starts[3]+1)])[order(heatmap.order,na.last=F),]
		#heatmap.table<-heatmap.table[order(heatmap.order,na.last=F),]
		# identify number of rows in each group for each sorting
		clusternums<-lapply(1:numsorts,function(h){ unlist(lapply(1:numgroups[h],function(x) length(which(group.table[,h]==x)))) } )
		grouprownumbers<-lapply(1:numgroups[1], function(x) which(group.table[,1]==x) )

		# set legend names and colors
		if(numgroups[1]>1){
			legendnames<-paste(c(1:numgroups[1],"all"),", n= ",c(clusternums,matrows),sep="")

			# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FIX THIS ######################################
			if(sort.settings[1,1]=="chrom"){
				legendnames<-paste(c(chroms,"all"),", n= ",c(clusternums,matrows),sep="")
			}

			legendcolors<-c(groupcolors,"gray50")

		} else{
			legendnames<-paste("all, n =",matrows)
			legendcolors<-"black"
		}
	}

	# #####################################
	# #               _       _
	# #  __   ___ __ | | ___ | |_ ___
	# #  \ \ / / '_ \| |/ _ \| __/ __|
	# #   \ V /| |_) | | (_) | |_\__ \
	# #    \_/ | .__/|_|\___/ \__|___/
	# #        |_|
	# #
	# #####################################

	#process each frag matrix
	if(is.null(fragmats) == FALSE){
		cat("starting vplots\n")
		for(v in 1:numfragmats){



			fragmat<-read.fmat(fragmats[v])

			fmatcols<-ncol(fragmat)
			fmatrows<-nrow(fragmat)

			if(align){
				curalignShifts <- alignShifts / fragwinsizes[v]
				alignedmat <- matrix(NA,nrow=fmatrows,ncol=3*fmatcols)
				for(i in 1:fmatrows){
					alignedstart <- fmatcols - curalignShifts[i]
					alignedmat[i,(alignedstart+1):(alignedstart+fmatcols)] <- fragmat[i,]
				}
				fragmat.rownames<-row.names(fragmat)
				fragmat<-alignedmat[,(fmatcols+1):(2*fmatcols)]
				row.names(fragmat)<-fragmat.rownames
				rm(alignedmat)
			}

			x<-((1:fmatcols)-(fmatcols/2))*fragwinsizes[v]

			fragsizes<-lapply(1:fmatcols,function(x) as.numeric(na.omit(as.numeric(unlist(strsplit(paste(fragmat[,x],collapse=","),","))))))

			fragsizes<-lapply(1:fmatcols,function(x) fragsizes[[x]][which(fragsizes[[x]] >= fragrange[1] & fragsizes[[x]] <= fragrange[2])])

			vmat<-matrix(0,ncol=fmatcols,nrow=fragrange[2]-fragrange[1])

			for(h in 1:fmatcols){
				vmat[,h]<-hist(fragsizes[[h]],breaks=fragrange[2]:fragrange[1],plot=F)$counts
			}

			#set vmat score limits
			vbot<-strsplit(vplot.lims[v],",")[[1]][1]
			vtop<-strsplit(vplot.lims[v],",")[[1]][2]
			if(grepl("%",vbot)){ vbot<-quantile(vmat , probs=as.numeric(gsub("%","",vbot)) /100,na.rm=T) }
			if(grepl("%",vtop)){ vtop<-quantile(vmat , probs=as.numeric(gsub("%","",vtop)) /100,na.rm=T) }
			vbot<-as.numeric(vbot)
			vtop<-as.numeric(vtop)
			if(vtop==vbot){vtop=vtop+1}

			# if(rpm==TRUE){
			# 	vscalar<-1000000/as.numeric(gsub("#","",readLines(pipe(paste("head -n 1",fragmats[v])))))
			# 	if(is.na(vscalar)==FALSE){vmat<-vmat*vscalar
			# 	} else{cat("WARNING: FRAGMENT COUNT NOT FOUND IN FRAGMAT, NOT NORMALIZING\n")}
			# }

			#vmat<-vmat[nrow(vmat):1,]

			cat("drawing vplot\n")

			brks=seq(vbot,vtop,by=(vtop-vbot)/100)



			# draw heatmap of all
			vplotname<-paste(pwd,"/",dname,"/","vplot_all_",fragmatnames[v],".png",sep="")
			png(width=1000,height=1000,file=vplotname)



			layout(matrix(c(4,2,1,3),nrow=2),heights=c(1,4),widths=c(4,1))
			#layout(matrix(1:2,nrow=2),heights=c(1,4))
			par(oma=c(2,2,2,2))
			par(mar=c(3,3,0,0))
			image(matrix(brks),col=vcolramps[[v]](100),breaks=brks,axes=F)
			axis(side=1,at=c(0,1),labels=c(brks[1],brks[length(brks)]))
			mtext("bin counts",side=1,line=2)

			par(mar=c(3,3,0,0))
			image(t(vmat),breaks=brks,col=vcolramps[[v]](100),xaxt='n',yaxt='n')
			axis(side=1,at=c(0,0.5,1),labels=c(min(x),0,max(x)))
			axis(side=2,at=(0:4)/4,labels=fragrange[1]+0:4*((fragrange[2]-fragrange[1])/4))
			mtext("Distance from feature",side=1,line=2)
			mtext("Fragment size (bp)",side=2,line=2)

			dx<-apply(vmat,2,mean,na.rm=T)
			dy<-apply(vmat,1,mean,na.rm=T)
			par(mar=c(3,0,0,0))
			plot(dy,1:(fragrange[2]-fragrange[1]),type="l",yaxs='i',axes=F,lwd=4)
			par(mar=c(0,3,1,0))
			plot(1:fmatcols,dx,type="l",xaxs='i',axes=F,lwd=4,main=fragmatnames[v])


			dev.off()

			if(is.null(mats)==FALSE){

				fragmat<-fragmat[order(heatmap.order),]

				for(q in 1:numgroups[1]){

					vplotname<-paste(pwd,"/",dname,"/","vplot_group",q,"_",fragmatnames[v],"_pergene.png",sep="")

					groupmat<-fragmat[grouprownumbers[[q]],]

					fragsizes<-lapply(1:fmatcols,function(x) as.numeric(na.omit(as.numeric(unlist(strsplit(paste(groupmat[,x],collapse=","),","))))))

					fragsizes<-lapply(1:fmatcols,function(x) fragsizes[[x]][which(fragsizes[[x]] >= fragrange[1] & fragsizes[[x]] <= fragrange[2])])
					vmat<-matrix(0,ncol=fmatcols,nrow=fragrange[2]-fragrange[1])

					for(h in 1:fmatcols){
						vmat[,h]<-hist(fragsizes[[h]],breaks=fragrange[2]:fragrange[1],plot=F)$counts
					}

					# if(rpm==TRUE){
					# 	vscalar<-1000000/as.numeric(gsub("#","",readLines(pipe(paste("head -n 1",fragmats[v])))))
					# 	if(is.na(vscalar)==FALSE){vmat<-vmat*vscalar
					# 	} else{cat("WARNING: FRAGMENT COUNT NOT FOUND IN FRAGMAT, NOT NORMALIZING\n")}
					# }

					vmat<-vmat[nrow(vmat):1,]

					pergenescale<-nrow(fragmat)/nrow(groupmat)

					brks=seq(vbot,vtop,by=(vtop-vbot)/100)

					png(width=1000,height=1000,file=vplotname)
					heatmap.2(vmat,Colv=NA,Rowv=NA,dendrogram="none",trace="none",breaks=brks,labRow=NA,labCol=NA,col=vcolramps[[v]](100),main=fragmatnames[v],cex.main=5)
					dev.off()
					png(width=1000,height=1000,file=vplotname)
					heatmap.2(vmat*pergenescale,Colv=NA,Rowv=NA,dendrogram="none",trace="none",breaks=brks,labRow=NA,labCol=NA,col=vcolramps[[v]](100),main=fragmatnames[v],cex.main=5)
					dev.off()
				}

			}

		}

		#montagecols<-numgroups[1]+1
		#system(paste("montage -geometry +2+2 -tile ",montagecols,"x",numfmats," ",paste(vplotnames1,collapse=" ", sep=" ")," ",pwd,"/",dname,"/vplot_montage_abs.png",sep=""))
		#system(paste("montage -geometry +2+2 -tile ",montagecols,"x",numfmats," ",paste(vplotnames2,collapse=" ", sep=" ")," ",pwd,"/",dname,"/vplot_montage_per.png",sep=""))
	}


	# ##############################################################################
	# TEMPORARY IMPLEMENTATION OF CROSS-MATRIX GROUP AVERAGES
	#grouplist=lapply(1:nummats, function(x){lapply(1:numgroups[1],function(x){})})
	#q=1
	# ##############################################################################

	#process each matrix
	hmnames<-paste(pwd,"/",dname,"/","heatmap_",basename(matnames),".png",sep="")

	#for(l in 1:nummats){
	for(l in 1:nummats){
		cat("starting mats\n")
		# load matrix and get dimensions
		curmat<-read.mat(mats[l])
		curmat[is.infinite(curmat)]<-NA

		matcols<-ncol(curmat)
		matrows<-nrow(curmat)

		if(align){
			cat("printing bp\n")
			curalignShifts <- alignShifts / winsizes[l]
			cat("printing numwins\n")
			alignedmat <- matrix(NA,nrow=matrows,ncol=3*matcols)
			for(i in 1:matrows){
				alignedstart <- matcols - curalignShifts[i]
				alignedmat[i,(alignedstart+1):(alignedstart+matcols)] <- curmat[i,]
			}
			curmat.rownames<-row.names(curmat)
			curmat<-as.matrix(alignedmat[,(matcols+1):(2*matcols)])
			row.names(curmat)<-curmat.rownames
			rm(alignedmat)
		}

		# sort matrix
		curmat<-as.matrix(curmat[order(heatmap.order),])

		# define x axis based on winsize and ncol assuming non-meta
		x<-((1:matcols)-(matcols/2))*winsizes[l]

		# set score limits for heatmap colors
		heatmap.mins<-unlist(lapply(strsplit(heatmap.lims,","),"[",1))
		heatmap.maxs<-unlist(lapply(strsplit(heatmap.lims,","),"[",2))
		if(grepl("%",heatmap.mins[l])){ heatmap.mins[l]<-quantile(curmat,probs=as.numeric(gsub("%","",heatmap.mins[l])) /100,na.rm=T) }
		if(grepl("%",heatmap.maxs[l])){ heatmap.maxs[l]<-quantile(curmat,probs=as.numeric(gsub("%","",heatmap.maxs[l])) /100,na.rm=T) }


		curmin<-as.numeric(heatmap.mins[l])
		curmax<-as.numeric(heatmap.maxs[l])
		if(curmin==curmax){curmax<-max(curmat,na.rm=TRUE)}
		if(curmin==curmax){curmax<-curmin+1}


		#calculate aggregate profiles


		groupmeans<-lapply(1:numgroups[1], function(y){
			colMeans(as.matrix(curmat[grouprownumbers[[y]],]),na.rm=TRUE)
		})
		groupmeans[[numgroups[1]+1]]<-colMeans(curmat,na.rm=TRUE)
		groupmeanmaxs<-unlist(lapply(groupmeans,max,na.rm=TRUE))
		groupmeanmins<-unlist(lapply(groupmeans,min,na.rm=TRUE))

		if(is.na(agg.lims) | is.null(agg.lims)){
			agg.min<-min(groupmeanmins)
			agg.max<-max(groupmeanmaxs)
		} else{
			agg.min<-as.numeric(unlist(lapply(strsplit(agg.lims[l],","),"[",1)))
			agg.max<-as.numeric(unlist(lapply(strsplit(agg.lims[l],","),"[",2)))
		}

		# ##############################################################################
		# TEMPORARY IMPLEMENTATION OF CROSS-MATRIX GROUP AVERAGES
		#if(is.null(crossmat)==FALSE & l %in% crossmat == TRUE){
		#	grouplist[[q]]=groupmeans
		#	q=q+1
		#}
		# ##############################################################################

		#make NAs 0 for aesthetics
		if(forcescore==TRUE){
			curmat[is.na(curmat)]<-0
			curmat[is.nan(curmat)]<-0
			curmat[is.infinite(curmat)]<-0
		}

		png(file=hmnames[l],height=5000,width=1000)

		# define color breaks
		brks<-c(seq(curmin,curmax,by=(curmax-curmin)/100),curmax)


		# copy matrix into new object
		mat<-curmat

		# create layout matrix
		m<-matrix(1:3,nrow=3)

		# set out-of-boundary scores to boundaries
		mat[which(mat>curmax)]<-curmax
		mat[which(mat<curmin)]<-curmin

		# create image layout and set margins
		layout(m,heights=c(1,20,5),widths=1)
		par(oma=c(10,0,10,0))
		par(mar=c(10,15,5,15))

		# draw color scale
		image(matrix(brks),col=scolramps[[l]](101),breaks=brks,axes=F)
		axis(side=1,at=c(0,1),labels=c(curmin,curmax),cex.axis=8,lwd=10,padj=2,line=0,tcl=-3)

		# title image
		mtext(matnames[l],side=3,cex=8,outer=T)

		# set margins and draw heatmap
		par(mar=c(5,15,5,15),xpd=TRUE)
		image(t(mat),breaks=brks,col=scolramps[[l]](101),axes=FALSE)


		# create group color bars next to heatmap
		for(g in 1:numsorts){
			for(i in 1:numgroups[g]){
				linevec <- group.table[,g]
				linevec[which(group.table[,g] != i)] <- NA
				linecoords <- (1:matrows)/matrows
				linecoords[which(is.na(linevec))] <- NA
				lines(rep(-0.01-g*0.05,matrows),linecoords[matrows:1],lwd=40,col=groupcolors[i],lend=1)
			}
		}
		axis(side=1,at=c(0,0.5,1),labels=F,lwd=10,padj=1,line=1,tcl=-3)

		# draw average plot
		par(mar=c(20,15,0,15))
		plot(0,type="n",xlim=c(min(x)-winsizes[l],max(x)),ylim=c(agg.min,agg.max),xaxs="i",axes=FALSE)
		#for(i in 1:(numgroups[1]+1)){
		for(i in 1:(numgroups[1])){
			lines(x,groupmeans[[i]],lwd=10,col=legendcolors[i])
		}
		axis(side=1,at=c(min(x)-winsizes[l],0,max(x)),labels=c(min(x)-winsizes[l],0,max(x)),cex.axis=8,lwd=10,padj=2,line=0,tcl=-3)
		axis(side=2,cex.axis=8,lwd=10,tcl=-3,padj=-1)
		mtext("Distance from feature (bp)",side=1,cex=5,outer=T,line=2)
		dev.off()


		#par(oldpar)
		# make pdf aggregate profile
		# pdf(file=paste0(removeext(hmnames[l]),".pdf"))
		# plot(0,type="n",xlim=c(min(x)-winsizes[l],max(x)),ylim=c(agg.min,agg.max),xaxs="i",axes=FALSE)
		# #for(i in 1:(numgroups[1]+1)){
		# for(i in 1:(numgroups[1])){
		# 	lines(x,groupmeans[[i]],lwd=10,col=legendcolors[i])
		# }
		# axis(side=1,at=c(min(x)-winsizes[l],0,max(x)),labels=c(min(x)-winsizes[l],0,max(x)),cex.axis=8,lwd=10,padj=2,line=0,tcl=-3)
		# axis(side=2,cex.axis=8,lwd=10,tcl=-3,padj=-1)
		# mtext("Distance from feature (bp)",side=1,cex=5,outer=T,line=2)
		# dev.off()


	}
	#},mc.cores=cores)

	# make montage
	#if(nummats < 10){
	#	system(paste("montage -geometry +2+2 -tile ",nummats,"x1 ",paste(hmnames,collapse=" ", sep=" ")," ",pwd,"/",dname,"/heatmap_montage.png",sep=""))
	#}

	# ##############################################################################
	# TEMPORARY IMPLEMENTATION OF CROSS-MATRIX GROUP AVERAGES
	# if(is.null(crossmat)==FALSE){
	# 	x<-((1:matcols)-(matcols/2))*winsizes[1]
	# 	pdf(file=paste(pwd,"/",dname,"/","group-averages",".pdf",sep=""))
	# 	for(i in 1:(numgroups[1])){
	# 		cat("i=",i,"\n")
	# 		groupscores=lapply(lapply(grouplist,"[",i),unlist)
	# 		allscores=unlist(groupscores)
	# 		print(allscores)
	# 		#plot(0,type="n",xlim=c(min(x),max(x)),ylim=c(min(allscores),max(allscores)),xlab="Distance from TSS",ylab="Average Score",main=paste("Group",i))
	# 		plot(0,type="n",xlim=c(min(x),max(x)),ylim=c(1,2.5),xlab="Distance from TSS",ylab="Average Score",main=paste("Group",i))
	# 		for(j in 1:length(crossmat)){
	# 			cat("j=",j,"\n")
	# 			lines(x,groupscores[[j]],lwd=3,col=matcolors[j])
	# 		}

	# 	}
	# 	dev.off()
	# }
	# ##############################################################################
}
