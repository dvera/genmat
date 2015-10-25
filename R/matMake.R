#' Create a matrix of data aligned at a set of features.
#'
#' \code{mat.make} creates a matrix of scores aligned at a set of genomic intervals, where each row is an interval and each column is a window surrounding each interval.
#'
#' @param scorefiles A character vector of paths to data files that are used to fill the matrix. Can be bed, bedGraph, wiggle, bigWig. Can also include URLs and/or gzipped forms of these files.
#' @param features A character vector of paths to genomic intervals corresponding to rows in the matrix that are used to align values in 'scorefiles'.
#' @param featurenames A character vector equal in length to 'features' that is used as a suffix in the output file names. If NULL, uses the file names of 'features'.
#' @param regionsize A positive integer specifying what size (bp) interval around each feature is used to generate the matrix. Must be divisible by 'windowsize'.
#' @param windowsize A positive integer specifying what size (bp) windows are used to generate the matrix.
#' @param strand Boolean indicating if the orientation of the scores around minus-stranded features should be reversed in the matrix.
#' @param featurecenter Boolean indicating if the matrix should be aligned to the midpoint of the intervals in 'features'. If FALSE, features will be aligned to the coordinate specified by 'start' (or 'stop' for minus-stranded features if 'strand' is TRUE).
#' @param start A positive integer specifying the column to be used in 'features' as the 5' end of features. Can be set to '7' for aligning at the start codons of genes in bed12 format, '3' for aligning at the end of genes, or '8' for aligning at stop codons.
#' @param stop A positive integer specifying the column to be used in 'features' as the 3' end of features.
#' @param prunefeaturesto A string specifying a bed file to filter intervals in 'features' by. Only invervals in 'features' that overlap with intervals in 'prunefeaturesto' will be used to construct the matrix. Useful if data in 'scorefiles' only exists for a proportion of the genome and wish to limit features used to construct the matrix to those regions.
#' @param narrowpeak Boolean indicating if the matrix are aligned to summits of peaks in a narrowpeak file that is present in 'features'. Only used if summits are specified in the narrowpeak file and if 'featurecenter' is TRUE.
#' @param maskbed A string specifying a bed file to filter values in the matrix by. If TRUE, only windows in the matrix that overlap with intervals in 'maskbed' will be assigned scores, and all other windows will contain 'NA'. Useful for preventing the FALSE assignement of zeroes in regions that have no value in 'scorefiles', such as in data from sequence-capture experiments.
#' @param bgfiller A value to be used for windows in the matrix that have no score in bedGraph files specified in 'scorefiles'. Can be a number or NA.
#' @param prunescores Boolean indicating if scores are intersected with the locations of windows in the matrix prior to generating the matrix. Does not affect or change the resulting matrix but does speed up the function for very large files in 'scoresfiles'.
#' @param rpm Boolean indicating if interval densities calculated from bed files specified in 'scorefiles' are normalized to reads-per-million (RPM). If TRUE, bed interval counts in each window is multiplied by 1000000/(# of intervals in bed file). Useful if using reads with varying amounts in 'scorefiles' and a comparison of resulting matrices is desired.
#' @param closest A string specifying a bed file to label features in the matrix. This will not affect the values in matrices but will change the names for each feature in the matrix based on the closest interval in 'closest'. Useful for assigning gene names to ChIP peaks used as 'features'.
#' @param cores A positive integer specifying how many bams to process simultaneously.
#' @param meta Boolean indicating if 'metafeatures' are used to generate the matrix. If FALSE, matrices are aligned to single coordinates, and all windows in the matrix represent a fixed-sized window surrounding the coordinates. If TRUE, matrices are aligned to both the 'start' and 'stop' coordinates, and the distances and windows between the start and stop are scaled to a fixed value specified in 'regionsize'. Useful for generating a matrix of values within and around varying-size genes.
#' @param metaflank A positive integeger specifying the distance from metafeature boundaries used to generate the matrix. Must be divisible by windowsize.
#' @param suffix A string that is appended to the output matrix names.
#' @param scoremat Boolean indicating if score matrices are generated relative to features.
#' @param fragmats Boolean indicating if a special type of matrix, a fragment-size matrix, is generated using the sizes of intervals for bed files in 'scorefiles'. Such matrices can be used to generate 'v-plots'.


matMake <-
function( scorefiles , features , featurenames = NULL , regionsize=2000 , windowsize=10, strand=FALSE , featurecenter=TRUE , start=2 , stop=3 , prunefeaturesto=NULL , narrowpeak=FALSE , maskbed=NULL , bgfiller=0 , prunescores=FALSE , rpm=TRUE , closest=NULL , cores="max" , meta=FALSE , metaflank=1000 , suffix=NULL , scoremat=TRUE , fragmats=FALSE ){

        # TO DO
        # #######
        # add bed info to mat rownames
        # move tmp files to a tmp directory

  options(scipen=99999)

  bed.split <- function( bed,regionsize,windowsize,covbedname ){
  	numwindows<-(regionsize/windowsize)

  	bedname<-basename(removeext(bed))
  	curbed<-read.tsv(bed)

  	bedrows<-nrow(curbed)
  	bedcols<-ncol(curbed)
  	#check parameters
  	if(ceiling(numwindows)!=floor(numwindows)){stop("regionsize is not a multiple of windowsize")}
  	if(ceiling(regionsize)!=floor(regionsize)){stop("regionsize must be an even number")}

  	#make covbed
  	flanks<- (  0:(numwindows-1)  ) * windowsize
  	winstarts<-as.numeric( unlist( lapply( curbed[,2] , function(x) x + flanks ) ) )
  	covbed<-data.frame(
  		"V1"=rep(curbed[,1],each=numwindows),
  		"V2"=format(winstarts,scientific=F,trim=T),
  		"V3"=format(winstarts+windowsize,scientific=F,trim=T),
  		"V4"=1:(bedrows*numwindows)
  	)
  	write.tsv(covbed,file=covbedname)
  	system(paste("sort -k1,1 -k2,2n -k3,3n",covbedname,"-o",covbedname))
  	covbedorder.call <- pipe(paste("cut -f 4",covbedname),open="r")
  	covbedorder<-as.numeric(readLines(covbedorder.call))
  	close(covbedorder.call)
  	return(covbedorder)
  }

  bed.split.meta <- function( bed, metasize, windowsize, flank, covbedname, start=2, stop=3 ){

  	options(scipen=99999)

  	curbed<-read.tsv(bed)
  	bedrows<-nrow(curbed)
  	bedcols<-ncol(curbed)

  	numwindows<-metasize/windowsize
  	numflankwindows<-flank/windowsize
  	leftwinstarts<-0:(numflankwindows-1) * windowsize - flank
  	leftwinends<-1:numflankwindows * windowsize - flank
  	rightwinstarts<-0:(numflankwindows-1) * windowsize
  	rightwinends<-1:numflankwindows * windowsize
  	genesizes<-curbed[,stop]-curbed[,start]
  	sizecos<-genesizes/metasize
  	genewinstarts<-mclapply(1:bedrows, function(x) curbed[x,start] + 0:(numwindows-1) * windowsize * sizecos[x], mc.cores=detectCores())
  	genewinends<-mclapply(1:bedrows, function(x) curbed[x,start] + 1:numwindows * windowsize * sizecos[x], mc.cores=detectCores())
  	genewinstarts<-mclapply(genewinstarts,round,mc.cores=detectCores())
  	genewinends<-mclapply(genewinends,round,mc.cores=detectCores())

  	covbed<-data.frame(
  		"V1"=rep(curbed[,1],each=numwindows+numflankwindows*2),
  		"V2"=format(unlist(lapply(1:bedrows,function(x){ c(leftwinstarts+curbed[x,start],genewinstarts[[x]],rightwinstarts+curbed[x,stop]) } )),scientific=F,trim=T),
  		"V3"=format(unlist(lapply(1:bedrows,function(x){ c(leftwinends+curbed[x,start], genewinends[[x]],rightwinends+curbed[x,stop]) } )),scientific=F,trim=T),
  		"V4"=1:(bedrows*(numwindows+numflankwindows*2))
  	)

  	write.tsv(covbed,file=covbedname)
  	system(paste("sort -k1,1 -k2,2n -k3,3n",covbedname,"-o",covbedname))
    covbedorder.call <- pipe(paste("cut -f 4",covbedname),open="r")
    covbedorder<-as.numeric(readLines(covbedorder.call))
    close(covbedorder.call)
  	return(covbedorder)


  }



  numscores <- length(scorefiles)
  numfeats <- length(features)
  if(is.null(featurenames)){
  	featurenames<-basename(removeext(features))
  } else{
  	if(length(featurenames) != numfeats ){stop("length of feature names should equal length of feature files")}
  }


	if(cores=="max"){cores=detectCores()-1}
	cores2<-floor(cores/numscores)
	if(cores2<1){cores2=1}
	if(cores2>numfeats){cores2 <- numfeats }

	numwindows<-regionsize/windowsize
	numfeats<-length(features)
	numscores<-length(scorefiles)

	#check parameters
	if(ceiling(numwindows)!=floor(numwindows)){stop("regionsize is not a multiple of windowsize")}
	if(ceiling(regionsize)!=floor(regionsize)){stop("regionsize must be an even number")}

	#make directory for saving files

	dnames<-paste0(featurenames,"_mat",windowsize)
	scorenames <- basename(removeext(scorefiles))
	snames<-lapply(1:numfeats, function(x) paste0(dnames[x],"/",scorenames,"_",featurenames[x],suffix,".mat",windowsize))
	names(snames) <- dnames

	#process features
	mclapply(1:numfeats,function(i){

  	featfile<-features[i]
		featname<-featurenames[i]
		dname <- dnames[i]
		system(paste("mkdir -p",dname))

		if(is.null(closest) == FALSE){
			featfile<-bed.closest(featfile,closest,strand=strand)
		}

		#recenter narrowpeaks
		if(file_ext(featfile) %in% c("narrowPeak","narrowpeak","np") & featurecenter & narrowpeak){
			cat("narrowPeak file detected\n")

			#check if peak location is in narrowPeak file and recenter features if OK
      headlines.call<-pipe(paste("head",featfile,"| awk '{print $10}'"),open="r")
      headlines<-as.numeric(readLines(headlines.call))
      close(headlines.call)

			if(length(which(headlines %in% (-1)>0))){
				cat("bad peak coordinate detected in narrowPeak file, treating as bed\n")
			} else{
				cat("peak coordinates detected in narrowPeak file, recentering bed features\n")
				outname<-paste(featname,"_recentered.np",sep="")
				system(paste("awk '{$2=$2+$10;$3=$2+1;print}' OFS='\t' ",featfile," | sort -k1,1 -k2,2n > ",outname,sep=""))
				featfile<-outname
			}
		}

		#create window bed for pruning
		if(meta==FALSE){
			cat("calculating matrix boundaries\n")
			featfile<-bed.recenter(featfile,regionsize,center=featurecenter,strand=strand,start=start,stop=stop )
		}

		#prune features to specified file
		if(is.null(prunefeaturesto) == FALSE){
			cat("pruning features\n")
			featfile<-bed.intersect(featfile,prunefeaturesto)
		}

		#copy processed features to output directory for later use
		system(paste("cp ",featfile," ",dname,"/",sep=""))

		#read in bed and get dimenscurions
		cat("reading in bed and getting info\n")
		curbed<-read.tsv(featfile)

		#discard regions beyond chromosome boundaries
    if(meta){
			curbed<-curbed[which(curbed[,2]>metaflank+1),]
		}
    write.tsv(curbed,file=featfile)

		bedcols<-ncol(curbed)
		bedrows<-nrow(curbed)

		# FEATURES SHOULD NOT BE MODIFIED AFTER THIS POINT #

		#find minus-stranded features
		if(bedcols < 6 & strand==TRUE){strand=FALSE;cat("WARNING: RUNNING WITH strand=FALSE BECAUSE NO STRAND COLUMN EXISTS\n")}
		if(strand){negrows<-which(curbed[,6]=="-")}

		#make windows to calculate scores for matrix
		cat("calculating coordinates for matrix windows\n")
		if(meta){
		      covbedname <- paste0(basename(removeext(featfile)),".covbed")
		      covbedorder<-bed.split.meta(featfile,regionsize,windowsize, metaflank, covbedname , start=start, stop=stop )
		      numwindows<-(metaflank*2+regionsize)/windowsize
		} else{
		      covbedname <- paste0(basename(removeext(featfile)),".covbed")
		      covbedorder <- bed.split( featfile , regionsize , windowsize , covbedname )
		}

		#make matrix of regions to (not) mask with NAs
		if(is.null(maskbed) == FALSE){

      cat("finding masking regions\n")
      maskmat.call <- pipe(paste("bedtools coverage -b",maskbed,"-a",covbedname," | sort -k4,4n | cut -f 5"),open="r")
			maskmat<-t(matrix(as.numeric(readLines(maskmat.call)),nrow=numwindows))
      close(maskmat.call)


    	if(strand){
        maskmat[negrows,1:numwindows]<-maskmat[negrows,numwindows:1]
		  }

		}

		#assign row names to matrix
		cat("naming features\n")
		if(bedcols<4){
			geneids<-1:bedrows
		} else{
    	  geneids<-curbed[,4]
		}

    if(length(unique(geneids)) != bedrows){
			geneids<-paste(geneids,1:bedrows,sep="-")
		}


    if(bedcols>5){
			strands<-curbed[,6]
		} else{
      strands<-rep("+",bedrows)
		}


		if(bedcols>12){
			symbs<-curbed[,13]
		} else{
			symbs<-geneids
		}

		matrownames <- paste(curbed[,1],curbed[,2],curbed[,3],geneids,strands,symbs,sep=";")




		#process each score file
		scores<-scorefiles
    outs<-mclapply(1:length(scores), function(j) {

    	scorefile<-scores[j]
			scorename<-scorenames[j]

			#DOWNLOAD, EXTRACT, AND/OR FORMAT CONVERSION
			if(grepl("tp://",scorefile)){
				cat(scorename,": URL detected, attempting to download file to current directory\n")
				download.file(scorefile,basename(scorefile))
				scorefile<-basename(scorefile)
			}
			if(file_ext(scorefile) == "gz"){
				cat(scorename,": extracting with gunzip to current directory\n")
				system(paste("gunzip -c",scorefile,">",scorename ) )
				scorefile<-paste(basename(removeext(scorefile)),sep="")
				scorename<-basename(removeext(scorefile))
			}
      if(file_ext(scorefile) %in% c("bam") == TRUE){
        scoretype="bam"
			}
			if(file_ext(scorefile) %in% c("wig","Wig") == TRUE){
				cat(scorename,": converting wig to bigWig\n")
				wig.2.bw(scorefile)
				scorefile<-paste(scorename,".bw",sep="")
			}
			if(file_ext(scorefile) %in% c("bb","bigbed","bigBed") == TRUE){
				cat(scorename,": converting bigBed to bed\n")
				bb.2.bed(scorefile)
				scorefile<-paste(scorename,".bed",sep="")
			}
			if(file_ext(scorefile) %in% c("bw","bigwig","bigWig") == TRUE){
				# cat(scorename,": converting bigWig to bedGraph\n")
				# bw.2.bg(scorefile)
				# scorefile<-paste(scorename,".bg",sep="")
        scoretype="bigWig"
			}
			if(file_ext(scorefile) %in% c("cbed","bed","cfbg","broadPeak","broadpeak","narrowPeak","narrowpeak") ==TRUE){
				scoretype="bed"
			}
			if(file_ext(scorefile) %in% c("bg","bedgraph","bedGraph") ==TRUE){
				scoretype="bedgraph"
			}

			scorename<-basename(removeext(scorefile))

			cat(scorename,": counting scores\n")
			numfrags<-filelines(scorefile)
			cat(scorename,":",numfrags,"fragments\n")

			#PRUNE READS/SCORES TO REGIONS AROUND FEATURES
			if(prunescores){
				cat(scorename,": pruning scores to regions of interest\n")
				scorefile<-bed.intersect(scorefile,featfile)
				cat(scorename,": counting pruned scores\n")
				prunedscorecount<-filelines(scorefile)
				cat(scorename,":",prunedscorecount,"scores after pruning\n")
			}

			if(scoremat){
				#get coverage/average
        if(scoretype=="bam"){
					cat(scorename,": finding coverage of",scorename,"on",featname,"\n")
          curmat.call <- pipe(paste("bedtools coverage -a",scorefile,"-b",covbedname," | sort -k4,4n | cut -f 5"),open="r")
					curmat<-t(matrix(as.numeric(readLines(curmat.call)),nrow=numwindows))
          close(curmat.call)
				}
				if(scoretype=="bed"){
					cat(scorename,": finding coverage of",scorename,"on",featname,"\n")
          curmat.call <- pipe(paste("bedtools coverage -a",scorefile,"-b",covbedname," | sort -k4,4n | cut -f 5"),open="r")
					curmat<-t(matrix(as.numeric(readLines(curmat.call)),nrow=numwindows))
          close(curmat.call)
				}
				if(scoretype=="bedgraph"){
					cat(scorename,": mapping scores in",scorename,"to",featname,"\n")
          curmat.call <- pipe(paste( "bedtools map -c 4 -o mean -null \"",bgfiller,"\" -a ",covbedname," -b ",scorefile," | sort -k4,4n | cut -f 5",sep="" ),open="r")
          curmat<-t(matrix(as.numeric(readLines(curmat.call)),nrow=numwindows))
          close(curmat.call)
		    }
        if(scoretype=="bigWig"){
					cat(scorename,": mapping scores in",scorename,"to",featname,"\n")
          curmat.call <- pipe(paste( "bigWigAverageOverBed",scorefile,covbedname,"/dev/stdout | sort -k1,1n | cut -f 6"),open="r")
          curmat<-t(matrix(as.numeric(readLines(curmat.call)),nrow=numwindows))
          close(curmat.call)
		    }

				#flip rows of negative-stranded features
				if(strand){
          curmat[negrows,1:numwindows]<-curmat[negrows,numwindows:1]
				}

				#add row names to matrix
        row.names(curmat)<-matrownames

		    #normalize by total fragments
				if(rpm & scoretype=="bed"){
					cat(scorename,": normalizing data\n")
					scalar<-1000000/numfrags
					curmat<-curmat*scalar
				}

				#mask uncovered regions with NA
				if(is.null(maskbed) == FALSE){
					cat(scorename,": masking data\n")
					curmat[maskmat==0]<-NA
				}

				#save matrix
				#outfilename<-paste(scorename,"_",featname,suffix,".mat",windowsize,sep="")
				write.mat(curmat,file=snames[[i]][j])
				cat(scorename,": matrix saved to",snames[[i]][j],"\n")

			}

			#make fragmats
			if(fragmats & scoretype == "bed"){

				cfbg<-cfbg.make(scorefile)

				cat(scorename,": creating fragmat\n")

        fmat.call <- pipe(paste("bedtools map -c 4 -o collapse -null \"NA\" -a ",covbedname," -b ",cfbg," | sort -k4,4n | cut -f 5",sep=""),open="r")
				fmat<-readLines(fmat.call)
        close(fmat.call)
				fmat<-t(matrix(fmat,nrow=numwindows))

				if(strand==TRUE){
					cat("adjusting for strand\n")
					fmat[negrows,1:numwindows]<-fmat[negrows,numwindows:1]
				}

				# save fmat
				row.names(fmat)<-matrownames
				fmatname<-paste(dname,"/",scorename,"_",featname,suffix,".fmat",windowsize,sep="")
				cat("saving matrix to",fmatname,"\n")
				write.mat(fmat,file=fmatname)
			}
		},mc.cores=cores, mc.preschedule=FALSE)

	}, mc.cores=cores2, mc.preschedule=FALSE)
	#})
	return(snames)
}
