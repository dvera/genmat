deploy<-function(message="no message"){
	library(travis)
	res <- system(paste0("git -C /lustre/maize/home/dlv04c/software/r/genmat/ add /lustre/maize/home/dlv04c/software/r/genmat/ &&\
	git -C /lustre/maize/home/dlv04c/software/r/genmat/ commit -a -m '",message,"' &&\
	git -C /lustre/maize/home/dlv04c/software/r/genmat/ push"))
	library(devtools)
	detach("package:genmat",unload=T)
	install_github("dvera/genmat")
	library(genmat)
}
