args <- commandArgs(trailingOnly = TRUE)


inpDir=args[1]
outDir=args[2]

library(RColorBrewer)
library(MASS)
#library(mining)
cov="2.5"
map = read.table(paste(inpDir,"/Scatter_Inverted_Repeats_Reverse_NT_ratio.txt",sep=""),header=T)

map<-map[,3:10]

dataInv<-as.matrix(map)

dataInv<-log10(dataInv+0.00001)


lmain=c("Mock-treated","5-AZA-CdR",expression(paste(ADAR1^'KD',' Mock-treated')),expression(paste(ADAR1^'KD',' 5-AZA-CdR')))

xlabels=c(expression('log'[10]*' (MDA5-protected/Total CytoRNA)'),expression('log'[10]*' (MDA5-protected/Total CytoRNA)'),expression('log'[10]*' (MDA5-protected/Total CytoRNA)'),expression('log'[10]*' (MDA5-protected/Total CytoRNA)'))

Ids=c("NT","AZA","shNT","shAZA")

cexaxis=1.5
cexlab=2
cextitle=2
cexmain=1.5


for (i in 1:2){

	id=Ids[i]


	png(paste(outDir,"/Inverted_Repeats_Reverse_NT_Scatter-",id,".png",sep=""),width=1500,height=1500,bg="white",res=200)

	par(mar=c(7,9,6,1),mgp=c(4, 1, 0))

	data<-dataInv[,c(i,i+4)]

	maxx=max(max(data[,1]))
	minx=min(min(data[,1]))

	maxy=max(max(data[,2]))
	miny=min(min(data[,2]))


	if (maxx>abs(minx)){
	        minx=-maxx
	}else{
	        maxx=abs(minx)

	}



	if (maxy>abs(miny)){
	        miny=-maxy
	}else{
	      	maxy=abs(miny)
	}




	axlimy=c(miny,maxy)


	axlimx=c(minx,maxx)

	cexaxis=1.5
	cexlab=2
	cextitle=2
	cexmain=1.5


	mycolRed <- rgb(255, 0, 0, max = 255, alpha = 100, names = "red50")
	mycolblue <- rgb(0, 0, 255, max = 255, alpha = 100, names = "blue50")


	k<- 11
	z <- kde2d(data[,1], data[,2], n=100)
	#print(z)
	my.cols <- rev(brewer.pal(k, "RdYlBu"))



	plot(data[,1],data[,2],pch=20,cex=1.5,ylab=xlabels[i],main="",xlab=xlabels[i],
		col="blue",ylim=axlimy,xlim=axlimx,cex.main=2,cex.lab=1.5,cex.axis=1.5,font=2)
	contour(z, drawlabels=FALSE, nlevels=k, lwd=3,col=my.cols, add=TRUE)

	mtext(lmain[i], line = 3,cex=2,font=2)

	mtext("IR Alus", line = 1,cex=1.5,font=2)




	mtext("Repeat2",side=1,line=5.5,adj=0.5,cex=1.5,font=2)


	mtext("Repeat1",side=2,line=7,adj=0.5,cex=1.5,font=2)


	lines(c(0,0),c(miny,maxy),col="black",lty=2,lwd=1.5)
	lines(c(minx,maxx),c(0,0),col="black",lty=2,lwd=1.5)
	box()


}


dev.off()



