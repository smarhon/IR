
args <- commandArgs(trailingOnly = TRUE)

inpDir=args[1]
outDir=args[2]

library(plotrix)
cov="2.5"
map = read.table(paste(inpDir,"/Inverted_Repeats_Reverse_NT_ratio_Histogram.txt",sep=""),header=T)


dataInv<-as.matrix(map)

map = read.table(paste(inpDir,"/nonInverted_Repeats_Reverse_NT_ratio_Histogram.txt",sep=""),header=T)

lmain=c("Mock-treated","5-AZA-CdR",expression(paste(ADAR1^'KD',' Mock-treated')),expression(paste(ADAR1^'KD',' 5-AZA-CdR')))

xlabels=c("log10 (MDA5-protected/Total CytoRNA)","log10 (MDA5-protected/Total CytoRNA)","log10 (MDA5-protected/Total CytoRNA)","log10 (MDA5-protected/Total CytoRNA)")



xlabels=c(expression('log'[10]*' (MDA5-protected/Total CytoRNA)'),expression('log'[10]*' (MDA5-protected/Total CytoRNA)'),expression('log'[10]*' (MDA5-protected/Total CytoRNA)'),expression('log'[10]*' (MDA5-protected/Total CytoRNA)'))

Ids=c("NT","AZA","shNT","shAZA")







datanonInv<-as.matrix(map)

dataInv<-dataInv[(dataInv[,1]>-3.0)&(dataInv[,1]<3.0),]

datanonInv<-datanonInv[(datanonInv[,1]>-3.0)&(datanonInv[,1]<3.0),]

xbklab<- seq(min(dataInv[,1]),max(dataInv[,1]),1)

xbk<-seq(1,dim(dataInv)[1],10)

cexaxis=1.5
cexlab=2
cextitle=2
cexmain=1.5



Ind<-c(2,3)
for (i in Ind){
	id=Ids[i-1]


	png(paste(outDir,"/Inverted_Repeats_Reverse_NT_ratio_Histogram-",id,".png",sep=""),width=1500,height=1500,bg="white",res=200)

	par(mar=c(6,6,4,1),mgp=c(4, 2, 0))

	maxy=max(max(max(dataInv[,2:5])),max(max(datanonInv[,2:5])))
	print (maxy)

	maxy=max(max(dataInv[,i]),max(datanonInv[,i]))

	print (maxy)
	maxy=(maxy%/%100+1)*100
	maxy=1200   
	axlimy=c(0,maxy)


	axlimy=c(0,maxy)

	cexaxis=1.5
	cexlab=2
	cextitle=2
	cexmain=1.5

	minx=min(dataInv[,1])

	maxx=max(dataInv[,1])

	mix=-3
	maxx=3

	axlimx=c(minx,maxx)

	print(maxy)

	mycolRed <- rgb(255, 0, 0, max = 255, alpha = 150, names = "red70")
	mycolblue <- rgb(0, 0, 255, max = 255, alpha = 200, names = "blue70")

	to=maxy-200
	from=600
	plot(datanonInv[,1],datanonInv[,i],lwd=3,xticlab=NA,yticlab=NA,type="h",ylim=axlimy,cex=1.5,border=mycolblue,main=lmain[i-1],font=2,
			ylab="Counts",xlab=xlabels[i-1],col=mycolblue,cex.main=2,cex.lab=1.8,cex.axis=1.5,axes=F,font.lab=2,main.font=2)


	par(new = TRUE)

	plot(dataInv[,1],dataInv[,i],lwd=3,xticlab=NA, yticlab=NA,type="h",ylim=axlimy,col=mycolRed,ylab="",xlab="",cex.lab=1.5,cex.axis=1.5,font=2,axes=F,main.font=2)


	axis(1,at=c(-3,-2,-1,0,1,2,3), labels=FALSE)
	text(x =c(-3,-2,-1,0,1,2,3),-100,labels = c(-3,-2,-1,0,1,2,3), ,srt = 0, pos=1,xpd = TRUE,cex=1.5,font=2)

	axis(2,at=c(0,400,800,1200), labels=FALSE)
	text(y =c(0,400,800,1200),-3.5,labels = c(0,400,800,1200), ,srt = 90, adj=0.5, xpd = TRUE,cex=1.5,font=2)

	legend(0.5,maxy-10,c("Non-IR Alus","IR Alus"), lty=c(1,1), cex=1.3, lwd=c(6,6),col=c(mycolblue,mycolRed),text.font=2)



	box()


}


dev.off()



