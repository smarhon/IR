
args <- commandArgs(trailingOnly = TRUE)


inpDir=args[1]
outDir=args[2]


#library("parallel")
#library(afex)

cov="2.5"
file1=paste(inpDir,"/LacZAZA_AverageProfile_Protected_Regions_Inverted_Repeats_Reverse_AZA.txt",sep="")

file2=paste(inpDir,"/LacZAZA_AverageProfile_Protected_Regions_nonInverted_Repeats_Reverse_AZA.txt",sep="")


file3=paste(inpDir,"/LacZNT_AverageProfile_Protected_Regions_Inverted_Repeats_Reverse_AZA.txt",sep="")

file4=paste(inpDir,"/LacZNT_AverageProfile_Protected_Regions_nonInverted_Repeats_Reverse_AZA.txt",sep="")

file5=paste(inpDir,"/LacZAZA_Average_Peaks_Protected_Regions_Inverted_Repeats_Reverse_AZA.txt",sep="")

file6=paste(inpDir,"/LacZNT_Average_Peaks_Protected_Regions_Inverted_Repeats_Reverse_AZA.txt",sep="")




map1=read.table(file1,header=T)

data1=as.matrix(map1)

map2=read.table(file2,header=T)

data2=as.matrix(map2)


data<-cbind(data1,data2)

map1=read.table(file3,header=T)

data1=as.matrix(map1)

map2=read.table(file4,header=T)

data2=as.matrix(map2)


dataNT<-cbind(data1,data2)

map1=read.table(file5,row.names=1,header=T)
map2=read.table(file6,row.names=1,header=T)

data1=as.matrix(map1)
data2=as.matrix(map2)
data2<-data2[rownames(data1),]
data2=as.matrix(data2)
dataP<-cbind(data1[,1],data2[,1])
pval=c()
for (i in 1:1000){
	rn<-rownames(dataP)
	rn<-sample(rn)
	dataP2<-dataP[rn[1:200],]



	y<-wilcox.test(dataP2[,1], dataP2[,2],paired=TRUE)

	#print(y[3])
	p=y[3]
	pval<-c(pval,as.numeric(p))
}

p=max(pval)
p=signif(p, digits = 2)
print(p)


dataPNon<-cbind(data1[,2],data2[,2])
pvalNon<-c()
for (i in 1:1000){
        rn<-rownames(dataPNon)
        rn<-sample(rn)
        dataP2<-dataPNon[rn[1:200],]



        y<-wilcox.test(dataP2[,1], dataP2[,2],paired=TRUE)

        #print(y[3])
        pNon=y[3]
        pvalNon<-c(pvalNon,as.numeric(pNon))
}



pNon=max(pvalNon)
pNon=signif(pNon, digits = 2)
print(pNon)
pNon=min(pvalNon)
pNon=signif(pNon, digits = 2)
print(pNon)


y<-ks.test(dataP[,1], dataP[,2],exact=TRUE)
print (y)
print("pTot=")
print(y[3])
pTot=as.numeric(y[3])
print(pTot)




y<-wilcox.test(dataP[,1], dataP[,2],paired=TRUE,exact=TRUE,correct=TRUE)
print (y)
print("pTot=")
print(y[3])
pTot=as.numeric(y[3])
print(pTot)
pTot=signif(pTot, digits = 2)
print(pTot)



y<-wilcox.test(data1[,2], data2[,2],paired=TRUE)
pTotNon=as.numeric(y[3])

pTotNon=signif(pTotNon, digits = 2)
print(pTotNon)

cexaxis=1.5
cexlab=1.5



maxx=dim(data)[1]
maxy=max(max(max(data)),max(max(dataNT)))+5
miny=0
bk=c(0,maxx/2,maxx)



pngFile=paste(outDir,"/AverageProfile_Protected_Regions_IR-nonIR_Reverse_AZA_profile.png",sep="")


cexmain=2

png(pngFile,width=1500, height=1500,res=200)
par(mar=c(5,6,5,2),mgp=c(4, 1, 0))


plot(data[,1],axes=F,type='l',lwd=3,ylim=c(0,maxy),cex.main=cexmain,cex.lab=cexlab,main="IR Transcripts",cex.axis=cexaxis,xlab="Genome Position",ylab="MDA5-protected RNA (RPM)",
col="#666666",font=2,font.lab=2,font.axis=2)
lines(dataNT[,1],col="#952A50",lwd=3,type='l')

axis(1,at=bk, cex.axis=1.5,labels=c("Start","Center","End"),font=2)
axis(2,cex.axis=1.5,font=2)


legend(400,maxy, c("5-AZA-CdR","Mock-treated"), lwd=c(3,3), col=c("#666666","#952A50"),lty=c(1,1), cex=1.5,text.font=2)
box()

print(dim(dataP))
mydata=list(dataP[,2],dataP[,1])
maxy=max(max(dataP))
miny=min(min(dataP))
y1=boxplot.stats(dataP[,1], coef =1.5, do.conf=FALSE, do.out=TRUE)
print (y1$stats)
maxy1=y1$stats[5]
miny1=y1$stats[1]
y2=boxplot.stats(dataP[,2], coef =1.5, do.conf=FALSE, do.out=TRUE)
print (y2$stats)
maxy2=y2$stats[5]
miny2=y2$stats[1]



maxy=max(maxy1,maxy2)+3
miny=min(miny1,miny2)


pngFile=paste(outDir,"/AverageProfile_Protected_Regions_IR-nonIR_Reverse_AZA_boxplot.png",sep="")


cexmain=2

png(pngFile,width=1500, height=1500,res=200)
par(mar=c(6,6,5,2),mgp=c(4, 1, 0))

boxplot(mydata, outline=FALSE,boxwex=0.25,main="IR Transcripts",ylim=c(miny,maxy),names=c("Mock-treated","5-AZA-CdR"),col=c("#952A50","#666666"),
cex.main=cexmain,cex.axis=cexaxis,cex.lab=cexlab,xlab="",ylab="MDA5-protected RNA (RPM)",font=2,par(font.lab=2,font=2,font.axis=2))
#mtext(paste("p = ",toString(pTot),sep="") ,adj=0.5,cex=1.5,line=-3,side=3,font=2)
mtext("p<2.2e-16" ,adj=0.5,cex=1.5,line=-3,side=3,font=2)


#mtext(paste("p = ",toString(pTotNon),sep="") ,adj=0.75,cex=1.5,line=-20,side=3)

dev.off()




