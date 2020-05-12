args <- commandArgs(trailingOnly = TRUE)


inpDir=args[1]
outDir=args[2]

map = read.table(paste(inpDir,"/PolyA_Counts_16Motifs_BamReads_AZA_AllMotifs_histogram.txt",sep=""),row.name=1,header=T)


data<-as.matrix(map)


o<-order(rowSums(data),decreasing=TRUE)

data<-data[o,]






x=seq(0,150)

cexaxis=1.5
cexlab=2
cextitle=2
cexmain=1.5



png(paste(outDir,"/PolyA_Counts_16Motifs_BamReads_AZA_AllMotifs_histogram.png",sep=""),width=1200,height=2400,bg="white", res=200)



par(mfrow=c(6,3),cex=0.9)
par(mar=c(3,3,2,2),mgp=c(2, 1, 0))

for (i in 1:dim(data)[1]){


maxy=max(data[i,])


#maxy=400
axlimy=c(0,maxy)


cexaxis=1
cexlab=1
cextitle=1
cexmain=1

minx=0



axlimx=c(0,150)

print(maxy)


plot(x,data[i,],lwd=1.5,type="l",cex=1,ylab="",xlab="",col="red",main=rownames(data)[i],xlim=axlimx,ylim=axlimy,cex.main=1.3,cex.lab=1,cex.axis=1,axes=F)

axis(1,at=c(0,150),labels=c("0","150bp"))

axis(2,at=c(0,maxy))


#legend(0.8,maxy,c("A--G","A--C","A--T","C--G","T--C","G--T"), lty=c(1,1,1,1,1,1), cex=1, lwd=c(3,3,3,3,3,3),col=c("red","blue","brown","black","darkgreen","purple"))

box()

}


dev.off()



