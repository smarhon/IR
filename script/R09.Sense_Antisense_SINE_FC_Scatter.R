
args <- commandArgs(trailingOnly = TRUE)


inpDir=args[1]
outDir=args[2]



#map = read.table("/cluster/home/smarhon/ADAR/hg19.RepeatElements.GC.CpG.txt",header=T)
map = read.table(paste(inpDir,"/LacZAZA_D5_vs_LacZNT_D5_Sense_AntiSense_FC.txt",sep=""),row.names=1,header=T)






dim(map)



data<-map[grepl("SINE",map[,9]),]


dim(data)

thr=1
thr2=1
dataRev<-data[(abs(data[,1])>=1)&(data[,3]<=0.05)&((data[,4]>0.05)|(abs(data[,2])<1))&((data[,8]>=thr)|(data[,7]>=thr2)),]

dataFwd<-data[(abs(data[,2])>=1)&(data[,4]<=0.05)&((data[,3]>0.05)|(abs(data[,1])<1))&((data[,6]>=thr)|(data[,5]>=thr2)),]


dataBoth<-data[(abs(data[,1])>=1)&(abs(data[,2])>=1)&(data[,4]<=0.05)&(data[,3]<=0.05),]


dataNone<-data[((abs(data[,1])<1)&(abs(data[,2])<1))|((data[,4]>0.05)&(data[,3]>0.05))|((data[,8]<thr)&(data[,7]<thr2))|((data[,6]<thr)&(data[,5]<thr2)),]

print(dim(dataBoth))

print(dim(dataRev))


print(dim(dataFwd))

TRboth=dataBoth[(dataBoth[,1]>0)&(dataBoth[,2]>0),]

print("TRboth=")
print (dim(TRboth)[1])

TRRev=dataRev[(dataRev[,1]>0)&(dataRev[,2]>0),]

print("TRRev=")
print (dim(TRRev)[1])



TRFwd=dataFwd[(dataFwd[,1]>0)&(dataFwd[,2]>0),]

print("TRFwd=")
print (dim(TRFwd)[1])





TLboth=dataBoth[(dataBoth[,1]>0)&(dataBoth[,2]<=0),]

print("TLboth=")
print (dim(TLboth)[1])

TLRev=dataRev[(dataRev[,1]>0)&(dataRev[,2]<=0),]

print("TLRev=")
print (dim(TLRev)[1])



TLFwd=dataFwd[(dataFwd[,1]>0)&(dataFwd[,2]<=0),]

print("TLFwd=")
print (dim(TLFwd)[1])




BRboth=dataBoth[(dataBoth[,1]<=0)&(dataBoth[,2]>0),]

print("BRboth=")
print (dim(BRboth)[1])

BRRev=dataRev[(dataRev[,1]<=0)&(dataRev[,2]>0),]

print("BRRev=")
print (dim(BRRev)[1])



BRFwd=dataFwd[(dataFwd[,1]<=0)&(dataFwd[,2]>0),]

print("BRFwd=")
print (dim(BRFwd)[1])




BLboth=dataBoth[(dataBoth[,1]<=0)&(dataBoth[,2]<=0),]

print("BLboth=")
print (dim(BLboth)[1])

BLRev=dataRev[(dataRev[,1]<=0)&(dataRev[,2]<=0),]

print("BLRev=")
print (dim(BLRev)[1])



BLFwd=dataFwd[(dataFwd[,1]<=0)&(dataFwd[,2]<=0),]

print("BLFwd=")
print (dim(BLFwd)[1])













maxx=max(max(dataRev[,1]),max(dataRev[,2]),max(dataFwd[,1]),max(dataFwd[,2]),max(dataBoth[,1]),max(dataBoth[,2]),max(dataNone[,1]),max(dataNone[,2]))




minx=min(min(dataRev[,1]),min(dataRev[,2]),min(dataFwd[,1]),min(dataFwd[,2]),min(dataBoth[,1]),min(dataBoth[,2]),min(dataNone[,1]),min(dataNone[,2]))


maxx
minx

maxx=10
minx=-maxx








png(paste(outDir,"/Sense_Antisense_SINE_FC_Scatter.png",sep=""),width=1800,height=1800,bg="white",res=200)

par(mfrow=c(1,1),cex=0.9)
par(mar=c(9,10,4,1),mgp=c(6, 2, 0))




axlimx=c(minx,maxx)
axlimy=c(minx,maxx)
 
cexaxis=1.5
cexlab=2
cextitle=2
cexmain=1.5




plot(dataNone[,2],dataNone[,1],pch=20,xlim=axlimx,ylim=axlimy,cex.lab=cexlab,cex.axis=cexaxis,xlab="log(5-AZA-CdR/Mock-treatd)\nSense",ylab="AntiSense \n log(5-AZA-CdR/Mock-treatd)",
	cex.main=cexmain,main="SINEs",col="gray",font=2,font.lab=2)

points(dataRev[,2],dataRev[,1],pch=20,col="blue")

points(dataFwd[,2],dataFwd[,1],pch=20,col="red")


points(dataBoth[,2],dataBoth[,1],pch=20,col="green")


lines(c(0,0),c(minx,maxx),col="black",lty=2,lwd=1.5)
lines(c(minx,maxx),c(0,0),col="black",lty=2,lwd=1.5)



legend(minx,maxx, c("Both NS","AntiSense Sig.","Sense Sig.","Both Sig."),pch=20,pt.cex=3,col=c("gray","blue","red","green"),cex=1.3,text.font=2)



dev.off()

