args <- commandArgs(trailingOnly = TRUE)


inpDir=args[1]
outDir=args[2]



map<- read.table(paste(inpDir,"/Variants_Counts.txt",sep=""),row.names=1,header=T)

data<-as.matrix(map)

maxy=max(max(data))
maxy=((maxy%/%1000)+1)*1000

maxy=15000

print (maxy)
cnames<-colnames(data)

pngfile=paste(outDir,"/Variants_Counts.png",sep="")




names <- c("LACZ_NT_0_1","LACZ_AZA_0_1","LACZ_NT_2_1","LACZ_AZA_2_1","Sh1_NT_0_1","Sh1_AZA_0_1","Sh1_NT_2_1","Sh1_AZA_2_1")




lnames <- c(expression(paste(ADAR1^'WT',' MT')),expression(paste(ADAR1^'WT',' 5-AZA-CdR')),expression(paste(ADAR1^'WT',' MT')),expression(paste(ADAR1^'WT',' 5-AZA-CdR')),expression(paste(ADAR1^'KD',' MT')),expression(paste(ADAR1^'KD',' 5-AZA-CdR')),expression(paste(ADAR1^'KD',' MT')),expression(paste(ADAR1^'KD',' 5-AZA-CdR')))









print(names)
print(rownames(data))

data<-data[names,]
data<-t(data)



miny=0

axlimy=c(0,maxy)
        cexaxis=1
        cexlab=1.5
        cextitle=1
        cexmain=1.5

axlimx=c(1,dim(data)[2]+1)
 png(pngfile,width=2000,height=1200,bg="white",res=200)
 par(mar=c(12,7,2,2),mgp=c(5, 1, 0))
pcode=23
cexpoint=c(1,1.5,1,1,1,1,1,1,1,1,1,1)
mycol=c("black","red","black","black","black","black","black","black","black","black","black","black")
plot(c(0,0),type="n",xlim=axlimx,ylim=axlimy,cex.lab=cexlab,cex.axis=cexaxis,xlab="",ylab="Loci (Position) Counts",srt=45,cex.main=cexmain,main="A-to-I editing counts in Alus",axes=F,font.lab=2,font=2)
points(rep(1,dim(data)[1]),data[,1], col =mycol, type = "p",pch=pcode,cex=cexpoint,bg=mycol)     

points(rep(2,dim(data)[1]),data[,2], col =mycol, type = "p",pch=pcode,cex=cexpoint,bg=mycol)

points(rep(3,dim(data)[1]),data[,3], col =mycol, type = "p",pch=pcode,cex=cexpoint,bg=mycol)
points(rep(4,dim(data)[1]),data[,4], col =mycol, type = "p",pch=pcode,cex=cexpoint,bg=mycol)
points(rep(5,dim(data)[1]),data[,5], col =mycol, type = "p",pch=pcode,cex=cexpoint,bg=mycol)
points(rep(6,dim(data)[1]),data[,6], col =mycol, type = "p",pch=pcode,cex=cexpoint,bg=mycol)
points(rep(7,dim(data)[1]),data[,7], col =mycol, type = "p",pch=pcode,cex=cexpoint,bg=mycol)
points(rep(8,dim(data)[1]),data[,8], col =mycol, type = "p",pch=pcode,cex=cexpoint,bg=mycol)




#axis(1,at=c(0,18),labels=NA)
#axis(2,at=c(0,5000,10000,13000), srt=45 ,cex.axis=1)

axis(2,at=c(0,5000,10000,15000), labels=FALSE)
text(y =c(0,5000,10000,15000),0.4,labels = c(0,5000,10000,15000), ,srt = 90, adj=0.5, xpd = TRUE,cex=1.4,font=2)


lines(c(0,9),c(0,0))
lines(c(0,9),c(maxy,maxy))
lines(c(9,9),c(0,maxy))


legend(7.5,14700,c("A to I","Others"), pch=c(23,23), cex=1.3,col=c("red","black"),pt.bg=c("red","black"),pt.cex=c(1.5,1),text.font=2)


	#print(par("usr")[3])
	#text(seq(1,8,by=1),-500, srt = 60, adj= 1, xpd = TRUE,labels =lnames, cex=1.5,font=2)
	
	dev.off()





