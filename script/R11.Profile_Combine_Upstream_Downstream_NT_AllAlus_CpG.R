args <- commandArgs(trailingOnly = TRUE)


inpDir=args[1]
outDir=args[2]


library("parallel")




fileAll=paste(inpDir,"/Profile_Upstream_Downstream_IRregion_Regions_AllAlus_Alu_IR_CpG.txt",sep="")




fileAZA=paste(inpDir,"/Profile_Upstream_Downstream_IRregion_Regions_NT_Alu_IR_CpG.txt",sep="")



map=read.table(fileAll,header=F)
	


dataAll=as.matrix(map)
print (dim(dataAll))

dataAll<-dataAll[,1]

dataAll<-c(dataAll[1:1000],dataAll[1201:2200])

print (length(dataAll))






map=read.table(fileAZA,header=F)


dataAZA=as.matrix(map)

dataAZA=as.matrix(dataAZA)



dataAZA<-dataAZA[,1]

dataAZA<-c(dataAZA[1:1000],dataAZA[1201:2200])



	cexaxis=1.5
        cexlab=2


	bk=c(0,1000,2000)

	ymax=max(max(max(dataAll)),max(max(dataAZA)))+0.005

	print (ymax)

	pngFile=paste(outDir,"/Profile_Upstream_Downstream_Regions_AllAlus_NT_Alu_IR_CpG.png",sep="")

	png(pngFile,width=2500, height=2500,res=300)

	par(mfrow=c(1,1))

	 par(mar=c(5,7,5,2),mgp=c(5, 2, 0))






	
        plot(dataAZA,axes=F,type='l',main="" ,xlim=c(0,2000),ylim=c(0,ymax),cex.lab=cexlab,cex.axis=cexaxis,lwd=2,xlab="",ylab="CpG  density",col="red",font=2,font.lab=2)

        par(new = TRUE)

        plot(dataAll,axes=F,type='l',main="" ,xlim=c(0,2000),ylim=c(0,ymax),cex.lab=cexlab,cex.axis=cexaxis,lwd=2,xlab="",ylab="",col="blue")



	axis(1,at=c(0,1000,2000), labels=FALSE)
        text(x =c(0,1000,2000),-0.0015,labels = c("-50Kb","","+50Kb"), ,srt = 0, pos=1,xpd = TRUE,cex=1.3,font=2)





        axis(2,at=c(0.0,0.005,0.01,0.015,0.02,0.025), labels=FALSE)
        text(y =c(0.0,0.005,0.01,0.015,0.02,0.025),-200,labels = c("0.000","0.005","0.010","0.015","0.020","0.025"), ,srt = 90, adj=0.5, xpd = TRUE,cex=1.3,font=2)



        legend(600,ymax, c("All IR-Alus (n=746,470)","Baseline IR-Alus (n=1,040)"), lwd=c(3,3), col=c("blue","red"),lty=c(1,1),cex=1.2,text.font=2)

        box()




