args <- commandArgs(trailingOnly = TRUE)


inpDir=args[1]
outDir=args[2]


library("parallel")


fileAll=paste(inpDir,"/Profile_RunThrough_Regions_AllAlus_Alu_IR_CpG_Island.txt",sep="")


fileAZA=paste(inpDir,"/Profile_RunThrough_Regions_NT_Alu_IR_CpG_Island.txt",sep="")





map=read.table(fileAll,header=F)
	


dataAll=as.matrix(map)





map=read.table(fileAZA,header=F)



dataAZA=as.matrix(map)





	cexaxis=1.5
        cexlab=2


	bk=c(0,50,100)

	ymax=max(max(max(dataAll)),max(max(dataAZA)))
	ymax=0.12
	print (ymax)

	pngFile=paste(outDir,"/Profile_RunThrough_Regions_NT_AllAlus_Alu_IR_CpG_Island.png",sep="")

	png(pngFile,width=2500, height=2500,res=300)

	par(mfrow=c(1,1))

	 par(mar=c(5,7,5,2),mgp=c(5, 2, 0))


	plot(dataAll,axes=F,type='l',main="" ,xlim=c(0,100),ylim=c(0,ymax),cex.lab=cexlab,cex.axis=cexaxis,lwd=2,xlab="",ylab="CpG Island Intersection density",col="blue",
	font=2,font.lab=2)
        
	par(new = TRUE)

	plot(dataAZA,axes=F,type='l',main="" ,xlim=c(0,100),ylim=c(0,ymax),cex.lab=cexlab,cex.axis=cexaxis,lwd=2,xlab="",ylab="",col="red")


	        axis(1,at=c(0,50,100), labels=FALSE)
        text(x =c(0,50,100),-0.01,labels = c("-50Kb","Center of IR-Alus","+50Kb"), ,srt = 0, pos=1,xpd = TRUE,cex=1.3,font=2)





        axis(2,at=c(0.0,0.02,0.04,0.06,0.08,0.1), labels=FALSE)
        text(y =c(0.0,0.02,0.04,0.06,0.08,0.1),-10,labels = c("0.00","0.02","0.04","0.06","0.08","0.10"), ,srt = 90, adj=0.5, xpd = TRUE,cex=1.3,font=2)



	legend(30,ymax, c("All IR-Alus (n=746,470)","Baseline IR-Alus (n=1,040)"), lwd=c(3,3), col=c("blue","red"),lty=c(1,1), cex=1.2,text.font=2)

	box()


