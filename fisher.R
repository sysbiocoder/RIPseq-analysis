#Fischer Exact test
require(stringr)

# Total number of reads
rm(list = ls())
# Read feature count 
	args = commandArgs() 
	ix = grep("--args", args, fixed = TRUE)
 	if (ix < length(args))
	{
    		for (i in (ix+1):length(args))
		 {
    	 		print(paste("argument #",i-ix,"-->",args[i], sep=" "), sep = "\n")  
		}
		#Parameters passed
		ctl<-args[8]
		ctl
		tmt <- args[9]
		tmt
		infile <- args[10]
		infile
		outfolder <- args[11]
		outfolder
	} 
library(edgeR)#
infile
ctl
tmt
data <- read.table(infile, header=TRUE,sep="\t")
outfolder
ctcol<-strtoi(ctl)
tmcol <- strtoi(tmt)	
f1<-names(data)[ctcol]
f2<-names(data)[tmcol]
len <- "Length"
ofile <- paste(outfolder,str_trim(f1),"_",str_trim(f2),sep="")
newdata <- data[,c(f1,f2,len)]
row.names(newdata) <- data$Geneid
head(newdata)

newdata$CTL <- newdata[[f1]] #CTL
newdata$TMT <- newdata[[f2]] #TMT

df <- data.frame(newdata$CTL,newdata$TMT)
row.names(df) <- row.names(newdata) 

counts <-  df[rowSums(df) != 0,]


dge <- DGEList(counts=newdata[,1:2],genes=data.frame(Length=newdata$Length))
dge <- calcNormFactors(dge)
RPKM <- rpkm(dge)
head(RPKM)


#------------------------------------Precalculations--------------------------------------------------------------------

ftres= c()
geneid= c()
for( i in 1: nrow(counts))
{
res <- matrix(c(counts[i,1], counts[i,2],sum(counts[,1]) -counts[i,1],sum(counts[,2]) -counts[i,2]),ncol=2,byrow=TRUE)
ft <- fisher.test(res)
#print(counts[i,])
ftres[i] <- ft$p.value
geneid[i] <- row.names(counts[i,])
}
head(counts)

#----------------------------------------Eqn-2----------------------------------------------------------------------------------

Pvalue <- ftres
FDR <-  p.adjust(Pvalue,method="BH")
#---------------------FC----calculation-------------------------------------------------------------------------------------------
m1 <- data.frame(geneid,Pvalue,FDR)
log2FC <- log2(RPKM[,2]/RPKM[,1])
log2FC <- data.frame(log2FC,RPKM[,1],RPKM[,2])
log2FC$geneid <- row.names(log2FC)
counts$geneid <- row.names(counts)
lfc <- merge(log2FC, counts,by="geneid")
#----------------------------------------Master file----------------------------------------------------------------------------------
acr <-merge(lfc,m1,by="geneid")
head
file1 <- paste(ofile,"Master.txt",sep="")
write.table(acr,file1,sep="\t",row.names=F,quote=F)
