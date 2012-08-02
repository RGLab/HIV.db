pathToData<-"~/workspace/HIVdb/Slide-Comparison/data"


lFolders<-list.files(pathToData)
#MAKE MAPPING
#lFolders<-c("2F5", "7B2", "gp1403 Pre and Post", "HIVIG and IVIG", "T141485 v0 and v7", "Mk7403 pre and post")
for(idx in 1:length(lFolders))
{
	curFolder<-paste(pathToData, lFolders[idx], sep="/")
	print(curFolder)

	
	fNameList<-c()
	ptidList<-c()
	visitList<-c()
	batchList<-c()
	pmtList<-c()
	slideList<-c()
	mappingFile<-paste(curFolder, "mapping.csv", sep="/")
	unlink(mappingFile)
	
	lFiles<-list.files(curFolder)
	print(length(lFiles))
	
	for(idx2 in 1:length(lFiles))
	{

		fName<-substr(lFiles[idx2], 0, nchar(lFiles[idx2])-4)
		fNameList<-c(fNameList, fName)	
		lSplit<-unlist(strsplit(fName, split=c("-|_")))
		ptidList<-c(ptidList, fName)
		if(lSplit[3]=="Pre" | lSplit[2]=="IVIG" | lSplit[4]=="pre" | lSplit[3]=="v0")
			visit<-"Pre"
		else
			visit<-"Post"
		visitList<-c(visitList, visit)
		batchList<-c(batchList, lSplit[1])
		pmtList<-c(pmtList, as.numeric(lSplit[6]))
		slideList<-c(slideList, lSplit[5])
		
	}
	#fName,ptid,visit,batch,pmt,slide
	df<-data.frame(filename=fNameList, ptid=ptidList, visit=visitList, batch=batchList, pmt=pmtList, slide=slideList)
	write.csv(df, file=mappingFile,row.names=FALSE)
	cat("Mapping file for ",lFolders[idx],"written in", mappingFile,"\n")
}



#MAKE PEP_
#read data
pathToDoc<-"~/workspace/HIVdb/Slide-Comparison/doc"
JPTxls<-paste(pathToDoc, "JPT_2105.xls", sep="/")
xls<-read.xls(JPTxls) #library(gdata)

#remove triplicate, remove control rows and order
nxls<-xls[1:(nrow(xls)/3),]
nxls<-nxls[!nxls$Peptide.align=="control",]
nxls<-nxls[with(nxls, order(Annotation.1)),]
#Extract info ffom annotation col
annoList<-strsplit(as.character(nxls$Annotation.1), split="_")
protList<-unlist(lapply(annoList, function(x){x[5]}))
Uprot<-unique(protList)
pep_DB<-RangedData(space=Uprot)
for(idx in 1:length(Uprot))
{
	prot<-Uprot[idx]
	print(prot)
	prot_idx<-which(protList==prot)
	peptides<-nxls[prot_idx,]$ID
	peptidePos<-unlist(lapply(annoList[prot_idx], function(x){x[3]}))
	peptideNb<-unlist(lapply(annoList[prot_idx], function(x){x[1]}))
	aligned<-nxls[prot_idx,]$Peptide.align.
	hxb2_equi<-nxls[prot_idx,]$HXB2.alig.
	#I made sure all peptides are actually of length 15: which(nchar(as.character(rNames))!=15)
	ranges<-IRanges(start=as.integer(peptidePos)+1, width=15)
	zMatrix<-.zSum(as.character(peptides))
	pep_prot<-RangedData(ranges, aligned, peptideNb, zMatrix, space=prot)
	rownames(pep_prot)<-peptides
	pep_DB<-rbind(pep_DB,pep_prot)
}




#for env
#Apparently: peptides of 15aa, overlapping by 11
#hxb2 MRVKEKYQHLWRWGWRWGTMLLGMLMICSATEKLWVTV
#      GIRKNYQHLWRWGTM
rNames<-nxls[ENV_idx,]$ID
peptidePos<-unlist(lapply(annoList[ENV_idx], function(x){x[3]}))
peptideNb<-unlist(lapply(annoList[ENV_idx], function(x){x[1]}))
aligned<-nxls[ENV_idx,]$Peptide.align.
trimmed<-nxls[ENV_idx,]$HXB2.alig.

ENVranges<-IRanges(start=as.integer(peptidePos)+1, width=15)
pep_machin<-RangedData(ENVranges, aligned, peptideNb, space="env")
pep_machin<-RangedData(ranges(pep_machin), values(pep_machin),aligned)
pep_machin<-RangedData(ranges(pep_machin), values(pep_machin),peptideNb)
