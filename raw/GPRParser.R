#To read xls files
library(gdata)



JPTdesign<-read.xls("~/Programs/git/HIV.db/inst/extdata/JPT_Annotation_gp160_slide_18Mar11.xls")
#########################################################
##1.pair wise alignment of consecutive peptide sequence to construct entire clade sequence
#######################################################
res<-vector("list",7)
names(res)<-c("M","A","B","C","D","CRF01","CRF02")
for(c in 1:7)
{
	clades<-as.character(JPTdesign[,6+c])
	pepClade<-as.character(JPTdesign$Name[!is.na(clades)])
	start<-rep(0,length(pepClade))
	end<-rep(0,length(pepClade))
	tmpSeq<-do.call(function(x,...){paste(x,sep="",...)},as.list(rep("-",860)))
	start[1]<-1
	end[1]<-start[1]+nchar(pepClade[1])-1
	substr(tmpSeq,start[1],end[1])<-pepClade[1]
	for(i in 2:length(pepClade))
	{
		tmp<-pairwiseAlignment(pepClade[i],pepClade[i-1],type="overlap")
		start[i]<-start[i-1]+start(tmp@subject@range)-1
		end[i]<-start[i]+nchar(pepClade[i])-1
		substr(tmpSeq,start[i],end[i])<-pepClade[i]
	}
	### I remove the unecessary characters at the end
	res[[c]]<-list(peptide=pepClade,start=start,end=end,consensus=substr(tmpSeq,start[1],end[length(pepClade)]))
}
