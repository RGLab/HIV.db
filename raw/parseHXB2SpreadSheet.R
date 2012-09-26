# TODO: Add comment
# 
# Author: mike
###############################################################################

#############################
#parse hxb2 first
##############################
library(gdata)
library(HIV.db)
path=system.file("raw",package="HIV.db")
hxb2<-read.xls(paste(path,"hxb2.xls",sep="/"))

gene<-subset(hxb2,!gene%in%c("",1:1000),select=c("gene"))
gene

write.csv(gene,file=paste(path,"output/gene.csv",sep=""))

protein<-subset(hxb2,!protein.1%in%c("",1:1000),select=c("protein.1"))
protein
write.csv(protein,file=paste(path,"output/protein.csv",sep=""))

RNA<-subset(hxb2,!RNA%in%c("",1:1000),select=c("RNA"))
RNA<-subset(RNA,!grepl("http",RNA))
RNA
write.csv(RNA,file=paste(path,"output/RNA.csv",sep=""))

protein.site<-subset(hxb2,!protein.site%in%c("",1:1000),select=c("protein.site"))
protein.site
write.csv(protein.site,file=paste(path,"output/protein.site.csv",sep=""))


#######################################################
#then antibody table based on manually parsed hxb2 table 
####################################################
path="~/rglab/workspace/HIV.db/inst/extdata/parsed"
HIV_db<-read.csv(paste(path,"hxb2Table.csv",sep="/"),as.is=TRUE)
Antibody_db<-read.csv(paste(path,"ab_summary_manually corrected.csv",sep="/"),as.is=TRUE)
#Antibody_db$Epitope<-as.character(Antibody_db$Epitope)

#fill empty name
Antibody_db[Antibody_db$MAb.Name=="",]$MAb.Name<-"unknown"
#separate records with two binding sites
multiSites<-!is.na(Antibody_db$HXB2.start.2)
singleSite<-Antibody_db[!multiSites,c(-6,-7)]

site1<-Antibody_db[multiSites,c(-6,-7)]
for(i in 1:nrow(site1))
{
	seqs<-strsplit(as.character(site1[i,]$Epitope),"\\+")[[1]]
	site1[i,]$Epitope<-seqs[1]
	
}

site2<-Antibody_db[multiSites,c(-4,-5)]
for(i in 1:nrow(site2))
{
	seqs<-strsplit(as.character(site2[i,]$Epitope),"\\+")[[1]]
	site2[i,]$Epitope<-seqs[2]
	
}

#combine them
colnames(site2)[4:5]<-colnames(site1)[4:5]
Antibody_db<-rbind(singleSite,rbind(site1,site2))
head(Antibody_db)
#convert peptide position to hxb2 coorndiates 
#(Tat,Rev have two seperate region, the minor one was marked by added 1to the name
# such as Tat1,Rev1 in both orignal antibody table and parsed hxb2Table)
Antibody_db<-merge(Antibody_db,HIV_db,by.x="Protein",by.y="t_name")
Antibody_db$HXB2.start<- Antibody_db$t_start+(Antibody_db$HXB2.start-1)*3
Antibody_db$HXB2.end<- Antibody_db$t_start+(Antibody_db$HXB2.end-1)*3
head(Antibody_db)

path="~/rglab/workspace/HIV.db/inst/extdata/parsed"
write.csv(Antibody_db[,c(1:7,13)],file=paste(path,"antibody.csv",sep="/"),row.names=FALSE)
##maually cal the coornidates for p24-p2p7p1p6 and add to parsed table

#######################################################
#convert the manually collected antibody table into the same format from LANL 
####################################################
path="~/rglab/workspace/HIV.db/inst/extdata/parsed"
path1="~/rglab/workspace/HIV.db/inst/extdata/raw"
hxb2Table<-read.csv(paste(path,"hxb2Table.csv",sep="/"),as.is=TRUE)
Antibody<-read.csv(paste(path,"antibody.csv",sep="/"),as.is=TRUE)
binding<-read.csv(paste(path1,"binding.csv",sep="/"),as.is=TRUE)

binding1<-merge(binding,HIV_db$"hxb2Table",by.x="binds_to",by.y="t_name")
env_start<-subset(hxb2Table,t_name=="env")$t_start
#convert relative AA position to the abs HXB2 coordinates
binding1$HXB2.start<-env_start+(binding1$start-1)*3
binding1$HXB2.end<-env_start+(binding1$stop-1)*3
#extract AA seq from binding protein feature
binding1$Epitope<-
		apply(binding1,1,function(x){
#					browser()
					substr(
					getAA(getFeature(HIV_db,name="env"))
					,x["start"],x["stop"])
				})
#change the column names and append the necessary columns
names(binding1)[which(names(binding1)=="antibody")]<-"MAb.Name"
names(binding1)[which(names(binding1)=="binds_to")]<-"Protein"
binding1$Subtype<-""
binding1$Species<-"mouse"

write.csv(binding1[,names(Antibody)],file=paste(path,"bindings.csv",sep="/"),row.names=FALSE)



###preprocess gp120additional_raw
path="~/rglab/workspace/HIV.db/inst/extdata/parsed"
path1="~/rglab/workspace/HIV.db/inst/extdata/raw"
gp<-read.csv(paste(path1,"gpadditional_raw.csv",sep="/"),as.is=TRUE)
hxb2Table<-read.csv(paste(path,"hxb2Table.csv",sep="/"),as.is=TRUE)
ID_start<-max(hxb2Table$t_ID)
gp$t_ID<-ID_start+(1:nrow(gp))
env_start<-subset(hxb2Table,t_name=="env")$t_start
#convert relative AA position to the abs HXB2 coordinates
gp$t_start<-env_start+(gp$t_start-1)*3
gp$t_end<-env_start+(gp$t_end-1)*3
write.csv(gp,file=paste(path,"gpadditional.csv",sep="/"),row.names=FALSE)



