# TODO: Add comment
# 
# Author: mike
###############################################################################


.readTblfromHIVdb<-function(tblname)
{
	tblname<-paste(tblname,".csv",sep="")
	txdb_file <- system.file("extdata/parsed",tblname ,package="HIV.db")
    read.table(txdb_file,header = TRUE,sep=",",as.is = TRUE)
	
}

##the ranges information for each sequence are not meaningful here 
##since it is simply calculated from the sequential positions of these sequences from fasta file
.readAASeq<-function(fileName)
{
	seqfile <- system.file("extdata/parsed",fileName ,package="HIV.db")
	ret<-read.AAStringSet(seqfile)
	names(ret)<-tolower(sub("HXB2_AA_","",names(ret)))#strip the prefix
	ret
	
}

.readDNASeq<-function(fileName)
{
	seqfile <- system.file("extdata/parsed",fileName ,package="HIV.db")
	read.DNAStringSet(seqfile)
}


loadFeatures<-function()
{
#	browser()
	HIV_db<-new.env(hash=TRUE, parent=emptyenv())
	ret<-.readTblfromHIVdb(tblname="hxb2Table")
	ret<-subset(ret,t_category!="internal")
	ret1<-.readTblfromHIVdb("gpadditional")
#	ret1$t_name<-tolower(ret1$t_name)
	ret<-rbind(ret,ret1)
	
	assign("hxb2Table",ret,HIV_db)
	ret<-.readTblfromHIVdb("antibody")
	ret<-rbind(ret,.readTblfromHIVdb("bindings"))
	ret$X<-rownames(ret)
	ret$MAb.Name<-sub("&","",ret$MAb.Name)
	assign("antibody",ret,HIV_db)
	assign("hxb2AA",.readAASeq("hxb2_AA.fasta"),HIV_db)
	assign("hxb2DNA",.readDNASeq("hxb2_DNA.fasta"),HIV_db)
	HIV_db
#	new("HivFeature",FeatureID=as.integer(0),name="ALL",
#			category="hxb2",start=as.integer(1),
#			end=as.integer(9719),HIV_db=HIV_db)
}

lsCategory<-function(HIV_db)
{
	tbl<-get("hxb2Table",HIV_db)
	unique(tbl$t_category)
}
#when range is provided,return all the features that have intersections with the range
.getFeature<-function(HIV_db,category=NULL,name=NULL,start=NULL,end=NULL,frame=NULL,...)
{
#	browser()
	###if category is epitope then call getEpitope method to query antibodybinding table
	if(!is.null(category)&&tolower(category)=="epitope")
		return(getEpitope(HIV_db,name=name,start=start,end=end,frame=frame,...))
	####otherwise, query hxb2Table
	tbl<-get("hxb2Table",HIV_db)
	ret<-tbl
	if(!is.null(frame))
		ret<-subset(ret,t_frame%in%frame)
	if(!is.null(category))
		ret<-subset(ret,t_category%in%category)
	if(!is.null(name))
		ret<-subset(ret,t_name%in%name)

	if(!is.null(start))
	{
		ret<-subset(ret,t_end>=start)
#		len<-nrow(ret)

	}
	if(!is.null(end))
	{
		ret<-subset(ret,t_start<=end)
#		len<-nrow(ret)
#		ret<-ret[order(ret$t_end)[1],]#get the closest one
#		len<-nrow(ret)
#		ret<-ret[order(ret$t_start)[len],]#get the closest one	
	}
#	browser()
	ret<-HivFeature(ret,HIV_db)
	
#	if(children)
#		ret<-lapply(ret,function(feature){	
#			c(feature,getChildren(feature,recursive=TRUE))
#		})
#	browser()

	#	if(length(ret)==1)
#		ret[[1]]
#	else
#		ret
	ret
	
}
