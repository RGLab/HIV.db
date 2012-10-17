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
.readAASeq<-function(fileName, genome="hxb2")
{
	seqfile <- system.file("extdata/parsed",fileName ,package="HIV.db")
	ret<-readAAStringSet(seqfile)
	prefix<-paste(toupper(genome), "_AA_", sep="")
	names(ret)<-tolower(sub(prefix,"",names(ret)))#strip the prefix
	ret
	
}

.readDNASeq<-function(fileName)
{
	seqfile <- system.file("extdata/parsed",fileName ,package="HIV.db")
	readDNAStringSet(seqfile)
}


loadFeatures<-function(ref="env",DNA=FALSE, refScale=NULL, genome="hxb2")	
{	
	HIV_db<-new.env(hash=TRUE, parent=emptyenv())
	
	tableName<-paste(genome,"Table", sep="")
	ret<-.readTblfromHIVdb(tblname=tableName)
	ret<-subset(ret,t_category!="internal")
	
	if(genome=="hxb2")
	{
	  ret1<-.readTblfromHIVdb("gpadditional")
	  ret<-rbind(ret,ret1)
    }
	
	#Selection of what is relevant considering the reference
	refFeature<-subset(ret,tolower(t_name)==tolower(ref))
	#ret<-subset(ret,t_name%in%ref)
	DNArefStart<-refFeature[["t_start"]]
	DNArefEnd<-refFeature[["t_end"]]
	ret<-subset(ret,t_start>=DNArefStart)
	ret<-subset(ret,t_end<=DNArefEnd)
	#Change the coordinates to AA relative to the ref
	if(!DNA)
	{
		ret<-subset(ret,t_frame==refFeature[["t_frame"]])
		ret[["t_start"]]<-sapply(ret[["t_start"]], function(x){ceiling((x-DNArefStart)/3)})
		ret[["t_end"]]<-sapply(ret[["t_end"]], function(x){ceiling((x-DNArefStart)/3)})
	}
	if(!is.null(refScale))
	{
		ret[["t_start"]]<-coord2ext(ret[["t_start"]],refScale)
		ret[["t_end"]]<-coord2ext(ret[["t_end"]], refScale)
	}
	
	assign("FeatureTable",ret,HIV_db)

	AbTableName<-paste("antibody", genome, sep="_")
	ret<-.readTblfromHIVdb(AbTableName)
	if(genome=="hxb2")
	{
	  ret<-rbind(ret,.readTblfromHIVdb("bindings"))
    }
	ret<-subset(ret,start>=DNArefStart)
	ret<-subset(ret,end<=DNArefEnd)
	#Change the coordinates to AA relative to the ref
	if(!DNA)
	{
		ret<-subset(ret,t_frame==refFeature[["t_frame"]])
		ret[["start"]]<-sapply(ret[["start"]], function(x){ceiling((x-DNArefStart)/3)})
		ret[["end"]]<-sapply(ret[["end"]], function(x){ceiling((x-DNArefStart)/3)})
	}	
	if(!is.null(refScale))
	{
		ret[["start"]]<-coord2ext(ret[["start"]], refScale)
		ret[["end"]]<-coord2ext(ret[["end"]], refScale)
	}
	
	ret$X<-rownames(ret)
	ret$MAb.Name<-sub("&","",ret$MAb.Name)

	
	assign("antibody",ret,HIV_db)


	AAfasta<-paste(genome,"_AA.fasta",sep="")
	DNAfasta<-paste(genome,"_DNA.fasta",sep="")
	assign("AA",.readAASeq(AAfasta, genome)[[tolower(ref)]],HIV_db)
	assign("DNA",subseq(.readDNASeq(DNAfasta), DNArefStart, DNArefEnd) ,HIV_db) #keep only the sequence coding for the ref
	assign("ref", ref, HIV_db) #keep track of the ref
	assign("genome", genome, HIV_db) #keep track of the genome

	HIV_db
}

lsCategory<-function(HIV_db)
{

	tbl<-get("FeatureTable",HIV_db)
	unique(tbl$t_category)
}
#when range is provided,return all the features that have intersections with the range
.getFeature<-function(HIV_db,category=NULL,name=NULL,start=NULL,end=NULL,frame=NULL,...)
{
	###if category is epitope then call getEpitope method to query antibodybinding table
	if(!is.null(category)&&tolower(category)=="epitope")
		return(getEpitope(HIV_db,name=name,start=start,end=end,frame=frame,...))
	####otherwise, query hxb2Table
	ret<-get("FeatureTable", HIV_db)
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
	ret<-HivFeature(ret,HIV_db)
	if(is.null(ret))
		message("No feature found, returning NULL.")
	ret
	
}
