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
	ret<-read.AAStringSet(seqfile)
	prefix<-paste(toupper(genome), "_AA_", sep="")
	names(ret)<-tolower(sub(prefix,"",names(ret)))#strip the prefix
	ret
	
}

.readDNASeq<-function(fileName)
{
	seqfile <- system.file("extdata/parsed",fileName ,package="HIV.db")
	read.DNAStringSet(seqfile)
}


loadFeatures<-function(ref="env",DNA=FALSE, refScale=NULL, genome="hxb2")	
{	
	HIV_db<-new.env(hash=TRUE, parent=emptyenv())
	
	tableName<-paste(genome,"Table", sep="")
	ret<-.readTblfromHIVdb(tblname=tableName)
#	ret<-.readTblfromHIVdb(tblname="hxb2Table")
	ret<-subset(ret,t_category!="internal")
	
	if(genome=="hxb2") ###
	{
	  ret1<-.readTblfromHIVdb("gpadditional")
	  ret<-rbind(ret,ret1)
    }
	
	#Selection of what is relevant considering the reference
	refFeature<-subset(ret,t_name==ref)
	#ret<-subset(ret,t_name%in%ref)
	ret<-subset(ret,t_start>=refFeature[["t_start"]])
	ret<-subset(ret,t_end<=refFeature[["t_end"]])
	#Change the coordinates to AA relative to the ref
	if(!DNA)
	{
		ret<-subset(ret,t_frame==refFeature[["t_frame"]])
		ret[["t_start"]]<-sapply(ret[["t_start"]], function(x){ceiling((x-refFeature[["t_start"]])/3)})
		ret[["t_end"]]<-sapply(ret[["t_end"]], function(x){ceiling((x-refFeature[["t_start"]])/3)})
	}
	if(!is.null(refScale))
	{
		ret[["t_start"]]<-coord2ext(ret[["t_start"]],refScale)
		ret[["t_end"]]<-coord2ext(ret[["t_end"]], refScale)
	}
	
	assign(tableName,ret,HIV_db) ###
#	assign("hxb2Table",ret,HIV_db)
	if(genome=="hxb2") ###
	{
	ret<-.readTblfromHIVdb("antibody")
	ret<-rbind(ret,.readTblfromHIVdb("bindings"))
	ret<-subset(ret,HXB2.start>=refFeature[["t_start"]])
	ret<-subset(ret,HXB2.end<=refFeature[["t_end"]])
	#Change the coordinates to AA relative to the ref
	if(!DNA)
	{
		ret<-subset(ret,t_frame==refFeature[["t_frame"]])
		ret[["HXB2.start"]]<-sapply(ret[["HXB2.start"]], function(x){ceiling((x-refFeature[["t_start"]])/3)})
		ret[["HXB2.end"]]<-sapply(ret[["HXB2.end"]], function(x){ceiling((x-refFeature[["t_start"]])/3)})
	}	
	if(!is.null(refScale))
	{
		ret[["HXB2.start"]]<-coord2ext(ret[["HXB2.start"]], refScale)
		ret[["HXB2.end"]]<-coord2ext(ret[["HXB2.end"]], refScale)
	}
	
	ret$X<-rownames(ret)
	ret$MAb.Name<-sub("&","",ret$MAb.Name)

	
	assign("antibody",ret,HIV_db)
    } ### 
	AAName<-paste(genome,"AA",sep="")
	DNAName<-paste(genome,"DNA",sep="")
	AAfasta<-paste(genome,"_AA.fasta",sep="")
	DNAfasta<-paste(genome,"_DNA.fasta",sep="")
	assign(AAName,.readAASeq(AAfasta, genome),HIV_db)
	assign(DNAName,.readAASeq(DNAfasta, genome),HIV_db)
	assign("genome", genome, HIV_db) #keep track of the genome


	HIV_db
}

lsCategory<-function(HIV_db)
{
	tbl<-get("hxb2Table",HIV_db)
	unique(tbl$t_category)
}
#when range is provided,return all the features that have intersections with the range
.getFeature<-function(HIV_db,category=NULL,name=NULL,start=NULL,end=NULL,frame=NULL,...)
{
	
	###if category is epitope then call getEpitope method to query antibodybinding table
	if(!is.null(category)&&tolower(category)=="epitope")
		return(getEpitope(HIV_db,name=name,start=start,end=end,frame=frame,...))
	####otherwise, query hxb2Table
	genome<-getGenome(HIV_db)
	tableName<-paste(genome,"Table", sep="") ###
	tbl<-get(tableName,HIV_db)
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




###
# Convert the coordinates of an object into the extended coordinate system
###
setGeneric("coord2ext", def=function(obj, refScale) standardGeneric("coord2ext"))

setMethod("coord2ext", signature=(obj="numeric"), function(obj, refScale)
		{
			extVec<-sapply(obj, function(x){
						min(
								if(x<0)
									x
								else if(x==0)
									which(refScale==1)
								else if(x>refScale[length(refScale)])
									which(refScale==refScale[length(refScale)])
								else if(length(which(refScale==x)))
									which(refScale==x)
								else
									NaN
						)})#sapply#function#min
			extVec<-extVec[!is.na(extVec)]
			return(extVec)
		})



setMethod("coord2ext", signature=(obj="RangedData"), function(obj, refScale)
		{
			if(start(obj)[[1]]==0) { start(obj)[[1]]=1 } #To avoid Inf values
			
			extStart<-coord2ext(start(obj),refScale)
			extEnd<-coord2ext(end(obj),refScale)
			#assign new start coordinates after end to avoid width<0 issues
			end(obj)<-extEnd	
			start(obj)<-extStart
			return(obj)
		})

