# TODO: Add comment
# 
# Author: mike
###############################################################################



setMethod("show",
		signature=signature(object="HivFeature"),
		definition=function(object){


			values(object)<-values(object)[,!colnames(object)%in%c("FeatureID","parentID")]
			
			
			callNextMethod()
#			cat("HIV Feature:",object@name,"\n")
#			cat("type:",object@category,"\n")
#			cat("frame:",object@frame,"\n")
#			cat("Coordinates:",paste(getHXB2Coordinates(object)[1],getHXB2Coordinates(object)[2],sep="-"), "\n")

		})

#to be deprecated 		
setMethod("getHXB2Coordinates",
		signature=signature(object="HivFeature"),
		definition=function(object){		
#		c(object@start,object@end)
#			browser()
			cbind(start(object),end(object))
			
	})

setMethod("start",
		signature=signature(x="HivFeature"),
		definition=function(x){		
#		c(object@start,object@end)
#			browser()
			start(ranges(x)[[1]])
		})

setMethod("end",
		signature=signature(x="HivFeature"),
		definition=function(x){		
#		c(object@start,object@end)
#			browser()
			end(ranges(x)[[1]])
		})

#to be deprecated 		
setMethod("getFrame",
		signature=signature(object="HivFeature"),
		definition=function(object){		
#			object@frame
#			browser()
			values(object)[[1]][,"frame"]
		})
#to be deprecated 		
setMethod("getFrame",
		signature=signature(object="list"),
		definition=function(object){		
			unlist(lapply(object,getFrame))
		})

#to be deprecated 		
setMethod("getName",
		signature=signature(object="HivFeature"),
		definition=function(object){		
#			object@name
#			browser()
#			names(ranges(object)[[1]])
			values(object)[[1]][,"name"]
		})

#setMethod("getName",
#		signature=signature(object="Epitope"),
#		definition=function(object){		
##			object@name
##			browser()
##			names(ranges(object)[[1]])
#			values(object)[[1]][,"MAb.Name"]
#			
#		})

#to be deprecated 		
setMethod("getName",
		signature=signature(object="list"),
		definition=function(object){		
			unlist(lapply(object,getName))
		})

setGeneric("FeatureID", 
		function(object)
			standardGeneric("FeatureID"))
setMethod("FeatureID",
		signature=signature(object="HivFeature"),
		definition=function(object){		
#			object@name
#			browser()
			values(object)[[1]][,"FeatureID"]
		})

setGeneric("parentID", 
		function(object)
			standardGeneric("parentID"))
setMethod("parentID",
		signature=signature(object="HivFeature"),
		definition=function(object){		
#			object@name
#			browser()
			values(object)[[1]][,"parentID"]
		})


setMethod("getHIVdb",
		signature=signature(object="HivFeature"),
		definition=function(object){		
			object@HIV_db
		})

setMethod("getAnnotationTable",
		signature=signature(object="HivFeature"),
		definition=function(object){
			genome<-getGenome(object)
			tableName<-paste(genome,"Table", sep="")
			get(tableName,getHIVdb(object))
		})

setMethod("getAntibodyTable",
		signature=signature(object="HivFeature"),
		definition=function(object){		
			get("antibody",getHIVdb(object))
		})


#setMethod("getChildren",
#		signature=signature(object="list"),
#		definition=function(object,...){		
##			browser()
#			lapply(object,getChildren,...)
#			
#		})

setMethod("getChildren",
		signature=signature(object="HivFeature"),
		definition=function(object,...){		


			category<-list(...)$category
			recursive=list(...)$recursive
			ret<-getAnnotationTable(object)
#			
			if(!is.null(category))
				ret<-subset(ret,t_category%in%category)
			
			if(is.null(recursive))
				recursive<-TRUE
						
			
			ret<-lapply(seq_along(1:nrow(object)),function(i){
								curF<-object[i,]
								if(recursive)
									subset(ret,t_start>=start(curF)&t_end<=end(curF)&t_ID!=FeatureID(curF))
								else
									subset(ret,t_parentID==FeatureID(curF))
								
							})
							
			ret<-do.call(rbind,ret)
			ret<-unique(ret)
			
			ret<-HivFeature(ret,getHIVdb(object))
			ret
	
	})


	
#setMethod("getParent",
#		signature=signature(object="list"),
#		definition=function(object,...){		
#			
#			lapply(object,getParent,...)
#			
#		})	
	
setMethod("getParent",
		signature=signature(object="HivFeature"),
		definition=function(object, ...){
#			db=list(...)$db

			db<-getAnnotationTable(object)
			recursive=list(...)$recursive
			if(is.null(recursive))
				recursive<-FALSE
			if(recursive)
				ret<-subset(db,t_start<=start(object)&t_end>=end(object)&t_ID!=FeatureID(object))
			else
			{
				parentID<-subset(db,t_ID==FeatureID(object))$t_parentID
				ret<-subset(db,t_ID==parentID)
			}
			ret<-HivFeature(ret,getHIVdb(object))
			
#			if(length(ret)==1)
#				ret[[1]]
#			else
#				ret
			
			ret
		})

setMethod("getEpitope",
		signature=signature(object="list"),
		definition=function(object,...){
			
#			do.call(rbind,lapply(object,getEpitope,...))
			ret<-lapply(object,getEpitope,...)
#			names(ret)<-getName(object)
			unlist(ret)
		})



setMethod("getEpitope",
		signature=signature(object="HivFeature"),
		definition=function(object,start=NULL,end=NULL,...){
				#filter by feature coordinates
				hxb2_range<-getHXB2Coordinates(object)
#				ret<-subset(ret,HXB2.end>=hxb2_range[1]&HXB2.start<=hxb2_range[2])
				start<-max(hxb2_range[1],start)
				end<-min(hxb2_range[1,][2],end)
#				name<-list(...)$name
#				species<-list(...)$species
#				frame<-list(...)$frame
#			browser()
				getEpitope(getHIVdb(object),start=start,end=end,...)
		})

setMethod("getEpitope",
		signature=signature(object="environment"),
		definition=function(object,start=NULL,end=NULL,name=NULL,species=NULL,frame=NULL){
							
		#filter by start and end position
		ret<-get("antibody",object)
#		start=list(...)$start
#		end=list(...)$end
		if(!is.null(start))
			ret<-subset(ret,end>=start)
		if(!is.null(end))
			ret<-subset(ret,start<=end)
	
#		name<-list(...)$name
		if(!is.null(name))
			ret<-subset(ret,MAb.Name%in%name)
		
		#filter by species
#		species<-list(...)$species
		if(!is.null(species))
			ret<-subset(ret,Species%in%species)
#		frame<-list(...)$frame
		if(!is.null(frame))
			ret<-subset(ret,t_frame%in%frame)
#		browser()
		if(nrow(ret)>0)
			Epitope(ret,object)
		else
			NULL
		
})

setMethod("getFeature",
		signature=signature(object="environment"),
		definition=function(object,...){
			
			.getFeature(object,...)			
		})


setMethod("getDNA",
		signature=signature(object="environment"),
		definition=function(object, ...){
			start=list(...)$start
			end=list(...)$end
			genome<-getGenome(object)
			DNAName<-paste(genome,"DNA",sep="")
			subseq(get(DNAName,object)[[1]],start,end)
#			ret<-substr(get("hxb2DNA",object),start,end)
#			names(ret)<-""
#			ret
			
		})




setMethod("getDNA",
		signature=signature(object="HivFeature"),
		definition=function(object,...){
	DNA_list<-lapply(seq_along(1:nrow(object)),function(i){
						curF<-object[i,]
						
			hxb2_range<-getHXB2Coordinates(curF)		
			start=list(...)$start
			end=list(...)$end
			start<-max(start,hxb2_range[1])
			end<-min(end,hxb2_range[2])
			
			getDNA(getHIVdb(object),start=start,end=end)
			
		})
	names(DNA_list)<-getName(object)
	DNA_list
})
#setMethod("getDNA",
#		signature=signature(object="list"),
#		definition=function(object,...){
#			
#			lapply(object,getDNA,...)
#			
#		})

setMethod("isRoot",
		signature=signature(object="HivFeature"),
		definition=function(object){
			ifelse(parentID(object)==0,TRUE,FALSE)
		})

setMethod("getAA",
		signature=signature(object="environment"),
		definition=function(object, ...){
			
			name=tolower(list(...)$name)
			##currently only query the Full AA sequence 
#			subseq(get("hxb2AA",object)[2],start,end)
			
#			ret<-substr(get("hxb2AA",object)[[name]],start,end)
#			names(ret)<-""
#			browser()
			genome<-getGenome(object)
			AAName<-paste(genome, "AA", sep="")
			seqSet<-get(AAName,object)
			if(name%in%names(seqSet))
#				return(as.character(get("hxb2AA",object)[[name]]))
				get(AAName,object)[[name]]
			else
				return(NULL)
		})


setMethod("getAA",
		signature=signature(object="HivFeature"),
		definition=function(object, ...){
			
	AA_list<-lapply(seq_along(1:nrow(object)),function(i){
				curF<-object[i,]
				
				###search for top node for AA sequence
				curNode<-curF
	#			browser()
				while(!isRoot(curNode))
				{
					curNode<-getParent(curNode)	
				}
				
				offset<-getHXB2Coordinates(curNode)[1]
				hxb2_range<-getHXB2Coordinates(curF)
				AA_range<-(hxb2_range-offset)/3+1
				
	#			##ask for posType when it is not provided
	#			start<-list(...)$start
	#			end<-list(...)$end
	#			posType=list(...)$posType
	#			if(!is.null(start)&&is.null(posType)){
	#				message("Please choose the type of position:\n");
	#				posType<-menu(c("relative AA position","absolute HXB2 DNA coordinates"),graphics=FALSE)
	#				
	#			}
	#			
	#			if(!is.null(posType)&&posType==2)
	#			{
	#				start=(start-offset)/3+1
	#				end=(end-offset)/3+1
	#			}
	#			start<-max(start,AA_range[1])
	#			end<-min(end,AA_range[2])
				start<-AA_range[1]
				end<-AA_range[2]
	#			getAA(getHIVdb(object),start=start,end=end,name=getName(curNode))
#				browser()
	
				parentSeq<-getAA(getHIVdb(curF),name=getName(curNode))
				if(!is.null(parentSeq))
					substr(parentSeq,start,end)
				else
					""
			})
#	browser()
	names(AA_list)<-getName(object)
	AA_list
})

setMethod("getAA",
		signature=signature(object="Epitope"),
		definition=function(object,...){
			
			AA_list<-lapply(seq_along(1:nrow(object)),function(i){
						curF<-object[i,]
						AAString(values(curF)[[1]][,"Epitope"])
			
					})
			names(AA_list)<-getName(object)
			AA_list			
})
#setMethod("getAA",
#		signature=signature(object="list"),
#		definition=function(object,...){
#			
#			lapply(object,getAA,...)
#		})

#accessor to the genome elt of an environment object
setMethod("getGenome",
		signature=signature(object="environment"),
		definition=function(object, ...){
			return(object$genome)
		})

#for an HivFeature
setMethod("getGenome",
		signature=signature(object="HivFeature"),
		definition=function(object, ...){
			return(getGenome(getHIVdb(object)))
		})


# Input: A DataFrame like pep_hxb2 or pep_mac239
#        The only requirements are petides as rownames and a `clade` column
# Output: matrix of logical with peptides as rownames and clades as colnames
setMethod("clade",
		signature=signature(object="RangedData"),
		definition=function(object)
{
	cladeList<-unique(unlist(strsplit(levels(as.factor(object$clade)),",")))
	len<-nrow(object)
	retMatrix<-c()
	s2<-system.time(
			{
				
	for(pepIdx in 1:len)
	{
		pepClades<-unlist(strsplit(object[pepIdx,]$clade, split=","))
		tmpList<-lapply(cladeList, function(x)
				{
					if(x %in% pepClades)
						TRUE
					else
						FALSE
				})
		retMatrix<-c(retMatrix, tmpList)
	}
	})#s2

	print(c(s2))

	dim(retMatrix)<-c(length(cladeList), len)
	retMatrix<-t(retMatrix)
	rownames(retMatrix)<-rownames(object)
	colnames(retMatrix)<-cladeList

	return(retMatrix)
})