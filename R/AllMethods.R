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
			get("FeatureTable",getHIVdb(object))
		})

setMethod("getAntibodyTable",
		signature=signature(object="HivFeature"),
		definition=function(object){		
			get("antibody",getHIVdb(object))
		})

setMethod("getChildren",
		signature=signature(object="HivFeature"),
		definition=function(object,...){		


			category<-list(...)$category
			recursive=list(...)$recursive
			ret<-getAnnotationTable(object)
			
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
				getEpitope(getHIVdb(object),start=start,end=end,...)
		})

setMethod("getEpitope",
		signature=signature(object="environment"),
		definition=function(object,start=NULL,end=NULL,name=NULL,species=NULL,frame=NULL){
							
		#filter by start and end position
		ret<-get("antibody",object)

		queryStart<-start
		queryEnd<-end
		if(!is.null(queryStart))
			ret<-subset(ret,end>=queryStart)
		if(!is.null(queryEnd))
			ret<-subset(ret,start<=queryEnd)
	
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
		definition=function(object, name=NULL){
			if(is.null(name))
			{
				return(object$DNA)
			}
			else
			{
				idxs<-which(object$FeatureTable$t_name==name)
				names<-object$FeatureTable[idxs,]$t_name
				starts<-object$FeatureTable[idxs,]$t_start
				ends<-object$FeatureTable[idxs,]$t_end
				seqs<-c()
				for(i in 1:length(idxs))
					seqs<-c(seqs, subseq(object$DNA, starts[i]*3, ends[i]*3))
				names(seqs)<-names
				return(seqs)
			}
#			start=list(...)$start
#			end=list(...)$end
#			if(is.null(start))
#				start<-1
#			if(is.null(end))
#				end<-length(get("DNA",object)[[1]])
#			subseq(get("DNA",object)[[1]],start,end)
##			ret<-substr(get("hxb2DNA",object),start,end)
##			names(ret)<-""
##			ret
			
		})




setMethod("getDNA",
		signature=signature(object="HivFeature"),
		definition=function(object, name=NULL){
#		DNA_list<-lapply(seq_along(1:nrow(object)),function(i){
#						curF<-object[i,]
#						
#			hxb2_range<-getHXB2Coordinates(curF)		
#			start=list(...)$start
#			end=list(...)$end
#			start<-max(start,hxb2_range[1])
#			end<-min(end,hxb2_range[2])
#			
#			getDNA(getHIVdb(object),start=start,end=end)
#			
#		})
#	names(DNA_list)<-getName(object)
#	DNA_list
			getDNA(getHIVdb(object), name)
})


setMethod("isRoot",
		signature=signature(object="HivFeature"),
		definition=function(object){
			ifelse(parentID(object)==0,TRUE,FALSE)
		})

#setMethod("getAA",
#		signature=signature(object="environment"),
#		definition=function(object, ...){
#			
#			name<-list(...)$name
#			if(is.null(name))
#			{
#				name<-object$ref
#			}
#			name<-tolower(name)
#			##currently only query the Full AA sequence 
##			subseq(get("hxb2AA",object)[2],start,end)
#			
##			ret<-substr(get("hxb2AA",object)[[name]],start,end)
##			names(ret)<-""
##			browser()
#			genome<-getGenome(object)
#			AAName<-paste(genome, "AA", sep="")
#			seqSet<-get(AAName,object)
#			if(name%in%names(seqSet))
##				return(as.character(get("hxb2AA",object)[[name]]))
#				get(AAName,object)[[name]]
#			else
#				return(NULL)
#		})

setMethod("getAA",
		signature=signature(object="environment"),
		definition=function(object, name=NULL)
		{
			if(is.null(name))
			{
				return(object$AA)
			}
			else
			{
				idxs<-which(object$FeatureTable$t_name==name)
				names<-object$FeatureTable[idxs,]$t_name
				starts<-object$FeatureTable[idxs,]$t_start
				ends<-object$FeatureTable[idxs,]$t_end
				seqs<-c()
				for(i in 1:length(idxs))
				  seqs<-c(seqs, subseq(object$AA, starts[i], ends[i]))
			  	names(seqs)<-names
				return(seqs)
			}
		})

setMethod("getAA",
		signature=signature(object="HivFeature"),
		definition=function(object, name=NULL)
		{
			getAA(getHIVdb(object), name)
		})
			




#setMethod("getAA",
#		signature=signature(object="HivFeature"),
#		definition=function(object, ...){
#			name<-list(...)$name
#
#	AA_list<-lapply(seq_along(1:nrow(object)),function(i){
#				curF<-object[i,]
#				
#				###search for top node for AA sequence
#				curNode<-curF
#				while(!isRoot(curNode))
#				{
#					curNode<-getParent(curNode)	
#				}
#				
#				offset<-getHXB2Coordinates(curNode)[1]
#				hxb2_range<-getHXB2Coordinates(curF)
#				AA_range<-(hxb2_range-offset)/3+1
#				
#	#			##ask for posType when it is not provided
#	#			start<-list(...)$start
#	#			end<-list(...)$end
#	#			posType=list(...)$posType
#	#			if(!is.null(start)&&is.null(posType)){
#	#				message("Please choose the type of position:\n");
#	#				posType<-menu(c("relative AA position","absolute HXB2 DNA coordinates"),graphics=FALSE)
#	#				
#	#			}
#	#			
#	#			if(!is.null(posType)&&posType==2)
#	#			{
#	#				start=(start-offset)/3+1
#	#				end=(end-offset)/3+1
#	#			}
#	#			start<-max(start,AA_range[1])
#	#			end<-min(end,AA_range[2])
#				start<-AA_range[1]
#				end<-AA_range[2]
#	
#				parentSeq<-getAA(getHIVdb(curF),name=getName(curNode))
#				if(!is.null(parentSeq))
#					substr(parentSeq,start,end)
#				else
#					""
#			})
#	names(AA_list)<-getName(object)
#	AA_list
#})

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

#accessor to the ref elt of an envir
setMethod("getRef",
		signature=signature(object="environment"),
		definition=function(object){
			return(object$ref)
		})

#for an HivFeature
setMethod("getRef",
		signature=signature(object="environment"),
		definition=function(object){
			return(getRef(getHIVdb(object)))
		})

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

setGeneric("coord2ext", def=function(obj, refScale){ standardGeneric("coord2ext") })

# Convert the coordinates of an object into the extended coordinate given a scale
# Input: An object and a reference scale
# Output: an object of the same type with coordinates in extended system
setMethod("coord2ext", signature=(obj="numeric"), function(obj, refScale){   
  extVec<-sapply(obj, function(x){
    min(
      if(x<0){
        x   
      }else if(x==0){
        which(refScale==1)
      }else if(x>refScale[length(refScale)]){
        which(refScale==refScale[length(refScale)])
      }else if(length(which(refScale==x))){
        which(refScale==x)
      }else{            
        NaN 
      })})
  extVec<-extVec[!is.na(extVec)]
  return(extVec)
})              

setMethod("coord2ext", signature=(obj="RangedData"), function(obj, refScale){
  if(start(obj)[[1]]==0) { start(obj)[[1]]=1 } #To avoid Inf values
          
  extStart<-coord2ext(start(obj),refScale)
  extEnd<-coord2ext(end(obj),refScale)
  #assign new start coordinates after end to avoid width<0 issues
  end(obj)<-extEnd
  start(obj)<-extStart
  return(obj)
})

setMethod("coord2ext", signature=(obj="IRanges"), function(obj, refScale){
  if(start(obj)[[1]]==0) { start(obj)[[1]]=1 } #To avoid Inf values

  extStart<-coord2ext(start(obj),refScale)
  extEnd<-coord2ext(end(obj),refScale)
  #assign new start coordinates after end to avoid width<0 issues
  end(obj)<-extEnd
  start(obj)<-extStart
  return(obj)
})

