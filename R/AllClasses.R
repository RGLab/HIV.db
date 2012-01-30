# TODO: Add comment
# 
# Author: mike
###############################################################################



#setClass("HivFeature",  
#		representation(FeatureID = "integer",
#				name="character",
#				category="character",
#				start= "integer",
#				end= "integer",
#				frame= "integer",
#				parentID="integer",
#				HIV_db="environment"),
#		prototype(
#				FeatureID= integer(0),
#				name=character(0),
#				category=character(0),
#				start= integer(0),
#				end=integer(0),
#				parentID=integer(0),
#				HIV_db=new.env(hash=TRUE,parent=emptyenv())
#				)
#);



#setClass("Epitope",  
#		representation(FeatureID = "integer",
#				Epitope="character",
#				Subtype="character",
#				Species="character"),
#		prototype(
#				Epitope=character(0),
#				Subtype=character(0),
#				Species=character(0)
#				),
#		contains="HivFeature"
#);

setClass("HivFeature",  
		representation(
				HIV_db="environment"),
		prototype(
				HIV_db=new.env(hash=TRUE,parent=emptyenv())
		)
		,contains="RangedData"
);

#constructor for data.frame
HivFeature<-function(object,HIV_db){	
			
#		lapply(1:nrow(object),function(i){
#			new("HivFeature",FeatureID=object[i,"t_ID"],name=object[i,"t_name"],
#					category=object[i,"t_category"],start=object[i,"t_start"],
#					end=object[i,"t_end"],parentID=object[i,"t_parentID"],frame=object[i,"t_frame"],HIV_db=HIV_db)
#		})
#			browser()
			if(nrow(object)>0)
			{
				ir<-IRanges(start=object[,"t_start"]
						,end=object[,"t_end"]
#						,name=object[,"t_name"]
						)
				df<-DataFrame(FeatureID=object[,"t_ID"]
						,name=object[,"t_name"]
						,category=object[,"t_category"]
						,parentID=object[,"t_parentID"]
						,frame=object[,"t_frame"])
#				browser()
				rownames(df) <- names(ir) ## ensure these are identical
				#			N <- sum(elementLengths(ir))
				space <- Rle(factor("1"))
				irs<-split(ir,space)
				dfs<-split(df, space)
				new("HivFeature",ranges=irs
						,values=dfs
						,HIV_db=HIV_db
				)
			}else
			{
				NULL
				
			}
			
			
			
		}
		
setClass("Epitope",  
#		representation(FeatureID = "integer",
#				Epitope="character",
#				Subtype="character",
#				Species="character"),
#		prototype(
#				Epitope=character(0),
#				Subtype=character(0),
#				Species=character(0)
#				),
		,contains="HivFeature"
);
#constructor for Epitope from data.frame
Epitope<-function(object,HIV_db){	
#			browser()
#			lapply(1:nrow(object),function(i){
#						new("Epitope",FeatureID=as.integer(object[i,"X"]),name=object[i,"MAb.Name"],
#								category="Epitope",start=object[i,"HXB2.start"],
#								end=object[i,"HXB2.end"],frame=object[i,"t_frame"],HIV_db=HIV_db,
#								Epitope=object[i,"Epitope"],Species=object[i,"Species"],Subtype=object[i,"Subtype"])
#					})
			if(nrow(object)>0)
			{
				ir<-IRanges(start=object[,"HXB2.start"]
						,end=object[,"HXB2.end"]
#						,name=object[,"MAb.Name"]
						)
				df<-DataFrame(FeatureID=object[,"X"]
#						,MAb.Name=object[,"MAb.Name"]
							,name=object[,"MAb.Name"]
						,category="Epitope"
						,frame=object[,"t_frame"]
						,Epitope=object[,"Epitope"]
						,Species=object[,"Species"]
						,Subtype=object[,"Subtype"]
				)
#				browser()
				rownames(df) <- names(ir) ## ensure these are identical
				#			N <- sum(elementLengths(ir))
				space <- Rle(factor("1"))
				irs<-split(ir,space)
				dfs<-split(df, space)
				new("Epitope",ranges=irs
						,values=dfs
						,HIV_db=HIV_db
				)
			}else
			{
				NULL
				
			}
		}