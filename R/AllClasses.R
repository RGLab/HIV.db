# TODO: Add comment
# 
# Author: mike
###############################################################################

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
  contains="HivFeature"
);
#constructor for Epitope from data.frame
Epitope<-function(object,HIV_db){
			if(nrow(object)>0)
			{
				ir<-IRanges(start=object[,"start"]
						,end=object[,"end"]
						)
				df<-DataFrame(FeatureID=object[,"X"]
						,name=object[,"MAb.Name"]
						,category="Epitope"
						,frame=object[,"t_frame"]
						,Epitope=object[,"Epitope"]
						,Species=object[,"Species"]
						,Subtype=object[,"Subtype"]
				)
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
