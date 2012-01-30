\name{getChildren-methods}
\docType{methods}
\alias{getChildren-methods}
\alias{getChildren,HivFeature-method}
\alias{getChildren,list-method}
\alias{getParent-methods}
\alias{getParent,HivFeature-method}
\alias{getParent,list-method}
\alias{getChildren}
\alias{getParent}
\title{query the children or parent HivFeature}
\description{
The children/parent relations are decided by the range of HXB2 coordinates of HivFeatures.
These methods provide the interface to query the HIV database based on the hierarchical structure of HivFeatures.  
}
\section{Methods}{
\describe{

	
\item{\code{signature(object = "HivFeature")}}{
 					A logical scalar \code{recursive} indicates if all the descendants or all the ancestors are returned besides the immediate ones. 
 					\code{category} can be used in \code{getChildren} to further narrow down the range the HivFeatures to be searched}
 					\emph{Usage:}
						
			      \code{ getChildren(feature1,recursive=FALSE,category=NULL)}
			      \code{ getParent(feature1,recursive=FALSE)}


\item{\code{signature(object = "list")}}{}
}
}



\value{
  
  A HivFeature object of a list of HivFeature objects
  
}
\seealso{
\code{\link{getParent}},\code{\link{getFeature}}}
\author{
Mike Jiang,Raphael Gottardo
}

\examples{
library(HIV.db)

## load Table from db
HIV_db<-loadFeatures()

## query by parent
feature1<-getFeature(HIV_db, name="gag")

getChildren(feature1)

## query by children
feature3<-getFeature(HIV_db, name="p24")
getParent(feature3[[1]])

feature4<-getFeature(HIV_db, name="p15")
getParent(feature4)

feature5<-getFeature(HIV_db, name="V1")
getParent(feature5)


#get all the Descendants or all the Ancestors including the immediate
getParent(feature5,recursive=TRUE)
getChildren(feature1,recursive=TRUE)


}

\keyword{methods}

