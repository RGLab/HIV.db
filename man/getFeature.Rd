\name{getFeature-methods}
\docType{methods}
\alias{getFeature}
\alias{getFeature,environment-method}

\title{
query the annotation feature from the HIV feature database
}
\description{
The method provides the simply query interface to the HIV feature database  
}


\section{Methods}{
\describe{

	\item{\code{signature(object = "environment")}}{

		\code{object} is an environment that stores the HIV feature database,
		\code{\dots}:} are the additional arguments:
						  \describe{
						  
								  	\item{\code{name}:}{a character that stores the annotation name }
									\item{\code{category}:}{a character that stores the annotation category}
									\item{\code{start,end}:}{integers that store the coordinates of the annotation on the HIV genome}
									\item{\code{frame}:}{an integer storing the frame number of the annotation}
									}
						}
	\emph{Usage:}

  \code{getFeature(object,start=NULL,end=NULL,name=NULL,species=NULL,frame=NULL)}
 }

\value{
  
  A HivFeature object or a list of HivFeature objects
  
}
\author{
Mike Jiang,Raphael Gottardo
}

\seealso{
\code{\link{getEpitope}}, \code{\linkS4class{HivFeature}}}

\examples{
library(HIV.db)

## load Table from db
HIV_db<-loadFeatures()
## query by name or category
getFeature(HIV_db, category=c("gene"))
getFeature(HIV_db, category=c("gene"),frame=3)
getFeature(HIV_db, name=c("gp120","gp41"))
getFeature(HIV_db, category="loop")

## query by hxb2 coordinates
getFeature(HIV_db,start=7100,end=7200)
getFeature(HIV_db,start=6225,end=8795)

}

