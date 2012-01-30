\name{getEpitope-methods}
\docType{methods}
\alias{getEpitope-methods}
\alias{getEpitope,environment-method}
\alias{getEpitope,HivFeature-method}
\alias{getEpitope,list-method}
\alias{getEpitope}

\title{query Epitopes from HIV feature database}
\description{
 The methods provide the interfaces to query the HIV feature database for Epitopes 
}
\section{Methods}{
\describe{

\item{\code{signature(object = "environment")}}{
					Query the Epitopes directly from the data source stored in the environment.
					arguments are:
					\describe{
						\item{\code{name}:}{a character that passes the feature name }
						\item{\code{species}:}{a character that passes the species information}
						\item{\code{start,end}:}{integers that passe the coordinates of the annotation on the HIV genome}
						\item{\code{frame}:}{an integer passes the frame number}
					}
					\emph{Usage:}
		
				      \code{getEpitope(object,start=NULL,end=NULL,name=NULL,species=NULL,frame=NULL)}
}

\item{\code{signature(object = "HivFeature")}}{
				This method perform the same kind of query as the first one and use the \code{start,end} coordinates of the \code{HivFeature} object
				as the extra query condition.
				\emph{Usage:}
		
			     \code{getEpitope(object,start=NULL,end=NULL,...)}
			     \describe{
						
						\item{\code{start,end}:}{the coordinates of the HXB2 genome that intersect with coordinates of the HiveFeature object.}
						\item{\code{\dots}:}{additional arguments passed to getEpitope(object = "HivFeature") method}
					}
}			

\item{\code{signature(object = "list")}}{
				This method takes a list of \code{HivFeature} objects as the input
				\emph{Usage:}
		
			     \code{getEpitope(object,...)}
				\code{\dots} are the arguments passed to getEpitope(object = "HivFeature") method
}
}}

\value{
\code{signature(object = "environment")} and \code{signature(object = "HivFeature")} returns a Epitope object,
\code{signature(object = "list")} returns a list of Epitope objects
}
\author{
Mike Jiang,Raphael Gottardo
}

\seealso{
\code{\link{getFeature}}
}
\examples{
library(HIV.db)

## load Table from db
HIV_db<-loadFeatures()
getEpitope(HIV_db,start=4240,end=5096,species="mouse",name="4D6",frame=3)
}
\keyword{methods}


