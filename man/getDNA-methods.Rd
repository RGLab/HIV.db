\name{getDNA-methods}
\docType{methods}
\alias{getDNA-methods}
\alias{getDNA,environment-method}
\alias{getDNA,HivFeature-method}
\alias{getDNA,list-method}
\alias{getDNA}

\title{ query DNA sequence}
\description{
  query DNA sequences of \code{\link{HivFeature}}s or \code{\link{Epitope}}s or 
 a certain range within HXB2 genome 
}
\section{Methods}{
\describe{

\item{\code{signature(object = "environment")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(object = "HivFeature")}}{
%%  ~~describe this method here~~
}

\item{\code{signature(object = "list")}}{
%%  ~~describe this method here~~
}
}}

\value{
a character of DNA sequence 
or a list of characters of DNA sequences
}

\seealso{
\code{\link{getAA}},\code{\link{getFeature}}
}
\author{
Mike Jiang,Raphael Gottardo
}


\examples{
library(HIV.db)

## load Table from db
HIV_db<-loadFeatures()
##########################
#query DNA sequence 
#######################

feature6<-getFeature(HIV_db, name="p7")

#by feature
getDNA(feature6[[1]])

#by feature list
feature2<-getFeature(HIV_db, name="pol")
getDNA(getChildren(feature2))

#by range
getDNA(HIV_db,start=1921,end=2085)

##if both range and feature are provided,the intersection of two ranges are used for query
getDNA(feature6[[1]],start=1921,end=2095)

#by Epitope
getDNA(getEpitope(HIV_db,start=4240,end=5096,species="mouse",name="4D6",frame=3))



}
\keyword{methods}
