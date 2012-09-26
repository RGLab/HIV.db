\name{getAA-methods}
\docType{methods}
\alias{getAA-methods}
\alias{getAA,environment-method}
\alias{getAA,Epitope-method}
\alias{getAA,HivFeature-method}
\alias{getAA,list-method}
\alias{getAA}

\title{query amino acids sequence}
\description{
 query amino acids sequences of \code{\link{HivFeature}}s or \code{\link{Epitope}}s 
}
\section{Methods}{
\describe{

\item{\code{signature(object = "environment")}}{

}

\item{\code{signature(object = "Epitope")}}{
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
a character of amino acids sequence 
or a list of characters of amino acids sequences
}
\seealso{
\code{\link{getDNA}},\code{\link{getFeature}}

}

\author{
Mike Jiang,Raphael Gottardo
}

\examples{
library(HIV.db)

## load Table from db
HIV_db<-loadFeatures(ref="env")
##########################
#query amino acid sequence 
#######################
#by environment
getAA(HIV_db)

#by feature
feature1<-getFeature(HIV_db)
getAA(feature1)
#by name
getAA(feature1, name="gp120")

#by Epitope
getAA(getEpitope(HIV_db, start=100, end=200))

}
\keyword{methods}

