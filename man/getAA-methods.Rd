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
HIV_db<-loadFeatures()
##########################
#query amino acid sequence 
#######################
#by name
getAA(HIV_db,name="pol")
getAA(HIV_db,name="gag")

feature6<-getFeature(HIV_db, name="p7")

#by feature
getAA(feature6[[1]])

#by feature list
feature2<-getFeature(HIV_db, name="pol")
getAA(getChildren(feature2))

#by Epitope
getAA(getEpitope(HIV_db,start=4240,end=5096,species="mouse",name="4D6",frame=3))

}
\keyword{methods}

