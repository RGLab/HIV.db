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
HIV_env<-loadFeatures(ref="env")
##########################
#query DNA sequence 
#######################

feature1<-getFeature(HIV_env)

#by feature
getDNA(feature1)

#by feature list
getDNA(feature1, name=c("gp41","gp120"))

#by Epitope
getDNA(getEpitope(HIV_env,name="2F5"))



}
\keyword{methods}
