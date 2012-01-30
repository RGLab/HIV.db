\name{loadFeatures}
\alias{loadFeatures}

\title{
Load the feature database
}
\description{
load all the data tables into an environment variable served as a data source for query
}
\usage{
loadFeatures()
}

\value{
  
  An environment that stores HIV feature database
  }
 
 \author{
Mike Jiang,Raphael Gottardo
}
  
\seealso{
\code{\linkS4class{HivFeature}},\code{\link{lsCategory}}
}




\examples{
library(HIV.db)

## load Table from db
HIV_db<-loadFeatures()
ls(HIV_db)
}

