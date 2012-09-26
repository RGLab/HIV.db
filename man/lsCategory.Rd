\name{lsCategory}
\alias{lsCategory}

\title{
list all the availabe category in HIV.db
}

\description{
list all the available categories of HIV features in HIV.db 

}

\usage{
	lsCategory(HIV_db)
}

\arguments{
	\item{HIV_db}{an environment that stores the HIV feature database}
}

\value{
  
  A character vector that list all the category names
  
}
\author{
Mike Jiang,Raphael Gottardo
}

\seealso{
\code{\linkS4class{HivFeature}}}



\examples{
library(HIV.db)

## load Table from db
HIV_db<-loadFeatures()

###list available category
lsCategory(HIV_db)
}

\keyword{functions}

