\name{loadFeatures}
\alias{loadFeatures}

\title{
Load the feature database
}
\description{
load all the data tables into an environment variable served as a data source for query
}
\usage{
loadFeatures(ref="env",DNA=FALSE, refScale=NULL, genome="hxb2")
}

\arguments{
	\item{ref}{A \code{character}. The name of the reference feature, for subsetting purpose. The loaded feature will all be children of the reference: they belong to the same frame and are located within the ranges of the reference.}
	\item{DNA}{A \code{logical}. If set to TRUE, the coordinates are nucleotide position relative to the start of the selected genome.}
	\item{refScale}{A \code{numeric} vector. Used to convert the table into extended coordinate system. See \code{\link[Pviz]{readAlign}} for information on coordinates system and how to create refScale.}
	\item{genome}{A \code{character}. The name of the genome. Limited to hxb2 and mac239 at the moment.}
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

