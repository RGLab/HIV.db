\name{HivFeature-class}
\Rdversion{1.1}
\docType{class}
\alias{HivFeature-class}
\alias{HivFeature}
\alias{getFrame,HivFeature-method}
\alias{getHXB2Coordinates,HivFeature-method}
\alias{getName,HivFeature-method}
\alias{show,HivFeature-method}
\alias{getFrame}
\alias{getFrame,list-method}
\alias{getHXB2Coordinates}
\alias{getName}
\alias{getName,list-method}

\alias{show}

\title{Class \code{"HivFeature"}}
\description{
The class that stores HIV features.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("HivFeature", ...)}.
 
}
\section{Slots}{
  \describe{
    \item{\code{FeatureID}:}{an integer that stores the unique annotation ID }
    \item{\code{name}:}{a character that stores the annotation name }
    \item{\code{category}:}{a character that stores the annotation category}
    \item{\code{start,end}:}{integers that store the coordinates of the annotation on the HIV genome}
    \item{\code{frame}:}{an integer storing the frame number of the annotation}
    \item{\code{parentID}:}{an integer storing the parent annotation ID}
    \item{\code{HIV_db}:}{an environment that stores HIV feature database}
  }
}
\section{Methods}{
  \describe{
    
    \item{getFrame}{\code{signature(object = "HivFeature")}: ... }
    \item{getHXB2Coordinates}{\code{signature(object = "HivFeature")}: ... }
    \item{getName}{\code{signature(object = "HivFeature")}: ... }
    \item{show}{\code{signature(object = "HivFeature")}: ... }
	 }
}

\author{
Mike Jiang,Raphael Gottardo
}

\seealso{
\code{\linkS4class{Epitope}},\code{\link{getAA}},\code{\link{getDNA}},\code{\link{getFeature}},\code{\link{getChildren}},
\code{\link{getParent}} 
}
\examples{
library(HIV.db)

## load Table from db with default parameters (ref="env").
HIV_db<-loadFeatures()

feature1<-getFeature(HIV_db)
show(feature1)

getHXB2Coordinates(feature1)
getName(feature1)
getFrame(feature1)
}
\keyword{classes}
