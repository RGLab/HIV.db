\name{Epitope-class}
\Rdversion{1.1}
\docType{class}
\alias{Epitope-class}
\alias{Epitope}

\title{Class \code{"Epitope"}}
\description{
A class storing the HIV Epitope information.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Epitope", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{Epitope}:}{Object of class \code{"character"},stores the epitope sequence }
    \item{\code{Subtype}:}{Object of class \code{"character"} ,stoers the subtype information of the epitope }
    \item{\code{Species}:}{Object of class \code{"character"},store the specie that the epitope represents }
    \item{\code{FeatureID}:}{an integer that stores the unique annotation ID }
    \item{\code{name}:}{a character that stores the annotation name }
    \item{\code{category}:}{a character that stores the annotation category}
    \item{\code{start,end}:}{integers that store the coordinates of the annotation on the HIV genome}
    \item{\code{frame}:}{an integer storing the frame number of the annotation}
    \item{\code{parentID}:}{an integer storing the parent anotation ID}
    \item{\code{HIV_db}:}{an environment that stores all the data tables}
  }
}
\section{Extends}{
Class \code{"\linkS4class{HivFeature}"}, directly.
}
\author{
Mike Jiang,Raphael Gottardo
}


\seealso{
\code{\linkS4class{HivFeature}}
,\code{\link{getEpitope}},\code{\link{getAA}}

}

\keyword{classes}
