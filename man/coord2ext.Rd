\name{coord2ext}
\alias{coord2ext}
\alias{coord2ext-methods}
\alias{coord2ext,numeric-method}
\alias{coord2ext,RangedData-method}

\title{
	Methods to convert a sequence coordinates into the extended coordinate system.
}

\description{
	Converts the given \code{"integer"} or \code{\linkS4class{RangedData}} object into the extended coordinates system using a reference scale.
	The positions to be converted must be inferior or equal to the length of the scale.
}

\section{Methods}{
	\describe{
		\item{\code{signature(obj = "numeric")}}{
			Convert a vector of \code{"numeric"} into the extended coordinates system.}

		\item{\code{signature(obj = "RangedData")}}{
			Convert the start and end  of a \code{\linkS4class{RangedData}} object into the extended coordinate system.}
 	}
 }

\author{
Renan Sauteraud
}

\keyword{methods}