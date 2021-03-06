% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotPCPInteractive.R
\name{plotPCPInteractive}
\alias{plotPCPInteractive}
\title{Highlight parallel coordinate plot lines inside selection box}
\usage{
plotPCPInteractive(pcpDat, option = "deleteInteger")
}
\arguments{
\item{pcpDat}{the data frame that contains the parallel coordinate plot values}

\item{option}{the interactivity option ("deleteInteger", "delete", "highlight"); default ("deleteInteger")}
}
\description{
Highlight parallel coordinate plot lines inside selection box
}
\examples{
set.seed(3)
f = function(){1.3*rnorm(500)}
pcpDat = data.frame(ID = paste0("ID", 1:500), A.1=f(), A.2=f(), A.3=f(), B.1=f(), B.2=f(), B.3=f())
pcpDat$ID = as.character(pcpDat$ID)
plotPCPInteractive(pcpDat = pcpDat)
}
