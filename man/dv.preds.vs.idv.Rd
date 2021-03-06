% Xpose 4
% An R-based population pharmacokinetic/
% pharmacodynamic model building aid for NONMEM.
% Copyright (C) 1998-2004 E. Niclas Jonsson and Mats Karlsson.
% Copyright (C) 2005-2008 Andrew C. Hooker, Justin J. Wilkins, 
% Mats O. Karlsson and E. Niclas Jonsson.
% Copyright (C) 2009-2010 Andrew C. Hooker, Mats O. Karlsson and 
% E. Niclas Jonsson.

% This file is a part of Xpose 4.
% Xpose 4 is free software; you can redistribute it and/or
% modify it under the terms of the GNU Lesser General Public License
% as published by the Free Software Foundation, either version 3
% of the License, or (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.

% You should have received a copy of the GNU Lesser General Public License
% along with this program.  A copy can be cound in the R installation
% directory under \share\licenses. If not, see http://www.gnu.org/licenses/.

\name{dv.preds.vs.idv}
\alias{dv.preds.vs.idv}
\title{Observations (DV), individual predictions (IPRED) and population 
  predictions (IPRED) plotted against the independent variable (IDV), 
  for Xpose 4}
\description{
  This is a compound plot consisting of plots of observations (DV),
  individual predictions (IPRED), and population predictions (PRED) 
  against the independent variable (IDV), a specific function in Xpose 4. 
  It is a wrapper encapsulating arguments to the \code{xpose.plot.default} 
  function. 
}
\usage{
dv.preds.vs.idv(object,
           ylb  = "Observations/Predictions",
           layout=c(3,1),
           smooth=TRUE,
           scales=list(),
           ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An xpose.data object. }
  \item{ylb}{A string giving the label for the y-axis. \code{NULL} if none.}
  \item{layout}{A list controlling the number of columns and rows in a compound
    plot. The default is 2 columns and 1 row.}  
  \item{smooth}{Logical value indicating whether an x-y smooth should be superimposed.
    The default is TRUE.}
  \item{scales}{A list to be used for the \code{scales} argument in
    \code{xyplot}.}
  \item{\dots}{Other arguments passed to \code{link[xpose4generic]{xpose.plot.default}}.}
}
\details{
  A wide array of extra options controlling \code{xyplots} are 
  available. See \code{\link[xpose4generic]{xpose.plot.default}} and 
  \code{\link[xpose4generic]{xpose.panel.default}} for details.
}
\value{
  Returns a compound plot comprising plots of observations (DV),
  individual predictions (IPRED), and population predictions (PRED) 
  against the independent variable (IDV).
}
\author{E. Niclas Jonsson, Mats Karlsson, Andrew Hooker & Justin Wilkins}
\seealso{
  \code{\link{dv.vs.idv}},
  \code{\link{ipred.vs.idv}},
  \code{\link{pred.vs.idv}},
  \code{\link[xpose4generic]{xpose.plot.default}}, 
  \code{\link[xpose4generic]{xpose.panel.default}},
  \code{\link[lattice]{xyplot}},
  \code{\link[xpose4data]{xpose.prefs-class}},
  \code{\link[xpose4data]{xpose.data-class}}
}
\examples{
\dontrun{
## We expect to find the required NONMEM run and table files for run
## 5 in the current working directory
xpdb5 <- xpose.data(5)
}

## Here we load the example xpose database from the xpose4data package
data(simpraz.xpdb)
xpdb <- simpraz.xpdb

## A vanilla plot
dv.preds.vs.idv(xpdb)

## Custom colours and symbols, IDs
dv.preds.vs.idv(xpdb, cex=0.6, pch=8, col=1, ids=TRUE)
}
\keyword{methods}
