# Xpose 4
# An R-based population pharmacokinetic/
# pharmacodynamic model building aid for NONMEM.
# Copyright (C) 1998-2004 E. Niclas Jonsson and Mats Karlsson.
# Copyright (C) 2005-2008 Andrew C. Hooker, Justin J. Wilkins, 
# Mats O. Karlsson and E. Niclas Jonsson.
# Copyright (C) 2009-2010 Andrew C. Hooker, Mats O. Karlsson and 
# E. Niclas Jonsson.

# This file is a part of Xpose 4.
# Xpose 4 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation, either version 3
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program.  A copy can be cound in the R installation
# directory under \share\licenses. If not, see http://www.gnu.org/licenses/.

"ind.plots.wres.qq" <-
  function(object,
           main = "Default",
                                        #xlb  = NULL,
                                        # ylb  = xlabel(xvardef("dv",object),object),
                                        #ylb = NULL,
           layout=c(4,4),
           inclZeroWRES=FALSE,
           subset=xsubset(object),
           type="o",
           pch=object@Prefs@Graph.prefs$pch,
           col=object@Prefs@Graph.prefs$col,
           cex=object@Prefs@Graph.prefs$cex,
           abllty = object@Prefs@Graph.prefs$abllty,
           abllwd = object@Prefs@Graph.prefs$abllwd,
           ablcol = object@Prefs@Graph.prefs$ablcol,
           prompt = FALSE,
           main.cex=0.9,
           mirror=NULL,
           max.plots.per.page=1,
           ...) {

    ## check for mirror
    if(!is.null(mirror)){
      cat("Mirror not currently implemented for individual plots\n")
      return()
    }

    ## Make sure we have the necessary variables defined in the ##
    ## object.                                                  ##
    if(is.null(check.vars(c("id","wres"),object))) {
      cat("The ID and WRES variables all need to be defined in the current database!\n")
      return(NULL)
    }

    ## Fix any main and/or axis titles
    default.plot.title <- "Individual plots of WRES"
    plotTitle <- xpose.multiple.plot.title(object=object,
                                           plot.text = default.plot.title,
                                           main=main,
                                           ...)
    
    
    ## Bin them
    length.id <- length(unique(object@Data[[xvardef("id",object)]]))
    list.id   <- unique(object@Data[[xvardef("id",object)]])
    plots.per.page <- layout[1] * layout[2]
    plots.cur <- 0
    pages <- 1
    page.breaks <- c(0)
    old.obj <- object
    new.obj <- object
    new.obj@Data <- NULL
    
    for (i in list.id) {
      plots.cur <- plots.cur + 1
      if (plots.cur == plots.per.page) {
        pages <- pages + 1
        plots.cur <- 0
        page.breaks <- c(page.breaks, i)
      }
    }
    if (max(page.breaks) < max(list.id)) {
      page.breaks <- c(page.breaks, max(list.id))
    }
    id.levels <- levels(cut(object@Data$ID, page.breaks, include.lowest=T))
    old.obj@Data$bin <- cut(object@Data$ID, page.breaks, include.lowest=T)
    

    plot.num <- 0
    plotList <- vector("list",length(id.levels))     
    for (i in id.levels) {    ## start loop

      new.obj@Data <- subset(old.obj@Data, bin == i)

      
      ## Set up the data ##
      ## Figure out what variables we have defined
      select <- xvardef("wres",object)
      
      numpans <- length(select)

      nobj <- new("xpose.data",
                  Runno=object@Runno,
                  Data = NULL 
                  )
      Data(nobj) <- Data(new.obj,inclZeroWRES=inclZeroWRES,
                         subset=subset)
      
      
      xplot <- xpose.plot.qq(xvardef("wres",nobj),
                             new.obj,
                                        #xlb = xlb,
                                        #ylb = ylb,
                             by=xvardef("id",nobj),
                             main=plotTitle,
                                        #group="ind",
                             layout=layout,
                             scales=list(cex=0.7,tck=0.5),
                             aspect="fill",
                             xvar = xvardef("wres",object),
                             force.by.factor=TRUE,
                             ids=F,
                             pch=pch,
                                        #col=col,
                             abllty=abllty,
                             abllwd=abllwd,
                             ablcol=ablcol,
                             subset=subset,
                             as.table=TRUE,
                             main.cex=main.cex,
                             ...)
      plot.num <- plot.num+1
      plotList[[plot.num]] <- xplot
      

    }

    obj <- xpose.multiple.plot(plotList,max.plots.per.page=max.plots.per.page,plotTitle=NULL,prompt=prompt,...)
    return(obj)

  }
