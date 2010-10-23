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

"ind.plots.cwres.hist" <-
  function(object,
           main = "Default",
           #xlb  = NULL,
                                        # ylb  = xlabel(xvardef("dv",object),object),
           ylb = NULL,
           layout=c(4,4),
           inclZeroWRES=FALSE,
           subset=xsubset(object),
           hicol = object@Prefs@Graph.prefs$hicol,
           hilty = object@Prefs@Graph.prefs$hilty,
           hilwd = object@Prefs@Graph.prefs$hilwd,
           hidcol = object@Prefs@Graph.prefs$hidcol,
           hidlty = object@Prefs@Graph.prefs$hidlty,
           hidlwd = object@Prefs@Graph.prefs$hidlwd,
           hiborder = object@Prefs@Graph.prefs$hiborder,
           prompt = FALSE,
           mirror=NULL,
           main.cex=0.9,
           max.plots.per.page=1,
                                        #lty=c(0,1,1),
                                        #pch=c(21,32,32),
                                        #type="o",
                                        #col=c(1,object@Prefs@Graph.prefs$smcol,object@Prefs@Graph.prefs$lmcol),
                                        #lwd=1,
           ...) {

    ## Make sure we have the necessary variables defined in the ##
    ## object.                                                  ##
    if(is.null(check.vars(c("id","cwres"),object))) {
      return(NULL)
    }
    

    ## check for mirror
    if(!is.null(mirror)){
      cat("Mirror not currently implemented for individual plots\n")
      return()
    }

    
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
      select <- xvardef("cwres",object)
      
      numpans <- length(select)

      nobj <- new("xpose.data",
                  Runno=object@Runno,
                  Data = NULL 
                  )
      Data(nobj) <- Data(new.obj,inclZeroWRES=inclZeroWRES,
                         subset=subset)
      
      ## Fix any main and/or axis titles
      default.plot.title <- "Individual plots of CWRES"
      plotTitle <- xpose.multiple.plot.title(object=object,
                                             plot.text = default.plot.title,
                                             main=main,
                                             ...)

      
                                        # Set y axis title
##       if (is.null(xlb)) {
##         xlb <- xlabel(xvardef("cwres",object),object)
##       }
##       if (is.null(ylb)) {
##         ylb <- "Proportion"
##       }
      
      
      xplot <- xpose.plot.histogram(xvardef("cwres",nobj),
                                    new.obj,
                                    #xlb = xlb,
                                    #ylb = ylb,
                                    by=xvardef("id",nobj),
                                    main=plotTitle,
                                        #group="ind",
                                    layout=layout,
                                    scales=list(cex=0.7,tck=0.5),
                                    aspect="fill",
                                    xvar = xvardef("cwres",object),
                                    force.by.factor=TRUE,
                                    ids=F,
                                    subset=subset,
                                    as.table=TRUE,
                                    hicol = hicol,
                                    hilty = hilty,
                                    hilwd = hilwd,
                                    hidcol = hidcol,
                                    hidlty = hidlty,
                                    hidlwd = hidlwd,
                                    hiborder = hiborder,
                                    main.cex=main.cex,
                                    ...)
      plot.num <- plot.num+1
      plotList[[plot.num]] <- xplot
      

    }

    obj <- xpose.multiple.plot(plotList,max.plots.per.page=max.plots.per.page,prompt=prompt,
                               plotTitle=NULL,...)
    return(obj)

  }
