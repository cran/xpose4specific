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

runsum.old <- function(object,
                   modfile=paste("run",object@Runno,".mod",sep=""),
                   listfile=paste("run",object@Runno,".lst",sep=""),
                   main=NULL,
                   subset=xsubset(object),
                   ...
                   ) {


  ## Read model file
  if(is.readable.file(modfile)) {
    modfile <- scan(modfile,sep="\n",what=character(),quiet=TRUE)
    mod.file.lines <- length(modfile)
  } else {
    return()
  }

  ## Global settings concerning number of lines, number of columns etc.
   
  txtnrow  <- 63       # Number of rows in each column
  minrow   <- 5        # Minimum number of rows in each column
  
  parameter.list <- create.parameter.list(listfile)
  attach(parameter.list,warn.conflicts=F)

  ## Set up screen
  ##oldpar <- par(mar=c(2.7,2.7,2,1.5)+0.1,
  ##		mgp=c(1.7,0.6,0))
  ##  if(!printer)
  ##    on.exit(par(oldpar),add=T)

  ## Need to set mar to something very small, maybe for all screens.
  ## Set up coordinates for original screen
  m      <- matrix(0,6,4)
  #m[1,]  <- c(0,    1,    0.97, 1   )         # Title
  #m[2,]  <- c(0,    0.25, 0.77, 0.97)         # First gof
  #m[3,]  <- c(0.25, 0.5,  0.77, 0.97)         # Second gof
  #m[4,]  <- c(0.5,  0.75, 0.77, 0.97)         # Third gof
  #m[5,]  <- c(0.75, 1,    0.77, 0.97)         # Fourth gof
  #m[6,]  <- c(0,    1,    0,    0.77)         # Text area
  
  m[1,]  <- c(0,    1,    0.93, 1   )         # Title
  m[2,]  <- c(0,    0.25, 0.73, 0.93)         # First gof
  m[3,]  <- c(0.25, 0.5,  0.73, 0.93)         # Second gof
  m[4,]  <- c(0.5,  0.75, 0.73, 0.93)         # Third gof
  m[5,]  <- c(0.75, 1,    0.73, 0.93)         # Fourth gof
  m[6,]  <- c(0,    1,    0,    0.71)         # Text area
  
  close.screen(all=T)
  
  #cat("Check 1\n")
  
  split.screen(m)
  
  ## Add the plots
  add.plots(object,main=main,
            subset=subset,
            idscex=0.6,
            ...
            )
            
  #cat("Check 1 fin\n")

  ## Split screen 6 for text
  screen(6)
  split.screen(figs=c(1,2))
  
  #cat("Check 2\n") 
  screen(7)
  par(mar=c(0,0.5,0,0.5)+0.3)
  
  plot(0,0,type="n",ylim=c(1,txtnrow),xlim=c(0,8),xaxt="n",yaxt="n",
       ylab="",xlab="",bty="n")

  ## Add the termination messages
  if(seenterm == 1) {
    termtxt <- term
    
    xval  <- rep(0,length(termtxt))
    yval  <- txtnrow:(txtnrow-length(termtxt)+1)
  
    text(xval,yval,termtxt,adj=0,cex=0.4,font=2)
    currow <- txtnrow-length(termtxt)
  } else {
    currow <- txtnrow
  }

  ## Add objective
  if(seenobj == 1) {
    text(0,currow-1,paste("Objective:",ofv),adj=0,cex=0.4,font=2)
    currow <- currow-1
  } 
  

  ## Now we should print the parameters and their names
  ## Decide how many rows in each column of parameters
  if(seenseth ==1 || seenseom==1 || seensesi==1) {
    
    if(npar <= minrow*3) {
      par.in.col <- minrow
    } else {
      par.in.col <- floor(npar/3+1)
    }
    
    have.ses     <- 1
    ncols.needed <- 3
  } else {
    
    if(npar <= minrow*4) {
      par.in.col <- minrow
    } else {
      par.in.col <- floor(npar/4+1)
      }
    
    have.ses     <- 0
    ncols.needed <- 4
  }
    
  
  ## To make it easier to print the parameters we add empty entries
  ## in parnam, parval, separnam and separval to fill all columns
  
  if(have.ses == 0) max.par      <- 4*par.in.col
  if(have.ses == 1) max.par      <- 3*par.in.col
  extra.pars   <- max.par-npar

  if(extra.pars > 0) {
    parnam[(length(parnam)+1):max.par]   <- ""
    separnam[(length(separnam)+1):max.par] <- ""
    parval[(length(parval)+1):max.par]   <- ""
    separval[(length(separval)+1):max.par] <- ""
  }

  
  ## Print the parameters
  start.col <- 0
  if(have.ses == 0){
    start.col[1:4] <- c(0,2,4,6)
  }

  if(have.ses == 1)  {
    start.col[1:3] <- c(0,3,6)
  }

  curr.col <- 0
  for( i in start.col) {
    curr.col     <- curr.col+1
    range.to.print <- 
      (1+length(parnam)-par.in.col*(ncols.needed-(curr.col-1))):
	(length(parnam)-par.in.col*(ncols.needed-curr.col))
    
    ## Print the parameter names
    xval <- rep(i,1)
    yval <- currow-2
    if(!all(parnam[range.to.print]=="")){
      text(xval,yval,"Par",adj=0,cex=0.4,font=2)
    }
    xval <- rep(i,par.in.col)  
    yval <- (currow-3):(currow-2-par.in.col)
    text(xval,yval,parnam[range.to.print],adj=0,cex=0.4,font=2)
    
    ## Print the parameter values
    xval <- rep(i+1,1)
    yval <- currow-2
    if(!all(parval[range.to.print]=="")){
      text(xval,yval,"Val",adj=0,cex=0.4,font=2)
    }
    xval <- rep(i+1,par.in.col)
    yval <- (currow-3):(currow-2-par.in.col)
    text(xval,yval,parval[range.to.print],adj=0,cex=0.4,font=2)
    
    if(have.ses == 0) next
    
#### To save space the SE names are not printed. The code for, however, are left
#### if we in the future want to add this featre.

    ## Print the SE names
    ##  xval <- rep(i+2,par.in.col)
    ##  yval <- (currow-2):(currow-1-par.in.col)
    ##  text(xval,yval,separnam[range.to.print],adj=0,cex=0.4,font=2)
    ## Print the parameter values
    xval <- rep(i+2,1)
    yval <- currow-2
    if(!all(separval[range.to.print]=="")){
      text(xval,yval,"CV",adj=0,cex=0.4,font=2)
    }
    xval <- rep(i+2,par.in.col)
    yval <- (currow-3):(currow-2-par.in.col)
    text(xval,yval,separval[range.to.print],adj=0,cex=0.4,font=2)
  }

  currow <- currow-2-par.in.col
  
  ## Now it is time to print the mdoel file

  ## Number of lines left in column

  num.lines.left <- currow - 2

  ## If there is room for model file in the current column
  if(mod.file.lines < num.lines.left) {

    xval <- rep(0,mod.file.lines)
    yval <- (currow-2):(currow-1-mod.file.lines)

    text(xval,yval,modfile,adj=0,cex=0.4,font=2)
    
  } else {
    
    xval <- rep(0,num.lines.left)
    yval <- (currow-2):(currow-1-num.lines.left)

    text(xval,yval,modfile[1:num.lines.left],adj=0,cex=0.4,font=2)

    ## If we need another column

    #cat("Check 3\n")
    screen(8)
    par(mar=c(0,0.5,0,0.5)+0.3)
    plot(0,0,type="n",ylim=c(1,txtnrow),xlim=c(0,8),xaxt="n",yaxt="n",
	 ylab="",xlab="",bty="n")

    lines.left.to.print <- mod.file.lines-num.lines.left
    
    currow <- txtnrow
    xval <- rep(0,lines.left.to.print)
    yval <- currow:(currow-lines.left.to.print+1)
    
    text(xval,yval,modfile[(num.lines.left+1):mod.file.lines],
	 adj=0,cex=0.4,font=2)
  }

  detach(parameter.list)
  invisible()
  
  #return()

}



"add.plots" <-
  function(object,
           main=NULL,
           subset=NULL,
           inclZeroWRES = FALSE,
           onlyfirst    = FALSE,
           col   = object@Prefs@Graph.prefs$col,
           pch   = object@Prefs@Graph.prefs$pch,
           lty   = object@Prefs@Graph.prefs$lty,
           smlwd = object@Prefs@Graph.prefs$smlwd,
           smlty = object@Prefs@Graph.prefs$smlty,
           smcol = object@Prefs@Graph.prefs$smcol,
           abllwd= object@Prefs@Graph.prefs$abllwd,
           abllty= object@Prefs@Graph.prefs$abllty,
           ablcol= object@Prefs@Graph.prefs$ablcol,
           ids   =object@Prefs@Graph.prefs$ids,
           idsmode=object@Prefs@Graph.prefs$idsmode,
           idsext =object@Prefs@Graph.prefs$idsext,
           idscex= object@Prefs@Graph.prefs$idscex,
           idsdir= object@Prefs@Graph.prefs$idsdir,
           ...
           )
{

  if(!is.null(subset)) {
    maintit <- paste("Summary of run ",object@Runno,
                     ", ",subset,sep="")
  } else {
    maintit <- paste("Summary of run ",object@Runno,sep="")
  }

  ## Get the data
  data <- Data(object,inclZeroWRES,onlyfirst=onlyfirst,subset=subset)

  if(is.null(xvardef("dv",object))){
        cat("Dependent variable (DV) is not defined in database\n")
        cat("and is required to display plots\n")
    return()
  }
  if(is.null(xvardef("ipred",object))){
        cat("IPRED is not defined in database\n")
        cat("and is required to display plots\n")
    return()
  }
  if(is.null(xvardef("iwres",object))){
        cat("IWRES is not defined in database\n")
        cat("and is required to display plots\n")
    return()
  }
  if(is.null(xvardef("wres",object))){
        cat("WRES is not defined in database\n")
        cat("and is required to display plots\n")
    return()
  }
  if(is.null(xvardef("idv",object))){
        cat("independent variable (IDV) is not defined in database\n")
        cat("and is required to display plots\n")
    return()
  }
  if(is.null(xvardef("pred",object))){
        cat("PRED is not defined in database\n")
        cat("and is required to display plots\n")
    return()
  }
  if(is.null(xvardef("id",object))){
        cat("ID is not defined in database\n")
        cat("and is required to display plots\n")
    return()
  }
  
  ids<-unique(data[,xvardef("id",object)])

  ## Set summary title
  screen(1)
  par(mar=c(0,0,0,0))
  
  plot(0,0,
       type="n",
       ylim=c(0,1),
       xlim=c(0,1),
       xaxt="n",
       yaxt="n",
       bty="n",
       ylab="",
       xlab="")
  text(0.5,0.5,maintit,cex=1.2)

  ## Add plots
  ## PRED vs DV
  screen(2)
  par(mar=c(0,0.5,0,0.5)+0.3)

  ylm<-range(data[,xvardef("dv",object)])
  xlm<-range(data[,xvardef("pred",object)])

  plot(0,0,
       type="n",
       bty="o",
       ylim=ylm,
       xlim=xlm,
       lty=1,
       lwd=0.5,
       ylab="",
       xlab="",
       cex=0.5,
       tcl=-0.1,
       mgp = c(3,0.2,0),
       cex.axis=0.7)

  IDS  <- data[,xvardef("id",object)]
  DV   <- data[,xvardef("dv",object)]
  IDV  <- data[,xvardef("idv",object)]
  PRED <- data[,xvardef("pred",object)]
  IPRED<- data[,xvardef("ipred",object)]
  IWRES<- data[,xvardef("iwres",object)]
  WRES <- data[,xvardef("wres",object)]

  for (i in ids) {
    seli <- IDS == i
    yplt <- DV[seli]
    xplt <- PRED[seli]
    xord <- sort.list(xplt)
    xplt <- xplt[xord]
    yplt <- yplt[xord]
    lines(xplt,yplt,lty=lty,pch=pch,col=col,lwd=0.2)
  }

  if(!is.null(ids)) {
    addid(PRED,DV,data[,xvardef("id",object)],
          idsmode=idsmode,
          idsext=idsext,
          idscex=idscex,
          idsdir=idsdir,
          gridmode=FALSE)
  }
  
  abline(0,1,lty=abllty,lwd=abllwd,col=ablcol)
  ##title <- paste(xlabel(xvardef("dv",object),object),"vs",
  ##               xlabel(xvardef("pred",object),object))
  title <- paste("DV","vs","PRED")
  mtext(title,side=3,line=0.2,cex=0.5)

  ## IPRED vs DV
  screen(3)
  par(mar=c(0,0.5,0,0.5)+0.3)

  ylm<-range(DV)
  xlm<-range(IPRED)

  plot(0,0,
       type="n",
       bty="o",
       ylim=ylm,
       xlim=xlm,
       lty=1,
       lwd=0.5,
       ylab="",
       xlab="",
       cex=0.5,
       tcl=-0.1,
       mgp = c(3,0.2,0),
       cex.axis=0.7)

  for (i in ids) {
    seli <- IDS == i
    yplt <- DV[seli]
    xplt <- IPRED[seli]
    xord <- sort.list(xplt)
    xplt <- xplt[xord]
    yplt <- yplt[xord]
    lines(xplt,yplt,lty=lty,pch=pch,col=col,lwd=0.2)
  }

  if(!is.null(ids)) {
    addid(IPRED,DV,data[,xvardef("id",object)],
          idsmode=idsmode,
          idsext=idsext,
          idscex=idscex,
          idsdir=idsdir,
          gridmode=FALSE)

  }
  
  abline(0,1,lty=abllty,lwd=abllwd,col=ablcol)
  ##title <- paste(xlabel(xvardef("dv",object),object),"vs",
  ##               xlabel(xvardef("ipred",object),object))
  title <- paste("DV","vs","IPRED")
  mtext(title,side=3,line=0.2,cex=0.5)

  ## ABS IWRES vs IPRED
  screen(4)
  par(mar=c(0,0.5,0,0.5)+0.3)
  
  plot(IPRED,
       abs(IWRES),
       type="p",
       bty="o",
       mkh=0.05,
       pch=pch,
       col=col,
       ylab="",
       xlab="",
       cex=0.5,
       tcl=-0.1,
       mgp = c(3,0.2,0),
       cex.axis=0.7)

  lines(lowess(IPRED,abs(IWRES)),col=smcol,lwd=smlwd,lty=smlty)
  ##title <- paste("|",
  ##               xlabel(xvardef("iwres",object),object),"| vs ",
  ##               xlabel(xvardef("ipred",object),object),sep=" ")
  title <- paste("|IWRES|","vs","IPRED")
  mtext(title,side=3,line=0.2,cex=0.5)

  ## WRES vs TIME 
  screen(5)
  par(mar=c(0,0.5,0,0.5)+0.3)
  
  xlm<-range(IDV)
  ylm<-range(WRES)
  
  plot(0,0,
       type="n",
       bty="o",
       lty=1,
       lwd=0.5,
       ylim=ylm,
       xlim=xlm,
       ylab="",
       xlab="",
       cex=0.5,
       tcl=-0.1,
       mgp = c(3,0.2,0),
       cex.axis=0.7)
  
  lines(lowess(IDV,WRES))
  
  for (i in ids) {
    seli<-IDS==i
    xplt<-IDV[seli]
    yplt<-WRES[seli]
    xord <- sort.list(xplt)
    xplt <- xplt[xord]
    yplt <- yplt[xord]
    lines(xplt,yplt,lty=lty,col=col,pch=pch,lwd=0.2)
  }

  if(!is.null(ids)) {
    addid(IDV,WRES,data[,xvardef("id",object)],
          idsmode=idsmode,
          idsext=idsext,
          idscex=idscex,
          idsdir=idsdir,
          gridmode=FALSE)
  }
  
  abline(0,0,lty=abllty,lwd=abllwd,col=ablcol)
  ##title <- paste(xlabel(xvardef("wres",object),object),"vs",
  ##               xlabel(xvardef("idv",object),object))
  title <- paste("WRES","vs",xlabel(xvardef("idv",object),object))
  mtext(title,side=3,line=0.2,cex=0.5)



}
