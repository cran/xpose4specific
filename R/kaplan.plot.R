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

kaplan.plot <-
  function(x="TIME",y="DV",id="ID",
           data=d,
           evid="EVID",
           by=NULL,
           xlab="Time",
           ylab="Survival (%)",
           object=NULL,
           events.to.plot="All",
           sim.data=NULL,
           sim.zip.file=NULL,
           VPC = FALSE,
           nsim.lab="simNumber",
           sim.evct.lab="counter",
           probs=c(0.025,0.975),
           add.baseline=T,
           add.last.area=T,
           subset=NULL,
           ##subset.real=NULL,
           main="Default",
           nbins=NULL,
           real.se= if(!is.null(sim.data)) F else T,
           real.se.type="l",
           real.type="l",
           real.lwd=1,
           real.se.lty=2,
           real.se.lwd=0.5,
           real.se.col="red",
           inclZeroWRES=TRUE,
           onlyfirst=FALSE,
           #RTTE=FALSE,
           samp=NULL,
           poly.alpha=1,
           poly.fill="lightgreen",
           poly.line.col="darkgreen",
           poly.lty=2,
           censor.lines=TRUE,
           ylim=c(-5,105),
           ...)
{

  ##Get data for xpose object
  if(!is.null(object)){
    if(!is.null(samp)) {
      data <- SData(object,inclZeroWRES,onlyfirst=onlyfirst,
                    subset=subset,samp=samp)
    } else {
      data <- Data(object,inclZeroWRES,onlyfirst=onlyfirst,subset=subset)
    }
  } else {
    if(!is.null(subset)){
      on.exit(detach(data))
      attach(data,warn.conflicts=F)
      data <- data[eval(parse(text=paste("data$", subset))),]
      if(dim(data)[1]==0) return(NULL)
    }
  }
  if (VPC && !is.null(object)) sim.zip.file <- paste("simtab",object@Runno,".zip",sep="")
  if (!is.null(sim.zip.file)){
      sim.data <- read.table(unz(sim.zip.file,
                                 sub(".zip", "",sim.zip.file)),
                             skip=0,header=T)
  }
  if(!is.null(sim.data)){
      if(!is.null(subset)){
          on.exit(detach(sim.data))
          attach(sim.data,warn.conflicts=F)
          sim.data <- sim.data[eval(parse(text=paste("sim.data$", subset))),]
          if(dim(sim.data)[1]==0) return(NULL)
      }
  }


  ## Check if this is RTTE or TTE data
  data$counter <- 0
  counter <- 0
  old.id <- 0
  for(i in 1:length(data$counter)){
      new.id <- data[[id]][i]
      if(new.id != old.id){counter <- 0}
      old.id <- new.id
      if(!is.null(data[[evid]])){ # filter on evid
          if (data[[evid]][i]==0) counter=counter+1
      } else { # all rows are events
          counter=counter+1
      }
      data$counter[i]=counter
  }
  events <- max(data$counter)
  RTTE = FALSE
  if(events>1) RTTE = TRUE
  max.events <- events

  if(!is.null(sim.data)){
      sim.events <- max(sim.data[[sim.evct.lab]])
      #if(sim.events>1) RTTE = TRUE
      #max.events <- max(sim.events,events)
  }

  full.data <- data
  y.true <- y
  if(!is.null(sim.data)){
    full.sim.data <- sim.data

    ## add sim.ID
    ## add a unique ID identifier to simulated data.
    rle.result <- rle(full.sim.data[[id]])
    rle.result$values <- 1:length(rle.result$values)
    new.id <- inverse.rle(rle.result)
    full.sim.data$sim.ID <- new.id
  }
  
  #browser()
  #by <- c("DOSE==0","DOSE!=0")
  by.no <- length(by)
  #by.no
  if(!all(is.na(match(events.to.plot,"All")))){
      num.of.plots <- max.events
      if(by.no>0) num.of.plots <- num.of.plots*by.no
      plotList <- vector("list",num.of.plots)
      event.list <- 1:max.events
  } else {
      num.of.plots <- length(events.to.plot)
      if(by.no>0) num.of.plots <- num.of.plots*by.no
      plotList <- vector("list",num.of.plots)
      event.list <- events.to.plot
  }
  plot.num <- 0
  for(event.no in event.list){
      if(is.null(by)) by=c(1)
      for(by.val in by){
                                        #i=1 # for testing
          if(by.val==1){
              by.val=NULL
              by=NULL
          }
          if(is.null(by.val)){
              data <- full.data
              if(!is.null(sim.data)) sim.data <- full.sim.data
          } else {
              data <- full.data[eval(parse(text=paste("full.data$", by.val))),]
              if(dim(data)[1]==0) return(NULL)
              if(!is.null(sim.data)){
                  sim.data <- full.sim.data[eval(parse(text=paste("full.sim.data$", by.val))),]
                  if(dim(sim.data)[1]==0) return(NULL)
              }
          }

          tmp.name <- paste("Event",event.no,sep=" ")
          if(!is.null(by.val)) tmp.name <- paste(tmp.name,by.val,sep=", ")
          tmp.cex <- 0.8

          data <- subset(data,counter<=event.no)
          if(!is.null(data[[evid]])){
              data <- data[data[evid]==0,]
          }

          if(!is.null(sim.data)) sim.data <- subset(sim.data,eval(parse(text=sim.evct.lab))<=event.no)
          
          
          data <- data[!duplicated(data[[id]], fromLast = TRUE),]
          if(!is.null(sim.data)) sim.data <- sim.data[!duplicated(sim.data$sim.ID, fromLast = TRUE),]

          ## account for censoring and values larger than 1
          data$tmp.event <- 0
          data$tmp.event[data[y.true]!=0] = 1 # values larger than 1 set to 1
          data$tmp.event[data$counter < event.no] = 0 # censored events set to zero
          y="tmp.event"


          ## need to account for simulated data as well
          if(!is.null(sim.data)){
              sim.data$tmp.event <- 0
              sim.data$tmp.event[sim.data[y.true]!=0] = 1 # values larger than 1 set to 1
              sim.data$tmp.event[sim.data[[sim.evct.lab]] < event.no] = 0 # censored events set to zero
          }

                                        #sim.data$tmp.event <- sim.data[[y]]


          ##d.sub[c("ID","new.DV","TIME")]

          ##d.sub$ETMP = 0
                                        #d.sub$ETMP[d.sub$TYPE!=3] = 1
                                        #names(d.sub)
                                        #kaplan.plot(x="TIME",y="DV",data=d.sub)


                                        #browser()

          ## make new column for when censoring (no event) occurs
          ## if(!RTTE){
          ##   data$tmp.event <- 0
          ##   data$tmp.event[data[y]!=0] = 1
          ##   sim.data$tmp.event <- sim.data[[y]]
          ##   y="tmp.event"
          ## }
          ##browser()

          S <- Surv(data[,x],data[,y])
          ##f.1 <- survfit(S)

          f.1 <- survfit(S~1)
          a.1 <- summary(f.1)
          ##plot(f.1)


          if(!is.null(sim.data)){
              times <- sort(unique(sim.data[,x]))
              m1 <- matrix(nrow=length(unique(sim.data[,nsim.lab])),#number of simulaitons
                           ncol=length(times)) # number of times
              for(i in unique(sim.data[,nsim.lab])){
                  ##tmp <- subset(sim.data[,nsim==i],nsim==i)
                  tmp <- sim.data[eval(parse(text=paste("sim.data$", nsim.lab,"==",i))),]
                  S.sim <- Surv(tmp[,x],tmp[,y])
                  ##f.sim <- survfit(S.sim~tmp$TRT)
                  f.1.sim <- survfit(S.sim~1)
                  ##a.sim <- summary(f.sim)
                  a.1.sim <- summary(f.1.sim)
                  tmp.times <- f.1.sim$time
                  col.index <- match(tmp.times,times)

                  ## plot the simulated k-m curves
                                        #browser()
                                        #plot(f.1.sim)

                  ## add the new information
                  m1[i,col.index] <- f.1.sim$surv
                  ## make sure the total prop is carried through in the times
                  if(is.na(m1[i,1])){
                      m1[i,1] <- 1
                  }
                  for(j in 2:length(times)){
                      if(is.na(m1[i,j])){
                          m1[i,j] <- m1[i,j-1]
                      }
                  }
              }


                                        #browser()
                                        #m1[,1] # this list should actally not contain any NA
              ## should be 1 or a value.
                                        #dim(m1)
                                        #m1[,129]

              if(is.null(nbins)){
                  quants <- matrix(nrow=2,ncol=dim(m1)[2])
                  time.bin <- times
                  for(j in 1:dim(quants)[2]){
                      quants[,j] <- quantile(m1[,j],probs=probs,na.rm=F)
                  }
              } else {
                  remainder.times <- 0
                  ncomb <- floor(dim(m1)[2]/nbins)

                  if(ncomb==0){
                      ncomb <- 1
                      quants <- matrix(nrow=2,ncol=dim(m1)[2])
                      time.bin <- c()
                  } else {
                      remainder.times <- dim(m1)[2]-ncomb*nbins
                      quants <- matrix(nrow=2,ncol=nbins)
                      time.bin <- c()
                  }

                  k.add.this <- 0
                  for(j in 1:dim(quants)[2]){
                      ##j <- 1
                      tmp.bin <- c()
                      if(remainder.times>0){
                          ncomb.tmp <- ncomb+1
                          remainder.times <- remainder.times - 1
                      } else {
                          ncomb.tmp <- ncomb
                      }
                      for(k in 1:ncomb.tmp){
                          ##k.add.this <- (j-1)*ncomb
                          if(k==1) time.bin <- c(time.bin,times[k+k.add.this])
                          if(k+k.add.this>dim(m1)[2]) next
                          tmp.bin <- c(tmp.bin,m1[,k+k.add.this])
                      }
                      k.add.this <- k.add.this+ncomb.tmp
                      quants[,j] <- quantile(tmp.bin,probs=probs,na.rm=F)
                  }
              }

              if(add.baseline){
                  tmp.mat <- matrix(nrow=dim(quants)[1],ncol=dim(quants)[2]+1)
                  tmp.mat[,1] <- 1
                  tmp.mat[,-1] <- quants
                  tmp.x <- c(0,time.bin)
              } else {
                  tmp.mat <- quants
                  tmp.x <- time.bin
              }

              PI.up <- c()
              PI.down <- c()
              PI.times <- c()
              n.times <- length(tmp.x)
              for(i in 1:n.times){
                  PI.reps=2
                  time.reps=2
                  if(i==1) time.reps=1
                  if(i==n.times) PI.reps=1
                  PI.up <- c(PI.up,rep(tmp.mat[2,i],PI.reps))
                  PI.down <- c(PI.down,rep(tmp.mat[1,i],PI.reps))
                  PI.times <- c(PI.times,rep(tmp.x[i],time.reps))
              }

              if(add.last.area){
                  if(tail(tail(PI.times,n=1))==tail(times,n=1)){
                      time.pt <- (tail(PI.times,n=1)-PI.times[1])*1.01
                  } else {
                      time.pt <- (tail(times,n=1)-PI.times[1])*1.01
                  }
                  PI.times=c(PI.times,time.pt)
                  PI.up <- c(PI.up,tail(PI.up,n=1))
                  PI.down <- c(PI.down,tail(PI.down,n=1))
              }

              PI.up <- PI.up*100
              PI.down <- PI.down*100
          } else{
              PI.up=NULL
              PI.down=NULL
              PI.times=NULL
          } # end simulation stuff

          if(add.baseline){
              tmp.y <- c(1,f.1$surv)
              tmp.x <- c(0,f.1$time)
              tmp.y.upper <- c(1,f.1$upper)
              tmp.y.lower <- c(1,f.1$lower)
          } else {
              tmp.y <- c(f.1$surv)
              tmp.x <- c(f.1$time)
              tmp.y.upper <- c(f.1$upper)
              tmp.y.lower <- c(f.1$lower)
          }

          real.data <- c()
          real.times <- c()
          real.data.upper <- c()
          real.data.lower <- c()
          for(i in 1:length(tmp.y)){
              y.reps=2
              x.reps=2
              if(i==1) x.reps=1
              if(i==length(tmp.y)) y.reps=1
              real.data <- c(real.data,rep(tmp.y[i],y.reps))
              real.times <- c(real.times,rep(tmp.x[i],x.reps))
              real.data.upper <- c(real.data.upper,rep(tmp.y.upper[i],y.reps))
              real.data.lower <- c(real.data.lower,rep(tmp.y.lower[i],y.reps))
          }

          real.data <- real.data*100
          real.data.upper <- real.data.upper*100
          real.data.lower <- real.data.lower*100

          xplot <-
              xyplot(real.data~real.times,
                     main=list(tmp.name,cex=tmp.cex),
                     ylim =ylim,
                     xlab=xlab,ylab=ylab,
                     real.type=real.type,
                     PI.up=PI.up,PI.down=PI.down,
                     PI.times=PI.times,
                     real.data.upper=real.data.upper,
                     real.data.lower=real.data.lower,
                     real.se.type=real.se.type,
                     real.se.lty=real.se.lty,
                     real.se.col=real.se.col,
                     f.1=f.1,
                     ...,
                     panel=function(x,y,PI.up,PI.down,PI.times,
                     real.data.upper,real.data.lower,
                     real.se.type=real.se.type,
                     real.type=real.type,
                     real.se.lty=real.se.lty,
                     real.se.col=real.se.col,
                     f.1=f.1,
                     ...){
                         if(!is.null(sim.data)){
                             grid.polygon(c(PI.times,rev(PI.times)),c(PI.up,rev(PI.down)),
                                          default.units="native",
                                          gp=gpar(fill=poly.fill,alpha=poly.alpha,col=poly.line.col,lty=poly.lty)
                                          )
                         }
                         panel.xyplot(x,y,type=real.type,lwd=real.lwd,...)
                         if(real.se) panel.xyplot(x,real.data.upper,type=real.se.type,
                                                  lty=real.se.lty,
                                                  col=real.se.col,
                                                  lwd=real.se.lwd,
                                                  ...)
                         if(real.se) panel.xyplot(x,real.data.lower,type=real.se.type,
                                                  lty=real.se.lty,
                                                  col=real.se.col,
                                                  lwd=real.se.lwd,
                                                  ...)

                         if(censor.lines){
                             if(any(f.1$n.censor>0)){
                                 panel.segments(x0=c(f.1$time[f.1$n.censor>0]),
                                                y0=c(f.1$surv[f.1$n.censor>0]*100-2),
                                                x1=c(f.1$time[f.1$n.censor>0]),
                                                y1=c(f.1$surv[f.1$n.censor>0]*100+2),
                                                ...)
                             }
                         }
                     }
                     )
          plot.num <- plot.num+1
          #if(plot.num>152) browser()
          plotList[[plot.num]] <- xplot
      }
  }
  default.plot.title <- "Kaplan-Meyer plots"
  if (num.of.plots==1) default.plot.title <- paste("Kaplan-Meyer plot of event",events.to.plot[1])
  if(!is.null(object)){
      plotTitle <- xpose.multiple.plot.title(object=object,
                                             plot.text = default.plot.title,
                                             main=main,
                                             ...)
  } else {
      if (is.null(main)){
          plotTitle <- NULL
      } else {
          if(!is.na(match(main,"Default"))) {
              plotTitle <- default.plot.title
              if (!is.null(subset)){
                  plotTitle <- paste(plotTitle,"\n[",subset,"]",sep="")
              }
          } else {
              plotTitle <- main
          }
      }
  }
  obj <- xpose.multiple.plot(plotList,plotTitle,...)
  return(obj)
}
