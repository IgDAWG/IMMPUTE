# Begin Functions #################################################################################
GetField <- function(x,Res) {
  Tmp <- unlist(strsplit(as.character(x),":"))
  if (length(Tmp)<2) {
    return(x)
  } else if (Res==1) {
    return(Tmp[1])
  } else if (Res > 1) {
    Out <- paste(Tmp[1:Res],collapse=":")
    return(Out)
  }
}

FixField <- function(x) {
  if(nchar(as.character(x))==1){
    return( paste("0",as.character(x),sep=""))
  } else { return(as.character(x)) }
}

percent <- function(x,y) {
  return(paste(round(x*100,y),"%",sep=""))
}

is.between <- function(x, a, b) {
  x > a & x <= b
}

colorBAR <- function(xS,xE,yS,yRange,Legend) {
  
  #Color Bar Legend
  ColmapL <- length(palette())
  Inc <- (0.04*yRange)
  yE <- yS + Inc
  Width <- (xE-xS)/ColmapL
  
  #Outline for color bar legend
  lines(x=c(xS,xE),y=c(yS,yS),lwd=3)
  lines(x=c(xS,xE),y=c(yE,yE),lwd=3)
  lines(x=c(xS,xS),y=c(yS,yE),lwd=3)
  lines(x=c(xE,xE),y=c(yS,yE),lwd=3)
  
  #Text for color bar legend
  text(xS,yE+(yRange*0.03),Legend[[2]],cex=1.0)
  text(xE,yE+(yRange*0.03),Legend[[1]],cex=1.0)
  text(mean(c(xS,xE) ),yS-(yRange*0.03),Legend[[3]],cex=0.75)
  
  for(z in 1:ColmapL) {
    xpts <- c(xS,xE,xE,xS)
    ypts <- c(yS,yS,yE,yE)
    polygon(xpts,ypts,col=1+(ColmapL-z),border=NA)
    xS <- xS + Width 
    xE <- xS + Width
  }
  
}

CheckMatch <- function(x,y) {
  tmp.lst <- list()
  for(z in 1:length(x)) {
    if(x %in% y) {
      tmp.lst[[z]] <- 1
    } else {
      tmp.lst[[z]] <- 0
    }
  }
  return(max(unlist(tmp.lst)))
}

CheckMatchv2 <- function(x,y) {
  tmp.lst <- list()
  for(z in 1:length(x)) {
    if(sum(x %in% y)>0) {
      tmp.lst[[z]] <- 1
    } else {
      tmp.lst[[z]] <- 0
    }
  }
  return(max(unlist(tmp.lst)))
}

CopyProb <- function(x) {
  Col <- grep("Prob",colnames(x))
  x[,Col+1] <- x[,Col]
  return(x)
}

Zipper <- function(x) {
  tmp.out <- list()
  for(y in 1:nrow(x)) {
    x.1 <- sort(unlist(strsplit(as.character(x[y,1]),"/")),decreasing=F)
    x.2 <- sort(unlist(strsplit(as.character(x[y,2]),"/")),decreasing=F)
    tmp <- expand.grid(x.1,x.2,stringsAsFactors=F)
    tmp <- matrix(data=paste(tmp[,1],tmp[,2],sep="^"),ncol=1)
    colnames(tmp) <- c("Genotype")
    tmp.out[[y]] <- tmp
  }
  tmp.out <- unique(do.call(rbind,tmp.out))
  return(tmp.out) 
}

HLA.Mirror <- function(x) {
  x.split <- sapply(x,strsplit,split="\\^")
  x.split <- do.call(rbind,x.split)
  rownames(x.split) <- NULL
  x.mirror <- matrix(paste(x.split[,2],x.split[,1],sep="^"),ncol=1)
  return(x.mirror)
}

# End Functions ###################################################################################