combineIndices <- function(...){

  args <- list(...)
  nargs <- length(args)
  if(nargs==1) return(args[[1]])
  largs <- sapply(args,length)
  if(length(unique(largs))>1) stop("arguments have different length")
  i <- args[[1]]
  for(ii in 2:nargs){

    #print(i)
    I <- max(i)
    j <- args[[ii]]
    J <- max(j)
    i <- i + I*(j-1)
  }
  i
}

interactIndices <- function(...){
  i <- combineIndices(...)
  match(i,sort(unique(i)))
}


rowsum2 <- function(x,...) UseMethod("rowsum2")

rowsum2.default <- function(x,...,default=NA){

  indices <- list(...)
  idims <- sapply(indices,max)
  i <- combineIndices(...)
  j <- match(i,unique(i))
  tmp <- rowsum(x,j)
  res <- array(default,idims)
  res[i] <- tmp[j]
  res
}

rowsum2.matrix <- function(x,...,default=NA){

  indices <- list(...)
  idims <- sapply(indices,max)
  i <- combineIndices(...)
  j <- match(i,unique(i))
  tmp <- rowsum(x,j)
  res <- matrix(default,nrow=prod(idims),ncol=ncol(x))
  res[i,] <- tmp[j,]
  dim(res) <- c(idims,ncol(x))
  res
}