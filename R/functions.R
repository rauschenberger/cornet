
#--- Workhorse function --------------------------------------------------------

#' @export
#' @aliases cornet-package
#' @title
#' Combined regression
#' 
#' @description
#' Implements lasso and ridge regression for dichotomised outcomes.
#' Such outcomes are not naturally but artificially binary.
#' They indicate whether an underlying measurement is greater than a threshold.
#'
#' @param y
#' continuous outcome\strong{:}
#' vector of length \eqn{n}
#' 
#' @param cutoff
#' cut-off point for dichotomising outcome into classes\strong{:}
#' \emph{meaningful} value between \code{min(y)} and \code{max(y)}
#' 
#' @param X
#' features\strong{:}
#' numeric matrix with \eqn{n} rows (samples)
#' and \eqn{p} columns (variables)
#' 
#' @param alpha
#' elastic net mixing parameter\strong{:}
#' numeric between \eqn{0} (ridge) and \eqn{1} (lasso)
#' 
#' @param foldid
#' fold identifiers\strong{:}
#' vector with entries between \eqn{1} and \code{nfolds};
#' or \code{NULL} (balance)
#' 
#' @param nfolds
#' number of folds\strong{:}
#' integer between \eqn{3} and \eqn{n}
#' 
#' @param type.measure
#' loss function for binary classification\strong{:}
#' character \code{"deviance"}, \code{"mse"}, \code{"mae"},
#' or \code{"class"} (see \code{\link[glmnet]{cv.glmnet}})
#' 
#' @param pi
#' pi sequence\strong{:}
#' vector of increasing values in the unit interval;
#' or \code{NULL} (default sequence)
#' 
#' @param npi
#' number of \code{pi} values (weighting)
#' 
#' @param sigma
#' sigma sequence\strong{:}
#' vector of increasing positive values;
#' or \code{NULL} (default sequence)
#' 
#' @param nsigma
#' number of \code{sigma} values (scaling)
#' 
#' @param ...
#' further arguments passed to \code{\link[glmnet]{glmnet}}
#' 
#' @details
#' The argument \code{family} is unavailable, because
#' this function fits a \emph{gaussian} model for the numeric response,
#' and a \emph{binomial} model for the binary response.
#' 
#' Linear regression uses the loss function \code{"deviance"} (or \code{"mse"}),
#' but the loss is incomparable between linear and logistic regression.
#' 
#' The loss function \code{"auc"} is unavailable for internal cross-validation.
#' If at all, use \code{"auc"} for external cross-validation only.
#' 
#' @return
#' Returns an object of class \code{cornet}, a list with multiple slots:
#' \itemize{
#'    \item \code{gaussian}: fitted linear model, class \code{glmnet}
#'    \item \code{binomial}: fitted logistic model, class \code{glmnet}
#'    \item \code{sigma}: scaling parameters \code{sigma},
#'           vector of length \code{nsigma}
#'    \item \code{pi}: weighting parameters \code{pi},
#'           vector of length \code{npi}
#'    \item \code{cvm}: evaluation loss,
#'           matrix with \code{nsigma} rows and \code{npi} columns
#'    \item \code{sigma.min}: optimal scaling parameter,
#'           positive scalar
#'    \item \code{pi.min}: optimal weighting parameter,
#'           scalar in unit interval
#'    \item \code{cutoff}: threshold for dichotomisation
#' }
#' 
#' @seealso
#' Methods for objects of class \code{cornet} include
#' \code{\link[=coef.cornet]{coef}} and
#' \code{\link[=predict.cornet]{predict}}.
#' 
#' @references 
#' Armin Rauschenberger and Enrico Glaab (2023).
#' "Predicting artificial binary outcomes from high-dimensional data in biomedicine".
#' \emph{Manuscript in preparation}.
#' 
#' @examples
#' n <- 100; p <- 200
#' y <- rnorm(n)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' net <- cornet(y=y,cutoff=0,X=X)
#' net
#' 
cornet <- function(y,cutoff,X,alpha=1,npi=101,pi=NULL,nsigma=99,sigma=NULL,nfolds=10,foldid=NULL,type.measure="deviance",...){
  
  #--- temporary ---
  # alpha <- 1; cutoff <- 0; npi <- 101; pi <- NULL; nsigma <- 99; sigma <- NULL; nfolds <- 10;  foldid <- NULL; type.measure <- "deviance"; logistic <- TRUE
  test <- list()
  test$combined <- TRUE
  test$switch <- TRUE
  
  #--- checks ---
  n <- length(y)
  .check(x=y,type="vector")
  if(all(y %in% c(0,1))){warning("Binary response.",call.=FALSE)}
  .check(x=cutoff,type="scalar",min=min(y),max=max(y))
  if(length(y)!=nrow(X)){stop("Contradictory sample size.",call.=FALSE)}
  .check(x=X,type="matrix",dim=c(n,NA))
  .check(x=alpha,type="scalar",min=0,max=1)
  .check(x=npi,type="scalar",min=1)
  .check(x=pi,type="vector",min=0,max=1,null=TRUE)
  .check(x=nsigma,type="scalar",min=1)
  .check(x=sigma,type="vector",min=.Machine$double.eps,null=TRUE)
  .check(x=nfolds,type="scalar",min=3,max=n)
  .check(x=foldid,type="vector",dim=n,values=seq_len(nfolds),null=TRUE)
  .check(x=type.measure,type="string",values=c("deviance","mse","mae","class"))
  # never auc (min/max confusion)!
  if(!is.null(list(...)$family)){stop("Reserved argument \"family\".",call.=FALSE)}
  
  # binarisation
  z <- 1*(y > cutoff)
  
  # fold identifiers
  if(is.null(foldid)){
    foldid <- palasso:::.folds(y=z,nfolds=nfolds)
  }
  
  #--- model fitting ---
  fit <- list()
  fit$gaussian <- glmnet::glmnet(y=y,x=X,family="gaussian",alpha=alpha,...)
  fit$binomial <- glmnet::glmnet(y=z,x=X,family="binomial",alpha=alpha,...)
   
  #--- sigma sequence ---
  if(is.null(sigma)){
    fit$sigma <- exp(seq(from=log(0.05*stats::sd(y)), # decrease 0.05?
                  to=log(10*stats::sd(y)),length.out=nsigma))
  } else {
    fit$sigma <- sigma
    nsigma <- length(sigma)
  }
  
  #--- pi sequence ---
  if(is.null(pi)){
    fit$pi <- seq(from=0,to=1,length.out=npi)
  } else {
    fit$pi <- pi
    npi <- length(pi)
  }
  
  #--- tuning parameters ---
  lab.pi <- paste0("pi",seq_len(npi))
  lab.sigma <- paste0("si",seq_len(nsigma))
  names(fit$sigma) <- lab.sigma
  names(fit$pi) <- lab.pi
  
  #--- cross-validation ---
  pred <- list()
  pred$y <- matrix(data=NA,nrow=n,ncol=length(fit$gaussian$lambda))
  pred$z <- matrix(data=NA,nrow=n,ncol=length(fit$binomial$lambda))

  if(test$combined){
    dimnames <- list(NULL,lab.sigma,lab.pi)
    pred$combined <- array(data=NA,dim=c(n,nsigma,npi),dimnames=dimnames)
  }

  for(k in seq_len(nfolds)){

    y0 <- y[foldid!=k]
    y1 <- y[foldid==k]
    z0 <- z[foldid!=k]
    z1 <- z[foldid==k]
    X0 <- X[foldid!=k,,drop=FALSE]
    X1 <- X[foldid==k,,drop=FALSE]
    
    # linear regression
    net <- glmnet::glmnet(y=y0,x=X0,family="gaussian",alpha=alpha,...)
    temp_y <- stats::predict(object=net,newx=X1,type="response",s=fit$gaussian$lambda)
    pred$y[foldid==k,seq_len(ncol(temp_y))] <- temp_y
    cvm <- palasso:::.loss(y=y1,fit=temp_y,family="gaussian",type.measure="deviance")[[1]]
    y_hat <- temp_y[,which.min(cvm)]
    
    # logistic regression
    net <- glmnet::glmnet(y=z0,x=X0,family="binomial",alpha=alpha,...)
    temp_z <- stats::predict(object=net,newx=X1,type="response",s=fit$binomial$lambda)
    pred$z[foldid==k,seq_len(ncol(temp_z))] <- temp_z
    cvm <- palasso:::.loss(y=z1,fit=temp_z,family="binomial",type.measure=type.measure)[[1]]
    z_hat <- temp_z[,which.min(cvm)]
    
    # combined regression
    if(test$combined){
      for(i in seq_along(fit$sigma)){
        for(j in seq_along(fit$pi)){
          cont <- stats::pnorm(q=y_hat,mean=cutoff,sd=fit$sigma[i])
          pred$combined[foldid==k,i,j] <- fit$pi[j]*cont + (1-fit$pi[j])*z_hat
        }
      }
    }
    
  }
  
  #--- evaluation ---
  
  # linear loss
  fit$gaussian$cvm <- palasso:::.loss(y=y,fit=pred$y,family="gaussian",type.measure="deviance")[[1]]
  fit$gaussian$lambda.min <- fit$gaussian$lambda[which.min(fit$gaussian$cvm)]
  
  # logistic loss
  fit$binomial$cvm <- palasso:::.loss(y=z,fit=pred$z,family="binomial",type.measure=type.measure)[[1]]
  fit$binomial$lambda.min <- fit$binomial$lambda[which.min(fit$binomial$cvm)]

  # combined loss
  if(test$combined){
    dimnames <- list(lab.sigma,lab.pi)
    fit$cvm <- matrix(data=NA,nrow=nsigma,ncol=npi,dimnames=dimnames)
    for(i in seq_len(nsigma)){
      for(j in seq_len(npi)){
        fit$cvm[i,j] <- palasso:::.loss(y=z,fit=pred$combined[,i,j],family="binomial",type.measure=type.measure)[[1]]
      }
    }
    temp <- which(fit$cvm==min(fit$cvm),arr.ind=TRUE,useNames=TRUE)
    if(nrow(temp)>1){message("Multiple minima.");temp <- temp[1,,drop=FALSE]}
    fit$sigma.min <- fit$sigma[temp[1]]
    fit$pi.min <- fit$pi[temp[2]]
    if(fit$cvm[names(fit$sigma.min),names(fit$pi.min)]!=min(fit$cvm)){stop("Internal error.")}
  }
  
  if(test$switch){
    fit$inf <- fit$cvm
    fit$inf[,!fit$pi %in% c(0,1)] <- Inf
    temp <- which(fit$inf==min(fit$inf),arr.ind=TRUE,useNames=TRUE)
    if(nrow(temp)>1){message("Multiple minima.");temp <- temp[1,,drop=FALSE]}
    fit$inf.sigma.min <- fit$sigma[temp[1]]
    fit$inf.pi.min <- fit$pi[temp[2]]
    if(fit$inf[names(fit$inf.sigma.min),names(fit$inf.pi.min)]!=min(fit$inf)){stop("Internal error.")}
  }
  
  #--- return ---
  fit$cutoff <- cutoff
  fit$info <- list(type.measure=type.measure,
                   sd.min=fit$sigma[which.min(fit$cvm[,fit$pi==1])], # continuous only
                   sd.y=stats::sd(y),
                   "+"=sum(z==1),
                   "-"=sum(z==0),
                   table=table(z),
                   n=n,p=ncol(X),
                   test=as.data.frame(test))

  class(fit) <- "cornet"
  return(fit)
}

#--- Methods -------------------------------------------------------------------

#' @export
#' @title
#' Combined regression
#'
#' @description
#' Prints summary of cornet object.
#'
#' @param x
#' \link[cornet]{cornet} object
#' 
#' @param ...
#' further arguments (not applicable)
#' 
#' @return
#' Returns sample size \eqn{n},
#' number of covariates \eqn{p},
#' information on dichotomisation,
#' tuned scaling parameter (sigma),
#' tuned weighting parameter (pi),
#' and corresponding loss.
#'
#' @examples
#' n <- 100; p <- 200
#' y <- rnorm(n)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' net <- cornet(y=y,cutoff=0,X=X)
#' print(net)
#' 
print.cornet <- function(x,...){
  cat("cornet object:\n")
  cat(paste0("n = ",x$info$n,","," p = ",x$info$p),"\n")
  cat(paste0("z = I(y > ",signif(x$cutoff,digits=2),"): "))
  cat(paste0(x$info$"+","+"," vs ",x$info$"-","-"),"\n")
  cat(paste0("sigma.min = ",signif(x$sigma.min,digits=1)),"\n")
  cat(paste0("pi.min = ",round(x$pi.min,digits=2)),"\n")
  type <- x$info$type.measure
  value <- signif(x$cvm[names(x$sigma.min),names(x$pi.min)],digits=2)
  cat(paste0(type," = ",value))
  return(invisible(NULL))
}

#' @export
#' @title
#' Extract estimated coefficients
#'
#' @description
#' Extracts estimated coefficients from linear and logistic regression,
#' under the penalty parameter that minimises the cross-validated loss.
#'
#' @param object
#' \link[cornet]{cornet} object
#' 
#' @param ...
#' further arguments (not applicable)
#' 
#' @return
#' This function returns a matrix with \eqn{n} rows and two columns,
#' where \eqn{n} is the sample size. It includes the estimated coefficients
#' from linear regression (1st column: \code{"beta"})
#' and logistic regression (2nd column: \code{"gamma"}).
#'
#' @examples
#' n <- 100; p <- 200
#' y <- rnorm(n)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' net <- cornet(y=y,cutoff=0,X=X)
#' coef(net)
#' 
coef.cornet <- function(object,...){
  
  if(length(list(...))!=0){warning("Ignoring arguments.")}
  
  s <- object$gaussian$lambda.min
  #beta <- glmnet::coef.glmnet(object=object$gaussian,s=s)
  beta <- stats::coef(object=object$gaussian,s=s)
  
  s <- object$binomial$lambda.min
  #gamma <- glmnet::coef.glmnet(object=object$binomial,s=s)
  gamma <- stats::coef(object=object$binomial,s=s)
  
  coef <- cbind(beta,gamma)
  colnames(coef) <- c("beta","gamma")
  return(coef)
}

#' @export
#' @title
#' Plot loss matrix
#'
#' @description
#' Plots the loss for different combinations of
#' scaling (sigma) and weighting (pi) parameters.
#'
#' @param x
#' \link[cornet]{cornet} object
#' 
#' @param ...
#' further arguments (not applicable)
#' 
#' @return
#' This function plots the evaluation loss (\code{cvm}).
#' Whereas the matrix has sigma in the rows, and pi in the columns,
#' the plot has sigma on the \eqn{x}-axis, and pi on the \eqn{y}-axis.
#' For all combinations of sigma and pi, the colour indicates the loss.
#' If the R package \code{RColorBrewer} is installed,
#' blue represents low. Otherwise, red represents low.
#' White always represents high.
#' 
#' @examples
#' n <- 100; p <- 200
#' y <- rnorm(n)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' net <- cornet(y=y,cutoff=0,X=X)
#' plot(net)
#' 
plot.cornet <- function(x,...){
  
  if(length(list(...))!=0){warning("Ignoring arguments.")}

  k <- 100
  levels <- stats::quantile(x$cvm,probs=seq(from=0,to=1,length.out=k+1))
  
  # colours
  if("RColorBrewer" %in% .packages(all.available=TRUE)){
    pal <- rev(c("white",RColorBrewer::brewer.pal(n=9,name="Greys")))
    col <- grDevices::colorRampPalette(colors=pal)(k)
  } else {
    col <- grDevices::heat.colors(n=k)
  }
  #col <- grDevices::grey(level=seq(from=0,to=1,length.out=k))
  
  nsigma <- length(x$sigma)
  npi <- length(x$pi)
  
  graphics::plot.new()
  graphics::par(xaxs="i",yaxs="i")
  graphics::plot.window(xlim=c(1-0.5,nsigma+0.5),ylim=c(1-0.5,npi+0.5))
  
  graphics::title(xlab=expression(sigma),ylab=expression(pi),cex.lab=1)
  #graphics::.filled.contour(x=seq_along(x$sigma),y=seq_along(x$pi),z=x$cvm,levels=levels,col=col)
  graphics::image(x=seq_along(x$sigma),y=seq_along(x$pi),z=x$cvm,breaks=levels,col=col,add=TRUE)
  graphics::box()
  
  ssigma <- which(x$sigma %in% x$sigma.min)
  spi <- which(x$pi %in% x$pi.min)
  
  if(length(ssigma)==1 & length(spi)==1){
    # axes with labels for tuned parameters
    graphics::axis(side=1,at=c(1,ssigma,nsigma),labels=signif(x$sigma[c(1,ssigma,nsigma)],digits=2))
    graphics::axis(side=2,at=c(1,spi,npi),labels=signif(x$pi[c(1,spi,npi)],digits=2))
    # point for tuned parameters
    graphics::points(x=ssigma,y=spi,pch=4,col="white",cex=1)
  } else {
    # axes with standard labels
    at <- seq(from=1,to=nsigma,length.out=5)
    graphics::axis(side=1,at=at,labels=signif(x$sigma,digits=2)[at])
    at <- seq(from=1,to=nsigma,length.out=5)
    graphics::axis(side=2,at=at,labels=signif(x$pi,digits=2)[at])
    # points for selected parameters
    isigma <- sapply(x$sigma.min,function(y) which(x$sigma==y))
    ipi <- sapply(x$pi.min,function(y) which(x$pi==y))
    graphics::points(x=isigma,y=ipi,pch=4,col="white",cex=1)
  }
  
}

#' @export
#' @title
#' Predict binary outcome
#'
#' @description
#' Predicts the binary outcome with linear, logistic, and combined regression.
#' 
#' @details
#' For linear regression, this function tentatively transforms
#' the predicted values to predicted probabilities,
#' using a Gaussian distribution with a fixed mean (threshold)
#' and a fixed variance (estimated variance of the numeric outcome).
#' 
#' @param object
#' \link[cornet]{cornet} object
#' 
#' @param newx
#' covariates\strong{:}
#' numeric matrix with \eqn{n} rows (samples)
#' and \eqn{p} columns (variables)
#' 
#' @param type
#' \code{"probability"}, \code{"odds"}, \code{"log-odds"}
#' 
#' @param ...
#' further arguments (not applicable)
#' 
#' @examples
#' n <- 100; p <- 200
#' y <- rnorm(n)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' net <- cornet(y=y,cutoff=0,X=X)
#' predict(net,newx=X)
#' 
predict.cornet <- function(object,newx,type="probability",...){
  
  if(length(list(...))!=0){warning("Ignoring arguments.")}
  
  x <- object; rm(object)
  
  test <- x$info$test
  
  .check(x=newx,type="matrix")
  .check(x=type,type="string",values=c("probability","odds","log-odds"))
  
  prob <- list()
  
  # intercept
  prob$intercept <- as.numeric(stats::predict(object=x$binomial,
                        newx=newx,type="response")[,"s0"])
  
  # linear and logistic
  link <- as.numeric(stats::predict(object=x$gaussian,
                  newx=newx,s=x$gaussian$lambda.min,type="response"))
  prob$gaussian <- stats::pnorm(q=link,mean=x$cutoff,sd=x$info$sd.min) # continuous only
  prob$binomial <- as.numeric(stats::predict(object=x$binomial,
                  newx=newx,s=x$binomial$lambda.min,type="response"))
  
  # combined
  if(test$combined){
    cont <- stats::pnorm(q=link,mean=x$cutoff,sd=x$sigma.min)
    prob$combined <- x$pi.min*cont + (1-x$pi.min)*prob$binomial
  }
  
  if(test$switch){
    cont <- stats::pnorm(q=link,mean=x$cutoff,sd=x$inf.sigma.min)
    prob$switch <- x$inf.pi.min*cont + (1-x$inf.pi.min)*prob$binomial
  }
  
  # consistency tests
  lapply(X=prob,FUN=function(p) .check(x=p,type="vector",min=0,max=1))
  .equal(link>x$cutoff,prob$gaussian>0.5)

  # transformation
  if(type=="probability"){
    frame <- prob
  } else if(type=="odds"){
    frame <- lapply(X=prob,FUN=function(x) x/(1-x))
  } else if(type=="log-odds"){
    frame <- lapply(X=prob,FUN=function(x) log(x/(1-x)))
  } else {
    stop("Invalid type.",call.=FALSE)
  }
  
  return(as.data.frame(frame))
}

#--- Internal functions --------------------------------------------------------

#' @title
#' Equality
#'
#' @description
#' Verifies whether two or more arguments are identical.
#'
#' @param ...
#' scalars, vectors, or matrices of equal dimensions
#' 
#' @param na.rm
#' remove missing values\strong{:}
#' logical
#' 
#' @examples 
#' cornet:::.equal(1,1,1)
#' 
.equal <- function(...,na.rm=FALSE){
  list <- list(...)
  cond <- vapply(X=list,
                 FUN=function(x) all(x==list[[1]],na.rm=na.rm),
                 FUN.VALUE=logical(1))
  if(any(!cond)){
    stop("Unequal elements.",call.=FALSE)
  }
  return(invisible(NULL))
}

#' @title
#' Arguments
#'
#' @description
#' Verifies whether an argument matches formal requirements.
#'
#' @param x
#' argument
#' 
#' @param type
#' character \code{"string"}, \code{"scalar"}, \code{"vector"}, \code{"matrix"}
#' 
#' @param dim
#' vector/matrix dimensionality\strong{:}
#' integer scalar/vector
#' 
#' @param miss
#' accept missing values\strong{:}
#' logical
#' 
#' @param min
#' lower limit\strong{:}
#' numeric
#' 
#' @param max
#' upper limit\strong{:}
#' numeric
#' 
#' @param values
#' only accept specific values\strong{:}
#' vector 
#' 
#' @param inf
#' accept infinite (\code{Inf} or \code{-Inf}) values\strong{:}
#' logical
#' 
#' @param null
#' accept \code{NULL}\strong{:}
#' logical
#'
#' @examples
#' cornet:::.check(0.5,type="scalar",min=0,max=1)
#' 
.check <- function(x,type,dim=NULL,miss=FALSE,min=NULL,max=NULL,values=NULL,inf=FALSE,null=FALSE){
  name <- deparse(substitute(x))
  if(null && is.null(x)){
    return(invisible(NULL)) 
  }
  if(type=="string"){
    cond <- is.vector(x) & is.character(x) & length(x)==1
  } else if(type=="scalar"){
    cond <- is.vector(x) & is.numeric(x) & length(x)==1
  } else if(type=="vector"){
    cond <- is.vector(x) & is.numeric(x)
  } else if(type=="matrix"){
    cond <- is.matrix(x) & is.numeric(x)
  } else {
    warning("Unknown type.")
  }
  if(!cond){
    stop(paste0("Argument \"",name,"\" does not match formal requirements."),call.=FALSE)
  }
  if(!is.null(dim) && length(dim)==1 && length(x)!=dim){
      stop(paste0("Argument \"",name,"\" has invalid length."),call.=FALSE)
  }
  if(!is.null(dim) && length(dim)>1 && any(dim(x)!=dim,na.rm=TRUE)){
      stop(paste0("Argument \"",name,"\" has invalid dimensions."),call.=FALSE)
  }
  
  #   } else if(length(dim)==2){
  #     if(!is.na(dim[1]) && nrow(x)!=dim[1]){
  #       stop(paste0("Argument \"",name,"\" has invalid row number."),call.=FALSE)
  #     }
  #     if(!is.na(dim[2]) && ncol(x)!=dim[2]){
  #       stop(paste0("Argument \"",name,"\" has invalid column number."),call.=FALSE)
  #     }
  #   } else {
  #     
  #   }
  # }
  
  # if(!is.null(length) && length(x)!=length){
  #   stop(paste0("Argument \"",name,"\" has invalid length."),call.=FALSE)
  # }
  # if(!is.null(nrow) && nrow(x)!=nrow){
  #   stop(paste0("Argument \"",name,"\" has invalid row number."),call.=FALSE)
  # }
  # if(!is.null(ncol) && ncol(x)!=ncol){
  #   stop(paste0("Argument \"",name,"\" has invalid column number."),call.=FALSE)
  # }
  
  if(!miss && any(is.na(x))){
    stop(paste0("Argument \"",name,"\" contains missing values."),call.=FALSE)
  }
  if(!is.null(min) && any(x<min,na.rm=TRUE)){
    stop(paste0("expecting ",name," >= ",min),call.=FALSE)
  }
  if(!is.null(max) && any(x>max,na.rm=TRUE)){
    stop(paste0("expecting ",name," <= ",max),call.=FALSE)
  }
  if(!is.null(values) && any(!x %in% values,na.rm=TRUE)){
    stop(paste0("Argument \"",name,"\" contains invalid values."),call.=FALSE)
  }
  if(!inf && any(is.infinite(values))){
    stop(paste0("Argument \"",name,"\ contains infinite values."),call.=FALSE)
  }
  return(invisible(NULL))
}

#--- Application ---------------------------------------------------------------

#' @export
#' @title
#' Performance measurement
#'
#' @description
#' Compares models for a continuous response with a cut-off value.
#' 
#' @details
#' Computes the cross-validated loss of logistic and combined regression.
#' 
#' @inheritParams  cornet
#' 
#' @param nfolds.ext
#' number of external folds
#' 
#' @param foldid.int
#' number of internal folds
#' 
#' @param foldid.ext
#' external fold identifiers\strong{:}
#' vector of length \eqn{n} with entries
#' between \eqn{1} and \code{nfolds.ext};
#' or \code{NULL}
#' 
#' @param nfolds.int
#' internal fold identifiers\strong{:}
#' vector of length \eqn{n} with entries
#' between \eqn{1} and \code{nfolds.int};
#' or \code{NULL}
#' 
#' @param rf
#' comparison with random forest\strong{:}
#' logical
#' 
#' @param svm
#' comparison with support vector machine\strong{:}
#' logical
#' 
#' @param ... further arguments passed to
#' \code{\link[cornet]{cornet}} or \code{\link[glmnet]{glmnet}}
#' 
#' @examples
#' \dontshow{#n <- 50; p <- 20
#' #y <- rnorm(n)
#' #X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' #loss <- cv.cornet(y=y,cutoff=0,X=X,nfolds.ext=2)
#' #loss}
#' \dontrun{n <- 100; p <- 200
#' y <- rnorm(n)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' loss <- cv.cornet(y=y,cutoff=0,X=X,rf=TRUE,svm=FALSE,knn=TRUE)
#' loss}
#' 
cv.cornet <- function(y,cutoff,X,alpha=1,nfolds.ext=5,nfolds.int=10,foldid.ext=NULL,foldid.int=NULL,type.measure="deviance",rf=FALSE,svm=FALSE,knn=FALSE,...){
  
  if(rf){
    if(!"randomForest" %in% rownames(utils::installed.packages())){
      stop("Requires package 'randomForest'.")
    }
  }
  
  if(svm){
    if(!"e1071" %in% rownames(utils::installed.packages())){
      stop("Requires package 'e1071'.")
    }
  }
  
  if(knn){
    if(!"class" %in% rownames(utils::installed.packages())){
      stop("Requires package 'class'.")
    }
  }
  
  z <- 1*(y > cutoff)
  if(is.null(foldid.ext)){
    foldid.ext <- palasso:::.folds(y=z,nfolds=nfolds.ext)
  } else {
    nfolds.ext <- max(foldid.ext)
  }
  
  #--- cross-validated loss ---
  
  cols <- c("intercept","binomial","gaussian","combined","switch","rf"[rf],"svm"[svm],"knn"[knn])
  pred <- matrix(data=NA,nrow=length(y),ncol=length(cols),
                 dimnames=list(NULL,cols))
  
  for(i in seq_len(nfolds.ext)){
    
    y0 <- y[foldid.ext!=i]
    z0 <- z[foldid.ext!=i]
    X0 <- X[foldid.ext!=i,]
    X1 <- X[foldid.ext==i,]
    
    if(is.null(foldid.int)){
      foldid <- palasso:::.folds(y=z0,nfolds=nfolds.int)
    } else {
      foldid <- foldid.int[foldid.ext!=i]
    }
    
    fit <- cornet(y=y0,cutoff=cutoff,X=X0,alpha=alpha,type.measure=type.measure,foldid=foldid,...)
    tryCatch(expr=plot.cornet(fit),error=function(x) NULL)
    temp <- predict.cornet(fit,newx=X1)
    if(any(temp<0|temp>1)){stop("Outside unit interval.",call.=FALSE)}
    model <- colnames(temp) # was "colnames(pred)"
    for(j in seq_along(model)){
      pred[foldid.ext==i,model[j]] <- temp[[model[j]]]
    }
    
    if(rf){
      object <- randomForest::randomForest(x=X0,y=as.factor(z0),norm.votes=TRUE)
      pred[foldid.ext==i,"rf"] <- stats::predict(object,newdata=X1,type="prob")[,2]
      #colnames(X0) <- paste0("x",1:ncol(X0))
      #rf <- caret::train(x=X0,y=as.factor(z0),method="rf")
    }
    
    if(svm){
      #object <- e1071::svm(x=X0,y=as.factor(z0),probability=TRUE)
      #object <- e1071::best.svm(x=X0,y=as.factor(z0),probability=TRUE,gamma=2^(-1:1),cost=2^(2:4))
      #pred[foldid.ext==i,"svm"] <- attributes(stats::predict(object,newdata=X1,probability=TRUE))$probabilities[,2]
    }
    
    if(knn){
      #object <- e1071::tune.knn(x=X0,y=as.factor(z0),probability=TRUE,k=1:10)
      #t <- e1071::knn(train=X0,test=X1,cl=as.factor(z0),k=object$best.parameters,prob=TRUE)
      #pred[foldid.ext==i,"svm"] <- attributes(stats::predict(object,newdata=X1,probability=TRUE))$probabilities[,2]
      #colnames(X0) <- colnames(X1) <- paste0("x",seq_len(p))
      #knn <- caret::train(x=X0,y=as.factor(ifelse(z0==1,"one","zero")),method="knn",trControl=caret::trainControl(classProbs=TRUE))
      #predict(knn,newdata=X1,classProbs=TRUE)
      cvm <- numeric()
      for(k in seq_len(50)){
        temp  <- class::knn.cv(train=X0,cl=as.factor(z0),k=k)
        cvm[k] <- mean(z0!=temp)
      }
      temp <- class::knn(train=X0,test=X1,cl=as.factor(z0),k=which.min(cvm),prob=TRUE)
      pred[foldid.ext==i,"knn"] <- ifelse(temp==1,attributes(temp)$prob,1-attributes(temp)$prob)
    }
    
  }
  
  type <- c("deviance","class","mse","mae","auc")
  loss <- lapply(X=type,FUN=function(x) palasso:::.loss(y=z[foldid.ext!=0],fit=pred[foldid.ext!=0,],family="binomial",type.measure=x,foldid=foldid.ext[foldid.ext!=0])[[1]])
  names(loss) <- type
  
  if("MLmetrics" %in% installed.packages()[,1]){
    loss$prauc <- apply(pred[foldid.ext!=0,],2,function(x) MLmetrics::PRAUC(y_pred=x,y_true=z[foldid.ext!=0]))
  }

  # if(FALSE){
  #   #--- deviance residuals --- (removed during revision)
  #   
  #   # squared deviance residuals
  #   limit <- 1e-05
  #   pred[pred < limit] <- limit
  #   pred[pred > 1 - limit] <- 1 - limit
  #   res <- -2 * (z * log(pred) + (1 - z) * log(1 - pred))
  #   rxs <- res[,"binomial"]
  #   rys <- res[,"combined"]
  #   
  #   # residual increase/decrease
  #   loss$resid.factor <- stats::median((rys-rxs)/rxs)
  #   
  #   # paired test for each fold
  #   loss$resid.pvalue <- numeric()
  #   for(i in seq_len(nfolds.ext)){
  #     cond <- foldid.ext==i
  #     loss$resid.pvalue[i] <- stats::wilcox.test(x=rxs[cond],y=rys[cond],
  #                                                paired=TRUE,alternative="greater")$p.value
  #   }
  # }
  
  loss <- lapply(loss,function(x) signif(x,digits=6)) # trial stability
  
  return(loss)
  
}

#' @title
#' Single-split test
#'
#' @description
#' Compares models for a continuous response with a cut-off value.
#' 
#' @details
#' Splits samples into \eqn{80} percent for training
#' and \eqn{20} percent for testing,
#' calculates squared deviance residuals of logistic and combined regression,
#' conducts the paired one-sided Wilcoxon signed rank test,
#' and returns the \eqn{p}-value.
#' For the multi-split test,
#' use the median \eqn{p}-value from \eqn{50} single-split tests
#' (van de Wiel 2009).
#' 
#' @inheritParams cornet
#' 
#' @examples
#' n <- 100; p <- 200
#' y <- rnorm(n)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' cornet:::.test(y=y,cutoff=0,X=X)
#' 
.test <- function(y,cutoff,X,alpha=1,type.measure="deviance"){
  
  z <- 1*(y > cutoff)
  fold <- palasso:::.folds(y=z,nfolds=5)
  fold <- ifelse(fold==1,1,0)
  
  fit <- cornet(y=y[fold==0],cutoff=cutoff,X=X[fold==0,],
                        alpha=alpha,type.measure=type.measure)
  tryCatch(expr=plot.cornet(fit),error=function(x) NULL)
  pred <- predict.cornet(fit,newx=X[fold==1,])
  if(any(pred<0|pred>1)){stop("Outside unit interval.",call.=FALSE)}
  
  limit <- 1e-05
  pred[pred < limit] <- limit
  pred[pred > 1 - limit] <- 1 - limit
  res <- -2 * (z[fold==1] * log(pred) + (1 - z[fold==1]) * log(1 - pred))
  #pvalue <- stats::wilcox.test(x=res[,"binomial"],y=res[,"combined"],paired=TRUE,alternative="greater")$p.value
  
  pvalue <- list()
  pvalue$log <- stats::wilcox.test(x=res[,"combined"],y=res[,"binomial"],paired=TRUE,alternative="less")$p.value
  pvalue$lin <- stats::wilcox.test(x=res[,"combined"],y=res[,"gaussian"],paired=TRUE,alternative="less")$p.value
  
  ##equality logistic deviance
  #colMeans(res)
  #palasso:::.loss(y=z[fold==1],fit=pred,family="binomial",type.measure="deviance")
  
  return(pvalue)
}

#' @title
#' Data simulation
#'
#' @description
#' Simulates data for unit tests
#' 
#' @param n
#' sample size\strong{:}
#' positive integer
#' 
#' @param p
#' covariate space\strong{:}
#' positive integer
#' 
#' @param cor
#' correlation coefficient \strong{:}
#' numeric between \eqn{0} and \eqn{1}
#' 
#' @param prob
#' effect proportion\strong{:}
#' numeric between \eqn{0} and \eqn{1}
#' 
#' @param sd
#' standard deviation\strong{:}
#' positive numeric
#' 
#' @param exp
#' exponent\strong{:}
#' positive numeric
#' 
#' @param frac
#' class proportion\strong{:}
#' numeric between \eqn{0} and \eqn{1}
#' 
#' @details
#' For simulating correlated features (\code{cor}\eqn{>0}),
#' this function requires the R package MASS
#' (see \code{\link[MASS]{mvrnorm}}).
#' 
#' @return
#' Returns invisible list with elements \code{y} and \code{X}.
#' 
#' @examples
#' data <- cornet:::.simulate(n=10,p=20)
#' names(data)
#' 
.simulate <- function(n , p, cor = 0, prob = 0.1, sd = 1, exp = 1, frac = 1) {
  
  .check(x=n,type="scalar",min=1)
  .check(x=p,type="scalar",min=1)
  .check(x=cor,type="scalar",min=0,max=1)
  .check(x=prob,type="scalar",min=0,max=1)
  .check(x=sd,type="scalar",min=0)
  .check(x=exp,type="scalar",min=0)
  .check(x=frac,type="scalar",min=0,max=1)
  
  mu <- rep(x = 0, times = p)
  Sigma <- matrix(data = NA, nrow = p, ncol = p)
  Sigma <- cor^abs(row(Sigma) - col(Sigma))
  diag(Sigma) <- 1
  
  #warning("DEACTIVE THIS!")
  #Sigma <- round(Sigma,digits=5)
  #X <- mvtnorm::rmvnorm(n = n, mean = mu, sigma = Sigma,method = "chol")
  #X <- round(X,digits=5)
  
  #if (requireNamespace("MASS", quietly = TRUE)) {
  #    X <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma) # not reproducible
  if (requireNamespace("mvtnorm", quietly = TRUE)) {
      X <- mvtnorm::rmvnorm(n = n, mean = mu, sigma = Sigma)    
  } else {
      if(cor!=0){stop("Correlation requires R package mvtnorm!",call.=FALSE)}
      X <- vapply(X = mu, FUN = function(x) stats::rnorm(n = n,mean = x),
                  FUN.VALUE = numeric(n))
  }
  
  beta <- stats::runif(n=p)*stats::rbinom(n = p,size = 1, prob = prob)
  if(FALSE){ # non-linear effects
    warning("temporary: non-linear effects")
    lambda <- sample(x=c(0,0.5,1,2),replace=TRUE,prob=c(0.25,0.25,0.25,0.25),size=p)
    tau <- sample(x=c(-2,-1,0,1,2),replace=TRUE,prob=c(0.2,0.2,0.2,0.2,0.2),size=p)
    Z <- sign(X-tau)*abs(X-tau)^rep(lambda,each=n)
  } else {
    Z <- X
  }
  mean <- Z %*% beta
  y <- stats::rnorm(n = n, mean = mean, sd = sd)
  y <- sign(y) * abs(y)^exp # departure from normality
  
  cutoff <- stats::quantile(y,probs=frac)
  
  list <- list(y = y, X = X, cutoff = cutoff)
  list <- lapply(list,function(x) signif(x,digits=6)) # trial stability
  # or use signif?
  return(invisible(list))
}

#--- Legacy --------------------------------------------------------------------

# multi.split <- function(...,n=50){
#  median(replicate(n=n,expr=cornet:::.test(...)))
# }

# @param prob
# (approximate) proportion of causal covariates\strong{:}
# numeric between \eqn{0} and \eqn{1}
# 
# @param fac
# noise factor\strong{:}
# positive real number
# 
#.simulate <- function(n,p,prob=0.2,fac=1){
#  beta <- stats::rnorm(n=p)
#  cond <- stats::rbinom(n=p,size=1,prob=prob)
#  beta[cond==0] <- 0
#  X <- matrix(stats::rnorm(n=n*p),nrow=n,ncol=p)
#  mean <- X %*% beta
#  y <- stats::rnorm(n=n,mean=mean,sd=fac*stats::sd(mean))
#  return(invisible(list(y=y,X=X)))
#}

# # Import this function from the palasso package.
# .loss <- function (y,fit,family,type.measure,foldid=NULL){
#   if (!is.list(fit)) {
#     fit <- list(fit)
#   }
#   loss <- list()
#   for (i in seq_along(fit)) {
#     if (is.vector(fit[[i]])) {
#       fit[[i]] <- as.matrix(fit[[i]])
#     }
#     if (is.null(foldid) & (family == "cox" | type.measure == 
#                            "auc")) {
#       stop("Missing foldid.", call. = FALSE)
#     }
#     if (family == "gaussian") {
#       if (type.measure %in% c("deviance", "mse")) {
#         loss[[i]] <- apply(X = fit[[i]], MARGIN = 2, 
#                            FUN = function(x) mean((x - y)^2))
#       }
#       else if (type.measure == "mae") {
#         loss[[i]] <- apply(X = fit[[i]], MARGIN = 2, 
#                            FUN = function(x) mean(abs(x - y)))
#       }
#       else {
#         stop("Invalid type measure.", call. = FALSE)
#       }
#     }
#     else if (family == "binomial") {
#       if (type.measure == "deviance") {
#         limit <- 1e-05
#         fit[[i]][fit[[i]] < limit] <- limit
#         fit[[i]][fit[[i]] > 1 - limit] <- 1 - limit
#         loss[[i]] <- apply(X = fit[[i]], MARGIN = 2, 
#                            FUN = function(x) mean(-2 * (y * log(x) + (1 - 
#                                                                         y) * log(1 - x))))
#       }
#       else if (type.measure == "mse") {
#         loss[[i]] <- apply(X = fit[[i]], MARGIN = 2, 
#                            FUN = function(x) 2 * mean((x - y)^2))
#       }
#       else if (type.measure == "mae") {
#         loss[[i]] <- apply(X = fit[[i]], MARGIN = 2, 
#                            FUN = function(x) 2 * mean(abs(x - y)))
#       }
#       else if (type.measure == "class") {
#         loss[[i]] <- apply(X = fit[[i]], MARGIN = 2, 
#                            FUN = function(x) mean(abs(round(x) - y)))
#       }
#       else if (type.measure == "auc") {
#         weights <- table(foldid)
#         cvraw <- matrix(data = NA, nrow = length(weights), 
#                         ncol = ncol(fit[[i]])) # typo in palasso package !
#         for (k in seq_along(weights)) {
#           cvraw[k, ] <- apply(X = fit[[i]], MARGIN = 2, 
#                               FUN = function(x) glmnet::auc(y = y[foldid == 
#                                                                     k], prob = x[foldid == k]))
#         }
#         loss[[i]] <- apply(X = cvraw, MARGIN = 2, FUN = function(x) stats::weighted.mean(x = x, 
#                                                                                          w = weights, na.rm = TRUE))
#         names(loss[[i]]) <- colnames(fit[[i]]) # typo in palasso package!
#       }
#       else {
#         stop("Invalid type measure.", call. = FALSE)
#       }
#     }
#     else if (family == "poisson") {
#       if (type.measure == "deviance") {
#         loss[[i]] <- apply(X = fit[[i]], MARGIN = 2, 
#                            FUN = function(x) mean(2 * (ifelse(y == 0, 
#                                                               0, y * log(y/x)) - y + x), na.rm = TRUE))
#       }
#       else if (type.measure == "mse") {
#         loss[[i]] <- apply(X = fit[[i]], MARGIN = 2, 
#                            FUN = function(x) mean((x - y)^2))
#       }
#       else if (type.measure == "mae") {
#         loss[[i]] <- apply(X = fit[[i]], MARGIN = 2, 
#                            FUN = function(x) mean(abs(x - y)))
#       }
#       else {
#         stop("Invalid type measure.", call. = FALSE)
#       }
#     }
#     else if (family == "cox") {
#       if (type.measure == "deviance") {
#         weights <- tapply(X = y[, "status"], INDEX = foldid, 
#                           FUN = sum)
#         loss[[i]] <- apply(X = fit[[i]], MARGIN = 2, 
#                            FUN = function(x) stats::weighted.mean(x = x/weights, 
#                                                                   w = weights, na.rm = TRUE))
#       }
#       else {
#         stop("Invalid type measure.", call. = FALSE)
#       }
#     }
#     else {
#       stop("Invalid family.", call. = FALSE)
#     }
#     if (sum(diff(is.na(loss[[i]]))) == 1) {
#       loss[[i]] <- loss[[i]][!is.na(loss[[i]])]
#     }
#   }
#   return(loss)
# }
# 
# # Import this function from the palasso package.
# .folds <- function (y, nfolds, foldid = NULL){
#   if(!is.null(foldid)){
#     return(foldid)
#   }
#   #if (survival::is.Surv(y)){ # active in palasso
#   #  y <- y[, "status"] # active in palasso
#   #} # active in palasso
#   if(all(y %in% c(0, 1))){
#     foldid <- rep(x = NA, times = length(y))
#     foldid[y == 0] <- sample(x = rep(x = seq_len(nfolds), 
#                                      length.out = sum(y == 0)))
#     foldid[y == 1] <- sample(x = rep(x = seq_len(nfolds), 
#                                      length.out = sum(y == 1)))
#   } else {
#     foldid <- sample(x = rep(x = seq_len(nfolds), length.out = length(y)))
#   }
#   return(foldid)
# }
