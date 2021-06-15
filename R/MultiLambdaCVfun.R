#####################################################################################################
#
#
#                                     SET UP PARALLEL COMPUTING + SVD based SINGLE lambda CV
#
#
#####################################################################################################

setupParallel <- function(ncpus=2, sourcefile=NULL,sourcelibraries=c("multiridge","survival","pROC","risksetROC")){
  if(ncpus < 2) {
    message("ERROR: use ncpus >= 2")
    return(NULL)
  }
  ex <- requireNamespace("snowfall")
  if(!ex) print("Please install snowfall package for allowing parallel computation")
  sfInit(parallel=TRUE,cpus=ncpus,slaveOutfile="test.txt")
  if(!is.null(sourcefile)) invisible(capture.output(sfSource(sourcefile)))
  sfLibrary("snowfall", character.only=TRUE)
  if(!is.null(sourcelibraries)) for(k in 1:length(sourcelibraries)) {
    #k<-1
    sl <- sourcelibraries[k]
    inst <- try(find.package(sl))
    if(class(sl) =="try-error"){
      message(paste("source library",sl,"not installed. Install when you wish to use it for parallel computing"))
    } else sfLibrary(sl, character.only=TRUE)
  }
  print("Set-up ready")
  #runsetup <<- TRUE
  #assign("runsetup", TRUE, envir = baseenv())
  #pkg.globals$runsetup <- TRUE
  #return()
}

#fast CV that does not depend on penalized package
fastCV2 <- function(XXblocks,Y,X1=NULL,kfold=10,
                    intercept=ifelse(class(Y)=="Surv", FALSE, TRUE),parallel=FALSE,fixedfolds = TRUE,
                    model=NULL, eps=1e-10, reltol=0.5, lambdamax = 10^6,traceCV=TRUE){
  # XXblocks=XXbl;Y=resp;kfold=10;fixedfolds = TRUE; intercept <-TRUE;
  # parallel=FALSE;model=NULL; eps=1e-10; lambdamax = 10^6
  if(is.null(model)){
    if(class(Y)=="Surv") model <- "cox" else {
      model <- ifelse(length(unique(Y)) ==2, "logistic","linear")
    }
  }

  if(parallel){if(sfCpus()==1) {
    message("ERROR: please run the setupParallel() function first")
    #message("Use sourcelibraries=c(\"penalized\")")
    return(NULL)
  }}
  if(class(Y) == "Surv" & model != "cox") {
    print("Response is of class survival. Switching to model=\'survival\' ")
    model <- "cox"
  }
  leftout <- CVfolds(Y=Y,model=model, kfold=kfold,nrepeat=1, fixedfolds=fixedfolds)
  #fold <- leftout
  fastCVX2 <- function(XX){
    #XX <- XXblocks[[1]]
    #if(class(kfold) == "list") {
    #converts leftout represented as list to fold represented as vector (as required by optL2() function)
    XXbl <- list(XX)
    ol1 <- try(optLambdas(penaltiesinit = NULL, XXblocks=XXbl, Y=Y, X1=X1, folds=leftout, intercept=intercept, model=model, reltol=reltol,traceCV=traceCV,
                          fixedseed=fixedfolds),silent = TRUE)
    if(class(ol1) == "try-error") ol1 <- list(lambda=lambdamax)
    return(ol1)
  }

  if(!parallel) ol1s <- lapply(XXblocks,fastCVX2) else {
    sfExport("Y","X1","kfold","intercept")
    ol1s <- sfLapply(XXblocks,fastCVX2)
  }
  lambdas0 <- unlist(lapply(ol1s,getElement,name="optpen"))
  lambdas <- sapply(lambdas0,function(x) min(x,lambdamax))
  lamol1s = c(list(lambdas),ol1s)
  nbl <- length(XXblocks)
  names(lamol1s) <- c("lambdas",paste("CVres",1:nbl,sep=""))
  return(lamol1s)
}


#####################################################################################################
#
#
#                                     FITTING WITH IWLS ALGORITHM
#
#
#####################################################################################################


#FITTING: efficient IWLS algorithm for logistic and linear (only one iteration)
IWLSridge <- function(XXT,Y,X1=NULL,intercept=TRUE,frac1=NULL,eps=1e-7,maxItr=25,trace=FALSE, model=NULL, E0=NULL){
  #Input:
  #XXT:  sample cross-product (nxn matrix) from penalized variables [computed by SigmaFromBlocks]
  #Y: nx1 response
  #X1: optional nxp1 matrix (p1 < n) with unpenalized covariates
  #intercept: use intercept? If yes, intercept is either estimated or set equal to logit(frac1)
  #frac1: fraction of 1's in the population; if non-null the intercept is estimated as logit(frac1); otherwise intercept is estimated from the data along with all other parameters
  #eps: numerical bound for convergence
  #maxItr: maximum number of iterations used in IWLS
  #trace: should output be traced?

  #Output: List object containing:Ypred=Ypred,convergence=convergence,nIt=it,Hres=Hres, linearized=linearized, unpen=unpen,eta0=eta0))
  #etas: linear predictors, vector (n)
  #Ypred: predicted responses, vector (n)
  #convergence (T/F): T if IWLS has converged under threshold eps, F otherwise
  #nIt: number of iterations needed for convergence
  #Hres: list containing matrices required for predictions on new samples
  #linearized: linearized response vector (n) of last IWLS iteration, required for predictions
  #unpen: does the model include unpenalized parameters, including the intercept?
  #eta0: scalar. Offset = logit(frac1), if frac1 is specified. The same for all samples. Required for predictions.


  #intercept<-FALSE
  #frac1 <- NULL
  #X1 <- matrix(rep(1,n),ncol=1)
  #XXblocks <- list(XXTmir);penalties <- c(mirlambda)  #CIN3=1, normal=1
  #XXblocks <- list(XXTmir,XXTmeth);penalties <- c(mirlambda,methlambda)  #CIN3=1, normal=1
  #maxItr <- 10
  #eps=1e-5

  #in terms of hat matrix
  #note: eta is lin pred without intercept; etatot incl intercept
  if(is.null(model)) model <- ifelse(length(unique(Y)) ==2, "logistic","linear")

  n <- length(Y)
  if(is.null(X1)) unpen <- FALSE else unpen <- TRUE #new 15/4/2021
  if(intercept & is.null(frac1)) X1 <- cbind(matrix(rep(1,n),ncol=1),X1) #if frac1 is not given intercept is estimated alongside all other parameters
  if(intercept & !is.null(frac1)) eta0 <- log((frac1)/(1-frac1)) else eta0 <-  0 #if frac1 is given it is used to estimate the intercept
  if(is.null(E0)) eta <- rep(0,n) else eta <- E0

  if(!is.null(X1)) XXT <- XXT + X1 %*% t(X1)*(1/0.001)  #0.001: penalty for unpenalized variables

  etatot <- eta0 + eta
  it <- 1; nrm <- Inf


  if(model=="logistic") Ypred<-1/(1+exp(-etatot))
  if(model=="linear") {Ypred <- etatot; maxItr<-1}
  while(nrm>eps & it<=maxItr){
    etaold <- eta
    Ypredold <- Ypred
    if(model=="logistic") WV <-Ypred*(1-Ypred) #Var(Y)
    if(model=="linear") WV <- rep(1,n)
    #Winv <- diag(1/(Ypred*(1-Ypred)))
    #XtXDinv <- solve(t(X) %*% W %*% X + 2*diag(lambda,p)) #inverting pxp matrix
    znew <-  matrix(Y - Ypred,ncol=1)
    # inv <- solve(Winv + XXT)
    # Hmat <- XXT -  XXT %*% inv %*% XXT
    linresp <- znew + diag(WV) %*% etaold
    Hresl <- .Hpenlin(WV,XXT,linresp);Hmatl <- Hresl$Hmatl  #NEW 15/4
    #else { Hres <- .Hunpen(WV,XXT,X1); Hmat <- Hres$Hmat; unpen <- TRUE}
    # if(is.null(X1)) Hmat <- .Hpen(WV,XXT) else {
    #    Hmat <- .Hunpen(WV,XXT,X1);Hres <- NULL
    #  }
    #eta <- Hmat %*% znew + XXT %*% inv %*% etaold
    #eta <- Hmat %*% (znew + diag(WV) %*% etaold)
    eta <- Hmatl
    etatot <- eta0 + eta
    if(model=="logistic") Ypred<-.minmaxf(1/(1+exp(-etatot))) #predicted probability
    #if(model=="logistic") Ypred<-1/(1+exp(-etatot)) #predicted probability
    if(model=="linear") Ypred <- etatot
    #if(trace) {print(Ypred); print(WV)}
    if(trace) {print(Ypred)}
    it <- it+1
    nrm <- max(abs(Ypred-Ypredold))
    if(trace) print(nrm)
    if(is.na(nrm)) {it <- maxItr + 1; nrm <- Inf}

    #alt verified for mir data
    # eta2 <- datamir %*% solve(t(datamir) %*% W %*% datamir + diag(rep(penalties[1],ncol(datamir)))) %*% t(datamir) %*% z
  }
  linearized <- znew + diag(WV) %*% etaold
  if(trace) print(it)
  convergence <- nrm<eps
  return(list(etas=eta,Ypred=Ypred,convergence=convergence,nIt=it-1,Hresl=Hresl, linearized=linearized,
              unpen=unpen,intercept=intercept, eta0=eta0, X1=X1)) #KW, MW components of Hres needed for pred, and so is the linearized response; added X1: needed for predictions
}

#IWLS for Cox ridge
IWLSCoxridge <- function(XXT,Y,X1=NULL,intercept=FALSE,eps=1e-7,maxItr=25,trace=FALSE, E0=NULL){
  #Input:
  #XXT:  sample cross-product (nxn matrix) from penalized variables [computed by SigmaFromBlocks]
  #Y: survival response
  #X1: optional nxp1 matrix (p1 < n) with unpenalized covariates
  #intercept: use intercept? If yes, intercept is either estimated or set equal to logit(frac1)
  #frac1: fraction of 1's in the population; if non-null the intercept is estimated as logit(frac1); otherwise intercept is estimated from the data along with all other parameters
  #eps: numerical bound for convergence
  #maxItr: maximum number of iterations used in IWLS
  #trace: should output be traced?

  #Output: List object containing:Ypred=Ypred,convergence=convergence,nIt=it,Hres=Hres, linearized=linearized, unpen=unpen,eta0=eta0))
  #etas: linear predictors, vector (n)
  #Ypred: predicted responses, vector (n)
  #convergence (T/F): T if IWLS has converged under threshold eps, F otherwise
  #nIt: number of iterations needed for convergence
  #Hres: list containing matrices required for predictions on new samples
  #linearized: linearized response vector (n) of last IWLS iteration, required for predictions
  #unpen: does the model include unpenalized parameters, including the intercept?
  #eta0: scalar. Offset = logit(frac1), if frac1 is specified. The same for all samples. Required for predictions.


  #intercept<-FALSE
  #frac1 <- NULL
  #X1 <- matrix(rep(1,n),ncol=1)
  #XXblocks <- list(XXTmir);penalties <- c(mirlambda)  #CIN3=1, normal=1
  #XXblocks <- list(XXTmir,XXTmeth);penalties <- c(mirlambda,methlambda)  #CIN3=1, normal=1
  #maxItr <- 10
  #eps=1e-5

  #in terms of hat matrix
  #note: eta is lin pred without intercept; etatot incl intercept
  #XXT = SigmaFromBlocks(XXomics[2],penalties=lambdas[2]);Y=resp;X1=NULL;intercept=FALSE;frac1=NULL;eps=1e-7;maxItr=100;trace=TRUE
  #XXT = XXT; Y=surv;X1=NULL;intercept=FALSE;frac1=NULL;eps=1e-7;maxItr=25;trace=TRUE; E0 <- NULL
  n <- length(Y)
  eta0 <- 0
  if(is.null(X1)) unpen <- FALSE else unpen <- TRUE #new 15/4/2021
  if(is.null(E0)) eta <- rep(0,n) else eta <- E0
  if(intercept) X1 <- cbind(matrix(rep(1,n),ncol=1),X1) #if frac1 is not given intercept is estimated alongside all other parameters
  if(!is.null(X1)) XXT <- XXT + X1 %*% t(X1)*(1/0.001)  #0.001: penalty for unpenalized variables

  it <- 1; nrm <- Inf
  etatot <- eta
  H0 <- .breslow(Y,etatot)[,2]
  Ypred<- as.numeric(H0 * exp(etatot))
  Yev <- Y[,2]
  while(nrm>eps & it<=maxItr){
    etaold <- eta
    Ypredold <- Ypred
    WV <- Ypred + 10^(-10)
    znew <-  matrix(Yev - Ypred,ncol=1)
    #Hres <- .Hpen(WV,XXT);Hmat <- Hres$Hmat; unpen <- FALSE

    linresp <- znew + diag(WV) %*% etaold
    Hresl <- .Hpenlin(WV,XXT,linresp);Hmatl <- Hresl$Hmatl  #NEW 15/4

    #} else { Hres <- .Hunpen(WV,XXT,X1); Hmat <- Hres$Hmat; unpen <- TRUE}
    #eta <- Hmat %*% (znew + diag(WV) %*% etaold)
    eta <- Hmatl
    etatot <- eta
    H0 <- .breslow(Y,etatot)[,2]
    Ypred<- as.numeric(H0 * exp(etatot))
    if(trace) {print(etatot)}
    #if(trace) {print(WV)}
    it <- it+1
    nrm <- max(abs(Ypred-Ypredold))
    if(trace) print(nrm)
    if(is.na(nrm)) {it <- maxItr + 1; nrm <- Inf}

  }
  linearized <- znew + diag(WV) %*% etaold
  if(trace) print(it)
  convergence <- nrm<eps
  return(list(etas=eta,Ypred=Ypred,convergence=convergence,nIt=it-1,Hresl=Hresl,
              linearized=linearized, unpen=unpen,intercept=intercept,eta0=eta0, X1=X1)) #KW, MW components of Hres needed for pred, and so is the linearized response
}

#####################################################################################################
#
#
#                                     DATA INPUT, PREDICTIONS, COEFFICIENTS (BETAS)
#
#
#####################################################################################################



#CREATES list OF nxn MATRICES REPRESENTING X_b %*% t(X_b) [for fitting], or X^new_b %*% t(X_b) [for prediction]
createXXblocks <- function(datablocks,datablocksnew=NULL,which2pair=NULL){
  #datablocks <- list(X[,1:50],X[,51:75],X[,76:100]); which2pair=c(1,2)
  #datablocks <- list(X[,1:50],X[,51:100])
  nbl <- length(datablocks)
  nfeats <- unlist(lapply(datablocks,nrow))
  if(nbl>1) {
    neq <- sum(sapply(2:nbl,function(b) nfeats[1] == nfeats[b]))
    if(neq != (nbl-1)){
      print("ERROR: Data blocks should have the same number of rows (samples)")
      return(NULL)
    }
  }
  if(is.null(datablocksnew)) res <- lapply(datablocks, function(block) block %*% t(block)) else
    res <- lapply(1:nbl, function(i) datablocksnew[[i]] %*% t(datablocks[[i]]))

  if(!is.null(which2pair)) {
    ncol1 <- ncol(datablocks[[which2pair[1]]])
    ncol2 <- ncol(datablocks[[which2pair[2]]])
    if(ncol1 != ncol2){
      print("ERROR: Paired data blocks should have the same number of columns (features)")
      return(NULL)
    } else {
      permute <- c(((ncol1+1):(2*ncol1)), (1:ncol1))
      dbpaired <- cbind(datablocks[[which2pair[1]]],datablocks[[which2pair[2]]])
      if(is.null(datablocksnew)) res <- c(res, list(dbpaired[,permute] %*% t(dbpaired))) else {
        dbpairednew <- cbind(datablocksnew[[which2pair[1]]],datablocksnew[[which2pair[2]]])
        res <- c(res, list(dbpairednew[,permute] %*% t(dbpaired)))
      }
      cat("Paired block appended to unpaired ones:\n")
      cat(paste("Use pairing = c(",which2pair[1],",",which2pair[2],",",nbl+1,") for further computations",sep=""))
    }
  }
  return(res)
}

#CREATES LIST OF PAIRED DATA BLOCKS
createXblocks <- function(datablocks,which2pair=NULL){
  #datablocks <- list(X[,1:50],X[,51:75],X[,76:100]); which2pair=c(1,2)
  #datablocks <- list(X[,1:50],X[,51:100])
  nbl <- length(datablocks)
  nfeats <- unlist(lapply(datablocks,nrow))
  if(nbl>1) {
    neq <- sum(sapply(2:nbl,function(b) nfeats[1] == nfeats[b]))
    if(neq != (nbl-1)){
      print("ERROR: Data blocks should have the same number of rows (samples)")
      return(NULL)
    }
  }
  res <- lapply(datablocks, function(block) block)
  if(!is.null(which2pair)) {
    ncol1 <- ncol(datablocks[[which2pair[1]]])
    ncol2 <- ncol(datablocks[[which2pair[2]]])
    if(ncol1 != ncol2){
      print("ERROR: Paired data blocks should have the same number of columns (features)")
      return(NULL)
    } else {
      permute <- c(((ncol1+1):(2*ncol1)), (1:ncol1))
      dbpaired <- cbind(datablocks[[which2pair[1]]],datablocks[[which2pair[2]]])
      res <- c(res, list(dbpaired[,permute]))
      cat("Paired block appended to unpaired ones:\n")
      cat(paste("Use pairing = c(",which2pair[1],",",which2pair[2],",",nbl+1,") for further computations",sep=""))
    }
  }
  return(res)
}

#augment data with zero's; to allow pairing of data on DIFFERENT samples
augment <- function(Xdata1,Xdata2){
  #Xdata1 <- dataW3; Xdata2 <- dataF3
  ncol1 <- ncol(Xdata1)
  ncol2 <- ncol(Xdata2)
  if(ncol1 != ncol2){
    print("ERROR: Paired data blocks should have the same number of columns (features)")
    return(NULL)
  } else {
    nsam1 <- nrow(Xdata1)
    nsam2 <- nrow(Xdata2)
    return(list(Xaug1 = rbind(Xdata1,matrix(rep(0,ncol1*nsam2),nrow=nsam2)),Xaug2 = rbind(matrix(rep(0,ncol2*nsam1),nrow=nsam1),Xdata2)))
  }
}

#Computation of overall sample cross-product (covariance) from covariate-block specific cross-products and penalties
SigmaFromBlocks <-function(XXblocks,penalties,pairing=NULL){
  #XXblocks: list of B matrices, dim nxn, each representing block-specific cross-product
  # [each XXTb = t(X_b) %*% X_b for b=1, ..., B data types]. When pairing applies, should be augmented with
  #a cross-block t(X) %*% Q %*% X, where Q = 1 - diag(2), X = [X_{1,.1}, ..., X_{1,.p}, X_{2,.1}, ..., X_{2,.p}].
  #Q is a permutation matrix, so product can be efficiently computed.
  #So paired columns (variables) are juxtaposed in X. See examples.
  #pairing: vector of length 3. First elements denotes the index of the first block matrix (and penalty) involved in the pair,
  #second the index of the 2; third: index of the cross-block (and corresponding pairing penalty).
  #XXblocks =XXblocks; pairing=c(1,2,3);penalties=c(penaltiesstart,4000)

  nblocks <- length(XXblocks)
  if(nblocks != length(penalties)){
    print("Error: Number of penalty parameters should equal number of blocks")
    return(NULL)
  } else {
    if(is.null(pairing)) Sigma<-Reduce('+', lapply(1:nblocks,function(i) XXblocks[[i]] * 1/penalties[i]))
    else {
      nonpaired <- setdiff(1:nblocks,pairing)
      if(length(nonpaired) > 0) Sigma<-Reduce('+', lapply(nonpaired,function(i) XXblocks[[i]] * 1/penalties[i]))  else Sigma <- 0
      invpair <- solve(matrix(c(penalties[pairing[1]],-penalties[pairing[3]],
                                -penalties[pairing[3]],penalties[pairing[2]]),nrow=2))
      invpenpair <- c(invpair[1,1],invpair[2,2],invpair[1,2])
      Sigma <- Sigma + Reduce('+', lapply(pairing,function(i) XXblocks[[i]] * invpenpair[i]))
    }
    return(Sigma)
  }
}


#coefficient estimates given (converged) IWLS weights
betasout <- function(IWLSfit,Xblocks,X1=NULL,penalties,pairing=NULL){
  #IWLSfit <- fit1; Xblocks <- Xomics;penalties <- penpref;pairing=NULL
  n <- length(IWLSfit$etas)
  nfeats <- unlist(lapply(Xblocks, ncol))
  if(IWLSfit$intercept) X1 <- cbind(matrix(rep(1,n),ncol=1),X1)
  if(!is.null(X1)){Xblocks <- c(list(X1),Xblocks); penalties <- c(0.001,penalties)}
  LtX <- .LambdaInvtXFromBlocks(Xblocks,penalties,pairing)
    Mmatl <- IWLSfit$Hresl$Mmatl
    betas <- LtX %*% Mmatl
    if(is.null(X1)){betaunpen <- NULL; betapen <- betas} else {
      nunpen <- ncol(X1); betaunpen <- betas[1:nunpen];betapen <- betas[-(1:nunpen)]
        }
  nbl <- length(nfeats)
  cumindex <- cumsum(c(0,nfeats))
  betaout <- lapply(1:nbl,function(bl) betapen[(cumindex[bl]+1):cumindex[bl+1]])
  betaout <- c(list(betaunpen),betaout)
  return(betaout)
}




#Predictions for new samples
predictIWLS <- function(IWLSfit,X1new=NULL, Sigmanew){
  #IWLSfit: list object; output from IWLSridge
  #X1new: design matrix with unpenalized variables, p1 x nnew; should not contain the intercept
  #Sigmanew: (penalized) sample cross-product between new and training samples (nnew x n); may be computed using SigmaFromBlocks function
  #Output: vector (nnew) of linear predictors

  unpen <- IWLSfit$unpen
  X1 <- IWLSfit$X1
  if(unpen & is.null(X1new)){print("Fit contains unpenalized variables. Please provide argument X1new.")
    return(NULL)
  }
  eta0 <- IWLSfit$eta0
  intercept <- IWLSfit$intercept
  n <- nrow(Sigmanew)
  if(intercept) X1new <- cbind(matrix(rep(1,n),ncol=1),X1new)
  if(!is.null(X1new)){Sigmanew <- Sigmanew +  X1new %*% t(X1)*1/0.001}
  Mmatl <- IWLSfit$Hresl$Mmatl
  Hnewl <- Sigmanew %*% Mmatl
  pred <- Hnewl + eta0
  return(pred)
}

#####################################################################################################
#
#
#                                     CROSS-VALIDATION & SCORING
#
#
#####################################################################################################



#Function for creating (balanced) folds
CVfolds <- function(Y,model=NULL,balance=TRUE,kfold=10,fixedfolds=TRUE,nrepeat=1){ #response is required for balanced CV
  #response: response vector, length n
  #model: "logistic", "cox", etc
  #balance: should the splits balance levels of the response?
  #kfold: scalar, the number of folds
  #fixedfolds: should the folds be fixed? (for reproducibility)
  #nrepeat: number of repeats of the CV
  #Output: list object with kfold elements containing the sample indices of the left-out samples per fold

  if(is.null(model)){
    if(class(Y)=="Surv") model <- "cox" else {
      model <- ifelse(length(unique(Y)) ==2, "logistic","linear")
    }
  }
  response <- Y
  if(model=="linear") balance <- FALSE
  CVfoldsrep <- function(rep){
    nsam <- length(response)
    if (fixedfolds) set.seed(3534+rep-1) #else set.seed(NULL)  #changed 19/4
    if (!balance) {
      rand <- sample(1:nsam)
      grs1 <- floor(nsam/kfold)
      grs2 <- grs1 + 1
      ngr1 <- kfold * grs2 - nsam
      folds <- lapply(1:kfold, function(xg) {
        if (xg <= ngr1)
          els <- rand[(1 + (xg - 1) * grs1):(xg * grs1)]
        else els <- rand[(ngr1 * grs1 + 1 + (xg - ngr1 -
                                               1) * grs2):(ngr1 * grs1 + (xg - ngr1) * grs2)]
        return(sort(els))
      })
    }
    else {
      if (model == "logistic")
        if (class(response) == "factor")
          nev <- which((as.numeric(response) - 1) == 1)
        else nev <- which(response == 1)
        if (model == "cox") nev <- which(response[, 2] == 1)
        nsamev <- length(nev)
        randev <- sample(nev)
        grs1 <- floor(nsamev/kfold)
        grs2 <- grs1 + 1
        ngr1 <- kfold * grs2 - nsamev
        foldsev <- lapply(1:kfold, function(xg) {
          if (xg <= ngr1)
            els <- randev[(1 + (xg - 1) * grs1):(xg * grs1)]
          else els <- randev[(ngr1 * grs1 + 1 + (xg - ngr1 -
                                                   1) * grs2):(ngr1 * grs1 + (xg - ngr1) * grs2)]
          return(els)
        })
        nonev <- setdiff(1:nsam, nev)
        nsamnonev <- length(nonev)
        randnonev <- sample(nonev)
        grs1 <- floor(nsamnonev/kfold)
        grs2 <- grs1 + 1
        ngr1 <- kfold * grs2 - nsamnonev
        foldsnonev <- lapply(1:kfold, function(xg) {
          if (xg <= ngr1)
            els <- randnonev[(1 + (xg - 1) * grs1):(xg *
                                                      grs1)]
          else els <- randnonev[(ngr1 * grs1 + 1 + (xg - ngr1 -
                                                      1) * grs2):(ngr1 * grs1 + (xg - ngr1) * grs2)]
          return(els)
        })
        folds <- lapply(1:kfold, function(i) sort(c(foldsev[[i]],
                                                    foldsnonev[[i]])))
    }
    return(folds)
  }
  return(unlist(lapply(1:nrepeat,CVfoldsrep),recursive=FALSE))
}

#Score function for evaluation of predictions
Scoring <- function(lp,Y,model=NULL,score=ifelse(model=="linear","mse","loglik"),print=TRUE){
  #lp: linear predictor, same size as response
  #response: response vector
  #score: score
  #model: either "logistic", "cox", "linear"
  #output: score (numeric)
  if(is.null(model)){
    if(class(Y)=="Surv") model <- "cox" else {
      model <- ifelse(length(unique(Y)) ==2, "logistic","linear")
    }
  }

  if(model=="linear" && score=="loglik"){
    score <-  "mse"
    print("Switching to MSE score for linear model")
  }
  response <- Y
  minus <- FALSE
  if (model == "linear"){
    if(score=="mse") {thescore <- mean(-(response-lp)^2);minus<-TRUE}
    if(score=="abserror") {thescore <- mean(-abs(response-lp));minus<-TRUE}
    if(score=="cor") thescore <- cor(Y,lp)
    if(score=="spearman") thescore <- cor(Y,lp,method="spearman")
    if(score=="kendall") thescore <- cor(Y,lp,method="kendall")
  }
  if (model == "logistic"){
    probab <- as.numeric(1/(1+exp(-lp)))
    if (class(response) == "factor") Y <- as.numeric(response) - 1 else Y <- response
    if(score=="loglik") thescore <- mean(dbinom(Y,size=1,probab,log=TRUE))
    if(score=="auc") {thescore <- try(as.numeric(pROC::auc(Y,probab,smooth=TRUE)),silent=TRUE);
    suppressWarnings(if(class(thescore)=="try-error") thescore <- as.numeric(pROC::auc(Y,probab,smooth=FALSE)))}
    if(score=="rankboost") {whpos <- which(Y==1); thescore <- mean(1/(1+exp(-2*(outer(probab[whpos],probab[-whpos],"-")))))}
    if(score=="brier") {thescore <- mean(-(Y-probab)^2); minus<-TRUE }
  }
  if (model == "cox"){
    if(score=="loglik") {
      #Y <- surv; lp <- myfit$etas
      hazards <- .breslow(response,lp);
      di <- response[,2]
      ht <- hazards[,1];Ht <- hazards[,2]
      thescore <- sum(-Ht*exp(lp))
      di1 <- which(di==1)
      if(length(di1)>0) thescore <- thescore + sum(di[di1]*(log(ht[di1])+lp[di1]))
      thescore <- thescore/length(lp)
    }
    if(score=="cindex"){
      tm <- max(response[,1])
      iAUCr <- risksetAUC(Stime=response[,1], status= response[,2], marker=lp, method="Cox",plot=FALSE,tmax=tm)
      thescore<- iAUCr$Cindex
      if(is.na(thescore)) thescore <- 0.5
    }
  }
  if(score=="loglik") score <- "mean loglik"
  if(print) {if(minus) {score <- paste("minus",score)}; print(paste(score,": ",round(thescore,3)))}
  return(thescore)
}


#Computes predictive score using CV
CVscore <- function(penalties, XXblocks,Y,X1=NULL, pairing=NULL,folds, intercept=ifelse(class(Y)=="Surv", FALSE, TRUE),frac1=NULL,score="loglik", model = NULL,
                    eps=1e-7,maxItr=100,trace=FALSE,printCV=TRUE,save=FALSE,parallel=FALSE){

  #XXblocks: list of B matrices, dim nxn, each representing block-specific cross-product [each XXTb = t(X_b) %*% X_b for b=1, ..., B data types].
  #penalties: vector with B strictly positive penalties
  #response: nx1 response
  #X1: optional nxp1 matrix (p1 < n) with unpenalized covariates
  #pairing: see SigmaFromBlocks
  #folds: list of k folds, containing left-out samples per fold (may be computed by CVfolds function).
  #intercept: use intercept? If yes, intercept is either estimated or set equal to logit(frac1)
  #frac1: fraction of 1's in the population; if non-null the intercept is estimated as logit(frac1); otherwise intercept is estimated from the data along with all other parameters
  #score: score
  #model: either "logistic", "cox"
  #eps: numerical bound for IWLS convergence
  #maxItr: maximum number of iterations used in IWLS
  #trace: should IWLSoutput be traced?
  #printCV: should the CV-score be printed on screen?
  #save: boolean. If TRUE appends the penalties and resulting CVscore to global variable 'allscores'

  #Output: numeric, CV-ed prediction score for given penalties

  # penalties<-c(mirlambda,methlambda); XXblocks <- list(XXTmir,XXTmeth); response<-Y; X1=NULL;folds <- leftout;
  # intercept=ifelse(class(Y)=="Surv", FALSE, TRUE); frac1=NULL;score="loglik"; model="logistic"; eps=1e-7;maxItr=100;trace=FALSE
  #XXblocks =list(XXTcn);penalties <- 1824;response=surv;X1=NULL;intercept=FALSE;frac1=NULL;eps=1e-7;maxItr=100;trace=FALSE
  #parallel <- F; save<-F; printCV <-F; pairing <- NULL
  #XXblocks =XXblocks[2];penalties <- 200;response=Y2;X1=NULL;intercept=FALSE;frac1=NULL;eps=1e-7;maxItr=100;trace=FALSE
  #parallel <- F; save<-F; printCV <-F; pairing <- NULL; model="linear";score="mse"
  # penalties = lambdas[2]; XXblocks=XXomics[2];Y=resp;X1=NULL; pairing=NULL;
  # penalties=lambdas2; XXblocks=XXomics;Y=resp;X1<- NULL; folds=leftout; intercept=FALSE;frac1=NULL;score="loglik";
  #model="cox";eps= 1e-7;maxItr=100;trace=FALSE;
  # printCV=TRUE; save=FALSE;parallel=FALSE;pairing <- NULL
  #

  if(is.null(model)){
    if(class(Y)=="Surv") model <- "cox" else {
      model <- ifelse(length(unique(Y)) ==2, "logistic","linear")
    }
  }

  response <- Y
  if(parallel){if(sfCpus()==1) {
    message("ERROR: please run the setupParallel() function first")
    #message("Use sourcelibraries=c(\"penalized\")")
    return(NULL)
  }}
  nf <- length(folds)
  n <- length(response)
  nblocks <- length(XXblocks)
  SigmaLambda <- SigmaFromBlocks(XXblocks,penalties,pairing)
  if(parallel) system.time(sfExport("folds","model", "SigmaLambda","response","X1","intercept","frac1","eps","maxItr","trace","n"))
  lpsoutfun <- function(i){
    #i<-1
    outs <- folds[[i]]
    ins <- (1:n)[-outs]
    respin <- response[ins]
    respout <- response[outs]
    SigmaIn <- SigmaLambda[ins,ins] #Sigma matrix allows subsetting; SigmaIn is required for fitting
    SigmaOutIn <- SigmaLambda[outs,ins] #SigmaOutIn is required for prediction
    X1in <- X1[ins,]
    X1out <- X1[outs,]
    if(model=="logistic" | model=="linear") fitin <- try(IWLSridge(XXT=SigmaIn,Y=respin,X1=X1in,intercept=intercept,frac1=frac1,eps=eps,maxItr=maxItr,trace=trace, model=model))
    if(model=="cox") fitin <- try(IWLSCoxridge(XXT=SigmaIn,Y=respin,X1=X1in,intercept=intercept,eps=eps,maxItr=maxItr,trace=trace))
    if(class(fitin) =="try-error") lpsout <- rep(0,length(respout)) else lpsout <- predictIWLS(fitin, X1out,SigmaOutIn)
    return(lpsout)
  }
  if(parallel) lpsoutlist <- sfLapply(1:nf,lpsoutfun) else lpsoutlist <- lapply(1:nf,lpsoutfun)
  all_lpsout <- unlist(lpsoutlist)
  if(model!="cox") all_respout <- unlist(lapply(1:nf,function(i){out <- folds[[i]];respout <- response[out];return(respout)})) else {  #cox model
    all_evout <- unlist(lapply(1:nf,function(i){out <- folds[[i]];respout <- response[out,2];return(respout)}))
    all_timeout <- unlist(lapply(1:nf,function(i){out <- folds[[i]];respout <- response[out,1];return(respout)}))
    all_respout <- Surv(all_timeout,all_evout)
  }
  cvscore <- Scoring(all_lpsout,all_respout,score=score,model=model,print=FALSE)
  #nrep <- sum(unlist(lapply(folds,length)))/n #nr of CV repeats
  #if(length(cvscores)!=1) cvscore <- sum(cvscores)/nrep  else cvscore <- cvscores #CV is averaged across CV-repeats [but not for auc type scores]
  if(printCV) {print("penalties:"); print(penalties); print(paste("CV-score:", round(cvscore,3)))}
  return(cvscore)
  #return(lpsoutlist)
}

#####################################################################################################
#
#
#                                     OPTIMIZATION
#
#
#####################################################################################################

#OPTIMIZATION OF CV CRITERION; ONE OPTIMIZER
optLambdas <- function(penaltiesinit=NULL, XXblocks,Y,X1=NULL, pairing=NULL, folds, intercept=ifelse(class(Y)=="Surv", FALSE, TRUE),frac1=NULL,score="loglik", model=NULL,
                       epsIWLS=1e-3,maxItrIWLS=25, traceCV=TRUE, reltol=1e-4, optmethod=ifelse(length(penaltiesinit)==1,"Brent", "Nelder-Mead"),
                       maxItropt=500, save=FALSE, parallel=FALSE, fixedpen=NULL, fixedseed=TRUE){
  #Arguments additional to those of CVscore:
  #penaltiesinit: initial penalties; if NULL (not recommended) then set to 100. Initial penalties can be found by running function "fastCV"
  #multinit: initial penalty multipliers; these are multiplied by penaltiesinit to opbtain penalites; optimizer optimizes the multipliers (log scale)
  #traceCV: boolean: should the CV scores be traced during optimization?
  #reltol: relative tolerance of the optimizer (see optim)
  #optmethod: optimization method
  #Note: optimization of the penalty multipliers is done on log-scale. This allows unconstrained optimization
  # fixedseed: some methods, in particular SANN use random searches; then fixing the seed may be desirable for reproducibility

  #penaltiesinit=1000; multinit = c(1);XXblocks=list(XXTcn);response=surv;X1=NULL;folds=leftout; intercept=ifelse(class(Y)=="Surv", FALSE, TRUE);frac1=NULL;score="loglik";
  # model="cox"; eps=1e-7;maxItr=100;trace=FALSE;reltol=10^(-4); optmethod="Brent";save=T
  #penaltiesinit=penaltiesstart; XXblocks=XXblocks; pairing = c(1,2,3)
  # penaltiesinit=lambdas; XXblocks=XXbl;Y=resp;X1=NULL;folds=leftout;intercept=ifelse(class(Y)=="Surv", FALSE, TRUE);frac1=NULL;score="loglik"; model="logistic";maxItrIWLS = 25;
  # parallel=FALSE;epsIWLS=1e-3;maxItropt=500;trace=FALSE;reltol=10^(-4);optmethod="Nelder-Mead";save=T;traceCV=FALSE;pairing <- NULL
  # penaltiesinit=lambdas; XXblocks=XXomics12;Y=resp;X1=NULL;folds=leftout;intercept=ifelse(class(Y)=="Surv", FALSE, TRUE);frac1=NULL;score="loglik"; model="logistic";maxItrIWLS = 25;
  # parallel=TRUE;epsIWLS=1e-3;maxItropt=500;trace=FALSE;reltol=10^(-4);optmethod="Nelder-Mead";save=T;traceCV=FALSE;pairing = c(1,2,3);
  # fixedseed=TRUE;fixedpen=NULL

  if(is.null(model)){
    if(class(Y)=="Surv") model <- "cox" else {
      model <- ifelse(length(unique(Y)) ==2, "logistic","linear")
    }
  }

  if(fixedseed) set.seed(12345)

  response <- Y

  if(parallel){if(sfCpus()==1) {
    message("ERROR: please run the setupParallel() function first")
    #message("Use sourcelibraries=c(\"penalized\")")
    return(NULL)
  }}
  npen <- length(XXblocks)
  if(!is.null(fixedpen)) {
    if(length(fixedpen) > (npen-1) | max(fixedpen) > npen){
      message("ERROR: incorrect input for argument \'fixedpen\'")
      return(NULL)
    }
  }


  if(is.null(penaltiesinit)) {penaltiesinit <- rep(100,npen); if(!is.null(pairing)) {penpair <- pairing[3]; penaltiesinit[penpair] <- 25}
  } else {if(!is.null(pairing)) {
    if(length(penaltiesinit) < npen)
    {
      penaltiesinit <- c(penaltiesinit,paired=NA)
      penpair1 <- pairing[1]; penpair2 <- pairing[2]; penpair3 <- pairing[3]
      penaltiesinit[penpair3] <- sqrt(penaltiesinit[penpair1])*sqrt(penaltiesinit[penpair2])*0.25
      penaltiesinit[penpair1] <- penaltiesinit[penpair1]*(1+.25)
      penaltiesinit[penpair2] <- penaltiesinit[penpair2]*(1+.25)
    }
  }
  }  #use smaller paired penalty by default

  print("Initial penalties:")
  print(penaltiesinit)
  multinit = rep(1,length(penaltiesinit))

  if(!is.null(fixedpen)){
    multinit <- multinit[-fixedpen]
  }

  if(length(multinit) == 1 & optmethod != "SANN") optmethod <- "Brent"
  print(paste("Using",optmethod,"for optimization"))

  CVS <- function(logmults){
    allpenalties <- rep(NA,length(penaltiesinit))
    if(!is.null(fixedpen)) {
      allpenalties[fixedpen] <- penaltiesinit[fixedpen]
      allpenalties[-fixedpen] <- exp(logmults)*penaltiesinit[-fixedpen]
    } else allpenalties <- exp(logmults)*penaltiesinit

    res <- -CVscore(allpenalties, XXblocks=XXblocks,Y=response,X1=X1, pairing=pairing,
                    folds=folds, intercept=intercept,frac1=frac1,score=score, model=model, eps=epsIWLS,
                    maxItr=maxItrIWLS,trace=FALSE, printCV=traceCV, save=save,parallel=parallel)
    #assign("allscores", rbind(allscores,c(res,allpenalties)), envir = baseenv())
    return(res)
  }

  #assign("allscores", c(), envir = baseenv())
  #minimal penalty
  lb <- log10(0.1/max(penaltiesinit))

  if(optmethod != "Brent") optres <- optim(par = log(multinit), CVS, control=list(reltol=reltol,maxit=maxItropt), method=optmethod) else
    optres <- try(optim(par = log(multinit), CVS, control=list(reltol=reltol,maxit=maxItropt), method=optmethod,lower=lb,upper=10))
  if(class(optres)=="try-error") optres <- optim(par = log(multinit), CVS, control=list(reltol=reltol,maxit=maxItropt), method="Nelder-Mead")
  #if(!is.null(ncol(allscores))) colnames(allscores) <- c("Score", paste("pen",1:npen,sep=""))

  optpen <- rep(NA,length(penaltiesinit))
  if(!is.null(fixedpen)) {
    optpen[fixedpen] <- penaltiesinit[fixedpen]
    optpen[-fixedpen] <- as.numeric(exp(optres$par))*penaltiesinit[-fixedpen]
  } else optpen <- as.numeric(exp(optres$par))*penaltiesinit

  #optres2 <- c(optres,optpen=list(optpen),allsc = list(allscores))
  optres2 <- c(optres,optpen=list(optpen))
  return(optres2)
}

#WRAPPER FOR OPTIMIZING IN TWO STEP: 1 GLOBAL SEARCH, 2 LOCAL SEARCH
optLambdasWrap <- function(penaltiesinit=NULL, XXblocks,Y,X1=NULL, pairing=NULL, folds, intercept=ifelse(class(Y)=="Surv", FALSE, TRUE),frac1=NULL,score="loglik",
                           model=NULL,
                           epsIWLS=1e-3,maxItrIWLS=25, traceCV=TRUE, reltol=1e-4, optmethod1= "SANN", optmethod2 =ifelse(length(penaltiesinit)==1,"Brent", "Nelder-Mead"),
                           maxItropt1=10,maxItropt2=25, save=FALSE, parallel=FALSE, pref=NULL,fixedpen=NULL){

  # penaltiesinit=NULL; X1=NULL; pairing=NULL; folds; intercept=ifelse(class(Y)=="Surv", FALSE, TRUE);frac1=NULL;score="loglik"; model="logistic";
  # epsIWLS=1e-3;maxItrIWLS=25; traceCV=TRUE; reltol=1e-4; optmethod1= "SANN"; optmethod2 =ifelse(length(penaltiesinit)==1,"Brent", "Nelder-Mead");
  # maxItropt1=10;maxItropt2=25; save=FALSE; parallel=FALSE;  pref=NULL;
  # penaltiesinit=lambdas; XXblocks=XXall;Y=resp;folds=leftout;score="loglik"; model="cox"; intercept=FALSE;
  # parallel=TRUE;

  if(is.null(model)){
    if(class(Y)=="Surv") model <- "cox" else {
      model <- ifelse(length(unique(Y)) ==2, "logistic","linear")
    }
  }

  nbl <- length(XXblocks)

  if(is.null(pref)){prefer <- 1:nbl; peni <- penaltiesinit} else {prefer <- pref;
  if(is.null(penaltiesinit)) peni<-NULL else peni <- penaltiesinit[prefer]}

  jl0 <- try(optLambdas(penaltiesinit=peni, XXblocks=XXblocks[prefer],Y=Y,X1=X1, pairing=pairing, folds=folds, intercept=intercept,frac1=frac1,score=score, model=model,
                        epsIWLS=epsIWLS,maxItrIWLS=maxItrIWLS, traceCV=traceCV, reltol=reltol,optmethod=optmethod1,maxItropt=maxItropt1,
                        save=save, parallel=parallel, fixedpen=fixedpen))

  if(class(jl0) != "try-error") lambdas0 <- jl0$optpen else lambdas0 <- peni

  jl1 <- optLambdas(penaltiesinit=lambdas0, XXblocks=XXblocks[prefer],Y=Y,X1=X1, pairing=pairing, folds=folds, intercept=intercept,frac1=frac1,score=score, model=model,
                    epsIWLS=epsIWLS,maxItrIWLS=maxItrIWLS, traceCV=traceCV, reltol=reltol,optmethod=optmethod2,maxItropt=maxItropt2,
                    save=save, parallel=parallel, fixedpen=fixedpen )
  lambdas1 <- jl1$optpen

  res <- c(res0=list(jl0),res1=list(jl1))
  lambdas <- c(list(lambdas0),list(lambdas1))

  if(is.null(pref)) {
    return(list(res=res,lambdas=lambdas,optpen = lambdas1))
  } else {
    pref <- c(fixedpen,pref)
    lambdasall <- rep(NA,nbl)
    lambdasall[pref] <- lambdas1
    if(!is.null(penaltiesinit)) lambdasall[-pref] <- penaltiesinit[-pref] else lambdasall[-pref] <- 100
    jl2 <- try(optLambdas(penaltiesinit=lambdasall, XXblocks=XXblocks,Y=Y,X1=X1, pairing=pairing, folds=folds, intercept=intercept,frac1=frac1,score=score, model=model,
                          epsIWLS=epsIWLS,maxItrIWLS=maxItrIWLS, traceCV=traceCV, reltol=reltol,optmethod=optmethod1,maxItropt=maxItropt1,
                          save=save, parallel=parallel,fixedpen= pref ))

    if(class(jl2) != "try-error") lambdas2 <- jl2$optpen else lambdas2 <- lambdasall


    jl3 <- optLambdas(penaltiesinit=lambdas2, XXblocks=XXblocks,Y=Y,X1=X1, pairing=pairing, folds=folds, intercept=intercept,frac1=frac1,score=score, model=model,
                      epsIWLS=epsIWLS,maxItrIWLS=maxItrIWLS, traceCV=traceCV, reltol=reltol,optmethod=optmethod2,maxItropt=maxItropt2,
                      save=save, parallel=parallel, fixedpen=pref)
    lambdas3<- jl3$optpen

    res <- c(res,res2=list(jl2),res3=list(jl3))
    lambdas <- c(lambdas,list(lambdas2),list(lambdas3))
    return(list(res=res,lambdas=lambdas,optpen = lambdas3))
  }
}


#####################################################################################################
#
#
#                                    DOUBLE CV FOR PERFORMANCE EVALUATION
#
#
#####################################################################################################


doubleCV <- function(penaltiesinit,XXblocks,Y,X1=NULL,pairing=NULL,outfold=5, infold=10, nrepeatout=1, nrepeatin=1, balance=TRUE, fixedfolds=TRUE, intercept=ifelse(class(Y)=="Surv", FALSE, TRUE),frac1=NULL,score="loglik", model=NULL,
                     eps=1e-7,maxItr=10,trace=FALSE,printCV=TRUE,reltol=1e-4, optmethod1= "SANN", optmethod2 =ifelse(length(penaltiesinit)==1,"Brent", "Nelder-Mead"),
                     maxItropt1=10,maxItropt2=25, save=FALSE, parallel=FALSE, pref=NULL,fixedpen=NULL){
  #1. produce outfolds
  #2. subset XXblocks and response and X1
  #3. opt lambda per out-fold
  #4. Scoring

  # XXblocks <- XXblocksmm; response<-Y; X1=NULL;outfold=5;infold=10;nrepeatout=2; nrepeatin=1;balance=TRUE;fixedfolds=TRUE;
  # intercept=ifelse(class(Y)=="Surv", FALSE, TRUE); frac1=NULL;score="loglik"; model="logistic"; eps=1e-7;maxItr=100;trace=FALSE;printCV=TRUE;save=FALSE;parallel=F;
  #reltol=1e-4; optmethod=ifelse(length(XXblocks)==1,"Brent", "Nelder-Mead")
  if(is.null(model)){
    if(class(Y)=="Surv") model <- "cox" else {
      model <- ifelse(length(unique(Y)) ==2, "logistic","linear")
    }
  }

  response <- Y
  leftoutout <- CVfolds(Y=response,kfold=outfold,nrepeat=nrepeatout,fixedfolds=fixedfolds,balance=balance,model=model)
  n <- length(response)
  nbl <- length(XXblocks)
  npen <- length(penaltiesinit)
  if((npen != nbl & is.null(pairing)) | (npen != (nbl-1) & !is.null(pairing))){
    message("ERROR: please make sure penaltiesinit is the same length as XXblocks. (In case of pairing: paired penalty should not be specified; will be initialized automatically)")
    return(NULL)
  }

  if(parallel) sfExport("leftoutout","XXblocks","response","X1","intercept","frac1","eps","maxItr","trace","n","nrepeatin",
                        "fixedfolds","balance", "model","score","reltol","optmethod1", "optmethod2","maxItropt1","maxItropt2","pref","fixedpen")

  lpsoutfun <- function(i){
    if(parallel) sfCat(paste("Outer fold:",i),sep="\n") else print(paste("Outer fold:",i))
    #i <- 1
    outs <- leftoutout[[i]]
    ins <- (1:n)[-outs]
    respin <- response[ins]
    respout <- response[outs]
    XXblocksin <- lapply(XXblocks, function(block,whichin) block[whichin,whichin],whichin=ins)
    X1in <- X1[ins,]
    X1out <- X1[outs,]
    foldsin <- CVfolds(Y=respin,kfold=infold,nrepeat=nrepeatin,fixedfolds=fixedfolds,balance=balance,model=model)
    # if(is.null(penaltiesinit)){
    # if(is.null(pairing)) penaltiesi <- sapply(1:nbl,function(bl) fastCV2(XXblocksin[[bl]],Y=respin, folds=foldsin, penaltyinit=NULL, X1=X1in,
    #                                                 intercept=intercept,frac1=frac1,score=score, model=model, epsIWLS=eps,
    #                                                 maxItrIWLS=maxItr, traceCV=trace, reltol=reltol, parallel=FALSE, NMinit=FALSE)[[1]])
    #   else { #initialization based on parametrization in the manuscript
    #        penaltiesi <- rep(NA,nbl)
    #        pair3 <- pairing[3]
    #        allbl <- setdiff(1:nbl,pair3)
    #        penaltiesi[allbl] <- sapply(allbl,function(bl) fastCV2(XXblocksin[[bl]],Y=respin, folds=foldsin, penaltyinit=NULL, X1=X1in,
    #                                                        intercept=intercept,frac1=frac1,score=score, model=model, epsIWLS=eps,
    #                                                        maxItrIWLS=maxItr, traceCV=trace, reltol=reltol, parallel=FALSE, NMinit=FALSE)[[1]])
    #        penpair1i <- penaltiesi[pairing[1]];penpair2i <- penaltiesi[pairing[2]]
    #        penpair3 <- sqrt(penpair1i)*sqrt(penpair2i)/4
    #        penaltiesi[pair3] <- penpair3
    #        penaltiesi[pairing[1]] <- penpair1i*(1+1/4);penaltiesi[pairing[2]] <- penpair2i*(1+1/4)
    #   }
    # } else
    #penaltiesinit <- c(4,5,6); pairing=c(1,2)
    if(!is.null(pairing)) {
      penaltiesi <- penaltiesinit
      penaltiesi <- c(penaltiesi, sqrt(penaltiesi[pairing[1]])*sqrt(penaltiesi[pairing[2]])/4)
      penaltiesi[pairing] <-  penaltiesi[pairing]*(1+1/4)
    } else penaltiesi <- penaltiesinit
    if(nbl>1 | !is.null(penaltiesi)) optlam <- optLambdasWrap(penaltiesinit=penaltiesi, XXblocks=XXblocksin,Y=respin,X1=X1in, pairing=pairing, folds=foldsin, intercept=intercept,frac1=frac1,
                                                              score=score, model=model, epsIWLS=eps,maxItrIWLS=maxItr, traceCV=trace, reltol=reltol,
                                                              optmethod1=optmethod1,optmethod2=optmethod2, maxItropt1=maxItropt1, maxItropt2=maxItropt2, save=FALSE, parallel=FALSE, pref=pref,fixedpen=fixedpen)$optpen else optlam = penaltiesi
    print("Opt penalties:");print(optlam)

    SigmaLambda <- SigmaFromBlocks(XXblocks,optlam, pairing=pairing)
    SigmaIn <- SigmaLambda[ins,ins] #Sigma matrix allows subsetting; SigmaIn is required for fitting
    SigmaOutIn <- SigmaLambda[outs,ins] #SigmaOutIn is required for prediction
    if(model!="cox") fitin <- try(IWLSridge(SigmaIn,respin,X1=X1in,intercept=intercept,frac1=frac1,eps=eps,maxItr=maxItr,trace=trace))
    if(model=="cox") fitin <- try(IWLSCoxridge(SigmaIn,respin,X1=X1in,intercept=intercept,eps=eps,maxItr=maxItr,trace=trace))
    if(class(fitin) =="try-error") lpsout <- rep(0,length(respout)) else lpsout <- predictIWLS(fitin, X1out,SigmaOutIn)
    return(lpsout)
  }
  nf <- length(leftoutout)
  if(parallel) lpsoutlist <- sfLapply(1:nf,lpsoutfun) else lpsoutlist <- lapply(1:nf,lpsoutfun)
  all_lpsout <- unlist(lpsoutlist)
  #all_respout <- unlist(lapply(1:nf,function(i){out <- leftoutout[[i]];respout <- response[out];return(respout)}))
  if(model!="cox") all_respout <- unlist(lapply(1:nf,function(i){out <- leftoutout[[i]];respout <- response[out];return(respout)})) else {  #cox model
    all_evout <- unlist(lapply(1:nf,function(i){out <- leftoutout[[i]];respout <- response[out,2];return(respout)}))
    all_timeout <- unlist(lapply(1:nf,function(i){out <- leftoutout[[i]];respout <- response[out,1];return(respout)}))
    all_respout <- Surv(all_timeout,all_evout)
  }
  all_samout <- unlist(lapply(1:nf,function(i){out <- leftoutout[[i]];return(out)}))
  whfold <- unlist(lapply(1:nf,function(i){out <- leftoutout[[i]];return(rep(i,length(out)))}))
  res0 <- data.frame(sampleindex=all_samout,true=all_respout, linpred = all_lpsout, whfold=whfold)
  od <- order(all_samout)
  res <- res0[od,]
  # rn <- rownames(XXblocks[[1]])
  # if(!is.null(rn)) rownames(res) <- sapply(all_samout,function(sam) rn[sam])
  return(res)
}






#####################################################################################################
#
#
#                                     AUXILIARY FUNCTIONS
#
#
#####################################################################################################


#Computation of weighted Hat matrix, penalized variables only; needed in IWLS functions
.Hpen <- function(WV,XXT){ #WV:weigths as vector (n); XXT: (penalized) sample cross-product (nxn)
  n <-length(WV)
  Winv <- diag(1/WV)
  inv <- solve(Winv + XXT)
  Mmat <- diag(n) - inv %*% XXT
  Hmat <- XXT %*% Mmat
  return(list(Hmat=Hmat,Mmat=Mmat))
}

.Hpenlin <- function(WV,XXT,lr){ #WV:weigths as vector (n); XXT: (penalized) sample cross-product (nxn); lr: linearized response
  n <-length(WV)
  Winv <- diag(1/WV)
  XXTlr <- XXT %*% lr
  invl <- solve(Winv + XXT,XXTlr)
  Mmatl <- lr - invl
  Hmatl <- XXT %*% Mmatl
  return(list(Hmatl=Hmatl,Mmatl=Mmatl))
}

# n <- 500
# p <- 5000
# X <- matrix(rnorm(n*p),nrow=n)
# Y <- matrix(rnorm(n),nrow=n)
# WV <- rep(1,n)
# XXT <- X %*% t(X)
# system.time(res1 <- .Hpen(WV,XXT)$Hmat %*% Y)
# system.time(res2 <- .Hpenlin(WV,XXT,Y)$Hmatl)

#Computation of weighted Hat matrix, accounting for unpenalized variables; needed in IWLS functions
.Hunpen <- function(WV,XXT,X1){
  #WV:weigths as vector (n)
  #XXT: sample cross-product (nxn) penalized variables
  #X1 (p1xn) design matrix unpenalized variables
  n <- length(WV)
  WsqrtV <- sqrt(WV)
  WinvsqrtV <- 1/WsqrtV
  X1W <- WsqrtV * X1
  X1aux <- solve(t(X1W) %*% X1W) %*% t(X1W)
  X1Wproj <- X1W %*% X1aux
  P1W <- diag(n) - X1Wproj
  GammaW <- t(t(WsqrtV * XXT) * WsqrtV) #faster
  P1GammaW <- P1W %*% GammaW
  #Mir: vanaf dit stukje
  invforH2<-try(solve(diag(n) + P1GammaW))
  if(class(invforH2)[1]=="try-error"){
    svdXXT <- svd(diag(n) + P1GammaW)
    svdd <- svdXXT$d
    #reci <- 1/svdd[1:n]
    reci <- c(1/svdd[1:n-1],0)
    invforH2 <- svdXXT$v %*% (reci * t(svdXXT$u))
  }
  #Mir: tot hier
  #H2 <- Winvsqrt %*% GammaW %*% (diag(n) -invforH2 %*% P1GammaW) %*% P1W %*% Winvsqrt
  MW <- WsqrtV * t(t((diag(n) - invforH2 %*% P1GammaW) %*% P1W) * WinvsqrtV)
  H20 <- GammaW %*% (WinvsqrtV * MW)
  H2 <- t(t(WinvsqrtV * H20)) #faster
  #Hboth <- Winvsqrt %*% X1Wproj %*% (diag(n) - Wsqrt %*% H2 %*% Wsqrt) %*% Winvsqrt + H2
  KW <- t(t(X1aux %*% (diag(n) - t(t(WsqrtV * H2) * WsqrtV))) * WinvsqrtV)
  Hmat <-(WinvsqrtV * X1W) %*% KW  + H2
  return(list(Hmat=Hmat,MW=MW,KW=KW))
}

#auxiliary function for betasout
.LambdaInvtXFromBlocks<-function(Xblocks,penalties,pairing=NULL){
  #Xblocks: list of B matrices, dim nxp, each representing block-specific design
  #  When pairing applies, should be augmented with
  #a block Q %*% X, where Q = 1 - diag(2), X = [X_{1,.1}, ..., X_{1,.p}, X_{2,.1}, ..., X_{2,.p}].
  #Q is a permutation matrix, so product can be efficiently computed.
  #pairing: vector of length 3. First elements denotes the index of the first block matrix (and penalty) involved in the pair,
  #second the index of the 2; third: index of the cross-block (and corresponding pairing penalty).
  #Xblocks =Xblocks; pairing=c(1,2,3);penalties=c(penaltiesstart,4000)
  #Xblocks <- list(X[,1:p1],X[,-(1:p1)]);penalties <- c(1,5);pairing=NULL
  #Xblocks = Xblocks;penalties = c(3,2,1); pairing=c(1,2,3)

  nblocks <- length(Xblocks)
  if(nblocks != length(penalties)){
    print("Error: Number of penalty parameters should equal number of blocks")
    return(NULL)
  } else {
    if(is.null(pairing)) LtX <-Reduce('rbind', lapply(1:nblocks,function(i) t(Xblocks[[i]]) * 1/penalties[i]))
    else {
      nonpaired <- setdiff(1:nblocks,pairing)
      if(length(nonpaired) > 0) LtX <- Reduce('rbind', lapply(nonpaired,function(i) t(Xblocks[[i]]) * 1/penalties[i]))  else LtX <- c()
      invpair <- solve(matrix(c(penalties[pairing[1]],-penalties[pairing[3]],
                                -penalties[pairing[3]],penalties[pairing[2]]),nrow=2))
      invpenpair <- c(invpair[1,1],invpair[2,2],invpair[1,2])
      LtX <- rbind(LtX,t(Xblocks[[3]]) * invpenpair[3] + Reduce('rbind', lapply(pairing[1:2],function(i) t(Xblocks[[i]]) * invpenpair[i]))) #
    }
    return(LtX)
  }
}

#auxiliary function; thresholds predictions to avoid numerical problems; in IWLSRidge
.minmaxf <- function(vec,thr=10^(-3)){
  sapply(vec,function(x) min(max(x,thr),1-thr))
}


#NEEDED FOR IWLSCoxRidge
.breslow <- function(Y,lp){ #checked to coincide with basehaz from penalized (which coincides with survfit)
  #Y survival response; lp: linear predictor (length n)
  #Returns hazard and cumulative hazard at all timepoints
  ord <- order(Y)
  invord <- order(ord) #to put back in original order
  Ys <- Y[ord]
  di <- Ys[,2] #event
  ti <- Ys[,1] #time
  lp <- lp[ord]
  htsort <- di/rev(cumsum(exp(rev(lp))))
  ht <- htsort[invord]
  Ht <- cumsum(htsort)[invord]
  return(data.frame(ht=ht,Ht=Ht))
}


#####################################################################################################
#
#
#                         FUNCTIONS FOR MAXIMIZING MARGINAL LIKELIHOOD
#
#
#####################################################################################################


mgcv_lambda <- function(penalties, XXblocks,Y, model=NULL, printscore=TRUE, pairing=NULL, sigmasq = 1,
                        opt.sigma=ifelse(model=="linear",TRUE, FALSE)){

  #penalties <- allpenalties; XXblocks=XXblocks; Y=response; model=model,printscore=TRUE; pairing =NULL;sigmasq=sigmasq;
  if(is.null(model)){
    if(class(Y)=="Surv") model <- "cox" else {
      model <- ifelse(length(unique(Y)) ==2, "logistic","linear")
    }
  }
  #penalties=lambdas; XXblocks = XXbl;Y=resp; model="logistic";printscore=TRUE; intercept= TRUE
  if(opt.sigma){ #Mir
    sigmasq <- penalties[1]#Mir
    penalties <- penalties[-1]#Mir
  }#Mir
  if(model=="linear") fam <- "gaussian"
  if(model=="logistic") fam <- "binomial"
  if(model=="cox") fam <- "cox.ph"

  XXT <- SigmaFromBlocks(XXblocks,penalties=penalties, pairing=pairing)
  n <- nrow(XXT)
  XXTi <- try(solve(XXT), silent=TRUE)
  if(class(XXTi)[1]=="try-error"){ #rank of matrix might be n-1 due to scaling... take Moore-Penrose inverse
    pmt <- proc.time()
    svdXXT <- svd(XXT)
    svdd <- svdXXT$d
    #reci <- 1/svdd[1:n]
    reci <- c(1/svdd[1:n-1],0)
    XXTi <- svdXXT$v %*% (reci * t(svdXXT$u))
    proc.time()-pmt
  }

  Xgam <- diag(n)
  PP = list(Xgam=list(XXTi,sp=1))

  #MML
  pred <- try(gam(Y ~ 0 + Xgam,family=fam, paraPen=PP,method="ML",scale=sigmasq)) #no intercept

  score <- pred$gcv.ubre
  if(printscore) {
    print(penalties)
    print(score)
  }
  return(score)
}

optLambdas_mgcv <- function(penaltiesinit=NULL, XXblocks,Y, pairing=NULL,model=NULL, reltol=1e-4, optmethod =ifelse(length(penaltiesinit)==1,"Brent", "Nelder-Mead"),
                            maxItropt=500,tracescore=TRUE,fixedpen=NULL, fixedseed = TRUE, sigmasq = 1,
                            opt.sigma=ifelse(model=="linear",TRUE, FALSE)){
  #penaltiesinit=lambdas; Xblocks=Xomics;Y=resp; model="cox";reltol=1e-4; optmethod= "Nelder-Mead";maxItropt=20
  # penaltiesinit=lambdas; XXblocks=XXblocks;Y=respnum; pairing<- NULL;model="linear";reltol=1e-4; optmethod= "SANN";maxItropt=2;
  # tracescore=TRUE;fixedseed=TRUE; pref=NULL;fixedpen=NULL; sigmasq = 1;opt.sigma=ifelse(model=="linear",TRUE, FALSE)
   # penaltiesinit=c(100); XXblocks=XXblocks[1];
   # Y=respnum;reltol=0.1;pairing=NULL;model=NULL; optmethod=ifelse(length(penaltiesinit)==1,"Brent", "Nelder-Mead");
   # maxItropt=500;tracescore=TRUE;fixedpen=NULL; fixedseed = TRUE; sigmasq = 1; opt.sigma=FALSE


  if(is.null(model)){
    if(class(Y)=="Surv") model <- "cox" else {
      model <- ifelse(length(unique(Y)) ==2, "logistic","linear")
    }
  }
  if(model!="linear") opt.sigma <- FALSE
  if(opt.sigma){#Mir
    sigmainit <- sigmasq#Mir
    penaltiesinit <- penaltiesinit#Mir
  }#Mir
  response <- Y
  npen <- length(XXblocks)
  if(opt.sigma & optmethod=="Brent") optmethod <- "Nelder-Mead"


  if(!is.null(fixedpen)) {
    if(length(fixedpen) > (npen-1) | max(fixedpen) > npen){
      message("ERROR: incorrect input for argument \'fixedpen\'")
      return(NULL)
    }
  }

  if(fixedseed) set.seed(12345)
  if(is.null(sigmasq)&model!="linear") sigmasq <- 1#Mir


  if(is.null(penaltiesinit)) {penaltiesinit <- rep(100,npen); if(!is.null(pairing)) {penpair <- pairing[3]; penaltiesinit[penpair] <- 25}
  } else {if(!is.null(pairing)) {
    if(length(penaltiesinit) < npen)
    {
      penaltiesinit <- c(penaltiesinit,paired=NA)
      penpair1 <- pairing[1]; penpair2 <- pairing[2]; penpair3 <- pairing[3]
      penaltiesinit[penpair3] <- sqrt(penaltiesinit[penpair1])*sqrt(penaltiesinit[penpair2])*0.25
      penaltiesinit[penpair1] <- penaltiesinit[penpair1]*(1+.25)
      penaltiesinit[penpair2] <- penaltiesinit[penpair2]*(1+.25)
    }
  }
  }  #use smaller paired penalty by default

  print("Initial penalties:")
  print(penaltiesinit)
  multinit = rep(1,length(penaltiesinit))

  if(!is.null(fixedpen)){
    multinit <- multinit[-fixedpen]
  }

  if((length(multinit)+opt.sigma) == 1 & optmethod != "SANN") optmethod <- "Brent"  #Mir print(paste("Using",optmethod,"for optimization"))
  print(paste("Using",optmethod,"for optimization"))


  if(opt.sigma){
    CVS <- function(logmults){
      # sigmasq <- logmults[1]
      # logmults <- logmults[-1]
      allpenalties <- rep(NA,length(penaltiesinit))
      if(!is.null(fixedpen)) {
        allpenalties[fixedpen] <- penaltiesinit[fixedpen]
        allpenalties[-fixedpen] <- exp(logmults[-1])*penaltiesinit[-fixedpen]
        allpenalties <-c(exp(logmults[1]), allpenalties) #add sigma
      } else allpenalties <-  c(exp(logmults[1]),exp(logmults[-1])*penaltiesinit)

      res <- mgcv_lambda(allpenalties, XXblocks=XXblocks,Y=response,
                         model=model,printscore=tracescore,opt.sigma=TRUE)
      #allscores_mgcv <<- c(allscores_mgcv,res)
      #assign("allscores_mgcv", rbind(allscores_mgcv,c(res,allpenalties)), envir = baseenv())
      return(res)
    }

    #assign("allscores_mgcv", c(), envir = baseenv())  #saves all evaluations
    #minimal penalty
    lb <- log10(0.1/max(penaltiesinit))

    if(optmethod != "Brent"){
      optres <- optim(par = c(log(sigmainit),log(multinit)), CVS, control=list(reltol=reltol,maxit=maxItropt), method=optmethod)
    }else{
      optres <- try(optim(par = log(multinit), CVS, control=list(reltol=reltol,maxit=maxItropt), method=optmethod,lower=lb,upper=10))
      #optres <- try(optimise(CVS, tol=reltol,interval=c(-10,10)))
      }

    if(class(optres)=="try-error") optres <- optim(par = c(log(sigmainit),log(multinit)), CVS, control=list(reltol=reltol,maxit=maxItropt), method="Nelder-Mead")
    #if(!is.null(ncol(allscores_mgcv))) colnames(allscores_mgcv) <- c("Score", "sigmasq", paste("pen",1:npen,sep=""))

    optpen <- c(as.numeric(exp(optres$par[1])),as.numeric(exp(optres$par[-1]))*penaltiesinit)

    #optres2 <- c(optres,optpen=list(optpen),allsc = list(allscores_mgcv))
    optres2 <- c(optres,optpen=list(optpen))
    return(optres2)
  } else{#Mir: tot hier, in deze else loop staat wat het anders was geweest
    CVS <- function(logmults){
      #logmults <- rep(1,4)
      allpenalties <- rep(NA,length(penaltiesinit))
      if(!is.null(fixedpen)) {
        allpenalties[fixedpen] <- penaltiesinit[fixedpen]
        allpenalties[-fixedpen] <- c(exp(logmults[1]),exp(logmults[-1])*penaltiesinit[-fixedpen])
      } else allpenalties <- exp(logmults)*penaltiesinit
      res <- mgcv_lambda(allpenalties, XXblocks=XXblocks,Y=response,
                         model=model,printscore=tracescore,sigmasq=sigmasq, opt.sigma=FALSE)
      #assign("allscores_mgcv", rbind(allscores_mgcv,c(res,allpenalties)), envir = baseenv())
      return(res)
    }

    #assign("allscores_mgcv", c(), envir = baseenv())  #saves all evaluations

    if(optmethod != "Brent") optres <- optim(par = log(multinit), CVS, control=list(reltol=reltol,maxit=maxItropt), method=optmethod) else
      #optres <- try(optimise(CVS, tol=1,interval=c(-10,10)))
      optres <- try(optim(par = log(multinit), CVS, control=list(reltol=reltol,maxit=maxItropt), method=optmethod,lower=-5,upper=10))

    if(class(optres)=="try-error") optres <- optim(par = log(multinit), CVS, control=list(reltol=reltol,maxit=maxItropt), method="Nelder-Mead")
    #if(!is.null(ncol(allscores_mgcv))) colnames(allscores_mgcv) <- c("Score", paste("pen",1:npen,sep=""))

    optpen <- as.numeric(exp(optres$par))*penaltiesinit

    #optres2 <- c(optres,optpen=list(optpen),allsc = list(allscores_mgcv))
    optres2 <- c(optres,optpen=list(optpen))
    return(optres2)
  }
}

optLambdas_mgcvWrap <- function(penaltiesinit=NULL, XXblocks,Y, pairing=NULL, model=NULL,
                                reltol=1e-4, optmethod1= "SANN", optmethod2 =ifelse(length(penaltiesinit)==1,"Brent", "Nelder-Mead"),
                                maxItropt1=10,maxItropt2=25,tracescore=TRUE,fixedseed=TRUE, pref=NULL,fixedpen=NULL, sigmasq = 1,
                                opt.sigma=ifelse(model=="linear",TRUE, FALSE)){

  # penaltiesinit=NULL; X1=NULL; pairing=NULL; folds; intercept=TRUE;frac1=NULL;score="loglik"; model="logistic";
  # epsIWLS=1e-3;maxItrIWLS=25; traceCV=TRUE; reltol=1e-4; optmethod1= "SANN"; optmethod2 =ifelse(length(penaltiesinit)==1,"Brent", "Nelder-Mead");
  # maxItropt1=10;maxItropt2=25; save=FALSE; parallel=FALSE; warmstart=FALSE; pref=NULL;
  # penaltiesinit=lambdas; XXblocks=XXall;Y=resp;folds=leftout;score="loglik"; model="cox"; intercept=FALSE;
  # parallel=TRUE;warmstart=FALSE

  # penaltiesinit=lambdas; XXblocks=XXblocks;Y=respnum; pairing<- NULL;model="linear";reltol=1e-4; optmethod1= "SANN";maxItropt1=10;optmethod2= "Nelder-Mead";maxItropt2=25;
  # tracescore=TRUE;fixedseed=TRUE; pref=NULL;fixedpen=NULL; sigmasq = 1;opt.sigma=ifelse(model=="linear",TRUE, FALSE)
   if(is.null(model)){
    if(class(Y)=="Surv") model <- "cox" else {
      model <- ifelse(length(unique(Y)) ==2, "logistic","linear")
    }
  }
  if(model!="linear") opt.sigma <- FALSE #Mir
  if(opt.sigma){#Mir
    sigmainit <- sigmasq#Mir
    penaltiesinit <- penaltiesinit#Mir
  }#Mir
  if(is.null(sigmasq)&model!="linear") sigmasq <- 1#Mir

  nbl <- length(XXblocks)

  if(is.null(pref)){prefer <- 1:nbl; peni <- penaltiesinit} else {prefer <- pref;
  if(is.null(penaltiesinit)) peni<-NULL else peni <- penaltiesinit[prefer]}

  jl0 <- try(optLambdas_mgcv(penaltiesinit=peni, XXblocks=XXblocks[prefer],Y=Y, pairing=pairing, model=model,
                             reltol=reltol,optmethod=optmethod1,maxItropt=maxItropt1,tracescore=tracescore,
                             fixedseed=fixedseed, fixedpen=fixedpen,sigmasq=sigmasq,opt.sigma=opt.sigma ))

  if(class(jl0) != "try-error"){
    if(opt.sigma){#Mir
      lambdas0 <- jl0$optpen[-1]#Mir
      sigma0 <- jl0$optpen[1]#Mir
    }else{#Mir
      lambdas0 <- jl0$optpen
      sigma0 <- sigmasq
    }#Mir
  }else lambdas0 <- peni

  jl1 <- optLambdas_mgcv(penaltiesinit=lambdas0, XXblocks=XXblocks[prefer],Y=Y, pairing=pairing, model=model,
                         reltol=reltol,optmethod=optmethod2,maxItropt=maxItropt2,tracescore=tracescore,
                         fixedseed=fixedseed, fixedpen=fixedpen,sigmasq=sigma0,opt.sigma=opt.sigma )

  lambdas1 <- jl1$optpen

  if(opt.sigma){#Mir
    lambdas1 <- jl1$optpen[-1]#Mir
    sigma1 <- jl1$optpen[1]#Mir
  }else{#Mir
    lambdas1 <- jl1$optpen
    sigma1 <- sigmasq
  }#Mir
  res <- c(res0=list(jl0),res1=list(jl1))
  lambdas <- c(list(lambdas0),list(lambdas1))
  sigmas <- c(sigma0,sigma1)

  if(is.null(pref)) {
    return(list(res=res,lambdas=lambdas,optpen = lambdas1,sigmas=sigmas))
  } else {
    pref <- c(fixedpen,pref)
    lambdasall <- rep(NA,nbl)
    lambdasall[pref] <- lambdas1
    if(!is.null(penaltiesinit)) lambdasall[-pref] <- penaltiesinit[-pref] else lambdasall[-pref] <- 100
    jl2 <- try(optLambdas_mgcv(penaltiesinit=lambdasall, XXblocks=XXblocks,Y=Y, pairing=pairing, model=model,
                               reltol=reltol,optmethod=optmethod1,maxItropt=maxItropt1,tracescore=tracescore,
                               fixedseed=fixedseed, fixedpen=pref,sigmasq=sigmasq,opt.sigma=opt.sigma  ) )

    if(class(jl2) != "try-error"){
      if(opt.sigma){#Mir
        lambdas2 <- jl2$optpen[-1]#Mir
        sigma02 <- jl2$optpen[1]#Mir
      }else{#Mir
        lambdas2 <- jl2$optpen
        sigma02 <- sigmasq
      }#Mir
    }else lambdas2 <- lambdasall


    jl3 <- optLambdas_mgcv(penaltiesinit=lambdas2, XXblocks=XXblocks,Y=Y, pairing=pairing, model=model,
                           reltol=reltol,optmethod=optmethod2,maxItropt=maxItropt2,tracescore=tracescore,
                           fixedseed=fixedseed, fixedpen=pref,sigmasq=sigma02,opt.sigma=opt.sigma )

    if(opt.sigma){#Mir
      lambdas3 <- jl3$optpen[-1]#Mir
      sigma3 <- jl3$optpen[1]#Mir
    }else{#Mir
      lambdas3 <- jl3$optpen
      sigma3 <- sigmasq
    }#Mir
    res <- c(res,res2=list(jl2),res3=list(jl3))
    lambdas <- c(lambdas,list(lambdas2),list(lambdas3))
    sigmas <- c(sigmas,sigma02,sigma3)
    return(list(res=res,lambdas=lambdas,optpen = lambdas3,sigmas=sigmas))
  }
}

mlikCV <- function(penaltiesinit,XXblocks,Y,pairing=NULL, outfold=5, nrepeatout=1,balance=TRUE, fixedfolds=TRUE, model=NULL,
                   intercept=ifelse(class(Y)=="Surv", FALSE, TRUE),
                   reltol=1e-4, trace=FALSE, optmethod1= "SANN", optmethod2 =ifelse(length(penaltiesinit)==1,"Brent", "Nelder-Mead"),
                   maxItropt1=10,maxItropt2=25,parallel=FALSE, pref=NULL,fixedpen=NULL, sigmasq = 1,
                   opt.sigma=ifelse(model=="linear",TRUE, FALSE)){
  #1. produce outfolds
  #2. subset XXblocks and response
  #3. opt lambda per out-fold
  #4. Scoring

  # XXblocks <- XXblocksmm; response<-Y; X1=NULL;outfold=5;infold=10;nrepeatout=2; nrepeatin=1;balance=TRUE;fixedfolds=TRUE;
  # intercept=ifelse(class(Y)=="Surv", FALSE, TRUE); frac1=NULL;score="loglik"; model="logistic"; eps=1e-7;maxItr=100;trace=FALSE;printCV=TRUE;save=FALSE;parallel=F;
  #reltol=1e-4; optmethod=ifelse(length(XXblocks)==1,"Brent", "Nelder-Mead")
  if(is.null(model)){
    if(class(Y)=="Surv") model <- "cox" else {
      model <- ifelse(length(unique(Y)) ==2, "logistic","linear")
    }
  }
  response <- Y
  leftoutout <- CVfolds(Y=response,kfold=outfold,nrepeat=nrepeatout,fixedfolds=fixedfolds,balance=balance,model=model)
  n <- length(response)
  nbl <- length(XXblocks)
  npen <- length(penaltiesinit)
  if(npen != nbl ){
    message("ERROR: please make sure penaltiesinit is the same length as XXblocks.")
    return(NULL)
  }

  if(parallel) sfExport("leftoutout","XXblocks","response","n","fixedfolds","balance", "model","reltol","trace","optmethod1", "optmethod2","maxItropt1","maxItropt2")

  lpsoutfun <- function(i){
    if(parallel) sfCat(paste("Outer fold:",i),sep="\n") else print(paste("Outer fold:",i))
    #i <- 1
    outs <- leftoutout[[i]]
    ins <- (1:n)[-outs]
    respin <- response[ins]
    respout <- response[outs]
    XXblocksin <- lapply(XXblocks, function(block,whichin) block[whichin,whichin],whichin=ins)
    #X1in <- X1[ins,]
    #X1out <- X1[outs,]
    # foldsin <- CVfolds(Y=respin,kfold=infold,nrepeat=nrepeatin,fixedfolds=fixedfolds,balance=balance,model=model)
    # if(is.null(penaltiesinit)){
    # if(is.null(pairing)) penaltiesi <- sapply(1:nbl,function(bl) fastCV2(XXblocksin[[bl]],Y=respin, folds=foldsin, penaltyinit=NULL, X1=X1in,
    #                                                 intercept=intercept,frac1=frac1,score=score, model=model, epsIWLS=eps,
    #                                                 maxItrIWLS=maxItr, traceCV=trace, reltol=reltol, parallel=FALSE, NMinit=FALSE)[[1]])
    #   else { #initialization based on parametrization in the manuscript
    #        penaltiesi <- rep(NA,nbl)
    #        pair3 <- pairing[3]
    #        allbl <- setdiff(1:nbl,pair3)
    #        penaltiesi[allbl] <- sapply(allbl,function(bl) fastCV2(XXblocksin[[bl]],Y=respin, folds=foldsin, penaltyinit=NULL, X1=X1in,
    #                                                        intercept=intercept,frac1=frac1,score=score, model=model, epsIWLS=eps,
    #                                                        maxItrIWLS=maxItr, traceCV=trace, reltol=reltol, parallel=FALSE, NMinit=FALSE)[[1]])
    #        penpair1i <- penaltiesi[pairing[1]];penpair2i <- penaltiesi[pairing[2]]
    #        penpair3 <- sqrt(penpair1i)*sqrt(penpair2i)/4
    #        penaltiesi[pair3] <- penpair3
    #        penaltiesi[pairing[1]] <- penpair1i*(1+1/4);penaltiesi[pairing[2]] <- penpair2i*(1+1/4)
    #   }
    # } else
    #penaltiesinit <- c(4,5,6); pairing=c(1,2)
    penaltiesi <- penaltiesinit
    if(!is.null(pairing)) {
      penaltiesi <- penaltiesinit
      penaltiesi <- c(penaltiesi, sqrt(penaltiesi[pairing[1]])*sqrt(penaltiesi[pairing[2]])/4)
      penaltiesi[pairing] <-  penaltiesi[pairing]*(1+1/4)
    } else penaltiesi <- penaltiesinit

    optlam <- optLambdas_mgcvWrap(penaltiesinit=penaltiesi, XXblocks=XXblocksin,Y=respin,pairing=pairing,
                                  model=model, tracescore=trace, reltol=reltol,
                                  optmethod1=optmethod1,optmethod2=optmethod2, maxItropt1=maxItropt1,
                                  maxItropt2=maxItropt2,pref=pref,fixedpen=fixedpen, sigmasq = sigmasq,
                                  opt.sigma=opt.sigma)$optpen
    print("Opt penalties:");print(optlam)

    SigmaLambda <- SigmaFromBlocks(XXblocks,optlam, pairing=pairing)
    SigmaIn <- SigmaLambda[ins,ins] #Sigma matrix allows subsetting; SigmaIn is required for fitting
    SigmaOutIn <- SigmaLambda[outs,ins] #SigmaOutIn is required for prediction
    if(model!="cox") fitin <- try(IWLSridge(SigmaIn,respin,X1=NULL,intercept=intercept,trace=FALSE))
    if(model=="cox") fitin <- try(IWLSCoxridge(SigmaIn,respin,X1=NULL,intercept=intercept,trace=FALSE))
    if(class(fitin) =="try-error") lpsout <- rep(0,length(respout)) else lpsout <- predictIWLS(fitin, X1new=NULL,SigmaOutIn)
    return(lpsout)
  }
  nf <- length(leftoutout)
  if(parallel) lpsoutlist <- sfLapply(1:nf,lpsoutfun) else lpsoutlist <- lapply(1:nf,lpsoutfun)
  all_lpsout <- unlist(lpsoutlist)
  #all_respout <- unlist(lapply(1:nf,function(i){out <- leftoutout[[i]];respout <- response[out];return(respout)}))
  if(model!="cox") all_respout <- unlist(lapply(1:nf,function(i){out <- leftoutout[[i]];respout <- response[out];return(respout)})) else {  #cox model
    all_evout <- unlist(lapply(1:nf,function(i){out <- leftoutout[[i]];respout <- response[out,2];return(respout)}))
    all_timeout <- unlist(lapply(1:nf,function(i){out <- leftoutout[[i]];respout <- response[out,1];return(respout)}))
    all_respout <- Surv(all_timeout,all_evout)
  }
  all_samout <- unlist(lapply(1:nf,function(i){out <- leftoutout[[i]];return(out)}))
  whfold <- unlist(lapply(1:nf,function(i){out <- leftoutout[[i]];return(rep(i,length(out)))}))
  res0 <- data.frame(sampleindex=all_samout,true=all_respout, linpred = all_lpsout, whfold=whfold)
  od <- order(all_samout)
  res <- res0[od,]
  # rn <- rownames(XXblocks[[1]])
  # if(!is.null(rn)) rownames(res) <- sapply(all_samout,function(sam) rn[sam])
  return(res)
}





