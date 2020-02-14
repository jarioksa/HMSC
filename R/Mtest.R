#' @title Matrix Method Tester for Posterior Sampling
#'
#' @description Tests the use of \CRANpkg{Matrix} for sampling the
#'     posterior with block-conditional Gibbs MCMC sampler. The
#'     function is based on \code{\link{sampleMcmc}}, but only calls
#'     once the specified samplers which were modified to use
#'     \CRANpkg{Matrix} functions instead of \pkg{base} and
#'     \pkg{stats} functions. The results must remain unchanged after
#'     moving to \CRANpkg{Matrix} and we see if they are any faster.
#'
#'
#' @param hM a fitted \code{Hmsc} model object
#' @param initPar a named list of parameter values used for initiation
#'     of MCMC states
#' @param dataParList a named list with pre-computed \code{Qg},
#'     \code{iQg}, \code{RQg}, \code{detQg}, \code{rLPar} parameters
#' @param updater a named list, specifying which conditional updaters
#'     should be ommitted
#' @param fromPrior whether prior (TRUE) or posterior (FALSE) is to be
#'     sampled
#'
#' @return Sampled list of updated parameters.
#'
#' @seealso \code{\link{sampleMcmc}}, \code{\link{Hmsc}}
#'
#' @examples
#' data(TD)
#' set.seed(4711)
#' system.time(out <- Mtest(TD$m))
#' out
#' rm(.Random.seed)
#'
#' @export

`Mtest` <-
    function(hM, initPar=NULL, dataParList=NULL, updater=list(),
             fromPrior = FALSE)
{
    ## We do only one iteration (and one chain): set sampleMcmc para here
    samples <- 1L
    nChains <- 1L
    chain <- 1L
    nParallel <- 1L
    transient <- 0L
    thin <- 1L
    verbose <- FALSE
    adaptNf <- rep(transient, hM$nr)
    alignPost <- FALSE

   X1 = hM$XScaled
   if (hM$ncsel>0){
      if(is.matrix(X1)){
         X2=X1
         X1=list()
         for (j in 1:hM$ns){
            X1[[j]] = X2
         }
      }
   }

   Tr = hM$TrScaled
   Y = hM$YScaled
   distr = hM$distr
   Pi = hM$Pi
   dfPi = hM$dfPi
   C = hM$C
   nr = hM$nr

   mGamma = hM$mGamma
   iUGamma = chol2inv(chol(hM$UGamma))
   V0 = hM$V0
   f0 = hM$f0
   aSigma = hM$aSigma
   bSigma = hM$bSigma
   rhopw = hM$rhopw

   if(is.null(dataParList))
      dataParList = computeDataParameters(hM)
   Qg = dataParList$Qg
   iQg = dataParList$iQg
   RQg = dataParList$RQg
   detQg = dataParList$detQg
   rLPar = dataParList$rLPar

   hM$postList = vector("list", nChains)
   hM$repList = vector("list", nChains)
   initSeed = sample.int(.Machine$integer.max, nChains)

   ## switching of the augmented updaters if the required conditions
   ## are not met
   EPS = 1e-6
   updaterWarningFlag = TRUE
   ## updater$Gamma2
   if(!identical(updater$Gamma2, FALSE) && any(abs(mGamma) > EPS)) {
      updater$Gamma2 = FALSE
      if(updaterWarningFlag)
         message("Setting updater$Gamma2=FALSE due to non-zero mGamma")
   }
   if(!identical(updater$Gamma2, FALSE) &&
      any(abs(iUGamma - kronecker(iUGamma[1:hM$nc,1:hM$nc],
                                  diag(hM$nt))) > EPS)){
      updater$Gamma2 = FALSE
      if(updaterWarningFlag)
         message("Setting updater$Gamma2=FALSE due to non-kronecker structure of UGamma matrix")
   }
   if(!identical(updater$Gamma2, FALSE) && (!is.null(C))){
      updater$Gamma2 = FALSE
      if(updaterWarningFlag)
         message("Setting updater$Gamma2=FALSE due to specified phylogeny matrix")
   }
   ## updater$GammaEta
   if(!identical(updater$GammaEta, FALSE) && any(abs(mGamma) > EPS)){
      updater$GammaEta = FALSE
      if(updaterWarningFlag)
          message("Setting updater$GammaEta=FALSE due to non-zero mGamma")
   }
   if(!identical(updater$GammaEta, FALSE) && hM$nr==0){
      updater$GammaEta = FALSE
      if(updaterWarningFlag)
         message("Setting updater$GammaEta=FALSE due to absence of random effects included to the model")
   }
    if(nChains>1)
         print(sprintf("Computing chain %d", chain))
    set.seed(initSeed[chain])
    parList = computeInitialParameters(hM, initPar)

    Gamma = parList$Gamma
    V = parList$V
    iV = chol2inv(chol(V))
    Beta = parList$Beta
    BetaSel = parList$BetaSel
    PsiRRR = parList$PsiRRR
    DeltaRRR = parList$DeltaRRR
    wRRR = parList$wRRR
    sigma = parList$sigma
    iSigma = 1 / sigma
    Lambda = parList$Lambda
    Eta = parList$Eta
    Alpha = parList$Alpha
    Psi = parList$Psi
    Delta = parList$Delta
    rho = parList$rho
    Z = parList$Z

    X1A = X1

    if(hM$ncsel>0){
        for (i in 1:hM$ncsel){
            XSel = hM$XSelect[[i]]
            for (spg in 1:length(XSel$q)){
                if(!BetaSel[[i]][spg]){
                    fsp = which(XSel$spGroup==spg)
                    for (j in fsp){
                        X1A[[j]][,XSel$covGroup]=0
                    }
                }
            }
        }
    }

    X = X1A
    if(hM$ncRRR>0){
        XB=hM$XRRRScaled%*%t(wRRR)
        if(is.matrix(X)){
            X=cbind(X,XB)
        } else {
            for (j in 1:hM$ns){
                X[[j]] = cbind(X[[j]],XB)
            }
        }
    }

    if(!identical(updater$Gamma2, FALSE) && (!is.matrix(X))){
        updater$Gamma2 = FALSE
        if(updaterWarningFlag)
            message("Setting updater$Gamma2=FALSE due to X is not a matrix")
    }
    if(!identical(updater$GammaEta, FALSE) && (!is.matrix(X))){
        updater$GammaEta = FALSE
        if(updaterWarningFlag)
            message("Setting updater$GammaEta=FALSE due to X is not a matrix")
    }

    postList = vector("list", samples)
    for(iter in seq_len(transient+samples*thin)){
        if(!identical(updater$Gamma2, FALSE))
            Gamma = updateGamma2(Z=Z,Gamma=Gamma,iV=iV,iSigma=iSigma,
                                 Eta=Eta,Lambda=Lambda, X=X,Pi=Pi,
                                 dfPi=dfPi,Tr=Tr,C=C,rL=hM$rL,
                                 iQg=iQg, mGamma=mGamma,
                                 iUGamma=iUGamma)

        if(!identical(updater$GammaEta, FALSE)){
            GammaEtaList <-
                updateGammaEta(Z=Z, Gamma=Gamma,
                               V=chol2inv(chol(iV)), iV=iV,
                               id=iSigma, Eta=Eta, Lambda=Lambda,
                               Alpha=Alpha, X=X, Pi=Pi, dfPi=dfPi,
                               Tr=Tr, rL=hM$rL, rLPar=rLPar,
                               Q=Qg[,,rho], iQ=iQg[,,rho],
                               RQ=RQg[,,rho], U=hM$UGamma,
                               iU=iUGamma)
             Gamma = GammaEtaList$Gamma
             Eta = GammaEtaList$Eta
        }
    }
    ## stop here!
    return(list(Gamma = Gamma, Eta = Eta))
    ## skip the rest
}
