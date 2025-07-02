# Functions for fitting a Hawkes discrete HMM

################################################################################
# Utils
################################################################################
MakeDiscreteHawkes1D <- function(contPath, N, Tmax=1){
  # contPath = list(times, states)
  # Simulation of a 1D discrete Hawkes process
  timesDisc <- table(ceiling(N*contPath$times/Tmax))
  y <- rep(0, N)
  y[as.numeric(names(timesDisc))] <- timesDisc
  jumpsDiffDisc <- diff(c(floor(N*contPath$hiddenPath$times/Tmax), N))
  z <- as.vector(sapply(1:length(contPath$hiddenPath$states), function(i){
    rep(contPath$hiddenPath$states[i], each=jumpsDiffDisc[i])
  }))
  return(list(y=y, z=z))
}

################################################################################
Parms2Vec <- function(parms){
  # Gathers all parms in a single vector, ordering the transition matrix by row
  tmp <- parms; tmp$pi <- t(tmp$pi)
  return(unlist(tmp))
}

################################################################################
Parms2VecEmission <- function(parms){
  return(c(parms$alpha, parms$beta, parms$mu))
}

################################################################################
Vec2ParmsEmission <- function(vecEmission, Q){
  # Splits the parameter vector into a list
  names(vecEmission) <- c()
  alpha <- vecEmission[1]
  beta <- vecEmission[2]
  mu <- vecEmission[3:(2+Q)]
  return(list(alpha=alpha, beta=beta, mu=mu))
}

################################################################################
ABM2AlphaBetaMu <- function(abm, delta){
  mu <- abm[3] * delta
  beta <- exp(-abm[2]*delta)
  alpha <- abm[1] * (1/beta - 1) / abm[2]
  alphaBetaMu <- c(alpha, beta, mu); names(alphaBetaMu) <- c('alpha', 'beta', 'mu')
  return(alphaBetaMu)
}

################################################################################
# Function for the homogeneous Hawkes in discrete time
################################################################################
GradLogLikDiscreteHawkes1D <- function(parms, data, w=rep(1, length(data$y)), logParms=FALSE){
  # Gradient of the log-lik of a discrete Hawkes process
  if(logParms){parms <- exp(parms)}
  alpha <- parms[1]; beta <- parms[2]; mu <- parms[3]
  N <- length(data$y) 
  u <- 0
  gradu <- rep(0, 3)
  a <- data$y[1]/(mu + beta*u) - 1
  grad <- w[1] * a * c(beta*gradu[1], u+beta*gradu[2], 1)
  for(t in 2:N){
    gradu <- c(data$y[t-1] + beta*gradu[1], u+beta*gradu[2], 0)
    u <- alpha*data$y[t-1] + beta*u
    a <- data$y[t]/(mu + beta*u) - 1
    grad <- grad + w[t] * a * c(beta*gradu[1], u+beta*gradu[2], 1)
  }
  names(grad) <- c('alpha', 'beta', 'mu')
  if(logParms){grad <- grad*parms}
  return(grad)
}

################################################################################
LogLikDiscreteHawkes1D <- function(parms, data, w=rep(1, length(data$y)), logParms=FALSE){
  # Log-likelihood of a 1D discrete Hawkes process
  if(logParms){parms <- exp(parms)}
  alpha <- parms[1]; beta <- parms[2]; mu <- parms[3]
  N <- length(data$y) 
  u <- 0; 
  logLik <- w[1]*dpois(data$y[1], mu, log=TRUE)
  for(t in 2:N){
    u <- alpha*data$y[t-1] + beta*u
    logLik <- logLik + w[t]*dpois(data$y[t], mu + beta*u, log=TRUE)
  }
  return(logLik)
}

################################################################################
# Initialization for the discrete-time Hawkes HMM (hawkesbow + Poisson HMM)
################################################################################
HawkesBow2ABM <- function(bow){
  abm <- c(bow[2]*bow[3], bow[3], bow[1])
  names(abm) <- c('a', 'b', 'm')
  return(abm)
}

########################################################
StatDist <- function(pi)
{
  Eig <- eigen(t(pi))
  q <- which.min(abs(Eig$values-1))
  nu <- Re(as.vector(Eig$vectors[,q]))
  nu <- nu / sum(nu)
  return(nu)
}

########################################################
PoissonHMMuniInit <- function(y, Q){
  n <- length(y)
  tau <- Mclust(sqrt(y)+sqrt(1+y), G=Q, modelNames='E')$z
  eta <- t(tau[1:(n-1), ])%*%tau[2:n, ]
  return(list(tau=tau, eta=eta))
}

########################################################
Forward <- function(pi, nu, logPhi)
{
  # Entropy computation: see Hernando & al (05)
  n <- nrow(logPhi); Q <- ncol(logPhi)
  phi <- exp(logPhi)
  Fw <- matrix(0, n, Q); H <- Fw; 
  
  logL <- 0
  Fw[1,] <- phi[1,] * nu
  logL <- logL + log(sum(Fw[1,]))
  Fw[1,] <- Fw[1,] / sum(Fw[1,])
  for (t in (2:n)){
    Ptmp <- Fw[t-1,] %*% pi
    Fw[t,] <- (Ptmp) * phi[t,]
    logL <- logL + log(sum(Fw[t,]))
    Fw[t,] <- Fw[t,] / sum(Fw[t,])
    Ptmp <- Ptmp / sum(Ptmp)
    H[t, ] <- (H[t-1, ] * Ptmp) - (Ptmp * log(Ptmp))
  }
  Ent <- sum(H[n,] * Fw[n]) - sum(Fw[n, ]*log(Fw[n,])) 
  return(list(Fw=Fw, logL=logL, Ent=Ent))
}

########################################################
Backward <- function(pi, Fw)
{
  n <- dim(Fw)[1]
  Q <- dim(Fw)[2]
  # tau(t,q) <- E(z(t,q) | X)
  tau <- matrix(0, n, Q)
  G <- tau
  tau[n,] <- Fw[n, ]
  # eta(q, l) <- sum_t E(z(t, q) z(t+1, l) | X)
  eta <- etaTmp <- matrix(0, Q, Q)
  condEntropy <- 0
  for (t in ((n-1):1)){
    G[t+1, ] <- Fw[t,] %*% pi
    B <- as.vector(tau[t+1, ] / G[t+1, ])
    tau[t,] <- Fw[t,] * (pi %*% B)
    tau[t,] <- tau[t,] / sum(tau[t,])
    etaTmp <- pi * outer(Fw[t,], B)
    piTmp <- etaTmp / tau[t,]
    condEntropy <- condEntropy - sum(etaTmp * log(piTmp + (piTmp==0)))
    eta <- eta + etaTmp
  }
  condEntropy <- condEntropy - sum(tau[1,] * log(tau[1, ] + (tau[1, ]==0)))
  return(list(tau=tau, eta=eta, G=G, condEntropy=condEntropy))
}

########################################################
HMMEstep <- function(pi, nu, logPhi, epsTau=1e-6)
{
  Fw <- Forward(pi, nu, logPhi)
  Bw <- Backward(pi, Fw$Fw)
  logL <- Fw$logL
  Ent <- Fw$Ent
  Fw <- Fw$Fw
  tau <- Bw$tau
  G <- Bw$G
  eta <- Bw$eta
  
  # Smoothing
  tau <- tau + epsTau;   tau <- tau / rowSums(tau)
  
  return(list(tau=tau, Fw=Fw, G=G, eta=eta, logL=logL, Ent=Ent, condEntropy=Bw$condEntropy))
}

########################################################
PoissonHMMuniMstep <- function(y, tau, eta, offset=rep(0, length(y))){
  Q <- ncol(tau); n <- length(y)
  pi <- eta/rowSums(as.matrix(eta))
  nu <- StatDist(pi)
  Gamma = t(y)%*%tau / colSums(tau)
  logPhi = matrix(0, n, Q); 
  for (q in 1:Q){logPhi[, q] = dpois(y, Gamma[q]+offset, log=TRUE)}
  return(list(pi=pi, nu=nu, Gamma=Gamma, logPhi=logPhi))
}

#######################################################
# Viterbi path
Viterbi <- function(phi, nu, pi)
{
  n <- dim(phi)[1]
  Q <- dim(phi)[2]
  # Prob of the most probable path
  lik <- matrix(0, n, Q)
  transMat <- lik
  lik[1,] <- phi[1,] * nu
  logLikTot <- log(sum(lik[1, ]))
  lik[1, ] <- lik[1, ] / sum(lik[1, ])
  for (t in (2:n)){
    likTrans <- matrix(rep(lik[t-1,], each=Q), Q, Q) * pi
    lik[t,] <- phi[t,] * apply(likTrans, 1, max)
    logLikTot <- logLikTot + log(sum(lik[t, ]))
    lik[t, ] <- lik[t, ] / sum(lik[t, ])
    transMat[t, ] <- apply(likTrans, 1, which.max)   
  }
  # Most probable path
  z <- rep(0, n)
  z[n] <- which.max(lik[n,])
  for (t in ((n-1):1)){
    z[t] <- transMat[(t+1), z[t+1]]
  }
  likTot <- exp(logLikTot)
  
  return(list(z=z, logLik=logLikTot))
}

########################################################
PoissonHMMuniEM <- function(y, Q, offset=rep(0, length(y)), epsTol=1e-6, maxIter=1e3){
  n <- length(y)
  if (Q == 1){
    tau <- matrix(1, length(y), 1); eta <- 1
  }else{
    z <- ceiling(Q*rank(y-offset)/n); tau <- matrix(0, n, Q)
    for(q in 1:Q){tau[which(z==q), q] <- 1}; tau <- tau + 1/Q; tau <- tau / rowSums(tau)
    eta <- t(tau[-n, ])%*%tau[-1, ]
    tau <- tau; eta <- eta
  }
  diff <- 2*epsTol; iter <- 0
  while((diff > epsTol) & (iter <= maxIter)){
    iter=iter+1
    mStep <- PoissonHMMuniMstep(y, tau, eta, offset)
    eStep <- HMMEstep(mStep$pi, mStep$nu, mStep$logPhi)
    diff <- max(abs(tau - eStep$tau))
    tau <- eStep$tau
  }
  viterbi <- Viterbi(phi=exp(mStep$logPhi), nu=mStep$nu, pi=mStep$pi)$z
  return(list(pi=mStep$pi, nu=mStep$nu, Gamma=mStep$Gamma, tau=eStep$tau, 
              logPhi=mStep$logPhi, logL=eStep$logL, Ent=eStep$Ent, 
              classif=apply(eStep$tau, 1, which.max), iter=iter, eStep=eStep, viterbi=viterbi))
}

################################################################################
# Init EM (could be with small Ninit < N)
InitEMHMMHawkes1D <- function(contPath, nInit, Q, Tmax=1, 
                              controlNloptr=list(xtol_rel=1e-9, maxeval=1e3)){
  deltaInit <- Tmax/nInit
  bow <- hawkesbow::mle(contPath$times, "Exponential", end=1, opts=controlNloptr)$par
  bowCont <- HawkesBow2ABM(bow)
  discPathInit <- MakeDiscreteHawkes1D(contPath=contPath, N=nInit, Tmax=Tmax)
  bowDiscInit <- ABM2AlphaBetaMu(bowCont, delta=deltaInit)
  if(Q==1){
    gamma <- mean(discPathInit$y)
    logPhi <- dpois(discPathInit$y, gamma, log=TRUE)
    hmm <- list(pi=matrix(1, 1, 1), nu=rep(1, 1), Gamma=gamma, 
                tau=matrix(1, length(discPathInit$y), 1),
                logPhi=logPhi, logL=sum(logPhi), Ent=0, iter=1)
  }else{hmm <- PoissonHMMuniEM(y=discPathInit$y, Q=Q)}
  hmm$mu <- as.vector(hmm$Gamma*(1 - bowCont[1]/bowCont[2]))
  adhoc <- list(nu=rep(1/Q, Q), pi=matrix(1, Q, Q) + (N/10)*diag(Q)); 
  adhoc$pi <- adhoc$pi/rowSums(adhoc$pi)
  parmsDiscInit <- list(alpha=bowDiscInit[1], beta=bowDiscInit[2], mu=hmm$mu, nu=adhoc$nu, pi=adhoc$pi, tau=hmm$tau, hmm=hmm)
  return(parmsDiscInit)
}

################################################################################
# E step for the discrete-time Hawkes HMM EM algorithm
################################################################################
ForwardHawkesHmm1D <- function(y, parms){
  # Initialization
  N <- length(y); Q <- length(parms$mu); 
  u <- rep(0, N); phi <- matrix(0, N, Q)
  for(t in 2:N){u[t] <- parms$beta*u[t-1] + parms$alpha*y[t-1]}
  ent <- fwd <- matrix(0, N, Q); 
  fwd[1, ] <- phi[1, ] <- dpois(y[1], parms$mu)
  logL <- log(sum(fwd[1, ]))
  fwd[1, ] <- fwd[1, ] / sum(fwd[1, ])
  # Forward recursion
  for(t in 2:N){
    pTmp <- fwd[t-1,] %*% parms$pi
    phi[t, ] <- dpois(y[t], parms$mu + parms$beta*u[t])
    fwd[t, ] <- (fwd[t-1, ]%*%parms$pi) * phi[t, ]
    pCondy <- sum(fwd[t, ])
    logL <- logL + log(pCondy)
    fwd[t, ] <- fwd[t, ] / pCondy
    pTmp <- pTmp / rowSums(pTmp)
    ent[t, ] <- (ent[t-1, ] * pTmp) - (pTmp * log(pTmp))
  }
  ent <- sum(ent[N,] * fwd[N]) - sum(fwd[N, ]*log(fwd[N,]))
  phi[which(phi < 1e-30)] <- 1e-30
  return(list(fwd=fwd, logL=logL, phi=phi, ent=ent))
}

################################################################################
BackwardHawkesHmm1D <- function(parms, forward){
  # Initialization
  N <- nrow(forward$fwd); Q <- length(parms$mu)
  tau <- G <- matrix(0, N, Q)
  tau[N, ] <- forward$fwd[N, ]
  # Backward recursion
  eta <- matrix(0, Q, Q)
  condEnt <- 0
  for (t in ((N-1):1)){
    G[t+1, ] <- forward$fwd[t, ] %*% parms$pi
    B <- as.vector(tau[t+1, ] / G[t+1, ])
    tau[t,] <- forward$fwd[t,] * (parms$pi %*% B)
    tau[t,] <- tau[t,] / sum(tau[t,])
    # eta <- eta + parms$pi * outer(forward$fwd[t,], B)
    etaTmp <- parms$pi * outer(forward$fwd[t,], B)
    piTmp <- etaTmp / tau[t,]
    condEnt <- condEnt - sum(etaTmp * log(piTmp + (piTmp==0)))
    eta <- eta + etaTmp
  }
  condEnt <- condEnt - sum(tau[1,] * log(tau[1, ] + (tau[1, ]==0)))
  return(list(tau=tau, eta=eta, condEnt=condEnt))
}
Estep <- function(y, parms){
  forward <- ForwardHawkesHmm1D(y, parms) 
  backward <- BackwardHawkesHmm1D(parms, forward)
  return(list(fwd=forward$fwd, logL=forward$logL, tau=backward$tau, eta=backward$eta, 
              phi=forward$phi, ent=forward$ent, condEnt=backward$condEnt))
}

################################################################################
# M step for the discrete-time Hawkes HMM EM algorithm
################################################################################
EspLogLikEmission <- function(data, parms, eStep){
  Q <- ncol(eStep$tau)
  espLogLik <- 0
  for(q in 1:Q){
    espLogLik <- espLogLik + 
      LogLikDiscreteHawkes1D(parms=c(alpha=parms$alpha, beta=parms$beta, mu=parms$mu[q]), 
                             data=data, w=eStep$tau[, q], logParms=FALSE)
  }
  return(espLogLik)
}
EspLogLikEmissionMu <- function(data, mu, eStep, alphaBeta){
  return(EspLogLikEmission(data=data, parms=Vec2ParmsEmission(c(alphaBeta, mu), Q=ncol(eStep$tau)), eStep=eStep))
}
EspLogLikEmissionParmsVec <- function(data, parmsVec, eStep){
  return(EspLogLikEmission(data=data, parms=Vec2ParmsEmission(parmsVec, Q=ncol(eStep$tau)), eStep=eStep))
}
EspLogLik <- function(data, parms, eStep){
  return(EspLogLikTransition(data, parms, eStep) + EspLogLikEmission(data, parms, eStep))
}

################################################################################
# Gradients
GradEspLogLikEmission <- function(data, parms, eStep){
  # Gradient for the emission parms: alpha, beta, mu[q]
  Q <- ncol(eStep$tau)
  gradAlpha <- gradBeta <- 0; gradMu <- rep(0, Q)
  for(q in 1:Q){
    gradTmp <- GradLogLikDiscreteHawkes1D(parms=c(alpha=parms$alpha, beta=parms$beta, mu=parms$mu[q]), 
                                          data=data, w=eStep$tau[, q], logParms=FALSE)
    gradAlpha <- gradAlpha + gradTmp[1]
    gradBeta <- gradBeta + gradTmp[2]
    gradMu[q] <- gradTmp[3]
  }
  grad <- c(gradAlpha, gradBeta, gradMu)
  names(grad) <- c('alpha', 'beta', paste0('mu', 1:Q))
  return(grad)
}
GradEspLogLikEmissionParmsVec <- function(data, parmsVec, eStep){
  return(GradEspLogLikEmission(data=data, parms=Vec2ParmsEmission(parmsVec, Q=ncol(eStep$tau)), eStep=eStep))
}
GradEspLogLikEmissionMu <- function(data, mu, eStep, alphaBeta){
  return(GradEspLogLikEmission(data=data, parms=Vec2ParmsEmission(c(alphaBeta, mu), Q=ncol(eStep$tau)), eStep=eStep)[-c(1:2)])
}

################################################################################
Mstep <- function(data, parms, eStep){
  nu <- eStep$tau[1, ]; Q <- length(nu)
  pi <- eStep$eta / rowSums(eStep$eta)
  
  parmsTransition <- c(nu, as.vector(t(pi)))
  logParmsEmission <- optim(par=log(parmsEmissionInit),
                            fn=EspLogLikLogParms, gr=GradEspLogLikLogParms,
                            parmsTransition=parmsTransition, y=y, eStep=eStep,
                            method='BFGS', control=list(fnscale=-1))$par
  # Re-ordering states
  logMu <- logParmsEmission[3:(2+Q)]; orderMu <- order(logMu)
  nu <- nu[orderMu]; pi <- pi[orderMu, orderMu]
  logParmsEmission <- c(logParmsEmission[1:2], logMu[orderMu])
  parmsEmission <- exp(logParmsEmission)
  parmsTransition <- c(nu, as.vector(t(pi)))
  return(Vec2Parms(c(parmsEmission, parmsTransition), Q))
}


################################################################################
# Fit discrete-time Hawkes HMM
################################################################################
# EM 
EMHMMHawkes1D <- function(discPath, discParms, Q, iterMax=1e3, tol=1e-4, print=FALSE, 
                          controlOptim=list(maxit=1e3, fnscale=-1, trace=0, factr=1e-9, pgtol=1e-9)){
  parmsDiscInit <- discParms
  eStep <- Estep(y=discPath$y, parms=parmsDiscInit)
  logLpath <- eStep$logL
  tau <- eStep$tau
  iter <- 0; diff <- 2*tol
  cat('iter = ')
  while((iter < iterMax) & (diff > tol)){
    iter <- iter+1
    # E step + transition
    eStep <- Estep(y=discPath$y, parms=discParms)
    parmsTransition <- list(nu=eStep$tau[1, ], pi=eStep$eta/rowSums(eStep$eta))
    # Path
    logL <- eStep$logL; logLpath <- c(logLpath, logL)
    # M step 
    fitEmission <- optim(par=Parms2VecEmission(discParms), fn=EspLogLikEmissionParmsVec, gr=GradEspLogLikEmissionParmsVec, 
                         data=discPath, eStep=eStep, method='L-BFGS-B', lower=rep(1e-6, 2+Q), upper=c(Inf, 1-1e-6, rep(Inf, Q)), 
                         control=controlOptim)
    parmsEmission <- Vec2ParmsEmission(fitEmission$par, Q=Q)
    parmsNew <- list(alpha=parmsEmission$alpha, beta=parmsEmission$beta, mu=parmsEmission$mu, 
                     nu=parmsTransition$nu, pi=parmsTransition$pi)
    # Update
    if(iter > 1){diff <- max(abs(eStep$tau - tau))}
    if(print){
      cat(iter, discParmsPath[iter, 1:(2+Q)], logLpath[iter], diff, '')
    }else{if(iter%%round(sqrt(iterMax))==0){cat(iter, '')}}
    discParms <- parmsNew
    tau <- eStep$tau
  }
  cat(iter, '\n')
  classif <- apply(eStep$tau, 1, which.max)
  viterbi <- Viterbi(eStep$phi, parmsTransition$nu, parmsTransition$pi)$z
  return(list(discParms=discParms, 
              parmsTransition=parmsTransition, parmsEmission=parmsEmission,
              logL=logL, logLpath=logLpath, classif=classif, viterbi=viterbi, 
              eStep=eStep))
}

