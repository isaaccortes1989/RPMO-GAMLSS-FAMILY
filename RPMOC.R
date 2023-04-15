#Comments:
#pRPMOC: cdf for the RPMOC model
#dRPMOC: pdf for the RPMOC model
#qRPMOC: quantile function for the RPMOC model 
#rRPMOC: function to generate random variables that follow a RPMOC distribution  
#mu is the quantile parameter
#sigma is the scale parameter
#nu is the skewness parameter
#tau is the parameter of the fixed probability

rm(list=ls())

pRPMOC = function(q,mu=1,sigma=1,nu=1,tau=0.5,lower.tail=TRUE,log.p=FALSE){
  if(any(sigma<=0)) stop(paste("sigma must be positive","\n",""))
  if(any(nu<=0)) stop(paste("nu must be positive","\n",""))
  if(any(tau<=0)|any(tau>=1)) stop(paste("tau must be between 0 and 1","\n",""))
  kappa_1 <- ((nu*tau)/(1+tau*(nu-1)))
  t       <- ((q-mu)/sigma)+qcauchy(kappa_1)
  lGy     <- log(pcauchy(t))-log(((1-nu)*pcauchy(t)+nu)) 
  if(log.p == FALSE)
    cdf <- exp(lGy)  
  else cdf <- lGy
  if(lower.tail == TRUE)
    cdf <- cdf
  else cdf <- 1-cdf
  return(cdf)
}
###############
dRPMOC = function(x,mu=1,sigma=1,nu=1,tau=0.5,log=FALSE){
  if(any(sigma<=0)) stop(paste("sigma must be positive","\n",""))
  if(any(nu<=0))    stop(paste("nu must be positive","\n",""))
  if(any(tau<=0)|any(tau>=1)) stop(paste("tau must be between 0 and 1","\n",""))
  kappa_1 <- (nu*tau)/(1+tau*(nu-1))
  t       <- ((x-mu)/sigma)+qcauchy(kappa_1)
  lfy     <- log(nu)+dcauchy(t,log=TRUE)-log(sigma)-2*log(nu+(1-nu)*pcauchy(t))
  if(log== FALSE)
    fy <- exp(lfy)
  else fy <- lfy
  return(fy)
}
###############
qRPMOC = function(p,mu=1,sigma=1,nu=1,tau=0.5,lower.tail=TRUE,log.p=FALSE){ 
  if(any(sigma<=0)) stop(paste("sigma must be positive","\n",""))
  if(any(nu<=0)) stop(paste("nu must be positive","\n",""))
  if(any(tau<=0)|any(tau>=1)) stop(paste("tau must be between 0 and 1","\n",""))
  if(log.p == TRUE)
    p <- exp(p)  
  else p <- p
  if(lower.tail == TRUE)
    p <- p
  else p <- 1-p
  if(any(p<0)|any(p>1)) stop(paste("p must be between 0 and 1","\n",""))
  kappa_1 <- (nu*tau)/(1+tau*(nu-1))
  Q       <- sigma*(qcauchy((nu*p)/(1+p*(nu-1)))-qcauchy(kappa_1))+mu
  return(Q)
}
###############
rRPMOC = function(n,mu=1,sigma=1,nu=1,tau=0.5){
  if(any(sigma<=0)) stop(paste("sigma must be positive","\n",""))
  if(any(nu<=0)) stop(paste("nu must be positive","\n",""))
  if(any(tau<=0)|any(tau>=1)) stop(paste("tau must be between 0 and 1","\n",""))
  U <- runif(n)
  kappa_1 <- (nu*tau)/(1+tau*(nu-1))
  X <- sigma*(qcauchy((nu*U)/(1+U*(nu-1)))-qcauchy(kappa_1))+mu
  return(X)
} 


#Comments:
#Here we generate the RPMOC model structure for the GAMLSS family
#For more details, see chapter 6 of Stasinopoulos et al. (2017)

RPMOC<-function(mu.link = "identity", sigma.link = "log",nu.link = "log",tau.link = "identity"){
  mstats <- checklink("mu.link","RPMOC",substitute(mu.link),c("inverse","log","identity","own"))
  dstats <- checklink("sigma.link","RPMOC",substitute(sigma.link),c("inverse","log","identity","own"))
  vstats <- checklink("nu.link","RPMOC",substitute(nu.link),c("log","identity"))
  tstats <- checklink("tau.link","RPMOC",substitute(tau.link),c("identity"))
  structure(list(family = c("RPMOC","Reparameterized Marshall Olkin Cauchy"), 
                 parameters = list(mu=TRUE,sigma=TRUE,nu=TRUE,tau=FALSE),nopar=4,type = "Continuous", 
                 mu.link = as.character(substitute(mu.link)), sigma.link = as.character(substitute(sigma.link)), 
                 nu.link = as.character(substitute(nu.link)), tau.link = as.character(substitute(tau.link)), 
                 mu.linkfun = mstats$linkfun, sigma.linkfun = dstats$linkfun, 
                 nu.linkfun = vstats$linkfun, tau.linkfun = tstats$linkfun, 
                 mu.linkinv = mstats$linkinv, sigma.linkinv = dstats$linkinv, 
                 nu.linkinv = vstats$linkinv, tau.linkinv = tstats$linkinv, 
                 mu.dr = mstats$mu.eta, sigma.dr = dstats$mu.eta, nu.dr = vstats$mu.eta, tau.dr = tstats$mu.eta,
                 dldm = function(y,mu,sigma,nu,tau){
                   kappa_1 <- (nu*tau)/(1+tau*(nu-1))
                   t <- ((y-mu)/sigma)+qcauchy(kappa_1)
                   rho <- nu+(1-nu)*pcauchy(t)
                   Lambda <- 2/rho
                   zeta <- -2*t/(1+t^2) 
                   t_m <- -1/sigma
                   dldm <- t_m*(zeta-Lambda*(1-nu)*dcauchy(t))
                   return(dldm)
                 },d2ldm2 = function(y,mu,sigma,nu,tau){
                   kappa_1 <- (nu*tau)/(1+tau*(nu-1))
                   t <- ((y-mu)/sigma)+qcauchy(kappa_1)
                   rho <- nu+(1-nu)*pcauchy(t)
                   t_m <- -1/sigma
                   f_p <- (-2*t/(1+t^2))*dcauchy(t) 
                   psi <- 2*((t^2-1-t-t^3)/(1+t^2)^2)*dcauchy(t)
                   Lambda <- 2/rho
                   Lambda_m <- -0.5*Lambda^2*(1-nu)*dcauchy(t)*t_m
                   zeta <- -2*t/(1+t^2) 
                   zeta_m <- t_m*((psi/dcauchy(t))-zeta^2)  
                   d2ldm2 <- t_m*(zeta_m-(1-nu)*f_p*t_m*Lambda-(1-nu)*dcauchy(t)*Lambda_m)
                   return(d2ldm2)
                 },dldd = function(y,mu,sigma,nu,tau){
                   kappa_1 <- (nu*tau)/(1+tau*(nu-1))
                   t <- ((y-mu)/sigma)+qcauchy(kappa_1)
                   rho <- nu+(1-nu)*pcauchy(t)
                   Lambda <- 2/rho
                   zeta <- -2*t/(1+t^2) 
                   t_m <- -1/sigma
                   t_d <- -(y-mu)/sigma^2
                   dldd <- t_d*(zeta-Lambda*(1-nu)*dcauchy(t))+t_m
                   return(dldd) 
                 },d2ldd2 = function(y,mu,sigma,nu,tau){
                   kappa_1 <- (nu*tau)/(1+tau*(nu-1)) 
                   t <- ((y-mu)/sigma)+qcauchy(kappa_1) 
                   rho <- nu+(1-nu)*pcauchy(t) 
                   t_md <- 1/sigma^2 
                   t_dd <- 2*(y-mu)/sigma^3  
                   t_d <- -(y-mu)/sigma^2  
                   f_p <- (-2*t/(1+t^2))*dcauchy(t) 
                   psi <- 2*((t^2-1-t-t^3)/(1+t^2)^2)*dcauchy(t) 
                   Lambda <- 2/rho 
                   Lambda_d <- -0.5*Lambda^2*(1-nu)*dcauchy(t)*t_d 
                   zeta <- -2*t/(1+t^2)  
                   zeta_d <- t_d*((psi/dcauchy(t))-zeta^2)  
                   d2ldd2 <- t_md+t_dd*(zeta-(1-nu)*Lambda*dcauchy(t))+t_d*(zeta_d-(1-nu)*Lambda_d*dcauchy(t)-(1-nu)*Lambda*f_p*t_d)
                   return(d2ldd2)
                 },d2ldmdd = function(y,mu,sigma,nu,tau){
                   kappa_1 <- (nu*tau)/(1+tau*(nu-1)) 
                   t <- ((y-mu)/sigma)+qcauchy(kappa_1) 
                   rho <- nu+(1-nu)*pcauchy(t) 
                   t_md <- 1/sigma^2 
                   t_d <- -(y-mu)/sigma^2 
                   t_m <- -1/sigma 
                   f_p <- (-2*t/(1+t^2))*dcauchy(t) 
                   psi <- 2*((t^2-1-t-t^3)/(1+t^2)^2)*dcauchy(t)
                   Lambda <- 2/rho 
                   Lambda_d <- -0.5*Lambda^2*(1-nu)*dcauchy(t)*t_d
                   zeta <- -2*t/(1+t^2) 
                   zeta_d <- t_d*((psi/dcauchy(t))-zeta^2)  
                   d2ldmdd <- t_md*(zeta-(1-nu)*dcauchy(t)*Lambda)+t_m*(zeta_d-(1-nu)*f_p*t_d*Lambda-(1-nu)*dcauchy(t)*Lambda_d)
                   return(d2ldmdd)
                 },dldv = function(y,mu,sigma,nu,tau){
                   kappa_1 <- (nu*tau)/(1+tau*(nu-1))
                   kappa_2 <- (tau*(1-tau))/((1+tau*(nu-1))^2)
                   t <- ((y-mu)/sigma)+qcauchy(kappa_1)
                   rho <- nu+(1-nu)*pcauchy(t)
                   Lambda <- 2/rho
                   zeta <- -2*t/(1+t^2) 
                   t_n <- kappa_2/dcauchy(qcauchy(kappa_1))
                   dldv <- (1/nu)+zeta*t_n-Lambda*(1-pcauchy(t)+(1-nu)*dcauchy(t)*t_n)
                   return(dldv)
                 },dldt = function(y){
                   rep(0,length(y))
                 },d2ldv2 = function(y,mu,sigma,nu,tau){
                   kappa_1 <- (nu*tau)/(1+tau*(nu-1)) 
                   kappa_2 <- (tau*(1-tau))/((1+tau*(nu-1))^2) 
                   kappa_3 <- 2*tau^2*(1-tau)/(1+tau*(nu-1))^3 
                   t <- ((y-mu)/sigma)+qcauchy(kappa_1) 
                   rho <- nu+(1-nu)*pcauchy(t) 
                   t_n  <- kappa_2/dcauchy(qcauchy(kappa_1))  
                   f_p <- (-2*t/(1+t^2))*dcauchy(t)
                   psi <- 2*((t^2-1-t-t^3)/(1+t^2)^2)*dcauchy(t) 
                   zeta<- -2*t/(1+t^2) 
                   zeta_n <- t_n*((psi/dcauchy(t))-zeta^2) 
                   Lambda <- 2/rho 
                   Lambda_n <- -0.5*Lambda^2*(1-pcauchy(t)+(1-nu)*dcauchy(t)*t_n) 
                   t_nn <- -(-2*qcauchy(kappa_1)/(1+qcauchy(kappa_1)^2))*kappa_2^2/(dcauchy(qcauchy(kappa_1)))^2-kappa_3/dcauchy(qcauchy(kappa_1))
                   d2ldv2_1 <-  -(1/nu^2)+zeta_n*t_n+zeta*t_nn-Lambda_n*(1-pcauchy(t)+(1-nu)*dcauchy(t)*t_n)
                   d2ldv2_2 <- -Lambda*(-2*dcauchy(t)*t_n+(1-nu)*f_p*t_n^2+(1-nu)*dcauchy(t)*t_nn) 
                   d2ldv2 <- d2ldv2_1+d2ldv2_2
                   return(d2ldv2)
                 },d2ldt2 = function(y){
                   rep(0,length(y))
                 },d2ldmdv = function(y,mu,sigma,nu,tau){
                   kappa_1 <- (nu*tau)/(1+tau*(nu-1))
                   kappa_2 <- (tau*(1-tau))/((1+tau*(nu-1))^2)
                   t <- ((y-mu)/sigma)+qcauchy(kappa_1)
                   t_m <- -1/sigma 
                   t_n <- kappa_2/dcauchy(qcauchy(kappa_1))
                   f_p <- (-2*t/(1+t^2))*dcauchy(t) 
                   psi <- 2*((t^2-1-t-t^3)/(1+t^2)^2)*dcauchy(t) 
                   rho <- nu+(1-nu)*pcauchy(t)
                   Lambda <- 2/rho
                   Lambda_n <- -0.5*Lambda^2*(1-pcauchy(t)+(1-nu)*dcauchy(t)*t_n)
                   zeta <- -2*t/(1+t^2) 
                   zeta_n <- t_n*((psi/dcauchy(t))-zeta^2) 
                   d2ldmdv <- t_m*(zeta_n+dcauchy(t)*Lambda-(1-nu)*f_p*t_n*Lambda-(1-nu)*dcauchy(t)*Lambda_n)
                   return(d2ldmdv)
                 },d2ldmdt = function(y){
                   rep(0,length(y))
                 },d2ldddv = function(y,mu,sigma,nu,tau){
                   kappa_1 <- (nu*tau)/(1+tau*(nu-1)) 
                   kappa_2 <- (tau*(1-tau))/((1+tau*(nu-1))^2) 
                   t <- ((y-mu)/sigma)+qcauchy(kappa_1) 
                   rho <- nu+(1-nu)*pcauchy(t) 
                   t_d  <- -(y-mu)/sigma^2  
                   t_n  <- kappa_2/dcauchy(qcauchy(kappa_1)) 
                   f_p <- (-2*t/(1+t^2))*dcauchy(t) 
                   psi <- 2*((t^2-1-t-t^3)/(1+t^2)^2)*dcauchy(t) 
                   zeta <- -2*t/(1+t^2) 
                   zeta_n <- t_n*((psi/dcauchy(t))-zeta^2) 
                   Lambda <- 2/rho 
                   Lambda_n <- -0.5*Lambda^2*(1-pcauchy(t)+(1-nu)*dcauchy(t)*t_n) 
                   d2ldddv <- t_d*(zeta_n+Lambda*dcauchy(t)-(1-nu)*Lambda_n*dcauchy(t)-(1-nu)*Lambda*f_p*t_n) 
                   return(d2ldddv)
                 },d2ldddt = function(y){
                   rep(0, length(y))
                 }, d2ldvdt = function(y){
                   rep(0, length(y))
                 },G.dev.incr = function(y,mu,sigma,nu,tau,...) -2*dRPMOC(y,mu,sigma,nu,tau,log = TRUE),
                 rqres = expression(rqres(pfun = "pRPMOC",type = "Continuous",y=y,mu=mu,sigma=sigma,nu=nu,tau=tau)), 
                 mu.initial = expression(mu <-rep(median(y),length(y))), sigma.initial = expression(sigma <- rep((sd(y)+0.001),length(y))), 
                 nu.initial = expression(nu<-rep(1,length(y))),tau.initial = expression(tau<-rep(0.1,length(y))), 
                 mu.valid = function(mu) TRUE, sigma.valid = function(sigma) all(sigma>0), nu.valid = function(nu) all(nu>0),
                 tau.valid = function(tau) all((tau>0)&(tau<1)),y.valid = function(y) TRUE),class = c("gamlss.family","family"))}

