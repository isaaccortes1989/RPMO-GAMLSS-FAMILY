#===Comments:
#===pRPMOT: cdf for the RPMOT model with df=15
#===dRPMOT: pdf for the RPMOT model df=15
#===qRPMOT: quantile function for the RPMOT model df=15
#===rRPMOT: function to generate random variables that follow a RPMOT distribution df=15  
#===mu is the quantile parameter 
#===sigma is the scale parameter 
#===nu is the skewness parameter
#===tau is the parameter of the fixed probability

rm(list=ls())
#-------------------------------------------------------------------------------
pRPMOT = function(q,mu=1,sigma=1,nu=1,tau=0.5,lower.tail=TRUE,log.p=FALSE){
  if(any(sigma<=0)) stop(paste("sigma must be positive","\n",""))
  if(any(nu<=0)) stop(paste("nu must be positive","\n",""))
  if(any(tau<=0)|any(tau>=1)) stop(paste("tau must be between 0 and 1","\n",""))
  df      <- 15
  kappa_1 <- ((nu*tau)/(1+tau*(nu-1)))
  t       <- ((q-mu)/sigma)+qt(kappa_1, df = df)
  lGy     <- log(pt(t,df = df))-log(((1-nu)*pt(t,df = df)+nu)) 
  if(log.p == FALSE)
    cdf    <- exp(lGy)  
  else cdf <- lGy
  if(lower.tail == TRUE)
    cdf <- cdf
  else cdf <- 1-cdf
  return(cdf)
}

#---------------------------------------------------------
dRPMOT = function(x,mu=1,sigma=1,nu=1,tau=0.5,log=FALSE){
  if(any(sigma<=0)) stop(paste("sigma must be positive","\n",""))
  if(any(nu<=0))    stop(paste("nu must be positive","\n",""))
  if(any(tau<=0)|any(tau>=1)) stop(paste("tau must be between 0 and 1","\n",""))
  df      <- 15
  kappa_1 <- (nu*tau)/(1+tau*(nu-1))
  t       <- ((x-mu)/sigma)+qt(kappa_1,df=df)
  lfy     <- log(nu)+dt(t,df=df,log=TRUE)-log(sigma)-2*log(nu+(1-nu)*pt(t,df=df))
  if(log== FALSE)
    fy <- exp(lfy)
  else fy <- lfy
  return(fy)
}#BIEN
#-------------------------------------------------------------------------------
qRPMOT = function(p,mu = 1,sigma = 1,nu = 1,tau = 0.5,lower.tail = TRUE,log.p = FALSE){ 
  if(any(sigma<=0)) stop(paste("sigma must be positive","\n",""))
  if(any(nu<=0)) stop(paste("nu must be positive","\n",""))
  if(any(tau<=0)|any(tau>=1)) stop(paste("tau must be between 0 and 1","\n",""))
  df <- 15
  if(log.p == TRUE)
    p  <- exp(p)  
  else p <- p
  if(lower.tail == TRUE)
    p  <- p
  else p <- 1-p
  if(any(p<0)|any(p>1)) stop(paste("p must be between 0 and 1","\n",""))
  kappa_1 <- (nu*tau)/(1+tau*(nu-1))
  Q       <- sigma*(qt((nu*p)/(1+p*(nu-1)),df = df)-qt(kappa_1,df = df))+mu
  return(Q)
}
#-------------------------------------------------------------------------------
rRPMOT = function(n,mu=1,sigma=1,nu=1,tau=0.5){
  if(any(sigma<=0)) stop(paste("sigma must be positive","\n",""))
  if(any(nu<=0))    stop(paste("nu must be positive","\n",""))
  if(any(tau<=0)|any(tau>=1)) stop(paste("tau must be between 0 and 1","\n",""))
  df <- 15
  U <- runif(n)
  kappa_1 <- (nu*tau)/(1+tau*(nu-1))
  X       <- sigma*(qt((nu*U)/(1+U*(nu-1)),df=df)-qt(kappa_1,df=df))+mu
  return(X)
}

#===Comments:
#===Here we generate the RPMOT model structure for the GAMLSS family
#===For more details, see chapter 6 of Stasinopoulos et al. (2017)
#-------------------------------------------------------------------------------
RPMOT<-function(mu.link = "identity", sigma.link = "log",nu.link = "log",tau.link = "identity"){
  mstats <- checklink("mu.link","RPMOT",substitute(mu.link),c("inverse","log","identity","own"))
  dstats <- checklink("sigma.link","RPMOT",substitute(sigma.link),c("inverse","log","identity","own"))
  vstats <- checklink("nu.link","RPMOT",substitute(nu.link),c("log","identity"))
  tstats <- checklink("tau.link","RPMOT",substitute(tau.link),c("identity"))
  structure(list(family = c("RPMOT","Reparameterized Marshall Olkin T"), 
                 parameters = list(mu=TRUE,sigma=TRUE,nu=TRUE,tau=FALSE),nopar=4,type = "Continuous", 
                 mu.link = as.character(substitute(mu.link)), sigma.link = as.character(substitute(sigma.link)), 
                 nu.link = as.character(substitute(nu.link)), tau.link = as.character(substitute(tau.link)), 
                 mu.linkfun = mstats$linkfun, sigma.linkfun = dstats$linkfun, 
                 nu.linkfun = vstats$linkfun, tau.linkfun = tstats$linkfun, 
                 mu.linkinv = mstats$linkinv, sigma.linkinv = dstats$linkinv, 
                 nu.linkinv = vstats$linkinv, tau.linkinv = tstats$linkinv, 
                 mu.dr = mstats$mu.eta, sigma.dr = dstats$mu.eta, nu.dr = vstats$mu.eta, tau.dr = tstats$mu.eta,
                 dldm = function(y,mu,sigma,nu,tau){
                   df      <- 15
                   kappa_1 <- (nu*tau)/(1+tau*(nu-1))
                   t       <- ((y-mu)/sigma)+qt(kappa_1, df = df)
                   rho     <- nu+(1-nu)*pt(t, df = df)
                   Lambda  <- 2/rho
                   zeta    <- -((df+1)/df)*(t/(1+t^2/df))
                   t_m     <- -1/sigma
                   dldm    <- t_m*(zeta-Lambda*(1-nu)*dt(t,df = df))
                   return(dldm)
                 },d2ldm2 = function(y,mu,sigma,nu,tau){
                   df       <- 15
                   kappa_1  <- (nu*tau)/(1+tau*(nu-1))
                   t        <- ((y-mu)/sigma)+qt(kappa_1,df = df)
                   rho      <- nu+(1-nu)*pt(t,df = df)
                   t_m      <- -1/sigma
                   f_p      <- -((df+1)/df)*(t/(1+t^2/df))*dt(t,df=df)  
                   psi      <- dt(t,df = df)*(((df+1)/df)^2*t^2/(1+t^2/df)^2-(df+1)*(df-t^2)/(df+t^2)^2)
                   Lambda   <- 2/rho
                   Lambda_m <- -0.5*Lambda^2*(1-nu)*dt(t,df=df)*t_m 
                   zeta     <- -((df+1)/df)*(t/(1+t^2/df)) 
                   zeta_m   <- t_m*((psi/dt(t,df=df))-zeta^2)   
                   d2ldm2   <- t_m*(zeta_m-(1-nu)*f_p*t_m*Lambda-(1-nu)*dt(t,df = df)*Lambda_m)
                   return(d2ldm2)
                 },dldd = function(y,mu,sigma,nu,tau){
                   df      <- 15
                   kappa_1 <- (nu*tau)/(1+tau*(nu-1))
                   t       <- ((y-mu)/sigma)+qt(kappa_1,df = df)
                   rho     <- nu+(1-nu)*pt(t,df = df)
                   Lambda  <- 2/rho
                   zeta    <- -((df+1)/df)*(t/(1+t^2/df)) 
                   t_m     <- -1/sigma
                   t_d     <- -(y-mu)/sigma^2
                   dldd    <- t_d*(zeta-Lambda*(1-nu)*dt(t,df = df))+t_m
                   return(dldd) 
                 },d2ldd2 = function(y,mu,sigma,nu,tau){
                   df       <- 15
                   kappa_1  <- (nu*tau)/(1+tau*(nu-1)) 
                   t        <- ((y-mu)/sigma)+qt(kappa_1,df=df) 
                   rho      <- nu+(1-nu)*pt(t,df=df) 
                   t_md     <- 1/sigma^2 
                   t_dd     <- 2*(y-mu)/sigma^3  
                   t_d      <- -(y-mu)/sigma^2 
                   f_p      <- -((df+1)/df)*(t/(1+t^2/df))*dt(t,df=df) 
                   psi      <- dt(t,df=df)*(((df+1)/df)^2*t^2/(1+t^2/df)^2-(df+1)*(df-t^2)/(df+t^2)^2)
                   Lambda   <- 2/rho
                   Lambda_d <- -0.5*Lambda^2*(1-nu)*dt(t,df=df)*t_d 
                   zeta     <- -((df+1)/df)*(t/(1+t^2/df)) 
                   zeta_d   <- t_d*((psi/dt(t,df=df))-zeta^2) 
                   d2ldd2   <- t_md+t_dd*(zeta-(1-nu)*Lambda*dt(t,df=df))+t_d*(zeta_d-(1-nu)*Lambda_d*dt(t,df=df)-(1-nu)*Lambda*f_p*t_d)
                   return(d2ldd2)
                 },d2ldmdd = function(y,mu,sigma,nu,tau){
                   df       <- 15
                   kappa_1  <- (nu*tau)/(1+tau*(nu-1)) 
                   t        <- ((y-mu)/sigma)+qt(kappa_1,df=df) 
                   rho      <- nu+(1-nu)*pt(t,df=df) 
                   t_md     <- 1/sigma^2 
                   t_d      <- -(y-mu)/sigma^2 
                   t_m      <- -1/sigma 
                   f_p      <- -((df+1)/df)*(t/(1+t^2/df))*dt(t,df=df) 
                   psi      <- dt(t,df=df)*(((df+1)/df)^2*t^2/(1+t^2/df)^2-(df+1)*(df-t^2)/(df+t^2)^2)
                   Lambda   <- 2/rho 
                   Lambda_d <- -0.5*Lambda^2*(1-nu)*dt(t,df=df)*t_d 
                   zeta     <- -((df+1)/df)*(t/(1+t^2/df)) 
                   zeta_d   <- t_d*((psi/dt(t,df=df))-zeta^2) 
                   d2ldmdd  <- t_md*(zeta-(1-nu)*dt(t,df=df)*Lambda)+t_m*(zeta_d-(1-nu)*f_p*t_d*Lambda-(1-nu)*dt(t,df=df)*Lambda_d)
                   return(d2ldmdd)
                 },dldv = function(y,mu,sigma,nu,tau){
                   df      <- 15
                   kappa_1 <- (nu*tau)/(1+tau*(nu-1))
                   kappa_2 <- (tau*(1-tau))/((1+tau*(nu-1))^2)
                   t       <- ((y-mu)/sigma)+qt(kappa_1,df = df)
                   rho     <- nu+(1-nu)*pt(t,df = df)
                   Lambda  <- 2/rho
                   zeta    <- -((df+1)/df)*(t/(1+t^2/df)) 
                   t_n     <- kappa_2/dt(qt(kappa_1,df = df),df = df)
                   dldv    <- (1/nu)+zeta*t_n-Lambda*(1-pt(t,df = df)+(1-nu)*dt(t,df = df)*t_n)
                   return(dldv)
                 },dldt = function(y){
                   rep(0,length(y))
                 },d2ldv2 = function(y,mu,sigma,nu,tau){
                   df       <- 15
                   kappa_1  <- (nu*tau)/(1+tau*(nu-1)) 
                   kappa_2  <- (tau*(1-tau))/((1+tau*(nu-1))^2)
                   kappa_3  <- 2*tau^2*(1-tau)/(1+tau*(nu-1))^3 
                   t        <- ((y-mu)/sigma)+qt(kappa_1,df=df) 
                   rho      <- nu+(1-nu)*pt(t,df=df) 
                   t_n      <- kappa_2/dt(qt(kappa_1,df=df),df=df)
                   f_p      <- -((df+1)/df)*(t/(1+t^2/df))*dt(t,df=df) 
                   psi      <- dt(t,df=df)*(((df+1)/df)^2*t^2/(1+t^2/df)^2-(df+1)*(df-t^2)/(df+t^2)^2)
                   zeta     <- -((df+1)/df)*(t/(1+t^2/df)) 
                   zeta_n   <- t_n*((psi/dt(t,df=df))-zeta^2) 
                   Lambda   <- 2/rho 
                   Lambda_n <- -0.5*Lambda^2*(1-pt(t,df=df)+(1-nu)*dt(t,df=df)*t_n)
                   t_m      <- -1/sigma 
                   t_nn     <- ((df+1)/df)*(qt(kappa_1,df=df)/(1+qt(kappa_1,df=df)^2/df))*kappa_2^2/(dt(qt(kappa_1,df=df),df=df))^2 -kappa_3/dt(qt(kappa_1,df=df),df=df)
                   d2ldv2_1 <- -(1/nu^2)+zeta_n*t_n+zeta*t_nn-Lambda_n*(1-pt(t,df=df)+(1-nu)*dt(t,df=df)*t_n) 
                   d2ldv2_2 <- -Lambda*(-2*dt(t,df=df)*t_n+(1-nu)*f_p*t_n^2+(1-nu)*dt(t,df=df)*t_nn)
                   d2ldv2   <- d2ldv2_1+d2ldv2_2
                   return(d2ldv2)
                 },d2ldt2 = function(y){
                   rep(0,length(y))
                 },d2ldmdv = function(y,mu,sigma,nu,tau){
                   df       <- 15
                   kappa_1  <- (nu*tau)/(1+tau*(nu-1)) 
                   kappa_2  <- (tau*(1-tau))/((1+tau*(nu-1))^2)
                   t        <- ((y-mu)/sigma)+qt(kappa_1,df=df) 
                   t_m      <- -1/sigma 
                   t_n      <- kappa_2/dt(qt(kappa_1,df=df),df=df)
                   f_p      <- -((df+1)/df)*(t/(1+t^2/df))*dt(t,df=df) 
                   psi      <- dt(t,df=df)*(((df+1)/df)^2*t^2/(1+t^2/df)^2-(df+1)*(df-t^2)/(df+t^2)^2)
                   rho      <- nu+(1-nu)*pt(t,df=df) 
                   Lambda   <- 2/rho 
                   Lambda_n <- -0.5*Lambda^2*(1-pt(t,df=df)+(1-nu)*dt(t,df=df)*t_n)
                   zeta     <- -((df+1)/df)*(t/(1+t^2/df)) 
                   zeta_n   <- t_n*((psi/dt(t,df=df))-zeta^2) 
                   d2ldmdv  <- t_m*(zeta_n+dt(t,df=df)*Lambda-(1-nu)*f_p*t_n*Lambda-(1-nu)*dt(t,df=df)*Lambda_n)
                   return(d2ldmdv)
                 },d2ldmdt = function(y){
                   rep(0,length(y))
                 },d2ldddv = function(y,mu,sigma,nu,tau){
                   df       <- 15
                   kappa_1  <- (nu*tau)/(1+tau*(nu-1)) 
                   kappa_2  <- (tau*(1-tau))/((1+tau*(nu-1))^2)
                   t        <- ((y-mu)/sigma)+qt(kappa_1,df=df) 
                   rho      <- nu+(1-nu)*pt(t,df=df) 
                   t_d      <- -(y-mu)/sigma^2 
                   t_n      <- kappa_2/dt(qt(kappa_1,df=df),df=df)
                   f_p      <- -((df+1)/df)*(t/(1+t^2/df))*dt(t,df=df) 
                   psi      <- dt(t,df=df)*(((df+1)/df)^2*t^2/(1+t^2/df)^2-(df+1)*(df-t^2)/(df+t^2)^2)
                   zeta     <- -((df+1)/df)*(t/(1+t^2/df)) 
                   zeta_n   <- t_n*((psi/dt(t,df=df))-zeta^2) 
                   Lambda   <- 2/rho 
                   Lambda_n <- -0.5*Lambda^2*(1-pt(t,df=df)+(1-nu)*dt(t,df=df)*t_n)
                   d2ldddv  <- t_d*(zeta_n+Lambda*dt(t,df=df)-(1-nu)*Lambda_n*dt(t,df=df)-(1-nu)*Lambda*f_p*t_n) 
                   return(d2ldddv)
                 },d2ldddt = function(y){
                   rep(0, length(y))
                 },d2ldvdt = function(y){
                   rep(0, length(y))
                 },G.dev.incr = function(y,mu,sigma,nu,tau,...) -2*dRPMOT(y,mu,sigma,nu,tau,log = TRUE),
                 rqres = expression(rqres(pfun = "pRPMOT",type = "Continuous",y=y,mu=mu,sigma=sigma,nu=nu,tau=tau)), 
                 mu.initial = expression(mu <-rep(median(y),length(y))), sigma.initial = expression(sigma <- rep((sd(y)+0.001),length(y))), 
                 nu.initial = expression(nu<-rep(1,length(y))),tau.initial = expression(tau<-rep(0.1,length(y))), 
                 mu.valid = function(mu) TRUE, sigma.valid = function(sigma) all(sigma>0), nu.valid = function(nu) all(nu>0),
                 tau.valid = function(tau) all((tau>0)&(tau<1)),y.valid = function(y) TRUE),class = c("gamlss.family","family"))}

