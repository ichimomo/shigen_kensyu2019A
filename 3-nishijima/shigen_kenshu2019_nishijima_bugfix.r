
# set working directory

# read library
library(tidyverse)
library(ggplot2)

### GLMM ###

# data simulation
N = 100
beta = c(1,0.5,1) # parameter setting
X = seq(-1,1,length=N) # explanatory variable
set.seed(12345)
epsilon = rnorm(n=N,mean=0,sd=beta[3]) # random variable
r = beta[1] + beta[2]*X + epsilon # linear predictor
Y = rpois(n=N, lambda = exp(r)) # response variable
id = factor(1:N)

plot(Y~X,pch=16,cex=1.5)

poisson_model = glm(Y~X, family=poisson(link="log"))
summary(poisson_model)


library(lme4)
poisson_model_RE = glmer(Y~X + (1|id), family=poisson(link="log"))
summary(poisson_model_RE)

AIC(poisson_model,poisson_model_RE)

coef(poisson_model)
coef(poisson_model_RE)
plot(coef(poisson_model_RE)$id[,1]-summary(poisson_model_RE)$coefficients[1,1]~epsilon,
     xlab="True",ylab="Estimated",main="Random effect")
abline(0,1)

# use nlminb (optim)
poisson_obj = function(beta) {
  r = beta[1] + beta[2]*X
  negative_loglik = -sum(dpois(Y, lambda = exp(r), log = TRUE))
  return(negative_loglik)
}

poisson_opt = nlminb(beta[1:2],poisson_obj)
poisson_opt$par #parameter
poisson_opt$objective #negative loglik

poissonRE_obj = function(beta) {
  negative_loglik = 0  
  for (i in 1:length(X)) {
    likelihood = function(epsilon) {
      r = beta[1] + beta[2]*X[i] + epsilon
      dpois(Y[i],exp(r))*dnorm(epsilon,0,beta[3])
    }
    res = integrate(likelihood,-Inf,Inf,abs.tol = 0L)
    negative_loglik = negative_loglik-log(res$value)
    }
  return(negative_loglik)
}

poissonRE_opt = nlminb(beta, poissonRE_obj)
poissonRE_opt$par # parameters
poissonRE_opt$objective #negative log-likelihood

## use TMB
library(TMB)
compile("glmm_poisson.cpp")
dyn.load(dynlib("glmm_poisson"))

data = list(Y=Y,X=X)
params = list(beta=beta,epsilon=rep(0,N))
obj = MakeADFun(data, params, random="epsilon",DLL="glmm_poisson")
opt = nlminb(obj$par, obj$fn, obj$gr)
sdrep = sdreport(obj) #random effect
obj$env$parList()$beta
obj$env$parList()$epsilon


### surplus production model using TMB ###

### surplus production model using TMB ###

compile("surplus_production.cpp")
dyn.load(dynlib("surplus_production"))

# Schaefer 
nyear = 50
r = 0.5
k = 100
(MSY = r*k/4)
sigma_pro = 0.4
sigma_obs = 0.2
B0 = 50
F0 = 0.3
q = 1
set.seed(124)

F_vec <- B_vec <- C_vec <- NULL
for(i in 1:nyear) {
  F_y = exp(log(F0)+rnorm(1,0,0.1))
  F_vec = c(F_vec,F_y)
  F0=F_y
  C_y = B0*(1-exp(-F_y))
  C_vec = c(C_vec,C_y)
  B_y = (B0+r*B0*(1-B0/k)-C_y)*exp(rnorm(1,-sigma_pro^2/2,sigma_pro))
  B_vec = c(B_vec,B_y)
  B0 = B_y
}

Catch = C_vec
cpue = exp(rnorm(nyear,log(q*B_vec),sigma_obs))

data = list(Catch=Catch,cpue=cpue,SP_type=0)
params = list(log_r=log(r),log_K=log(k),log_sigma_pro=log(sigma_pro),
              log_q=log(q),log_sigma_obs=log(sigma_obs),log_B=rep(log(0.5*k),nyear))
obj = MakeADFun(data, params, random="log_B",DLL="surplus_production")
opt = nlminb(obj$par, obj$fn, obj$gr)
sdrep = sdreport(obj,bias.correct=TRUE) #random effect
(r_est = exp(obj$env$parList()$log_r))
(K_est =exp(obj$env$parList()$log_K))
r_est*K_est/4  # MSY estimate

exp(obj$env$parList()$log_sigma_pro)
exp(obj$env$parList()$log_sigma_obs)
B_est = sdrep$unbiased$value

matplot(cbind(B_vec,B_est),type="l",lwd=3, main = "Biomass",ylab="Biomass")
matplot(cbind(B_vec/k,B_est/K_est),type="l",lwd=3,main = "Depletion rate",ylab="B/K")

# Fox model

data2 = list(Catch=Catch,cpue=cpue,SP_type=1)
# params = list(log_r=log(r),log_K=log(k),log_sigma_pro=log(sigma_pro),
#               log_q=log(q),log_sigma_obs=log(sigma_obs),log_B=rep(log(0.5*k),nyear))
obj2 = MakeADFun(data2, params, random="log_B",DLL="surplus_production")
opt2 = nlminb(obj2$par, obj2$fn, obj2$gr)
sdrep2 = sdreport(obj2,bias.correct = TRUE) #random effect

(r_est2 = exp(obj2$env$parList()$log_r))
(K_est2 =exp(obj2$env$parList()$log_K))
r_est2*K_est2/(exp(1)*log(K_est2)) # MSY estimates

B_est2 = sdrep2$unbiased$value

matplot(cbind(B_vec,B_est2),type="l",lwd=2)
matplot(cbind(B_vec/k,B_est2/K_est2),type="l",lwd=2)


### delay-difference model ### 

sim_dd = function(F,l0,beta,sigma,x0=0.5,nyear=200,nsim=100,seed=1,rec_dev=NULL) {
  argname <- ls() 
  arglist <- lapply(argname,function(xx) eval(parse(text=xx)))
  names(arglist) <- argname
  
  set.seed(seed)
  if(is.null(rec_dev)){
    rec_dev = matrix(rnorm(nsim*nyear,-sigma^2/2,sigma),nrow=nsim)
  }
  M = -log(l0)
  l_f = exp(-F)
  x0 = rep_len(x0,nsim)
  r<-X<-Catch<-c()
  for(t in 1:nyear){
    r_t = ifelse(x0<beta/l_f,l_f*x0/beta,1)*(1-l0)*exp(rec_dev[,t])
    r = cbind(r,r_t)
    X_t = r_t+x0*exp(-M-F)
    X = cbind(X,X_t)
    Catch = cbind(Catch,X_t*(1-exp(-M-F))*F/(F+M))
    x0 = X_t
  }
  list(input=arglist,r=r,X=X,Catch=Catch,M=M,rec_dev=rec_dev)
}

# sim_dd_res = sim_dd(F=0.2,l0=0.6,beta=0.3,sigma=0,nyear=200,nsim=1000)

est_MSY_dd = function(sim_dd_res,interval=c(-10,10)) {
  input_tmp = sim_dd_res$input
  input_tmp$rec_dev = sim_dd_res$rec_dev
  obj_f = function(x) {
    input_tmp$F = exp(x)
    sim_dd_test = do.call(sim_dd, input_tmp)
    -log(mean(sim_dd_test$Catch[,ncol(sim_dd_test$Catch)]))
  }
  opt=optimize(obj_f,interval)
  msy=exp(-opt$objective)
  Fmsy=exp(opt$minimum)
  input_tmp$F=Fmsy
  sim_dd_msy=do.call(sim_dd,input_tmp)
  Xmsy=mean(sim_dd_msy$X[,ncol(sim_dd_msy$X)])
  Smsy=Xmsy*exp(-Fmsy)
  return(c(msy=msy,Fmsy=Fmsy,Xmsy=Xmsy,Smsy=Smsy))
}


sim_dd_res = sim_dd(F=0.2,l0=0.6,beta=0.9,sigma=0,nyear=200,nsim=1000)
test_msy = est_MSY_dd(sim_dd_res)
test_msy

sim_dd_res = sim_dd(F=0.2,l0=0.1,beta=0.3,sigma=0.9,nyear=200,nsim=1000)
test_msy = est_MSY_dd(sim_dd_res)
test_msy

par_grid = expand.grid(l0=seq(0.1,0.9,by=0.4),
                       beta=seq(0.05,0.95,by=0.05),
                       sigma=seq(0,1.2,by=0.4))

msy_res = c()
for (i in 1:nrow(par_grid)) {
  sim_dd_res = sim_dd(F=0.2,l0=par_grid[i,1],beta=par_grid[i,2],sigma=par_grid[i,3],nyear=200,nsim=1000)
  test_msy = est_MSY_dd(sim_dd_res)
  msy_res = rbind(msy_res,test_msy)
  cat(i,"\n")
}

msy_res = cbind(par_grid,msy_res) %>% as_tibble()
msy_res = mutate(msy_res,id=1:nrow(msy_res))

write.csv(msy_res, file="msy_res.csv",row.names = FALSE)

msy_res_gg = filter(msy_res,l0 %in% c(0.1,0.5,0.9)) %>%
  filter(beta<0.95) %>% 
  mutate(sigma_rec = factor(sigma))

(g1 = ggplot(msy_res_gg,aes(x=beta,y=Fmsy,group=sigma_rec,colour=sigma_rec))+
    theme_bw(base_size=16)+
    geom_path(size=2)+
    facet_wrap(~l0,ncol=3)+
    # ylim(0,NA)+
    scale_y_continuous(expand=c(0,0),limits=c(0,max(msy_res_gg$Fmsy)*1.05))
)

(g2 = ggplot(msy_res_gg,aes(x=beta,y=msy,group=sigma_rec,colour=sigma_rec))+
    theme_bw(base_size=16)+
    geom_path(size=2)+
    facet_wrap(~l0,ncol=3)+
    # ylim(0,NA)+
    scale_y_continuous(expand=c(0,0),limits=c(0,max(msy_res_gg$msy)*1.05))
)

(g3 = ggplot(msy_res_gg,aes(x=beta,y=Smsy,group=sigma_rec,colour=sigma_rec))+
    theme_bw(base_size=16)+
    geom_path(size=2)+
    facet_wrap(~l0,ncol=3)+
    # ylim(0,NA)+
    scale_y_continuous(expand=c(0,0),limits=c(0,max(msy_res_gg$Smsy)*1.05))
  
)

ggsave("delay-difference_Fmsy_g1.png",g1, width=200, height=120, dpi=600, units="mm")
ggsave("delay-difference_msy_g1.png",g2, width=200, height=120, dpi=600, units="mm")
ggsave("delay-difference_Smsy_g1.png",g3, width=200, height=120, dpi=600, units="mm")


### VPA and SAM ### 

source("rvpa1.9.4.r")
compile("rvpa_tmb4.cpp")
dyn.load(dynlib("rvpa_tmb4"))

source("rsam0.5.r")
compile("jsam9.cpp")
dyn.load(dynlib("jsam9"))

### data handling
caa <- read.csv("caa0_2018.csv",row.names=1)
waa <- read.csv("waa1_2018.csv",row.names=1) #2018年の6+を変更
maa <- read.csv("maa1_2018.csv",row.names=1)
cpue <- read.csv("cpue4_bic2018.csv",row.names=1)

dat <- data.handler(caa, waa, maa, cpue, M=0.4)

vpa_base <- vpa(dat,
                tune = TRUE,
                term.F = "all",
                alpha = 1,
                abund = c("N","N","SSB","SSB"),
                min.age = c(0,0,0,0),
                max.age = c(0,0,6,6),
                est.method="ml", #重み推定（最尤法）
                b.est = TRUE, #b推定
                sel.def="max",
                b.fix = c(NA,NA,1,1), #北上期と秋季はb推定し、たもすくいと産卵量はb=1にfix
                lambda = 0.6754305, 
                last.catch.zero=TRUE,
                rec.new.index = 1:2, #最新年の加入量を予測する際に使うindex
                # rec.new=24302,
                fc.year=2015:2017,
                plot = FALSE,
                TMB=TRUE,
                p.init=0.2,
                use.index =1:4,
                sigma.constraint = c(1,1,2,3) #北上期と秋季のsigmaを一定
)

sam_base <- sam(dat, #rvpaと同じデータ
                last.catch.zero = TRUE,
                abund = c("N","N","SSB","SSB"),
                cpp.file.name = "jsam9",
                index.age=c(0,0,0,0),
                b.est=TRUE,
                b.fix=c(NA,NA,1,1),
                SR = "BH",
                varC = c(0,1,2,2,3,4,4),
                varF = c(0,0,1,1,1,1,1),
                varN = c(0,1,1,1,1,1,1),
                varN.fix=c(NA,1e-4), 
                rho.mode=3,
                bias.correct = TRUE,
                get.random.vcov = FALSE
)

rbind(colSums(vpa_base$baa),colSums(sam_base$baa)) # biomass
rbind(colSums(vpa_base$ssb),colSums(sam_base$ssb)) # SSB
