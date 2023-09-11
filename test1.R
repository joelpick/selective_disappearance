library(lme4)
n=1000
AFR <- rpois(n,3)+1
repro_lifespans <- rpois(n,3)+1#rgeom(n,0.3)+1
LR <- AFR + repro_lifespans -1


mean_age <- (LR + AFR )/ 2
beta_l <- -1
beta_f <- 2
scaled_mean_age <-(beta_l*LR + beta_f*AFR )/ 2

beta<-(beta_l*var(LR)/var(LR + AFR) + beta_f*var(AFR)/var(LR + AFR)  + 2*beta_l*beta_f*cov(AFR,LR)/var(LR + AFR))

var(mean_age)
var(scaled_mean_age)
cov(mean_age,scaled_mean_age)

1/4*(beta_l*(cov(AFR,LR)+var(LR)) + beta_f*(cov(AFR,LR)+var(AFR)))

(beta_l*(cov(AFR,LR)+var(LR)) + beta_f*(cov(AFR,LR)+var(AFR))) / (4*(var(LR) + var(AFR) + 2*cov(AFR,LR)))

beta <- sqrt(var(scaled_mean_age)/var(mean_age))/2
plot(mean_age*beta,scaled_mean_age);abline(0,1)
  
var(LR)
var(AFR)
cov(AFR,LR)

(beta_l^2*var(LR) + beta_f^2*var(AFR) + 2*beta_l*beta_f*cov(AFR,LR))/4

(var(LR) + var(AFR) + 2*cov(AFR,LR))/4

(beta_l^2*var(LR) + beta_f^2*var(AFR) + 2*beta_l*beta_f*cov(AFR,LR))/(var(LR) + var(AFR) + 2*cov(AFR,LR))
beta_l^2*var(LR)/var(LR + AFR) + beta_f^2*var(AFR)/var(LR + AFR)  + 2*beta_l*beta_f*cov(AFR,LR)/var(LR + AFR)

## maybe its to do with the covariance between mean age and (beta_l*l+beta_f*f)/2

beta_a <- -0.1
beta_l <- 0.2
beta_f <- 0.5

dat<-data.frame(
	id=rep(1:n,repro_lifespans),
	age=c(sapply(1:n, function(x) AFR[x]:LR[x]),recursive=TRUE),
	FR = rep(AFR,repro_lifespans),
	LR = rep(LR,repro_lifespans)
)
dat$y <- dat$age*-0.1 + dat$LR*beta_l + dat$FR*beta_f + rnorm(nrow(dat),0,0.5)
age_bar <- tapply(dat$age,dat$id,mean)
dat$age_bar <- age_bar[match(dat$id,names(age_bar))]
dat$age_dev <- dat$age - dat$age_bar

mod_func <- function(formula) c(coef(lm(formula,dat)),Vr=summary(lm(formula,dat))$sigma^2)
list(
	mod_a = mod_func(y~age),
	mod_l = mod_func(y~age+LR),
	mod_f = mod_func(y~age+FR),
	mod_fl = mod_func(y~age+FR+LR),
	mod_c = mod_func(y~age_dev+age_bar),
	mod_c2 = mod_func(y~age+age_bar)
	)

	mod_fl = mod_func(y~age+FR+LR)

E_beta_ab <-(beta_l*(cov(dat$FR,dat$LR)+var(dat$LR)) + beta_f*(cov(dat$FR,dat$LR)+var(dat$FR))) / (var(dat$LR) + var(dat$FR) + 2*cov(dat$FR,dat$LR))

2* (mod_fl["LR"]*(cov(dat$FR,dat$LR)+var(dat$LR)) + mod_fl["FR"]*(cov(dat$FR,dat$LR)+var(dat$FR))) / (var(dat$LR) + var(dat$FR) + 2*cov(dat$FR,dat$LR))


# var(dat$age_bar)

mod_c[1] + (mod_c[3]-mod_c[2])/2 
(mod_c[3]-mod_c[2])/2 



var(dat$y)
var(dat$LR)
var(dat$FR)
sum(cov(dat[,c("FR","LR")]))

## AFR and LR likely to be correlated? because LR is AFR plus repro lifespan?
 
sqrt(
	(
		(mod_fl["LR"]-mod_fl["age"])^2*var2(dat$LR) + 
		(mod_fl["FR"]-mod_fl["age"])^2*var2(dat$FR) + 
		2*(mod_fl["FR"]-mod_fl["age"])*(mod_fl["LR"]-mod_fl["age"])*cov2(dat$FR,dat$LR)
	)/(
		sum(cov3(dat[,c("FR","LR")]))
	))

sqrt(
	(
		mod_fl["LR"]^2*var2(dat$LR) + 
		mod_fl["FR"]^2*var2(dat$FR) + 
		2*mod_fl["FR"]*mod_fl["LR"]*cov2(dat$FR,dat$LR)
	)/(
		sum(cov3(dat[,c("FR","LR")]))
	))*2 + 



sqrt((mod_fl["LR"]^2*var(dat$LR) + mod_fl["FR"]^2*var(dat$FR) + 2*mod_fl["FR"]*mod_fl["LR"]*cov(dat$FR,dat$LR))/sum(cov(dat[,c("FR","LR")])))

var2 <- function(y) mean((y-mean(y))^2)
cov2 <- function(x,y) mean((x-mean(x))*(y-mean(y)))
cov3 <- function(x) cov(x)*(n-1)/n
