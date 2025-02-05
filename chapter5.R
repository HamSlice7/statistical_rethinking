##5.1
#load data and copy
library(rethinking)
data("WaffleDivorce")
d <- WaffleDivorce

#standardize variables - zero centered, sd 1
d$D <- standardize(d$Divorce)
d$M <- standardize(d$Marriage)
d$A <- standardize(d$MedianAgeMarriage)

plot(D ~ M, d)
plot(D ~ A, d)

##5.2
sd(d$MedianAgeMarriage)

##5.3 - approximate posterior 
m5.1 <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + bA * A,
    a ~ dnorm(0, 0.2),
    bA ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data = d
)

##5.4 simulate from the priors
set.seed(10)
prior <- extract.prior(m5.1) #extracting random values from prior distributions
mu <- link(m5.1, post = prior, data = list(A=c(-2,2))) #finding mean D using randomly generated parameters from prior distribution
plot( NULL , xlim=c(-2,2) , ylim=c(-2,2) )
for (i in 1:50) {
  lines(c(-2,2), mu[i,], col = col.alpha("black", 0.4))
}

##5.5 Posterior Predictions
#compute percentile interval of mean
A_seq <- seq(from = -3, to = 3.2, length.out=30)
mu <- link(m5.1, data = list(A=A_seq)) #extracting mean D from fitted parameters
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

#plot it all
plot(D ~A, data = d, col=rangi2)
lines(A_seq, mu.mean, lwd=2)
shade(mu.PI, A_seq)
precis(m5.1)

##5.6
#Fitting model for marriage rate and divorce rate
m5.2 <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + bM * M, 
    a ~ dnorm(0,0.2),
    bM ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data = d
)

#compute percentile interval of mean
M_seq <- seq(from = -3, to = 3.2, length.out=30)
mu <- link(m5.2, data = list(M=M_seq)) #extracting mean D from fitted parameters
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

#plot it all
plot(D ~M, data = d, col=rangi2)
lines(M_seq, mu.mean, lwd=2)
shade(mu.PI, M_seq)
precis(m5.2)

##5.7 - Constructing a DAG
library(dagitty)
dag5.1 <- dagitty("dag{A -> D; A -> M; M ->D}")
coordinates(dag5.1) <- list(x=c(A=0, D=1, M=2), y=c(A=0, D=1, M=0))
drawdag(dag5.1)

##5.8
DMA_dag2 <- dagitty('dag{D <- A -> M}')
impliedConditionalIndependencies(DMA_dag2)

##5.9
DMA_dag1 <- dagitty('dag{D <- A -> M -> D}')
impliedConditionalIndependencies(DMA_dag1)

##5.10 multiple linear regression
m5.3 <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + bM*M + bA*A,
    a ~ dnorm(0,0.2),
    bM ~ dnorm(0,0.5),
    bA ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data = d
)
precis(m5.3)

##5.11
plot( coeftab(m5.1,m5.2,m5.3), par=c("bA","bM") )

##5.13
m5.4 <- quap(
  alist(
    M ~ dnorm(mu, sigma),
    mu <- a + bAM * A,
    a ~ dnorm(0,0.2),
    bAM ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data = d
)

##5.14
mu <- link(m5.4) #getting predicted value for each state
mu_mean <- apply(mu, 2, mean)
mu_resid <- d$M - mu_mean

##5.15
#call link without specifying new data so it uses original data
mu <- link(m5.3)

#summarize samples across cases
mu_mean <- apply(mu, 2, mean)
mu_PI <- apply(mu,2,PI)

#simulate observations
#again no new data, so uses original data
D_sim <- sim( m5.3 , n=1e4 )
D_PI <- apply(D_sim, 2, PI)

##5.16
plot(mu_mean ~ d$D)
abline(a=0, b=1, lty=2)

for (i in 1:nrow(d)) {
  lines(rep(d$D[i], 2), mu_PI[,i], col=rangi2)
}

##5.17
identify( x=d$D , y=mu_mean , labels=d$Loc )

##5.19
m5.3_A <- quap(
  alist(
    ## A -> D <- M
    D ~ dnorm(mu, sigma),
    mu <- a + bM*M + bA*A,
    a ~ dnorm(0, 0.2),
    bM ~ dnorm(0, 0.5),
    bA ~ dnorm(0, 0.5),
    sigma ~ dexp(1),
    ## A -> M
    M ~ dnorm(mu_M, sigma_M),
    mu_M <- aM + bAM*A,
    aM ~ dnorm(0, 0.2),
    bAM ~ dnorm(0, 0.5),
    sigma_M ~ dexp(1)
  ), data = d
)

precis(m5.3_A)

##5.20
A_seq <- seq(from=-2, to=2, length.out=30)

##5.21
#prep data
sim_dat <- data.frame(A=A_seq)

#simulate M and then D using A_seq
s <- sim(m5.3_A, data = sim_dat, vars=c("M", "D"))

##5.22
plot(sim_dat$A, colMeans(s$D), ylim=c(-2, 2), type = "l", xlab="manipulated A", ylab = "counterfactual D")
shade(apply(s$D, 2, PI), sim_dat$A)
mtext( "Total counterfactual effect of A on D" )

##5.23
#new data frame, standardized to mean 26.1 and std dev 1.24
sim2_dat <- data.frame( A=(c(20,30)-26.1)/1.24 ) #z-score standardization
s2 <- sim(m5.3_A, data = sim2_dat, vars=c("M", "D"))
mean(s2$D[,2] - s2$D[,1]) #expected causal effect of increasing median age at marriage from 20 to 30 on divorce rate

##5.24
sim_dat <- data.frame( M=seq(from=-2,to=2,length.out=30) , A=0 )
s <- sim(m5.3_A, data = sim_dat, vars="D")

plot( sim_dat$M , colMeans(s) , ylim=c(-2,2) , type="l" ,
      xlab="manipulated M" , ylab="counterfactual D" )
shade( apply(s,2,PI) , sim_dat$M )
mtext( "Total counterfactual effect of M on D" )

###Overthinking: Simulating counterfactuals
##5.25
#define a range of values that we want to assign to A
A_seq <- seq( from=-2 , to=2 , length.out=30 )

##5.26
#extract the posterior samples, because we’ll simulate observations for each set of samples
post <- extract.samples(m5.3_A)
#I used the with function, which saves us having to type post$ in front of every parameter name
#For each value in A_seq, generate 1 sample M per mean and sd (10,000 each as there are 10,000 models in post)
M_sim <- with(post, sapply(1:30, 
          function(i) rnorm( 1e3 , aM + bAM*A_seq[i] , sigma_M ) ) )

D_sim <- with( post , sapply( 1:30 ,
         function(i) rnorm( 1e3 , a + bA*A_seq[i] + bM*M_sim[,i] , sigma ) ) )                       
