library(rethinking)
library(rstan)
data("chimpanzees")
d <- chimpanzees

#11.2 - create index variables to represent the different combinations of predictors
d$treatment <- 1 + d$prosoc_left + 2*d$condition

#11.3 verifying index variable creation worked as intended
xtabs(~ treatment + prosoc_left + condition, d)

#11.4 Constructing a single variable logistic regression with a flat prior for a
m11.1 <- quap(
  alist(
    pulled_left ~ dbinom(1,p),
    logit(p) <- a,
    a ~ dnorm(0,10)
  ), data = d
)

#11.5 - sample from the prior
set.seed(1999)
prior <- extract.prior(m11.1, n=1e4)

#11.6 using inverse-link function, inverse of logit is inv_logit, to get p
p <- inv_logit(prior$a)
dens(p, adj=0.1)

#11.7 - using flat prior for beta
m11.2 <- quap(
  alist(
    pulled_left ~ dbinom(1,p),
    logit(p) <- a + b[treatment],
    a ~ dnorm(0,1.5),
    b[treatment] ~ dnorm(0,10)
  ), data = d
)
set.seed(1999)
prior <- extract.prior(m11.2, n=1e4)
#probability of pulling left for each treatment based on flat prior for beta
p <- sapply( 1:4 , function(k) inv_logit( prior$a + prior$b[,k] ) )

#plot absolute prior difference between the first two treatments
dens(abs(p[,1] - p[,2]), adj=0.1)


#11.9 changing prior for b
m11.3 <- quap(
  alist(
    pulled_left ~ dbinom(1, p),
    logit(p) <- a + b[treatment],
    a ~ dnorm(0,1.5),
    b[treatment] ~ dnorm(0,0.5)
  ), data = d
)
set.seed(1999)
prior <- extract.prior(m11.3, n = 1e4)
p <- sapply( 1:4 , function(k) inv_logit( prior$a + prior$b[,k] ) )
mean( abs( p[,1] - p[,2] ) )

##11.10
#trimmed data list
dat_list <- list(
  pulled_left = d$pulled_left,
  actor = d$actor,
  treatment = as.integer(d$treatment)
)

##11.11
#Using Markov chain to estimate posterior
m11.4 <- ulam(
  alist(
    pulled_left ~ dbinom( 1 , p ) ,
    logit(p) <- a[actor] + b[treatment] ,
    a[actor] ~ dnorm( 0 , 1.5 ),
    b[treatment] ~ dnorm( 0 , 0.5 )
  ) , data=dat_list , chains=4 , log_lik=TRUE )
precis( m11.4 , depth=2 )

##11.12 Looking at the intercepts on the outcome scale
post <- extract.samples(m11.4)
p_left <- inv_logit(post$a)
plot(precis(as.data.frame(p_left)))

##11.3
labs <- c("R/N","L/N","R/P","L/P")
plot( precis( m11.4 , depth=2 , pars="b" ) , labels=labs )

##11.14
diffs <- list(
  db13 = inv_logit(post$b[,1]) - inv_logit(post$b[,3]),
  db24 = inv_logit(post$b[,2]) - inv_logit(post$b[,4])
)
plot(precis(diffs))

#11.15 - proportion of left pulls in each combination of actor and treatment
pl <- by(d$pulled_left, list(d$actor, d$treatment), mean)
pl[1,]


##11.17
#using fitted models to make predictions
dat <- list(actor = rep(1:7, each = 4), treatment=rep(1:4, times = 7))
p_post <- link(m11.4, data = dat)
p_mu <- apply(p_post, 2, mean)
p_ci <- apply(p_post, 2, PI)

##11.23
post <- extract.samples(m11.4)
mean(exp(post$b[,4] - post$b[,2]))

##11.24
data(chimpanzees)
d <- chimpanzees
d$treatment <- 1 + d$prosoc_left + 2*d$condition
d$side <- d$prosoc_left + 1 # right 1, left 2
d$cond <- d$condition + 1 # no partner 1, partner 2

d_aggregated <- aggregate(
  d$pulled_left ,
  list( treatment=d$treatment , actor=d$actor ,
        side=d$side , cond=d$cond ) ,
  sum )
colnames(d_aggregated)[5] <- "left_pulls"

##11.25
dat <- with( d_aggregated , list(
  left_pulls = left_pulls,
  treatment = treatment,
  actor = actor,
  side = side,
  cond = cond ) )

#now there are 18 trials on each row
m11.6 <- ulam(
  alist(
    left_pulls ~ dbinom( 18 , p ) ,
    logit(p) <- a[actor] + b[treatment] ,
    a[actor] ~ dnorm( 0 , 1.5 ) ,
    b[treatment] ~ dnorm( 0 , 0.5 )
  ), data = dat, chains = 4, log_lik = TRUE
)
precis(m11.6, depth = 2) # same as m11.4

##11.28
data(UCBadmit)
d <- UCBadmit

##11.29
dat_list <- list(
  admit = d$admit,
  applications = d$applications,
  gid = ifelse( d$applicant.gender=="male" , 1 , 2 ))

m11.7 <- ulam(
  alist(
    admit ~ dbinom(applications, p),
    logit(p) <- a[gid],
    a[gid] ~ dnorm(0,1.5)
  ), data = dat_list, chains = 4
)
precis(m11.7, depth = 2)
       
##11.30
#Calculating log-odds difference and probability difference for males vs females
post <- extract.samples(m11.7)
diff_a <- post$a[,1] - post$a[,2]
diff_p <- inv_logit(post$a[,1]) - inv_logit(post$a[,2])
precis( list( diff_a=diff_a , diff_p=diff_p ) )

##11.31
#checking the posterior predictions for the model
postcheck( m11.7 ) #see that women admissions are higher for most os the department which seems wrong 
# draw lines connecting points from same dept
for ( i in 1:6 ) {
  x <- 1 + 2*(i-1)
  y1 <- d$admit[x]/d$applications[x]
  y2 <- d$admit[x+1]/d$applications[x+1]
  lines( c(x,x+1) , c(y1,y2) , col=rangi2 , lwd=2 )
  text( x+0.5 , (y1+y2)/2 + 0.05 , d$dept[x] , cex=0.8 , col=rangi2 )
}

##11.32
#Constructing a model where department is a predictor along with gender
dat_list$dept_id <- rep(1:6, each=2)
m11.8 <- ulam(
  alist(
    admit ~ dbinom(applications, p),
    logit(p) <- a[gid] + delta[dept_id],
    a[gid] ~ dnorm(0,1.5),
    delta[dept_id] ~ dnorm(0,1.5)
  ), data = dat_list, chains=4, iter = 4000
)

precis(m11.8, depth = 2)

##11.33
#examining the difference in log-odds and probability of acceptance between men and women across departments
post <- extract.samples(m11.8)
diff_a <- post$a[,1] - post$a[,2]
diff_p <- inv_logit(post$a[,1]) - inv_logit(post$a[,2])
precis(list(diff_a = diff_a, diff_p = diff_p))

##11.34
pg <- with(dat_list, sapply(1:6, function(k)
  applications[dept_id==k]/sum(applications[dept_id==k])))
rownames(pg) <- c("male", "female")
colnames(pg) <- unique(d$dept)
round(pg,2)

postcheck(m11.8)
for ( i in 1:6 ) {
  x <- 1 + 2*(i-1)
  y1 <- d$admit[x]/d$applications[x]
  y2 <- d$admit[x+1]/d$applications[x+1]
  lines( c(x,x+1) , c(y1,y2) , col=rangi2 , lwd=2 )
  text( x+0.5 , (y1+y2)/2 + 0.05 , d$dept[x] , cex=0.8 , col=rangi2 )
}


##11.35
y <- rbinom(1e5, 1000, 1/1000)
c(mean(y), var(y))

##11.36
data(Kline)
d <- Kline
d

##11.37
d$P <- scale(log(d$population))
d$contact_id <- ifelse(d$contact=="high", 2, 1)

#11.38 plotting prior of norm(0,10) on log scale for the intercept


curve(dlnorm(x,0,10), from = 0, to = 100, n = 200)

##11.40
#more informative prior on the log norm scale
curve(dlnorm(x, 3, 0.5), from = 0, to = 100, n=200)

a <- rnorm(100, 3, 0.5)
for (i in 1:N) {
  print(exp(a)[i])
} 


##11.41
N <- 100
a <- rnorm(N, 3, 0.5)
b <- rnorm(N, 0, 10)
plot(NULL, xlim=c(-2,2), ylim=c(0,100))
for ( i in 1:N ) curve( exp( a[i] + b[i]*x ) , add=TRUE , col=grau() )


##11.42 - using tighter priors for population
set.seed(10)
N <- 100
a <- rnorm( N , 3 , 0.5 )
b <- rnorm( N , 0 , 0.2 )
plot( NULL , xlim=c(-2,2) , ylim=c(0,100) )
for ( i in 1:N ) curve( exp( a[i] + b[i]*x ) , add=TRUE , col=grau() )

##11.43
x_seq <- seq( from=log(100) , to=log(200000) , length.out=100 )
lambda <- sapply( x_seq , function(x) exp( a + b*x ) )
plot( NULL , xlim=range(x_seq) , ylim=c(0,500) , xlab="log population" , ylab="total tools" )
for ( i in 1:N ) lines( x_seq , lambda[i,] , col=grau() , lwd=1.5 )

##11.44
plot( NULL , xlim=range(exp(x_seq)) , ylim=c(0,500) , xlab="population" ,  ylab="total tools" )
for ( i in 1:N ) lines( exp(x_seq) , lambda[i,] , col=grau() , lwd=1.5 )

##11.45
dat <- list(
  T = d$total_tools,
  P = d$P,
  cid = d$contact_id
)

#intercept only model
m11.9 <- ulam(
  alist(
    T ~ dpois(lambda),
    log(lambda) <- a,
    a ~ dnorm(3, 0.5)
  ), data = dat, chains = 4, log_lik = TRUE
)

#interaction model
m11.10 <- ulam(
  alist(
    T ~ dpois(lambda),
    log(lambda) <- a[cid] + b[cid] * P,
    a[cid] ~ dnorm(3, 0.5),
    b[cid] ~ dnorm(0,0.2)
  ), data = dat, chains = 4, log_lik = TRUE
)

##11.46
compare( m11.9 , m11.10 , func=PSIS ) # some points are overfitting

##11.47
k <- PSIS( m11.10 , pointwise=TRUE )$k
plot( dat$P , dat$T , xlab="log population (std)" , ylab="total tools" ,
      col=rangi2 , pch=ifelse( dat$cid==1 , 1 , 16 ) , lwd=2 ,
      ylim=c(0,75) , cex=1+normalize(k) )

# set up the horizontal axis values to compute predictions at
ns <- 100
P_seq <- seq(from = -1.4, to = 3, length.out = ns)

# predictions for cid = 1 (low contact)
lambda <- link(m11.10, data = data.frame(P=P_seq, cid = 1))
lmu <- apply(lambda, 2, mean)
lci <- apply(lambda, 2, PI)
lines(P_seq, lmu, lty=2, lwd=1.5)
shade(lci, P_seq, xpd=TRUE)

#predictions for cid=2 (high contact)
lambda <- link(m11.10, data = data.frame(P=P_seq, cid = 2))
lmu <- apply(lambda, 2, mean)
lci <- apply(lambda, 2, PI)
lines(P_seq, lmu, lty = 1, lwd=1.5)
shade(lci, P_seq, xpd=TRUE)

##11.48
plot(d$population, d$total_tools, xlab = "population", ylab = "total tools", 
     col=rangi2, pch=ifelse(dat$cid==1, 1, 16), lwd = 2, ylim=c(0,75), cex=1+normalize(k))

ns <- 100 
P_seq <- seq(from = -5, to = 3, length.out = ns)
# 1.53 is sd of log(population)
# 9 is mean of log(population)
pop_seq <- exp(P_seq*1.53 + 9)

lambda <- link( m11.10 , data=data.frame( P=P_seq , cid=1 ) )
lmu <- apply( lambda , 2 , mean )
lci <- apply( lambda , 2 , PI )
lines( pop_seq , lmu , lty=2 , lwd=1.5 )
shade( lci , pop_seq , xpd=TRUE )

lambda <- link( m11.10 , data=data.frame( P=P_seq , cid=2 ) )
lmu <- apply( lambda , 2 , mean )
lci <- apply( lambda , 2 , PI )
lines( pop_seq , lmu , lty=1 , lwd=1.5 )
shade( lci , pop_seq , xpd=TRUE )

##11.49 - scientific model - come nback to this 
dat2 <- list( T=d$total_tools, P=d$population, cid=d$contact_id )
m11.11 <- ulam(
  alist(
    T ~ dpois( lambda ),
    lambda <- exp(a[cid])*P^b[cid]/g,
    a[cid] ~ dnorm(1,1),
    b[cid] ~ dexp(1),
    g ~ dexp(1)
  ), data=dat2 , chains=4 , log_lik=TRUE )


plot(d$population, d$total_tools, xlab = "population", ylab = "total tools", 
     col=rangi2, pch=ifelse(dat$cid==1, 1, 16), lwd = 2, ylim=c(0,75), cex=1+normalize(k))

ns <- 100 
P_seq <- seq(from = -5, to = 3, length.out = ns)
# 1.53 is sd of log(population)
# 9 is mean of log(population)
pop_seq <- exp(P_seq*1.53 + 9)

lambda <- link( m11.11 , data=data.frame( P=P_seq , cid=1 ) )
lmu <- apply( lambda , 2 , mean )
lci <- apply( lambda , 2 , PI )
lines( pop_seq , lmu , lty=2 , lwd=1.5 )
shade( lci , pop_seq , xpd=TRUE )

lambda <- link( m11.11 , data=data.frame( P=P_seq , cid=2 ) )
lmu <- apply( lambda , 2 , mean )
lci <- apply( lambda , 2 , PI )
lines( pop_seq , lmu , lty=1 , lwd=1.5 )
shade( lci , pop_seq , xpd=TRUE )


##11.50
num_day <- 30
y <- rpois(num_day, 1.5)

##11.51
num_weeks <- 4
y_new <- rpois(num_weeks, 0.5*7)

##11.52
y_all <- c(y, y_new)
exposure <- c( rep(1,30) , rep(7,4) )
monastery <- c( rep(0,30) , rep(1,4) )
d <- data.frame( y=y_all , days=exposure , monastery=monastery )


##11.53
#compute the offset
d$log_days <- log( d$days )

# fit the model
m11.12 <- quap(
  alist(
    y ~ dpois( lambda ),
    log(lambda) <- log_days + a + b*monastery,
    a ~ dnorm( 0 , 1 ),
    b ~ dnorm( 0 , 1 )
  ), data=d )


##11.54
post <- extract.samples(m11.12)
lambda_old <- exp(post$a)
lambda_new <- exp(post$a + post$b)
precis(data.frame(lambda_old, lambda_new))

##11.55
#simulate career choices among 500 individuals
N <- 500 # number of individuals
income <- c(1,2,5) #expected income of each career
score <- 0.5*income #score for each income
#next line converts scores to probabilities
p <- softmax(income[1], income[2], income[3])

#now simulate choice
#outcome career holds event type values, not counts
career <- rep(NA, N) #empty vector of choices for each individual
# sample chosen career for each individual
set.seed(34302)
for (i in 1:N) career[i] <- sample(1:3, size = 1, prob = p)

##11.56
code_m11.13 <- "
data{
int N; // number of individuals
int K; // number of possible careers
int career[N]; // outcome
vector[K] career_income;
}
parameters{
vector[K-1] a; // intercepts
real<lower=0> b; // association of income with choice
}
model{
vector[K] p;
vector[K] s;
a ~ normal( 0 , 1 );
b ~ normal( 0 , 0.5 );
s[1] = a[1] + b*career_income[1];
s[2] = a[2] + b*career_income[2];
s[3] = 0; // pivot
p = softmax( s );
career ~ categorical( p );
}"


##11.57
dat_list <- list( N=N , K=3 , career=career , career_income=income )
m11.13 <- stan( model_code=code_m11.13 , data=dat_list , chains=4 )
precis( m11.13 , 2 )



##11.58
post <- extract.samples( m11.13 )

#set up logit scores
s1 <- with(post, a[,1] + b*income[1])
s2_orig <- with( post , a[,2] + b*income[2] )
s2_new <- with( post , a[,2] + b*income[2]*2 )

#compute probabilities for original and counterfactual
p_orig <- sapply(1:length(post$b), function(i) 
  softmax(c(s1[i], s2_orig[i], 0)))
p_new <- sapply( 1:length(post$b) , function(i)
  softmax( c(s1[i],s2_new[i],0) ) )

# summarize
p_diff <- p_new[2,] - p_orig[2,]
precis( p_diff )

##11.59
N <- 500
#simulate family incomes for each individual
family_income <- runif(N)
#assign a unique coefficient for each type of event
b <- c(-2,0,2)
career <- rep(NA,N) #empty vector of choices for each individual
for (i in 1:N){
  score <- 0.5*(1:3) + b*family_income[i]
  p <- softmax(score[1],score[2],score[3])
  career[i] <- sample(1:3, size = 1, prob = p)
}


##11.60
