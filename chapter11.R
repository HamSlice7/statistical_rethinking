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
