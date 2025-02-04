#4.2
prod( 1 + runif(12,0,0.1) )

#4.3
growth <- replicate( 10000 , prod( 1 + runif(12,0,0.1) ) )  
dens( growth , norm.comp=TRUE )

#4.4
big <- replicate( 10000 , prod( 1 + runif(12,0,0.5) ) )  
small <- replicate( 10000 , prod( 1 + runif(12,0,0.01) ) )
dens(big, norm.comp = TRUE)
dens(small, nrom.comp = TRUE)

#4.5
log.big <- replicate( 10000 , log(prod(1 + runif(12,0,0.5))) )
dens(log.big)

#4.6
w <- 6; n <- 9;  
p_grid <- seq(from=0,to=1,length.out=100) 
posterior <- dbinom(w,n,p_grid)*dunif(p_grid,0,1) 
posterior <- posterior/sum(posterior)


#4.7
library(rethinking) 
data(Howell1) 
d <- Howell1

#4.8
str(d)

#4.9
precis(d)

#4.10
d$height

#4.11
d2 <- d[ d$age >= 18 , ]
dens(d2$height)

#4.12
curve( dnorm( x , 178 , 20 ) , from=100 , to=250 )

#4.13
curve( dunif( x , 0 , 50 ) , from=-10 , to=60 )

#4.14
sample_mu <- rnorm( 1e4 , 178 , 20 )  
dens(sample_mu)
sample_sigma <- runif( 1e4 , 0 , 50 ) 
dens(sample_sigma)
prior_h <- rnorm( 1e4 , sample_mu , sample_sigma ) 
dens( prior_h )

#4.15
sample_mu <- rnorm( 1e4 , 178 , 100 ) 
prior_h <- rnorm( 1e4 , sample_mu , sample_sigma ) 
dens( prior_h )

#4.16
mu.list <- seq( from=150, to=160 , length.out=100 ) 
sigma.list <- seq( from=7 , to=9 , length.out=100 ) 
post <- expand.grid( mu=mu.list , sigma=sigma.list ) 
post$LL <- sapply( 1:nrow(post) , function(i) sum( dnorm( d2$height , post$mu[i] , post$sigma[i] , log=TRUE ) ) ) 
post$prod <- post$LL + dnorm( post$mu , 178 , 20 , TRUE ) + dunif( post$sigma , 0 , 50 , TRUE ) 
post$prob <- exp( post$prod - max(post$prod) )

#4.17
contour_xyz( as.vector(post$mu) , post$sigma , post$prob )

#4.18
image_xyz( post$mu , post$sigma , post$prob )

#4.19
sample.rows <- sample( 1:nrow(post) , size=1e4 , replace=TRUE , prob=post$prob )
sample.mu <- post$mu[ sample.rows ]
sample.sigma <- post$sigma[ sample.rows ]

#4.20
plot( sample.mu , sample.sigma , cex=0.5 , pch=16 , col=col.alpha(rangi2,0.1) )

#4.21
dens( sample.mu )  
dens( sample.sigma )

#4.22
PI( sample.mu )  
PI( sample.sigma )


#4.26
library(rethinking)  
data(Howell1) 
d <- Howell1 
d2 <- d[ d$age >= 18 , ]

#4.27
flist <- alist(  height ~ dnorm( mu , sigma ) , 
                 mu ~ dnorm( 178 , 20 ) , 
                 sigma ~ dunif( 0 , 50 ) )
#4.28
m4.1 <- quap( flist , data=d2 )

#4.29
precis( m4.1 )


#4.31
m4.2 <- quap( alist( height ~ dnorm( mu , sigma ) , 
                     mu ~ dnorm( 178 , 0.1 ) , 
                     sigma ~ dunif( 0 , 50 ) ) , data=d2 ) 
precis( m4.2 )


#4.34
library(rethinking)  
post <- extract.samples( m4.1 , n=1e4 )
head(post)

#4.37
library(rethinking)  
data(Howell1); 
d <- Howell1; 
d2 <- d[ d$age >= 18 , ] 
plot( d2$height ~ d2$weight )

#4.38
set.seed(2971) 
N <- 100 # 100 lines 
a <- rnorm( N , 178 , 20 )
b <- rnorm( N , 0 , 10 )

#4.39
plot( NULL , xlim=range(d2$weight) , ylim=c(-100,400) , xlab="weight" , ylab="height" ) 
abline( h=0 , lty=2 ) 
abline( h=272 , lty=1 , lwd=0.5 ) 
mtext( "b ~ dnorm(0,10)" ) 
xbar <- mean(d2$weight) 
for ( i in 1:N ) curve( a[i] + b[i]*(x - xbar) , 
                        from=min(d2$weight) , 
                        to=max(d2$weight) , 
                        add=TRUE , 
                        col=col.alpha("black",0.2) )
#4.40
b <- rlnorm( 1e4 , 0 , 1 )  
dens( b , xlim=c(0,5) , adj=0.1 )

#4.41
set.seed(2971)  
N <- 100 # 100 lines 
a <- rnorm( N , 178 , 20 ) 
b <- rlnorm( N , 0 , 1 )

plot( NULL , xlim=range(d2$weight) , ylim=c(-100,400) , xlab="weight" , ylab="height" ) 
abline( h=0 , lty=2 ) 
abline( h=272 , lty=1 , lwd=0.5 ) 
mtext( "b ~ dnorm(0,10)" ) 
xbar <- mean(d2$weight) 
for ( i in 1:N ) curve( a[i] + b[i]*(x - xbar) , 
                        from=min(d2$weight) , 
                        to=max(d2$weight) , 
                        add=TRUE , 
                        col=col.alpha("black",0.2) )

#4.42
# load data again, since it's a long way back library(rethinking) 
data(Howell1) 
d <- Howell1 
d2 <- d[ d$age >= 18 , ]
# define the average weight, x-bar 
xbar <- mean(d2$weight)
# fit model
m4.3 <- quap( 
  alist( 
    height ~ dnorm( mu , sigma ) , 
    mu <- a + b*( weight - xbar ) , 
    a ~ dnorm( 178 , 20 ) , 
    b ~ dlnorm( 0 , 1 ) , 
    sigma ~ dunif( 0 , 50 ) 
    ) , data=d2 )

#4.43
m4.3 <- quap( alist( 
  height ~ dnorm( mu , sigma ) , 
  mu <- a + exp(log_b)*( weight - xbar ) , 
  a ~ dnorm( 178 , 20 ) , 
  log_b ~ dnorm( 0 , 1 ) , 
  sigma ~ dunif( 0 , 50 ) ) , 
  data=d2 )

#4.44
precis( m4.3 )

#4.45
round( vcov( m4.3 ) , 3 )

#4.46
plot( height ~ weight , data=d2 , col=rangi2 )
post <- extract.samples( m4.3 )
a_map <- mean(post$a)
b_map <- mean(post$b)
curve( a_map + b_map*(x - xbar) , add=TRUE )

#4.47
post <- extract.samples( m4.3 )
post[1:5,]

#4.48
N <- 10
dN <- d2[ 1:N , ]
mN <- quap( 
  alist( 
    height ~ dnorm( mu , sigma ) ,
    mu <- a + b*( weight - mean(weight) ) , 
    a ~ dnorm( 178 , 20 ) , 
    b ~ dlnorm( 0 , 1 ) , 
    sigma ~ dunif( 0 , 50 )
    ) , data=dN )

#4.49
# extract 20 samples from the posterior
post <- extract.samples( mN , n=20 )
# display raw data and sample size
plot( dN$weight , dN$height , 
      xlim=range(d2$weight) , ylim=range(d2$height) , 
      col=rangi2 , xlab="weight" , ylab="height" ) 
mtext(concat("N = ",N))
# plot the lines, with transparency
for ( i in 1:20 )
  curve( post$a[i] + post$b[i]*(x-mean(dN$weight)) ,
         col=col.alpha("black",0.3) , add=TRUE )

#4.50
post <- extract.samples( m4.3 )
mu_at_50 <- post$a + post$b * ( 50 - xbar )

#4.51
dens( mu_at_50 , col=rangi2 , lwd=2 , xlab="mu|weight=50" )

#4.52
PI( mu_at_50 , prob=0.89 )

#4.53
mu <- link( m4.3 )
str(mu)

#4.54
# define sequence of weights to compute predictions for
# these values will be on the horizontal axis
weight.seq <- seq( from=25 , to=70 , by=1 )

# use link to compute mu
# for each sample from posterior
# and for each weight in weight.seq
mu <- link( m4.3 , data=data.frame(weight=weight.seq) )
str(mu)

#4.55
# use type="n" to hide raw data
plot(height ~ weight, d2, type="n")

#loop over samples and plot each mu value
for (i in 1:100) {
  points(weight.seq, mu[i,], pch=16, col=col.alpha(rangi2, 0.1)) }


#4.56
#summarize the distribution of mu
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI, prob=0.89)

#4.57
#plot raw data
#fading out points to make line and interval more visible
plot(height ~ weight, data = d2, col = col.alpha(rangi2, 0.5))

#plot the MAP line, aka the mean mu for each weight
lines(weight.seq, mu.mean)

#plot a shaded region for 89% PI
shade(mu.PI, weight.seq)

#4.58 under the hood of link()
post <- extract.samples(m4.3) 
mu.link <- function(weight) post$a + post$b*(weight - xbar) 
weight.seq <- seq(from=25, to=70, by=1)
mu <- sapply(weight.seq, mu.link) #plugging in desired weights into the linear functions from the posterior distribution
mu.mean <- apply(mu, 2, mean)
mu.CI <- apply(mu, 2, PI, prob=0.89)

#4.59
sim.height <- sim(m4.3, data = list(weight=weight.seq), n=1e4) #sampling 10000 heights from Gaussian distribution for the particular expected height corresponding to the weight

#4.60
height.PI <- apply(sim.height, 2, PI, prob = 0.95) # 89% posterior prediction interval of observable (according to the model) heights, across the values of weight in weight.seq

#4.61
#plot raw data
plot(height ~ weight, d2, col=col.alpha(rangi2, 0.5))

#Draw MAP line
lines(weight.seq, mu.mean)

#Draw HPDI region for line
shade(mu.PI, weight.seq)

#Draw PI region for simulated heights
shade(height.PI, weight.seq)


#4.63 - Under the hood for sim()
post <- extract.samples(m4.3)
weight.seq <- 25:70
#For each expected heights in the posterior distribution, sample 10,000 heights from a normal distribution and do this for each input weight
sim.height <- sapply(weight.seq, function(weight)
  rnorm(
    n=nrow(post),
    mean=post$a + post$b*(weight - xbar),
    sd=post$sigma ))

height.PI <- apply(sim.height, 2, PI, prob = 0.89)

#4.64
d <- Howell1
plot(height ~ weight, d)

#4.65 - adding a quadratic term
d$weight_s <- (d$weight - mean(d$weight))/sd(d$weight)
d$weight_s2 <- d$weight_s^2
m4.5 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b1*weight_s + b2*weight_s2,
    a ~ dnorm(178, 20),
    b1 ~ dlnorm(0,1),
    b2 ~ dnorm(0,1),
    sigma ~ dunif(0,50)
  ), data = d
)

#4.66
precis(m4.5)

#4.67
weight.seq <- seq(from = -2.2, to = 2, length.out=30)
pred_dat <- list(weight_s=weight.seq, weight_s2=weight.seq^2)
mu <- link(m4.5, data = pred_dat)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI, prob=0.89)
sim.height <- sim(m4.5, data = pred_dat)
height.PI <- apply(sim.height, 2, PI, prob=0.89)

#4.68
plot(height ~ weight_s, d, col = col.alpha(rangi2,0.5))
lines(weight.seq, mu.mean)
shade(mu.PI, weight.seq)
shade(height.PI, weight.seq)
