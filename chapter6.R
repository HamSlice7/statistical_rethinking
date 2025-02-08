##6.2
N <- 100 #number of individuals
set.seed(909)
height <- rnorm(N, 10, 2) #sim total height of each
leg_prop <- runif(N, 0.4, 0.5) #leg as proportion of height
leg_left <- leg_prop*height + rnorm(N, 0, 0.02) # sim left leg as proportion + error
leg_right <- leg_prop*height + rnorm(N, 0, 0.02) # sim right leg as proportion + error
d <- data.frame(height, leg_left, leg_right)

##6.3
m6.1 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + bl*leg_left + br*leg_right,
    a <- dnorm(10, 100),
    bl ~ dnorm(2, 10),
    br ~ dnorm(2,10),
    sigma ~ dexp(1)
  ), data = d
)
precis(m6.1)

##6.4
plot(precis(m6.1))

##6.5
post <- extract.samples(m6.1)
plot(bl ~ br, post, col=col.alpha(rangi2,0.1) , pch=16)

##6.6
sum_blbr <- post$bl + post$br
dens(sum_blbr, col=rangi2, lwd=2, xlab="sum of bl and br")

##6.7
m6.2 <- quap(
  alist(
    height ~ dnorm( mu , sigma ) ,
    mu <- a + bl*leg_left,
    a ~ dnorm( 10 , 100 ) ,
    bl ~ dnorm( 2 , 10 ) ,
    sigma ~ dexp( 1 )
  ) , data=d )
precis(m6.2)

##6.8
library(rethinking)
data(milk)
d <- milk
d$K <- standardize(d$kcal.per.g)
d$F <- standardize(d$perc.fat) 
d$L <- standardize(d$perc.lactose)

##6.9
# kcal.per.g regressed on perc.fat
m6.3 <- quap(
  alist(
    K ~ dnorm( mu , sigma ) ,
    mu <- a + bF*F ,
    a ~ dnorm( 0 , 0.2 ) ,
    bF ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data=d )
# kcal.per.g regressed on perc.lactose
m6.4 <- quap(
  alist(
    K ~ dnorm( mu , sigma ) ,
    mu <- a + bL*L ,
    a ~ dnorm( 0 , 0.2 ) ,
    bL ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data=d )

precis(m6.3)
precis(m6.4)

##6.10
m6.5 <- quap(
  alist(
    K ~ dnorm( mu , sigma ) ,
    mu <- a + bF*F + bL*L ,
    a ~ dnorm( 0 , 0.2 ) ,
    bF ~ dnorm( 0 , 0.5 ) ,
    bL ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ), data=d )
precis(m6.5)
    
    
##6.11
pairs( ~ kcal.per.g + perc.fat + perc.lactose , data=d , col=rangi2 )

##6.12
d <- milk
sim.coll <- function(r=0.9) {
  d$x <- rnorm(nrow(d), mean=r*d$perc.fat, sd=sqrt((1-r^2)*var(d$perc.fat)))
  m <- lm(kcal.per.g ~ perc.fat + x, data = d)
  sqrt( diag( vcov(m) ) )[2] # stddev of parameter
}
rep.sim.coll <- function(r=0.9, n=100) {
  stddev <- replicate(n, sim.coll(r))
  mean(stddev)
}

r.seq <- seq(from=0,to=0.99,by=0.01)
stddev <- sapply( r.seq , function(z) rep.sim.coll(r=z,n=100) )
plot( stddev ~ r.seq , type="l" , col=rangi2, lwd=2 , xlab="correlation" )


#6.13
set.seed(71)
N <- 100

#simulate initial heights
h0 <- rnorm(N, 10, 2)

#assign treatments and simulate fungus and growth
treatment <- rep(0:1, each=N/2)
fungus <- rbinom(N, size = 1, prob = 0.5-treatment*0.4)
h1 <- h0 + rnorm(N, 5 - 3*fungus)

# compose a clean data frame
d <- data.frame( h0=h0 , h1=h1 , treatment=treatment , fungus=fungus )
precis(d)

##6.14
sim_p <- rlnorm(1e4, 0, 0.25)
precis(sim_p)

##6.15
m6.6 <- quap(
  alist(
    h1 ~ dnorm(mu, sigma),
    mu <- h0*p,
    p ~ dlnorm(0, 0.25),
    sigma ~ dexp(1)
  ), data = d
)

precis(m6.6)

##6.16
m6.7 <- quap(
  alist(
    h1 ~ dnorm(mu, sigma),
    mu <- h0 * p,
    p <- a + bt*treatment + bf*fungus,
    a ~ dlnorm(0, 0.2),
    bt ~ dnorm( 0 , 0.5 ),
    bf ~ dnorm(0,0.5),
    sigma ~ dexp( 1 )
  ), data = d2
)
precis(m6.7)

##6.17
m6.8 <- quap(
  alist(
    h1 ~ dnorm( mu , sigma ),
    mu <- h0 * p,
    p <- a + bt*treatment,
    a ~ dlnorm( 0 , 0.2 ),
    bt ~ dnorm( 0 , 0.5 ),
    sigma ~ dexp( 1 )
  ), data = d2)
precis(m6.8)

##6.18
library(dagitty)
plant_dag <- dagitty( "dag {
  H_0 -> H_1
  F -> H_1
  T -> F
}")

coordinates( plant_dag ) <- list( x=c(H_0=0,T=2,F=1.5,H_1=1) ,
                                  y=c(H_0=0,T=0,F=0,H_1=0) )
drawdag( plant_dag )

##6.19
impliedConditionalIndependencies(plant_dag)

##6.20
set.seed(71)
N <- 1000
h0 <- rnorm(N,10,2)
treatment <- rep( 0:1 , each=N/2 )
M <- rbern(N)
fungus <- rbinom( N , size=1 , prob=0.5 - treatment*0.4 + 0.4*M )
h1 <- h0 + rnorm( N , 5 + 3*M )
d2 <- data.frame( h0=h0 , h1=h1 , treatment=treatment , fungus=fungus )

##6.21
d <- sim_happiness(seed=1977, N_years = 1000)
precis(d)

##6.22
d2 <- d[ d$age>17 , ] # only adults
d2$A <- ( d2$age - 18 ) / ( 65 - 18 )

##6.23
d2$mid <- d2$married + 1
m6.9 <- quap(
  alist(
    happiness ~ dnorm(mu, sigma),
    mu <- a[mid] + bA*A, 
    a[mid] ~ dnorm(0,1),
    bA ~ dnorm(0,2),
    sigma ~ dexp(1)
  ), data = d2
)
precis(m6.9, depth = 2)
