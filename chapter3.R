library(rethinking)
data("Howell1")

#Restrict to individuals above 18
d <- Howell1[Howell1$age>=18,]

#function to simulate weights of individuals from weight 
sim_weight <- function(H, b, sd) {
  U <- rnorm(length(H), 0, sd)
  W <- b*H + U
  return(W)
}

#simulate height data - uniform distribution of heights
H <- runif(200, min = 130, max = 170)

W <- sim_weight(H, b=0.5, sd=5)

plot(W ~ H, col=2, lwd=3)


mu.list <- seq( from=150, to=160 , length.out=100 )
sigma.list <- seq( from=7 , to=9 , length.out=100 )
post <- expand.grid( mu=mu.list , sigma=sigma.list )
post$LL <- sapply( 1:nrow(post) , function(i) sum( dnorm( d$height , post$mu[i] , post$sigma[i] , log=TRUE ) ) )
post$prod <- post$LL + dnorm( post$mu , 178 , 20 , TRUE ) +
  dunif( post$sigma , 0 , 50 , TRUE )
dunif( post$sigma , 0 , 50 , TRUE )

dunif( post$sigma , 0 , 50 , TRUE )

image_xyz( post$mu , post$sigma , post$prob )
