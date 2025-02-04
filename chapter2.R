#Using the binomial distribution to compute the likelihood of 6 successes in 9 trials where the proababilty of sucess in each trial is 0.5
dbinom(6, size=9, prob = 0.5)


##Grid approximation
#define grid --> different proportions of water
p_grid <- seq(from = 0, to=1, length.out=10)

#define prior
prior <- rep(1,10)

#compute likelihood at each value in the grid. Using binomial distribution as likelihood
likelihood <- dbinom(6, size = 9, prob = p_grid)

#compute product of likelihood and prior
unstd.posterior <- likelihood * prior

#standardize the posterior so it sums to 1
posterior <- unstd.posterior / sum(unstd.posterior)

plot(p_grid, posterior, type="b", xlab = "probability of water", ylab = "posterior probability")
mtext("20 points")



##2M1

p_grid <- seq(from = 0, to =1, length.out = 25)
prior <- rep(1,25)
likelihood <- dbinom(3, size = 3, prob = p_grid)
unstd.posterior <-  likelihood * prior
posterior <- unstd.posterior/sum(unstd.posterior)
plot(p_grid, posterior, type="b", xlab = "probability of water", ylab = "posterior probability")

##2M2
p_grid <- seq(from = 0, to =1, length.out = 25)
prior <- ifelse(p_grid<0.5, 0, 1)
likelihood <- dbinom(5, size = 7, prob = p_grid)
unstd.posterior <-  likelihood * prior
posterior <- unstd.posterior/sum(unstd.posterior)
plot(p_grid, posterior, type="b", xlab = "probability of water", ylab = "posterior probability")

