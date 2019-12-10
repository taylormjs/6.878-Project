require(glm2, quietly = TRUE)
data(crabs)
crabs.boot <- crabs[crabs$Rep1, -c(5:6)]

# This algorithm doesn't converge.
t.glm <- system.time(
  fit.glm <- glm(Satellites ~ Width + Dark + GoodSpine, data = crabs.boot, family = poisson(identity),
                 start = rep(1, 4), maxit = 500)
)

# Make some addreg object.
t.cem <- system.time(
  fit.cem <- addreg(Satellites ~ Width + Dark + GoodSpine, data = crabs.boot, family = poisson,
                    start = rep(1, 4), method="cem")
)

fit.em <- addreg(Satellites ~ Width + Dark + GoodSpine, data = crabs.boot, family = poisson,
                 start = rep(1, 4), method="em")

# Run the EM algorithm to fit it.
t.cem.acc <- system.time(fit.cem.acc <- update(fit.cem, accelerate = "squarem"))
t.em.acc <- system.time(fit.em.acc <- update(fit.em, accelerate = "squarem"))