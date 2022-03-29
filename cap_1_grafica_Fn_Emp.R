Fn <- function(x,m) {sapply(x,function(x) mean( m <= x))}
m <- 10
ma.obs10 <- rnorm(m,0,1)
rang_obs10 <- seq(min(ma.obs10), max(ma.obs10), length.out = m)
evFn10 <- Fn(rang_obs10, ma.obs10)
plot(rang_obs10, evFn10, type = "s", xlim = c(min(ma.obs10), max(ma.obs10)), ylim = c(0,1),
    ylab = "Fn. empírica.", xlab = "observaciones", main = "Fn. con muestra tamaño 10")
points(rang_obs10, evFn10)
aux <- seq(min(ma.obs10), max(ma.obs10), length.out = 1000)
lines(aux, pnorm(aux), col = "red", lwd = 3)

dev.new()
m2 <- 100
ma.obs100 <- rnorm(m2,0,1)
rang_obs100 <- seq(min(ma.obs100), max(ma.obs100), length.out = m2)
evFn100 <- Fn(rang_obs100, ma.obs100)
plot(rang_obs, evFn100, xlim = c(min(ma.obs), max(ma.obs)), ylim = c(0,1),
    ylab = "Fn. empírica.", xlab = "observaciones", main = "Fn. con muestra tamaño 100")
aux2 <- seq(min(ma.obs), max(ma.obs), length.out = 1000)
lines(aux2, pnorm(aux), col = "red", lwd = 3)
dev.new()
m3 <- 1000
ma.obs1000 <- rnorm(m3,0,1)
rang_obs1000 <- seq(min(ma.obs1000), max(ma.obs1000), length.out = m3)
evFn1000 <- Fn(rang_obs1000, ma.obs1000)
plot(rang_obs1000, evFn1000, xlim = c(min(ma.obs1000), max(ma.obs1000)), ylim = c(0,1),
    ylab = "Fn. empírica.", xlab = "observaciones", main = "Fn. con muestra tamaño 1000")
aux3 <- seq(min(ma.obs), max(ma.obs), length.out = 1000)
lines(aux3, pnorm(aux), col = "red", lwd = 3)

dev.new()
m4 <- 50
ma.obs50 <- rnorm(m4,0,1)
rang_obs50 <- seq(min(ma.obs50), max(ma.obs50), length.out = m4)
evFn50 <- Fn(rang_obs50, ma.obs50)
plot(rang_obs50, evFn50, xlim = c(min(ma.obs50), max(ma.obs50)), ylim = c(0,1),
    ylab = "Fn. empírica.", xlab = "observaciones", main = "Fn. con muestra tamaño 50")
aux4 <- seq(min(ma.obs50), max(ma.obs50), length.out = 1000)
lines(aux4, pnorm(aux), col = "red", lwd = 3)
####================

dev.new(); par(mfrow = c(2,2))
plot(rang_obs10, evFn10, type = "s", xlim = c(min(ma.obs10), max(ma.obs10)), ylim = c(0,1),
    ylab = "Fn. empírica.", xlab = "observaciones", main = "Fn. con muestra tamaño 10")
points(rang_obs10, evFn10)
aux <- seq(min(ma.obs10), max(ma.obs10), length.out = 1000)
lines(aux, pnorm(aux), col = "red", lwd = 3)

plot(rang_obs50, evFn50, xlim = c(min(ma.obs50), max(ma.obs50)), ylim = c(0,1),
    ylab = "Fn. empírica.", xlab = "observaciones", main = "Fn. con muestra tamaño 50")
aux4 <- seq(min(ma.obs50), max(ma.obs50), length.out = 1000)
lines(aux4, pnorm(aux), col = "red", lwd = 3)


plot(rang_obs100, evFn100, xlim = c(min(ma.obs100), max(ma.obs100)), ylim = c(0,1),
    ylab = "Fn. empírica.", xlab = "observaciones", main = "Fn. con muestra tamaño 100")
aux <- seq(min(ma.obs100), max(ma.obs100), length.out = 1000)
lines(aux, pnorm(aux), col = "red", lwd = 3)


plot(rang_obs1000, evFn1000, xlim = c(min(ma.obs1000), max(ma.obs1000)), ylim = c(0,1),
    ylab = "Fn. empírica.", xlab = "observaciones", main = "Fn. con muestra tamaño 1000")
aux <- seq(min(ma.obs1000), max(ma.obs1000), length.out = 1000)
lines(aux, pnorm(aux), col = "red", lwd = 3)

