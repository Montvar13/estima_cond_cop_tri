library(copula)

gumbel.cop <- gumbelCopula(2)
dev.new();par(mfrow = c(1,2))
persp (gumbel.cop, dCopula, xlab = "u", ylab = "v", zlab = "z", col = "yellow")
contour (gumbel.cop, dCopula, xlab = "u", ylab = "v", col = "red2", lwd = 2)
image(gumbel.cop, dCopula, xlab = "u", ylab = "v", col=rainbow(4))


clayton.cop <- claytonCopula(2)
dev.new();par(mfrow = c(1,2))
persp (clayton.cop, dCopula, xlab = "u", ylab = "v", zlab = "z",
       col = "royalblue")
contour (clayton.cop, dCopula, xlab = "u", ylab = "v", col = "darkorange2", lwd=2)

frank.cop <- frankCopula(-10)
dev.new();par(mfrow = c(1,2))
persp (frank.cop, dCopula, xlab = "u", ylab = "v", zlab = "z", col = "darkolivegreen2",
       theta = 10)
contour (frank.cop, dCopula, xlab = "u", ylab = "v", col = "goldenrod4", lwd = 2)


####
frank1.cop <- frankCopula(10)
par(mfrow = c(1,2))
persp (frank1.cop, dCopula, xlab = "u", ylab = "v", zlab = "z", col = "chocolate1")
contour (frank1.cop, dCopula, xlab = "u", ylab = "v", col = "green")


dens.clay <- function(u,v,t) {(u^{-t-1}*v^{-t-1}*(t+1)*(u^{-t} + v^{-t} - 1)^{-1/t - 2})*(u^{-t} + v^{-t} > 1)}