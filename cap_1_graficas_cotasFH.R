###CURVAS DE NIVEL DE LAS COTAS DE FRÉCHET-HOEFFDING

#FUNCIÓN SIMILAR A OUTER.
rat <- function(x,y,f) { 
 mat <- matrix( ,length(x), length(y))
 for(i in 1:length(x)) { 
 a <- rep(x[i], length(x))
 mat[ i, ] <- mapply(f, a, y)
 }
 mat
}

m <- 25
u <- seq(0, 1, length = m)
v <- seq(0, 1, length = m)
W <- function(u,v) { max(u + v - 1, 0) }
z <- rat(u, v, W)

image(u, v, z, main = "Cota inferior de Fréchet Hoeffding")
filled.contour(u, v, z, color = function(n) hcl.colors(n, "Oranges"),
       plot.title = title(main = "Cota inferior de Fréchet-Hoeffding",
       xlab = "u", ylab = "v"),
       plot.axes = {axis(1, seq(0, 1, by = 0.2))
                    axis(2, seq(0, 1, by = 0.2)) }
)
dev.new()
persp(u, v, z, col = "gold", theta = 355, phi = 20,
     main = "Cota inferior de Fréchet-Hoeffding",
     xlab = "u", ylab = "v", zlab = "W(u,v)", ) 




#=====
m <- 25
u <- seq(0, 1, length = m)
v <- seq(0, 1, length = m)
zp <- rat(u, v, min)

filled.contour(u, v, zp, color = function(n) hcl.colors(n, "Oranges"),
       plot.title = title(main = "Cota superior de Fréchet-Hoeffding",
       xlab = "u", ylab = "v"),
       plot.axes = {axis(1, seq(0, 1, by = 0.2))
                    axis(2, seq(0, 1, by = 0.2)) }
)

dev.new()
persp(u, v, zp, col = "gold", theta = 355, phi = 20,
     main = "Cota superior de Fréchet-Hoeffding",
     xlab = "u", ylab = "v", zlab = "M(u,v)") 



#============================================================================

diagW <- function(t) { max(2*t-1,0) }  
t <- seq(0, 1, length.out = 10000)
plot(t, sapply(t, diagW), col = "blue", lwd = 1, xlab = "t", ylab ="t", 
main = "Secciones diagonales") 
lines(t, t, col = "red", lwd = 3)  
legend(0.1, 0.99, legend = c("Cópula M", "Cópula W"), col = c("red", "blue"),
       lty = 1, cex = 1)
legend(0.4, 0.38, legend = "C", col = "white", box.lty=0)
#===========================================================================
u <- seq(0, 1, length = 25)
v <- seq(0, 1, length = 25)
pi <- function(u,v) {u*v}
z <- outer(u,v,pi)
dev.new()
filled.contour(u, v, z, color = function(n) hcl.colors(n, "Oranges"),
       plot.title = title(main = "Cópula producto.",
       xlab = "u", ylab = "v"),
       plot.axes = {axis(1, seq(0, 1, by = 0.2))
                    axis(2, seq(0, 1, by = 0.2)) }
)
dev.new()
persp(u, v, z, col = "gold", theta = 355, phi = 20,
     main = "Cópula producto",
     xlab = "u", ylab = "v", zlab = "M(u,v)") 
dev.new()
diagpi <- function(t) {t^2}
t <- seq(0, 1, length.out = 10000)
plot(t, sapply(t, diagpi), col = "green", lwd = 1, xlab = "t", ylab = "t", 
    main = "Sección diagonal de la cópula producto") 

dev.new()
plot(t, sapply(t, diagW), col = "blue", lwd = 1, xlab = "t", ylab ="t", 
main = "Secciones diagonales") 
lines(t, t, col = "red", lwd = 3)
lines(t, sapply(t, diagpi), col = "green", lwd = 2)  
legend(0.1, 0.99, legend = c("Cópula M", "Cópula producto", "Cópula W"), col = c("red", "green", "blue"),
       lty = 1, cex = 1)

#======================00
fx <- function(x){dexp(x, 1)}
fy <- function(y){dgamma(y, 2, 1)}
x <- seq(0, 5, length.out = 1000)
y <- seq(0, 10, length.out = 1000)
dev.new(); par(mfrow = c(1,2))
plot(x, fx(x), main = "Densidad de X", col = "red", lwd = 3, type = "l")
plot(y, fy(y), main = "Densidad de Y", col = "blue", lwd = 2, type = "l")

