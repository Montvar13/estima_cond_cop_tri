# Código para gráficas y tablas del capítulo 3.

#PAQUETES REQUERIDOS

library(ADGofTest)
library(subcopem2D)
library(scatterplot3d)

#FUNCIONES NECESARIAS.

#X
fx <- function(x) dunif(x,0,1)
invFx <- function(u)qunif(u,0,1)

#Y
fy <- function(y) {-log(y)*((0<y)&(y<1))}
invFylx <- function(v,x)qunif(v,0,x)

#Z
#Teta = 0
fz_0 <- function(z) {-log(z)*((0<z)&(z<1))}
invFz1xy_0 <- function(w,x,y) qunif(w, 0, x)

#Teta = 1
fz_1 <- function(z) {(log(z) + 2/sqrt(z) - 2)*((0 < z)&(z < 1))}
invFz1xy_1 <- function(w,x,y) qunif(w, 0, x*y)

#Teta>0 y theta != 1
fz_t <- function(z,teta) {(z^(-teta/(teta + 1)) - 1)/(teta - teta^2) + 
                           (teta*z^((1 - teta)/teta)*(1 - 
                                                      1/(z^(1/(teta^2 + teta)))))/(1 - 
                                                                                   teta)}
invFz1xy_t <- function(w,x,y,t) qunif(w, 0, x*y^t)


# SIMULACIÓN DE (X,Y,Z) PARA TETA = 0.

m <- 3000
u <- runif(m,0,1)
x_0 <- invFx(u)
v <- runif(m,0,1)
y_0 <- mapply(invFylx, v, x_0)
w <- runif(m,0,1)
z_0 <- mapply(invFz1xy_0, w, x_0, y_0)
ma.obs_0 <- cbind(x_0,y_0,z_0)
head(ma.obs_0)

# SIMULACIONES PARA TETA = 0.1
m <- 3000
t <- 0.1
u <- runif(m,0,1)
x_01 <- invFx(u)
v <- runif(m,0,1)
y_01 <- mapply(invFylx, v, x_01)
w <- runif(m,0,1)
z_01 <- mapply(invFz1xy_t, w, x_01, y_01, rep(t, m))
ma.obs_01 <- cbind(x_01,y_01,z_01)
head(ma.obs_01)

# ANÁLISIS DE LAS OBSERVACIONES

summary(x_0)
summary(y_0)
summary(z_0)		

#Figura 3.7

dev.new(); par(mfrow = c(2,3))

hist(x_0, prob = T, col = "yellow", main = "Histograma de X", xlab = "x", breaks = 10)
rx_0 <- seq(min(x_0), max(x_0), length.out = 10000)
lines(rx_0, fx(rx_0), col = "red")

hist(y_0, prob = T, col = "yellow", main = "Histograma de Y", xlab = "y", breaks = 10)
ry_0 <- seq(min(y_0), max(y_0), length.out = 10000)
lines(ry_0, fy(ry_0), col = "red")

hist(z_0, prob = T, col = "yellow", main = "Histograma de Z", xlab = "z", breaks = 10)
rz_0 <- seq(min(z_0), max(z_0), length.out = 10000)
lines(rz_0, fz_0(rz_0), col = "red")

hist(x_01, prob = T, col = "yellow", main = "Histograma de X", xlab = "x", breaks = 10)
rx_01 <- seq(min(x_01), max(x_01), length.out = 10000)
lines(rx_01, fx(rx_01), col = "red")

hist(y_01, prob = T, col = "yellow", main = "Histograma de Y", xlab = "y", breaks = 10)
ry_01 <- seq(min(y_01), max(y_01), length.out = 10000)
lines(ry_01, fy(ry_01), col = "red")

hist(z_01, prob = T, col = "yellow", main = "Histograma de Z", xlab = "z", breaks = 10)
rz_01 <- seq(min(z_01), max(z_01), length.out = 10000)
lines(rz_01, mapply(fz_t, rz_01, rep(t,m)), col = "red")


# SIMULACIONES PARA TETA = 0.5

m <- 3000
t <- 0.5
u <- runif(m,0,1)
x_05 <- invFx(u)
v <- runif(m,0,1)
y_05 <- mapply(invFylx, v, x_05)
w <- runif(m,0,1)
z_05 <- mapply(invFz1xy_t, w, x_05, y_05, rep(t, m))
ma.obs_05 <- cbind(x_05,y_05,z_05)
head(ma.obs_05)


# SIMULACIONES PARA TETA = 0.7

m <- 3000
t <- 0.7
u <- runif(m,0,1)
x_07 <- invFx(u)
v <- runif(m,0,1)
y_07 <- mapply(invFylx, v, x_07)
w <- runif(m,0,1)
z_07 <- mapply(invFz1xy_t, w, x_07, y_07, rep(t, m))
ma.obs_07 <- cbind(x_07,y_07,z_07)
head(ma.obs_07)

###ANÁLISIS DE LAS OBSERVACIONES

summary(x_05)
summary(y_05)
summary(z_05)

# Figuras 3.8

dev.new(); par(mfrow = c(2,3))

hist(x_05, prob = T, col = "yellow", main = "Histograma de X", xlab = "x", breaks = 10)
rx_05 <- seq(min(x_05), max(x_05), length.out = 10000)
lines(rx_05, fx(rx_05), col = "red")

hist(y_05, prob = T, col = "yellow", main = "Histograma de Y", xlab = "y", breaks = 10)
ry_05 <- seq(min(y_05), max(y_05), length.out = 10000)
lines(ry_05, fy(ry_05), col = "red")

hist(z_05, prob = T, col = "yellow", main = "Histograma de Z", xlab = "z", breaks = 10)
rz_05 <- seq(min(z_05), max(z_05), length.out = 10000)
lines(rz_05, mapply(fz_t, rz_05, rep(t,m)), col = "red")

hist(x_07, prob = T, col = "yellow", main = "Histograma de X", xlab = "x", breaks = 10)
rx_07 <- seq(min(x_07), max(x_07), length.out = 10000)
lines(rx_07, fx(rx_07), col = "red")

hist(y_07, prob = T, col = "yellow", main = "Histograma de Y", xlab = "y", breaks = 10)
ry_07 <- seq(min(y_07), max(y_07), length.out = 10000)
lines(ry_07, fy(ry_07), col = "red")

hist(z_07, prob = T, col = "yellow", main = "Histograma de Z", xlab = "z", breaks = 10)
rz_07 <- seq(min(z_07), max(z_07), length.out = 10000)
lines(rz_07, mapply(fz_t, rz_07, rep(t,m)), col = "red")



# SIMULACIONES PARA TETA = 1

m <- 3000
u <- runif(m,0,1)
x_1 <- invFx(u)
v <- runif(m,0,1)
y_1 <- mapply(invFylx, v, x_1)
w <- runif(m,0,1)
z_1 <- mapply(invFz1xy_1, w, x_1, y_1)
ma.obs_1 <- cbind(x_1,y_1,z_1)
head(ma.obs_1)

# SIMULACIONES PARA TETA = 1.1

m <- 3000
t <- 1.1
u <- runif(m,0,1)
x_11 <- invFx(u)
v <- runif(m,0,1)
y_11 <- mapply(invFylx, v, x_11)
w <- runif(m,0,1)
z_11 <- mapply(invFz1xy_t, w, x_11, y_11, rep(t, m))
ma.obs_11 <- cbind(x_11,y_11,z_11)
head(ma.obs_11)

# ANÁLISIS DE LAS OBSERVACIONES

summary(x_1)
summary(y_1)
summary(z_1)

#Figura 3.9

dev.new(); par(mfrow = c(2,3))

hist(x_1, prob = T, col = "yellow", main = "Histograma de X", xlab = "x", breaks = 10)
rx_1 <- seq(min(x_1), max(x_1), length.out = 10000)
lines(rx_1, fx(rx_1), col = "red")

hist(y_1, prob = T, col = "yellow", main = "Histograma de Y", xlab = "y", breaks = 10)
ry_1 <- seq(min(y_1), max(y_1), length.out = 10000)
lines(ry_1, fy(ry_1), col = "red")

hist(z_1, prob = T, col = "yellow", main = "Histograma de Z", xlab = "z", breaks = 10)
rz_1 <- seq(min(z_1), max(z_1), length.out = 10000)
lines(rz_1, fz_1(rz_1), col = "red")

hist(x_11, prob = T, col = "yellow", main = "Histograma de X", xlab = "x", breaks = 10)
rx_11 <- seq(min(x_11), max(x_11), length.out = 10000)
lines(rx_11, fx(rx_11), col = "red")

hist(y_11, prob = T, col = "yellow", main = "Histograma de Y", xlab = "y", breaks = 10)
ry_11 <- seq(min(y_11), max(y_11), length.out = 10000)
lines(ry_11, fy(ry_11), col = "red")

hist(z_11, prob = T, col = "yellow", main = "Histograma de Z", xlab = "z", breaks = 10)
rz_11 <- seq(min(z_11), max(z_11), length.out = 10000)
lines(rz_11, mapply(fz_t, rz_11, rep(t,m)), col = "red")


# SIMULACIONES PARA TETA = 1.5

m <- 3000
t <- 1.5
u <- runif(m,0,1)
x_15 <- invFx(u)
v <- runif(m,0,1)
y_15 <- mapply(invFylx, v, x_15)
w <- runif(m,0,1)
z_15 <- mapply(invFz1xy_t, w, x_15, y_15, rep(t, m))
ma.obs_15 <- cbind(x_15,y_15,z_15)
head(ma.obs_15)

# SIMULACIONES PARA TETA = 2

m <- 3000
t <- 2
u <- runif(m,0,1)
x_2 <- invFx(u)
v <- runif(m,0,1)
y_2 <- mapply(invFylx, v, x_2)
w <- runif(m,0,1)
z_2 <- mapply(invFz1xy_t, w, x_2, y_2, rep(t, m))
ma.obs_2 <- cbind(x_2,y_2,z_2)
head(ma.obs_2)

# ANÁLISIS DE LAS OBSERVACIONES

summary(x_15)
summary(y_15)
summary(z_15)

#Figura 3.10

dev.new(); par(mfrow = c(2,3))

hist(x_15, prob = T, col = "yellow", main = "Histograma de X", xlab = "x", breaks = 10)
rx_15 <- seq(min(x_15), max(x_15), length.out = 10000)
lines(rx_15, fx(rx_15), col = "red")

hist(y_15, prob = T, col = "yellow", main = "Histograma de Y", xlab = "y", breaks = 10)
ry_15 <- seq(min(y_15), max(y_15), length.out = 10000)
lines(ry_15, fy(ry_15), col = "red")

hist(z_15, prob = T, col = "yellow", main = "Histograma de Z", xlab = "z", breaks = 10)
rz_15 <- seq(min(z_15), max(z_15), length.out = 10000)
lines(rz_15, mapply(fz_t, rz_15, rep(t,m)), col = "red")

hist(x_2, prob = T, col = "yellow", main = "Histograma de X", xlab = "x", breaks = 10)
rx_2 <- seq(min(x_2), max(x_2), length.out = 10000)
lines(rx_2, fx(rx_2), col = "red")

hist(y_2, prob = T, col = "yellow", main = "Histograma de Y", xlab = "y", breaks = 10)
ry_2 <- seq(min(y_2), max(y_2), length.out = 10000)
lines(ry_2, fy(ry_2), col = "red")

hist(z_2, prob = T, col = "yellow", main = "Histograma de Z", xlab = "z", breaks = 10)
rz_2 <- seq(min(z_2), max(z_2), length.out = 10000)
lines(rz_2, mapply(fz_t, rz_2, rep(t,m)), col = "red")



# Gráficas de dispersión 3D:

# OBSERVACIONES (X,Y,Z)

# Teta = 0
# Figura 3.11
dev.new()
scatterplot3d(ma.obs_0[,1:3], angle = 45, color = "steelblue", pch = 16,
              main = "Observaciones de (X,Y,Z)",
              xlab = "x", ylab = "y", zlab = "z")

# Teta = 0.5
# Figura 3.12
dev.new()
scatterplot3d(ma.obs_05[,1:3], angle = 45, color = "steelblue", pch = 16,
              main = "Observaciones de (X,Y,Z)",
              xlab = "x", ylab = "y", zlab = "z")

# Teta = 1
# Figura 3.13
dev.new()
scatterplot3d(ma.obs_1[,1:3], angle = 45, color = "steelblue", pch = 16,
              main = "Observaciones de (X,Y,Z)",
              xlab = "x", ylab = "y", zlab = "z")

# Teta = 1.5
# Figura 3.14
dev.new()
scatterplot3d(ma.obs_15[,1:3], angle = 45, color = "steelblue", pch = 16,
              main = "Observaciones de (X,Y,Z)",
              xlab = "x", ylab = "y", zlab = "z")

# Teta = 2
# Figura 3.15
dev.new()
scatterplot3d(ma.obs_2[,1:3], angle = 45, color = "steelblue", pch = 16,
              main = "Observaciones de (X,Y,Z)",
              xlab = "x", ylab = "y", zlab = "z")


# GRÁFICAS DE LAS SIMULACIONES DE (X,Y|Z), (Y,Z|X) Y (X,Z|Y) PARA 
# z = (0.25, 0.75) 
# teta = (0, 0.5)

# Funciones:

#Z
#teta = 0
Fz_0 <- function(z) {(z*(1-log(z)))*((0 < z)&(z < 1)) + 1*(z>=1)}
#teta = 1 
Fz_1 <- function(z) {(z*log(z)+4*sqrt(z)-3*z)*((0 < z)&(z < 1)) + 1*(z>=1)}
#teta > 0 & teta != 1
Fz_t <- function(z,t) {((1/(1-t))*((((t+1)*(1-t^2))/(t))*z^(1/(t+1)) - 
                       z/t + t^2*z^(1/t)))*((0 <= z)&(z < 1)) + 1*(z>=1)}
invFz_0 <- function(w) {uniroot(function(z) Fz_0(z) - w, 
                        interval = c(0.00000000001,0.999))$root}
invFz_1 <- function(w) {uniroot(function(z) Fz_1(z) - w, 
                        interval = c(0.00000000001,0.999))$root}
invFz_t <- function(w,t){uniroot(function(z) Fz_t(z,t) - w, 
                         interval = c(0.00000000000001,0.999))$root}
#X|Z
#teta = 0
Fx1z_0 <- function(x,z) {(1-log(x)/log(z))*((z <= x)&(x < 1)) + 1*(x >= 1)} 
#teta =1
Fx1z_1 <- function(x,z) {((1/(log(z)+2/sqrt(z)-2))*(2/sqrt(z) + 
                         (log(z)-2*log(x)-2)/x))*((sqrt(z) <= x)&(x < 1)) + 1*(x >= 1)}
#teta > 0 & teta != 1
Fx1z_t <- function(x,z,t){((z^(-t/(t+1)) - x^(-t) + 
                            t^2*z^((1-t)/t)*(x^(-1/t)-z^(-1/(t^2+t))))/(z^(-t/(t+1)) - 1 + 
                            t^2*z^((1-t)/t)*(1-z^(-1/(t^2+t)))))*((z^(1/(t+1)) <= x)&(x < 1))
                            + 1*(x >= 1)}
invFx1z_0 <- function(u,z){z^(1-u)}
invFx1z_1 <- function(u,z){uniroot(function(x) Fx1z_1(x,z) - u, 
                                   interval = c(0.0001,0.999))$root}
invFx1z_t <- function(u,z,t){uniroot(function(x) Fx1z_t(x,z,t) - u, 
                             interval = c(0.0001,0.999), 
                             tol = .Machine$double.eps^0.25)$root}
#Y|Z
#teta = 0
Fy1z_0 <- function(y,z){((y/log(z))*(1-1/z))*((0 < y)&(y <= z)) + 
                        ((1/log(z))*(log(z)-log(y)+y-1))*((z<y)&(y<1)) + 1*(y>=1)}
#teta = 1
Fy1z_1 <- function(y,z) {(1/(log(z)+2/sqrt(z)-2))*((y/z - 1 + log(z) - 
                         log(y))*((z <= y)&(y <= sqrt(z))) + (2/sqrt(z) - 1/y - 1 + 
                         log(z)-log(y))*((sqrt(z) < y)&(y < 1))) + 1*(y>=1)}
#teta > 0 & teta !=1 
Fy1z_t <- function(y,z,t) {((t-t^2)/(z^(-t/(t+1))-1+t^2*z^(1/t - 
                            1)*(1-z^(-1/(t^2+t)))))*(((y - z^(1/t))/z + 
                            (z^((1-t)/t)-y^(1-t))/(1-t))*((z^(1/t) <= y)&(y <= z^(1/(t+1))))
                            + (z^(-t/(t+1))*(1+1/t) - (y^(-t))/t - (1/(1-t))*(y^(1-t) - 
                            t*z^(1/t - 1)))*((z^(1/(t+1)) < y)&(y < 1))) + 1*(y>=1)}

invFy1z_0 <- function(v,z){((z*v*log(z))/(z-1))*((0 < v)&(v <= (z-1)/log(z))) + 
                           (uniroot(function(y) (1/log(z))*(log(z)-log(y)+y-1) - v, 
                           interval = c(0.000001,0.999))$root)*(((z-1)/log(z)<v)&(v<1))}

invFy1z_1 <- function(v,z){uniroot(function(y) Fy1z_1(y,z) - v, 
                           interval = c(0.00000001,0.999))$root}

invFy1z_t <- function(v,z,t){uniroot(function(y) Fy1z_t(y,z,t) - v, 
                             interval = c(0.000000001,0.999), 
                             tol = .Machine$double.eps^0.25)$root}

#X|YZ
#teta = 0
invFx1yz_0 <- function(u,y,z) {((z/(1-u*(1-z)))*((0 <= u)&(u <= 1)))*(y <= z) + 
                               ((u*(1-y)+y)*((0 <= u)&(u <= 1)))*(z < y)}
#teta = 1
invFx1yz_1 <- function(u,y,z) {((z/(y-u*(y-z)))*((0 <= u)&(u <= 1)))*((y^2 <= z)&(z < y)) +
                               ((u*(1-y)+y)*((0 <= u)&(u <= 1)))*(z < y^2)}
#teta > 0 & teta != 1
invFx1yz_t <- function(u,y,z,t){((z/(y^{t} - u*(y^{t} - z)))*((0 <= u)&(u <= 1)))*(
                                 y^{1+t} <= z) + ((u*(1-y)+y)*((0 <= u)&(u <= 1)))*(
                                 z < y^{1+t})}
#X|Y
invFx1y <- function(u,y) {(y^{1-u})*((0 <= u)&(u <= 1))}

#Z|XY
invFz1xy <- function(w,x,y,t) qunif(w,0,x*y^{t})


# (X,Y|Z)

# Figura 3.16
# SI TETA = 0

# z = 0.25

m <- 4000
(selz25_0 <- 0.25)
v <- runif(m, 0, 1)
y1selz25_0 <- mapply(invFy1z_0, v, rep(selz25_0, m))
u <- runif(m, 0, 1)
x1selz25_0 <- mapply(invFx1yz_0, u, y1selz25_0, rep(selz25_0, m))
xy1selecz25_0 <- cbind(x1selz25_0, y1selz25_0)
summary(xy1selecz25_0)
matUV1z25_0 <- apply(xy1selecz25_0, 2, rank)/nrow(xy1selecz25_0)

## z = 0.75

(selz75_0 <- 0.75)
v <- runif(m, 0, 1)
y1selz75_0 <- mapply(invFy1z_0, v, rep(selz75_0, m))
u <- runif(m, 0, 1)
x1selz75_0 <- mapply(invFx1yz_0, u, y1selz75_0, rep(selz75_0, m))
xy1selecz75_0 <- cbind(x1selz75_0, y1selz75_0)
summary(xy1selecz75_0)
matUV1z75_0 <- apply(xy1selecz75_0, 2, rank)/nrow(xy1selecz75_0)

dev.new(); par(mfrow = c(2,2))
plot(matUV1z25_0, main = "Pseudo-observaciones de (X,Y|Z)", xlab = "u|z", ylab = "v|z")
plot(matUV1z75_0, main = "Pseudo-observaciones de (X,Y|Z)", xlab = "u|z", ylab = "v|z")
plot(xy1selecz25_0, main = "Observaciones de (X,Y|Z)", xlab = "x|z", ylab = "y|z", 
     xlim = c(0,1))
plot(xy1selecz75_0, main = "Observaciones de (X,Y|Z)", xlab = "x|z", ylab = "y|z",
     xlim = c(0,1))

# Figura 3.17

#SI TETA = 0.5
m <- 3000
t <- 0.5

#Z = 0.25

(selz25_05 <- 0.25)
v <- runif(m, 0, 1)
y1selz25_05 <- mapply(invFy1z_t, v, rep(selz25_05, m), rep(t,m))
u <- runif(m, 0, 1)
x1selz25_05 <- mapply(invFx1yz_t, u, y1selz25_05, rep(selz25_05, m), rep(t,m))
xy1selecz25_05 <- cbind(x1selz25_05, y1selz25_05)
summary(xy1selecz25_05)
matUV1z25_05 <- apply(xy1selecz25_05, 2, rank)/nrow(xy1selecz25_05)

#Z = 0.75

(selz75_05 <- 0.75)
v <- runif(m, 0, 1)
y1selz75_05 <- mapply(invFy1z_t, v, rep(selz75_05, m), rep(t,m))
u <- runif(m, 0, 1)
x1selz75_05 <- mapply(invFx1yz_t, u, y1selz75_05, rep(selz75_05, m), rep(t,m))
xy1selecz75_05 <- cbind(x1selz75_05, y1selz75_05)
summary(xy1selecz75_05)
matUV1z75_05 <- apply(xy1selecz75_05, 2, rank)/nrow(xy1selecz75_05)

dev.new(); par(mfrow = c(2,2))
plot(matUV1z25_05, main = "Pseudo-observaciones de (X,Y|Z)", xlab = "u|z", ylab = "v|z")
plot(matUV1z75_05, main = "Pseudo-observaciones de (X,Y|Z)", xlab = "u|z", ylab = "v|z")
plot(xy1selecz25_05, main = "Observaciones de (X,Y|Z)", xlab = "x|z", ylab = "y|z", 
     xlim = c(0,1))
plot(xy1selecz75_05, main = "Observaciones de (X,Y|Z)", xlab = "x|z", ylab = "y|z", 
     xlim = c(0,1))

#PARA (X,Z|Y)

#Figura 3.18
#Teta = 0

m <- 3000
t <- 0
#Y = 0.25
(ysel25_0 <- 0.25)
u <- runif(m)
x1ysel25_0 <- mapply(invFx1y, u, rep(ysel25_0, m))
w <- runif(m)
z1xy_sel25_0 <- mapply(invFz1xy, w, x1ysel25_0, rep(ysel25_0, m), rep(t, m))
xz1ysel25_0 <- cbind(x1ysel25_0, z1xy_sel25_0)
summary(xz1ysel25_0)
matUV1y25_0 <- apply(xz1ysel25_0, 2, rank)/nrow(xz1ysel25_0)

#Y = 0.75
(ysel75_0 <- 0.75)
u <- runif(m)
x1ysel75_0 <- mapply(invFx1y, u, rep(ysel75_0, m))
w <- runif(m)
z1xy_sel75_0 <- mapply(invFz1xy, w, x1ysel75_0, rep(ysel75_0, m), rep(t, m))
xz1ysel75_0 <- cbind(x1ysel75_0, z1xy_sel75_0)
summary(xz1ysel75_0)
matUV1y75_0 <- apply(xz1ysel75_0, 2, rank)/nrow(xz1ysel75_0)

dev.new(); par(mfrow = c(2,2))
plot(matUV1y25_0, main = "Pseudo-observaciones de (X,Z|Y)", xlab = "u|y", ylab = "w|y")
plot(matUV1y75_0, main = "Pseudo-observaciones de (X,Z|Y)", xlab = "u|y", ylab = "w|y")
plot(xz1ysel25_0, main = "Observaciones de (X,Z|Y)", xlab = "x|z", ylab = "y|z",
      xlim = c(0,1), ylim = c(0,1))
plot(xz1ysel75_0, main = "Observaciones de (X,Z|Y)", xlab = "x|z", ylab = "y|z",
      xlim = c(0,1), ylim = c(0,1))


# Figura 3.19

#Teta = 0.5
m <- 3000
t <- 0.5
#Y = 0.25
(ysel25_05 <- 0.25)
u <- runif(m)
x1ysel25_05 <- mapply(invFx1y, u, rep(ysel25_05, m))
w <- runif(m)
z1xy_sel25_05 <- mapply(invFz1xy, w, x1ysel25_05, rep(ysel25_05, m), rep(t, m))
xz1ysel25_05 <- cbind(x1ysel25_05, z1xy_sel25_05)
summary(xz1ysel25_05)
matUV1y25_05 <- apply(xz1ysel25_05, 2, rank)/nrow(xz1ysel25_05)

#Y = 0.75
(ysel75_05 <- 0.75)
u <- runif(m)
x1ysel75_05 <- mapply(invFx1y, u, rep(ysel75_05, m))
w <- runif(m)
z1xy_sel75_05 <- mapply(invFz1xy, w, x1ysel75_05, rep(ysel75_05, m), rep(t, m))
xz1ysel75_05 <- cbind(x1ysel75_05, z1xy_sel75_05)
summary(xz1ysel75_05)
matUV1y75_05 <- apply(xz1ysel75_05, 2, rank)/nrow(xz1ysel75_05)

dev.new(); par(mfrow = c(2,2))
plot(matUV1y25_05, main = "Pseudo-observaciones de (X,Z|Y)", xlab = "u|y", ylab = "w|y")
plot(matUV1y75_05, main = "Pseudo-observaciones de (X,Z|Y)", xlab = "u|y", ylab = "w|y")
plot(xz1ysel25_05, main = "Observaciones de (X,Z|Y)", xlab = "x|z", ylab = "y|z",
      xlim = c(0,1), ylim = c(0,1))
plot(xz1ysel75_05, main = "Observaciones de (X,Z|Y)", xlab = "x|z", ylab = "y|z",
      xlim = c(0,1), ylim = c(0,1))

# (Y,Z|X)

#X
invFx <- function(u) qunif(u)
#Y|X
invFy1x <- function(v,x) qunif(v,0,x)
#Z|XY
invFz1xy <- function(w,x,y,t) qunif(w,0,x*y^{t})

# Figura 3.20

#PARA TETA = 0
m <- 3000
t <- 0
# Z = 0.25
(selx25_0 <- 0.25 )
v <- runif(m)
y1selx25_0 <- mapply(invFy1x, v, rep(selx25_0, m))
w <- runif(m)
z1selx25_0 <- mapply(invFz1xy, w, rep(selx25_0, m), y1selx25_0, rep(t, m))
yz1selx25_0 <- cbind(y1selx25_0, z1selx25_0)
summary(yz1selx25_0)
matUV1x25_0 <- apply(yz1selx25_0, 2, rank)/nrow(yz1selx25_0)

# Z = 0.75
(selx75_0 <- 0.75)
v <- runif(m)
y1selx75_0 <- mapply(invFy1x, v, rep(selx75_0, m))
w <- runif(m)
z1selx75_0 <- mapply(invFz1xy, w, rep(selx75_0, m), y1selx75_0, rep(t, m))
yz1selx75_0 <- cbind(y1selx75_0, z1selx75_0)
summary(yz1selx75_0)
matUV1x75_0 <- apply(yz1selx75_0, 2, rank)/nrow(yz1selx75_0)

dev.new(); par(mfrow = c(2,2))
plot(matUV1x25_0, main = "Pseudo-observaciones de (Y,Z|X)", xlab = "v|x", ylab = "w|x")
plot(matUV1x75_0, main = "Pseudo-observaciones de (Y,Z|X)", xlab = "v|x", ylab = "w|x")
plot(yz1selx25_0, main = "Observaciones de (Y,Z|X)", xlab = "x|z", ylab = "y|z",
      xlim = c(0,1), ylim = c(0,1))
plot(yz1selx75_0, main = "Observaciones de (Y,Z|X)", xlab = "x|z", ylab = "y|z",
      xlim = c(0,1), ylim = c(0,1))

# Figura 3.21

##PARA TETA = 0.5
m <- 3000
t <- 0.5

#Z = 0.25
(selx25_05 <- 0.25)
v <- runif(m)
y1selx25_05 <- mapply(invFy1x, v, rep(selx25_05, m))
w <- runif(m)
z1selx25_05 <- mapply(invFz1xy, w, rep(selx25_05, m), y1selx25_05, rep(t, m))
yz1selx25_05 <- cbind(y1selx25_05, z1selx25_05)
summary(yz1selx25_05)
matUV1x25_05 <- apply(yz1selx25_05, 2, rank)/nrow(yz1selx25_05)

#Z = 0.75
(selx75_05 <- 0.75)
v <- runif(m)
y1selx75_05 <- mapply(invFy1x, v, rep(selx75_05, m))
w <- runif(m)
z1selx75_05 <- mapply(invFz1xy, w, rep(selx75_05, m), y1selx75_05, rep(t, m))
yz1selx75_05 <- cbind(y1selx75_05, z1selx75_05)
summary(yz1selx75_05)
matUV1x75_05 <- apply(yz1selx75_05, 2, rank)/nrow(yz1selx75_05)

dev.new(); par(mfrow = c(2,2))
plot(matUV1x25_05, main = "Pseudo-observaciones de (Y,Z|X)", xlab = "v|x", ylab = "w|x")
plot(yz1selx25_05, main = "Observaciones de (Y,Z|X)", xlab = "x|z", ylab = "y|z",
      xlim = c(0,1), ylim = c(0,1))
plot(matUV1x75_05, main = "Pseudo-observaciones de (Y,Z|X)", xlab = "v|x", ylab = "w|x")
plot(yz1selx75_05, main = "Observaciones de (Y,Z|X)", xlab = "x|z", ylab = "y|z",
      xlim = c(0,1), ylim = c(0,1))

