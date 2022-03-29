library(TwoCop)
library(subcopem2D)
dif.copula <- function(g1, g2,m){
  copempg1 <- subcopemc(g1,50)$matrix
  copempg2 <- subcopemc(g2,50)$matrix
  T <- sum((copempg1 - copempg2)^{2})
  n <- nrow(g1)
  aux.index <- seq(1, 2*n, by = 1)
  aux <- rbind(g1,g2)
  T_est <- numeric(m)
  for(i in 1:m){
              mg1.index <- sample(aux.index, n)
              mg1 <- aux[mg1.index, ]
              mg2 <- aux[-mg1.index,]
              copempmg1 <- subcopemc(mg1,51)$matrix
              copempmg2 <- subcopemc(mg2,51)$matrix
              T_est[i] <- sum((copempmg1 - copempmg2)^{2})
              }
  pvalue <- sum(T <= T_est)/m
 return(pvalue)
}

#FUNCIONES NECESARIAS.
#Y|X
Fy1x <- function(y,x) punif(y, 0, x) 

#Z|X
#teta = 0
Fz1x_0 <- function(z,x) punif(z, 0, x)
#teta = 1
Fz1x_1 <- function(z,x) {( (2*z*log(x)-z*log(z)+z)/(x^2))*((0<z)&(z<x^2)) + 1*(x^2 <= z)}
#teta > 0 & teta !=1
Fz1x_t <- function(z,x,t) {((1/(1-t))*(z/x^(t+1) - (t*z^(1/t))/x^(1/t +1)))*((0 <= z)&(z < x^(t+1))) + 1*(x^(t+1)<=z)} 

#X|Y
Fx1y <- function(x,y) {(1-log(x)/log(y))*((y <= x)&(x < 1)) + 1*(x >= 1)}

#Z|Y
Fz1y <- function(z,y,t) {(-(z/log(y))*(1/y^{t+1} - 1/y^{t}))*((0 <= z)&(z <= y^{1+t})) + (z/(y^{t}*log(y)) - 1/log(y) - log(z)/log(y) + t+1)*((y^{t+1} < z)&(z < y^{t})) + 1*(z >= y^{t})}

#X|Z
#teta = 0
Fx1z_0 <- function(x,z) {(1-log(x)/log(z))*((z <= x)&(x < 1)) + 1*(x >= 1)}
#teta = 1
Fx1z_1 <- function(x,z) {((1/(log(z)+2/sqrt(z)-2))*(2/sqrt(z) + (log(z)-2*log(x)-2)/x))*((sqrt(z) <= x)&(x < 1)) + 1*(x >= 1)}
#teta > 0 & teta != 1
Fx1z_t <- function(x,z,t){((z^(-t/(t+1)) - x^(-t) + t^2*z^((1-t)/t)*(x^(-1/t)-z^(-1/(t^2+t))))/(z^(-t/(t+1)) - 1 + t^2*z^((1-t)/t)*(1-z^(-1/(t^2+t)))))*((z^(1/(t+1)) <= x)&(x < 1)) + 1*(x >= 1)}


#Y|Z
#teta = 0
Fy1z_0 <- function(y,z){((y/log(z))*(1-1/z))*((0 < y)&(y <= z)) + ((1/log(z))*(log(z)-log(y)+y-1))*((z<y)&(y<1)) + 1*(y>=1)}
#teta = 1
Fy1z_1 <- function(y,z) {(1/(log(z)+2/sqrt(z)-2))*((y/z - 1 + log(z) - log(y))*((z <= y)&(y <= sqrt(z))) + (2/sqrt(z) - 1/y - 1 + log(z)-log(y))*((sqrt(z) < y)&(y < 1))) + 1*(y>=1)}
#teta > 0 & teta !=1 
Fy1z_t <- function(y,z,t) {((t-t^2)/(z^(-t/(t+1))-1+t^2*z^(1/t -1)*(1-z^(-1/(t^2+t)))))*(((y - z^(1/t))/z + (z^((1-t)/t)-y^(1-t))/(1-t))*((z^(1/t) <= y)&(y <= z^(1/(t+1)))) + (z^(-t/(t+1))*(1+1/t) - (y^(-t))/t - (1/(1-t))*(y^(1-t)-t*z^(1/t - 1)))*((z^(1/(t+1)) < y)&(y < 1))) + 1*(y>=1)}

invFx <- function(u)qunif(u,0,1)
invFylx <- function(v,x)qunif(v,0,x)
#Teta = 0
invFz1xy_0 <- function(w,x,y) qunif(w, 0, x)
#Teta = 1
invFz1xy_1 <- function(w,x,y) qunif(w, 0, x*y)
#Teta>0 y theta != 1
invFz1xy_t <- function(w,x,y,t) qunif(w, 0, x*y^t)

########## PARA TETA = 0.
#####SIMULACIÓN DE (X,Y,Z)
m <- 5000
u <- runif(m,0,1)
x_0 <- invFx(u)
v <- runif(m,0,1)
y_0 <- mapply(invFylx, v, x_0)
w <- runif(m,0,1)
z_0 <- mapply(invFz1xy_0, w, x_0, y_0)
ma.obs_0 <- cbind(x_0,y_0,z_0)
head(ma.obs_0)

###(Y,Z|X)
(mediana.x_0 <- median(x_0))
vi.dx_0 <- mapply(Fy1x, y_0, x_0)
wi.dx_0 <- mapply(Fz1x_0, z_0, x_0)
viwi.dx_0 <- cbind(vi.dx_0,wi.dx_0)
head(viwi.dx_0)
aux1.dx_0 <- cbind(x_0, viwi.dx_0)
head(aux1.dx_0)
part_AM.aux.dx_0 <- aux1.dx_0[aux1.dx_0[,1] <= mediana.x_0, ]  
part_DM.aux.dx_0 <- aux1.dx_0[aux1.dx_0[,1] >  mediana.x_0, ]
part_AM.dx_0 <- part_AM.aux.dx_0[,2:3] 
part_DM.dx_0 <- part_DM.aux.dx_0[,2:3] 
dev.new(); par(mfrow = c(1,3))
plot(part_AM.dx_0, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (Y,Z|X)", 
     sub = "X <= Mediana{X}")
plot(part_DM.dx_0, pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (Y,Z|X)",
     sub = "X > Mediana{X}")
plot(part_AM.dx_0, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (Y,Z|X)")
points(part_DM.dx_0, pch = 20, lwd = 2)
system.time(
pvalue.dx_0 <- dif.copula(part_AM.dx_0, part_DM.dx_0,10000)) 

###(X,Z|Y)
t <- 0
(mediana.y_0 <- median(y_0))
ui.dy_0 <- mapply(Fx1y, x_0, y_0)
wi.dy_0 <- mapply(Fz1y, z_0, y_0, rep(t, m))
uiwi.dy_0 <- cbind(ui.dy_0, wi.dy_0)
head(uiwi.dy_0)
aux1.dy_0 <- cbind(y_0, uiwi.dy_0)
head(aux1.dy_0)
part_AM.aux.dy_0 <- aux1.dy_0[aux1.dy_0[,1] <= mediana.y_0, ]  
part_DM.aux.dy_0 <- aux1.dy_0[aux1.dy_0[,1] >  mediana.y_0, ]
part_AM.dy_0 <- part_AM.aux.dy_0[,2:3] 
part_DM.dy_0 <- part_DM.aux.dy_0[,2:3] 
dev.new(); par(mfrow = c(1,3))
plot(part_AM.dy_0, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Z|Y)", 
     sub = "Y <= Mediana{Y}")
plot(part_DM.dy_0, pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Z|Y)",
     sub = "Y > Mediana{Y}")
plot(part_AM.dy_0, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Z|Y)")
points(part_DM.dy_0, pch = 20, lwd = 2)
pvalue.dy_0 <- dif.copula(part_AM.dy_0, part_DM.dy_0,10000)
 
###(X,Y|Z)
(mediana.z_0 <- median(z_0))
ui.dz_0 <- mapply(Fx1z_0, x_0, z_0)
vi.dz_0 <- mapply(Fy1z_0, y_0, z_0)
uivi.dz_0 <- cbind(ui.dz_0, vi.dz_0)
head(uivi.dz_0)
aux1.dz_0 <- cbind(z_0, uivi.dz_0)
head(aux1.dz_0)
part_AM.aux.dz_0 <- aux1.dz_0[aux1.dz_0[,1] <= mediana.z_0, ]  
part_DM.aux.dz_0 <- aux1.dz_0[aux1.dz_0[,1] >  mediana.z_0, ]
part_AM.dz_0 <- part_AM.aux.dz_0[,2:3] 
part_DM.dz_0 <- part_DM.aux.dz_0[,2:3] 
dev.new(); par(mfrow = c(1,3))
plot(part_AM.dz_0, col = "blue", pch = 20, lwd =2,
     main = "Pseudo-observaciones de la cópula de (X,Y|Z)", 
     sub = "Z <= Mediana{Z}")
plot(part_DM.dz_0, pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Y|Z)",
     sub = "Z > Mediana{Z}")
plot(part_AM.dz_0, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Y|Z)")
points(part_DM.dz_0, pch = 20, lwd = 2)
pvalue.dz_0 <- dif.copula(part_AM.dz_0, part_DM.dz_0,10000)
 

########## PARA TETA = 0.1
#####SIMULACIÓN DE (X,Y,Z)
m <- 1000
t <- 0.1
u <- runif(m,0,1)
x_01 <- invFx(u)
v <- runif(m,0,1)
y_01 <- mapply(invFylx, v, x_01)
w <- runif(m,0,1)
z_01 <- mapply(invFz1xy_t, w, x_01, y_01, rep(t, m))
ma.obs_01 <- cbind(x_01,y_01,z_01)
head(ma.obs_01)

###(Y,Z|X)
(mediana.x_01 <- median(x_01))
vi.dx_01 <- mapply(Fy1x, y_01, x_01)
wi.dx_01 <- mapply(Fz1x_t, z_01, x_01, rep(t,m))
viwi.dx_01 <- cbind(vi.dx_01,wi.dx_01)
head(viwi.dx_01)
aux1.dx_01 <- cbind(x_01, viwi.dx_01)
head(aux1.dx_01)
part_AM.aux.dx_01 <- aux1.dx_01[aux1.dx_01[,1] <= mediana.x_01, ]  
part_DM.aux.dx_01 <- aux1.dx_01[aux1.dx_01[,1] >  mediana.x_01, ]
part_AM.dx_01 <- part_AM.aux.dx_01[,2:3] 
part_DM.dx_01 <- part_DM.aux.dx_01[,2:3] 
dev.new(); par(mfrow = c(1,3))
plot(part_AM.dx_01, col = "blue", pch = 20, lwd = 2, 
     main = "Pseudo-observaciones de la cópula de (Y,Z|X)", 
     sub = "X <= Mediana{X}")
plot(part_DM.dx_01, pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (Y,Z|X)",
     sub = "X > Mediana{X}")
plot(part_AM.dx_01, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (Y,Z|X)")
points(part_DM.dx_01, pch = 20, lwd = 2)
dif.copula(part_AM.dx_01, part_DM.dx_01,10000)

###(X,Z|Y)
t <- 0.1
(mediana.y_01 <- median(y_01))
ui.dy_01 <- mapply(Fx1y, x_01, y_01)
wi.dy_01 <- mapply(Fz1y, z_01, y_01, rep(t, m))
uiwi.dy_01 <- cbind(ui.dy_01, wi.dy_01)
head(uiwi.dy_01)
aux1.dy_01 <- cbind(y_01, uiwi.dy_01)
head(aux1.dy_01)
part_AM.aux.dy_01 <- aux1.dy_01[aux1.dy_01[,1] <= mediana.y_01, ]  
part_DM.aux.dy_01 <- aux1.dy_01[aux1.dy_01[,1] >  mediana.y_01, ]
part_AM.dy_01 <- part_AM.aux.dy_01[,2:3] 
part_DM.dy_01 <- part_DM.aux.dy_01[,2:3] 
dev.new(); par(mfrow = c(1,3))
plot(part_AM.dy_01, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Z|Y)", 
     sub = "Y <= Mediana{Y}")
plot(part_DM.dy_01, pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Z|Y)",
     sub = "Y > Mediana{Y}")
plot(part_AM.dy_01, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Z|Y)")
points(part_DM.dy_01, pch = 20, lwd = 2)
dif.copula(part_AM.dy_01, part_DM.dy_01,10000)

###(X,Y|Z)
(mediana.z_01 <- median(z_01))
ui.dz_01 <- mapply(Fx1z_t, x_01, z_01, rep(t,m))
vi.dz_01 <- mapply(Fy1z_t, y_01, z_01, rep(t,m))
uivi.dz_01 <- cbind(ui.dz_01, vi.dz_01)
head(uivi.dz_01)
aux1.dz_01 <- cbind(z_01, uivi.dz_01)
head(aux1.dz_01)
part_AM.aux.dz_01 <- aux1.dz_01[aux1.dz_01[,1] <= mediana.z_01, ]  
part_DM.aux.dz_01 <- aux1.dz_01[aux1.dz_01[,1] >  mediana.z_01, ]
part_AM.dz_01 <- part_AM.aux.dz_01[,2:3] 
part_DM.dz_01 <- part_DM.aux.dz_01[,2:3] 
dev.new(); par(mfrow = c(1,3))
plot(part_AM.dz_01, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Y|Z)", 
     sub = "Z <= Mediana{Z}")
plot(part_DM.dz_01, pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Y|Z)",
     sub = "Z > Mediana{Z}")
plot(part_AM.dz_01, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Y|Z)")
points(part_DM.dz_01, pch = 20, lwd = 2)
dif.copula(part_AM.dz_01, part_DM.dz_01,10000)


########## PARA TETA = 0.5
#####SIMULACIÓN DE (X,Y,Z)
m <- 5000
t <- 0.5
u <- runif(m,0,1)
x_05 <- invFx(u)
v <- runif(m,0,1)
y_05 <- mapply(invFylx, v, x_05)
w <- runif(m,0,1)
z_05 <- mapply(invFz1xy_t, w, x_05, y_05, rep(t, m))
ma.obs_05 <- cbind(x_05,y_05,z_05)
head(ma.obs_05)

###(Y,Z|X)
(mediana.x_05 <- median(x_05))
vi.dx_05 <- mapply(Fy1x, y_05, x_05)
wi.dx_05 <- mapply(Fz1x_t, z_05, x_05, rep(t,m))
viwi.dx_05 <- cbind(vi.dx_05,wi.dx_05)
head(viwi.dx_05)
aux1.dx_05 <- cbind(x_05, viwi.dx_05)
head(aux1.dx_05)
part_AM.aux.dx_05 <- aux1.dx_05[aux1.dx_05[,1] <= mediana.x_05, ]  
part_DM.aux.dx_05 <- aux1.dx_05[aux1.dx_05[,1] >  mediana.x_05, ]
part_AM.dx_05 <- part_AM.aux.dx_05[,2:3] 
part_DM.dx_05 <- part_DM.aux.dx_05[,2:3] 
dev.new(); par(mfrow = c(1,3))
plot(part_AM.dx_05, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (Y,Z|X)", 
     sub = "X <= Mediana{X}")
plot(part_DM.dx_05, pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (Y,Z|X)",
     sub = "X > Mediana{X}")
plot(part_AM.dx_05, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (Y,Z|X)")
points(part_DM.dx_05, pch = 20, lwd = 2)
dif.copula(part_AM.dx_05, part_DM.dx_05,10000)
 
###(X,Z|Y)
t <- 0.5
(mediana.y_05 <- median(y_05))
ui.dy_05 <- mapply(Fx1y, x_05, y_05)
wi.dy_05 <- mapply(Fz1y, z_05, y_05, rep(t, m))
uiwi.dy_05 <- cbind(ui.dy_05, wi.dy_05)
head(uiwi.dy_05)
aux1.dy_05 <- cbind(y_05, uiwi.dy_05)
head(aux1.dy_05)
part_AM.aux.dy_05 <- aux1.dy_05[aux1.dy_05[,1] <= mediana.y_05, ]  
part_DM.aux.dy_05 <- aux1.dy_05[aux1.dy_05[,1] >  mediana.y_05, ]
part_AM.dy_05 <- part_AM.aux.dy_05[,2:3] 
part_DM.dy_05 <- part_DM.aux.dy_05[,2:3] 
dev.new(); par(mfrow = c(1,3))
plot(part_AM.dy_05, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Z|Y)", 
     sub = "Y <= Mediana{Y}")
plot(part_DM.dy_05, pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Z|Y)",
     sub = "Y > Mediana{Y}")
plot(part_AM.dy_05, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Z|Y)")
points(part_DM.dy_05, pch = 20, lwd = 2)
dif.copula(part_AM.dy_05, part_DM.dy_05,10000)

###(X,Y|Z)
(mediana.z_05 <- median(z_05))
ui.dz_05 <- mapply(Fx1z_t, x_05, z_05, rep(t,m))
vi.dz_05 <- mapply(Fy1z_t, y_05, z_05, rep(t,m))
uivi.dz_05 <- cbind(ui.dz_05, vi.dz_05)
head(uivi.dz_05)
aux1.dz_05 <- cbind(z_05, uivi.dz_05)
head(aux1.dz_05)
part_AM.aux.dz_05 <- aux1.dz_05[aux1.dz_05[,1] <= mediana.z_05, ]  
part_DM.aux.dz_05 <- aux1.dz_05[aux1.dz_05[,1] >  mediana.z_05, ]
part_AM.dz_05 <- part_AM.aux.dz_05[,2:3] 
part_DM.dz_05 <- part_DM.aux.dz_05[,2:3] 
dev.new(); par(mfrow = c(1,3))
plot(part_AM.dz_05, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Y|Z)", 
     sub = "Z <= Mediana{Z}")
plot(part_DM.dz_05, pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Y|Z)",
     sub = "Z > Mediana{Z}")
plot(part_AM.dz_05, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Y|Z)")
points(part_DM.dz_05, pch = 20, lwd = 2)
dif.copula(part_AM.dz_05, part_DM.dz_05,10000)


########## PARA TETA = 0.7
#####SIMULACIÓN DE (X,Y,Z)
m <- 5000
t <- 0.7
u <- runif(m,0,1)
x_07 <- invFx(u)
v <- runif(m,0,1)
y_07 <- mapply(invFylx, v, x_07)
w <- runif(m,0,1)
z_07 <- mapply(invFz1xy_t, w, x_07, y_07, rep(t, m))
ma.obs_07 <- cbind(x_07,y_07,z_07)
head(ma.obs_07)

###(Y,Z|X)
(mediana.x_07 <- median(x_07))
vi.dx_07 <- mapply(Fy1x, y_07, x_07)
wi.dx_07 <- mapply(Fz1x_t, z_07, x_07, rep(t,m))
viwi.dx_07 <- cbind(vi.dx_07,wi.dx_07)
head(viwi.dx_07)
aux1.dx_07 <- cbind(x_07, viwi.dx_07)
head(aux1.dx_07)
part_AM.aux.dx_07 <- aux1.dx_07[aux1.dx_07[,1] <= mediana.x_07, ]  
part_DM.aux.dx_07 <- aux1.dx_07[aux1.dx_07[,1] >  mediana.x_07, ]
part_AM.dx_07 <- part_AM.aux.dx_07[,2:3] 
part_DM.dx_07 <- part_DM.aux.dx_07[,2:3] 
dev.new(); par(mfrow = c(1,3))
plot(part_AM.dx_07, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (Y,Z|X)", 
     sub = "X <= Mediana{X}")
plot(part_DM.dx_07, pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (Y,Z|X)",
     sub = "X > Mediana{X}")
plot(part_AM.dx_07, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (Y,Z|X)")
points(part_DM.dx_07, pch = 20, lwd = 2)
dif.copula(part_AM.dx_07, part_DM.dx_07,10000)
  
###(X,Z|Y)
t <- 0.7
(mediana.y_07 <- median(y_07))
ui.dy_07 <- mapply(Fx1y, x_07, y_07)
wi.dy_07 <- mapply(Fz1y, z_07, y_07, rep(t, m))
uiwi.dy_07 <- cbind(ui.dy_07, wi.dy_07)
head(uiwi.dy_07)
aux1.dy_07 <- cbind(y_07, uiwi.dy_07)
head(aux1.dy_07)
part_AM.aux.dy_07 <- aux1.dy_07[aux1.dy_07[,1] <= mediana.y_07, ]  
part_DM.aux.dy_07 <- aux1.dy_07[aux1.dy_07[,1] >  mediana.y_07, ]
part_AM.dy_07 <- part_AM.aux.dy_07[,2:3] 
part_DM.dy_07 <- part_DM.aux.dy_07[,2:3] 
dev.new(); par(mfrow = c(1,3))
plot(part_AM.dy_07, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Z|Y)", 
     sub = "Y <= Mediana{Y}")
plot(part_DM.dy_07, pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Z|Y)",
     sub = "Y > Mediana{Y}")
plot(part_AM.dy_07, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Z|Y)")
points(part_DM.dy_07, pch = 20, lwd = 2)
dif.copula(part_AM.dy_07, part_DM.dy_07,10000)
  

###(X,Y|Z)
(mediana.z_07 <- median(z_07))
ui.dz_07 <- mapply(Fx1z_t, x_07, z_07, rep(t,m))
vi.dz_07 <- mapply(Fy1z_t, y_07, z_07, rep(t,m))
uivi.dz_07 <- cbind(ui.dz_07, vi.dz_07)
head(uivi.dz_07)
aux1.dz_07 <- cbind(z_07, uivi.dz_07)
head(aux1.dz_07)
part_AM.aux.dz_07 <- aux1.dz_07[aux1.dz_07[,1] <= mediana.z_07, ]  
part_DM.aux.dz_07 <- aux1.dz_07[aux1.dz_07[,1] >  mediana.z_07, ]
part_AM.dz_07 <- part_AM.aux.dz_07[,2:3] 
part_DM.dz_07 <- part_DM.aux.dz_07[,2:3] 
dev.new(); par(mfrow = c(1,3))
plot(part_AM.dz_07, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Y|Z)", 
     sub = "Z <= Mediana{Z}")
plot(part_DM.dz_07, pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Y|Z)",
     sub = "Z > Mediana{Z}")
plot(part_AM.dz_07, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Y|Z)")
points(part_DM.dz_07, pch = 20, lwd = 2)
dif.copula(part_AM.dz_07, part_DM.dz_07,10000)
  


########## PARA TETA = 1
#####SIMULACIÓN DE (X,Y,Z)
m <- 5000
u <- runif(m,0,1)
x_1 <- invFx(u)
v <- runif(m,0,1)
y_1 <- mapply(invFylx, v, x_1)
w <- runif(m,0,1)
z_1 <- mapply(invFz1xy_1, w, x_1, y_1)
ma.obs_1 <- cbind(x_1,y_1,z_1)
head(ma.obs_1)


###(Y,Z|X)
(mediana.x_1 <- median(x_1))
vi.dx_1 <- mapply(Fy1x, y_1, x_1)
wi.dx_1 <- mapply(Fz1x_1, z_1, x_1)
viwi.dx_1 <- cbind(vi.dx_1,wi.dx_1)
head(viwi.dx_1)
aux1.dx_1 <- cbind(x_1, viwi.dx_1)
head(aux1.dx_1)
part_AM.aux.dx_1 <- aux1.dx_1[aux1.dx_1[,1] <= mediana.x_1, ]  
part_DM.aux.dx_1 <- aux1.dx_1[aux1.dx_1[,1] >  mediana.x_1, ]
part_AM.dx_1 <- part_AM.aux.dx_1[,2:3] 
part_DM.dx_1 <- part_DM.aux.dx_1[,2:3] 
dev.new(); par(mfrow = c(1,3))
plot(part_AM.dx_1, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (Y,Z|X)", 
     sub = "X <= Mediana{X}")
plot(part_DM.dx_1, pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (Y,Z|X)",
     sub = "X > Mediana{X}")
plot(part_AM.dx_1, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (Y,Z|X)")
points(part_DM.dx_1, pch = 20, lwd = 2)
dif.copula(part_AM.dx_1, part_DM.dx_1,10000)


###(X,Z|Y)
t <- 1
(mediana.y_1 <- median(y_1))
ui.dy_1 <- mapply(Fx1y, x_1, y_1)
wi.dy_1 <- mapply(Fz1y, z_1, y_1, rep(t, m))
uiwi.dy_1 <- cbind(ui.dy_1, wi.dy_1)
head(uiwi.dy_1)
aux1.dy_1 <- cbind(y_1, uiwi.dy_1)
head(aux1.dy_1)
part_AM.aux.dy_1 <- aux1.dy_1[aux1.dy_1[,1] <= mediana.y_1, ]  
part_DM.aux.dy_1 <- aux1.dy_1[aux1.dy_1[,1] >  mediana.y_1, ]
part_AM.dy_1 <- part_AM.aux.dy_1[,2:3] 
part_DM.dy_1 <- part_DM.aux.dy_1[,2:3] 
dev.new(); par(mfrow = c(1,3))
plot(part_AM.dy_1, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Z|Y)", 
     sub = "Y <= Mediana{Y}")
plot(part_DM.dy_1, pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Z|Y)",
     sub = "Y > Mediana{Y}")
plot(part_AM.dy_1, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Z|Y)")
points(part_DM.dy_1, pch = 20, lwd = 2)
dif.copula(part_AM.dy_1, part_DM.dy_1,10000)

 
###(X,Y|Z)
(mediana.z_1 <- median(z_1))
ui.dz_1 <- mapply(Fx1z_1, x_1, z_1)
vi.dz_1 <- mapply(Fy1z_1, y_1, z_1)
uivi.dz_1 <- cbind(ui.dz_1, vi.dz_1)
head(uivi.dz_1)
aux1.dz_1 <- cbind(z_1, uivi.dz_1)
head(aux1.dz_1)
part_AM.aux.dz_1 <- aux1.dz_1[aux1.dz_1[,1] <= mediana.z_1, ]  
part_DM.aux.dz_1 <- aux1.dz_1[aux1.dz_1[,1] >  mediana.z_1, ]
part_AM.dz_1 <- part_AM.aux.dz_1[,2:3] 
part_DM.dz_1 <- part_DM.aux.dz_1[,2:3] 
dev.new(); par(mfrow = c(1,3))
plot(part_AM.dz_1, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Y|Z)", 
     sub = "Z <= Mediana{Z}")
plot(part_DM.dz_1, pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Y|Z)",
     sub = "Z > Mediana{Z}")
plot(part_AM.dz_1, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Y|Z)")
points(part_DM.dz_1, pch = 20, lwd = 2)
dif.copula(part_AM.dz_1, part_DM.dz_1,10000)
 

########## PARA TETA = 1.1
#####SIMULACIÓN DE (X,Y,Z)
m <- 5000
t <- 1.1
u <- runif(m,0,1)
x_11 <- invFx(u)
v <- runif(m,0,1)
y_11 <- mapply(invFylx, v, x_11)
w <- runif(m,0,1)
z_11 <- mapply(invFz1xy_t, w, x_11, y_11, rep(t, m))
ma.obs_11 <- cbind(x_11,y_11,z_11)
head(ma.obs_11)

###(Y,Z|X)
(mediana.x_11 <- median(x_11))
vi.dx_11 <- mapply(Fy1x, y_11, x_11)
wi.dx_11 <- mapply(Fz1x_t, z_11, x_11, rep(t,m))
viwi.dx_11 <- cbind(vi.dx_11,wi.dx_11)
head(viwi.dx_11)
aux1.dx_11 <- cbind(x_11, viwi.dx_11)
head(aux1.dx_11)
part_AM.aux.dx_11 <- aux1.dx_11[aux1.dx_11[,1] <= mediana.x_11, ]  
part_DM.aux.dx_11 <- aux1.dx_11[aux1.dx_11[,1] >  mediana.x_11, ]
part_AM.dx_11 <- part_AM.aux.dx_11[,2:3] 
part_DM.dx_11 <- part_DM.aux.dx_11[,2:3] 
dev.new(); par(mfrow = c(1,3))
plot(part_AM.dx_11, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (Y,Z|X)", 
     sub = "X <= Mediana{X}")
plot(part_DM.dx_11, pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (Y,Z|X)",
     sub = "X > Mediana{X}")
plot(part_AM.dx_11, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (Y,Z|X)")
points(part_DM.dx_11, pch = 20, lwd = 2)
dif.copula(part_AM.dx_11, part_DM.dx_11,10000)


###(X,Z|Y)
(mediana.y_11 <- median(y_11))
ui.dy_11 <- mapply(Fx1y, x_11, y_11)
wi.dy_11 <- mapply(Fz1y, z_11, y_11, rep(t, m))
uiwi.dy_11 <- cbind(ui.dy_11, wi.dy_11)
head(uiwi.dy_11)
aux1.dy_11 <- cbind(y_11, uiwi.dy_11)
head(aux1.dy_11)
part_AM.aux.dy_11 <- aux1.dy_11[aux1.dy_11[,1] <= mediana.y_11, ]  
part_DM.aux.dy_11 <- aux1.dy_11[aux1.dy_11[,1] >  mediana.y_11, ]
part_AM.dy_11 <- part_AM.aux.dy_11[,2:3] 
part_DM.dy_11 <- part_DM.aux.dy_11[,2:3] 
dev.new(); par(mfrow = c(1,3))
plot(part_AM.dy_11, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Z|Y)", 
     sub = "Y <= Mediana{Y}")
plot(part_DM.dy_11, pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Z|Y)",
     sub = "Y > Mediana{Y}")
plot(part_AM.dy_11, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Z|Y)")
points(part_DM.dy_11, pch = 20, lwd = 2)
dif.copula(part_AM.dy_11, part_DM.dy_11,10000)


###(X,Y|Z)
(mediana.z_11 <- median(z_11))
ui.dz_11 <- mapply(Fx1z_t, x_11, z_11, rep(t,m))
vi.dz_11 <- mapply(Fy1z_t, y_11, z_11, rep(t,m))
uivi.dz_11 <- cbind(ui.dz_11, vi.dz_11)
head(uivi.dz_11)
aux1.dz_11 <- cbind(z_11, uivi.dz_11)
head(aux1.dz_11)
part_AM.aux.dz_11 <- aux1.dz_11[aux1.dz_11[,1] <= mediana.z_11, ]  
part_DM.aux.dz_11 <- aux1.dz_11[aux1.dz_11[,1] >  mediana.z_11, ]
part_AM.dz_11 <- part_AM.aux.dz_11[,2:3] 
part_DM.dz_11 <- part_DM.aux.dz_11[,2:3] 
dev.new(); par(mfrow = c(1,3))
plot(part_AM.dz_11, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Y|Z)", 
     sub = "Z <= Mediana{Z}")
plot(part_DM.dz_11, pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Y|Z)",
     sub = "Z > Mediana{Z}")
plot(part_AM.dz_11, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Y|Z)")
points(part_DM.dz_11, pch = 20, lwd = 2)
dif.copula(part_AM.dz_11, part_DM.dz_11,10000)
 

########## PARA TETA = 1.5
#####SIMULACIÓN DE (X,Y,Z)

m <- 5000
t <- 1.5
u <- runif(m,0,1)
x_15 <- invFx(u)
v <- runif(m,0,1)
y_15 <- mapply(invFylx, v, x_15)
w <- runif(m,0,1)
z_15 <- mapply(invFz1xy_t, w, x_15, y_15, rep(t, m))
ma.obs_15 <- cbind(x_15,y_15,z_15)
head(ma.obs_15)

###(Y,Z|X)
(mediana.x_15 <- median(x_15))
vi.dx_15 <- mapply(Fy1x, y_15, x_15)
wi.dx_15 <- mapply(Fz1x_t, z_15, x_15, rep(t,m))
viwi.dx_15 <- cbind(vi.dx_15,wi.dx_15)
head(viwi.dx_15)
aux1.dx_15 <- cbind(x_15, viwi.dx_15)
head(aux1.dx_15)
part_AM.aux.dx_15 <- aux1.dx_15[aux1.dx_15[,1] <= mediana.x_15, ]  
part_DM.aux.dx_15 <- aux1.dx_15[aux1.dx_15[,1] >  mediana.x_15, ]
part_AM.dx_15 <- part_AM.aux.dx_15[,2:3] 
part_DM.dx_15 <- part_DM.aux.dx_15[,2:3] 
dev.new(); par(mfrow = c(1,3))
plot(part_AM.dx_15, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (Y,Z|X)", 
     sub = "X <= Mediana{X}")
plot(part_DM.dx_15, pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (Y,Z|X)",
     sub = "X > Mediana{X}")
plot(part_AM.dx_15, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (Y,Z|X)")
points(part_DM.dx_15, pch = 20, lwd = 2)
dif.copula(part_AM.dx_15, part_DM.dx_15,10000)


###(X,Z|Y)
(mediana.y_15 <- median(y_15))
ui.dy_15 <- mapply(Fx1y, x_15, y_15)
wi.dy_15 <- mapply(Fz1y, z_15, y_15, rep(t, m))
uiwi.dy_15 <- cbind(ui.dy_15, wi.dy_15)
head(uiwi.dy_15)
aux1.dy_15 <- cbind(y_15, uiwi.dy_15)
head(aux1.dy_15)
part_AM.aux.dy_15 <- aux1.dy_15[aux1.dy_15[,1] <= mediana.y_15, ]  
part_DM.aux.dy_15 <- aux1.dy_15[aux1.dy_15[,1] >  mediana.y_15, ]
part_AM.dy_15 <- part_AM.aux.dy_15[,2:3] 
part_DM.dy_15 <- part_DM.aux.dy_15[,2:3] 
dev.new(); par(mfrow = c(1,3))
plot(part_AM.dy_15, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Z|Y)", 
     sub = "Y <= Mediana{Y}")
plot(part_DM.dy_15, pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Z|Y)",
     sub = "Y > Mediana{Y}")
plot(part_AM.dy_15, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Z|Y)")
points(part_DM.dy_15, pch = 20, lwd = 2)
dif.copula(part_AM.dy_15, part_DM.dy_15,10000)
   

###(X,Y|Z)
(mediana.z_15 <- median(z_15))
ui.dz_15 <- mapply(Fx1z_t, x_15, z_15, rep(t,m))
vi.dz_15 <- mapply(Fy1z_t, y_15, z_15, rep(t,m))
uivi.dz_15 <- cbind(ui.dz_15, vi.dz_15)
head(uivi.dz_15)
aux1.dz_15 <- cbind(z_15, uivi.dz_15)
head(aux1.dz_15)
part_AM.aux.dz_15 <- aux1.dz_15[aux1.dz_15[,1] <= mediana.z_15, ]  
part_DM.aux.dz_15 <- aux1.dz_15[aux1.dz_15[,1] >  mediana.z_15, ]
part_AM.dz_15 <- part_AM.aux.dz_15[,2:3] 
part_DM.dz_15 <- part_DM.aux.dz_15[,2:3] 
dev.new(); par(mfrow = c(1,3))
plot(part_AM.dz_15, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Y|Z)", 
     sub = "Z <= Mediana{Z}")
plot(part_DM.dz_15, pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Y|Z)",
     sub = "Z > Mediana{Z}")
plot(part_AM.dz_15, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Y|Z)")
points(part_DM.dz_15, pch = 20, lwd = 2)
dif.copula(part_AM.dz_15, part_DM.dz_15,10000)


########## PARA TETA = 2
#####SIMULACIÓN DE (X,Y,Z)
m <- 5000
t <- 2
u <- runif(m,0,1)
x_2 <- invFx(u)
v <- runif(m,0,1)
y_2 <- mapply(invFylx, v, x_2)
w <- runif(m,0,1)
z_2 <- mapply(invFz1xy_t, w, x_2, y_2, rep(t, m))
ma.obs_2 <- cbind(x_2,y_2,z_2)
head(ma.obs_2)

###(Y,Z|X)
(mediana.x_2 <- median(x_2))
vi.dx_2 <- mapply(Fy1x, y_2, x_2)
wi.dx_2 <- mapply(Fz1x_t, z_2, x_2, rep(t,m))
viwi.dx_2 <- cbind(vi.dx_2,wi.dx_2)
head(viwi.dx_2)
aux1.dx_2 <- cbind(x_2, viwi.dx_2)
head(aux1.dx_2)
part_AM.aux.dx_2 <- aux1.dx_2[aux1.dx_2[,1] <= mediana.x_2, ]  
part_DM.aux.dx_2 <- aux1.dx_2[aux1.dx_2[,1] >  mediana.x_2, ]
part_AM.dx_2 <- part_AM.aux.dx_2[,2:3] 
part_DM.dx_2 <- part_DM.aux.dx_2[,2:3] 
dev.new(); par(mfrow = c(1,3))
plot(part_AM.dx_2, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (Y,Z|X)", 
     sub = "X <= Mediana{X}")
plot(part_DM.dx_2, pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (Y,Z|X)",
     sub = "X > Mediana{X}")
plot(part_AM.dx_2, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (Y,Z|X)")
points(part_DM.dx_2, pch = 20, lwd = 2)
dif.copula(part_AM.dx_2, part_DM.dx_2,10000)


###(X,Z|Y)
(mediana.y_2 <- median(y_2))
ui.dy_2 <- mapply(Fx1y, x_2, y_2)
wi.dy_2 <- mapply(Fz1y, z_2, y_2, rep(t, m))
uiwi.dy_2 <- cbind(ui.dy_2, wi.dy_2)
head(uiwi.dy_2)
aux1.dy_2 <- cbind(y_2, uiwi.dy_2)
head(aux1.dy_2)
part_AM.aux.dy_2 <- aux1.dy_2[aux1.dy_2[,1] <= mediana.y_2, ]  
part_DM.aux.dy_2 <- aux1.dy_2[aux1.dy_2[,1] >  mediana.y_2, ]
part_AM.dy_2 <- part_AM.aux.dy_2[,2:3] 
part_DM.dy_2 <- part_DM.aux.dy_2[,2:3] 
dev.new(); par(mfrow = c(1,3))
plot(part_AM.dy_2, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Z|Y)", 
     sub = "Y <= Mediana{Y}")
plot(part_DM.dy_2, pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Z|Y)",
     sub = "Y > Mediana{Y}")
plot(part_AM.dy_2, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Z|Y)")
points(part_DM.dy_2, pch = 20, lwd = 2)
dif.copula(part_AM.dy_2, part_DM.dy_2,10000)


###(X,Y|Z)
(mediana.z_2 <- median(z_2))
ui.dz_2 <- mapply(Fx1z_t, x_2, z_2, rep(t,m))
vi.dz_2 <- mapply(Fy1z_t, y_2, z_2, rep(t,m))
uivi.dz_2 <- cbind(ui.dz_2, vi.dz_2)
head(uivi.dz_2)
aux1.dz_2 <- cbind(z_2, uivi.dz_2)
head(aux1.dz_2)
part_AM.aux.dz_2 <- aux1.dz_2[aux1.dz_2[,1] <= mediana.z_2, ]  
part_DM.aux.dz_2 <- aux1.dz_2[aux1.dz_2[,1] >  mediana.z_2, ]
part_AM.dz_2 <- part_AM.aux.dz_2[,2:3] 
part_DM.dz_2 <- part_DM.aux.dz_2[,2:3] 
dev.new(); par(mfrow = c(1,3))
plot(part_AM.dz_2, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Y|Z)", 
     sub = "Z <= Mediana{Z}")
plot(part_DM.dz_2, pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Y|Z)",
     sub = "Z > Mediana{Z}")
plot(part_AM.dz_2, col = "blue", pch = 20, lwd = 2,
     main = "Pseudo-observaciones de la cópula de (X,Y|Z)")
points(part_DM.dz_2, pch = 20, lwd = 2)
dif.copula(part_AM.dz_2, part_DM.dz_2,10000)
