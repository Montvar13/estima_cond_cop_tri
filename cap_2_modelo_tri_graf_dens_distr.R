  #EJERCICIO 2
#GRÁFICAS DE MARGINALES

#X
fx <- function(x) { dunif(x) }
Fx <- function(x) {punif(x)}

u <- runif(10000, 0, 1)
x <- qunif(u)
hist(x, prob = T, col = "aquamarine", main = "X~Unif(0,1)", xlab = "x", breaks = 15)
lines(x, fx(x), col = "red", lwd =2)
aux <- seq(0.0001, 0.999, length.out = 10000)
plot(aux, fx(aux), main = "Función de densidad de X", col = "red",
     ylim = c(0,1.1), ylab = "fx(x)", xlab = "x")
#=================================================================================
#Y
fy <- function(y) { (-log(y))*((0 < y)&(y < 1))}
Fy <- function(y) { (y - y*log(y))*((0 < y)&(y < 1)) }
aux2 <- seq(0.0001, 0.999, length.out = 10000)
plot(aux2, fy(aux2), main = "Función de densidad de Y", col = "red",
     ylab = "fy(y)", xlab = "y", type = "l", lwd = 3)
m <- 10000
u <- runif(m, 0.00001, 0.9999)
invFy <- function(u) {uniroot(function(y) Fy(y) - u, interval = c(0.000001,0.99999))}$root
y <- numeric(m)
for(k in 1:m){
y[k] <- invFy(u[k])
}
dev.new()
hist(y, freq = FALSE, col = "aquamarine", main = "Y~fy(y)", xlab = "y", breaks = 15)
lines(sort(y), sapply(sort(y), fy), col = "red", lwd = 2)
legend(0.7, 3, legend = "fy(y)", col = "red", lty = 1, cex = 1)

#=================================================================================
#Z

fz <- function(z,teta) {
            if(teta == 0) {fz <- function(z){(-log(z))*((0 < z)&(z < 1))}}
            if(teta == 1) {fz <- function(z) {(log(z) + 2/sqrt(z) - 2)*((0 < z)&(z < 1))}}
            if(teta > 0 & teta != 1) {fz <- function(z) {(z^(-teta/(teta + 1)) - 1)/(teta - teta^2) + (teta*z^((1 - teta)/teta)*(1 - 1/(z^(1/(teta^2 + teta)))))/(1 - teta)}

}
fz(z)
}

Fz <- function(z, teta) {
              if(teta == 0) {Fz <- function(z){z - z*log(z)}}
              if(teta == 1) {Fz <- function(z){z*log(z) + 4*sqrt(z) - 3*z}}
              if(teta > 0 & teta !=1) {Fz <- function(z) {(1/(1-teta))*(-((teta-1)*(teta+1)^2*z^(1/(teta+1))+z-teta^3*z^(1/teta))/(teta))}}
Fz(z)
}

u <- runif(m, 0.001, 0.999)

invFz <- function(u, teta){ 
        invFz <- function(u){uniroot(function(z) Fz(z,teta) - u, interval = c(0.000000001, 0.999999))$root}
invFz(u)
}
z <- numeric(m)
#TETA = 0
for(k in 1:m) {
z[k] <- invFz(u[k],0)
}

#TETA = 1
u <- runif(m, 0.1,0.999)
z1 <- numeric(m)
for(k in 1:m) {
z1[k] <- invFz(u[k],1)
}

#TETA = 4
u <- runif(m, 0.10,0.999)
z4 <- numeric(m)
for(k in 1:m) {
z4[k] <- invFz(u[k],4)
}

dev.new(); par(mfrow = c(2,2))
hist(z, prob = TRUE, main = "Z~Fz(z), teta = 0", col = "aquamarine")
lines(sort(z), mapply(fz, sort(z), rep(0, m)), col = "red", lwd = 2)
legend(0.7 , 3.5, legend = "fz(z)", col = "red", lty = 1, cex = 0.7)

hist(z1, prob = TRUE, main = "Z~Fz(z), teta = 1", col = "aquamarine")
lines(sort(z1), mapply(fz, sort(z1), rep(1, m)), col = "red", lwd = 2)
legend(0.5 , 10, legend = "fz(z)", col = "red", lty = 1, cex = 0.7)

hist(z4, prob = TRUE, main = "Z~Fz(z), teta = 4", col = "aquamarine")
lines(sort(z4), mapply(fz, sort(z4), rep(4, m)), col = "red", lwd = 2)
legend(0.35 , 13, legend = "fz(z)", col = "red", lty = 1, cex = 0.7)
par(mfrow = c(2,2))
aux3 <- seq(0.000001,0.999999, length.out = 100000)
plot(aux3, mapply(fz, aux3, rep(0,100000)), main = "Densidad de Z, teta = 0",
     col = "red", type = "l", lwd = 3, ylab = "fz(z)", xlab = "z", ylim = c(0,10))
aux4 <- seq(0.000001,0.999999, length.out = 100000)
plot(aux4, mapply(fz, aux4, rep(1,100000)), main = "Densidad de Z, teta = 1",
     col = "red", type = "l", lwd = 3, ylab = "fz(z)", xlab = "z", ylim = c(0,10))
aux5 <- seq(0.000001,0.999999, length.out = 100000)
plot(aux5, mapply(fz, aux5, rep(0.5,100000)), main = "Densidad de Z, teta = 0.5",
     col = "red", type = "l", lwd = 3, ylab = "fz(z)", xlab = "z", ylim = c(0,10))
aux6 <- seq(0.000001,0.999999, length.out = 100000)
plot(aux6, mapply(fz, aux6, rep(1.5,100000)), main = "Densidad de Z, teta = 1.5",
     col = "red", type = "l", lwd = 3, ylab = "fz(z)", xlab = "z", ylim = c(0,10))



#===============================================================================

###CONJUNTA BIVARIADA
##XY

fxy <- function(x,y) {(1/x)*((0 < y)&(y < x)&(x < 1))} 
x <- seq(0.1, 0.99, length.out = 70)
y <- seq(0.1, 0.99, length.out = 70)
fxy_ev <- outer(x, y, fxy)
dev.new(); par(mfrow = c(1,2))
persp(x,y,fxy_ev, col = "yellow", main = "Densidad del vector aleatorio (X,Y)",
      zlab = "fxy(x,y)", theta = 90)
image(x, y, fxy_ev, col = heat.colors(20), main = "Conjuntos de nivel de fxy")

##YZ

fyz <- function(y,z,teta) {(((1/z - 1/y^teta)*(z >= y^(teta+1)))+((1/y^(teta+1) - 1/y^teta)*(z < y^(teta+1))))*(((0 < z)&(z < y^teta))&((0 < y)&(y < 1)))}
z <- seq(0.1, 0.99, length.out = 60)
y <- seq(0.1, 0.99, length.out = 60)

#TETA = 0*
fyz_ev <- outer(y, z, function(y,z) fyz(y,z,0))
dev.new(); par(mfrow = c(1,2))
persp(y,z,fyz_ev, col = "yellow", main = "Densidad del vector aleatorio (Y,Z), teta = 0",
      zlab = "fyz(y,z)", theta = 90)
image(y, z, fyz_ev, col = heat.colors(20), main = "Conjuntos de nivel de fyz, teta = 0")

#TETA = 0.1
fyz_ev.1 <- outer(y, z, function(y,z) fyz(y,z,0.1))
dev.new(); par(mfrow = c(1,2))
persp(y,z,fyz_ev.1, col = "yellow", main = "Densidad del vector aleatorio (Y,Z), teta = 0.1",
      zlab = "fyz(y,z)", theta = 90)
image(y, z, fyz_ev.1, col = heat.colors(20), main = "Conjuntos de nivel de fyz, teta = 0.1")

#TETA = 0.50*

fyz_ev.5 <- outer(y, z, function(y,z) fyz(y,z,0.5))
dev.new(); par(mfrow = c(1,2))
persp(y,z,fyz_ev.5, col = "yellow", main = "Densidad del vector aleatorio (Y,Z), teta = 0.5",
      zlab = "fyz(y,z)", theta = 100)
image(y, z, fyz_ev.5, col = heat.colors(20), main = "Conjuntos de nivel de fyz, teta = 0.5")



#TETA = 0.95
fyz_ev.95 <- outer(y, z, function(y,z) fyz(y,z,0.95))
dev.new(); par(mfrow = c(1,2))
persp(y,z,fyz_ev.95, col = "yellow", main = "Densidad del vector aleatorio (Y,Z), teta = 0.95",
      zlab = "fyz(y,z)", theta = 130, phi = 30)
image(y, z, fyz_ev.95, col = heat.colors(20), main = "Conjuntos de nivel de fyz, teta = 0.95")

#TETA = 1

fyz_ev1 <- outer(y, z, function(y,z) fyz(y,z,1))
dev.new(); par(mfrow = c(1,2))
persp(y,z,fyz_ev1, col = "yellow", main = "Densidad del vector aleatorio (Y,Z), teta = 1",
      zlab = "fyz(y,z)", theta = 110)
image(y, z, fyz_ev1, col = heat.colors(20), main = "Conjuntos de nivel de fyz, teta = 1")

#TETA = 1.1

fyz_ev1.1 <- outer(y, z, function(y,z) fyz(y,z,1.1))
dev.new(); par(mfrow = c(1,2))
persp(y,z,fyz_ev1.1, col = "yellow", main = "Densidad del vector aleatorio (Y,Z), teta = 1.1",
      zlab = "fyz(y,z)", theta = 0)
image(y, z, fyz_ev1.1, col = heat.colors(20), main = "Conjuntos de nivel de fyz, teta = 1.1")


#TETA = 1.5

fyz_ev1.5 <- outer(y, z, function(y,z) fyz(y,z,1.5))
dev.new(); par(mfrow = c(1,2))
persp(y,z,fyz_ev1.5, col = "yellow", main = "Densidad del vector aleatorio (Y,Z), teta = 1.5",
      zlab = "fyz(y,z)", theta = 140)
image(y, z, fyz_ev1.5, col = heat.colors(20), main = "Conjuntos de nivel de fyz, teta = 1.5")


#XZ

fxz <- function(x,z,teta) {
            if(teta == 0) {
                   fxz <- function(x,z) {(1/x)*((0 < z)&(z < x)&(x < 1))} }
            if(teta == 1) { 
                   fxz <- function(x,z) {((2*log(x) - log(z))/(x^2))*(((0 < z)&(z < x^2))&((0 < x)&(x < 1)))} }
            if(teta > 0 & teta !=1) {
                   fxz <- function(x,z) {(1/(1-teta))*(1/x^(teta+1) - z^((1-teta)/teta)/x^((1+teta)/teta))*(((0 < z)&(z < x^(teta+1)))&((0 < x)&(x < 1)))}}
fxz(x,z)
}	
x <- seq(0.1, 0.99, length.out = 60)
z <- seq(0.1, 0.99, length.out = 60)

#TETA = 0

fxz_ev0 <- outer(x,z,function(x,z) fxz(x,z,0))
dev.new(); par(mfrow = c(1,2))
persp(x,z,fxz_ev0, col = "yellow", main = "Densidad del vector aleatorio (X,Z), teta = 0",
      zlab = "fxz(x,z)", theta = 90)
image(x, z, fxz_ev0, col = heat.colors(20), main = "Conjuntos de nivel de fxz, teta = 0")

#TETA = 0.1
fxz_ev.1 <- outer(x, z, function(x,z) fxz(x,z,0.1))
dev.new(); par(mfrow = c(1,2))
persp(x,z,fxz_ev.1, col = "yellow", main = "Densidad del vector aleatorio (X,Z), teta = 0.1",
      zlab = "fxz(x,z)", phi = 30)
image(x, z, fxz_ev.1, col = heat.colors(20), main = "Conjuntos de nivel de fxz, teta = 0.1")

#TETA = 0.5
fxz_ev.5 <- outer(x, z, function(x,z) fxz(x,z,0.5))
dev.new(); par(mfrow = c(1,2))
persp(x,z,fxz_ev.5, col = "yellow", main = "Densidad del vector aleatorio (X,Z), teta = 0.5",
      zlab = "fxz(x,z)", phi = 10)
image(x, z, fxz_ev.5, col = heat.colors(20), main = "Conjuntos de nivel de fxz, teta = 0.5")


#TETA = 1
fxz_ev1 <- outer(x,z,function(x,z) fxz(x,z,1))
dev.new(); par(mfrow = c(1,2))
persp(x,z,fxz_ev1, col = "yellow", main = "Densidad del vector aleatorio (X,Z), teta = 1",
      zlab = "fxz(x,z)", theta = 0)
image(x, z, fxz_ev1, col = heat.colors(20), main = "Conjuntos de nivel de fxz, teta = 1")


#TETA = 1.5
fxz_ev1.5 <- outer(x,z,function(x,z) fxz(x,z,1.5))
dev.new(); par(mfrow = c(1,2))
persp(x,z,fxz_ev1.5, col = "yellow", main = "Densidad del vector aleatorio (X,Z), teta = 1.5",
      zlab = "fxz(x,z)", theta = 0)
image(x, z, fxz_ev1.5, col = heat.colors(20), main = "Conjuntos de nivel de fxz, teta = 1.5")

#TETA = 2
fxz_ev2 <- outer(x,z,function(x,z) fxz(x,z,2))
dev.new(); par(mfrow = c(1,2))
persp(x,z,fxz_ev2, col = "yellow", main = "Densidad del vector aleatorio (X,Z), teta = 2",
      zlab = "fxz(x,z)", theta = 0)
image(x, z, fxz_ev2, col = heat.colors(20), main = "Conjuntos de nivel de fxz, teta = 2")
#===========================================================================================

###CONDICIONALES(1 VARIABLES)

##XY|Z 

fxylz <- function(x,y,z,teta) {
                 if(teta == 0) {fxy1z <- function(x,y,z){(-1/(log(z)*x^2))*((((abs(z-y)+z+y)/2 < x)&(x < 1))&((0 < y)&(y < 1)))}} 
                 if(teta == 1) {fxy1z <- function(x,y,z){(1/(x^2*y*(log(z)+2/sqrt(z)-2)))*(((z/x < y)&(y < 1))&((y < x)&(x < 1)))}}
                 if(teta > 0 & teta!= 1) {fxy1z <- function(x,y,z) {((teta-teta^2)/(x^2*y^teta*(z^(-teta/(teta+1))-1+teta^2*z^((1-teta)/teta)*(1-1/z^(1/(teta^2+teta))))))*((((z/x)^(1/teta) < y)&(y < x))&((z^(1/(teta+1)) < x)&(x < 1)))}}
fxy1z(x,y,z)
}			

x <- seq(0.01, 0.99, length.out = 45)
y <- seq(0.01, 0.99, length.out = 45)

#TETA = 0 |	 z = 0.25

fxy1z_0.25_ev0 <- outer(x,y,function(x,y) fxylz(x,y,0.25,0))
dev.new(); par(mfrow = c(1,2))
persp(x, y, fxy1z_0.25_ev0, col = "yellow", main = "Densidad del vector aleatorio (X,Y |Z = 0.25) teta = 0",
      zlab = "fxy1z(x,y)", phi = 25, teta = 80)
image(x, y, fxy1z_0.25_ev0, col = heat.colors(20), main = "Conjuntos de nivel de fxy|z = 0.25 teta = 0")
                 


#TETA = 0 | z = 0.5

fxy1z_0.5_ev0 <- outer(x,y,function(x,y) fxylz(x,y,0.5,0))
dev.new(); par(mfrow = c(1,2))
persp(x, y, fxy1z_0.5_ev0, col = "yellow", main = "Densidad del vector aleatorio (X,Y |Z = 0.5) teta = 0",
      zlab = "fxy1z(x,y)", phi = 40, teta = 80)
image(x, y, fxy1z_0.5_ev0, col = heat.colors(20), main = "Conjuntos de nivel de fxy|z = 0.5 teta = 0")
                 

#TETA = 0 |	 z = 0.75

fxy1z_0.75_ev0 <- outer(x,y,function(x,y) fxylz(x,y,0.75,0))
dev.new(); par(mfrow = c(1,2))
persp(x, y, fxy1z_0.75_ev0, col = "yellow", main = "Densidad del vector aleatorio (X,Y |Z = 0.75) teta = 0",
      zlab = "fxy1z(x,y)", teta = 80, phi = 25)
image(x, y, fxy1z_0.75_ev0, col = heat.colors(20), main = "Conjuntos de nivel de fxy|z = 0.75 teta = 0")
                 
#TETA = 1 | z = 0.25

fxy1z_0.25_ev1 <- outer(x,y,function(x,y) fxylz(x,y,0.25,1))
dev.new(); par(mfrow = c(1,2))
persp(x, y, fxy1z_0.25_ev1, col = "yellow", main = "Densidad del vector aleatorio (X,Y |Z = 0.25) teta = 1",
      zlab = "fxy1z(x,y)", theta = 25, phi = 35)
image(x, y, fxy1z_0.25_ev1, col = heat.colors(20), main = "Conjuntos de nivel de fxy|z = 0.25 teta = 1")

#TETA = 1 | z = 0.5

fxy1z_0.5_ev1 <- outer(x,y,function(x,y) fxylz(x,y,0.5,1))
dev.new(); par(mfrow = c(1,2))
persp(x, y, fxy1z_0.5_ev1, col = "yellow", main = "Densidad del vector aleatorio (X,Y |Z = 0.5) teta = 1",
      zlab = "fxy1z(x,y)", theta = 70, phi = 40)
image(x, y, fxy1z_0.5_ev1, col = heat.colors(20), main = "Conjuntos de nivel de fxy|z = 0.5 teta = 1")

x <- seq(0.01, 0.99, length.out = 65)
y <- seq(0.01, 0.99, length.out = 65)


#TETA = 0.5 | z = 0.5

fxy1z_0.5_ev0.5 <- outer(x,y,function(x,y) fxylz(x,y,0.5,0.5))
dev.new(); par(mfrow = c(1,2))
persp(x, y, fxy1z_0.5_ev0.5, col = "yellow", main = "Densidad del vector aleatorio (X,Y |Z = 0.5) teta = 0.5",
      zlab = "fxy1z(x,y)", theta = 40, phi = 35)
image(x, y, fxy1z_0.5_ev0.5, col = heat.colors(20), main = "Conjuntos de nivel de fxy|z = 0.5 teta = 0.5")

#TETA = 2 | z = 0.25

fxy1z_0.25_ev2 <- outer(x,y,function(x,y) fxylz(x,y,0.25,2))
dev.new(); par(mfrow = c(1,2))
persp(x, y, fxy1z_0.25_ev2, col = "yellow", main = "Densidad del vector aleatorio (X,Y |Z = 0.25) teta = 2",
      zlab = "fxy1z(x,y)", theta = 90, phi = 20)
image(x, y, fxy1z_0.25_ev2, col = heat.colors(20), main = "Conjuntos de nivel de fxy|z = 0.25 teta = 2")


#TETA = 0.15 | z = 0.70

fxy1z_0.70_ev0.15 <- outer(x,y,function(x,y) fxylz(x,y,0.70,0.15))
dev.new(); par(mfrow = c(1,2))
persp(x, y, fxy1z_0.70_ev0.15, col = "yellow", main = "Densidad del vector aleatorio (X,Y |Z = 0.70) teta = 0.15",
      zlab = "fxy1z(x,y)", theta = 40)
image(x, y, fxy1z_0.70_ev0.15, col = heat.colors(20), main = "Conjuntos de nivel de fxy|z = 0.70 teta = 0.15")

#TETA = 2 | z = 0.70

fxy1z_0.70_ev2 <- outer(x,y,function(x,y) fxylz(x,y,0.70,2))
dev.new(); par(mfrow = c(1,2))
persp(x, y, fxy1z_0.70_ev2, col = "yellow", main = "Densidad del vector aleatorio (X,Y |Z = 0.70) teta = 2",
      zlab = "fxy1z(x,y)", theta = 110, phi = 30)
image(x, y, fxy1z_0.70_ev2, col = heat.colors(20), main = "Conjuntos de nivel de fxy|z = 0.70 teta = 2")
 

##YZ|X

fyz1x <- function(y,z,x,teta) {(1/(x^2*y^teta))*(((0 < z)&(z < x*y^teta))&((0 < y)&(y < x)))}                

y <- seq(0.1, 0.9, length.out = 50)
z <- seq(0.1, 0.9, length.out = 50)

#TETA = 0, x = 0.25
fyz1x_0.25_ev0 <- outer(y, z, function(y,z) fyz1x(y,z,0.25,0))
dev.new(); par(mfrow = c(1,2))
persp(y, z, fyz1x_0.25_ev0, col = "yellow", main = "Densidad del vector aleatorio (Y,Z |X = 0.25) teta = 0",
      zlab = "fyz1x(y,z)", theta = 70, phi = 45)
image(y, z, fyz1x_0.25_ev0, col = heat.colors(20), main = "Conjuntos de nivel de fyz|x = 0.25 teta = 0")

#TETA = 0, x = 0.75
fyz1x_0.75_ev0 <- outer(y, z, function(y,z) fyz1x(y,z,0.75,0))
dev.new(); par(mfrow = c(1,2))
persp(y, z, fyz1x_0.75_ev0, col = "yellow", main = "Densidad del vector aleatorio (Y,Z |X = 0.75) teta = 0",
      zlab = "fyz1x(y,z)", theta = 70, phi = 45)
image(y, z, fyz1x_0.75_ev0, col = heat.colors(20), main = "Conjuntos de nivel de fyz|x = 0.75 teta = 0")



#TETA = 0.5, x = 0.45
fyz1x_0.45_ev0.5 <- outer(y, z, function(y,z) fyz1x(y,z,0.45,0.5))
dev.new(); par(mfrow = c(1,2))
persp(y, z, fyz1x_0.45_ev0.5, col = "yellow", main = "Densidad del vector aleatorio (Y,Z |X = 0.45) teta = 0.5",
      zlab = "fyz1x(y,z)", theta = 70, phi = 15)
image(y, z, fyz1x_0.45_ev0.5, col = heat.colors(20), main = "Conjuntos de nivel de fyz|x = 0.45 teta = 0.5")

#TETA = 0.5, x = 0.75
fyz1x_0.75_ev0.5 <- outer(y, z, function(y,z) fyz1x(y,z,0.75,0.5))
dev.new(); par(mfrow = c(1,2))
persp(y, z, fyz1x_0.75_ev0.5, col = "yellow", main = "Densidad del vector aleatorio (Y,Z |X = 0.75) teta = 0.5",
      zlab = "fyz1x(y,z)", theta = 70, phi = 45)
image(y, z, fyz1x_0.75_ev0.5, col = heat.colors(20), main = "Conjuntos de nivel de fyz|x = 0.75 teta = 0.5")



#TETA = 1, x = 0.45
fyz1x_0.45_ev1 <- outer(y, z, function(y,z) fyz1x(y,z,0.45,1))
dev.new(); par(mfrow = c(1,2))
persp(y, z, fyz1x_0.45_ev1, col = "yellow", main = "Densidad del vector aleatorio (Y,Z |X = 0.45) teta = 1",
      zlab = "fyz1x(y,z)", theta = 80, phi = 15)
image(y, z, fyz1x_0.45_ev1, col = heat.colors(20), main = "Conjuntos de nivel de fyz|x = 0.45 teta = 1")

#TETA = 1, x = 0.75
fyz1x_0.75_ev1 <- outer(y, z, function(y,z) fyz1x(y,z,0.75,1))
dev.new(); par(mfrow = c(1,2))
persp(y, z, fyz1x_0.75_ev1, col = "yellow", main = "Densidad del vector aleatorio (Y,Z |X = 0.75) teta = 1",
      zlab = "fyz1x(y,z)", theta = 80, phi = 15)
image(y, z, fyz1x_0.75_ev1, col = heat.colors(20), main = "Conjuntos de nivel de fyz|x = 0.75 teta = 1")


#TETA = 1.5, x = 0.45
fyz1x_0.45_ev1.5 <- outer(y, z, function(y,z) fyz1x(y,z,0.45,1.5))
dev.new(); par(mfrow = c(1,2))
persp(y, z, fyz1x_0.45_ev1.5, col = "yellow", main = "Densidad del vector aleatorio (Y,Z |X = 0.45) teta = 1.5",
      zlab = "fyz1x(y,z)", theta = 50, phi = 15)
image(y, z, fyz1x_0.45_ev1.5, col = heat.colors(20), main = "Conjuntos de nivel de fyz|x = 0.45 teta = 1.5")

#TETA = 1.5, x = 0.75
fyz1x_0.75_ev1.5 <- outer(y, z, function(y,z) fyz1x(y,z,0.75,1.5))
dev.new(); par(mfrow = c(1,2))
persp(y, z, fyz1x_0.75_ev1.5, col = "yellow", main = "Densidad del vector aleatorio (Y,Z |X = 0.75) teta = 1.5",
      zlab = "fyz1x(y,z)", theta = 50, phi = 15)
image(y, z, fyz1x_0.75_ev1.5, col = heat.colors(20), main = "Conjuntos de nivel de fyz|x = 0.75 teta = 1.5")


##XZ|Y

fxz1y <- function(x,z,y,teta) {(-1/(x^2*y^teta*log(y)))*(((0 < z)&(z < x*y^teta))&((y < x)&(x < 1)))}
x <- seq(0.1,0.9, length.out = 60)
z <- seq(0.1,0.9, length.out = 60)

#TETA = 0, y = 0.25

fxz1y_0.25_ev0 <- outer(x,z,function(x,z) fxz1y(x,z,0.25,0))
dev.new(); par(mfrow = c(1,2))
persp(x, z, fxz1y_0.25_ev0, col = "yellow", main = "Densidad del vector aleatorio (X,Z |Y = 0.25) teta = 0",
      zlab = "fxz1y(x,z)", theta = 58, phi = 15)
image(x, z, fxz1y_0.25_ev0, col = heat.colors(20), main = "Conjuntos de nivel de fxz|y = 0.25 teta = 0")

#TETA = 0, y = 0.75

fxz1y_0.75_ev0 <- outer(x,z,function(x,z) fxz1y(x,z,0.75,0))
dev.new(); par(mfrow = c(1,2))
persp(x, z, fxz1y_0.75_ev0, col = "yellow", main = "Densidad del vector aleatorio (X,Z |Y = 0.75) teta = 0",
      zlab = "fxz1y(x,z)", theta = 50, phi = 15)
image(x, z, fxz1y_0.75_ev0, col = heat.colors(20), main = "Conjuntos de nivel de fxz|y = 0.75 teta = 0")


#TETA = 1, y = 0.25

fxz1y_0.25_ev1 <- outer(x,z,function(x,z) fxz1y(x,z,0.25,1))
dev.new(); par(mfrow = c(1,2))
persp(x, z, fxz1y_0.25_ev1, col = "yellow", main = "Densidad del vector aleatorio (X,Z |Y = 0.25) teta = 1",
      zlab = "fxz1y(x,z)", theta = 50, phi = 15)
image(x, z, fxz1y_0.25_ev1, col = heat.colors(20), main = "Conjuntos de nivel de fxz|y = 0.25 teta = 1")

#TETA = 1, y = 0.75

fxz1y_0.75_ev1 <- outer(x,z,function(x,z) fxz1y(x,z,0.75,1))
dev.new(); par(mfrow = c(1,2))
persp(x, z, fxz1y_0.75_ev1, col = "yellow", main = "Densidad del vector aleatorio (X,Z |Y = 0.75) teta = 1",
      zlab = "fxz1y(x,z)", theta = 50, phi = 15)
image(x, z, fxz1y_0.75_ev1, col = heat.colors(20), main = "Conjuntos de nivel de fxz|y = 0.75 teta = 1")




#TETA = 0.5, y = 0.25

fxz1y_0.25_ev0.5 <- outer(x,z,function(x,z) fxz1y(x,z,0.25,0.5))
dev.new(); par(mfrow = c(1,2))
persp(x, z, fxz1y_0.25_ev0.5, col = "yellow", main = "Densidad del vector aleatorio (X,Z |Y = 0.25) teta = 0.5",
      zlab = "fxz1y(x,z)", theta = 50, phi = 15)
image(x, z, fxz1y_0.25_ev0.5, col = heat.colors(20), main = "Conjuntos de nivel de fxz|y = 0.25 teta = 0.5")

#TETA = 0.5, y = 0.75

fxz1y_0.75_ev0.5 <- outer(x,z,function(x,z) fxz1y(x,z,0.75,0.5))
dev.new(); par(mfrow = c(1,2))
persp(x, z, fxz1y_0.75_ev0.5, col = "yellow", main = "Densidad del vector aleatorio (X,Z |Y = 0.75) teta = 0.5",
      zlab = "fxz1y(x,z)", theta = 50, phi = 15)
image(x, z, fxz1y_0.75_ev0.5, col = heat.colors(20), main = "Conjuntos de nivel de fxz|y = 0.75 teta = 0.5")


#TETA = 1.5, y = 0.25

fxz1y_0.25_ev1.5 <- outer(x,z,function(x,z) fxz1y(x,z,0.25,1.5))
dev.new(); par(mfrow = c(1,2))
persp(x, z, fxz1y_0.25_ev1.5, col = "yellow", main = "Densidad del vector aleatorio (X,Z |Y = 0.25) teta = 1.5",
      zlab = "fxz1y(x,z)", theta = 50, phi = 15)
image(x, z, fxz1y_0.25_ev1.5, col = heat.colors(20), main = "Conjuntos de nivel de fxz|y = 0.25 teta = 1.5")

#TETA = 1.5, y = 0.75

fxz1y_0.75_ev1.5 <- outer(x,z,function(x,z) fxz1y(x,z,0.75,1.5))
dev.new(); par(mfrow = c(1,2))
persp(x, z, fxz1y_0.75_ev1.5, col = "yellow", main = "Densidad del vector aleatorio (X,Z |Y = 0.75) teta = 1.5",
      zlab = "fxz1y(x,z)", theta = 50, phi = 15)
image(x, z, fxz1y_0.75_ev1.5, col = heat.colors(20), main = "Conjuntos de nivel de fxz|y = 0.75 teta = 1.5")

#=================================================================================

###CONDICIONAL EN 2 VARIABLES.

##X|YZ


invFx1yz <- function(u,y,z,teta) {(z/(y^teta-u*(y^teta-z)))*(z/y^teta >= y)+(y/(1-u*(1-y)))*(z/y^teta < y)}
fx1yz <- function(x,y,z,teta) {((z/(x^2*(y^(teta)-z)))*((y <= z/y^teta) & (z/y^teta < x) & (x < 1)))+((y/(x^2*(1-y)))*((z/y^teta < y) & (y < x) & (x < 1)))}
m <- 10000
u <- runif(m,0,1)
x <- seq(0.1,0.9,length.out = m)
xs <- seq(0.01,0.99,length.out = m)

#TETA = 0; Y = 0.2; Z = 0.2
y <- 0.2; z <- 0.2; teta <- 0
x_0.2_0.2_0 <- mapply(invFx1yz,u, rep(y,m), rep(z,m), rep(teta,m))
fx1yz_0.2_0.2_ev0 <- mapply(fx1yz, xs, rep(y,m), rep(z,m), rep(teta,m))

#TETA = 0; Y = 0.2; Z = 0.5
y <- 0.2; z <- 0.5; teta <- 0
x_0.2_0.5_0 <- mapply(invFx1yz,u, rep(y,m), rep(z,m), rep(teta,m))
fx1yz_0.2_0.5_ev0 <- mapply(fx1yz, xs, rep(y,m), rep(z,m), rep(teta,m))


#TETA = 1; Y = 0.5; Z = 0.2
y <- 0.5; z <- 0.2; teta <- 1
x_0.5_0.2_1 <- mapply(invFx1yz,u, rep(y,m), rep(z,m), rep(teta,m))
fx1yz_0.5_0.2_ev1 <- mapply(fx1yz, xs, rep(y,m), rep(z,m), rep(teta,m))

#TETA = 0.5; Y = 0.5; Z = 0.2
y <- 0.5; z <- 0.2; teta <- 0.5
x_0.5_0.2_0.5 <- mapply(invFx1yz,u, rep(y,m), rep(z,m), rep(teta,m))
fx1yz_0.5_0.2_ev0.5 <- mapply(fx1yz, xs, rep(y,m), rep(z,m), rep(teta,m))

dev.new(); par(mfrow = c(1,2))
plot(xs, fx1yz_0.2_0.2_ev0, col = "red", lwd = 1,
     main = "Densidad de (X|Y = 0.2, Z = 0.2)", xlab = "x", ylab = "fx|yz(x|0.2,0.2)")
hist(x_0.2_0.2_0, prob = TRUE, main = "Histograma (X|Y=0.2, Z=0.2), teta = 0", xlim = c(0,1),
    col = "yellow", xlab = "x")
lines(xs, fx1yz_0.2_0.2_ev0, col = "red", lwd = 2)
legend(0.5,4, legend = "fx|yz(x)", col = "red", lty = 1, cex = 0.8, box.lty = 0) 

hist(x_0.2_0.5_0, prob = TRUE, main = "Histograma (X|Y=0.2, Z=0.5), teta = 0", xlim = c(0,1),
    col = "yellow", xlab = "x")
lines(xs, fx1yz_0.2_0.5_ev0, col = "red", lwd = 2)
legend(0.6, 3.5, legend = "fx|yz(x)", col = "red", lty = 1, cex = 0.8, box.lty = 0) 

dev.new(); par(mfrow = c(1,2))
plot(xs, fx1yz_0.5_0.2_ev1, col = "red", lwd = 2,
     main = "Densidad de (X|Y = 0.5, Z = 0.1)", xlab = "x", 
     ylab = "fx|yz(x|0.5,0.1)")
hist(x_0.5_0.2_1, prob = TRUE, main = "Histograma (X|Y=0.5, Z=0.1), teta = 1", xlim = c(0,1),
    col = "yellow", xlab = "x")
lines(xs, fx1yz_0.5_0.2_ev1, col = "red", lwd = 2)
legend(0.6, 3, legend = "fx|yz(x)", col = "red", lty = 1, cex = 0.8, box.lty = 0) 

hist(x_0.5_0.2_0.5, prob = TRUE, main = "Histograma (X|Y=0.5, Z=0.1), teta = 0.5", xlim = c(0,1),
    col = "yellow", xlab = "x")
lines(xs, fx1yz_0.5_0.2_ev0.5, col = "red", lwd = 2)
legend(0.6, 3, legend = "fx|yz(x)", col = "red", lty = 1, cex = 0.8, box.lty = 0) 
 

##Y|XZ 

invFy1xz <- function(u,x,z,teta){
                  if(teta == 0){invFy1xz <- function(u,x,z){(u*x)*(z <= x)}} 
                  if(teta == 1){invFy1xz <- function(u,x,z){(x^(2*u-1)/z^(u-1))*(z <= x^2)}}
       if(teta > 0 & teta != 1){invFy1xz <- function(u,x,z){(u*(x^(1-teta)-z^((1-teta)/teta)*x^((teta-1)/teta))+ (z/x)^((1-teta)/teta))^(1/(1-teta))*(z <= x^(teta+1))}}
invFy1xz(u,x,z)
}
fy1xz <- function(y,x,z,teta) {
            if(teta == 0){fy1xz <- function(y,x,z){(1/x)*(((0 <= y)&(y <= x))&(z <= x))}} 
            if(teta == 1){fy1xz <- function(y,x,z){(1/(y*(2*log(x)-log(z))))*(((z/x < y)&(y < x))&(z < x^2))}}
            if(teta > 0 & teta!=1) {fy1xz <- function(y,x,z){((1-teta)/(y^teta*(x^(1-teta)-z^((1-teta)/teta)*x^((teta-1)/teta))))*((((z/x)^(1/teta) < y)&(y < x))&(z < x^(teta+1)))}}
fy1xz(y,x,z)
}
m <- 10000
u <- runif(m,0,1)
ys <- seq(0.01,0.99,length.out = m)

#TETA = 0, x  = 0.5, z = 0.25
teta <- 0; x <- 0.5; z <- 0.25
y_0.5_0.25_0 <- mapply(invFy1xz, u, rep(x,m), rep(z,m), rep(teta,m))
fy1xz_0.5_0.25_ev0 <- mapply(fy1xz, sort(y_0.5_0.25_0), rep(x,m), rep(z,m), rep(teta,m))
hist(y_0.5_0.25_0, prob = TRUE, main = "Histograma (Y|X=0.5, Z=0.1)~Unif(0,0.5), teta = 0",
    col = "yellow", xlab = "y") 
lines(sort(y_0.5_0.25_0), fy1xz_0.5_0.25_ev0, col = "red", lwd = 2)

#TETA = 1, x  = 0.7, z = 0.3
teta <- 1; x <- 0.7; z <- 0.3
y_0.7_0.3_1 <- mapply(invFy1xz, u, rep(x,m), rep(z,m), rep(teta,m))
fy1xz_0.7_0.3_ev1 <- mapply(fy1xz, ys, rep(x,m), rep(z,m), rep(teta,m))

#TETA = 1, x  = 0.5, z = 0.2
teta <- 1; x <- 0.5; z <- 0.2
y_0.5_0.2_1 <- mapply(invFy1xz, u, rep(x,m), rep(z,m), rep(teta,m))
fy1xz_0.5_0.2_ev1 <- mapply(fy1xz, ys, rep(x,m), rep(z,m), rep(teta,m))

dev.new();par(mfrow = c(1,2))
hist(y_0.7_0.3_1, prob = TRUE, main = "Histograma (Y|X=0.7, Z=0.3), teta = 1",
    col = "yellow", xlab = "y", xlim = c(0,1)) 
lines(ys, fy1xz_0.7_0.3_ev1, col = "red", lwd = 2)
legend(0.6, 4, legend = "fy|xz(y)", col = "red", lty = 1, cex = 0.8, box.lty = 0)

hist(y_0.5_0.2_1, prob = TRUE, main = "Histograma (Y|X=0.5, Z=0.2), teta = 1",
    col = "yellow", xlab = "y", xlim = c(0,1)) 
lines(ys, fy1xz_0.5_0.2_ev1, col = "red", lwd = 2)
legend(0.6, 4, legend = "fy|xz(y)", col = "red", lty = 1, cex = 0.8, box.lty = 0)

#========================================================================

#TETA = 0.5, x  = 0.7, z = 0.3
teta <- 0.5; x <- 0.7; z <- 0.3

fy1xz_0.7_0.3_ev0.5 <- mapply(fy1xz, y, rep(x,m), rep(z,m), rep(teta,m))
dev.new()
plot(y,fy1xz_0.7_0.3_ev0.5, type = "l", col = "red", lwd = 2,
    main = "Y|X = 0.7, Z = 0.3 Teta = 0.5", ylab = "fy1xz(y|0.7,0.3)")

#TETA = 1.5, x  = 0.7, z = 0.3
teta <- 1.5; x <- 0.7; z <- 0.3

fy1xz_0.7_0.3_ev1.5 <- mapply(fy1xz, y, rep(x,m), rep(z,m), rep(teta,m))
dev.new()
plot(y,fy1xz_0.7_0.3_ev1.5, type = "l", col = "red", lwd = 2,
    main = "Y|X = 0.7, Z = 0.3 Teta = 1.5", ylab = "fy1xz(y|0.7,0.3)")

##Z|XY

fz1xy <- function(z,x,y,teta){1/(x*y^teta)*((0 < z )&(z < x*y^teta))}

m <- 1000
z <- seq(0.1,0.9, length.out = m)


#TETA = 0, x = 0.5, y = 0.25
teta <- 0; x <- 0.5; y <- 0.25
fz1xy_0.5_0.25_ev0 <- mapply(fz1xy, z, rep(x,m), rep(y,m), rep(teta,m))
dev.new()
plot(z,fz1xy_0.5_0.25_ev0, type = "l", col = "red", lwd = 2,
    main = "(Z|X = 0.5, Y = 0.25).Teta = 0~Unif(0,0.5)", ylab = "fz1xy(z|0.5,0.25)")

#TETA = 1, x = 0.5, y = 0.8
teta <- 1; x <- 0.5; y <- 0.8
fz1xy_0.5_0.8_ev1 <- mapply(fz1xy, z, rep(x,m), rep(y,m), rep(teta,m))
dev.new()
plot(z,fz1xy_0.5_0.8_ev1, type = "l", col = "red", lwd = 2,
    main = "(Z|X = 0.5, Y = 0.8).Teta = 1", ylab = "fz1xy(z|0.5,0.8)")


#TETA = 1.5, x = 0.7, y = 0.8
teta <- 1.5; x <- 0.7; y <- 0.8
fz1xy_0.7_0.8_ev1.5 <- mapply(fz1xy, z, rep(x,m), rep(y,m), rep(teta,m))
dev.new()
plot(z,fz1xy_0.7_0.8_ev1.5, type = "l", col = "red", lwd = 2,
    main = "(Z|X = 0.7, Y = 0.8).Teta = 1.5", ylab = "fz1xy(z|0.7,0.8)")

#============================================================================
#SIMULACIONES

fy1x <- function(y,x) dunif(y,0,x)
invFy1x <- function(v,x) qunif(v,0,x)
#(X,Y)
u <- runif(5000, 0, 1)
x <- qunif(u, 0, 1)
v <- runif(5000, 0, 1)
y <- mapply(invFy1x, v, x)
x.y <- cbind(x,y)
dev.new()
plot(x.y, main = "Simulaciones de (X,Y)", xlab = "x", ylab = "y")

#(X,Y|Z)

#Teta = 0
invFx1z_0 <- function(u,z){z^(1-u)*((0 < u) & (u < 1))}
Fy1xz_0 <- function(y,x,z) {punif(y, 0, x)*(z <= x)}
invFy1xz_0 <- function(v,x,z) {qunif(v, 0, x)*(z <= x)}
z <- 0.25
u <- runif(5000, 0, 1)
x1z_0 <- mapply(invFx1z_0, u, rep(z, 5000))
v <- runif(5000, 0, 1)
y1z_0 <- mapply(invFy1xz_0, v, x1z_0, rep(z, 5000))
dev.new()
plot(x1z_0, y1z_0, main = "Simulaciones de (X,Y|Z = 0.5), teta = 0", xlab = "x|z", ylab = "y|z",
     xlim = c(0,1))
lines(rep(z,1000), seq(0,1,length.out = 1000), col = "red", lwd = 3)
legend(0.3, 0.8, legend = "z", col ="red", lty = 1, cex = 0.8, box.lty = 0)

#Teta = 1
Fx1z_1 <- function(x,z) {((1/(log(z)+2/sqrt(z)-2))*(2/sqrt(z) + (log(z)-2*log(x)-2)/x))*((sqrt(z) <= x)&(x < 1)) + 1*(x >= 1)}
invFx1z_1 <- function(u,z){uniroot(function(x) Fx1z_1(x,z) - u, interval = c(0.00001,0.99999))$root}
Fy1xz_1 <- function(y,x,z) {(((log(y) + log(x) - log(z))/(2*log(x) - log(z)))*((z/x <= y) & (y < x)) + 1*(y >= x))*(z <= x^2)}
invFy1xz_1 <- function(v,x,z){uniroot(function(y) Fy1xz_1(y,x,z) - v, interval = c(0.00001, 0.9999))$root}
z <- 0.25
m <- 2000
u <- runif(m, 0, 1)
x1z_1 <- mapply(invFx1z_1, u, rep(z, m))
v <- runif(m, 0, 1)
y1z_1 <- mapply(invFy1xz_1, v, x1z_1, rep(z, m))
dev.new()
plot(x1z_1, y1z_1, main = "Simulaciones de (X,Y|Z = 0.25), teta = 1", xlab = "x|z", ylab = "y|z",
     xlim = c(0,1), ylim = c(0,1))

#Teta > 0 & teta != 1

Fx1z_t <- function(x,z,t){((z^(-t/(t+1)) - x^(-t) + t^2*z^((1-t)/t)*(x^(-1/t)-z^(-1/(t^2+t))))/(z^(-t/(t+1)) - 1 + t^2*z^((1-t)/t)*(1-z^(-1/(t^2+t)))))*((z^(1/(t+1)) <= x)&(x < 1)) + 1*(x >= 1)}
invFx1z_t <- function(u,z,t){uniroot(function(x) Fx1z_t(x,z,t) - u, interval = c(0.000001,0.99999), tol = .Machine$double.eps^0.25)$root}
Fy1xz_t <- function(y,x,z,t){(((y^(1 - t) - (z/x)^(1/t - 1))/(x^(1 - t) - (z/x)^(1/t - 1)))*(((z/x)^(1/t) < y) & (y < x)) + 1*(y >= x))*(z <= x^{t + 1})}
invFy1xz_t <- function(v,x,z,t){(v*(x^(1 - t) - (z/x)^(1/t - 1)) + (z/x)^(1/t - 1))^(1/(1 - t))}
z <- 0.5
t <- 1.5
m <- 2500
u <- runif(m, 0, 1)
x1z_0.25_0.15 <- mapply(invFx1z_t, u, rep(z,m), rep(t,m))
v <- runif(m, 0, 1)
y1z_0.25_0.15 <- mapply(invFy1xz_t, v, x1z_0.25_0.15, rep(z,m), rep(t,m))
dev.new()
plot(x1z_0.25_0.15,y1z_0.25_0.15, main = "Simulaciones de (X,Y|Z = 0.5), teta = 1.5", xlab = "x|z", ylab = "y|z",
     xlim = c(0,1), ylim = c(0,1))

z <- 0.25
t <- 2
m <- 2500
u <- runif(m, 0, 1)
x1z_0.25_0.15 <- mapply(invFx1z_t, u, rep(z,m), rep(t,m))
v <- runif(m, 0, 1)
y1z_0.25_0.15 <- mapply(invFy1xz_t, v, x1z_0.25_0.15, rep(z,m), rep(t,m))
plot(x1z_0.25_0.15,y1z_0.25_0.15, main = "Simulaciones de (X,Y|Z = 0.25), teta = 2", xlab = "x|z", ylab = "y|z",
     xlim = c(0,1), ylim = c(0,1))

#(X,Z)

#Teta = 0
Fz1x_0 <- function(z,x) {punif(z, 0, x)}
invFz1x_0 <- function(w,x){qunif(w, 0, x)}
m <- 10000
u <- runif(m, 0, 1)
x <- qunif(u, 0, 1)
w <- runif(m, 0, 1)
z <- mapply(invFz1x_0, w, x)
plot(x, z, main = "Simulaciones de (X,Z), teta = 0", xlab = "x", ylab = "z",
     xlim = c(0,1), ylim = c(0,1))

#Teta = 1
Fz1x_1 <- function(z,x) {( (2*z*log(x)-z*log(z)+z)/(x^2))*((0<z)&(z<x^2)) + 1*(x^2 <= z)}
invFz1x_1 <- function(w,x){uniroot(function(z) Fz1x_1(z,x) - w, interval = c(0.00000000001,0.9999999999))$root}
m <- 10000
u <- runif(m, 0, 1)
x <- qunif(u, 0, 1)
w <- runif(m, 0, 1)
z <- mapply(invFz1x_1, w, x)
plot(x, z, main = "Simulaciones de (X,Z), teta = 1", xlab = "x", ylab = "z",
     xlim = c(0,1), ylim = c(0,1))

#Teta > 0 & Teta != 1
Fz1x_t <- function(z,x,t) {((1/(1-t))*(z/x^(t+1) - (t*z^(1/t))/x^(1/t +1)))*((0 <= z)&(z < x^(t+1))) + 1*(x^(t+1)<=z)}
invFz1x_t <- function(w,x,t){uniroot(function(z) Fz1x_t(z,x,t) - w, interval = c(0.0000000000001,0.9999999999))$root}
m <- 10000
t <- 0.1
u <- runif(m, 0, 1)
x <- qunif(u, 0, 1)
w <- runif(m, 0, 1)
z <- mapply(invFz1x_t, w, x, rep(t,m))
plot(x, z, main = "Simulaciones de (X,Z), teta = 0.1", xlab = "x", ylab = "z",
     xlim = c(0,1), ylim = c(0,1))


m <- 10000
t <- 0.5
u <- runif(m, 0, 1)
x <- qunif(u, 0, 1)
w <- runif(m, 0, 1)
z <- mapply(invFz1x_t, w, x, rep(t,m))
plot(x, z, main = "Simulaciones de (X,Z), teta = 0.5", xlab = "x", ylab = "z",
     xlim = c(0,1), ylim = c(0,1))

m <- 10000
t <- 1.5
u <- runif(m, 0, 1)
x_1 <- qunif(u, 0, 1)
w <- runif(m, 0, 1)
z_1 <- mapply(invFz1x_t, w, x_1, rep(t,m))
plot(x_1, z_1, main = "Simulaciones de (X,Z), teta = 1.5", xlab = "x", ylab = "z",
     xlim = c(0,1), ylim = c(0,1), add = TRUE)



#(X,Z|Y)
Fx1y <- function(x,y){(1-log(x)/log(y))*((y <= x)&(x < 1)) + 1*(x >= 1)}
invFx1y <- function(u,y){(y^(1 - u))*((0 <= y) & (y < 1)) + 1*(u == 1)}
Fz1xy <- function(z,x,y,t){dunif(z, 0, x*y^t)}
invFz1xy <- function(w, x, y, t){qunif(w, 0, x*y^t)}
m <- 5000
t <- 0
y <- 0.75
u <- runif(m, 0, 1)
x1y <- mapply(invFx1y, u, rep(y,m))
w <- runif(m, 0, 1)
z1y <- mapply(invFz1xy, w, x1y, rep(y,m), rep(t,m)) 
plot(x1y, z1y, main = "Simulaciones de (X,Z|Y = 0.75), teta = 0", xlab = "x|y", ylab = "z|y",
     xlim = c(0,1), ylim = c(0,1))

t <- 1
y <- 0.25
u <- runif(m, 0, 1)
x1y <- mapply(invFx1y, u, rep(y,m))
w <- runif(m, 0, 1)
z1y <- mapply(invFz1xy, w, x1y, rep(y,m), rep(t,m)) 
plot(x1y, z1y, main = "Simulaciones de (X,Z|Y = 0.25), teta = 1", xlab = "x|y", ylab = "z|y",
     xlim = c(0,1), ylim = c(0,1))

t <- 1
y <- 0.75
u <- runif(m, 0, 1)
x1y <- mapply(invFx1y, u, rep(y,m))
w <- runif(m, 0, 1)
z1y <- mapply(invFz1xy, w, x1y, rep(y,m), rep(t,m)) 
plot(x1y, z1y, main = "Simulaciones de (X,Z|Y = 0.75), teta = 1", xlab = "x|y", ylab = "z|y",
     xlim = c(0,1), ylim = c(0,1))

t <- 0.5
y <- 0.25
u <- runif(m, 0, 1)
x1y <- mapply(invFx1y, u, rep(y,m))
w <- runif(m, 0, 1)
z1y <- mapply(invFz1xy, w, x1y, rep(y,m), rep(t,m)) 
plot(x1y, z1y, main = "Simulaciones de (X,Z|Y = 0.25), teta = 0.5", xlab = "x|y", ylab = "z|y",
     xlim = c(0,1), ylim = c(0,1))
lines(rep(y,1000), seq(0,1,length.out = 1000))
az <- seq(0,1,length.out = 1000)
faux <- function(x) x/2
bz <- sapply(az,faux)
lines(az,bz)
t <- 1.5
y <- 0.3
u <- runif(m, 0, 1)
x1y <- mapply(invFx1y, u, rep(y,m))
w <- runif(m, 0, 1)
z1y <- mapply(invFz1xy, w, x1y, rep(y,m), rep(t,m)) 
plot(x1y, z1y, main = "Simulaciones de (X,Z|Y = 0.25), teta = 0", xlab = "x|y", ylab = "z|y",
     xlim = c(0,1), ylim = c(0,1))

#(Y,Z)
Fy <- function(y) {(y*(1-log(y)))*((0 < y)&(y < 1)) + 1*(y>=1)}
invFy <- function(v) {uniroot(function(y) Fy(y) - v, interval = c(0.00000001,0.9999999))$root}
Fz1y_t <- function(z,y,t) {((-z/log(y))*(1/y^(t+1) - 1/y^t))*((0 <= z)&(z <= y^(t+1))) + (z/(y^t*log(y)) - 1/log(y) - log(z)/log(y) +(t+1))*((y^(t+1) < z)&(z < y^t)) + 1*(z>=y^t)}
invFz1y_t <- function(w,y,t){((log(y)*y^(t + 1)*w)/(y - 1))*((0 <= w) & (w <= (y - 1)/log(y))) + 
(uniroot(function(z) Fz1y_t(z,y,t) - w, interval = c(0.00000000000000001,0.9999999999999))$root)*(((y - 1)/log(y) < w) & (w < (y - 1)/(y*log(y))))}

m <- 5000
t <- 0
v <- runif(m, 0, 1)
y <- mapply(invFy, v)
w <- runif(m, 0, 1)
z <- mapply(invFz1y_t, w, y, rep(t,m))
plot(y, z, main = "Simulaciones de (Y,Z), teta = 0", xlab = "y", ylab = "z",
     xlim = c(0,1), ylim = c(0,1))

m <- 5000
t <- 0.5
v <- runif(m, 0, 1)
y <- mapply(invFy, v)
w <- runif(m, 0, 1)
z <- mapply(invFz1y_t, w, y, rep(t,m))
plot(y, z, main = "Simulaciones de (Y,Z), teta = 0.5", xlab = "y", ylab = "z",
     xlim = c(0,1), ylim = c(0,1))


m <- 5000
t <- 1
v <- runif(m, 0, 1)
y <- mapply(invFy, v)
w <- runif(m, 0, 1)
z <- mapply(invFz1y_t, w, y, rep(t,m))
plot(y, z, main = "Simulaciones de (Y,Z), teta = 1", xlab = "y", ylab = "z",
     xlim = c(0,1), ylim = c(0,1))


m <- 5000
t <- 1.5
v <- runif(m, 0, 1)
y <- mapply(invFy, v)
w <- runif(m, 0, 1)
z <- mapply(invFz1y_t, w, y, rep(t,m))
plot(y, z, main = "Simulaciones de (Y,Z), teta = 1.5", xlab = "y", ylab = "z",
     xlim = c(0,1), ylim = c(0,1))


#(Y,Z|X)
Fy1x <- function(y,x){(y/x)*((0 <= y)&(y < x)) + 1*(y>=x)}
invFy1x <- function(v,x){qunif(v, 0, x)}
Fz1xy <- function(z,x,y,t){dunif(z, 0, x*y^t)}
invFz1xy <- function(w, x, y, t){qunif(w, 0, x*y^t)}

m <- 1000
t <- 0
x <- 0.25
v <- runif(m, 0, 1)
y1x <- mapply(invFy1x, v, rep(x,m))
w <- runif(m, 0, 1)
z1x <- mapply(invFz1xy, w, rep(x,m), y1x, rep(t,m))
plot(y1x,z1x, main = "Simulaciones de (Y,Z|X = 0.25), teta = 0", xlab = "y", ylab = "z",
     xlim = c(0,1), ylim = c(0,1))

m <- 1000
t <- 0.5
x <- 0.45
v <- runif(m, 0, 1)
y1x <- mapply(invFy1x, v, rep(x,m))
w <- runif(m, 0, 1)
z1x <- mapply(invFz1xy, w, rep(x,m), y1x, rep(t,m))
plot(y1x,z1x, main = "Simulaciones de (Y,Z|X = 0.45), teta = 0.5", xlab = "y", ylab = "z",
     xlim = c(0,1), ylim = c(0,1))

m <- 3000
t <- 0.5
x <- 0.75
v <- runif(m, 0, 1)
y1x <- mapply(invFy1x, v, rep(x,m))
w <- runif(m, 0, 1)
z1x <- mapply(invFz1xy, w, rep(x,m), y1x, rep(t,m))
plot(y1x,z1x, main = "Simulaciones de (Y,Z|X = 0.75), teta = 0.5", xlab = "y", ylab = "z",
     xlim = c(0,1), ylim = c(0,1))

m <- 1000
t <- 1
x <- 0.45
v <- runif(m, 0, 1)
y1x <- mapply(invFy1x, v, rep(x,m))
w <- runif(m, 0, 1)
z1x <- mapply(invFz1xy, w, rep(x,m), y1x, rep(t,m))
plot(y1x,z1x, main = "Simulaciones de (Y,Z|X = 0.45), teta = 1", xlab = "y", ylab = "z",
     xlim = c(0,1), ylim = c(0,1))

m <- 3000
t <- 1
x <- 0.75
v <- runif(m, 0, 1)
y1x <- mapply(invFy1x, v, rep(x,m))
w <- runif(m, 0, 1)
z1x <- mapply(invFz1xy, w, rep(x,m), y1x, rep(t,m))
plot(y1x,z1x, main = "Simulaciones de (Y,Z|X = 0.75), teta = 1", xlab = "y", ylab = "z",
     xlim = c(0,1), ylim = c(0,1))

m <- 1000
t <- 1.5
x <- 0.45
v <- runif(m, 0, 1)
y1x <- mapply(invFy1x, v, rep(x,m))
w <- runif(m, 0, 1)
z1x <- mapply(invFz1xy, w, rep(x,m), y1x, rep(t,m))
plot(y1x,z1x, main = "Simulaciones de (Y,Z|X = 0.45), teta = 1.5", xlab = "y", ylab = "z",
     xlim = c(0,1), ylim = c(0,1))