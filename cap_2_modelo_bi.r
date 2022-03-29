                           ###EJERCICIO 1

#PAQUETES REQUERIDOS.
library(subcopem2D)
library(copula)
library(ADGofTest)
#==========================================================================================
#Desarrollo del modelo te�rico bivariado...

Psi_x <- function(y,x) {
         ((1-exp(-x)-x*exp(-y))/(1-exp(-x)))*(x < y) + 
         ((1-exp(-y)-y*exp(-y))/(1-exp(-x)))*((0 < y)&(y <= x))}

Fi_y <- function(x,y) {
        ((1-exp(-x)-x*exp(-y))/(1-exp(-y)*(y+1)))*((0 < x )&(x < y)) + 1*(y <= x)}

derivPsi_x <- function(y,x) {
              ((x*exp(-y))/(1-exp(-x)) )*(x < y) + ( (y*exp(-y))/(1-exp(-x)) )*(
              (0 < x)&(y <= x))}
derivFi_y <- function(x,y) {
             (( exp(-y)*(x-y) )/( 1-exp(-y)*(y+1) ))*((0 < x)&(x < y)) + 0*(y <= x)}

y <- seq(0.0001, 10, length.out = 1000)
dev.new()
plot(y, Psi_x(y,1), type = "l", lwd = 2, main = "Psi_x(y)", ylab = "Psi_x(y)")
lines(y, Psi_x(y,1.5), col = "brown", lwd = 2)
lines(y, Psi_x(y,2), col = "green", lwd = 2)
lines(y, Psi_x(y,5), col = "red", lwd = 2)
lines(y, Psi_x(y,10), col = "blue", lwd = 2)
legend(6, 0.8, legend = c("Psi_x(y)|x = 1", "Psi_x(y)|x = 1.5", "Psi_x(y)|x = 2",
       "Psi_x(y)|x = 5", "Psi_x(y)|x = 10"), col = c("black", "brown", "green", 
       "red", "blue"), lty = 1, cex = 0.8)
dev.new()
plot(y, derivPsi_x(y,1), type = "l", lwd = 2, main = "Derivada de Psi_x(y)", 
     ylab = "Psi_x(y)")
lines(y, derivPsi_x(y,1.5), col = "brown", lwd = 2)
lines(y, derivPsi_x(y,2), col = "green", lwd = 2)
lines(y, derivPsi_x(y,5), col = "red", lwd = 2)
lines(y, derivPsi_x(y,10), col = "blue", lwd = 2)
legend(6, 0.4, legend = c("d/dy(Psi_x(y)) |x = 1", "d/dy(Psi_x(y)) |x = 1.5", 
       "d/dy(Psi_x(y)) |x = 2", "d/dy(Psi_x(y)) |x = 5", "d/dy(Psi_x(y)) |x = 10"),
       col = c("black", "brown", "green", "red", "blue"), lty = 1, cex = 0.8)


x <- seq(0.001, 10, length.out = 1000)
dev.new()
plot(x, Fi_y(x,1), type = "l", lwd = 2, main = "Fi_y(x)", ylab = "Fi_y(x)")
lines(x, Fi_y(x,1.5), col = "brown", lwd = 2)
lines(x, Fi_y(x,2), col = "green", lwd = 2)
lines(x, Fi_y(x,5), col = "red", lwd = 2)
lines(x, Fi_y(x,10), col = "blue", lwd = 2)
legend(6, 0.8, legend = c("Fi_y(x)|y = 1", "Fi_y(x)|y = 1.5", "Fi_y(x)|y = 2",
       "Fi_y(x)|y = 5", "Fi_y(x)|y = 10"), col = c("black", "brown", "green", 
       "red", "blue"), lty = 1, cex = 0.8)

dev.new()
plot(x, derivFi_y(x,1), type = "l", lwd = 2, main = "Derivada de Fi_y(x)", 
     ylab = "(Fi_y(x))'")
lines(x, derivFi_y(x,1.5), col = "brown", lwd = 2)
lines(x, derivFi_y(x,2), col = "green", lwd = 2)
lines(x, derivFi_y(x,5), col = "red", lwd = 2)
lines(x, derivFi_y(x,10), col = "blue", lwd = 2)
legend(6, -0.4, legend = c("d/dy(Fi_y(x)) |y = 1", "d/dy(Fi_y(x)) |y = 1.5", 
       "d/dy(Fi_y(x)) |y = 2", "d/dy(Fi_y(x)) |y = 5", "d/dy(Fi_y(x)) |y = 10"),
       col = c("black", "brown", "green", "red", "blue"), lty = 1, cex = 0.8)


#==========================================================================================

# FUNCIONES EXTRA.
# EN ESTA PATE �NICAMENTE DEFINO ALGUNAS FUNCIONES QUE VOY A UTLIZAR MAS ADELANTE.

# FUNCI�N DE DISTRIBUCI�N EMP�RICA.
Fn <- function(x, muestra) sapply(x, function(z) mean(muestra <= z))

# FUNCI�N EMP�RICA DE CUANTILES.
cuantil.ep <- function(p, muestra){
  aviso <- ""
  n <- length(muestra)
  if (p < 1/(n+1)){
    cuantil <- min(muestra)
    aviso <- "cuantil fuera de rango muestral"
  }
  if (p > n/(n+1)){
    cuantil <- max(muestra)
    aviso <- "cuantil fuera de rango muestral"
  }
  if ((1/(n+1) <= p)&(p <= n/(n+1))){
    if(((n+1)*p)%in%(1:n)){
      j <- as.integer((n+1)*p)
      cuantil <- sort(muestra)[j]
    }
    else {
      j <- floor((n+1)*p)
      x.ord <- sort(muestra)[j:(j+1)]
      cuantil <- (j+1-(n+1)*p)*x.ord[1] + ((n+1)*p - j)*x.ord[2]
    }
  }
  if (aviso!="") warning(aviso)
  return(cuantil)
}
cuantil <- function(p, x.obs) sapply(p, cuantil.ep, muestra = x.obs)

#FUNCI�N SIMILAR A OUTER.
rat <- function(x,y,f) { 
 mat <- matrix( ,length(x), length(y))
   for(i in 1:length(x)) { 
     a <- rep(x[i], length(x))
     mat[ i, ] <- mapply(f, a, y)
   }
 mat
}
#==========================================================================================
# MODELO TE�RICO.
# EN ESTA PARTE PROGRAMO LAS FUNCIONES RELACIONADAS AL MODELO TE�RICO DEL EJERCICIO

# DISTRIBUCIONES MARGINALES.
Fx <- function(x) { pexp(x, 1)}
Fy <- function(y) { pgamma(y, 2, 1) }

# INVERSAS DISTRIBUCIONES MARGINALES. 
invFx <- function(u) { qexp(u, 1) }
invFy <- function(v) { qgamma(v, 2, 1) }

# P(Y <= y | X = x) Y SU INVERSA.
Fylx <- function(y,x) { 1- exp(x-y) }
invFylx <- function(z,x) { x - log(1-z) }

# DENSIDAD CONJUNTA
fxy <- function(x,y) { exp(-y)*(0<x & x<y) }

# DISTRIBUCI�N CONJUNTA.
Fxy <- function(x,y) {
       (1-exp(-x)-x*exp(-y))*(0<x & x<y) + (1-exp(-y)-y*exp(-y))*(0<y & y<=x) }

# AQU� PROGRAM� LA C�PULA SUBYACENTE V�A COROLARIO DE SKLAR.
# DESPU�S SIMPLEMENTE REALIC� UN GR�FICO DE LAS SUPERFICIES, CORROBORANDO LO QUE 
# DICE EL COROLARIO DE SKLAR 

#C�PULA TE�RICA.
m <- 51
C_XY <- function(u,v) {
        (u + (log(1-u))*exp(-invFy(v)))*(0 < u & u < 1 - exp(-invFy(v))) + 
         v*( (0 < 1-exp(-invFy(v)))&(1-exp(-invFy(v)) <= u))
  }
us <- seq(0, 1, length.out = m)
vs <- seq(0, 1, length.out = m)
z <- outer(us,vs,C_XY)
zf <- outer(invFx(us), invFy(vs), Fxy)
dev.new()
image(us,vs,z, col = heat.colors(20), main = "C�pula te�rica",
      xlab = "u", ylab ="v")


# DENSIDAD DE LA C�PULA TE�RICA. 
densC <- function(u,v) {(1/(invFy(v)*(1-u)))*((0 < u)&(u < 1-exp(-invFy(v))))}
ud <- seq(0.01, 0.95, length.out = 100)
vd <- seq(0.01, 0.95, length.out = 100)
zd <- outer(ud,vd,densC)
image(ud, vd, zd, col = heat.colors(20), main = "Densidad de la c�pula", xlab = "u", 
      ylab = "v")
#=========================================================================================

# SIMULACIONES DEL VECTOR (X,Y).
# LAS SIGUIENTES LINEAS SON PARA GENERAR UNA MUESTRA ALEATORIA CON LA DISTRIBUCI�N
# DEL EJERCICIO.
m <- 10000
u <- runif(m, 0, 1)
x <- invFx(u)
v <- runif(m, 0, 1)
y <- mapply(invFylx,v, x)

# EN EL SIGUIENTE GR�FICO APARECEN LA DENSIDAD DE PUNTOS (X,Y) Y A LADO LA DENSIDAD
# CONJUNTA DEL MODELO.
dev.new()
par(mfrow = c(1,2))
plot(x, y, main = "Observaciones del vector aleatorio (X,Y)")
xd <- seq(0, max(max(x), max(y)), length.out = 50)
yd <- seq(0, max(max(x), max(y)), length.out = 50)
zd <- outer(xd, yd, fxy)
zd[lower.tri(zd, diag = F)] <- NA
dev.new()
persp(xd, yd, zd, theta = 180, phi = 25, main = "Densidad conjunta de (X,Y)", col = "gold", 
      xlab = "x", ylab = "y", zlab = "f(x,y)")
filled.contour(xd, yd, zd, color = function(n) hcl.colors(n, "Oranges"), 
               plot.title = title("Densidad conjunta de (X,Y)", xlab = "x", ylab = "y"))

##ESTE GR�FICO REPRESENTA LOS HISTOGRAMAS MARGINALES DE X, Y. CORROBORANDO QUE 
#COINCIDAN CON LAS TE�RICAS.
dev.new()
par(mfrow = c(1,2))
hist(x, prob = T, main = "Marginal X", col = "yellow")
xi <- seq(0, max(x) + sd(x), length.out = 1000)
lines(xi, dexp(xi), col = "red", lwd = 2) 
hist(y, prob = T, main = "Marginal Y", col = "green", ylim = c(0,0.4))
yi <- seq(0, max(y) + sd(y), length.out = 1000)
lines(yi, dgamma(yi, 2, 1), col = "red", lwd = 2)
ad.test(x, pexp, 1)
ad.test(y, pgamma, 2, 1)

#DEFINIMOS LA MUESTRA CONJUNTA.
ma.obs <- cbind(x,y)
head(ma.obs)
#==========================================================================================

##AN�LISIS DE LAS OBSERVACIONES.
#ACONTINUACI�N RELIZO UN POCO DE AN�LISIS DE LA MUESTRA CONJUNTA

matUV <- apply(ma.obs, 2, rank)/nrow(ma.obs) #Pseudo-observaciones
dev.new()
par(mfrow = c(1,2))
plot(matUV, main = "Pseudo-observaciones", xlab = "u", ylab = "v")
plot(ma.obs, main = "Observaciones del vector (X,Y)") 
copemp <- subcopemc(ma.obs, m = 50)
(Pearson <- cor(ma.obs))
(Spearman <- cor(ma.obs, method = "spearman"))
(Kendal <- cor(ma.obs, method = "kendal"))
(depMono <- copemp$depMon)
(depSup <- copemp$depSup)
#==========================================================================================

##COMPARACI�N; C�PULA TE�RICA VS C�PULA EMP�RICA.
#EN ESTA PARTE HICE LA COMPARATIVA DE LAS CURVAS DE NIVEL DE LA C�PULA
#TE�RICA (V�A COROLARIO DE SKLAR) VS LA C�PULA EMPIRICA (LA QUE RESULTA DEL 
#PAQUETE SUBCOPEM2D).

vcopemp <- copemp$matrix
ucop <- seq(0, 1, length.out = 51)
vcop <- seq(0, 1, length.out = 51)
dev.new()
contour(us, vs, z, col = "red", main = "C�pula te�rica vs C�pula emp�rica")
contour(ucop, vcop, vcopemp, add = T, col = "blue")
legend("bottomright", legend = c("C�pula te�rica", "C�pula emp�rica"), 
       fill = c("red", "blue"))
#==========================================================================================

## COMPARACI�N; C�PULA TE�RICA VS C�PULA "CONDICIONAL (M�TODO 1, 2 Y 3)".

## M�TODO 1 (Fxy(x,y) = P(Y <= y | X <= x)*Fx(x)).

#  P(Y <= y | X <= x) MUESTRAL.

Aylx <- function(y,x,m) {
    cond <- m[which(m[ ,1] <= x), ]
 if(nrow(cond) == 0 && class(cond) == "matrix"){
  (Aylx <- 0)
 } 
 else if(class(cond) == "numeric") { 
  as.numeric(cond[2] <= y)
 }
 else {
  length(which(cond[ ,2] <= y))/nrow(cond)
 }
} 

## Fxy DE LA FORMA P(Y <= y | X <= x)*Fx(x).
Fxy_1 <- function(a,b) { Aylx(b, a, ma.obs)*Fn(a,x) }

# CON BASE EN EL COROLARIO DE SKLAR Y EL CAMBIO DE VARIABLE u = Fx(x), v = Fy(y)
# NOS QUEDA QUE LA C�PULA RESULTANTE POR ESA V�A QUEDA COMO: 
# C(u,v) = u*P(Y <= Fy^-1(v) | X <= Fx^-1(u))

copcond <- function(u,v) { u*Aylx(cuantil(v,y), cuantil(u,x), ma.obs) }

# LO SIGUIENTE QUE HAGO ES REALIZAR UNA COMPARACI�N ENTRE LAS CURVAS DE NIVEL DE 
# LA C�PULA TE�RICA VS LA C�PULA QUE RESULTA POR ESTE M�TODO.

zcc <- rat(ucop,vcop,copcond)
dev.new()
contour(us, vs, z, col = "red", main = "C�pula te�rica vs C�pula condicional 1")
contour(ucop, vcop, zcc, col = "blue", add = T)
legend("bottomright", legend = c("C�pula te�rica", "C�pula condicional 1"), 
       fill = c("red", "blue"))

# M�TODO 2 (Fxy(x,y) = Fy(y) - (1-Fx(x))*P(Y <= y | X > x))
# LA SIGUIENTE FUNCI�N REPRESENTA P(Y <= y | X > x) MUESTRAL.

Mylx <- function(y,x,m) { 
   cond1 <- m[which(m[,1] > x), ]
 if(nrow(cond1) == 0 && class(cond1) == "matrix" ){
  (Mylx <- 0)
 } 
 else if(class(cond1) == "numeric") { 
  as.numeric(cond1[2] <= y)
 }
 else {
  length(which(cond1[ ,2] <= y))/nrow(cond1)
 }
} 

##PROGRAMO Fxy DE LA FORMA Fy(y) - (1-Fx(x))*P(Y <= y | X > x)
Fxy_2 <- function(a,b) { Fn(b,y) - (1 - Fn(a,x))*Mylx(b, a, ma.obs) }

# CON BASE EN EL COROLARIO DE SKLAR Y EL CAMBIO DE VARIABLE u = Fx(x), v = Fy(y)
# NOS QUEDA QUE LA C�PULA RESULTANTE POR ESA V�A QUEDA COMO: 
# C(u,v) = v - (1-u)*P(Y <= Fy^-1(v) | X > Fx^-1(u))

copcond2 <- function(u,v) { v - (1 - u)*Mylx(cuantil(v,y), cuantil(u,x), ma.obs) }

# LO SIGUIENTE QUE HAGO ES REALIZAR UNA COMPARACI�N ENTRE LAS CURVAS DE NIVEL DE 
# LA C�PULA TE�RICA VS LA C�PULA QUE RESULTA POR ESTE M�TODO.

zcc2 <- rat(ucop,vcop,copcond2)
dev.new()
contour(us, vs, z, col = "red", main = "C�pula te�rica VS C�pula condicional 2")
contour(ucop, vcop, zcc2, col = "blue", add = T)
legend("bottomright", legend = c("C�pula te�rica", "C�pula condicional 2"), 
       fill = c("red", "blue"))


# M�TODO 3, COMB. LINEAL CONVEXA
# EN ESTE M�TODO UNO UNA COMBINACI�N LINEAL CONVEXA DE AMBOS M�TODOS, CON EL 
# PROPOSITO DE TRATAR DE USAR MAS DATOS PARA LOS C�LCULOS. EL PAR�METRO PARA
# LA COMBINACI�N LINEAL ALFA = Fx(x)

Fxy_3 <- function(a,b) { Fn(a,x)*Fxy_1(a,b) + (1-Fn(a,x))*Fxy_2(a,b) }

# V�A COROLARIO DE SKLAR, LA C�PULA RESULTANTE ES LA SIGUIENTE:
copcond3 <- function(u,v) { (u)*copcond(u, v) + (1-u)*copcond2(u,v) }

# LO SIGUIENTE QUE HAGO ES REALIZAR UNA COMPARACI�N ENTRE LAS CURVAS DE NIVEL DE 
# LA C�PULA TE�RICA VS LA C�PULA QUE RESULTA POR ESTE M�TODO.

zcc3 <- rat(ucop,vcop,copcond3)
dev.new()
contour(us, vs, z, col = "red", main = "C�pula te�rica vs C�pula condicional 3")
contour(ucop, vcop, zcc3, col = "blue", add = T)
legend("bottomright", legend = c("C�pula te�rica", "C�pula condicional 3"), 
       fill = c("red", "blue"))


#==========================================================================================

# GR�FICAS C�PULA EMP�RICA VS C�PULA "CONDICIONAL"
# EN ESTA PARTE SE GENERAN LAS GR�FICAS DE LAS CURVAS DE NIVEL DE LA C�PULA EMP�RICA
# VS LAS C�PULAS GENERADAS CON LOS 3 M�TODOS PROPUESTOS.
# M�TODO 1

dev.new()
contour(ucop, vcop, vcopemp, col = "red", main = "C�pula emp�rica vs C�pula condicional 1")
contour(ucop, vcop, zcc, col = "blue", add = T)
legend("bottomright", legend = c("C�pula emp�rica", "C�pula condicional 1"), 
       fill = c("red", "blue"))

#M�TODO 2

dev.new()
contour(ucop, vcop, vcopemp, col = "red", main = "C�pula emp�rica vs C�pula condicional 2")
contour(ucop, vcop, zcc2, col = "blue", add = T)
legend("bottomright", legend = c("C�pula emp�rica", "C�pula condicional 2"), 
       fill = c("red", "blue"))

#M�TODO 3

dev.new()
contour(ucop, vcop, vcopemp, col = "red", main = "C�pula emp�rica vs C�pula condicional 3")
contour(ucop, vcop, zcc3, col = "blue", add = T)
legend("bottomright", legend = c("C�pula emp�rica", "C�pula condicional 3"), 
       fill = c("red", "blue"))

#==========================================================================================
# C�LCULO EMP�RICO DE LA DISTANCIA ENTRE PAR DE C�PULAS.

# SE DEFINE EL SIGUIENTE CONJUNTO DONDE SER�N EVALUADAS LA C�PULA TE�RICA Y LAS 
# C�PULAS POR LOS 3 M�TODOS. DADO QUE LA C�PULA TE�RICA RESULT� TENER DENTRO DE SU
# ESTRUCTURA EL T�RMINO "1/1-u", PRESENTA PROBLEMAS SI UN ELEMENTO DEL CONJUNTO ES 
# 1, RAZ�N POR LA CUAL SE DEFINIERON DE ESA MANERA LOS CONJUNTOS. 

ucop1 <- seq(0, 0.9999999, length.out = 51)
vcop1 <- seq(0, 0.9999999, length.out = 51)

####dist(TE�RICA,COND1)

#EVALUACIONES DE LA C�PULA TE�RICA
zt <- outer(ucop1, vcop1, C_XY)

#EVALUACI�N DE LA C�PULA COND1
zcc <- rat(ucop1, vcop1, copcond)
(dtc1 <- sum(abs(zt - zcc))/(51^2))

####d(TE�RICA, COND2)

#EVALUACI�N DE LA C�PULA COND2
zcc2 <- rat(ucop1, vcop1, copcond2)
(dtc2 <- sum(abs(zt - zcc2))/(51^2))

####d(TE�RICA, COND3)

#EVALUACI�N DE LA C�PULA COND3
zcc3 <- rat(ucop1, vcop1, copcond3)
(dtc3 <- sum(abs(zt - zcc3))/(51^2))

####d(EMP�RICA, COND1)

#EVALUACIONES DE LA C�PULA EMP�RICA
vcopemp <- copemp$matrix
#EVALUACI�N DE LA C�PULA COND1
zcc <- rat(ucop,vcop,copcond)
(dec1 <- sum(abs(vcopemp - zcc))/(51^2))


####d(EMP�RICA, COND2)
#EVALUACI�N DE LA C�PULA COND2
zcc2 <- rat(ucop,vcop,copcond2)
(dec2 <- sum(abs(vcopemp - zcc2))/(51^2))


####d(EMP�RICA, COND3)
#EVALUACI�N DE LA C�PULA COND3
zcc3 <- rat(ucop,vcop,copcond3)
(dec3 <- sum(abs(vcopemp - zcc3))/(51^2))
#==========================================================================================

#AUTOMATIZACI�N DEL PROCESO ANTERIOR PARA 1000 VECES ESTE PROCESO

ucop1 <- seq(0, 0.9999999, length.out = 51)
vcop1 <- seq(0, 0.9999999, length.out = 51)
ucop <- seq(0, 1, by = 1/50)
vcop <- seq(0, 1, by = 1/50)
mat <- matrix( , 10000, 6)
colnames(mat) <- c("d(teo,cond1)", "d(teo,cond2)", "d(teo,cond3)", "d(emp,cond1)", 
                   "d(emp,cond2)", "d(emp,cond3)")  

system.time(
for(k in 1:1000) {
u <- runif(10000, 0, 1)
x <- invFx(u)
v <- runif(10000, 0, 1)
y <- mapply(invFylx,v, x)
ma.obs <- cbind(x,y)
zt <- outer(ucop1, vcop1, C_XY)
zcc2 <- rat(ucop1, vcop1, copcond2)
zcc <- rat(ucop1, vcop1, copcond)
zcc3 <- rat(ucop1, vcop1, copcond3)
mat[k, 1] <- sum(abs(zt - zcc))/(50^2) 
mat[k, 2] <- sum(abs(zt - zcc2))/(50^2)
mat[k, 3] <- sum(abs(zt - zcc3))/(50^2)
}) 

   

system.time(
for(k in 1:1000) {
u <- runif(10000, 0, 1)
x <- invFx(u)
v <- runif(10000, 0, 1)
y <- mapply(invFylx,v, x)
ma.obs <- cbind(x,y)
copemp <- subcopemc(ma.obs, m = 50)
vcopemp <- copemp$matrix
zcc <- rat(ucop, vcop, copcond)
zcc2 <- rat(ucop, vcop, copcond2)
zcc3 <- rat(ucop, vcop, copcond3)
mat[k, 4] <- sum(abs(vcopemp - zcc))/(50^2) 
mat[k, 5] <- sum(abs(vcopemp - zcc2))/(50^2)
mat[k, 6] <- sum(abs(vcopemp - zcc3))/(50^2)
}) 
 