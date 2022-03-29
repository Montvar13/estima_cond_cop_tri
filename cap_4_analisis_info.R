library(fitdistrplus)
library(ADGofTest)
datos <- read.csv(file = file.choose())
str(datos)
colnames(datos)[3] <- "municipio"
str(datos)
n <- nrow(datos)

Pnat <- datos$Pnat
idh <- datos$IDH
homP1000 <- datos$HOMP1000

#================================================================================
#=================== VISUALIZACIÓN CONJUNTA DE LOS DATOS ========================
#================================================================================
 
library(scatterplot3d)
dev.new()
scatterplot3d(datos[ ,4:6], angle = 10, color = "steelblue", pch = 1,
              main = "Datos observados",
              xlab = "%Nativos", ylab = "HomP1000", zlab = "IDH",
              zlim = c(0,1), xlim = c(0,1))

dev.new()
scatterplot3d(datos[ ,c(5,6,4)], angle = 40, color = "steelblue", pch = 1,
              main = "Datos observados",
              xlab = "HomP1000", ylab = "IDH", zlab = "%Nativos",
              zlim = c(0,1), ylim = c(0,1))

dev.new()
scatterplot3d(datos[ ,c(5,4,6)], angle = 50, color = "steelblue", pch = 1,
              main = "Datos observados",
              xlab = "HomP1000", ylab = "%Nativos", zlab = "IDH",
              zlim = c(0,1), ylim = c(0,1))

dev.new()
scatterplot3d(datos[ ,c(6,4,5)], angle = 20, color = "steelblue", pch = 1,
              main = "Datos observados",
              xlab = "IDH", ylab = "%Nativos", zlab = "HomP1000",
              xlim = c(0,1), ylim = c(0,1))

dev.new()
scatterplot3d(datos[ ,c(6,5,4)], angle = 20, color = "steelblue", pch = 1,
              main = "Datos observados",
              xlab = "IDH", ylab = "HomP1000", zlab = "%Nativos",
              zlim = c(0,1), xlim = c(0,1))


#============================================================================================

####ANÁLISIS MARGINAL DE LA VARIABLE "PORCENTAJE DE NATIVOS".
str(Pnat)
sum(duplicated(Pnat))
dev.new()
hist(Pnat, prob = T, col = "yellow", main = "Porcentaje de nativos", 
     xlab = "%", xlim = c(0,1), breaks = 15)
summary(Pnat)
#-----------------------------------------------------------------------------
#BETA

para_beta <- fitdist(Pnat, "beta")$estimate
sims <- rbeta(length(Pnat), para_beta[1], para_beta[2])
summary(sims)
dev.new()
hist(sims, prob = TRUE, xlim = c(0,1), col = "yellow")
ad.test(Pnat, pbeta, para_beta[1], para_beta[2])

#MODELACIÓN CON DISTRIBUCIÓN KUMARASWAMY
library(extraDistr)
dkum <- function(x, a, b) {a*b*x^{a-1}*(1 - x^{a})^{b - 1}}
pkum <- function(q, a, b) {1 - (1 - q^{a})^{b}}
qkum <- function(p, a, b) {(1 - (1 - p)^{1/b})^{1/a}}
(para_kumar <- fitdist(Pnat, "kum", start=list(a=0.1, b=0.1))$estimate)
sims <- rkumar(1905, para_kumar[1], para_kumar[2])
dev.new()
hist(sims, prob = TRUE, col = "green",
     main = "Simulaciones de una Kumaraswamy", xlim = c(0,1))
summary(sims)
ad.test(Pnat, pkumar, para_kumar[1], para_kumar[2])
################################################################################
################################################################################
################################################################################
################################################################################
#MODELACIÓN PARETO AJUSTADA
dpare <- function(x, b) {b*(1 - x)^{b-1}}
ppare <- function(q, b) {1 - (1 - q)^{b}}
qpare <- function(p, b) {1 - (1 - p)^{1/b}}
rpare <- function(n, b) {
         u <- runif(n)
         rpare <- mapply(qpare, u, rep(b, n))
         return(rpare)}
(para_pare <- fitdist(Pnat, "pare", start=list( b=0.1))$estimate)
sims <- rpare(1905, para_pare)
dev.new()
hist(sims, prob = TRUE, col = "green",
     main = "Simulaciones de una Pareto ajustada", xlim = c(0,1))
summary(sims)
ad.test(Pnat, ppare, para_pare)

#------------------------------------------------------------------------------
##SUAVUZAMIENTO POR POLINOMIOS DE BERNSTEIN
# Función de distribución empírica
Fn <- function(x, muestra) sapply(x, function(z) mean(muestra <= z)) 
Fn.inv.Bernstein <- function(u, muestra){
#
# Entrada:
#                 u = probabilidad para calcular cuantil
#           muestra = vector de obervaciones univariadas x1,...,xn
#
    x <- sort(muestra)
    n <- length(x)
    xm <- rep(0, n + 1)
    xm[1] <- x[1]
    xm[n + 1] <- x[n]
    xm[2:n] <- (x[1:(n-1)] + x[2:n])/2
    return(sum(xm * dbinom(0:n, n, u)))
}
Fn.Bernstein.aux <- function(u, xdada, muestra) Fn.inv.Bernstein(u, muestra) - xdada
Fn.Bernstein <- function(x, muestra) uniroot(Fn.Bernstein.aux, interval = c(0, 1),
                         xdada = x, muestra = muestra, tol = 0.00001)$root

sims <- sapply(runif(n), Fn.inv.Bernstein, Pnat)
par(mfrow = c(1,2))
hist(Pnat, prob = T, col = "yellow", main = "Porcentaje de nativos", 
     xlab = "%", xlim = c(0,1), breaks = 15)
hist(sims, prob = TRUE, col = "yellow", main = "Simulación vía Bernstein", 
     breaks = 15, xlab = "%", xlim = c(0,1))
summary(Pnat)
summary(sims)
x <- seq(min(Pnat), max(Pnat), length.out = 10000)
dist <- sapply(x, Fn.Bernstein, Pnat) 
dev.new()
plot(x, dist, type = "l", lwd = 3, main = "Función de distribución", xlim = c(0,1))
lines(seq(max(Pnat), 1, length.out = 1000), rep(1, 1000), lwd = 3)

#============================
z.t <- (log(Pnat - min(Pnat) + 1))^(1)
z.t[which(z.t == min(z.t))] <- 0.00001
dev.new()
hist(z.t, prob = TRUE, xlim = c(0,1))

#MODELACIÓN PARETO AJUSTADA
dpare <- function(x, b) {b*(1 - x)^{b-1}}
ppare <- function(q, b) {1 - (1 - q)^{b}}
qpare <- function(p, b) {1 - (1 - p)^{1/b}}
rpare <- function(n, b) {
         u <- runif(n)
         rpare <- mapply(qpare, u, rep(b, n))
         return(rpare)}
(param_z.t.p <- fitdist(z.t, "pare", start=list( b=0.1))$estimate)
sims_z.t <- rpare(n, param_z.t.p)
dev.new(); par(mfrow = c(1,2))
hist(sims_z.t, prob = TRUE, col = "green",
     main = "Simulaciones de una Pareto ajustada", xlim = c(0,1))
hist(z.t, prob = TRUE, xlim = c(0,1))
summary(sims_z.t)
summary(z.t)
ad.test(z.t, ppare, param_z.t.p)

#=======================================================================================================

 
####ANÁLISIS MARGINAL DE LA VARIABLE "IDH".
str(idh)
sum(duplicated(idh))
###DEBIDO A LA NATURALEZA DE LA VARIABLE ALEATORIA, SE CONSIDERARÁ DICHA VARIABLE
###COMO CONTINUA Y A LOS VALORES REPETIDOS SE LES APLICARÁ JITTER.
idh <- jitter(idh)
summary(idh)
hist(idh, prob = TRUE, main = "Histograma del IDH", col = "green", xlim = c(0,1))
param_idh <- fitdist(idh, "beta")$estimate
sims_idh <- rbeta(n, param_idh[1], param_idh[2])
summary(sims_idh)
x <- seq(0,1, length.out = 10000)
y <- mapply(dbeta, x, 28.12193, 14.18486)
dev.new(); par(mfrow = c(1,2))
hist(idh, prob = TRUE, main = "Histograma del IDH", col = "green", xlim = c(0,1))
lines(x,y, col = "red", lwd = 2)
hist(sims_idh, prob = TRUE, main = "Simulaciones del IDH", xlim = c(0,1),
    col = "green", xlab = "sims IDH")
lines(x,y, col = "red", lwd = 2)
ad.test(idh, pbeta, param_idh[1], param_idh[2])

#===================================================================================

####ANÁLISIS MARGINAL DE LA VARIABLE "homP1000".
str(homP1000)
(num.rep <- sum(duplicated(homP1000)))
#EL % DE DATOS REPETIDOS ES DE: 
porc_rep <- 100*num.rep/n
cat(porc_rep, "%")

pos_rep <- which(duplicated(homP1000) == TRUE)
val_rep <- homP1000[pos_rep]
table(val_rep)
num_rep_0 <- length(which(val_rep == 0))
num_rep_no0 <- length(which(val_rep != 0))
#% de registros con valor 0.
num_rep_0/n

###APROXIMADAMENTE 1/3 DE LOS REGISTROS TIENEN EL VALOR 0, RAZÓN POR LO CUAL, LA
###VARIABLE "homP1000" TENDRÁ MASA EN 0. SÓLO 5 VALORES DISTINTOS DE 0 SE REPITEN
###POR LO QUE LA VARIABLE "homP1000" SERÁ MIXTA, CON MASA EN 0 Y CONTINUA EN EL 
###INTERVALO (0,1].

(masa_0 <- length(which(val_rep == 0))/n)
reg_0 <- homP1000[which(homP1000 == 0)]
reg_NO_0 <- homP1000[-which(homP1000 == 0)]
sum(duplicated(reg_NO_0))
pos_rep_reg_NO_0 <- which(duplicated(reg_NO_0) == T)
dat_rep_jit <- jitter(reg_NO_0[pos_rep_reg_NO_0])
reg_NO_0[pos_rep_reg_NO_0] <- dat_rep_jit
sum(duplicated(reg_NO_0))
summary(reg_NO_0)
dev.new()
hist(reg_NO_0, prob = T, col = "blue",
     xlab = "Homocidios", main = "Homocidios por 1000 habitantes")

homP1000[which(homP1000 != 0)] <-  reg_NO_0
pos_rep <- which(duplicated(homP1000) == TRUE)
val_rep <- homP1000[pos_rep]
table(val_rep)
summary(homP1000)

###TEOREMA DE DESCOMPOSICIÓN
##F_{HomP1000}(x) = alfa*F_{D}(x) + (1 - alfa)*F_{C}(x)
F.D <- function(x) {1*(x >= 0)}
#Haría falta modelar la parte continua de la variable.
summary(reg_NO_0)
dev.new()
hist(reg_NO_0, prob = T, col = "blue",
     xlab = "Homocidios", main = "Homocidios por 1000 habitantes")

#MODELACIÓN CAUCHY Z:= |X|
#dcauajust <- function(x, b) {(2*b)/(pi*(x^{2} + b^2))}
#pcauajust <- function(q, b) {(2/pi)*atan(q/b)}
#qcauajust <- function(p, b) {b*tan(pi/2 * p)}
#rcauajust <- function(n, b) {
#           u <- runif(n)
#           rcauajust <- mapply(qcauajust, u, rep(b, n))
#           return(rcauajust)}
#(param_reg_no_0 <- fitdist(reg_NO_0, "cauajust", start=list( b=0.005))$estimate)
#sims.NO_0 <- rcauajust(n, param_reg_no_0)
#hist(sims.NO_0)
 
#MODELACIÓN PARETO AJUSTADA Z := X - a
#dpareajust <- function(x, a, b) {(b*a^{b})/(x + a)^{b + 1}}
#ppareajust <- function(q, a, b) {1 - (a/(q + a))^{b}}
#qpareajust <- function(p, a, b) {a/(1 - p)^{1/b} - a}
#rpareajust <- function(n, a, b) {
#            u <- runif(n)
#            rpareajust <- mapply(qpareajust, u, rep(a, n), rep(b, n))
#            return(rpareajust)}
#(param_reg_no_0 <- fitdist(reg_NO_0, "pareajust", start = list(a = 0.001, b = 0.005))$estimate)
#sims.NO_0 <- rpareajust(n, param_reg_no_0[1], param_reg_no_0[2])
#hist(sims.NO_0, prob = TRUE)
#summary(sims.NO_0); summary(reg_NO_0)
#ad.test(reg_NO_0, ppareajust, param_reg_no_0[1], param_reg_no_0[2])

#GAMMA
#(param_reg_no_0 <- fitdist(reg_NO_0, "gamma")$estimate)
#sims.NO_0 <- rgamma(n, param_reg_no_0[1], param_reg_no_0[2])
#hist(sims.NO_0, prob = TRUE)
#summary(sims.NO_0); summary(reg_NO_0)
#ad.test(reg_NO_0, pgamma, param_reg_no_0[1], param_reg_no_0[2])

#LOGIT-NORM
#(param_reg_no_0 <- fitdist(reg_NO_0, "logitnorm", start = list(mu = 0.001, sigma = 0.001))$estimate)
#sims.NO_0 <- rgamma(n, param_reg_no_0[1], param_reg_no_0[2])
#hist(sims.NO_0, prob = TRUE)
#summary(sims.NO_0); summary(reg_NO_0)
#ad.test(reg_NO_0, pgamma, param_reg_no_0[1], param_reg_no_0[2])


#MODELACIÓN CON BERNSTEIN
sims.homP1000_NO_0 <- sapply(runif(n), Fn.inv.Bernstein, reg_NO_0)
dev.new(); par(mfrow = c(1,2))
hist(reg_NO_0, prob = T, col = "blue", main = "Número de homicidios distinto de 0",
     xlim = c(0,8), xlab = "Homicidios")
hist(sims.homP1000_NO_0, prob = TRUE, col = "blue", 
     main = "Simulación vía Bernstein", xlim = c(0,8), xlab = "Homicidios")
summary(reg_NO_0)
summary(sims.homP1000_NO_0)

#Con lo anterior se tiene que F.C(x) = Fn.Bernstein(x, reg_NO_0) 
F.C <- function(x) Fn.Bernstein(x, reg_NO_0)
x1 <- seq(-2, 7, length.out = 1000)
x2 <- seq(min(homP1000), max(homP1000), length.out = 1000)
y2 <- sapply(x2, Fn.Bernstein, homP1000)
dev.new(); par(mfrow = c(1,2))
plot(x1, sapply(x1, F.D), lwd = 2, main = "F.D(x)", xlab = "x", ylab = "F.D(x)",
     ylim = c(0,1))
abline(v = 0, lty = 5, col = "red", lwd = 3)
plot(x2, y2, type = "l", lwd = 4, main = "F.C(x)", xlab = "x", ylab = "F.C(x)") 

#DEFINIENDO alfa = P(HomP1000 = 0)
F_HomP1000 <- function(x){
                       if(x <= 0) {res <- masa_0*F.D(x)}
                       else if(x > max(homP1000)){res <- 1} 
                       else {res <- masa_0*F.D(x) + (1 - masa_0)*F.C(x)}
   return(res)
} 
x3 <- seq(-1, 0, length.out = 200)
x3 <- c(x3, seq(min(homP1000), 10, length.out = 800))
system.time(y3 <- sapply(x3, F_HomP1000))
dev.new()
plot(x3, y3, main = "Función de distrubución de homP1000", xlab = "homP1000",
     ylab = "F", ylim = c(0,1))
points(0, masa_0, pch = 19, col = "red")
x4 <- seq(min(reg_NO_0), max(reg_NO_0), length.out = 10000)
y4 <- sapply(x4, function(x) Fn.Bernstein(x, homP1000))
dev.new()
plot(x4, y4, main = "F_HomP1000", xlab = "x",
     ylab = "F_{HomP1000}(x)", ylim = c(0,1), type = "l", lwd = 3, xlim = c(-1,max(homP1000)))
points(0, masa_0, pch = 19, col = "red")
lines(seq(-1, 0, length.out = 1000), rep(0, 1000), lwd = 3)


#================================================================================
#================================================================================
#======================= ANÁLISIS DEL MODELO A PARES ============================
#================================================================================
#================================================================================
library(subcopem2D)
# VECTOR ALEATORIO 1 = (IDH, Pnat)
V1 <- cbind(idh, Pnat)
head(V1)
dev.new(); par(mfrow = c(1,2))
plot(V1, main = "Observaciones de V1", ylab ="% Nativos", ylim = c(0,1))
pseudo_V1 <- 1/n * apply(V1, 2, rank)
head(pseudo_V1)
plot(pseudo_V1, main = "Pseudo-observaciones de V1", xlab = "u", ylab = "v",
     ylim = c(0,1))
V1_subcopem <- subcopemc(V1, m  = 50, display = TRUE)
V1_subcopem$depMon
V1_subcopem$depSup
rho <- cor(V1, method = "spearman")[1,2]
cor(V1, method = "kendall")
cor(V1, method = "pearson")


#VECTOR ALEATORIO 2 = (IDH, HomP1000)
V2 <- cbind(idh, homP1000)
head(V2)
str(V2)
dev.new(); par(mfrow = c(1,2))
plot(V2, main = "Observaciones de V2", xlab = "IDH", ylab = "HomP1000", 
     xlim = c(0,1), ylim = c(0,1))
pseudo_V2 <- 1/length(V2[,2]) * apply(V2, 2, rank)
plot(pseudo_V2)

V2_subcopem <- subcopem(V2, display = TRUE)
V2_subcopem$depMon
V2_subcopem$depSup


#DIST IDH|HomP1000 = 0
V2_0 <- V2[which(V2[ , 2] == 0), ]
dev.new()
plot(V2_0, ylim = c(0,1), xlim = c(0,1), pch = 20)
summary(V2_0[, 1])
dev.new()
hist(V2_0[, 1], prob = T, xlim = c(0,1), main = "IDH con HomP1000 = 0",
     xlab = "IDH", col = "yellow")

V2_No_0 <- V2[which(V2[ , 2] != 0), ]
V2_N_sub <- subcopem(V2_No_0, display = TRUE)
summary(V2_No_0[ , 1])
dev.new()
hist(V2_No_0[, 1], prob = T, xlim = c(0,1), main = "IDH con HomP1000 > 0",
     col = "yellow", xlab = "IDH")
dev.new(); par(mfrow = c(1,2))
plot(V2_0, ylim = c(0,1), xlim = c(min(idh),1), pch = 20)
plot(V2_No_0, xlim = c(0,1))
pseudo_V2 <- 1/length(V2_No_0[,2]) * apply(V2_No_0, 2, rank)
dev.new()
plot(pseudo_V2)

#VECTOR ALEATORIO 3 = (Pnat, HomP1000)
V3 <- cbind(Pnat, homP1000)
head(V3)
dev.new()
plot(V3, xlim = c(0,1), xlab = "% Mig", ylab = "HomP1000", 
     main = "Observaciones de V3")

V3_subcopem <- subcopem(V3, display = TRUE)
V3_0 <- V3[which(V3[ , 2] == 0), ]
head(V3_0[,1])
dev.new(); par(mfrow = c(1,2))
hist(V3_0[, 1], prob = TRUE, xlim = c(0,1), main = "% Nativos con HomP1000 = 0",
     xlab = "% Nativos", col = "green")
summary(V3_0[ , 1])

V3_No_0 <- V3[which(V3[ , 2] != 0), ]
head(V3_No_0)
dev.new()
hist(V3_No_0[, 1], prob = TRUE, xlim = c(0,1), main = "% Nativos con HomP1000 > 0",
     xlab = "% Nativos", col = "green")
summary(V3_No_0[, 1])
plot(V3_No_0)
pseudo.V3_no_0 <- 1/nrow(V3_No_0) * apply(V3_No_0, 2, rank)
dev.new()
plot(pseudo.V3_no_0)
dev.new()
hist(V3_No_0[ , 1], prob = TRUE, xlim = c(0,1))
summary(V3_No_0[ , 1])

V3_no_0_subcopem <- subcopem(V3_No_0, display = TRUE)

############################################################################
library(tidyverse)
datos_a <- as.data.frame(datos)
datos_a <- datos_a %>% 
           mutate(homP1000 = if_else(HOMP1000 == 0, "0", "Mayor a 0"))

ggplot(data = datos_a) +
geom_point(aes(POR_NAT, IDH, colour = homP1000)) +
xlab("% Nativos") +
ylab("IDH")
