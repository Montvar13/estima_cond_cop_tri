#==========================================================================================
#================= APLICACIÓN AL MODELO TEÓRICO, CAPÍTULO 3.2 =============================
#==========================================================================================

# Se simularán muestras del modelo para valores de teta = {0, 0.5, 1, 1.5},
# con dichas simulaciones se estimará la cópula de manera condicional y se 
# va a comparar con la cópula teórica. Las muestras será tamaño 100,000 y se
# irá decidiendo con base en los tiempos de ejecucuón de los códigos el número
# de evaluaciones. Después de obtener la matriz de evaluaciones se obtendrá el
# promedio de las diferencias en valor absoluto.

#                         P R I M E R A   P A R T E 

# A continuación se definirán las funciones teóricas del modelo y otras 
# que son necesarias para realizar el trabajo.

# 1.- FUNCIÓN SINILAR A OUTER EN 3D.

outer3d <- function(x, y, z, f){ 

    mat <- matrix( , 0, length(y))

    for(i.z in 1:length(z)){
                e.z <- rep(z[i.z], length(x)*length(x))
                e.x <- rep(x, length(z))
                e.y <- rep(y, each = length(x))

                mat.aux <- matrix(mapply(f, e.x, e.y, e.z), length(x), length(y))
                mat.0 <- rep(0 , length(y))
                mat.aux <- rbind(mat.0, mat.aux)

                mat <- rbind(mat, mat.aux)
    }
    mat.0 <- rep(0, nrow(mat))
    mat <- cbind(mat.0, mat)
    mat.0 <- matrix(0, length(x) + 1, length(y) + 1)
    mat <- rbind(mat.0, mat)
    mat <- unname(mat)
    return(mat)
}

# 2.- FUNCIONES DEL MODELO TEÓRICO:

# X
Fx <- function(x) {punif(x)}
fx <- function(x) { dunif(x) }

# Y
Fy <- function(y) {(y*(1-log(y)))*((0 < y)&(y < 1)) + 1*(y>=1)}
fy <- function(y) { (-log(y))*((0 < y)&(y < 1))}

# Z
# teta = 0
Fz_0 <- function(z) {(z*(1-log(z)))*((0 < z)&(z < 1)) + 1*(z>=1)}
# teta = 1 
Fz_1 <- function(z) {(z*log(z)+4*sqrt(z)-3*z)*((0 < z)&(z < 1)) + 1*(z>=1)}
# teta > 0 & teta != 1
Fz_t <- function(z,t) {((1/(1-t))*((((t+1)*(1-t^2))/(t))*z^(1/(t+1)) - z/t + 
                       t^2*z^(1/t)))*((0 <= z)&(z < 1)) + 1*(z>=1)}

# Z|X
# teta = 0
Fz1x_0 <- function(z,x) {punif(z, 0, x)}
# teta = 1
Fz1x_1 <- function(z,x) {((2*z*log(x)-z*log(z)+z)/(x^2))*((0<z)&(z<x^2)) + 1*(x^2 <= z)}
# teta > 0 & teta !=1
Fz1x_t <- function(z,x,t) {((1/(1-t))*(z/x^(t+1) - 
                           (t*z^(1/t))/x^(1/t + 1)))*((0 <= z)&(z < x^(t+1))) + 
                           1*(x^(t+1)<=z)} 
# Y|X
Fy1x <- function(y,x){(y/x)*((0 <= y)&(y < x)) + 1*(y>=x)}
#Z|Y
Fz1y_t <- function(z,y,t) {((-z/log(y))*(1/y^(t+1) - 1/y^t))*((0 <= z)&(z <= y^(t+1))) + 
                           (z/(y^t*log(y)) - 1/log(y) - log(z)/log(y) + 
                           (t+1))*((y^(t+1) < z)&(z < y^t)) + 1*(z>=y^t)}

# INVERSAS
# X
invFx <- function(u){qunif(u,0,1)}
# Y
invFy <- function(v) {uniroot(function(y) Fy(y) - v, interval = c(0.00000000001, 1))$root}
# Z
invFz_0 <- function(w) {uniroot(function(z) Fz_0(z) - w, 
                        interval = c(0.000000000001, 1))$root}

invFz_1 <- function(w) {
                         if(w == 0){
                              0
                         }else if(w == 1){ 
                              1
                         } else{
                          uniroot(function(z) Fz_1(z) - w, 
                          interval = c(0.000000000001, 0.999999999))$root
                         }
}

invFz_t <- function(w,t){ uniroot(function(z) Fz_t(z,t) - w, 
                          interval = c(0.00000000000001, 1))$root}
# Y|X
invFy1x <- function(v,x)qunif(v,0,x)
# Z|XY
# teta = 0
invFz1xy_0 <- function(w,x,y) qunif(w, 0, x)
# teta = 1
invFz1xy_1 <- function(w,x,y) qunif(w, 0, x*y)
# teta > 0 & teta != 1
invFz1xy_t <- function(w,x,y,t) qunif(w, 0, x*y^t)
# Z|X
invFz1x_0 <- function(w,x){qunif(w, 0, x)}
invFz1x_1 <- function(w,x){uniroot(function(z) Fz1x_1(z,x) - w, 
                           interval = c(0.0000000000001,0.999))$root}
invFz1x_t <- function(w,x,t){uniroot(function(z) Fz1x_t(z,x,t) - w, 
                             interval = c(0.00000000001,0.999))$root}
#Z|Y
invFz1y_t <- function(w,y,t){
                ((log(y)*y^(t + 1)*w)/(y - 1))*((0 <= w) & 
                (w <= (y - 1)/log(y))) + (uniroot(function(z) Fz1y_t(z,y,t) - w, 
                interval = c(0.0000000001,0.9999999999))$root)*(((y - 1)/log(y) < w) & 
                (w < (y - 1)/(y*log(y))))}



## CÓPULA TRIVARIADA TEÓRICA

# Para teta = 0
C_XYZ_0 <- function(u,v,w){
           u*((u <= invFy(v))&(u <= invFz_0(w))*(u < 1)) +
           (invFz_0(w)*(1 + log(u) - log(invFz_0(w))))*((invFz_0(w) < u)&(u <= invFy(v))) +
           (invFy(v)*(1 + log(u) - log(invFy(v))))*((invFy(v) < u)&(u <= invFz_0(w))) +
           (invFy(v)*(log(invFz_0(w)) - log(invFy(v)) - invFz_0(w)/u + 2))*((invFy(v) < 
           invFz_0(w))&(invFz_0(w) < u)) +
           (invFz_0(w)*(log(invFy(v)) - log(invFz_0(w)) - invFy(v)/u + 2))*((invFz_0(w) < 
           invFy(v))&(invFy(v) < u)) + 
           1*(u>=1 & invFy(v) >= 1 & invFz_0(w) >= 1)
           }

densC_XYZ_0 <- function(u,v,w){
               (1/(u^{2}*log(invFy(v))*log(invFz_0(w))))*((invFy(v) <= u)&(
               invFz_0(w) <= u))}

# Para teta = 1
C_XYZ_1 <- function(u,v,w){
           u*((u < invFy(v))&(u^{2} <= invFz_1(w))) +
           (4*sqrt(invFz_1(w)) + (invFz_1(w)*(log(invFz_1(w)) - 2*log(u) - 3))/u)*((u 
            <= invFy(v))*((0 < invFz_1(w))&(invFz_1(w) < u^{2})))+
           (invFy(v)*(1 + log(u) - log(invFy(v))))*((invFy(v) < u)*(u*invFy(v) <= 
           invFz_1(w))) +
           (log(invFz_1(w))*(invFy(v) + invFz_1(w)/u) - log(invFy(v))*(2*invFy(v) + 
           invFz_1(w)/u) - invFz_1(w)/u*(log(u) + 2) + 3*invFy(v))*((invFy(v) < u)&(
           (invFy(v)^2 < invFz_1(w))&(invFz_1(w) < u*invFy(v)))) +
           (4*sqrt(invFz_1(w)) - (invFz_1(w)/u)*(log(u) + log(invFy(v)) - log(invFz_1(w)) +
           2) - invFz_1(w)/invFy(v))*((invFy(v) < u)&((0 < invFz_1(w))&(invFz_1(w) <= 
           (invFy(v))^{2}))) +
           1*(u>=1 & invFy(v)>=1 & invFz_1(w) >= 1)
           }

densC_XYZ_1 <- function(u,v,w) {
               (-1/(u^{2}*log(invFy(v))*invFy(v)*(log(invFz_1(w)) + 2/sqrt(invFz_1(w)) - 
               2)))*((invFy(v) < u)&(invFz_1(w) < u*invFy(v)))}

# Para teta > 0 & teta != 1

C_XYZ_t <- function(u,v,w,t){
           u*((invFy(v) > u)&(invFz_t(w,t) >= u^{t+1})) +
           (invFy(v)*(1 + log(u) - log(invFy(v))))*((invFy(v) < u)&(invFz_t(w,t) >= 
           u*(invFy(v))^{t})) +
           ((invFz_t(w,t))^{1/(t+1)} + (1/(1-t))*(invFz_t(w,t)/t*((invFz_t(w,t))^{-t/(t+1)}
           - u^{-t}) - t^2*(invFz_t(w,t))^{1/t}*((invFz_t(w,t))^{-1/(t^2+t)} - 
           u^{-1/t}) ))*((invFy(v) >= u)&(u^{t+1} > invFz_t(w,t))) +
           (invFy(v)*(log(invFz_t(w,t)) - (t+1)*log(invFy(v)) + 1) + (1/(1-t))*(invFy(v) - 
           invFz_t(w,t)*(invFy(v))^{1 - t}/u - t^{2}*(invFy(v) - 
           (invFz_t(w,t)/u)^{1/t})))*((invFy(v) < u)&(((invFy(v))^{t+1} < invFz_t(w,t))&(
           invFz_t(w,t) < u*(invFy(v))^{t}))) +
           ((invFz_t(w,t))^{1/(t+1)} + ((invFz_t(w,t))^{1/(t+1)} - 
           invFz_t(w,t)*invFy(v)^{-t})/t + t*((invFz_t(w,t))^{1/(t+1)} - 
           (invFz_t(w,t)/u)^{1/t}) + (1/(1-t))*((invFz_t(w,t))^{1/(t+1)} - 
           t*(invFz_t(w,t))^{1/t}*((invFz_t(w,t))^{-1/(t^2+t)} - u^{-1/t}) - 
           invFz_t(w,t)*(invFy(v))^{1-t}/u))*((invFy(v) < u)&((0 < invFz_t(w,t))&(
           invFz_t(w,t) < (invFy(v))^{t+1}))) +
           1*(u>=1 & v >=1 & w>=1)}

densC_XYZ_t <- function(u,v,w,t) {
               (-1/((u^2)*((invFy(v))^{t})*(log(invFy(v)))*((invFz_t(w,t))^{-t/(t+1)} - 1 +
               t^{2}*(invFz_t(w,t))^{1/t - 1}*(1 - (invFz_t(w,t))^{-1/(t^2+t)}))))*((
               invFy(v) < u)&(invFz_t(w,t) < u*(invFy(v))^{t}))}

# FUNCIÓN PARA ESTIMAR LA CÓPULA DE MANERA CONDICIONAL
# Debido a que las variables para este modelo son todas continuas, se usará una
# función específica para este tipo de variables.

subcop_condc <- function(mat.XYZ, m = nrow(mat.XYZ)){

     n <- nrow(mat.XYZ)
     mat.UVW <- apply(mat.XYZ, 2, rank)

     kp <- round(seq(1, n, length.out = m), 0)
 
     mat.UVW <- mat.UVW[order(mat.UVW[ , 3]), ]
     mat.UVW <- as.data.frame(mat.UVW)
     colnames(mat.UVW) <- c("U", "V", "W")

     suma.aux <- function(kpU, kpV){sum(mat.UV_W$U <= kpU & mat.UV_W$V <= kpV)}
     suma <- function(kpU, kpV){mapply(suma.aux, kpU, kpV)}

     subcopula <- matrix( , m*(m+1), m+1)
     subcopula[ , 1] <- 0

     for(k in 1:length(kp)){
                   mat.UV_W <- mat.UVW[1:kp[k], 1:2]
                   n_k <- nrow(mat.UV_W)
                   subcopula[((k-1)*m + k), ] <- 0
                   subcopula[((k-1)*m + k + 1):(k*m + k), 2:(m+1)] <- 1/n_k * outer(kp, kp, suma) * kp[k]/n
     }

     sub_0 <- matrix(0, m+1, m+1)

     subcopula <- rbind(sub_0, subcopula)
     subcopula <- unname(subcopula)
     return(subcopula)
}

#                        S E G U N D A    P A R T E

# En esta sección se simularán las obervaciones del modelo teórico para cada uno
# de los valores que se van a estudiar, {0, 0.5, 1}, posteriormente se van a
# obtener las evaluaciones de la cópula teórica para cada modelo, obteniendo una
# matríz de k^2*k, para compararla con la cópula resultante de la estimación 
# condicional. Se calculará la media de las distancias absulutas. Este proceso se
# va a repetir 100 veces para cada modelo.

# TETA = 0

 n <- 27000
 uk <- 50
 kp <- round(seq(1, n, length.out = uk), 0)/n
 system.time(cop_teo_0 <- outer3d(kp, kp, kp, C_XYZ_0))
 dif_0 <- numeric()

for( i in 1:100){

 u <- runif(n,0,1)
 x_0 <- invFx(u)
 v <- runif(n,0,1)
 y_0 <- mapply(invFy1x, v, x_0)
 w <- runif(n,0,1)
 z_0 <- mapply(invFz1xy_0, w, x_0, y_0)
 ma.obs_0 <- cbind(x_0,y_0,z_0)
# head(ma.obs_0)
# summary(ma.obs_0)

subcond_0 <- subcop_condc(ma.obs_0, uk)

dif_0 <- c(dif_0, sum(abs(cop_teo_0 - subcond_0))/((uk+1)^3))
}

hist(dif_0, col = "yellow", main = "Dif Cop Cond VS Cop Teo", freq = FALSE)
summary(dif_0)

# TETA = 0.5

 n <- 27000
 uk <- 50
 kp <- round(seq(1, n, length.out = uk), 0)/n
 system.time(cop_teo_05 <- outer3d(kp, kp, kp, function(x,y,z) C_XYZ_t(x,y,z,0.5)))
 dif_05 <- numeric()

for( i in 1:100){

t <- 0.5
u <- runif(n,0,1)
x_05 <- invFx(u)
v <- runif(n,0,1)
y_05 <- mapply(invFy1x, v, x_05)
w <- runif(n,0,1)
z_05 <- mapply(invFz1xy_t, w, x_05, y_05, rep(t, n))
ma.obs_05 <- cbind(x_05,y_05,z_05)
#head(ma.obs_05)
#summary(ma.obs_05)
#rank.UVW_05 <- apply(ma.obs_05, 2, rank)

subcond_05 <- subcop_condc(ma.obs_05, uk)

dif_05 <- c(dif_05, sum(abs(cop_teo_05 - subcond_05))/((uk+1)^3))
}

hist(dif_05, col = "yellow", main = "Dif Cop Cond VS Cop Teo", freq = FALSE)
summary(dif_05)


# TETA = 1

 n <- 27000
 uk <- 50
 kp <- round(seq(1, n, length.out = uk), 0)/n
 system.time(cop_teo_1 <- outer3d(kp, kp, kp, C_XYZ_1))
 dif_1 <- numeric()

for( i in 1:100){

u <- runif(n,0,1)
x_1 <- invFx(u)
v <- runif(n,0,1)
y_1 <- mapply(invFy1x, v, x_1)
w <- runif(n,0,1)
z_1 <- mapply(invFz1xy_1, w, x_1, y_1)
ma.obs_1 <- cbind(x_1,y_1,z_1)

subcond_1 <- subcop_condc(ma.obs_1, uk)

dif_1 <- c(dif_1, sum(abs(cop_teo_1 - subcond_1))/((uk+1)^3))
}

hist(dif_1, col = "yellow", main = "Dif Cop Cond VS Cop Teo", freq = FALSE)
summary(dif_1)


dif <- cbind(dif_0, dif_05, dif_1)

#write.csv(dif, "C:/Users/DELL/Dropbox/Tesis/AVANCES/Ejercicio2/tesis_evs.csv")




# VINE CÓPULAS

# De acuerdo al análiis que se realizó, se determinó que se trabajaría con 
# teta = 1 y las cópulas C_YZ|X y CXY|Z.

#Descomposiciones resultantes y sugeridas por vine copulas:
# 1.- c_XYZ = c_YZ|X*c_XY*c_XZ
# 2.- c_XYZ = c_XY|Z*c_XZ*c_YZ

# MEDANTE EL PROCESO DE SIMULACIÓN DERIVADO DE VINE CÓPULAS
#PARA (1) OBTENEMOS LAS SIGUIENTES FUNCIONES:

hx1y <- function(u,v){
         (1 - (log(u)/log(invFy(v))))*((0 < invFy(v))&(invFy(v) < u)&(u < 1)) + 
         1*(u == 1)}

hz1y <- function(w,v){
         1*((0 < Fz_1(invFy(v)))&(Fz_1(invFy(v)) <= w)&(w < 1)) +
         ((invFz_1(w) - invFy(v)*(log(invFz_1(w)) + 1))/(invFy(v)*log(invFy(v))) + 
         2)*((0 < Fz_1((invFy(v))^2))&(Fz_1((invFy(v))^2) <= w)&(w < Fz_1(invFy(v)))) +
         ((invFz_1(w)*(invFy(v) - 1))/(log(invFy(v))*(invFy(v))^2))*((0 < w)&(w < 
         Fz_1((invFy(v))^2))&(Fz_1((invFy(v))^2) < 1))}

hz1xy <- function(w, u, y){
         ((w*y^u*log(y))/(y-1))*(((0 < u)&(u < 1))&((0 < w)&(w < (y-1)/(log(y)))))+
         (invFz1y_t(w,y,1)*y^(u-2))*(((0 < u)&(u < 1))&(((y-1)/(log(y)) < w)&(w < 
         (y^(1-u)-1)/(log(y))+u))) +
         1*((((0 < u)&(u < 1))&((y^(1-u)-1)/(log(y)) + u < w)) | (w == 1))}

# Y NECESITAMOS LAS SIGUIENTES INVERSAS

invhx1y <- function(g,v){
           ((invFy(v))^(1-g))*((0 < g)&(g < 1)) + 1*(g == 1)
}

invhz1x_y <- function(a,u,v){
             ((a*(v-1))/(log(v)*v^u))*((0 < a)&(a <= v^u)) + 
             ((a*v^(1 - u) - log(a*v^(2-u)) - 1)/log(v) + 2)*((v^u < a)&(a <= 1))}

hz1y_aux <- function(w,v){
        ((invFz_1(w) - invFy(v)*(log(invFz_1(w)) + 1))/(invFy(v)*log(invFy(v))) + 2)*((0 < 
        Fz_1((invFy(v))^2))&(Fz_1((invFy(v))^2) <= w)&(w <= Fz_1(invFy(v))))}

invhz1y <- function(b,v){
           if((0 < b)&(b < ((invFy(v) - 1)/log(invFy(v))))){
                (Fz_1((b*log(invFy(v))*invFy(v)^2)/(invFy(v) - 1)))
           } else if( ((invFy(v) - 1)/log(invFy(v)) <= b) & (b < 1) ){
                (uniroot(function(w) hz1y_aux(w, v) - b, interval = c(Fz_1(invFy(v)^2), 
                Fz_1(invFy(v))))$root) 
           } else if(b == 1){
                 1 
           } else { 
                 0
           }
}

## PROCESO DE SIMULACIÓN
# Paso 1

m <- 27000
p <- round(runif(m), 4)
r <- round(runif(m), 4)
s <- round(runif(m), 4)

V <- p

#Paso 2

U <- numeric()
for(i in 1:m){
  U[i] <- invhx1y(r[i], V[i])
}

#Paso 3

W <- numeric()
for(i in 1:m){
  W[i] <- invhz1y(invhz1x_y(s[i], hx1y(U[i], V[i]), V[i]), V[i])
  print(i)
}

cop_obs <- cbind(U, V, W)
write.xlsx(cop_obs, "C:/Users/DELL/Dropbox/Tesis/AVANCES/Ejercicio2/tesis_vine_cop_obs.csv")

cop_obs_resp <- cop_obs
cop_obs <- cop_obs[!is.na(cop_obs[ , "W"]), ]
cop_obs <- cop_obs[1:27000 , ]

 m <- 27000
 uk <- 50
 kp <- round(seq(1, m, length.out = uk), 0)/m

 cop_emp <- function(mat, kp1, kp2, kp3){
           (1/nrow(mat))* nrow(mat[mat[ , "U"] <= kp1 & 
                                   mat[ , "V"] <= kp2 & 
                                   mat[ , "W"] <= kp3, ])
  }

emp <- outer3d(kp, kp, kp, function(x,y,z) cop_emp(cop_obs, x,y,z))
teo <- outer3d(kp, kp, kp, C_XYZ_1)
dif <- (1/50^3)*sum(abs(emp - teo))

#POR OTRA PARTE, PARA (2) OBTENEMOS LAS SIGUIENTES FUNCIONES:

hz1x <- function(w,u){
        ((invFz_1(w)*(1 + 2*log(u) - 
        log(invFz_1(w))))/(u^2))*((0 < w)&(w < Fz_1(u^2))&(Fz_1(u^2) < 1)) + 
        1*((0 < Fz_1(u^2))&(Fz_1(u^2) <= w)&(w < 1))}

hz1x_aux <- function(w,u){
        ((invFz_1(w)*(1 + 2*log(u) - log(invFz_1(w))))/(u^2))*((0 < w)&(w <= 
         Fz_1(u^2))&(Fz_1(u^2) < 1))}

invhz1x <- function(a,u){
           if((0 < a) & (a < 1)){
               (uniroot(function(w) hz1x_aux(w, u) - a, 
               interval = c(0.00001, Fz_1(u^2)))$root) 
           } else if(a == 1){
             1
           } else {0}}

hy1x <- function(v,u){
        (invFy(v)/u)*((0 < v)&(v < invFy(u))) + 
        1*(((0 < Fy(u))&(Fy(u) <= v) & (v < 1)) | (v == 1))
}

invhy1x <- function(a,u){
           Fy(a*u)*((0 < a) & (a < 1)) + 
           1*(a == 1)
}

hy1zx <- function(v,w,u){
         ((log(v) + 2*log(u) - log(invFz1x_1(w,u)))/(2*log(u) - log(invFz1x_1(w,u))))*((0 <
         invFy(w))&(invFy(w) < v)&(v < 1)) + 
         1*(((0 < invFy(w))&(invFy(w) < v)&(v < 1)) & v == 1)}

invhy1zx <- function(b,w,u){
            (exp((b-1)*(2*log(u) - log(invFz1x_1(w,u)))))*((hy1zx(invFy(w),w,u) < b) & (
             b < 1)) +
             1*(b == 1)
}


## PROCESO DE SIMULACIÓN
# Paso 1

m <- 27000
p <- round(runif(m), 4)
r <- round(runif(m), 4)
s <- round(runif(m), 4)

U <- p

#PASO 2 

W <- numeric()

for(i in 1:m){
   W[i] <- invhz1x(r[i], U[i])
   print(i)
}

pos_na <- which(is.na(W))
U <- U[-pos_na]
W <- W[-pos_na]


#Paso 3

V <- numeric()

for(i in 1:m){
  V[i] <- invhy1x(invhy1zx(s[i], hz1x(W[i], U[i]), W[i]), U[i])
  print(i)
}

pos_na <- which(is.na(V))
U <- U[-pos_na]
W <- W[-pos_na]
V <- V[-pos_na]


cop_obs2 <- cbind(U, V, W)
cop_obs <- rbind(cop_obs, cop_obs2)
cop_obs_resp <- cop_obs
cop_obs <- cop_obs[1:27000, ]


# write.csv(cop_obs, "C:/Users/DELL/Dropbox/Tesis/AVANCES/Ejercicio2/tesis_vine_cop_obs2.csv")

m <- 27000
uk <- 50
kp <- round(seq(1, m, length.out = uk), 0)/m

cop_emp <- function(mat, kp1, kp2, kp3){
           (1/nrow(mat))* nrow(mat[mat[ , "U"] <= kp1 & 
                                   mat[ , "V"] <= kp2 & 
                                   mat[ , "W"] <= kp3, ])
}

emp <- outer3d(kp, kp, kp, function(x,y,z) cop_emp(cop_obs, x,y,z))
teo <- outer3d(kp, kp, kp, C_XYZ_1)
dif <- (1/50^3)*sum(abs(emp - teo))

