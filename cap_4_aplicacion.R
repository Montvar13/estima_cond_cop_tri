#m y DATOS
#CONTINUA

subcop_condc <- function(mat.XYZ, m = nrow(mat.XYZ)){

     n <- nrow(mat.XYZ)

     kpU <- round(seq(1, n, length.out = m), 0)
     kpV <- kpU
     kpW <- kpU
 
     mat.XYZ <- mat.XYZ[order(mat.XYZ[ , 3]), ]
     mat.XYZ <- as.data.frame(mat.XYZ)
     colnames(mat.XYZ) <- c("X", "Y", "Z")

     suma.aux <- function(kpU, kpV){sum(mat.XY_Z$X <= kpU & mat.XY_Z$Y <= kpV)}
     suma <- function(kpU, kpV){mapply(suma.aux, kpU, kpV)}

     subcopula <- matrix( , m*(m+1), m+1)
     subcopula[ , 1] <- 0

     for(k in 1:length(kpW)){
                   mat.XY_Z <- mat.XYZ[1:kpW[k], 1:2]
                   n_k <- nrow(mat.XY_Z)
                   subcopula[((k-1)*m + k), ] <- 0
                   subcopula[((k-1)*m + k + 1):(k*m + k), 2:(m+1)] <- 1/n_k * outer(kpU, kpV, suma) * kpW[k]/n
     }

     sub_0 <- matrix( , m+1, m+1)
     sub_0 <- 0

     subcopula <- rbind(sub_0, subcopula)
     return(subcopula)
}


#NO CONTINUAS

 datos <- read.csv("C:/Users/DELL/Dropbox/Tesis/NUEVOS DATOS/TRES_VARIABLES.csv", sep = ",")
 datos$IDH <- jitter(datos$IDH)
 mat.XYZ <- datos[ , c(6, 4, 5)]
 mat.XYZ <- as.data.frame(mat.XYZ)
 colnames(mat.XYZ) <- c("idh", "Pnat", "homP1000")
 n <- nrow(mat.XYZ)

 suma.aux <- function(kpU, kpV){sum(mat.XY_Z$idh <= kpU & mat.XY_Z$Pnat <= kpV)}
 suma <- function(kpU, kpV){mapply(suma.aux, kpU, kpV)}

 M <- 50

 X.u <- sort(unique(mat.XYZ$idh))
 m1 <- length(X.u)
 ind_x <- round(seq(1, m1, length.out = M), 0)
 X.u <- X.u[ind_x]

 Y.u <- sort(unique(mat.XYZ$Pnat))
 m2 <- length(Y.u)
 ind_y <- round(seq(1, m2, length.out = M), 0)
 Y.u <- Y.u[ind_y]

 Z.u <- sort(unique(mat.XYZ$homP1000))
 m3 <- length(Z.u)
 ind_z <- round(seq(1, m3, length.out = M), 0)
 Z.u <- Z.u[ind_z]

 ranZ <- cumsum(as.vector(table(mat.XYZ$homP1000))/n)
 ranZ.u <- ranZ[ind_z]

 mat.XYZ <- mat.XYZ[order(mat.XYZ$homP1000), ]

 subcopula <- matrix( , M*(M+1), M + 1)
 subcopula[ , 1] <- 0

 system.time(
   for(k in 1:length(Z.u)){

   mat.XY_Z <- mat.XYZ[mat.XYZ[ , 3] <= Z.u[k], 1:2]
   n_k <- nrow(mat.XY_Z)

   subcopula[((k-1)*M + k), ] <- 0
   subcopula[((k-1)*M + k + 1):(k*M + k), 2:(M+1)] <- 1/n_k * outer(X.u, Y.u, suma) * ranZ.u[k]

   }
 )
 sub_0 <- matrix( 0, M+1, M+1)
 subcopula <- rbind(sub_0, subcopula)

 write.xlsx(subcopula, "C:/Users/DELL/Dropbox/Tesis/NUEVOS DATOS/subcopula.xlsx")


#=================================================================================
#ANÁLISIS DE LA SUBCÓPULA DEL MODELO...

indU <- c(0, ind_x/m1)
indV <- c(0, ind_y/m2)
indW <- c(0, ranZ.u)

## CÓPULA PI PARA EL DOMINIO DE LA SUBCÓPULA

pi <- matrix( , 0, 51)
for(i in 1:length(indW)){
    pi <- rbind(pi, outer(indU, indV, function(x,y) x*y*indW[i]))
}

d <- max(subcopula - pi) - max(pi - subcopula) 


## CÓPULA W

# 1.- FUNCIÓN SINILAR A OUTER EN 3D.

cop_max_fn <- function(x, y, z){ 

    mat <- matrix( , 0, 3)

    e.z <- rep(rep(z, each = length(x)), length(y))
    e.x <- rep(x, length(y)*length(z))
    e.y <- rep(y, each = length(x)*length(y))
    mat <- rbind(mat, cbind(e.x, e.y, e.z))

    outer3d_aux <- matrix(apply(mat, 1, function(x) max(sum(x) - 2, 0)), length(x)* length(z), length(y))
    return(outer3d_aux)
}

cop_max <- cop_max_fn(indU, indV, indW)

d_max <- max(cop_max - pi) - max(pi - cop_max) 

dep_mon <- (-1)*d/d_max

sum(subcopula > pi)/132651



