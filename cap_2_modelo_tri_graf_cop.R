
  ###DISTRIBUCIONES MARGINALES X, Y, Z, X|Y, X|Z, Y|X, Y|Z, Z|X, Z|Y

#X
Fx <- function(x) {x*((0 <= x)&(x < 1)) + 1*(x>=1)}
fx <- function(x) { dunif(x) }

#Y
Fy <- function(y) {(y*(1-log(y)))*((0 < y)&(y < 1)) + 1*(y>=1)}
fy <- function(y) { (-log(y))*((0 < y)&(y < 1))}

#Z
#teta = 0
Fz_0 <- function(z) {(z*(1-log(z)))*((0 < z)&(z < 1)) + 1*(z>=1)}
#teta = 1 
Fz_1 <- function(z) {(z*log(z)+4*sqrt(z)-3*z)*((0 < z)&(z < 1)) + 1*(z>=1)}
#teta > 0 & teta != 1
Fz_t <- function(z,t) {((1/(1-t))*((((t+1)*(1-t^2))/(t))*z^(1/(t+1)) - z/t + t^2*z^(1/t)))*((0 <= z)&(z < 1)) + 1*(z>=1)}

#X|Y
Fx1y <- function(x,y){(1-log(x)/log(y))*((y <= x)&(x < 1)) + 1*(x >= 1)}

#X|Z
#teta = 0
Fx1z_0 <- function(x,z) {(1-log(x)/log(z))*((z <= x)&(x < 1)) + 1*(x >= 1)} 
#teta =1
Fx1z_1 <- function(x,z) {((1/(log(z)+2/sqrt(z)-2))*(2/sqrt(z) + (log(z)-2*log(x)-2)/x))*((sqrt(z) <= x)&(x < 1)) + 1*(x >= 1)}
#teta > 0 & teta != 1
Fx1z_t <- function(x,z,t){((z^(-t/(t+1)) - x^(-t) + t^2*z^((1-t)/t)*(x^(-1/t)-z^(-1/(t^2+t))))/(z^(-t/(t+1)) - 1 + t^2*z^((1-t)/t)*(1-z^(-1/(t^2+t)))))*((z^(1/(t+1)) <= x)&(x < 1)) + 1*(x >= 1)}

#Y|X
Fy1x <- function(y,x){(y/x)*((0 <= y)&(y < x)) + 1*(y>=x)}

#Y|Z
#teta = 0
Fy1z_0 <- function(y,z){((y/log(z))*(1-1/z))*((0 < y)&(y <= z)) + ((1/log(z))*(log(z)-log(y)+y-1))*((z<y)&(y<1)) + 1*(y>=1)}
#teta = 1
Fy1z_1 <- function(y,z) {(1/(log(z)+2/sqrt(z)-2))*((y/z - 1 + log(z) - log(y))*((z <= y)&(y <= sqrt(z))) + (2/sqrt(z) - 1/y - 1 + log(z)-log(y))*((sqrt(z) < y)&(y < 1))) + 1*(y>=1)}
#teta > 0 & teta !=1 
Fy1z_t <- function(y,z,t) {((t-t^2)/(z^(-t/(t+1))-1+t^2*z^(1/t -1)*(1-z^(-1/(t^2+t)))))*(((y - z^(1/t))/z + (z^((1-t)/t)-y^(1-t))/(1-t))*((z^(1/t) <= y)&(y <= z^(1/(t+1)))) + (z^(-t/(t+1))*(1+1/t) - (y^(-t))/t - (1/(1-t))*(y^(1-t)-t*z^(1/t - 1)))*((z^(1/(t+1)) < y)&(y < 1))) + 1*(y>=1)}

#Z|X
#teta = 0
Fz1x_0 <- function(z,x) {punif(z, 0, x)}
#teta = 1
Fz1x_1 <- function(z,x) {( (2*z*log(x)-z*log(z)+z)/(x^2))*((0<z)&(z<x^2)) + 1*(x^2 <= z)}
#teta > 0 & teta !=1
Fz1x_t <- function(z,x,t) {((1/(1-t))*(z/x^(t+1) - (t*z^(1/t))/x^(1/t +1)))*((0 <= z)&(z < x^(t+1))) + 1*(x^(t+1)<=z)} 

#Z|Y
Fz1y_t <- function(z,y,t) {((-z/log(y))*(1/y^(t+1) - 1/y^t))*((0 <= z)&(z <= y^(t+1))) + (z/(y^t*log(y)) - 1/log(y) - log(z)/log(y) +(t+1))*((y^(t+1) < z)&(z < y^t)) + 1*(z>=y^t)}



#INVERSAS
#X
invFx <- function(u){qunif(u,0,1)}
#Y
invFy <- function(v) {uniroot(function(y) Fy(y) - v, interval = c(0.0000001,0.999))$root}
#Z
invFz_0 <- function(w) {uniroot(function(z) Fz_0(z) - w, interval = c(0.0001,0.999))$root}
invFz_1 <- function(w) {uniroot(function(z) Fz_1(z) - w, interval = c(0.000000000001,0.999))$root}
invFz_t <- function(w,t){uniroot(function(z) Fz_t(z,t) - w, interval = c(0.0001,0.999))$root}
#X|Y
invFx1y <- function(u,y){(y^(1 - u))*((0 <= y) & (y < 1)) + 1*(u == 1)}
#X|Z
invFx1z_0 <- function(u,z){z^(1-u)}
invFx1z_1 <- function(u,z){uniroot(function(x) Fx1z_1(x,z) - u, interval = c(0.0001,0.999))$root}
invFx1z_t <- function(u,z,t){uniroot(function(x) Fx1z_t(x,z,t) - u, interval = c(0.0001,0.999), tol = .Machine$double.eps^0.25)$root}
#Y|X
invFy1x <- function(v,x){qunif(v, 0, x)}
#Y|Z
invFy1z_0 <- function(v,z){((z*v*log(z))/(z-1))*((0 < v)&(v <= (z-1)/log(z))) + (uniroot(function(y) (1/log(z))*(log(z)-log(y)+y-1) - v, interval = c(0.0001,0.999))$root)*(((z-1)/log(z)<v)&(v<1))}
invFy1z_1 <- function(v,z){uniroot(function(y) Fy1z_1(y,z) - v, interval = c(0.0001,0.999))$root}
invFy1z_t <- function(v,z,t){uniroot(function(y) Fy1z_t(y,z,t) - v, interval = c(0.0001,0.999), tol = .Machine$double.eps^0.25)$root}
#Z|X
invFz1x_0 <- function(w,x){qunif(w, 0, x)}
invFz1x_1 <- function(w,x){uniroot(function(z) Fz1x_1(z,x) - w, interval = c(0.0000000000001,0.999))$root}
invFz1x_t <- function(w,x,t){uniroot(function(z) Fz1x_t(z,x,t) - w, interval = c(0.00000000001,0.999))$root}
#Z|Y
invFz1y_t <- function(w,y,t){((log(y)*y^(t + 1)*w)/(y - 1))*((0 <= w) & (w <= (y - 1)/log(y))) + 
(uniroot(function(z) Fz1y_t(z,y,t) - w, interval = c(0.0000000001,0.9999999999))$root)*(((y - 1)/log(y) < w) & (w < (y - 1)/(y*log(y))))}


#(X,Y)

Fxy <- function(x,y){x*(((0<=x)&(x<1))&((x <= y))) +
                    (y*(1+log(x)-log(y)))*(((0<x)&(x<1))&((0<y)&(y<x))) +
                    (Fy(y))*((x>=1)&((0<y)&(y<1))) + 
                    (1)*(x>=1&y>=1)}
fxy <- function(x,y) {(1/x)*((0 < y)&(y < x)&(x < 1))} 


#(X,Z)

Fxz_0 <- function(x,z){x*(((0<=x)&(x<1))&((x <= z)&(z<1))) + (z*(1+log(x)-log(z)))*(((0<=x)&(x<1))&((0<z)&(z<x))) + (Fz_0(z))*((x>=1)&((0<z)&(z<1))) + (1)*(x>=1&z>=1)}

Fxz_1 <- function(x,z){x*(((0<=x)&(x<1))&((x^2<=z)&(z<1))) +
                      (4*sqrt(z)+(z*(log(z)-2*log(x))-3*z)/(x))*(((0<=x)&(x<1))&((0<z)&(z<x^2))) +
                      (z*log(z)+4*sqrt(z)-3*z)*((x>=1)&((0<z)&(z<x^2))) +
                      (1)*(x>=1&z>=1)}

Fxz_t <- function(x,z,t){x*(((0<=x)&(x<1))&((x^(t+1)<=z)&(z<1))) +
                        ((1/(1-t))*(((z/t)*(z^(-t/(t+1)) - x^(-t))) - ((t^2*z^(1/t))*(z^(-1/(t^2+t)) - x^(-1/t)))) + z^(1/(t+1)))*(((0<=x)&(x<1))&((0<z)&(z<x^(t+1)))) +
                        Fz_t(z,t)*((x>=1)&((0<z)&(z<1))) +
                        (1)*(x>=1&z>=1)}

#(Y,Z)


Fyz_0 <- function(y,z)  {(y*(log(z)-log(y)+2-z))*(((0<y)&(y<1))&((y<z)&(z<1))) +
                         (Fy(y))*(((0<y)&(y<1))&(z>=1)) +
                         (z*(log(y)-log(z)+2-y))*(((0<y)&(y<1))&((0<z)&(z<=y))) +
                         (Fz_0(z))*((y>=1)&((0<z)&(z<1))) +
                         1*(y>=1&z>=1)}

Fyz_1 <- function(y,z)  {(Fy(y))*(((0<y)&(y<1))&(z>=y)) +
                         (log(z)*(y+z)-log(y)*(z+2*y)+3*y-2*z)*(((0<y)&(y<1))&((y^2<=z)&(z<y))) + 
                         (z*(log(z)-log(y))+4*sqrt(z)-2*z-z/y)*(((0<y)&(y<1))&((0<z)&(z<y^2)))+
                         Fz_1(z)*((y>=1)&((0<z)&(z<1))) +
                         1*(y>=1&z>=1)}

Fyz_t <- function(y,z,t){(Fy(y))*(((0<y)&(y<1))&(z>=y^t)) +
                         (y*log(z) - (t+1)*(y*log(y) - y + z^(1/t)) - ((z*y^(1-t) - z^(1/t))/(1-t)) + y)*(((0<y)&(y<1))&((y^(t+1)<=z)&(z<y^t))) +
                         (t*(z^(1/(t+1)) -z^(1/t)) + (z^(1/t) - z*y^(1-t))/(1-t) + (z^(1/(t+1)) - z*y^(-t))/(t) + 2*z^(1/(t+1)) - z^(1/t))*(((0<y)&(y<1))&((0<z)&(z<y^(t+1)))) + 
                         Fz_t(z,t)*((y>=1)&((0<z)&(z<1))) +
                         1*(y>=1&z>=1)}


#(Y,Z|X=x)

Fyz1x_0 <- function(y,z,x) {(z*y/x^2)*(((0<z)&(z<x))&((0<y)&(y<x))) +
                             Fy1x(y,x)*((x<=z)&((0<y)&(y<x))) +
                             Fz1x_0(z,x)*(((0<z)&(z<x))&(x<=y)) +
                             1*(z>=x&y>=x)}
Fyz1x_1 <- function(y,z,x) {(z*log(y)-z*(log(z)-1)+z*log(x))*(((0<y)&(y<x))&((0<z)&(z<x*y))) +
                             Fy1x(y,x)*((x*y<=z)&((0<y)&(y<x))) + 
                             Fz1x_1(z,x)*(((0<z)&(z<x*y))&(x<=y)) +
                             1*(z>=x*y&y>=x)}
Fyz1x_t <- function(y,z,x,t){((1/(1-t))*((z*y^(1-t))/(x^2) - (t*z^(1/t))/(x^(1/t+1)) ))*(((0<y)&(y<x))&((0<z)&(z<x*y^t))) +
                             Fy1x(y,x)*((x*y^t<=z)&((0<y)&(y<x))) +
                             Fz1x_t(z,x,t)*(((0<z)&(z<x*y^t))&(x<=y)) +
                             1*(z>=x*y^t&y>=x)}

#(X,Z|Y)

Fxz1y <- function(x,z,y,t) {((-z/log(y))*(1/y^(t+1) - 1/(x*y^t)))*(((y<=x)&(x<1))&((0<=z)&(z<=y^(t+1)))) + 
                            (z/(x*y^t*log(y)) - 1/log(y) - log(z)/log(y) + (t+1))*(((y<=x)&(x<1))&((y^(t+1)<z)&(z<=x*y^t))) +
                            Fx1y(x,y)*(((y<=x)&(x<1))&(x*y^t<z)) +
                            Fz1y_t(z,y,t)*(x>=1) +
                            1*((x>=1)&(z>=y^t))}

#(X,Y|Z)

Fxy1z_0 <- function(x,y,z)  {((y/log(z))*(1/x - 1/z))*(((z<=x)&(x<=1))&((0<y)&(y<=z))) +
                             ((1/log(z))*(log(z)-log(y)+y/x-1))*(((z<=x)&(x<1))&((z<y)&(y<x))) +
                             Fy1z_0(y,z)*(1<=x) +
                             Fx1z_0(x,z)*(y>=x) +
                             1*(x>=1&y>=z)}
Fxy1z_1 <- function(x,y,z)  {((1/(log(z)+2/sqrt(z) -2))*((log(z)-log(y)-log(x)-1)/x + y/z))*(((sqrt(z)<=x)&(x<1))&((z/x<=y)&(y<sqrt(z)))) +
                             ((1/(log(z)+2/sqrt(z) -2))*((log(z)-log(y)-log(x)-1)/x + 2/sqrt(z)-1/y))*(((sqrt(z)<=x)&(x<1))&((sqrt(z)<=y)&(y<x))) +
                             Fy1z_1(y,z)*(x>=1) +
                             Fx1z_1(x,z)*(y>=x) +
                             1*(x>=1&y>=x)}
Fxy1z_t <- function(x,y,z,t){aux <- function(z,t)((t-t^2)/(z^(-t/(t+1)) - 1 + t^2*z^(1/t-1)*(1 - z^(-1/(t^2+t)))))
                            (aux(z,t)*(y/z - (1/(1-t))*(y^(1-t)/x - t*z^(1/t-1)/x^(1/t))))*(((z^(1/(t+1))<=x)&(x<1))&((z^1/t <=y)&(y<z^(1/(t+1))))) +
                            (aux(z,t)*(z^(-t/(t+1))*(1+1/t)-y^(-t)/t-(1/(1-t))*(y^(1-t)/x - t*z^(1/t-1)/x^(1/t))))*(((z^(1/(t+1))<=x)&(x<1))&((z^(1/(t+1))<=y)&(y<x))) +
                             Fy1z_t(y,z,t)*(x>=1) +
                             Fx1z_t(x,z,t)*((x<=y)&(y<1)) +
                             1*(x>=1&y>=1)} 
 
                            
##CÓPULAS!!
#FUNCIÓN SIMILAR A OUTER.
rat <- function(x,y,f) { 
 mat <- matrix( ,length(x), length(y))
 for(i in 1:length(x)) { 
 a <- rep(x[i], length(x))
 mat[ i, ] <- mapply(f, a, y)
 }
 mat
}

m <- 100
##(X,Y)
u <- seq(0.05, 0.99, length.out = m)
v <- seq(0.05, 0.99, length.out = m)


Cxy <- function(u,v){u*(((0<u)&(u<=invFy(v))&(invFy(v)<1))|(v == 1)) + 
                     (invFy(v)*(1+log(u)-log(invFy(v))))*((0 < invFy(v))&(invFy(v)<u)&(u<1)) +
                     v*(u == 1) +
                     1*(u>=1&v>=1) +
                     0*(u == 0 | v == 0)}
densCxy <- function(u,v){((-1/u)*(1/log(invFy(v))))*((0 < invFy(v))&(invFy(v)<u)&(u<1))}
dev.new()
simsCxy <- rat(u, v, Cxy)
image(u, v, simsCxy, main = "Cópula de (X,Y)", col = heat.colors(20))
simsdensCxy <- rat(u, v, densCxy)
dev.new()
image(u, v, simsdensCxy, main = "Densidad cópula de (X,Y)", col = heat.colors(20))

##(X,Z)
m <- 100
u <- seq(0.12, 0.99, length.out = m)
w <- seq(0.12, 0.99, length.out = m)


#Teta = 0

Cxz_0 <- Cxy
densCxz_0 <- densCxy
dev.new()
simsCxz_0 <- rat(u, v, Cxz_0)
image(u, v, simsCxz_0, main = "Cópula de (X,Z), teta = 0", col = heat.colors(20))
simsdensCxz_0 <- rat(u, v, densCxz_0)
dev.new()
image(u, v, simsdensCxz_0, main = "Densidad cópula de (X,Z), teta = 0", col = heat.colors(20))

#Teta = 1

Cxz_1 <- function(u,w){u*(((0<u^2)&(u^2<invFz_1(w))&(invFz_1(w)<1)) | (w == 1)) +
                      (4*sqrt(invFz_1(w)) + (invFz_1(w)*(log(invFz_1(w)) - 2*log(u)) - 3*invFz_1(w))/u)*((0 < invFz_1(w))&(invFz_1(w)<u^2)&(u^2<1)) +
                       w*(u == 1) +
                       1*(u>=1 & w>=1)}
densCxz_1 <- function(u,w){((2*log(u) - log(invFz_1(w)))/(u^2*(log(invFz_1(w)) + 2/sqrt(invFz_1(w)) - 2)))*((0<invFz_1(w))&(invFz_1(w)<u^2)&(u^2<1))}
dev.new()
simsCxz_1 <- rat(u, w, Cxz_1)
image(u, w, simsCxz_1, main = "Cópula de (X,Z), teta = 1", col = heat.colors(20))
simsdensCxz_1 <- rat(u, w, densCxz_1)
dev.new()
image(u, w, simsdensCxz_1, main = "Densidad cópula de (X,Z), teta = 1", col = heat.colors(20))

#Teta > 0 & teta != 1 

Cxz_t <- function(u,w,t){u*(((0<u^(t+1))&(u^(t+1)<invFz_t(w,t))&(invFz_t(w,t)<1)) | (w == 1)) +
                        ((1/(1-t))*((invFz_t(w,t)/t)*((invFz_t(w,t))^(-t/(t+1)) - u^(-t)) - (t^2*(invFz_t(w,t))^(1/t))*((invFz_t(w,t))^(-1/(t^2+t)) - u^(-1/t))) + (invFz_t(w,t))^(1/(t+1)))*((0 < invFz_t(w,t))&(invFz_t(w,t)<u^(t+1))&(u^(t+1)<1)) +
                         w*(u == 1) +
                         1*(u>=1 & w>=1)}
densCxz_t <- function(u,w,t){( (t*(1 - u^(t-1/t)*(invFz_t(w,t))^(1/t-1)))/(u^(t+1)*((invFz_t(w,t))^(-t/(t+1)) - 1 + t^2*(invFz_t(w,t))^(1/t-1)*(1 - (invFz_t(w,t))^(-1/(t^2+t))) )) )*((0 < invFz_t(w,t))&(invFz_t(w,t)<u^(t+1))&(u^(t+1)<1))}
simsCxz_.3 <- rat(u, w,function(u,w) Cxz_t(u,w,0.3))
dev.new()
image(u, w, simsCxz_.3, main = "Cópula de (X,Z), teta = 0.3", col = heat.colors(20))
simsdensCxz_.3 <- rat(u, w,function(u,w) densCxz_t(u,w,0.3))
dev.new()
image(u, w, simsdensCxz_.3, main = "Densidad cópula de (X,Z), teta = 0.3", col = heat.colors(20))

u <- seq(0.19, 0.99, length.out = m)
w <- seq(0.19, 0.99, length.out = m)


simsCxz_1.5 <- rat(u, w,function(u,w) Cxz_t(u,w,1.5))
dev.new()
image(u, w, simsCxz_1.5, main = "Cópula de (X,Z), teta = 1.5", col = heat.colors(20))
simsdensCxz_1.5 <- rat(u, w,function(u,w) densCxz_t(u,w,1.5))
dev.new()
image(u, w, simsdensCxz_1.5, main = "Densidad cópula de (X,Z), teta = 1.5", col = heat.colors(20))

##(Y,Z)

#teta = 0
m <- 100
v <- seq(0.1, 0.99, length.out = m)
w <- seq(0.1, 0.99, length.out = m)

Cyz_0 <- function(v,w){(invFy(v)*(log(invFz_0(w)) - log(invFy(v)) + 2 - invFz_0(w)))*( ((0<invFy(v))&(invFy(v)<1)) & ((invFy(v)<invFz_0(w))&(invFz_0(w)<1)) ) +
                       v*(((0<invFy(v))&(invFy(v)<1))&(invFz_0(w)>=1)) +
                       (invFz_0(w)*(log(invFy(v)) - log(invFz_0(w)) + 2 - invFy(v)))*( ((0<invFy(v))&(invFy(v)<1)) & ((0<invFz_0(w))&(invFz_0(w)<=invFy(v)))) +
                       w*((invFy(v)>=1)&((0<invFz_0(w))&(invFz_0(w)<1))) +
                       1*(invFz_0(w)>=1&invFy(v)>=1)}
densCyz_0 <- function(v,w){(((1/log(invFz_0(w)))*(1/invFz_0(w) - 1))*(1/log(invFy(v))))*(((0 < invFy(v))&(invFy(v) < 1)) & ((invFy(v) < invFz_0(w))&(invFz_0(w) < 1))) + 
                           (((1/log(invFz_0(w)))*(1/invFy(v)   - 1))*(1/log(invFy(v))))*(((0 < invFy(v))&(invFy(v) < 1)) & ((invFz_0(w) <= invFy(v))))}
simsCyz_0 <- rat(v, w, Cyz_0)
dev.new()
image(v, w, simsCyz_0, main = "Cópula de (Y,Z), teta = 0", col = heat.colors(20))
simsdensCyz_0 <- rat(v, w,densCyz_0)
dev.new()
image(v, w, simsdensCyz_0, main = "Densidad cópula de (Y,Z), teta = 0", col = heat.colors(20))

#teta = 1
m <- 100
v <- seq(0.1, 0.99, length.out = m)
w <- seq(0.1, 0.99, length.out = m)

Cyz_1 <- function(v,w){(v)*(((0 < invFy(v))&(invFy(v) < 1)) & ( invFz_1(w)>=invFy(v) )) +
                       (log(invFz_1(w))*(invFy(v) + invFz_1(w)) - log(invFy(v))*(invFz_1(w) + 2*invFy(v)) + 3*invFy(v) - 2*invFz_1(w) )*(((0 < invFy(v))&(invFy(v) < 1)) & ((invFy(v)^2 < invFz_1(w))&(invFz_1(w) < invFy(v)))) +
                       (invFz_1(w)*(log(invFz_1(w)) - log(invFy(v))) + 4*sqrt(invFz_1(w)) - 2*invFz_1(w) - invFz_1(w)/invFy(v) )*(((0 < invFy(v))&(invFy(v) < 1)) & ((0 < invFz_1(w))&(invFz_1(w) < invFy(v)^2))) +
                       (w)*(((0 < invFz_1(w))&(invFz_1(w) < 1)) & (invFy(v)>=1)) +
                       (1)*(invFy(v)>=1 & invFz_1(w)>=1)}
densCyz_1 <- function(v,w) ((1/(log(invFy(v))*(log(invFz_1(w)) + 2/sqrt(invFz_1(w)) - 2)))*(1/invFy(v) - 1/invFz_1(w)))*(((0<invFy(v))&(invFy(v)<1)) & ((invFy(v)^2<invFz_1(w))&(invFz_1(w)<invFy(v)))) + ((1/(log(invFy(v))*(log(invFz_1(w)) + 2/sqrt(invFz_1(w)) - 2)))*(1/invFy(v) - 1/invFy(v)^2))*(((0<invFy(v))&(invFy(v)<1)) & (invFz_1(w) < invFy(v)^2)) 
simsCyz_1 <- rat(v, w, Cyz_1)
dev.new()
image(v, w, simsCyz_1, main = "Cópula de (Y,Z), teta = 1", col = heat.colors(20))
simsdensCyz_1 <- rat(v, w,densCyz_1)
dev.new()
image(v, w, simsdensCyz_1, main = "Densidad cópula de (Y,Z), teta = 1", col = heat.colors(20))

#Teta > 0 & teta != 1

Cyz_t <- function(v,w,t){(v)*(((0 < invFy(v))&(invFy(v) < 1)) & ((invFy(v))^(t) <= invFz_t(w,t))) + 
                         (invFy(v)*log(invFz_t(w,t)) - (t+1)*(invFy(v)*log(invFy(v)) - invFy(v) + (invFz_t(w,t))^(1/t))  - ((invFz_t(w,t)*(invFy(v))^(1-t) -(invFz_t(w,t))^(1/t) )/(1-t)) + invFy(v))*(((0 < invFy(v))&(invFy(v) < 1)) & ( ((invFy(v))^(t+1) < invFz_t(w,t)) & (invFz_t(w,t) < (invFy(v))^(t)) )) + 
                         (t*((invFz_t(w,t))^(1/(t+1)) - (invFz_t(w,t))^(1/t)) + ((invFz_t(w,t))^(1/t) - invFz_t(w,t)*(invFy(v))^(1-t))/(1-t) + ((invFz_t(w,t))^(1/(t+1)) - invFz_t(w,t)*(invFy(v))^(-t))/t + 2*(invFz_t(w,t))^(1/(t+1)) - (invFz_t(w,t))^(1/t))*(((0 < invFy(v))&(invFy(v) < 1)) & (invFz_t(w,t) < (invFy(v))^(t+1))) +
                         (w)*((invFy(v) >= 1) & ((0 < invFz_t(w,t)) & (invFz_t(w,t) < 1))) +
                         (1)*(invFz_t(w,t)>=1 & invFy(v)>=1)}

densCyz_t <- function(v,w,t){(1/((log(invFy(v))*( (1/(1-t))*(((1-t^2)/t)*(invFz_t(w,t))^(-t/(t+1)) - 1/t + t*(invFz_t(w,t))^(1/t-1)) ))))*(1/invFy(v)^t - 1/invFz_t(w,t))*(((0 < invFy(v))&(invFy(v) < 1)) & (((invFy(v))^(t+1) <= invFz_t(w,t))&(invFz_t(w,t) < invFy(v)^t))) +
                             (1/((log(invFy(v))*( (1/(1-t))*(((1-t^2)/t)*(invFz_t(w,t))^(-t/(t+1)) - 1/t + t*(invFz_t(w,t))^(1/t-1))))))*(1/invFy(v)^t - 1/(invFy(v))^(t+1))*(((0 < invFy(v))&(invFy(v) < 1)) & (invFz_t(w,t) < (invFy(v))^(t+1)))}
simsCyz_0.5 <- rat(v, w, function(v,w)Cyz_t(v,w,0.5))
dev.new()
image(v, w, simsCyz_0.5, main = "Cópula de (Y,Z), teta = 0.5", col = heat.colors(20))
simsdensCyz_0.5 <- rat(v, w,function(v,w) densCyz_t(v,w,0.5))
dev.new()
image(v, w, simsdensCyz_0.5, main = "Densidad cópula de (Y,Z), teta = 0.5", col = heat.colors(20))

m <- 100
v <- seq(0.109, 0.99, length.out = m)
w <- seq(0.109, 0.99, length.out = m)


simsCyz_1.5 <- rat(v, w, function(v,w)Cyz_t(v,w,1.5))
dev.new()
image(v, w, simsCyz_1.5, main = "Cópula de (Y,Z), teta = 1.5", col = heat.colors(20))
simsdensCyz_1.5 <- rat(v, w,function(v,w) densCyz_t(v,w,1.5))
dev.new()
image(v, w, simsdensCyz_1.5, main = "Densidad cópula de (Y,Z), teta = 1.5", col = heat.colors(20))


#(X,Y|Z)
#teta = 0
m <- 100
u <- seq(0.01, 0.99, length.out = m)
v <- seq(0.01, 0.99, length.out = m)

Cxy1z_0 <- function(u,v,z){(0)*(u==0 | v==0) +
                          ((v/(z-1))*(z^u-1))*(((0 < u)&(u <=1 )) & ((0 < v)&(v < (z-1)/log(z)))) +
                          ( (1/log(z))*(log(z) - log(invFy1z_0(v,z)) + invFy1z_0(v,z)/(z^(1-u)) - 1) )*(((0 < u)&(u < 1))&(((z-1)/log(z) < v)&(v <u + (z^(1-u)-1)/log(z)))) +
                          (v)*((u == 1) & ((0 < v)&(v < 1)))+                       
                          (u)*(((0 < u)&(u < 1))&(v >= u + (z^(1-u)-1)/log(z))) +
                          (1)*(u==1&v==1)}
densCxy1z_0 <- function(u,v,z){(log(z)*z^u/(z-1))*(((0 < u)&(u < 1)) & ((0 < v)&(v <= (z-1)/log(z)))) +
                              ((z^(u-1)*log(z)*invFy1z_0(v,z))/(invFy1z_0(v,z)-1))*(((0 < u)&(u < 1)) & (((z-1)/log(z) < v)&(v < u+(z^(1-u)-1)/log(z))))}
simsCxy1z_0_0.3 <- rat(u, v, function(u,v)Cxy1z_0(u,v,0.3))
dev.new()
image(u, v, simsCxy1z_0_0.3, main = "Cópula de (X,Y|Z), teta = 0, z = 0.3", col = heat.colors(20))
simsdensCxy1z_0_0.3 <- rat(u, v,function(u,v) densCxy1z_0(u,v,0.3))
dev.new()
image(u, v, simsdensCxy1z_0_0.3, main = "Densidad cópula de (X,Y|Z), teta = 0, z = 0.3", col = heat.colors(20))

#teta = 1
m <- 100
u <- seq(0.05, 1, length.out = m)
v <- seq(0.05, 1, length.out = m)

Cxy1z_1 <- function(u,v,z){(0)*(u == 0|v == 0) +
                           ((1/(log(z) + 2/sqrt(z) - 2)) * ((log(z) - log(invFy1z_1(v,z)) - log(invFx1z_1(u,z)) -1)/invFx1z_1(u,z) + invFy1z_1(v,z)/z))*(((sqrt(z) < invFx1z_1(u,z))&(invFx1z_1(u,z) < 1)) & ((z/invFx1z_1(u,z) < invFy1z_1(v,z))&(invFy1z_1(v,z) < sqrt(z)))) +
                           (v)*(invFx1z_1(u,z) == 1 & ((z < invFy1z_1(v,z))&(invFy1z_1(v,z) < 1))) +
                           ((1/(log(z) + 2/sqrt(z) - 2)) * ((log(z) - log(invFy1z_1(v,z)) - log(invFx1z_1(u,z)) -1)/invFx1z_1(u,z) + 2/sqrt(z) - 1/invFy1z_1(v,z)))*(((sqrt(z) < invFx1z_1(u,z))&(invFx1z_1(u,z) < 1)) & ((sqrt(z) < invFy1z_1(v,z))&(invFy1z_1(v,z) < invFx1z_1(u,z)))) +
                           (u)*(((sqrt(z) < invFx1z_1(u,z))&(invFx1z_1(u,z) < 1)) & (invFx1z_1(u,z) <= invFy1z_1(v,z))) +
                           (1)*(u == 1& v == 1)}
densCxy1z_1 <- function(u,v,z){( ((1/(log(z) + 2/sqrt(z) - 2))*z*invFy1z_1(v,z))/(invFy1z_1(v,z)*(2*log(invFx1z_1(u,z)) - log(z))*(invFy1z_1(v,z) - z)) )*(((sqrt(z) < invFx1z_1(u,z))&(invFx1z_1(u,z) < 1)) & ((z/invFx1z_1(u,z) < invFy1z_1(v,z))&(invFy1z_1(v,z) < sqrt(z)))) +
                             ( ((1/(log(z) + 2/sqrt(z) - 2))*invFy1z_1(v,z))/((2*log(invFx1z_1(u,z)) - log(z))*(1 - invFy1z_1(v,z))) )*(((sqrt(z) < invFx1z_1(u,z))&(invFx1z_1(u,z) < 1)) & ((sqrt(z) < invFy1z_1(v,z))&(invFy1z_1(v,z) < invFx1z_1(u,z))))}
simsCxy1z_1_0.2 <- rat(u, v, function(u,v)Cxy1z_1(u,v,0.2))
dev.new()
image(u, v, simsCxy1z_1_0.2, main = "Cópula de (X,Y|Z), teta = 1, z = 0.2", col = heat.colors(20))

simsdensCxy1z_1_0.2 <- rat(u, v, function(u,v) densCxy1z_1(u,v,0.2))
dev.new()
image(u, v, simsdensCxy1z_1_0.2, main = "Densidad cópula de (X,Y|Z), teta = 1, z = 0.2", col = heat.colors(20))

#Z = 0.5
m <- 100
u <- seq(0.01, 0.990, length.out = m)
v <- seq(0.01, 0.990, length.out = m)

simsCxy1z_1_0.5 <- rat(u, v, function(u,v)Cxy1z_1(u,v,0.5))
dev.new()
image(u, v, simsCxy1z_1_0.5, main = "Cópula de (X,Y|Z), teta = 1, z = 0.5", col = heat.colors(20))

simsdensCxy1z_1_0.5 <- rat(u, v, function(u,v) densCxy1z_1(u,v,0.5))
dev.new()
image(u, v, simsdensCxy1z_1_0.5, main = "Densidad cópula de (X,Y|Z), teta = 1, z = 0.5", col = heat.colors(20))

#Z = 0.8
m <- 100
u <- seq(0.022, 0.98, length.out = m)
v <- seq(0.022, 0.98, length.out = m)

simsCxy1z_1_0.8 <- rat(u, v, function(u,v)Cxy1z_1(u,v,0.8))
dev.new()
image(u, v, simsCxy1z_1_0.8, main = "Cópula de (X,Y|Z), teta = 1, z = 0.8", col = heat.colors(20))

simsdensCxy1z_1_0.8 <- rat(u, v, function(u,v) densCxy1z_1(u,v,0.8))
dev.new()
image(u, v, simsdensCxy1z_1_0.8, main = "Densidad cópula de (X,Y|Z), teta = 1, z = 0.8", col = heat.colors(20))

# teta > 0 & teta != 1
m <- 100
u <- seq(0.001, 0.990, length.out = m)
v <- seq(0.001, 0.990, length.out = m)

Cxy1z_t <- function(u,v,z,t){(0)*(u == 0 | v == 0) +
                             (((t - t^2)/(z^(-t/(t+1)) - 1 + t^2*z^(1/t-1)*(1 - z^(-1/(t^2+t)))))*(invFy1z_t(v,z,t)/z - (1/(1-t))*(invFy1z_t(v,z,t)^(1-t)/invFx1z_t(u,z,t) - t*z^(1/t-1)/invFx1z_t(u,z,t)^(1/t))))*(((z^(1/(t+1)) < invFx1z_t(u,z,t))&(invFx1z_t(u,z,t) < 1)) & (((z/invFx1z_t(u,z,t))^(1/t) < invFy1z_t(v,z,t))&(invFy1z_t(v,z,t) < z^(1/(t+1))))) +
                             (((t - t^2)/(z^(-t/(t+1)) - 1 + t^2*z^(1/t-1)*(1 - z^(-1/(t^2+t)))))*(z^(-t/(t+1))*(1+1/t) - invFy1z_t(v,z,t)^(-t)/t - (1/(1-t))*(invFy1z_t(v,z,t)^(1-t)/invFx1z_t(u,z,t) - t*z^(1/t-1)/invFx1z_t(u,z,t)^(1/t))))*(((z^(1/(t+1)) < invFx1z_t(u,z,t))&(invFx1z_t(u,z,t) < 1)) & ((z^(1/(t+1)) < invFy1z_t(v,z,t))&(invFy1z_t(v,z,t) < invFx1z_t(u,z,t)))) +
                             (v)*(u == 1) +
                             (u)*((invFx1z_t(u,z,t) < invFy1z_t(v,z,t)) & (invFy1z_t(v,z,t) < 1)) +
                             (1)*(u == 1 & v == 1)}
densCxy1z_t <- function(u,v,z,t){((z^(-t/(t + 1)) - 1 + t^2*z^(1/t - 1)*(1 - z^(-1/(t^2 + t))))*(invFx1z_t(u,z,t))^(1/t - 1)*(z))/((t)*(invFx1z_t(u,z,t)^(1/t - t) - z^(1/t - 1))*(invFy1z_t(v,z,t)^t - z))*(((z^(1/(t + 1)) < invFx1z_t(u,z,t)) & (invFx1z_t(u,z,t) < 1)) & (((z/invFx1z_t(u,z,t))^(1/t) < invFy1z_t(v,z,t)) & (invFy1z_t(v,z,t) < z^(1/(t + 1)))))+ ((z^(-t/(t + 1)) - 1 + t^2*z^(1/t - 1)*(1 - z^(-1/(t^2 + t))))*(invFx1z_t(u,z,t)^(1/t - 1))*(invFy1z_t(v,z,t)))/((t)*(invFx1z_t(u,z,t)^(1/t - t) - z^(1/t - 1))*(1 - invFy1z_t(v,z,t)))*(((z^(1/(t + 1)) < invFx1z_t(u,z,t)) & (invFx1z_t(u,z,t) < 1)) & ((z^(1/(t + 1)) < invFy1z_t(v,z,t)) & (invFy1z_t(v,z,t) < invFx1z_t(u,z,t))))}
simsCxy1z_1.5_0.5 <- rat(u, v, function(u,v)Cxy1z_t(u,v,0.5,1.5))
dev.new()
image(u, v, simsCxy1z_1.5_0.5, main = "Cópula de (X,Y|Z), teta = 1.5, z = 0.5", col = heat.colors(20))

simsdensCxy1z_1.5_0.5 <- rat(u, v, function(u,v) densCxy1z_t(u,v,0.5,1.5))
dev.new()
image(u, v, simsdensCxy1z_1.5_0.5, main = "Densidad cópula de (X,Y|Z), teta = 1.5, z = 0.5", col = heat.colors(20))


#(X,Z|Y)
m <- 100
u <- seq(0.001, 0.990, length.out = m)
w <- seq(0.001, 0.990, length.out = m)

Cxz1y_t <- function(u,w,y,t) {(0)*(u == 0 | w == 0) +
                              (((y^u - 1)/(y - 1))*w)*(((0 < u)&(u < 1)) & ((0 < w)&(w < (y - 1)/log(y)))) +
                              ((((y^(u - t - 1))*invFz1y_t(w,y,t) - log(invFz1y_t(w,y,t)) - 1)/log(y)) + t + 1)*(((0 < u)&(u < 1)) & (((y - 1)/log(y) < w)&(w <= (u + (y^(1 - u) - 1)/log(y))))) +
                              (u)*((((0 < u)&(u < 1)) & (((y^(1-u) - 1)/log(y)+ u) < w)) | (w == 1)) +
                              (w)*(u == 1) +
                              (1)*(u == 1 & w == 0)}
densCxz1y_t <- function(u,w,y,t){((log(y)*y^u)/(y - 1))*(((0 < u)&(u < 1)) & ((0 < w)&(w < (y - 1)/log(y)))) +
                                 ((log(y)*invFz1y_t(w,y,t)*y^(u-1))/(invFz1y_t(w,y,t) - y^t))*(((0 < u)&(u < 1)) & (((y - 1)/log(y) < w)&(w <= (u + (y^(1 - u) - 1)/log(y)))))}
simsCxz1y_0_0.2 <- rat(u, w, function(u,w) Cxz1y_t(u,w,0.2,0))
dev.new()
image(u, w, simsCxz1y_0_0.2, main = "Cópula de (X,Z|Y), teta = 0, y = 0.2", col = heat.colors(10))


simsdensCxz1y_0_0.2 <- rat(u, w, function(u,w) densCxz1y_t(u,w,0.2,0))
dev.new()
image(u, w, simsdensCxz1y_0_0.2, main = "Cópula de (X,Z|Y), teta = 0, y = 0.2", col = heat.colors(10))


simsCxz1y_1_0.2 <- rat(u, w, function(u,w) Cxz1y_t(u,w,0.2,1.5))
dev.new()
image(u, w, simsCxz1y_1_0.2, main = "Cópula de (X,Z|Y), teta = 1, y = 0.2", col = heat.colors(10))

simsdensCxz1y_1_0.2 <- rat(u, w, function(u,w) densCxz1y_t(u,w,0.2,1.5))
dev.new()
image(u, w, simsdensCxz1y_1_0.2, main = "Cópula de (X,Z|Y), teta = 1, y = 0.2", col = heat.colors(10))


#(Y,Z|X)
#teta = 0
m <- 100
v <- seq(0.001, 0.990, length.out = m)
w <- seq(0.001, 0.990, length.out = m)

Cyz1x_0 <- function(v,w,x){(0)*(v == 0 | w == 0) +
                           (w*v)*(((0 < v)&(v < 1)) & ((0 < w)&(w < 1))) +
                           (v)*(((0 < v)&(v < 1)) & (w == 1)) +
                           (w)*((v == 1) & ((0 < w)&(w < 1))) +
                           (1)*(v == 1 & w == 1)}
densCyz1x_0 <- function(v,w,x){1}
simsCxz1y_0_0.2 <- rat(v, w, function(v,w) Cyz1x_0(v,w,0.2))
dev.new()
image(v, w, simsCxz1y_0_0.2, main = "Cópula de (Y,Z|X), teta = 0, x = 0.2", col = heat.colors(10))


simsdensCyz1x_0_0.2 <- rat(v, w, function(v,w) densCyz1x_0(v,w,0.2))
dev.new()
image(v, w, simsdensCyz1x_0_0.2, main = "Densidad cópula de (Y,Z|X), teta = 0, x = 0.2", col = heat.colors(10))

#teta = 1
m <- 300
v <- seq(0.03, 0.990, length.out = m)
w <- seq(0.03, 0.990, length.out = m)


Cyz1x_1 <- function(v,w,x){(0)*(v == 0 | w == 0) +
                           ((invFz1x_1(w,x)/x^2)*(log(v) + 2*log(x) - log(invFz1x_1(w,x)) + 1))*(((0 < v)&(v < 1)) & ((0 < w)&(w < (invFy1x(v,x)/x)*(log(x) + 1 - log(invFy1x(v,x)))))) +
                           (v)*(((0 < v)&(v < 1)) & ((invFy1x(v,x)/x)*(log(x) + 1 - log(invFy1x(v,x))) <= w)) +
                           (w)*(v == 1 & ((0 < w)&(w < ((invFy1x(v,x)/x)*(log(x) + 1 - log(invFy1x(v,x))))))) +
                           (1)*(v == 1 & w == 1)}
densCyz1x_1 <- function(v,w,x){(1/(v*(2*log(x) - log(invFz1x_1(w,x)))))*(((0 < v)&(v < 1)) & ((0 < w)&(w < (invFy1x(v,x)/x)*(log(x) + 1 - log(invFy1x(v,x))))))}
simsCxz1y_1_0.2 <- rat(v, w, function(v,w) Cyz1x_1(v,w,0.2))
dev.new()
image(v, w, simsCxz1y_1_0.2, main = "Cópula de (Y,Z|X), teta = 1, x = 0.2", col = heat.colors(10))


simsdensCyz1x_1_0.2 <- rat(v, w, function(v,w) densCyz1x_1(v,w,0.2))
dev.new()
image(v, w, simsdensCyz1x_1_0.2, main = "Densidad cópula de (Y,Z|X), teta = 1, x = 0.2", col = heat.colors(10))

#teta > 0 & teta != 1

Cyz1x_t <- function(v,w,x,t){(0)*(v == 0 | w == 0) +
                             ((1/(1 - t))*((invFz1x_t(w,x,t)*v^(1 - t))/x^(t + 1) - (t*invFz1x_t(w,x,t)^(1/t))/x^(1/t + 1)))*(((0 < v)&(v < 1)) & ((0 < w)&(w < ((v^t - t*v)/(1 - t))))) +
                             (v)*(((0 < v)&(v < 1)) & (((v^t - t*v)/(1 - t)) <= w)) +
                             (w)*(v == 1 & ((0 < w)&(w < 1))) +
                             (1)*(v == 1 & w == 1)}
densCyz1x_t <- function(v,w,x,t){( ((1 - t)*x^(1/t - t))/(v^t*(x^(1/t - t) - invFz1x_t(w,x,t)^(1/t - 1))) )*(((0 < v)&(v < 1)) & ((0 < w)&(w < ((v^t - t*v)/(1 - t)))))}

simsCxz1y_1.5_0.2 <- rat(v, w, function(v,w) Cyz1x_t(v,w,0.2,1.5))
dev.new()
image(v, w, simsCxz1y_1.5_0.2, main = "Cópula de (Y,Z|X), teta = 1.5, x = 0.2", col = heat.colors(10))


simsdensCyz1x_1.5_0.2 <- rat(v, w, function(v,w) densCyz1x_t(v,w,0.2,1.5))
dev.new()
image(v, w, simsdensCyz1x_1.5_0.2, main = "Densidad cópula de (Y,Z|X), teta = 1.5, x = 0.2", col = heat.colors(10))



simsCxz1y_0.01_0.2 <- rat(v, w, function(v,w) Cyz1x_t(v,w,0.2,0.01))
dev.new()
image(v, w, simsCxz1y_0.5_0.2, main = "Cópula de (Y,Z|X), teta = 0.5, x = 0.2", col = heat.colors(10))


simsdensCyz1x_0.1_0.2 <- rat(v, w, function(v,w) densCyz1x_t(v,w,0.2,0.1))
dev.new()
image(v, w, simsdensCyz1x_0.1_0.2, main = "Densidad cópula de (Y,Z|X), teta = 0.5, x = 0.2", col = heat.colors(10))

densCopTri1 <- function(u,v,w){densCyz1x_1(v,w,0.3)*densCxy(u,v)*densCxz_1(u,w)}
densCxyz_1 <- function(u,v,w) {(-1/(u^{2}*log(invFy(v))*invFy(v)*(log(invFz_1(w)) + 2/sqrt(invFz_1(w)) - 2)))}

rat3d <- function(x,y,z,f){
 mat1 <- matrix( , length(x)^2, length(x))
 i <- 1
 while(i < length(x)^2 + 1){
      for(j in 1:length(x)){
          a <- rep(x[j], length(x))
          for(k in 1:length(y)){
              b <- rep(y[k], length(x))
              mat1[i,] <- mapply(f, a,b,z)
              i <- i+1
           }
       }
 }
 mat1
}
