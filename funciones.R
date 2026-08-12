library(geobr)
library(ggplot2)
library(sf)
library(dplyr)
library(writexl)
library(maps)
library(fields)
library(leaflet)
library(RColorBrewer)
library(spdep)
library(latex2exp)
library(ggcorrplot)
library(invgamma)
library(mvtnorm)

library(MASS)
library(GGally)
library(network)
library(sna)
library(ggplot2)
library(invgamma)
library(bnlearn)
library(BNSL)
library(xtable)
library(corrplot)
library(expm)
library(truncnorm)
library(matrixcalc)
library(pracma)
library(MCMCpack)
library(directlabels)
library(reshape)
library(ppcor)
library(LaplacesDemon)
library(VGAM)


CorrConBX<-function(A, Sigma) #genera boxplot de la correlaciones condicionales de la muestra segun el grafo A y la matriz de covarianza
{
n=dim(A)[1]	
QQ<-ginv(Sigma)

i=1
j=1
Qij<-matrix(ncol=n, nrow=n)
for(i in 1:n){
  for(j in 1:n){ 		
    Qij[i, j]<--QQ[i,j]/(sqrt(QQ[i,i]*QQ[j,j]))    
  }
}

CormtvS1orden<- Qij[which(A == 1, arr.ind = TRUE)]

CormtvS2orden<- Qij[which((A%^%2) == 1, arr.ind = TRUE)]

CormtvS3orden<- Qij[which((A%^%3) == 1, arr.ind = TRUE)]

CormtvS4orden<- Qij[which((A%^%4) == 1, arr.ind = TRUE)]

CormtvS5orden<- Qij[which((A%^%5) == 1, arr.ind = TRUE)]

CormtvS6orden<- Qij[which((A%^%6) == 1, arr.ind = TRUE)]

CormtvS7orden<- Qij[which((A%^%7) == 1, arr.ind = TRUE)]

CormtvS8orden<- Qij[which((A%^%8) == 1, arr.ind = TRUE)]


bxplot<-boxplot(CormtvS1orden, CormtvS2orden,CormtvS3orden, CormtvS4orden,CormtvS5orden, CormtvS6orden,CormtvS7orden, CormtvS8orden,
ylab = "Conditional correlation", xlab ="
Neighbor order", 
at = c(1,2,3,4,5,6,7,8),
names = c("1","2", "3","4","5","6","7","8"),
las = 1,
col = c('powderblue', 'mistyrose', 'powderblue', 'mistyrose','powderblue', 'mistyrose','powderblue', 'mistyrose')
)

return(bxplot)

}

CorrBX<-function(x, A) #genera boxplot de la correlaciones marginales de la muestra x y segun el grafo A
{
CORR<-cor(x)
	
CormtvS1orden<- CORR[which(A == 1, arr.ind = TRUE)]

CormtvS2orden<-  CORR[which((A%^%2) == 1, arr.ind = TRUE)]

CormtvS3orden<-  CORR[which((A%^%3) == 1, arr.ind = TRUE)]

CormtvS4orden<-  CORR[which((A%^%4) == 1, arr.ind = TRUE)]

CormtvS5orden<- CORR[which((A%^%5) == 1, arr.ind = TRUE)]

CormtvS6orden<- CORR[which((A%^%6) == 1, arr.ind = TRUE)]

CormtvS7orden<- CORR[which((A%^%7) == 1, arr.ind = TRUE)]

CormtvS8orden<- CORR[which((A%^%8) == 1, arr.ind = TRUE)]


bxplot<-boxplot(CormtvS1orden, CormtvS2orden,CormtvS3orden, CormtvS4orden,CormtvS5orden, CormtvS6orden,CormtvS7orden, CormtvS8orden,
ylab = "Correlacion", xlab ="Orden vecino", 
at = c(1,2,3,4,5,6,7,8),
names = c("1","2", "3","4","5","6","7","8"),
las = 1,
col = c('powderblue', 'mistyrose', 'powderblue', 'mistyrose','powderblue', 'mistyrose','powderblue', 'mistyrose')
)
}


ComCarMTVv1<-function(x, A, gam)
{	
CORR_mtv<-cov2cor(x)	
	
d<-rowSums(A)
D<-diag(d)
Cov_CAr<-solve(D-gam*A)

theta_CAR <-rmvnorm(9000, mean_t, sigma= 10*Cov_CAr)

CORR_CAR<-cor(theta_CAR)


CormtvS1orden<-CORR_mtv[which(A == 1, arr.ind = TRUE)]

CormtvS2orden<-CORR_mtv[which((A%^%2) == 1, arr.ind = TRUE)]


CorcarS1orden<-CORR_CAR[which(A == 1, arr.ind = TRUE)]

CorcarS2orden<-CORR_CAR[which((A%^%2) == 1, arr.ind = TRUE)]

 
bx1<-boxplot(CormtvS1orden, CorcarS1orden,CormtvS2orden, CorcarS2orden, ylab = "Conditional correlation", xlab ="
Neighbor order", 
at = c(1,2,3,4),
names = c("1","1", "2","2"),
las = 1,
col = c('powderblue', 'mistyrose', 'powderblue', 'mistyrose')
)


 return(bx1)
} #genera boxplot de la correlaciones marginales de la muestra x y segun el grafo A compara con CAR, ESTO para vecinos de primeros ordenes


ComCarMTVv2<-function(x, A, gam)
{	
CORR_mtv<-cov2cor(x)	
d<-rowSums(A)
D<-diag(d)
Cov_CAr<-solve(D-gam*A)

theta_CAR <-rmvnorm(9000, mean_t, sigma= 10*Cov_CAr)

CORR_CAR<-cor(theta_CAR)


CormtvS3orden<-CORR_mtv[which((A%^%3) == 1, arr.ind = TRUE)]

CormtvS4orden<-CORR_mtv[which((A%^%4) == 1, arr.ind = TRUE)]


CorcarS3orden<-CORR_CAR[which((A%^%3) == 1, arr.ind = TRUE)]

CorcarS4orden<-CORR_CAR[which((A%^%4) == 1, arr.ind = TRUE)]



bx2<-boxplot(CormtvS3orden, CorcarS3orden , CormtvS4orden, CorcarS4orden, ylab = "Correlation", xlab ="Neighbor order", 
at = c(1,2,3,4),
names = c("3", "3", "4","4"),
las = 1,
col = c('powderblue', 'mistyrose', 'powderblue', 'mistyrose')
)

return(bx2)
} #genera boxplot de la correlaciones marginales de la muestra x y segun el grafo A compara con CAR, ESTO para vecinos de ultmos  ordenes


comprPosi<-function(Sigma){
for(i in 1:15){
while(is.positive.definite(round(Sigma,i))==TRUE){
return(i)
}
}
}
#verifica la aproximacion para que la matriz sea definida positiva, lo que pasa es que no da simetrica, a pesar de lo sea 


CorrDCovBX<-function(x, A, ymax) #genera boxplot de la correlaciones marginales de la matrix de covarianzax y segun el grafo A
{
CORR<-cov2cor(x)
	
CormtvS1orden<- CORR[which(A == 1, arr.ind = TRUE)]

CormtvS2orden<-  CORR[which((A%^%2) == 1, arr.ind = TRUE)]

CormtvS3orden<-  CORR[which((A%^%3) == 1, arr.ind = TRUE)]

CormtvS4orden<-  CORR[which((A%^%4) == 1, arr.ind = TRUE)]

CormtvS5orden<- CORR[which((A%^%5) == 1, arr.ind = TRUE)]

CormtvS6orden<- CORR[which((A%^%6) == 1, arr.ind = TRUE)]

CormtvS7orden<- CORR[which((A%^%7) == 1, arr.ind = TRUE)]

CormtvS8orden<- CORR[which((A%^%8) == 1, arr.ind = TRUE)]


bxplot<-boxplot(CormtvS1orden, CormtvS2orden,CormtvS3orden, CormtvS4orden,CormtvS5orden, CormtvS6orden,CormtvS7orden, CormtvS8orden, ylab = "Correlation", xlab ="Neighbor order", 
at = c(1,2,3,4,5,6,7,8),
names = c("1","2", "3","4","5","6","7","8"),
las = 1,
col = c('powderblue', 'mistyrose', 'powderblue', 'mistyrose','powderblue', 'mistyrose','powderblue', 'mistyrose')
)
}


#graficas para datos perdidos
lineasCondSombras<-function(theta_mtv, i,j, Covr, mean_t, nu, nlines){
	
n=dim(Covr)[1]
posi<-comprPosi(Covr)
QQ<-solve(round(Covr,posi))
theta_ij<-theta_mtv
theta_ij[i]<-NA
theta_ij[j]<-NA

mean_tcond<-(Covr[c(i,j), -c(i,j)]%*%ginv(Covr[-c(i,j), -c(i,j)])%*%(theta_mtv[-c(i,j)]-mean_t[-c(i,j)]))+mean_t[c(i,j)]

  beta<-t(theta_mtv[-c(i,j)]-mean_t[-c(i,j)])%*%QQ[-c(i,j), -c(i,j)]%*%(theta_mtv[-c(i,j)]-mean_t[-c(i,j)])
     
    CovCondi<-(as.numeric((nu+beta-2)/(nu+n-4)))*(ginv(QQ[c(i,j), c(i,j)])) 

nsim=nlines

theta_mtv_cond1<-matrix(nrow=nsim, ncol=2)

t=1
for(t in 1:nsim){
theta_mtv_cond<-rmvt(10000, mean_tcond, sig*round(CovCondi,2), n)
theta_mtv_cond1[t,]<-colMeans(theta_mtv_cond)
}

theta_ijCon<-matrix(ncol=nsim, nrow=n)

t=1

for(t in 1:nsim){
	theta_ijCon[,t]<-theta_ij
	theta_ijCon[i,t]<-theta_mtv_cond1[t,1]
	theta_ijCon[j,t]<-theta_mtv_cond1[t,2]
}


theta_ijConarri<-NULL
theta_ijConabajo<-NULL
k=1
for(k in 1:nsim){	if(theta_ijCon[i,k]>mean_tcond[1]){
		theta_ijConarri<-cbind(theta_ijConarri, theta_ijCon[,k])
	}else{
		theta_ijConabajo<-cbind(theta_ijConabajo, theta_ijCon[,k])
	}	
}


dfarriba <- data.frame(cbind(seq(1:n), theta_ijConarri))

dfabajo <- data.frame(cbind(seq(1:n), theta_ijConabajo))

darriba <- melt(dfarriba, id="X1")
dabajo <- melt(dfabajo, id="X1")


dOrig <- data.frame(cbind(seq(1:n), theta_ij))

dorig <- melt(dOrig , id="V1")

sp_sombras<-ggplot(dorig, aes(x = V1, y = value))+ geom_line()+ geom_line(data = darriba, aes(x=X1, y=value, color=variable), linetype = "dashed", alpha=0.1)+ geom_line(data = dabajo, aes(x=X1, y=value, color=variable), linetype = "dashed", alpha=0.4)+ labs(x = "Nodo", y = "Resultado")+ geom_point(aes(x=i, y=mean_tcond[1]), colour="black", shape=4)+geom_point(aes(x=j, y=mean_tcond[2]), colour="black", shape=4)+scale_colour_discrete(guide = 'none')


return(sp_sombras)

}



#graficas para datos perdidos
lineasCondColor<-function(theta_mtv, i,j, Covr, mean_t, nu, nlines){
	
	
n=dim(Covr)[1]
posi<-comprPosi(Covr)
QQ<-solve(round(Covr,posi))
theta_ij<-theta_mtv
theta_ij[i]<-NA
theta_ij[j]<-NA

mean_tcond<-(Covr[c(i,j), -c(i,j)]%*%ginv(Covr[-c(i,j), -c(i,j)])%*%(theta_mtv[-c(i,j)]-mean_t[-c(i,j)]))+mean_t[c(i,j)]

  beta<-t(theta_mtv[-c(i,j)]-mean_t[-c(i,j)])%*%QQ[-c(i,j), -c(i,j)]%*%(theta_mtv[-c(i,j)]-mean_t[-c(i,j)])
     
    CovCondi<-(as.numeric((nu+beta-2)/(nu+n-4)))*(ginv(QQ[c(i,j), c(i,j)])) 

nsim=nlines

theta_mtv_cond1<-matrix(nrow=nsim, ncol=2)

t=1
for(t in 1:nsim){
theta_mtv_cond<-rmvt(10000, mean_tcond, sig*round(CovCondi,2), n)
theta_mtv_cond1[t,]<-colMeans(theta_mtv_cond)
}

theta_ijCon<-matrix(ncol=nsim, nrow=n)

t=1

for(t in 1:nsim){
	theta_ijCon[,t]<-theta_ij
	theta_ijCon[i,t]<-theta_mtv_cond1[t,1]
	theta_ijCon[j,t]<-theta_mtv_cond1[t,2]
}


theta_ijConarri<-NULL
theta_ijConabajo<-NULL
k=1
for(k in 1:nsim){	if(theta_ijCon[i,k]>mean_tcond[1]){
		theta_ijConarri<-cbind(theta_ijConarri, theta_ijCon[,k])
	}else{
		theta_ijConabajo<-cbind(theta_ijConabajo, theta_ijCon[,k])
	}	
}


dfarriba <- data.frame(cbind(seq(1:n), theta_ijConarri))

dfabajo <- data.frame(cbind(seq(1:n), theta_ijConabajo))

darriba <- melt(dfarriba, id="X1")
dabajo <- melt(dfabajo, id="X1")


dOrig <- data.frame(cbind(seq(1:n), theta_ij))

dorig <- melt(dOrig , id="V1")

sp_color<-ggplot(dorig, aes(x = V1, y = value))+ geom_line()+ geom_line(data = darriba, aes(x=X1, y=value, group=variable),colour=c("#CC79A7"), alpha=0.5, linetype = "dashed")+geom_line(data = dabajo, aes(x=X1, y=value, group=variable), colour=c("#56B4E9"), linetype = "dashed")+ labs(x = "t", y = TeX("$\\theta_t$"))+ geom_point(aes(x=i, y=mean_tcond[1]), colour="black", shape=4)+ geom_point(aes(x=j, y=mean_tcond[2]), colour="black", shape=4)+ scale_colour_discrete(guide = 'none')


return(sp_color)

}



ComCarMTVG<-function(x, A, gam){	
CORR_mtv<-cov2cor(x)	
	
d<-rowSums(A)
D<-diag(d)
Cov_CAr<-solve(D-gam*A)

theta_CAR <-rmvnorm(9000, mean_t, sigma= 10*Cov_CAr)

CORR_CAR<-cor(theta_CAR)

CormtvS1orden<-CORR_mtv[which(A == 1, arr.ind = TRUE)]

CormtvS2orden<-CORR_mtv[which((A%^%2) == 1, arr.ind = TRUE)]

CorcarS1orden<-CORR_CAR[which(A == 1, arr.ind = TRUE)]

CorcarS2orden<-CORR_CAR[which((A%^%2) == 1, arr.ind = TRUE)]

CormtvS3orden<-CORR_mtv[which((A%^%3) == 1, arr.ind = TRUE)]

CormtvS4orden<-CORR_mtv[which((A%^%4) == 1, arr.ind = TRUE)]

CorcarS3orden<-CORR_CAR[which((A%^%3) == 1, arr.ind = TRUE)]

CorcarS4orden<-CORR_CAR[which((A%^%4) == 1, arr.ind = TRUE)]

bx1<-boxplot(CormtvS1orden, CorcarS1orden,CormtvS2orden, CorcarS2orden,CormtvS3orden, CorcarS3orden , CormtvS4orden, CorcarS4orden, ylab = "Correlation", xlab ="Neighbor order", 
at = c(1,2,3,4,5,6,7,8),
names = c("1","1", "2","2","3", "3", "4","4"),
las = 1,
col = c('powderblue', 'mistyrose', 'powderblue', 'mistyrose','powderblue', 'mistyrose', 'powderblue', 'mistyrose')
)
 return(bx1)
}




CAR<-function(Y,A,iters,burn){
	

 THETA<-NULL
 Y_pre<-NULL
   tick   <- proc.time()[3]
   # Bookkeeping
    n     <- length(Y)  
    m     <- rowSums(A)
    adj   <- apply(A==1,1,which)
    nei1  <- row(A)[A==1]
    nei2  <- col(A)[A==1]

   # Initial values
    sig2e <- var(Y)/2
    sig2s <- var(Y)/2
    rho   <- 0.99
    mu    <- mean(Y)
    theta <- 0.5*(Y-mu)

   # Pre-compute the determinant for all rho of interest

    C        <- diag(1/sqrt(m))%*%A%*%diag(1/sqrt(m))
    lambda   <- eigen(C)$values
    rho.grid <- seq(0.5,0.999,0.001)
    logd     <- rho.grid
    for(j in 1:length(logd)){
      logd[j] <- sum(log(1-rho.grid[j]*lambda))
    }   


   # Keep track of stuff

    keepers <- matrix(0,iters,5)
    colnames(keepers) <- c("sig2e","sig2s","mu","rho", "Ds")
    theta1 <- theta2 <- MSE1<- 0
 
   # GO!!!

   for(iter in 1:iters){
	cat("\n Iteration ", iter,"of ", iters)
      # THETA
       for(j in 1:n){
         neig     <- adj[[j]]
         ybar     <- mean(theta[neig])
         VVV      <- m[j]/sig2s + 1/sig2e
         MMM      <- rho*ybar*m[j]/sig2s + (Y[j]-mu)/sig2e
         theta[j] <- rnorm(1,MMM/VVV,1/sqrt(VVV))
       }
THETA<-cbind(THETA, theta)
      # VARIANCES

       TAT   <- sum(theta[nei1]*theta[nei2])
       TMT   <- sum(m*theta^2)
       sig2s <- 1/rgamma(1,n/2+.1,(TMT-rho*TAT)/2+.1) #sigmatheta
            sig2e <- 1/rgamma(1,n/2+.1,sum((Y-mu-theta)^2)+ .1)#sigmay

      # MEAN

       VVV  <- n/sig2e + 0.001
       MMM  <- sum(Y-theta)/sig2e
       mu   <- rnorm(1,MMM/VVV,1/sqrt(VVV))

      # CAR DEPENDENCE PARAMETER

       R    <- 0.5*logd + 0.5*rho.grid*TAT/sig2s
       rho  <- sample(rho.grid,1,prob=exp(R-max(R)))

      # KEEP TRACK OF STUFF

   
   Ds<- -2*sum(log(dmvnorm(Y, rep(mu, n)+theta, sig2e*diag(n))))
   
   MSE<-(Y-(rep(mu, n)+theta))^2
   
       keepers[iter,] <- c(sig2e,sig2s,mu,rho, Ds)
       
       if(iter>burn){
         theta1 <- theta1 + theta/(iters-burn)
         theta2 <- theta2 + theta*theta/(iters-burn)
           MSE1<-MSE1+MSE/(iters-burn)
       }

   }
   
   
   tock   <- proc.time()[3]
   output <- list(samps=keepers,
                  theta.mn=theta1, theta.var=theta2-theta1^2,
                  minutes=(tock-tick)/60, MSEES=MSE1, THETAES=THETA) 
 }

 
 
MTV<-function(Y,A, Ar, iters,burn){
	THETA<-NULL
 tick   <- proc.time()[3]	
 # Bookkeeping
mr<-rowSums(Ar)
Mr<-diag(rowSums(Ar))
	
l<-0.9
nu=n+2
CovR<-solve(Mr-l*Ar)
Covr<-C%*%CovR%*%t(C)
posi=1
posi<-comprPosi(Covr)
SgPrior<-round(Covr,posi)

   # Initial values
    sig2y <- var(Y)/2
    sig2t<-var(Y)/2
    gamma   <- 0.9
    mu    <- mean(Y)
    theta <- 0.5*(Y-mu)
    Sigma<-rinvwishart(nu, SgPrior)

 
   # Keep track of stuff

    keepers <- matrix(0,iters,5)
    colnames(keepers) <- c("sig2y","sig2t","mu","gamma", "Ds")
    
    theta1 <- theta2 <- MSE1<- 0

   # GO!!!
iter=1


b<-sig2t*diag(Sigma)*(nu-n-1)
phi<-mu

for(iter in 1:iters){
	cat("\n Iteration ", iter,"of ", iters)
      # THETA
      
       for(j in 1:n){
         MMM      <- (Y[j]-mu)/sig2y+phi/b[j]
         VVV<-1/sig2y+1/b[j]
         theta[j] <- rnorm(1,MMM/VVV,1/sqrt(VVV))
       }

###Sigma
K_pos<- SgPrior+(theta)%*%t(theta)/(sig2t*(nu-n-1))
nu_pos<-nu+n
Sigma<-rinvwishart(nu_pos, K_pos)
b<-sig2t*diag(Sigma)*(nu-n-1)


      # VARIANCE de y da mejor con varainza pequeña 

#sig2y <-var(Y)/2

sig2y <- 1/rgamma(1,n/2+.1,sum((Y-mu-theta)^2)+.1) #varianza de Y

sig2t <- 1/rgamma(1, 0.1+(n*nu)/2, tr(SgPrior%*%ginv(Sigma))/2+0.1)

      # MEAN de Y

VVV  <- n/sig2y + 0.001
MMM  <- sum(Y-theta)/sig2y
mu   <- rnorm(1,MMM/VVV,1/sqrt(VVV))

#phi<-mean(theta)

phi<-0
      # CAR DEPENDENCE PARAMETER

      # KEEP TRACK OF STUFF
      Ds<--2*sum(log(dmvnorm(Y, rep(mu, n)+theta, sig2y*diag(n))))
      
      MSE<-(Y-(rep(mu, n)+theta))^2

THETA<-cbind(THETA, theta)
      keepers[iter,] <- c(sig2y, sig2t, mu, gamma, Ds)
              
       if(iter>burn){
         theta1 <- theta1 + theta/(iters-burn)
         theta2 <- theta2 + theta*theta/(iters-burn)
         MSE1   <- MSE1 + MSE/(iters-burn)
       }
   }
      
tock   <- proc.time()[3]

output <- list(samps=keepers,
                  theta.mn=theta1, theta.var=theta2-theta1^2, 
                  minutes=(tock-tick)/60, MSEES=MSE1, THETAes=THETA) 
}




GraficaMicroMinasMSE<-function(Tabla, Variable, Resul){
		
Tabla<-cbind(Tabla, as.numeric(Resul))

colnames(Tabla)[length(Tabla)]<-c("Variable")

brasileiro <- merge(mun, Tabla, by.x = "ID", by.y = "CODE")

proj4string(brasileiro) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

brasileiro$Variable[is.na(brasileiro$Variable)] <- 0

name<-as.character(Variable)

namesEs<-as.character(mun$NM_MICRO)

Encoding(namesEs) <- "latin1"


pal <- colorBin("Blues",domain = NULL, n=10)

#pal <- colorBin("Blues", domain = 0:2)

popup <- paste0("<strong>Estado: </strong>", namesEs,  "<br><strong>Pontos: </strong>", brasileiro$Variable)

mapa<-leaflet(brasileiro) %>%
  # Opcion para anadir imagenes o mapas de fondo (tiles)
  addProviderTiles("Stamen.TonerHybrid") %>%
  # Funcion para agregar poligonos
  addPolygons(color = "#444444" ,
              weight = 1, 
              smoothFactor = 0.5,
              opacity = 1.0,
              fillOpacity = 0.5,
              fillColor = ~pal(brasileiro$Variable),    # Color de llenado
              layerId = ~brasileiro$Variable,                 
              highlightOptions = highlightOptions(color = "white", weight = 2,
                                                  bringToFront = TRUE), #highlight cuando pasas el cursor
              label = ~namesEs,                                  # etiqueta cuando pasas el cursor
              labelOptions = labelOptions(direction = "auto"),
              popup = popup) %>%                                        # mostrar el popup  
  addLegend(position = "topright", pal = pal, values = ~brasileiro$Variable,
            title = name) 
return(mapa)
}









GraficaMapaBrasil<-function(Tabla, Variable, Resul){
		
Tabla<-cbind(Tabla, as.numeric(Resul))

names(Tabla)[length(Tabla)]<-c("Variable")

brasileiro <- merge(shp, Tabla, by.x = "CD_GEOCUF", by.y = "CODIGO")

proj4string(brasileiro) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

brasileiro$NM_ESTADO<-as.character(brasileiro$NM_ESTADO)

Encoding(brasileiro$NM_ESTADO) <- "UTF-8"

brasileiro$Variable[is.na(brasileiro$Variable)] <- 0

pal <- colorBin("Blues",domain = NULL,n=5)

state_popup <- paste0("<strong>Estado: </strong>", brasileiro$NM_ESTADO,  "<br><strong>Pontos: </strong>", brasileiro$Variable)

name<-as.character(Variable)

mapa<-leaflet(data = brasileiro) %>% addProviderTiles("CartoDB.Positron") %>% addPolygons(fillColor = ~pal(brasileiro$Variable),  fillOpacity = 0.6,  color = "#BDBDC3", weight = 1,  popup = state_popup) %>% addLegend("bottomright", pal = pal, values = ~brasileiro$Variable, title = paste(name), opacity = 1)
return(mapa)
} #GraficaMapaBrasil(IDHMDAT, c("thetMTV"), theta1) Grfica 






GraficaMapaMinas<-function(Tabla, Variable, Resul){
	
Tabla<-cbind(Tabla, as.numeric(Resul))

names(Tabla)[length(Tabla)]<-c("Variable")

brasileiro <- merge(shp, Tabla, by.x = "CODMUN6", by.y = "code")

proj4string(brasileiro) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

brasileiro$SEM_ACENTO<-as.character(brasileiro$SEM_ACENTO)

Encoding(brasileiro$SEM_ACENTO) <- "UTF-8"

brasileiro$Variable[is.na(brasileiro$Variable)] <- 0

name<-as.character(Variable)

palBin <- colorBin("viridis", domain = brasileiro$Variable)

popup <- paste0("<strong>Estado: </strong>", brasileiro$SEM_ACENTO,  "<br><strong>Pontos: </strong>", brasileiro$Variable)

mapa<-leaflet(brasileiro) %>%
  # Opcion para anadir imagenes o mapas de fondo (tiles)
  addProviderTiles("Stamen.TonerHybrid") %>%
  # Funcion para agregar poligonos
  addPolygons(color = "#444444" ,
              weight = 1, 
              smoothFactor = 0.5,
              opacity = 1.0,
              fillOpacity = 0.5,
              fillColor = ~palBin(brasileiro$Variable),    # Color de llenado
              layerId = ~brasileiro$Variable,                 
              highlightOptions = highlightOptions(color = "white", weight = 2,
                                                  bringToFront = TRUE), #highlight cuando pasas el cursor
              label = ~brasileiro$SEM_ACENTO,                                  # etiqueta cuando pasas el cursor
              labelOptions = labelOptions(direction = "auto"),
              popup = popup) %>%                                        # mostrar el popup
  
  addLegend(position = "topright", pal = palBin, values = ~brasileiro$Variable,
            title = name) 
return(mapa)
} #Grafica minas (IDHMDAT, c("thetMTV"), theta1) Grficaa






GraficaMicroMinas<-function(Tabla, Variable, Resul){
		
Tabla<-cbind(Tabla, as.numeric(Resul))

colnames(Tabla)[length(Tabla)]<-c("Variable")

brasileiro <- merge(mun, Tabla, by.x = "ID", by.y = "CODE")

proj4string(brasileiro) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

brasileiro$Variable[is.na(brasileiro$Variable)] <- 0

name<-as.character(Variable)

namesEs<-as.character(mun$NM_MICRO)

Encoding(namesEs) <- "latin1"


#pal <- colorBin("Greens",domain = NULL, n=5)

pal <- colorBin("viridis", domain = brasileiro$Variable)

popup <- paste0("<strong>Estado: </strong>", namesEs,  "<br><strong>Pontos: </strong>", brasileiro$Variable)

mapa<-leaflet(brasileiro) %>%
  # Opcion para anadir imagenes o mapas de fondo (tiles)
  addProviderTiles("Stamen.TonerHybrid") %>%
  # Funcion para agregar poligonos
  addPolygons(color = "#444444" ,
              weight = 1, 
              smoothFactor = 0.5,
              opacity = 1.0,
              fillOpacity = 0.5,
              fillColor = ~pal(brasileiro$Variable),    # Color de llenado
              layerId = ~brasileiro$Variable,                 
              highlightOptions = highlightOptions(color = "white", weight = 2,
                                                  bringToFront = TRUE), #highlight cuando pasas el cursor
              label = ~namesEs,                                  # etiqueta cuando pasas el cursor
              labelOptions = labelOptions(direction = "auto"),
              popup = popup) %>%                                        # mostrar el popup  
  addLegend(position = "topright", pal = pal, values = ~brasileiro$Variable,
            title = name) 
return(mapa)
} #grafica para microregiones de Brasil





##MODEL1111

Model1<-function(Y,A, C, iters,burn){
	 tick   <- proc.time()[3]
set.seed(1)  
SIGMA<-DS<-RHO<-THETAEST<-RHOEST<-THETA<-alfa<-RHO_S<-NULL

Atri<-A 
Atri[lower.tri(A)] <- 0
veci<-which(Atri == 1, arr.ind = TRUE)

p=length(C[1,])
n=length(C[,1])
Ar<-t(C)%*%C-2*diag(rep(1,p))
Mr<-diag(rowSums(Ar))

 #es mejor un ICAR
 l=1
sigma<-1
Sigma_r<-sigma*ginv(Mr-l*Ar)
mean_t<-rep(0, n)
mean_r<-rep(0, p)
beta<-9000

RHO_sim<-rmvnorm(1, mean_r, Sigma_r) 

adj   <- apply(Ar==1,1,which)
nei1  <- row(Ar)[Ar==1]
nei2  <- col(Ar)[Ar==1]

mr<-rowSums(Ar)
sig2y <- var(Y)/2   
 tau_y<-1/sig2y
mu    <- mean(Y)
keepers <- matrix(0,iters,3)
colnames(keepers) <- c("sig2y","mu", "Ds")

t=1
ALFA<-NULL
for(t in 1:n){
ALFA[t]<- mean(rinvgamma(100, shape=1/2, scale = 1/9000)) # rinvgamma(n, shape, rate = 1, scale = 1/rate)
}

CrC<-C%*%diag(as.numeric(RHO_sim))%*%t(C)
CrCsD<-CrC-diag(diag(CrC))
SumVa<-diag(rowSums(abs(CrCsD)))

Q_pos<-(diag(ALFA)+ SumVa)-(CrCsD)

Spos<-ginv(tau_y*diag(n)+Q_pos)

###theta

mean_pos<-tau_y*Spos%*%(Y-mu)

theta_pos<-rmvnorm(1,mean_pos, Spos)


s=1
for(s in 1:iters){
cat("\n Iteration ", s,"of ", iters)
t=1
c<-NULL
for(t in 1:p){
prob<-(1+exp(-(abs(theta_pos[veci[t,1]]-theta_pos[veci[t,2]])-abs(theta_pos[veci[t,1]]+theta_pos[veci[t,2]]))/(2*sqrt(sigma))))^(-1)

c[t]<-rbern(1, prob)*(-1)
}

t=1
rho<-NULL

for(t in 1:p){
neig <- adj[[t]]

mean1<-2*(mean(RHO_sim[neig])+mean_r[t])

mean2<-(theta_pos[veci[t,1]]+c[t]*theta_pos[veci[t,2]])^2

rho[t]<-rnorm(1, (mean1+mean2)/(mr[t]), sigma/mr[t])
}


RHO_sim<-rho

t=1 
ALFA<-NULL
for(t in 1:n){
ALFA[t]<- rinv.gaussian(1, sqrt(beta*sigma)*abs(theta_pos[t])^(-1),1/beta)
}

CrC<-C%*%diag(rho)%*%t(C)
CrCsD<-CrC-diag(diag(CrC))
SumVa<-diag(rowSums(abs(CrCsD)))

Q_pos<-(diag(ALFA)+ SumVa)-(CrCsD)


      # MEAN varianza de Y

VVV  <- n/sig2y + 0.001
MMM  <- sum(Y-theta_pos)/sig2y
mu   <- rnorm(1,MMM/VVV,1/sqrt(VVV))
sig2y <- 1/rgamma(1,n/2+.1,sum((Y-mu-theta_pos)^2)+.1)

tau_y<-1/sig2y

Spos<-ginv(tau_y*diag(n)+Q_pos)

mean_pos<-tau_y*Spos%*%(Y-mu)

theta_pos<-rmvnorm(1, mean_pos, Spos)

Ds<--2*sum(log(dmvnorm(Y, rep(mu, n)+theta_pos, sig2y*diag(n))))

keepers[s,] <- c(sig2y, mu, Ds)

RHO<-cbind(RHO, rho)

THETA<-cbind(THETA, as.numeric(theta_pos))
alfa<-cbind(alfa, as.numeric(ALFA))
RHO_S<-cbind(RHO_S, as.numeric(RHO_sim))
} 
RHOEST<-rowMeans(RHO[,burn:iters])
THETAEST<-rowMeans(THETA[,burn:iters])
tock   <- proc.time()[3]
 output <- list(samps=keepers,
                  Rho=RHOEST, thetamn=THETAEST, minutes=(tock-tick)/60, alfa=alfa, RHO_S=RHO_S)
                  
                  }
                  
                  
                  
 #######model1-2
 
 

Model1c<-function(Y,A, C, iters,burn){
	 tick   <- proc.time()[3]
set.seed(1)  
SIGMA<-DS<-RHO<-THETAEST<-RHOEST<-THETA<-alfa<-RHO_S<-NULL

Atri<-A 
Atri[lower.tri(A)] <- 0
veci<-which(Atri == 1, arr.ind = TRUE)

p=length(C[1,])
n=length(C[,1])
Ar<-t(C)%*%C-2*diag(rep(1,p))
Mr<-diag(rowSums(Ar))

 #es mejor un ICAR
 l=1
sigma<-1
Sigma_r<-sigma*ginv(Mr-l*Ar)
mean_t<-rep(0, n)
mean_r<-rep(0, p)
beta<-9000

RHO_sim<-rmvnorm(1, mean_r, Sigma_r) 

adj   <- apply(Ar==1,1,which)
nei1  <- row(Ar)[Ar==1]
nei2  <- col(Ar)[Ar==1]

mr<-rowSums(Ar)
sig2y <- var(Y)/2   
 tau_y<-1/sig2y
mu    <- mean(Y)
keepers <- matrix(0,iters,3)
colnames(keepers) <- c("sig2y","mu", "Ds")

t=1
ALFA<-NULL
for(t in 1:n){
ALFA[t]<- mean(rinvgamma(100, shape=1/2, scale = 1/9000)) # rinvgamma(n, shape, rate = 1, scale = 1/rate)
}

CrC<-C%*%diag(as.numeric(RHO_sim))%*%t(C)
CrCsD<-CrC-diag(diag(CrC))
SumVa<-diag(rowSums(abs(CrCsD)))

Q_pos<-(diag(ALFA)+ SumVa)-(CrCsD)

Spos<-ginv(tau_y*diag(n)+Q_pos)

###theta

mean_pos<-tau_y*Spos%*%(Y-mu)

theta_pos<-rmvnorm(1,mean_pos, Spos)


s=1
for(s in 1:iters){
cat("\n Iteration ", s,"of ", iters)
t=1
c<-NULL
for(t in 1:p){
prob<-(1+exp(-(abs(theta_pos[veci[t,1]]-theta_pos[veci[t,2]])-abs(theta_pos[veci[t,1]]+theta_pos[veci[t,2]]))/(2*sqrt(sigma))))^(-1)


c[t]<-ifelse(rbern(1, prob)==1, 1, -1)
}

t=1
rho<-NULL

for(t in 1:p){
neig <- adj[[t]]

mean1<-2*(mean(RHO_sim[neig])+mean_r[t])

mean2<-(theta_pos[veci[t,1]]+c[t]*theta_pos[veci[t,2]])^2

rho[t]<-rnorm(1, (mean1+mean2)/(mr[t]), sigma/mr[t])

RHO_sim[t]<-rho[t]
rho[t]<-c[t]*rho[t]
}



t=1 
ALFA<-NULL
for(t in 1:n){
ALFA[t]<- rinv.gaussian(1, sqrt(beta*sigma)*abs(theta_pos[t])^(-1),1/beta)
}

CrC<-C%*%diag(rho)%*%t(C)
CrCsD<-CrC-diag(diag(CrC))
SumVa<-diag(rowSums(abs(CrCsD)))

Q_pos<-(diag(ALFA)+ SumVa)-(CrCsD)


      # MEAN varianza de Y

VVV  <- n/sig2y + 0.001
MMM  <- sum(Y-theta_pos)/sig2y
mu   <- rnorm(1,MMM/VVV,1/sqrt(VVV))
sig2y <- 1/rgamma(1,n/2+.1,sum((Y-mu-theta_pos)^2)+.1)

tau_y<-1/sig2y

Spos<-ginv(tau_y*diag(n)+Q_pos)

mean_pos<-tau_y*Spos%*%(Y-mu)

theta_pos<-rmvnorm(1, mean_pos, Spos)

Ds<--2*sum(log(dmvnorm(Y, rep(mu, n)+theta_pos, sig2y*diag(n))))

keepers[s,] <- c(sig2y, mu, Ds)

RHO<-cbind(RHO, rho)

THETA<-cbind(THETA, as.numeric(theta_pos))
alfa<-cbind(alfa, as.numeric(ALFA))
RHO_S<-cbind(RHO_S, as.numeric(RHO_sim))
} 
RHOEST<-rowMeans(RHO[,burn:iters])
THETAEST<-rowMeans(THETA[,burn:iters])
tock   <- proc.time()[3]
 output <- list(samps=keepers,
                  Rho=RHOEST, thetamn=THETAEST, minutes=(tock-tick)/60, alfa=alfa, RHO_S=RHO_S)
                  
                  }                 