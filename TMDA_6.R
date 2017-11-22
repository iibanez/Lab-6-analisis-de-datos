
require(e1071)
require(doParallel)

#root <- "C:\\Users\\Usuario\\Documents\\1. Universidad\\Nivel 10\\Topico II - Mineria de Datos Avanzados\\Lab6TMDA\\"
#datos <- read.csv(paste0(root,"\\G4_001.csv"))

#root <- "/Volumes/HDD/Google Drive/2017 - 2/Miner??a de Datos Avanzada/V.- Support Vector Regression y Recurrencia Externa/"
datos <- read.csv("G5_002.csv")



retrasos<- function(data, lags){
  largo <- length(data$PAM)
  retrasos <- data.frame(matrix(ncol=0, nrow=largo-lags))
  
  for(i in 1:lags){
    PAM <- data$PAM[i:(largo-lags+i-1)]
    retrasos<-data.frame(retrasos,PAM)   
  }
  retrasos[['VFSC']] <- data$VFSC[(lags+1):largo]
  return(retrasos)   
}

#Paralelizaci?n 
Ts=0.2
registerDoParallel(cores = 8)
cost <- 2^seq(-5, 15, 1)
nu <- seq(0.1, 0.9,0.1)
gamma<-2^seq(-15, 3, 1)

#Normalizaci?n de datos
PAM<-(datos$PAM-min(datos$PAM))/(max(datos$PAM)-min(datos$PAM))
VFSC<-(datos$VFSC-min(datos$VFSC))/(max(datos$VFSC)-min(datos$VFSC))
data.normalized <- data.frame(PAM,VFSC)

#Generaci?n de conjuntos de test
train.index <- nrow(data.normalized)/2.0
data.train <- data.normalized[1:train.index,]
data.test <- data.normalized[train.index:(nrow(data.normalized)-1),]


for(l in 3:8){

retDatos.train <- retrasos(data.train, l)
x.train=subset(retDatos.train, select = -VFSC)
y.train=retDatos.train$VFSC
  
retDatos.test <- retrasos(data.test, l)
x.test=subset(retDatos.test, select = -VFSC)
y.test=retDatos.test$VFSC


parms <- expand.grid(cost = cost, nu = nu, gamma=gamma)

salida <- (c( foreach(i = 1:nrow(parms), combine = rbind, .inorder = FALSE)
              %dopar% {
                c <- parms[i, ]$cost
                n <- parms[i, ]$nu
                g <- parms[i, ]$gamma
                modelo <- e1071::svm(x.train, y.train, type = "nu-regression", kernel = "radial", cost =
                                       c, nu = n, gamma=g)
                pred <- predict(modelo, x.test)
                corr_pred<-cor(pred,y.test,method = "pearson")
                c(l,c, n, g, corr_pred)
              }
              ))

output <- matrix(unlist(salida), ncol = 5, byrow = TRUE)
mejoresModelos<-output[order(output[,5], decreasing = TRUE),]
colnames(mejoresModelos) <- c("Lags","Costos","NU","gamma","Correlaci??n")
write.csv(mejoresModelos, file = paste0("Mejores Modelos L=",l," AB.csv"),row.names = FALSE,sep = ";",dec =",")

salida <- (c( foreach(i = 1:nrow(parms), combine = rbind, .inorder = FALSE)
              %dopar% {
                c <- parms[i, ]$cost
                n <- parms[i, ]$nu
                g <- parms[i, ]$gamma
                modelo <- e1071::svm(x.test, y.test, type = "nu-regression", kernel = "radial", cost =
                                       c, nu = n, gamma=g)
                pred <- predict(modelo, x.train)
                corr_pred<-cor(pred,y.train,method = "pearson")
                c(l,c, n, g, corr_pred)
              }
))

output <- matrix(unlist(salida), ncol = 5, byrow = TRUE)
mejoresModelos<-output[order(output[,5], decreasing = TRUE),]
colnames(mejoresModelos) <- c("Lags","Costos","NU","gamma","Correlaci??n")
write.csv(mejoresModelos, file = paste0("Mejores Modelos L=",l," BA.csv"),row.names = FALSE,sep = ";",dec =",")

}

mejoresModelos <- read.csv("Mejores Modelos L=6 AB_002.csv")[1:2,]

#output <- matrix(unlist(salida), ncol = 5, byrow = TRUE)
#mejoresModelos<-output[order(output[,5], decreasing = TRUE),]

inverseStep=matrix(1,180/Ts,1)
inverseStep[(90/Ts):(180/Ts),1]=0

for (i in 1:length(mejoresModelos[,1])){
  lag<-list(PAMn = mejoresModelos[i,1],VFSCn = 0)
  
  retDatos.train <- retrasos(data.train, lag$PAMn)
  x.train=subset(retDatos.train, select = -VFSC)
  y.train=retDatos.train$VFSC
  
  retDatos.test <- retrasos(data.test, lag$PAMn)
  x.test=subset(retDatos.test, select = -VFSC)
  y.test=retDatos.test$VFSC
  
  mejorModelo <- svm(x.train, y.train, kernel = "radial",type = "nu-regression", cost = mejoresModelos[i,2], nu = mejoresModelos[i,3], gamma=mejoresModelos[i,4])
  PAMn=inverseStep
  VFSCn=inverseStep
  data <- data.frame(PAMn,VFSCn)
  
  lag<-list(PAMn = mejoresModelos[i,1],VFSCn = 0)
  retDatos <- retrasos(data, lag$PAMn) 
  x=subset(retDatos, select = -VFSCn) 
  y=retDatos$VFSCn 
  stepTime=seq(Ts,(length(retDatos$PAM))*Ts,Ts)
  stepResponse <- predict(mejorModelo, x )
  plot(stepTime,retDatos$PAM,type="l", col="red") 
  lines(stepTime,stepResponse, col = "blue")
  legend("topright", c("Escalon de presi??n", "respuesta al escalon"), title = "autorregulacion", pch = 1, col=c("red","blue"),lty=c(1,1),inset = 0.01)
  print(paste("corr=",mejoresModelos[5]))
  readline(prompt="Press [enter] to continue")
}
