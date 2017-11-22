
require(doParallel) 
require(e1071) 
registerDoParallel(cores = 8) 
cost <- 2^seq(15, 15, 1)
nu <- seq(0.1, 0.1, 0.1)
gamma<-2^seq(1, 1, 1)
lagsList<-seq(6)
datos=read.csv("G5_001.csv") 
PAMn<-(datos$PAM-min(datos$PAM))/(max(datos$PAM)-min(datos$PAM)) 
VFSCn<-(datos$VFSC-min(datos$VFSC))/(max(datos$VFSC)-min(datos$VFSC)) 
data <- data.frame(PAMn,VFSCn)
Ts=0.2
parms <- expand.grid(lagsList=lagsList, cost = cost, nu = nu, gamma=gamma)

retardos_multi <- function( signalData, lags){ 
  signal.uni <- signalData
  max.lag <- max(unlist(lags)) + 1
  indices <- 1:nrow(signal.uni)
  lag.mat <- embed(indices, max.lag) 
  col.names <- list("PAMn","VFSCn") 
  columns <- NULL 
  lagged.columns.names <- c() 
  for(colname in col.names){
    lag.order <- lags[[colname]]
    columns[[colname]] <- signal.uni[lag.mat[, 1], colname] 
    if(!is.null(lag.order) && lag.order > 0)
      for(i in 1:lag.order){
        new.colname <- paste(colname, paste0("lag", i), sep = ".") 
        lagged.columns.names <- c(lagged.columns.names, new.colname) 
        columns[[new.colname]] <- signal.uni[lag.mat[, i+1], colname]
      }
  }
  folded.signal <- data.frame(columns)
  sorting <- order(lag.mat[, 1])
  folded.signal <- folded.signal[sorting, ]
  list(folded.signal = folded.signal, lagged.columns.names = lagged.columns.names)
}

salida <- (c( foreach(i = 1:nrow(parms), combine = rbind, .inorder = FALSE) %dopar% {
  c <- parms[i, ]$cost
  n <- parms[i, ]$nu
  g <- parms[i, ]$gamma
  l <- parms[i, ]$lagsList
l = 5
lag<-list(PAMn = l,VFSCn = 0)
signal.train <- retardos_multi(data, lag)
retDatos=signal.train$folded.signal
  x=subset(retDatos, select = -VFSCn)
  y=retDatos$VFSCn
  modelo <- e1071::svm(x, y, type = "nu-regression", kernel = "radial", cost = c, nu = n, gamma=g)
  pred <- predict(modelo, x) 
  corr_pred <- cor(pred,y,method = "pearson") 
  c(l, c, n, g, corr_pred)
}))

output <- matrix(unlist(salida), ncol = 5, by = TRUE)
mejoresModelos<-output[order(output[,5], decreasing = TRUE),]

inverseStep=matrix(1,180/Ts,1)
inverseStep[(90/Ts):(180/Ts),1]=0

for (i in 1:length(mejoresModelos[,1])){
i = 1
PAMn<-(datos$PAM-min(datos$PAM))/(max(datos$PAM)-min(datos$PAM)) 
VFSCn<-(datos$VFSC-min(datos$VFSC))/(max(datos$VFSC)-min(datos$VFSC)) 
data <- data.frame(PAMn,VFSCn)
lag<-list(PAMn = mejoresModelos[1],VFSCn = 0)
signal.train <- retardos_multi(data, lag)
retDatos=signal.train$folded.signal
# signal.train
# $folded.signal
# PAMn   PAMn.lag1        VFSCn
# 1    0.369455645 0.367439516 0.5163566389
# 2    0.368951613 0.369455645 0.5003207184
x=subset(retDatos, select = -VFSCn)
y=retDatos$VFSCn
mejorModelo <- svm(x, y, kernel = "radial",type = "nu-regression", cost = mejoresModelos[2], nu = mejoresModelos[3], gamma=mejoresModelos[4])
PAMn=inverseStep
VFSCn=inverseStep
data <- data.frame(PAMn,VFSCn)
lag<-list(PAMn = mejoresModelos[1],VFSCn = 0)
signal.train <- retardos_multi(data, lag) 
retDatos=signal.train$folded.signal
# PAMn   PAMn.lag1   PAMn.lag2   PAMn.lag3   PAMn.lag4   PAMn.lag5   PAMn.lag6        VFSCn
# 1    0.356350806 0.359879032 0.363911290 0.366935484 0.368951613 0.369455645 0.367439516 0.4823604875
# 2    0.352318548 0.356350806 0.359879032 0.363911290 0.366935484 0.368951613 0.369455645 0.4881334189
# 3    0.349798387 0.352318548 0.356350806 0.359879032 0.363911290 0.366935484 0.368951613 0.4964720975
# 4    0.348286290 0.349798387 0.352318548 0.356350806 0.359879032 0.363911290 0.366935484 0.5067350866
# 5    0.348286290 0.348286290 0.349798387 0.352318548 0.356350806 0.359879032 0.363911290 0.5195638230
# 6    0.353326613 0.348286290 0.348286290 0.349798387 0.352318548 0.356350806 0.359879032 0.5336754330
# 7    0.363407258 0.353326613 0.348286290 0.348286290 0.349798387 0.352318548 0.356350806 0.5484284798
x=subset(retDatos, select = -VFSCn) 
y=retDatos$VFSCn 
stepTime=seq(Ts,(length(retDatos$PAMn))*Ts,Ts)
stepResponse <- predict(mejorModelo, x )
plot(stepTime,retDatos$PAMn,type="l", col="red") 
lines(stepTime,stepResponse, col = "blue")
legend("topright", c("Escalon de presi??n", "respuesta al escalon"), title = "autorregulacion", pch = 1, col=c("red","blue"),lty=c(1,1),inset = 0.01)
print(paste("corr=",mejoresModelos[5]))
readline(prompt="Press [enter] to continue")
}
