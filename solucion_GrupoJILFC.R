# solucion trivial ejemplo
remove(list=ls())
setwd("/home/emiliano/Documents/estadistica/estadistica_y_optimizacion_master/proyecto/")
dir()
# cargamos funciones con soluciones triviales. 
# Cada equipo debera definir en su archivo de solucion_teamName.R todas las
# funciones incluidas abajo

#############################################################################################

teamName <- "GrupoJILFC"
# integrante 1: Aguilar Marin, Jeremy Joel
# integrante 2: Barba La Orden, Irene
# integrante 3: Benages Guijarro, Lucía
# integrante 4: Calvo Castillo, Fabián
# integrante 5: Torrico Castellon, Cristhian

##########################
# librerías utilizadas
##########################
library(forecast)  
library(FinTS)  
library(quadprog)

##########################
# seccion 2 - predicciones
##########################

# Funcion para predecir UN activo para el periodo test en el momento t
getPred <- function(x_train, x_test_past){
  # INPUTS: x_train serie de tiempo de todo periodo train para UN activo
  #       : x_test_past serie de tiempo de periodo test que va desde 0 hast t-1
  
  # OUTPUTS: mu_hat: prediccion del activo para el momento t del periodo test
  #          se_hat: desviación estándar de la predicción. 
  
  mod <- auto.arima(x_train, stepwise = FALSE, approximation = FALSE)
  x_past = c(x_train, x_test_past)
  fit_upd <- Arima(x_past, model = mod)
  fc <- forecast::forecast(fit_upd, h = 1)
  mu_hat  <- as.numeric(fc$mean[1])
  
  z80 <- qnorm(0.8)
  se_from_up <- (fc$upper[,"80%"][1] - fc$mean[1]) / z80
  se_from_lo <- (fc$mean[1]  - fc$lower[,"80%"][1]) / z80
  se_hat_aux <- pmax(se_from_up, se_from_lo) 
  se_hat  <- as.numeric(se_hat_aux)
  
    return(list(mu_hat=mu_hat, se_hat=se_hat))
}

#####################################
# seccion 3 - utilidad media-varianza
###########################################################
# 3.1 utilidad media-varianza, alfa_i positiva o negativa
###########################################################
gammaMV <- 20

# Funcion para estimar la matriz de covarianzas entre los rendimientos de los 5
# activos a partir de las desviaciones estándares (que vendran de su modelo Arima)
# y el historico de rendimientos (que pueden usar para estimar correlaciones) 
getSigmaMV <- function(sig, Xpast){
  # INPUT: sig, vector en R^5 con desviaciones estandar de 5 activos al momento t
  #        Xpast, matriz en R^(T x 5) con rendimientos de 5 activos desde t=0 hasta t-1
  # OUTPUT: Sigma: matriz en 5 x 5 con covarianzas entre activos
  
  R <- cor(Xpast, use="pairwise.complete.obs") 
  D <- diag(sig)  
  Sigma <- D %*% R %*% D
  return(Sigma)
}

# Funcion para optimizar la asignacion de cada activo dentro del portafolio
# tomando en cuenta la utilidad Media-Varianza CON posiciones cortas
getAlphaMV <- function(mu,Sigma, gamma){
  # INPUT: mu: vector en R^5 con los rendimientos esperados de los 5 activos para el periodo t
  #        Sigma: matriz en R^{5 x 5} con las covarianzas de los 5 activos
  #        gamma: valor en R que representa el apetito de riesgo. Entre mas bajo mas riesgo.
  # OUTPUT: alpha: vector en R^5 con la asignación elegida para los 5 activos resultante de optimizar
  # U-MV sujeto a que alpha sume a 1.
  
  vec <- rep(1, length(mu))
  invSigma <- solve(Sigma)
  
  # Cálculo de Lambda y Alpha óptimo
  lambda <- (t(vec) %*% invSigma %*% mu - gamma) / (t(vec) %*% invSigma %*% vec)
  alpha <- (1/gamma) * invSigma %*% (mu - as.numeric(lambda) * vec)
  
  return(as.numeric(alpha))
} 

############################################
# 3.2 utilidad media-varianza, alfa_i positiva 
############################################

gammaMVPos <- 20

# Funcion para estimar la matriz de covarianzas entre los rendimientos de los 5
# activos a partir de las desviaciones estándares (que vendran de su modelo Arima)
# y el historico de rendimientos (que pueden usar para estimar correlaciones) 
getSigmaMVPos <- function(sig, Xpast){
  # INPUT: sig, vector en R^5 con desviaciones estandar de 5 activos al momento t
  #        Xpast, matriz en R^(T x 5) con rendimientos de 5 activos desde t=0 hasta t-1
  # OUTPUT: Sigma: matriz en 5 x 5 con covarianzas entre activos
  
  # Para el modelo MV Positivo (usualmente es la misma lógica de riesgo)
  return(getSigmaMV(sig, Xpast))
}

# Funcion para optimizar la asignacion de cada activo dentro del portafolio
# tomando en cuenta la utilidad Media-Varianza SIN posiciones cortas
getAlphaMVPos <- function(mu,Sigma, gamma){
  # INPUT: mu: vector en R^5 con los rendimientos esperados de los 5 activos para el periodo t
  #        Sigma: matriz en R^{5 x 5} con las covarianzas de los 5 activos
  #        gamma: valor en R que representa el apetito de riesgo. Entre mas bajo mas riesgo.
  # OUTPUT: alpha: vector en R^5 con la asignación elegida para los 5 activos resultante de optimizar
  # U-MV sujeto a que alpha sume a 1 y alpha_i>0.
  
  n <- length(mu)
  
  # Configuración para solve.QP
  Dmat <- gamma * Sigma + diag(1e-8, n) # Ajuste numérico
  dvec <- mu
  Amat <- cbind(rep(1, n), diag(n))     # Restricciones: Suma=1, I>=0
  bvec <- c(1, rep(0, n))               # Valores: 1, 0, 0, 0...
  
  # Resolución numérica
  res <- tryCatch({
    sol<- solve.QP(Dmat, dvec, Amat, bvec, meq=1)
    sol$solution
  }, error = function(e) {
    return(rep(1/n, n)) # Si falla, pesos iguales
  })
  res[res < 0] <- 0 
  res <- res / sum(res)
  
  return(res)
}


################################################
# seccion 4 - 
# utilidad log, alfa_i positiva o negativa
################################################
gammaLog <- 20

# Funcion para estimar la matriz de covarianzas entre los rendimientos de los 5
# activos a partir de las desviaciones estándares (que vendran de su modelo Arima)
# y el historico de rendimientos (que pueden usar para estimar correlaciones) 
getSigmaLog <- function(sig, Xpast){
  # INPUT: sig, vector en R^5 con desviaciones estandar de 5 activos al momento t
  #        Xpast, matriz en R^(T x 5) con rendimientos de 5 activos desde t=0 hasta t-1
  # OUTPUT: Sigma: matriz en 5 x 5 con covarianzas entre activos
  
  Sigma <- getSigmaMV(sig, Xpast)
  return(Sigma)
}

# Funcion para optimizar la asignacion de cada activo dentro del portafolio
# tomando en cuenta la utilidad log CON posiciones cortas
getAlphaLog <- function(mu,Sigma, gamma){
  # INPUT: mu: vector en R^5 con los rendimientos esperados de los 5 activos para el periodo t
  #        Sigma: matriz en R^{5 x 5} con las covarianzas de los 5 activos
  #        gamma: valor en R que representa el apetito de riesgo. Entre mas bajo mas riesgo.
  # OUTPUT: alpha: vector en R^5 con la asignación elegida para los 5 activos resultante de optimizar
  # U-log sujeto a que alpha sume a 1.
  
  # Repetiremos la optimización hasta verificar que el alpha suma 1
  comprobar_alpha <- TRUE
  # Penalización para la restricción
  penal <- 5
  while(comprobar_alpha){
    # Definición de la función a maximizar, con un parámetro de penalización
    Ulog_penal <- function(alpha){
      ulog <- Ulog(alpha, gamma, mu, Sigma)
      penal <- penal*(sum(alpha)-1)^2
      
      return(as.numeric(ulog - penal))
    }
    # Partimos de un alpha cualquiera y hacemos primero una optimización con el método Nelder-Mead para encontrar 
    # un punto de inicio para optimizar
    # Podemos obtener a veces warnings si el logaritmo tiene argumento negativo y genera NaN
    alpha_0 <- rep(1/5,5)
    res_nm <- suppressWarnings(optim(par=alpha_0, fn =function(alpha) -Ulog_penal(alpha), method = "Nelder-Mead"))
    alpha_nm <- res_nm$par
    # Encontramos un mejor óptimo con el algoritmo BFGS
    # Pero no es alarmante ya que optim simplemente considera este punto inviable y sigue iterando
    res <- suppressWarnings(optim(par=alpha_nm, fn =function(alpha) -Ulog_penal(alpha), method = "BFGS"))
    summary(res)
    alpha <- res$par
    
    # Comprobamos que suma 1 y si no aumentamos la penalización
    if(sum(alpha)-1 < 1e-6){
      # Si suma 1 se acaba aquí el bucle y devolveremos este alpha
      comprobar_alpha <- FALSE
    } else{
      penal <- 10*penal
    }
  }
  return(alpha)
}


# ... HASTA AQUI
#############################################################################################


###############################################################
# Evaluación de soluciones
###############################################################
source("funciones/eval_funcs.R")
#source("eval_funcss.R")

setwd("/home/emiliano/Documents/estadistica/estadistica_y_optimizacion_master/proyecto/")
X <- read.csv("stock_returns_train_2.csv")
X <- ts(X)/100

# Validation mode - para que se evaluen asi mismos con el 
Xtrain <- window(X, start=1,end=8*12) # el start-end es un ejemplo, pueden cambiarlo
Xtest <- window(X, start=8*12+1,end=10*12)

# Test mode - no tendran el archivo stock_returns_test.csv asi que esto lo 
# ejecutaremos una vez entreguen soluciones
#Xtrain <- X
#Xtest <- ts(read.csv("stock_returns_test.csv"))/100


#seccion 2 - predicciones
set.seed(43)
res <- getPred_ts(Xtrain, Xtest, getPred)
mu_hat = res$mu_hat
se_hat = res$se_hat

# MAPE
i <- 2
plot(as.data.frame(Xtest)[,i], ty="l")
lines(mu_hat[,i], col="blue", ty="l")



rmse <- sqrt(mean((Xtest-mu_hat)^2))
evals <- c(rmse=rmse)
evals

# seccion 3 - utilidad media varianza
# utilidad media-varianza, alfa_i positiva o negativa

alpha_hat <- getAlpha_ts(mu_hat, se_hat, gammaMV, getSigmaMV, getAlphaMV, Xtrain, Xtest)
passChecks <- getChecks(alpha_hat, mode="sum1")
ret <- getRet(alpha_hat, Xtest, passChecks)
evals <- c(evals, retMV=ret)
Umv_rel <- getUEval(alpha_hat, mu_hat, se_hat, Xtrain, Xtest, gammaMV, getSigmaMV, passChecks, Umv)
evals <- c(evals,  Umv=Umv_rel)

# utilidad media-varianza, alfa_i positiva 

alpha_hat <- getAlpha_ts(mu_hat, se_hat, gammaMVPos, getSigmaMVPos, getAlphaMVPos, Xtrain, Xtest)
passChecks <- getChecks(alpha_hat, mode=c("sum1","pos"))
ret <- getRet(alpha_hat, Xtest, passChecks)
evals <- c(evals, retMVPos=ret)
Umv_rel <- getUEval(alpha_hat, mu_hat, se_hat, Xtrain, Xtest, gammaMVPos, getSigmaMVPos, passChecks, Umv)
evals <- c(evals,  UmvPos=Umv_rel)


# seccion 4 - 
# utilidad log, alfa_i positiva o negativa

alpha_hat <- getAlpha_ts(mu_hat, se_hat, gammaLog, getSigmaLog, getAlphaLog, Xtrain, Xtest)
passChecks <- getChecks(alpha_hat, mode=c("sum1"))
ret <- getRet(alpha_hat, Xtest, passChecks)
evals <- c(evals, retLog=ret)
Umv_rel <- getUEval(alpha_hat, mu_hat, se_hat, Xtrain, Xtest, gammaLog, getSigmaLog, passChecks, Umv)
evals <- c(evals,  Ulog=Umv_rel)

evals


