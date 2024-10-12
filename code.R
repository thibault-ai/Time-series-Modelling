#######################################################
### Definition du repertoire- importation base
###########################################################
# importation des librairies
require(zoo) #format de serie temporelle pratique et facile d'utilisation (mais plus volumineux)
require(tseries) #diverses fonctions sur les series temporelles


## importation de la base 
data <- read.csv("valeurs_mensuelles.csv",sep=";") # base de donnees originale
df <- data[-c(1:2), ]                             # on enleve les deux valeurs les plus recentes
MOIS <- as.yearmon(seq(from=2023+11/12, to=1990+0/12, by=-1/12)) # formatage des dates
spread <- zoo(df$spread, order.by=MOIS)

################################################################################
#  PARTIE 1 : LES DONNEES 
################################################################################

## 1- Representation de la serie 
plot(spread, xlab = "Année", ylab = "Indice")

## 1-b Justification de la non stationnarite

# Etape 1 : On fait une regression de spread sur dates.! 

dates <- as.yearmon(seq(from=1990+0/12, to =  2023+11/12 , by=1/12))
# on trouve : tendance et intercept apparemment non nuls ( regard des p_valeurs)
summary(lm(spread ~ dates))

# Implication: usage du test de Dickey fuller avec tendance et constante. 

require(fUnitRoots) 
adf <- adfTest(spread, lag=0, type="ct") #
#  test de racine unitaire est effectué sans retard; ct signifie constante et trend inclue. 
adf

# conclusion : la série n'est pas stationnaire 

# Le test n'est valide que si les residus sont des bruits blancs faible. 
Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}


Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))


# On rajoute des lags jusqu'à obtenirdes résidus décorrelés
# si au bout de 24 on a que des NA, c'est qu'on a que des garbages
series <- spread; kmax <- 24; adftype="ct"
adfTest_valid <- function(series, kmax, adftype){
  k <- 0    # lag
  noautocorr <- 0    ## est ce que les residus sont des bruits blancs ?
  while (noautocorr==0){
    # Un message d'avancelent pour l'utilisateur
    cat(paste0("ADF with ",k," lags: residuals OK? "))
    # on effectue le test ( en tout cas les calculs)
    # rappel : garbage in & garbage out
    # tant que le lag ne suffit pas à décorreler les résidus
    # vous recolterez du garbage en sortie. 
    adf <- adfTest(series, lags=k, type=adftype)
    pvals <- Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))[,2]
    if (sum(pvals<0.05,na.rm=T)==0) {
      noautocorr <- 1; cat("OK \n")
    } else cat("nope \n")
    k <- k+1
  }
  return(adf)
}

adf <- adfTest_valid(spread,24,adftype="ct")
# H0 : racine / H1: pas de racine. 
#
Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))

# H0 : racine / H1: pas de racine. 
# Les résidus à l'ordre 21 sont décorrelés ; La p_value du test peut etre interpreter aisement

adf
# conclusion : p_value = 0.70 ** on ne rejette pas H0 : presence racine unitaire

### Faire le test KPSS

require(urca)
kpss.test(spread) ### 

# CONCLUSION : spread est non stationnaire / il faut la stationnariser

## 2- STATIONNARISATION DE LA SERIE 

##  A .Difference premiere de "spread" et stocke le résultat dans la variable "dspread".
dspread <- diff(spread,1)

##  B. La serie est stationnaire : Test lyung box

# Application de la fonction adfTest_valid ;  
# à l'ordre 0 , les résidus semblent décorrelés. 
adf <- adfTest_valid(dspread,24,"nc")

## H0 : racine unitaire / les résultats du test 
#  nous indique que le modèle ne contient pas de racine unitaire
Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))
adf

## 3-  REPRESENTATION DE LA SERIE AVANT ET APRES LA DIFFERENCIATION

par(mfrow=c(1,2))
plot(spread, xlab="(a) Série brute", ylab="Indice");
plot(dspread, xlab=" (b) Série différenciée ",ylab="dindice")


############################################################
# PARTIE 2- LE MODELE ARMA
############################################################

### 4 - CHOIX SELECTION ET ESTIMATION DES PARAMETRES 

## ETAPE 1 - VERIFICATION DE LA STATIONNARITE  

# Constante et tendance présente. 
summary(lm(dspread ~ dates[-1]))

# constante pas significative. 

##
# Application de la fonction adfTest_valid ;  
# à l'ordre 0 , les résidus semblent décorrelés. 
adf <- adfTest_valid(dspread,24,"nc")

## H0 : racine unitaire / les résultats du test 
#  nous indique que le modèle ne contient pas de racine unitaire
Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))
adf

# KPSS test 
kpss.test(dspread)

# test de philipp peron
pp.test(dspread)

## CONCLUSION : la serie dspread est stationnaire 

## ETAPE 2 - CHOIX DU MODELE 

##  recherche des ordres p et q
par(mfrow=c(1,2))
acf(dspread,24, main="", xlab="Retard", ylab="Autocorrelation")
pacf(dspread,24,main="", xlab="Retard",ylab="Autocorrelation partielle");     
#on regarde jusqu à deux ans de retard

#on a pmax=4;qmax=4

##  ETAPE 3 SELECTION DU MEILLEUR MODELE

pqs <- expand.grid(0:pmax,0:qmax) #combinaisons possibles de p<=p* et q<=q*
mat <- matrix(NA, nrow=pmax+1, ncol=pmax+1)
rownames(mat) <- paste0("p=",0:pmax) #renomme les lignes
colnames(mat) <- paste0("q=",0:pmax) #renomme les colonnes
AICs <- mat #matrice ou assigner les AIC
BICs <- mat #matrice ou assigner les BIC
for (row in 1:dim(pqs)[1]){
  p <- pqs[row,1]
  q <- pqs[row,2]
  estim <- try(arima(dspread,c(p,0,q), include.mean=F)) #tente d'estimer l'ARIMA
  AICs[p+1,q+1] <- if (class(estim)=="try-error") NA else estim$aic
  BICs[p+1,q+1] <- if (class(estim)=="try-error") NA else BIC(estim)
}
AICs
BICs
AICs==min(AICs) # AIC min p=4 et q=3

BICs==min(BICs) # BIC min p=0 et q=1


# CONLUSION : arima (0,1,1) et arima(4,1,3) sont candidats
arima011 <- arima(spread,c(0,1,1),include.mean=F)
arima011
arima413<- arima(spread,c(4,1,3),include.mean=F)
arima413

#on évalue quel modèle a les meilleurs capacités prédictives
# sur la base du R carré ajusté
adj_r2 <- function(model){
  ss_res <- sum(model$residuals^2)
  ss_tot <- sum(dspread[-c(1:max(p,q))]^2)
  p <- model$arma[1]
  q <- model$arma[2]
  n <- model$nobs-max(p,q)
  adj_r2 <- 1-(ss_res/(n-p-q-1))/(ss_tot/(n-1))
  return(adj_r2)
}

adj_r2(arima011)
adj_r2(arima413)

# l'écart de prédiction n'est pas très important
# nous prenons le arima(0,1,1) car plus parcimonieux

####################################################################
# PARTIE 3: Prevision 
####################################################################

forecast <- forecast(arima011, h = 2)  # Predire les 2 prochaines valeurs

# Afficher la prediction avec la region de confiance
png("Prevision.png")
plot(forecast, main="Prediction avec région de confiance")
dev.off()

