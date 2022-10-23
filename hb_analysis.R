
#-----------------------------prepare data--------------------------------------
#install.packages("RSGHB")
#install.packages("dplyr")
#install.packages("CRAN")
library(readr)
library(RSGHB)
library(dplyr)
library(tidyfst)
library(glue)
library(stringr) 

# set working directory 
# repo root
repo_root <- getwd()

# get cleaned data 
df <- read.csv(list.files("data_clean/apple/", full.names = T))
# manipulate data to feed the HB model
set_id = rep(1:6, times = nrow(df)/6)
df$S_id <- set_id
df = subset(df, select = -c(1))

# rename columns
colnames(df)[1] <- 'ID'
colnames(df)[2] <- 'Choice'

# set working dir to output folder
setwd(file.path(getwd(), "output/"))

#----------------------------select treatment-----------------------------------
df<-df[!(df$ID==16 ),]  
df<-df[(df$treatment_1==1 ),]  

#-----------------------set GMO and in-store as -1------------------------------

x_1 <- c("organic_1", "nongmo_1", "CRISPR_1", "conventional_1")
x_2 <- c("organic_2", "nongmo_2", "CRISPR_2", "conventional_2")
x_3 <- c("organic_3", "nongmo_3", "CRISPR_3", "conventional_3")
x_4 <- c("organic_4", "nongmo_4", "CRISPR_4", "conventional_4")

df2 <- df %>% mutate_at(.vars = x_1, funs(ifelse(GMO_1 == 1, -1, .)))
df3 <- df2 %>% mutate_at(.vars = x_2, funs(ifelse(GMO_2 == 1, -1, .)))
df4 <- df3 %>% mutate_at(.vars = x_3, funs(ifelse(GMO_3 == 1, -1, .)))
df5 <- df4 %>% mutate_at(.vars = x_4, funs(ifelse(GMO_4 == 1, -1, .)))

x=1:4
online_all = sapply(x, function(xi) {glue("online_{xi}")}) 
df5[, online_all][df5[, online_all] == 0] <- -1



df <- df5 #[(df5$treatment_1==1),] 


#--------------------------start analysis---------------------------------------

# SETTINGS for RSGHB 
modelname <- "hb_likelihood"       # used for output
gNCREP <- 10000    # Number of burn-in iterations 
gNEREP <- 1000  	# Number of iterations to keep for averaging after convergence has been reached
gNSKIP <- 5  	# Number of iterations to do in between retaining draws for averaging
gINFOSKIP <- 250    	# How frequently to print info about the iteration process
degreesOfFreedom <- 5 # Additional degrees of freedom for the prior covariance matrix (not including the number of parameters. (Defaults to 5)
#priorVariance <- c(1,1,1,1,1,1, 0.5) # last param is for price (Defaults to 2.0)
priorVariance <- 1
writeModel <- T # Write model files?
gStoreDraws <- TRUE # Store random draws?

# The choice vectors: Dummy coding the choice vector allows for 
# easier coding of the the likelihood calculations.
df$choice1  <- (df$Choice==1)
df$choice2  <- (df$Choice==2)
df$choice3  <- (df$Choice==3)
df$choice4  <- (df$Choice==4)
df$choice5  <- (df$Choice==5)


#df <- df %>%
#  mutate(price_1 = recode(price_1, "1"= 1.99 , "2" = 2.99, "3" = 3.99, "4" = 4.99) 
#  )
#df <- df %>%
#  mutate(price_2 = recode(price_2, "1"= 1.99 , "2" = 2.99, "3" = 3.99, "4" = 4.99) 
#  )
#df <- df %>%
#  mutate(price_3 = recode(price_3, "1"= 1.99 , "2" = 2.99, "3" = 3.99, "4" = 4.99) 
#  )
#df <- df %>%
#  mutate(price_4 = recode(price_4, "1"= 1.99 , "2" = 2.99, "3" = 3.99, "4" = 4.99) 
#  )


# Fixed effects for images and incentive alignment (ia) conditions and for validation choices in these conditions
#gVarNamesFixed <- c( "btreatment_1" ,"btreatment_2")
# random effects
gVarNamesNormal <- c("none", "organic", "nongmo", "CRISPR", "conventional",'online' ,"price")

# Assuming normal distributions for all parameters, except price, which is one-sided triangular distribution
gDIST <- c(rep(1,length(gVarNamesNormal)-1), 4)	

# STARTING VALUES
#FC <- c(rep(0, length(gVarNamesFixed))) # for the fixed coefficients
svN <- c(rep(1, length(gVarNamesNormal) - 1), 1)  # for the random coefficients

# CONTROL LIST TO PASS TO doHB
control <- list(
  modelname=modelname,
  #gVarNamesFixed = gVarNamesFixed,
 # FC = FC,
  gVarNamesNormal=gVarNamesNormal,
  gDIST=gDIST,
#  svN=svN,
  gNCREP=gNCREP,
  gNEREP=gNEREP,
  gNSKIP=gNSKIP,
  gINFOSKIP=gINFOSKIP,
  gStoreDraws=gStoreDraws,  
  writeModel= writeModel,
  degreesOfFreedom=degreesOfFreedom,
  priorVariance=priorVariance
#  constraintsNorm=constraintsNorm
)


# likelihood function

likelihood<- function(fc ,b)
{
  # fixed effects for scale adjustment parameters
# cc <- 1
 # btreatment_1  <- fc[cc]; cc <- cc + 1
 # btreatment_2  <- fc[cc]; cc <- cc + 1

  # random effects
  cc = 1 # internal pointer
  bnone <- b[,cc];cc=cc+1
  borganic  <- b[,cc];cc=cc+1
  bnongmo <- b[,cc];cc=cc+1 
  bCRISPR <- b[,cc];cc=cc+1 
  bconventional <- b[,cc];cc=cc+1 
  bonline  <- b[,cc];cc=cc+1 
  bprice  <- b[cc]; cc <- cc + 1
  # in this stylized model bprice is normalized to 1.
  # for the Allenby et al. (2014) model, set bgamma to 1 and estimate a bprice parameter instead (adjust the v_i and h_i formulations below accordingly).
  # for  the Sonnier, Allenby, and Otter (2007) specification estimate bgamma as 1/b[,cc]
  
  
  
  # scale adjustment factors for the different experimental conditions
  #treatment_estimation <- exp(btreatment_1 * df$treatment_1 + btreatment_2 * df$treatment_2 )
  
  # utility functions for the regular choice sets
  # Utility functions 
  # "none", "organic", "nongmo", "CRISPR", "conventional",'online' ,"gamma" 
 
  v1 <- (borganic * df$organic_1 + bnongmo * df$nongmo_1 + bCRISPR * df$CRISPR_1 
             + bconventional * df$conventional_1 + bonline * df$online_1 - bprice * df$price_1)
          
  v2 <- (borganic * df$organic_2 + bnongmo * df$nongmo_2 + bCRISPR * df$CRISPR_2 
             + bconventional * df$conventional_2 + bonline * df$online_2 - bprice * df$price_2)
  v3 <- (borganic * df$organic_3 + bnongmo * df$nongmo_3 + bCRISPR * df$CRISPR_3 
             + bconventional * df$conventional_3 + bonline * df$online_3 - bprice * df$price_3)
  v4 <- (borganic * df$organic_4 + bnongmo * df$nongmo_4 + bCRISPR * df$CRISPR_4 
             + bconventional * df$conventional_4 + bonline * df$online_4 - bprice * df$price_4)
  
  v5 <- bprice * bnone

#  exp(v3)*df$choice3  
  p  <- (exp(v1)*df$choice1 + exp(v2) * df$choice2  + exp(v4) * df$choice4 + exp(v3) * df$choice3
         + exp(v5) * df$choice5) / (exp(v1) + exp(v2)+ exp(v3) + exp(v4) + exp(v5))
  
  return(p)
}

set.seed(2017)

doHB(likelihood, df, control)


#---------------Sample-level main effect coefficient means----------------------

model <- doHB(likelihood, df, control)

dfA <- read_csv("hb_likelihood_A.csv")

dfA_select <- dfA %>% select(2:8)

sapply(dfA_select, function(dfA_select) c( 
    "Stand dev" = sd(dfA_select), 
    "Mean"= mean(dfA_select,na.rm=TRUE),
    "n" = length(dfA_select),
    "Median" = median(dfA_select),
    "CoeffofVariation" = sd(dfA_select)/mean(dfA_select,na.rm=TRUE),
    "Minimum" = min(dfA_select),
    "Maximun" = max(dfA_select),
    "Upper Quantile" = quantile(dfA_select,1),
    "LowerQuartile" = quantile(dfA_select,0)
  )
)



# B Average individual level draws 
head(model[["B"]]) # not use in this project


# end 







