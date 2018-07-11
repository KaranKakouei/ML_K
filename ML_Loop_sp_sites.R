
rm(list=ls(all.names=TRUE))

if (!require("eHOF")) { install.packages("eHOF", dependencies = TRUE) ; library(eHOF)}
if (!require("foreach")) { install.packages("foreach", dependencies = TRUE) ; library(foreach)}
if (!require("doParallel")) { install.packages("doParallel", dependencies = TRUE) ; library(doParallel)}
if (!require("reshape2")) { install.packages("reshape2", dependencies = TRUE) ; library(reshape2)}
if (!require("splines")) { install.packages("splines", dependencies = TRUE) ; library(splines)}
if (!require("graphics")) { install.packages("graphics", dependencies = TRUE) ; library(graphics)}
if (!require("ggplot2")) { install.packages("ggplot2", dependencies = TRUE) ; library(ggplot2)}
if (!require("gridExtra")) { install.packages("gridExtra", dependencies = TRUE) ; library(gridExtra)}
if (!require("data.table")) { install.packages("data.table", dependencies = TRUE) ; library(data.table)}
if (!require("stringr")) { install.packages("stringr", dependencies = TRUE) ; library(stringr)}

# define the predict function
"predict.HOF" <-
  function (object, model, newdata, ...) {
    if(missing(model)) model <- pick.model(object, ...)
    p <- coef(object, model, ...)
    xrange <- object$range
    if (missing(newdata)) x <- object$x else x <- newdata
    M = max(Dat$y)
    fv <- HOF.fun(x=x, model=as.character(model), p=as.numeric(p), M=M, xrange)
    fv
  }

## Ref: 
# eHOF models
# Jansen F, Oksanen J (2013) How to model species responses along ecological
# gradients - Huisman-Olff-Fresco models revisited.  Journal of Vegetation Science
# 24, 1108-1117.

#################################################################################
####---- loop for model fits. Kinzig_dh4 ---- Treene_dh4
#################################################################################

# Read all files / dataframes
setwd("~/Directory/")
files <- list.files(pattern = '*.csv')
files

# Activate 16 CPUs for parallel analysis and time consuming
cl <- makePSOCKcluster(16, outfile="")
registerDoParallel(cl)
getDoParWorkers()

# Loop to set up models (statistical relationships) which allow predicting future values based on new data (Machine learning)
# Parallel loops to run the script on 16 dataframes
foreach (f=1:length(files), .packages=c("zoo", "eHOF", "data.table", "stringr")) %dopar% {
  
  # Reading the dataframe
  data <- fread(files[f], header=T)
  data <- as.data.frame(data)
  cols <- colnames(data)[313:447]
  
  # foreach (c=1:length(cols), .packages=c("zoo", "eHOF", "ggplot2")) %dopar% {
  for (c in 1:length(cols)) {
    
    # Define Environmental variable
    Env.var <- c("Env_flow")
    
    myHOF = NULL
    
    # Define columns
    colResponse=which(names(data)==cols[c])
    colEnvVAR=which(names(data)==Env.var)
    data[is.na(data)] <- 0
    newdata=data
    
    # Remove NAs
    if(any(is.na(data[, colResponses])) == TRUE) {
        newdata <- data[ -which(is.na(data[, colResponse])), ]
      } else {
        newdata=data
      }

      site_id <- newdata$id
      
      # Define x-axis values for three periods of baseline, horizon 2050 and horizon 2090
      Dat <- NULL
      Dat$EnvVar_1_2009 <- newdata$EnvVar_1_2009
      Dat$EnvVar_1_2010 <- newdata$EnvVar_1_2010
      Dat$EnvVar_1_2011 <- newdata$EnvVar_1_2011
      Dat$EnvVar_1_2012 <- newdata$EnvVar_1_2012
      Dat$EnvVar_1_2013 <- newdata$EnvVar_1_2013
      Dat$EnvVar_1_2014 <- newdata$EnvVar_1_2014
      Dat$EnvVar_1_2015 <- newdata$EnvVar_1_2015
      Dat$EnvVar_1_2016 <- newdata$EnvVar_1_2016
      Dat$EnvVar_1_2017 <- newdata$EnvVar_1_2017
      Dat$EnvVar_1_2018 <- newdata$EnvVar_1_2018
      Dat$EnvVar_1_2046 <- newdata$EnvVar_1_2046
      Dat$EnvVar_1_2047 <- newdata$EnvVar_1_2047
      Dat$EnvVar_1_2048 <- newdata$EnvVar_1_2048
      Dat$EnvVar_1_2049 <- newdata$EnvVar_1_2049
      Dat$EnvVar_1_2050 <- newdata$EnvVar_1_2050
      Dat$EnvVar_1_2051 <- newdata$EnvVar_1_2051
      Dat$EnvVar_1_2052 <- newdata$EnvVar_1_2052
      Dat$EnvVar_1_2053 <- newdata$EnvVar_1_2053
      Dat$EnvVar_1_2054 <- newdata$EnvVar_1_2054
      Dat$EnvVar_1_2055 <- newdata$EnvVar_1_2055
      Dat$EnvVar_1_2090 <- newdata$EnvVar_1_2090
      Dat$EnvVar_1_2091 <- newdata$EnvVar_1_2091
      Dat$EnvVar_1_2092 <- newdata$EnvVar_1_2092
      Dat$EnvVar_1_2093 <- newdata$EnvVar_1_2093
      Dat$EnvVar_1_2094 <- newdata$EnvVar_1_2094
      Dat$EnvVar_1_2095 <- newdata$EnvVar_1_2095
      Dat$EnvVar_1_2096 <- newdata$EnvVar_1_2096
      Dat$EnvVar_1_2097 <- newdata$EnvVar_1_2097
      Dat$EnvVar_1_2098 <- newdata$EnvVar_1_2098
      Dat$EnvVar_1_2099 <- newdata$EnvVar_1_2099
      
      # Set up the statistical/predictive relationships
      Dat$x_current <- newdata[,colEnvVAR]
      Dat$x_current <- round(Dat$x_current, digits = 2)
      Dat$y_response <- round(newdata[,colResponse], digits=1)
      Dat$y_response <- Dat$y_response
      motyp <- eHOF.modelnames[1:7] # Define model types (7 algorythms / models from low to high complexity)
      myHOF <- HOF(Dat$y_response, Dat$x_current, M = max(Dat$y_response),
                   family = poisson, bootstrap = 100, test = 'AICc', modeltypes=motyp, freq.limit=10) 
  
      # plot(myHOF, marginal='point', para=TRUE, onlybest=T, newdata=seq(min(myHOF$range), max(myHOF$range),length.out=10000))
      
      # Get the best selected model
      para <- Para(myHOF, newdata=seq(min(myHOF$range), max(myHOF$range),length.out=10000))
      info <- as.data.frame(unlist(para))
      info
      m <- as.character(droplevels(info[which(rownames(info)=="model"),]))
      
      # Get fitted values
      fit <- as.data.frame(fitted.values(myHOF))
      fit <- fit[paste(m)]
      colnames(fit) = "fit_current"
      
      # Add x-axis values
      fit$x_current <- Dat$x_current
      
      # Add IDs
      fit$id <- site_id
      
      # Predictions for baseline period
      fit$pred_EnvVar_1_2009 <- as.data.frame(predict.HOF(myHOF, model=paste(m), newdata=Dat$EnvVar_1_2009))
      fit$pred_EnvVar_1_2010 <- as.data.frame(predict.HOF(myHOF, model=paste(m), newdata=Dat$EnvVar_1_2010))
      fit$pred_EnvVar_1_2011 <- as.data.frame(predict.HOF(myHOF, model=paste(m), newdata=Dat$EnvVar_1_2011))
      fit$pred_EnvVar_1_2012 <- as.data.frame(predict.HOF(myHOF, model=paste(m), newdata=Dat$EnvVar_1_2012))
      fit$pred_EnvVar_1_2013 <- as.data.frame(predict.HOF(myHOF, model=paste(m), newdata=Dat$EnvVar_1_2013))
      fit$pred_EnvVar_1_2014 <- as.data.frame(predict.HOF(myHOF, model=paste(m), newdata=Dat$EnvVar_1_2014))
      fit$pred_EnvVar_1_2015 <- as.data.frame(predict.HOF(myHOF, model=paste(m), newdata=Dat$EnvVar_1_2015))
      fit$pred_EnvVar_1_2016 <- as.data.frame(predict.HOF(myHOF, model=paste(m), newdata=Dat$EnvVar_1_2016))
      fit$pred_EnvVar_1_2017 <- as.data.frame(predict.HOF(myHOF, model=paste(m), newdata=Dat$EnvVar_1_2017))
      fit$pred_EnvVar_1_2018 <- as.data.frame(predict.HOF(myHOF, model=paste(m), newdata=Dat$EnvVar_1_2018))

      # Predictions for horizon 2050
      fit$pred_EnvVar_1_2046 <- as.data.frame(predict.HOF(myHOF, model=paste(m), newdata=Dat$EnvVar_1_2046))
      fit$pred_EnvVar_1_2047 <- as.data.frame(predict.HOF(myHOF, model=paste(m), newdata=Dat$EnvVar_1_2047))
      fit$pred_EnvVar_1_2048 <- as.data.frame(predict.HOF(myHOF, model=paste(m), newdata=Dat$EnvVar_1_2048))
      fit$pred_EnvVar_1_2049 <- as.data.frame(predict.HOF(myHOF, model=paste(m), newdata=Dat$EnvVar_1_2049))
      fit$pred_EnvVar_1_2050 <- as.data.frame(predict.HOF(myHOF, model=paste(m), newdata=Dat$EnvVar_1_2050))
      fit$pred_EnvVar_1_2051 <- as.data.frame(predict.HOF(myHOF, model=paste(m), newdata=Dat$EnvVar_1_2051))
      fit$pred_EnvVar_1_2052 <- as.data.frame(predict.HOF(myHOF, model=paste(m), newdata=Dat$EnvVar_1_2052))
      fit$pred_EnvVar_1_2053 <- as.data.frame(predict.HOF(myHOF, model=paste(m), newdata=Dat$EnvVar_1_2053))
      fit$pred_EnvVar_1_2054 <- as.data.frame(predict.HOF(myHOF, model=paste(m), newdata=Dat$EnvVar_1_2054))
      fit$pred_EnvVar_1_2055 <- as.data.frame(predict.HOF(myHOF, model=paste(m), newdata=Dat$EnvVar_1_2055))

      # Predictions for horizon 2090
      fit$pred_EnvVar_1_2090 <- as.data.frame(predict.HOF(myHOF, model=paste(m), newdata=Dat$EnvVar_1_2090))
      fit$pred_EnvVar_1_2091 <- as.data.frame(predict.HOF(myHOF, model=paste(m), newdata=Dat$EnvVar_1_2091))
      fit$pred_EnvVar_1_2092 <- as.data.frame(predict.HOF(myHOF, model=paste(m), newdata=Dat$EnvVar_1_2092))
      fit$pred_EnvVar_1_2093 <- as.data.frame(predict.HOF(myHOF, model=paste(m), newdata=Dat$EnvVar_1_2093))
      fit$pred_EnvVar_1_2094 <- as.data.frame(predict.HOF(myHOF, model=paste(m), newdata=Dat$EnvVar_1_2094))
      fit$pred_EnvVar_1_2095 <- as.data.frame(predict.HOF(myHOF, model=paste(m), newdata=Dat$EnvVar_1_2095))
      fit$pred_EnvVar_1_2096 <- as.data.frame(predict.HOF(myHOF, model=paste(m), newdata=Dat$EnvVar_1_2096))
      fit$pred_EnvVar_1_2097 <- as.data.frame(predict.HOF(myHOF, model=paste(m), newdata=Dat$EnvVar_1_2097))
      fit$pred_EnvVar_1_2098 <- as.data.frame(predict.HOF(myHOF, model=paste(m), newdata=Dat$EnvVar_1_2098))
      fit$pred_EnvVar_1_2099 <- as.data.frame(predict.HOF(myHOF, model=paste(m), newdata=Dat$EnvVar_1_2099))
      
      # Add colnames
      colnames(fit[,4]) <- "pred_EnvVar_1_2009"
      colnames(fit[,5]) <- "pred_EnvVar_1_2010"
      colnames(fit[,6]) <- "pred_EnvVar_1_2011"
      colnames(fit[,7]) <- "pred_EnvVar_1_2012"
      colnames(fit[,8]) <- "pred_EnvVar_1_2013"
      colnames(fit[,9]) <- "pred_EnvVar_1_2014"
      colnames(fit[,10]) <- "pred_EnvVar_1_2015"
      colnames(fit[,11]) <- "pred_EnvVar_1_2016"
      colnames(fit[,12]) <- "pred_EnvVar_1_2017"
      colnames(fit[,13]) <- "pred_EnvVar_1_2018"
      colnames(fit[,14]) <- "pred_EnvVar_1_2046"
      colnames(fit[,15]) <- "pred_EnvVar_1_2047"
      colnames(fit[,16]) <- "pred_EnvVar_1_2048"
      colnames(fit[,17]) <- "pred_EnvVar_1_2049"
      colnames(fit[,18]) <- "pred_EnvVar_1_2050"
      colnames(fit[,19]) <- "pred_EnvVar_1_2051"
      colnames(fit[,20]) <- "pred_EnvVar_1_2052"
      colnames(fit[,21]) <- "pred_EnvVar_1_2053"
      colnames(fit[,22]) <- "pred_EnvVar_1_2054"
      colnames(fit[,23]) <- "pred_EnvVar_1_2055"
      colnames(fit[,24]) <- "pred_EnvVar_1_2090"
      colnames(fit[,25]) <- "pred_EnvVar_1_2091"
      colnames(fit[,26]) <- "pred_EnvVar_1_2092"
      colnames(fit[,27]) <- "pred_EnvVar_1_2093"
      colnames(fit[,28]) <- "pred_EnvVar_1_2094"
      colnames(fit[,29]) <- "pred_EnvVar_1_2095"
      colnames(fit[,30]) <- "pred_EnvVar_1_2096"
      colnames(fit[,31]) <- "pred_EnvVar_1_2097"
      colnames(fit[,32]) <- "pred_EnvVar_1_2098"
      colnames(fit[,33]) <- "pred_EnvVar_1_2099"
      head(fit)
      
      # get specific part of the file name
      files[f]
      name <- str_sub(files[f], start=27, end=-5)
      name
      
      # Save the predicted values
      write.csv(fit, paste0("~/Directory/", f, "/","EnvVar_1_", name, "_", cols[c], ".csv"))

  }
}

# Stop clustering process
stopCluster(cl)

###############################################################################################################
####---- merging all tables produced in the previous step
###############################################################################################################

CC_model <- c(1:16)

for(s in 1:length(CC_model)) {
  cat("Running site", s, "\n")
  
  # load the main table with IDs
  rm(table)
  setwd("~/Directory/")
  table <- read.csv2("Main_table.csv")
  
  # Read all tables
  setwd(paste0("~/Directory/EnvVar_1_/", s, "/"))
  files <- list.files(pattern = '*.csv')
  files
  
  # loop to merge tables based on IDs
  for(f in 1:length(files)) {
    cat("Running site", f, "\n")
    s <- read.csv(files[[f]])
    
    s2=NULL
    s2$id <- s[,4]
    s2$baseline <- rowMeans(s[, 6:25])
    s2$H2050 <- rowMeans(s[, 26:45])
    s2$H2090 <- rowMeans(s[, 46:65])
    s2 <- as.data.frame(s2)
    head(s2)
    
    name <- str_sub(files[f], start=79, end=-5)

    colnames(s2)[which(names(s2) == "baseline")] <- paste0("th3", "_", name, "_baseline")
    colnames(s2)[which(names(s2) == "H2050")] <- paste0("th3", "_", name, "_H2050")
    colnames(s2)[which(names(s2) == "H2090")] <- paste0("th3", "_", name, "_H2090")
    head(s2)
    
    table <- merge(table, s2, by.x="id", by.y="id", all.x=T)
    head(table)
    
    rm(s2)
  }
  
  # get specific part of file name
  files[f]
  name2 <- str_sub(files[f], start=5, end=-19)

  # save the merged dataframe
  write.csv(table, paste0("~/Directory/", "EnvVar_1_", name2, ".csv"))
  
}



###############################################################################################################

# End


