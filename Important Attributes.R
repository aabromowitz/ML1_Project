# The purpose of this script is to figure out which variables are the "5 best" to get plots for.
# This project is supposed to be in python though, so the plotting will be in python.

# Packages
library(caret)
library(pROC)
library(car)

# Pull in Data
data<-read.csv("https://raw.githubusercontent.com/aabromowitz/ML1_Project/main/AIDS_Classification_50000.csv",stringsAsFactors = T)

# Calculate AUC scores for the 22 columns
data$trt[data$trt == 0] <- 'ZDV'
data$trt[data$trt == 1] <- 'ZDV + ddl'
data$trt[data$trt == 2] <- 'ZDV + Zal'
data$trt[data$trt == 3] <- 'ddl only'
data$trt <- as.factor(data$trt)
data$hemo[data$hemo == 0] <- 'no'
data$hemo[data$hemo == 1] <- 'yes'
data$hemo <- as.factor(data$hemo)
data$homo[data$homo == 0] <- 'no'
data$homo[data$homo == 1] <- 'yes'
data$homo <- as.factor(data$homo)
data$drugs[data$drugs == 0] <- 'no'
data$drugs[data$drugs == 1] <- 'yes'
data$drugs <- as.factor(data$drugs)
data$karnof <- as.factor(data$karnof)
data$oprior[data$oprior == 0] <- 'no'
data$oprior[data$oprior == 1] <- 'yes'
data$oprior <- as.factor(data$oprior)
data$z30[data$z30 == 0] <- 'no'
data$z30[data$z30 == 1] <- 'yes'
data$z30 <- as.factor(data$z30)
data$race[data$race == 0] <- 'White'
data$race[data$race == 1] <- 'non-white'
data$race <- as.factor(data$race)
data$gender[data$gender == 0] <- 'F'
data$gender[data$gender == 1] <- 'M'
data$gender <- as.factor(data$gender)
data$str2[data$str2 == 0] <- 'naive'
data$str2[data$str2 == 1] <- 'experienced'
data$str2 <- as.factor(data$str2)
data$strat[data$strat == 1] <- 'Antiretroviral Naive'
data$strat[data$strat == 2] <- '> 1 but <= 52 weeks of prior antiretroviral therapy'
data$strat[data$strat == 3] <- '> 52 weeks'
data$strat <- as.factor(data$strat)
data$symptom[data$symptom == 0] <- 'asymp'
data$symptom[data$symptom == 1] <- 'symp'
data$symptom <- as.factor(data$symptom)
data$treat[data$treat == 0] <- 'ZDV'
data$treat[data$treat == 1] <- 'others'
data$treat <- as.factor(data$treat)
data$offtrt[data$offtrt == 0] <- 'no'
data$offtrt[data$offtrt == 1] <- 'yes'
data$offtrt <- as.factor(data$offtrt)
data$infected[data$infected == 0] <- 'No'
data$infected[data$infected == 1] <- 'Yes'
data$infected <- as.factor(data$infected)

# Calculate AUC scores for the 22 columns
set.seed(1)
vars <- names(data)
vars <- vars[vars!="infected"]
num_vars <- length(vars)
var_aucs <- data.frame("vars" = vars)
num_folds <- 10
for (j in 1:num_vars) {
  var <- vars[j]
  print(paste(j,'/',num_vars,': ',var,sep=''))
  folds <- createFolds(data$infected, k = num_folds)
  auc_scores <- numeric(num_folds)
  for (i in 1:num_folds) {
    tryCatch({
      train_indices <- unlist(folds[-i])
      test_indices <- unlist(folds[i])
      train <- data[train_indices, ]
      test <- data[test_indices, ]
      form <- as.formula(paste("infected ~ ",var,sep=""))
      model <- glm(form, data = train, family = "binomial")
      predictions <- predict(model, newdata = test, type = "response")
      roc <- roc(response=test$infected,predictor=predictions,levels=c("No", "Yes"),direction = "<")
      auc_scores[i] <- auc(roc)
    }, error = function(e) { 
      auc_scores <- auc_scores[-i]
    })
  }
  var_aucs$auc[var_aucs$var == var] <- mean(auc_scores)
}
# strat = 0.6464335

# Add to strat
set.seed(2)
vars <- names(data)
vars <- vars[vars!="infected"]
vars <- vars[vars!="strat"]
num_vars <- length(vars)
var_aucs <- data.frame("vars" = vars)
num_folds <- 10
for (j in 1:num_vars) {
  var <- vars[j]
  print(paste(j,'/',num_vars,': ',var,sep=''))
  folds <- createFolds(data$infected, k = num_folds)
  auc_scores <- numeric(num_folds)
  for (i in 1:num_folds) {
    tryCatch({
      train_indices <- unlist(folds[-i])
      test_indices <- unlist(folds[i])
      train <- data[train_indices, ]
      test <- data[test_indices, ]
      form <- as.formula(paste("infected ~ strat + ",var,sep=""))
      model <- glm(form, data = train, family = "binomial")
      predictions <- predict(model, newdata = test, type = "response")
      roc <- roc(response=test$infected,predictor=predictions,levels=c("No", "Yes"),direction = "<")
      auc_scores[i] <- auc(roc)
    }, error = function(e) { 
      auc_scores <- auc_scores[-i]
    })
  }
  var_aucs$auc[var_aucs$var == var] <- mean(auc_scores)
}
# cd420 = 0.6744314

# Add to strat + cd420
set.seed(3)
vars <- names(data)
vars <- vars[vars!="infected"]
vars <- vars[vars!="strat"]
vars <- vars[vars!="cd420"]
num_vars <- length(vars)
var_aucs <- data.frame("vars" = vars)
num_folds <- 10
for (j in 1:num_vars) {
  var <- vars[j]
  print(paste(j,'/',num_vars,': ',var,sep=''))
  folds <- createFolds(data$infected, k = num_folds)
  auc_scores <- numeric(num_folds)
  for (i in 1:num_folds) {
    tryCatch({
      train_indices <- unlist(folds[-i])
      test_indices <- unlist(folds[i])
      train <- data[train_indices, ]
      test <- data[test_indices, ]
      form <- as.formula(paste("infected ~ strat + cd420 + ",var,sep=""))
      model <- glm(form, data = train, family = "binomial")
      predictions <- predict(model, newdata = test, type = "response")
      roc <- roc(response=test$infected,predictor=predictions,levels=c("No", "Yes"),direction = "<")
      auc_scores[i] <- auc(roc)
    }, error = function(e) { 
      auc_scores <- auc_scores[-i]
    })
  }
  var_aucs$auc[var_aucs$var == var] <- mean(auc_scores)
}
# z30 = 0.6850368

# Try removing a variable
set.seed(4)
start_form_str <- 'infected ~ strat + cd420 + z30'
vars <- unlist(strsplit(start_form_str, "\\+"))
vars <- trimws(vars)
vars[1] <- substr(vars[1], 12, nchar(vars[1]))
num_vars <- length(vars)
var_aucs <- data.frame("vars" = vars)
num_folds <- 10
for (j in 1:num_vars) {
  var <- vars[j]
  print(var)
  folds <- createFolds(data$infected, k = num_folds)
  auc_scores <- numeric(num_folds)
  for (i in 1:num_folds) {
    train_indices <- unlist(folds[-i])
    test_indices <- unlist(folds[i])
    train <- data[train_indices, ]
    test <- data[test_indices, ]
    form <- as.formula(paste(start_form_str," - ",var,sep=""))
    model <- glm(form, data = train, family = "binomial")
    predictions <- predict(model, newdata = test, type = "response")
    roc <- roc(response=test$infected,predictor=predictions,levels=c("No", "Yes"),direction = "<")
    auc_scores[i] <- auc(roc)
  }
  var_aucs$auc[var_aucs$var == var] <- mean(auc_scores)
} 
# Didn't remove anything

# Add to strat + cd420 + z30
set.seed(5)
vars <- names(data)
vars <- vars[vars!="infected"]
vars <- vars[vars!="strat"]
vars <- vars[vars!="cd420"]
vars <- vars[vars!="z30"]
num_vars <- length(vars)
var_aucs <- data.frame("vars" = vars)
num_folds <- 10
for (j in 1:num_vars) {
  var <- vars[j]
  print(paste(j,'/',num_vars,': ',var,sep=''))
  folds <- createFolds(data$infected, k = num_folds)
  auc_scores <- numeric(num_folds)
  for (i in 1:num_folds) {
    tryCatch({
      train_indices <- unlist(folds[-i])
      test_indices <- unlist(folds[i])
      train <- data[train_indices, ]
      test <- data[test_indices, ]
      form <- as.formula(paste("infected ~ strat + cd420 + z30 + ",var,sep=""))
      model <- glm(form, data = train, family = "binomial")
      predictions <- predict(model, newdata = test, type = "response")
      roc <- roc(response=test$infected,predictor=predictions,levels=c("No", "Yes"),direction = "<")
      auc_scores[i] <- auc(roc)
    }, error = function(e) { 
      auc_scores <- auc_scores[-i]
    })
  }
  var_aucs$auc[var_aucs$var == var] <- mean(auc_scores)
}
# time = 0.6902066

# Try removing a variable
set.seed(6)
start_form_str <- 'infected ~ strat + cd420 + z30 + time'
vars <- unlist(strsplit(start_form_str, "\\+"))
vars <- trimws(vars)
vars[1] <- substr(vars[1], 12, nchar(vars[1]))
num_vars <- length(vars)
var_aucs <- data.frame("vars" = vars)
num_folds <- 10
for (j in 1:num_vars) {
  var <- vars[j]
  print(var)
  folds <- createFolds(data$infected, k = num_folds)
  auc_scores <- numeric(num_folds)
  for (i in 1:num_folds) {
    train_indices <- unlist(folds[-i])
    test_indices <- unlist(folds[i])
    train <- data[train_indices, ]
    test <- data[test_indices, ]
    form <- as.formula(paste(start_form_str," - ",var,sep=""))
    model <- glm(form, data = train, family = "binomial")
    predictions <- predict(model, newdata = test, type = "response")
    roc <- roc(response=test$infected,predictor=predictions,levels=c("No", "Yes"),direction = "<")
    auc_scores[i] <- auc(roc)
  }
  var_aucs$auc[var_aucs$var == var] <- mean(auc_scores)
} 
# Didn't remove anything

# Add to strat + cd420 + z30 + time
set.seed(7)
vars <- names(data)
vars <- vars[vars!="infected"]
vars <- vars[vars!="strat"]
vars <- vars[vars!="cd420"]
vars <- vars[vars!="z30"]
vars <- vars[vars!="time"]
num_vars <- length(vars)
var_aucs <- data.frame("vars" = vars)
num_folds <- 10
for (j in 1:num_vars) {
  var <- vars[j]
  print(paste(j,'/',num_vars,': ',var,sep=''))
  folds <- createFolds(data$infected, k = num_folds)
  auc_scores <- numeric(num_folds)
  for (i in 1:num_folds) {
    tryCatch({
      train_indices <- unlist(folds[-i])
      test_indices <- unlist(folds[i])
      train <- data[train_indices, ]
      test <- data[test_indices, ]
      form <- as.formula(paste("infected ~ strat + cd420 + z30 + time + ",var,sep=""))
      model <- glm(form, data = train, family = "binomial")
      predictions <- predict(model, newdata = test, type = "response")
      roc <- roc(response=test$infected,predictor=predictions,levels=c("No", "Yes"),direction = "<")
      auc_scores[i] <- auc(roc)
    }, error = function(e) { 
      auc_scores <- auc_scores[-i]
    })
  }
  var_aucs$auc[var_aucs$var == var] <- mean(auc_scores)
}
# cd40 = 0.6944189

# Try removing a variable
set.seed(8)
start_form_str <- 'infected ~ strat + cd420 + z30 + time + cd40'
vars <- unlist(strsplit(start_form_str, "\\+"))
vars <- trimws(vars)
vars[1] <- substr(vars[1], 12, nchar(vars[1]))
num_vars <- length(vars)
var_aucs <- data.frame("vars" = vars)
num_folds <- 10
for (j in 1:num_vars) {
  var <- vars[j]
  print(var)
  folds <- createFolds(data$infected, k = num_folds)
  auc_scores <- numeric(num_folds)
  for (i in 1:num_folds) {
    train_indices <- unlist(folds[-i])
    test_indices <- unlist(folds[i])
    train <- data[train_indices, ]
    test <- data[test_indices, ]
    form <- as.formula(paste(start_form_str," - ",var,sep=""))
    model <- glm(form, data = train, family = "binomial")
    predictions <- predict(model, newdata = test, type = "response")
    roc <- roc(response=test$infected,predictor=predictions,levels=c("No", "Yes"),direction = "<")
    auc_scores[i] <- auc(roc)
  }
  var_aucs$auc[var_aucs$var == var] <- mean(auc_scores)
} 
# Didn't remove anything

# Add to strat + cd420 + z30 + time + cd40
set.seed(9)
vars <- names(data)
vars <- vars[vars!="infected"]
vars <- vars[vars!="strat"]
vars <- vars[vars!="cd420"]
vars <- vars[vars!="z30"]
vars <- vars[vars!="time"]
vars <- vars[vars!="cd40"]
num_vars <- length(vars)
var_aucs <- data.frame("vars" = vars)
num_folds <- 10
for (j in 1:num_vars) {
  var <- vars[j]
  print(paste(j,'/',num_vars,': ',var,sep=''))
  folds <- createFolds(data$infected, k = num_folds)
  auc_scores <- numeric(num_folds)
  for (i in 1:num_folds) {
    tryCatch({
      train_indices <- unlist(folds[-i])
      test_indices <- unlist(folds[i])
      train <- data[train_indices, ]
      test <- data[test_indices, ]
      form <- as.formula(paste("infected ~ strat + cd420 + z30 + time + ",var,sep=""))
      model <- glm(form, data = train, family = "binomial")
      predictions <- predict(model, newdata = test, type = "response")
      roc <- roc(response=test$infected,predictor=predictions,levels=c("No", "Yes"),direction = "<")
      auc_scores[i] <- auc(roc)
    }, error = function(e) { 
      auc_scores <- auc_scores[-i]
    })
  }
  var_aucs$auc[var_aucs$var == var] <- mean(auc_scores)
}
# cd40 = 0.6944189

# Looking at the model
model <- glm('infected ~ strat + cd420 + z30 + time + cd40', data = data, family = "binomial")
summary(model) # all low p values
vif(model) # low vif values

model <- glm('infected ~ strat + cd420 + z30 + time + cd40 + str2', data = data, family = "binomial")
summary(model) # all low p values
vif(model) # low vif values

model <- glm('infected ~ strat + cd420 + z30 + time + cd40 + str2 + trt', data = data, family = "binomial")
summary(model) # not good
vif(model) # low vif values

model <- glm('infected ~ strat + cd420 + z30 + time + cd40 + str2 + preanti', data = data, family = "binomial")
summary(model) # all low p values
vif(model) # low vif values

################################################################################
# Plotting some of the variables
plot(data$time)
hist(data$time) # interesting, violin plot, or histogram

plot(data$z30) 

plot(data$preanti)
hist(data$preanti)
boxplot(data$preanti) # interesting, boxplot or histogram

plot(data$str2)

plot(data$strat) # could be an example of a categorical plot

hist(data$cd40)
boxplot(data$cd40) # interesting

hist(data$cd420)
boxplot(data$cd420)

################################################################################
# Try figuring out if there are any useful interactions
set.seed(10)
vars <- names(data)
vars <- vars[vars!="infected"]
allVars <- vars
for (i in 1:length(vars)){
  for (j in 1:i){
    if(vars[i]!=vars[j]) {
      allVars <- c(allVars,paste(vars[i],'*',vars[j],sep=''))
    }
  }
}
num_vars <- length(allVars)
var_aucs <- data.frame("vars" = allVars)
num_folds <- 10
for (j in 1:num_vars) {
  var <- allVars[j]
  print(paste(j,'/',num_vars,': ',var,sep=''))
  folds <- createFolds(data$infected, k = num_folds)
  auc_scores <- numeric(num_folds)
  for (i in 1:num_folds) {
    tryCatch({
      train_indices <- unlist(folds[-i])
      test_indices <- unlist(folds[i])
      train <- data[train_indices, ]
      test <- data[test_indices, ]
      form <- as.formula(paste("infected ~ ",var,sep=""))
      model <- glm(form, data = train, family = "binomial")
      predictions <- predict(model, newdata = test, type = "response")
      roc <- roc(response=test$infected,predictor=predictions,levels=c("No", "Yes"),direction = "<")
      auc_scores[i] <- auc(roc)
    }, error = function(e) { 
      auc_scores <- auc_scores[-i]
    })
  }
  var_aucs$auc[var_aucs$var == var] <- mean(auc_scores)
}
# cd420*strat = 0.6741977

# Add to strat + cd420
set.seed(11)
vars <- names(data)
vars <- vars[vars!="infected"]
allVars <- vars
for (i in 1:length(vars)){
  for (j in 1:i){
    if(vars[i]!=vars[j]) {
      allVars <- c(allVars,paste(vars[i],'*',vars[j],sep=''))
    }
  }
}
allVars <- allVars[allVars!="strat"]
allVars <- allVars[allVars!="cd420"]
num_vars <- length(allVars)
var_aucs <- data.frame("vars" = allVars)
num_folds <- 10
for (j in 1:num_vars) {
  var <- allVars[j]
  print(paste(j,'/',num_vars,': ',var,sep=''))
  folds <- createFolds(data$infected, k = num_folds)
  auc_scores <- numeric(num_folds)
  for (i in 1:num_folds) {
    tryCatch({
      train_indices <- unlist(folds[-i])
      test_indices <- unlist(folds[i])
      train <- data[train_indices, ]
      test <- data[test_indices, ]
      form <- as.formula(paste("infected ~ strat + cd420 + ",var,sep=""))
      model <- glm(form, data = train, family = "binomial")
      predictions <- predict(model, newdata = test, type = "response")
      roc <- roc(response=test$infected,predictor=predictions,levels=c("No", "Yes"),direction = "<")
      auc_scores[i] <- auc(roc)
    }, error = function(e) { 
      auc_scores <- auc_scores[-i]
    })
  }
  var_aucs$auc[var_aucs$var == var] <- mean(auc_scores)
}
# z30*time = 0.6901645

# Add to strat + cd420 + z30 + time
set.seed(12)
vars <- names(data)
vars <- vars[vars!="infected"]
allVars <- vars
for (i in 1:length(vars)){
  for (j in 1:i){
    if(vars[i]!=vars[j]) {
      allVars <- c(allVars,paste(vars[i],'*',vars[j],sep=''))
    }
  }
}
allVars <- allVars[allVars!="strat"]
allVars <- allVars[allVars!="cd420"]
allVars <- allVars[allVars!="z30"]
allVars <- allVars[allVars!="time"]
num_vars <- length(allVars)
var_aucs <- data.frame("vars" = allVars)
num_folds <- 10
for (j in 1:num_vars) {
  var <- allVars[j]
  print(paste(j,'/',num_vars,': ',var,sep=''))
  folds <- createFolds(data$infected, k = num_folds)
  auc_scores <- numeric(num_folds)
  for (i in 1:num_folds) {
    tryCatch({
      train_indices <- unlist(folds[-i])
      test_indices <- unlist(folds[i])
      train <- data[train_indices, ]
      test <- data[test_indices, ]
      form <- as.formula(paste("infected ~ strat + cd420 + z30 + time + ",var,sep=""))
      model <- glm(form, data = train, family = "binomial")
      predictions <- predict(model, newdata = test, type = "response")
      roc <- roc(response=test$infected,predictor=predictions,levels=c("No", "Yes"),direction = "<")
      auc_scores[i] <- auc(roc)
    }, error = function(e) { 
      auc_scores <- auc_scores[-i]
    })
  }
  var_aucs$auc[var_aucs$var == var] <- mean(auc_scores)
}
# cd40*trt = 0.6969354

# Add to strat + cd420 + z30 + time + cd40*trt
set.seed(13)
vars <- names(data)
vars <- vars[vars!="infected"]
allVars <- vars
for (i in 1:length(vars)){
  for (j in 1:i){
    if(vars[i]!=vars[j]) {
      allVars <- c(allVars,paste(vars[i],'*',vars[j],sep=''))
    }
  }
}
allVars <- allVars[allVars!="strat"]
allVars <- allVars[allVars!="cd420"]
allVars <- allVars[allVars!="z30"]
allVars <- allVars[allVars!="time"]
allVars <- allVars[allVars!="cd40*trt"]
allVars <- allVars[allVars!="cd40"]
allVars <- allVars[allVars!="trt"]
num_vars <- length(allVars)
var_aucs <- data.frame("vars" = allVars)
num_folds <- 10
for (j in 1:num_vars) {
  var <- allVars[j]
  print(paste(j,'/',num_vars,': ',var,sep=''))
  folds <- createFolds(data$infected, k = num_folds)
  auc_scores <- numeric(num_folds)
  for (i in 1:num_folds) {
    tryCatch({
      train_indices <- unlist(folds[-i])
      test_indices <- unlist(folds[i])
      train <- data[train_indices, ]
      test <- data[test_indices, ]
      form <- as.formula(paste("infected ~ strat + cd420 + z30 + time + cd40*trt + ",var,sep=""))
      model <- glm(form, data = train, family = "binomial")
      predictions <- predict(model, newdata = test, type = "response")
      roc <- roc(response=test$infected,predictor=predictions,levels=c("No", "Yes"),direction = "<")
      auc_scores[i] <- auc(roc)
    }, error = function(e) { 
      auc_scores <- auc_scores[-i]
    })
  }
  var_aucs$auc[var_aucs$var == var] <- mean(auc_scores)
}
# treat*str2 = 0.7005556

# Try removing a variable
set.seed(14)
start_form_str <- 'infected ~ strat + cd420 + z30 + time + cd40 + trt+ cd40*trt + treat*str2 + treat + str2'
vars <- unlist(strsplit(start_form_str, "\\+"))
vars <- trimws(vars)
vars[1] <- substr(vars[1], 12, nchar(vars[1]))
num_vars <- length(vars)
var_aucs <- data.frame("vars" = vars)
num_folds <- 10
for (j in 1:num_vars) {
  var <- vars[j]
  print(var)
  folds <- createFolds(data$infected, k = num_folds)
  auc_scores <- numeric(num_folds)
  for (i in 1:num_folds) {
    train_indices <- unlist(folds[-i])
    test_indices <- unlist(folds[i])
    train <- data[train_indices, ]
    test <- data[test_indices, ]
    form <- as.formula(paste(start_form_str," - ",var,sep=""))
    model <- glm(form, data = train, family = "binomial")
    predictions <- predict(model, newdata = test, type = "response")
    roc <- roc(response=test$infected,predictor=predictions,levels=c("No", "Yes"),direction = "<")
    auc_scores[i] <- auc(roc)
  }
  var_aucs$auc[var_aucs$var == var] <- mean(auc_scores)
} 
# str2 = 0.7006131

model <- glm('infected ~ strat + cd420 + z30 + time + cd40 + trt+ cd40*trt + treat*str2 + treat + str2 - str2', data = data, family = "binomial")
summary(model)

# Try removing a variable
set.seed(15)
start_form_str <- 'infected ~ strat + cd420 + z30 + time + cd40 + trt+ cd40*trt + treat*str2 + treat + str2 - str2'
vars <- unlist(strsplit(start_form_str, "\\+"))
vars <- trimws(vars)
vars[1] <- substr(vars[1], 12, nchar(vars[1]))
num_vars <- length(vars)
var_aucs <- data.frame("vars" = vars)
num_folds <- 10
for (j in 1:num_vars) {
  var <- vars[j]
  print(var)
  folds <- createFolds(data$infected, k = num_folds)
  auc_scores <- numeric(num_folds)
  for (i in 1:num_folds) {
    train_indices <- unlist(folds[-i])
    test_indices <- unlist(folds[i])
    train <- data[train_indices, ]
    test <- data[test_indices, ]
    form <- as.formula(paste(start_form_str," - ",var,sep=""))
    model <- glm(form, data = train, family = "binomial")
    predictions <- predict(model, newdata = test, type = "response")
    roc <- roc(response=test$infected,predictor=predictions,levels=c("No", "Yes"),direction = "<")
    auc_scores[i] <- auc(roc)
  }
  var_aucs$auc[var_aucs$var == var] <- mean(auc_scores)
} 
# Didn't remove anything

# Add to strat + cd420 + z30 + time + cd40 + trt + cd40*trt + treat*str2 - str2
set.seed(16)
vars <- names(data)
vars <- vars[vars!="infected"]
allVars <- vars
for (i in 1:length(vars)){
  for (j in 1:i){
    if(vars[i]!=vars[j]) {
      allVars <- c(allVars,paste(vars[i],'*',vars[j],sep=''))
    }
  }
}
allVars <- allVars[allVars!="strat"]
allVars <- allVars[allVars!="cd420"]
allVars <- allVars[allVars!="z30"]
allVars <- allVars[allVars!="time"]
allVars <- allVars[allVars!="cd40*trt"]
allVars <- allVars[allVars!="cd40"]
allVars <- allVars[allVars!="trt"]
allVars <- allVars[allVars!="treat*str2"]
allVars <- allVars[allVars!="treat"]
num_vars <- length(allVars)
var_aucs <- data.frame("vars" = allVars)
num_folds <- 10
for (j in 1:num_vars) {
  var <- allVars[j]
  print(paste(j,'/',num_vars,': ',var,sep=''))
  folds <- createFolds(data$infected, k = num_folds)
  auc_scores <- numeric(num_folds)
  for (i in 1:num_folds) {
    tryCatch({
      train_indices <- unlist(folds[-i])
      test_indices <- unlist(folds[i])
      train <- data[train_indices, ]
      test <- data[test_indices, ]
      form <- as.formula(paste("infected ~ strat + cd420 + z30 + time + cd40*trt + treat*str2 - str2 + ",var,sep=""))
      model <- glm(form, data = train, family = "binomial")
      predictions <- predict(model, newdata = test, type = "response")
      roc <- roc(response=test$infected,predictor=predictions,levels=c("No", "Yes"),direction = "<")
      auc_scores[i] <- auc(roc)
    }, error = function(e) { 
      auc_scores <- auc_scores[-i]
    })
  }
  var_aucs$auc[var_aucs$var == var] <- mean(auc_scores)
}
# offtrt*preanti = 0.7018407

# Try removing a variable
set.seed(17)
start_form_str <- 'infected ~ strat + cd420 + z30 + time + cd40 + trt+ cd40*trt + treat*str2 + treat - str2 + offtrt*preanti + offtrt + preanti'
vars <- unlist(strsplit(start_form_str, "[\\+\\-]"))
vars <- trimws(vars)
vars[1] <- substr(vars[1], 12, nchar(vars[1]))
num_vars <- length(vars)
var_aucs <- data.frame("vars" = vars)
num_folds <- 10
for (j in 1:num_vars) {
  var <- vars[j]
  print(var)
  folds <- createFolds(data$infected, k = num_folds)
  auc_scores <- numeric(num_folds)
  for (i in 1:num_folds) {
    train_indices <- unlist(folds[-i])
    test_indices <- unlist(folds[i])
    train <- data[train_indices, ]
    test <- data[test_indices, ]
    form <- as.formula(paste(start_form_str," - ",var,sep=""))
    model <- glm(form, data = train, family = "binomial")
    predictions <- predict(model, newdata = test, type = "response")
    roc <- roc(response=test$infected,predictor=predictions,levels=c("No", "Yes"),direction = "<")
    auc_scores[i] <- auc(roc)
  }
  var_aucs$auc[var_aucs$var == var] <- mean(auc_scores)
} 
# preanti = 0.7018475

# Try removing a variable
set.seed(18)
start_form_str <- 'infected ~ strat + cd420 + z30 + time + cd40 + trt+ cd40*trt + treat*str2 + treat - str2 + offtrt*preanti + offtrt - preanti'
vars <- unlist(strsplit(start_form_str, "[\\+\\-]"))
vars <- trimws(vars)
vars[1] <- substr(vars[1], 12, nchar(vars[1]))
num_vars <- length(vars)
var_aucs <- data.frame("vars" = vars)
num_folds <- 10
for (j in 1:num_vars) {
  var <- vars[j]
  print(var)
  folds <- createFolds(data$infected, k = num_folds)
  auc_scores <- numeric(num_folds)
  for (i in 1:num_folds) {
    train_indices <- unlist(folds[-i])
    test_indices <- unlist(folds[i])
    train <- data[train_indices, ]
    test <- data[test_indices, ]
    form <- as.formula(paste(start_form_str," - ",var,sep=""))
    model <- glm(form, data = train, family = "binomial")
    predictions <- predict(model, newdata = test, type = "response")
    roc <- roc(response=test$infected,predictor=predictions,levels=c("No", "Yes"),direction = "<")
    auc_scores[i] <- auc(roc)
  }
  var_aucs$auc[var_aucs$var == var] <- mean(auc_scores)
} 
# Didn't remove anything

# Add to strat + cd420 + z30 + time + cd40 + trt + cd40*trt + treat*str2 - str2
set.seed(19)
vars <- names(data)
vars <- vars[vars!="infected"]
allVars <- vars
for (i in 1:length(vars)){
  for (j in 1:i){
    if(vars[i]!=vars[j]) {
      allVars <- c(allVars,paste(vars[i],'*',vars[j],sep=''))
    }
  }
}
allVars <- allVars[allVars!="strat"]
allVars <- allVars[allVars!="cd420"]
allVars <- allVars[allVars!="z30"]
allVars <- allVars[allVars!="time"]
allVars <- allVars[allVars!="cd40*trt"]
allVars <- allVars[allVars!="cd40"]
allVars <- allVars[allVars!="trt"]
allVars <- allVars[allVars!="treat*str2"]
allVars <- allVars[allVars!="treat"]
allVars <- allVars[allVars!="offtrt*preanti"]
allVars <- allVars[allVars!="offtrt"]
num_vars <- length(allVars)
var_aucs <- data.frame("vars" = allVars)
num_folds <- 10
for (j in 1:num_vars) {
  var <- allVars[j]
  print(paste(j,'/',num_vars,': ',var,sep=''))
  folds <- createFolds(data$infected, k = num_folds)
  auc_scores <- numeric(num_folds)
  for (i in 1:num_folds) {
    tryCatch({
      train_indices <- unlist(folds[-i])
      test_indices <- unlist(folds[i])
      train <- data[train_indices, ]
      test <- data[test_indices, ]
      form <- as.formula(paste("infected ~ strat + cd420 + z30 + time + cd40*trt + treat*str2 - str2 + offtrt*preanti - preanti + ",var,sep=""))
      model <- glm(form, data = train, family = "binomial")
      predictions <- predict(model, newdata = test, type = "response")
      roc <- roc(response=test$infected,predictor=predictions,levels=c("No", "Yes"),direction = "<")
      auc_scores[i] <- auc(roc)
    }, error = function(e) { 
      auc_scores <- auc_scores[-i]
    })
  }
  var_aucs$auc[var_aucs$var == var] <- mean(auc_scores)
}
# cd820*wtkg = 0.7021376

# Try removing a variable
set.seed(20)
start_form_str <- 'infected ~ strat + cd420 + z30 + time + cd40 + trt+ cd40*trt + treat*str2 + treat - str2 + offtrt*preanti + offtrt - preanti + cd820*wtkg + cd820 + wtkg'
vars <- unlist(strsplit(start_form_str, "[\\+\\-]"))
vars <- trimws(vars)
vars[1] <- substr(vars[1], 12, nchar(vars[1]))
num_vars <- length(vars)
var_aucs <- data.frame("vars" = vars)
num_folds <- 10
for (j in 1:num_vars) {
  var <- vars[j]
  print(var)
  folds <- createFolds(data$infected, k = num_folds)
  auc_scores <- numeric(num_folds)
  for (i in 1:num_folds) {
    train_indices <- unlist(folds[-i])
    test_indices <- unlist(folds[i])
    train <- data[train_indices, ]
    test <- data[test_indices, ]
    form <- as.formula(paste(start_form_str," - ",var,sep=""))
    model <- glm(form, data = train, family = "binomial")
    predictions <- predict(model, newdata = test, type = "response")
    roc <- roc(response=test$infected,predictor=predictions,levels=c("No", "Yes"),direction = "<")
    auc_scores[i] <- auc(roc)
  }
  var_aucs$auc[var_aucs$var == var] <- mean(auc_scores)
} 
# cd40 = 0.7021947

# Try removing a variable
set.seed(20)
start_form_str <- 'infected ~ strat + cd420 + z30 + time + trt+ cd40*trt - cd40 + treat*str2 + treat - str2 + offtrt*preanti + offtrt - preanti + cd820*wtkg + cd820 + wtkg'
vars <- unlist(strsplit(start_form_str, "[\\+\\-]"))
vars <- trimws(vars)
vars[1] <- substr(vars[1], 12, nchar(vars[1]))
num_vars <- length(vars)
var_aucs <- data.frame("vars" = vars)
num_folds <- 10
for (j in 1:num_vars) {
  var <- vars[j]
  print(var)
  folds <- createFolds(data$infected, k = num_folds)
  auc_scores <- numeric(num_folds)
  for (i in 1:num_folds) {
    train_indices <- unlist(folds[-i])
    test_indices <- unlist(folds[i])
    train <- data[train_indices, ]
    test <- data[test_indices, ]
    form <- as.formula(paste(start_form_str," - ",var,sep=""))
    model <- glm(form, data = train, family = "binomial")
    predictions <- predict(model, newdata = test, type = "response")
    roc <- roc(response=test$infected,predictor=predictions,levels=c("No", "Yes"),direction = "<")
    auc_scores[i] <- auc(roc)
  }
  var_aucs$auc[var_aucs$var == var] <- mean(auc_scores)
} 
# Didn't remove anything

# Add to strat + cd420 + z30 + time + cd40*trt - cd40 + treat*str2 - str2 + offtrt*preanti - preanti + cd820*wtkg
set.seed(21)
vars <- names(data)
vars <- vars[vars!="infected"]
allVars <- vars
for (i in 1:length(vars)){
  for (j in 1:i){
    if(vars[i]!=vars[j]) {
      allVars <- c(allVars,paste(vars[i],'*',vars[j],sep=''))
    }
  }
}
allVars <- allVars[allVars!="strat"]
allVars <- allVars[allVars!="cd420"]
allVars <- allVars[allVars!="z30"]
allVars <- allVars[allVars!="time"]
allVars <- allVars[allVars!="cd40*trt"]
allVars <- allVars[allVars!="trt"]
allVars <- allVars[allVars!="treat*str2"]
allVars <- allVars[allVars!="treat"]
allVars <- allVars[allVars!="offtrt*preanti"]
allVars <- allVars[allVars!="offtrt"]
allVars <- allVars[allVars!="cd820*wtkg"]
allVars <- allVars[allVars!="cd820"]
allVars <- allVars[allVars!="wtkg"]
num_vars <- length(allVars)
var_aucs <- data.frame("vars" = allVars)
num_folds <- 10
for (j in 1:num_vars) {
  var <- allVars[j]
  print(paste(j,'/',num_vars,': ',var,sep=''))
  folds <- createFolds(data$infected, k = num_folds)
  auc_scores <- numeric(num_folds)
  for (i in 1:num_folds) {
    tryCatch({
      train_indices <- unlist(folds[-i])
      test_indices <- unlist(folds[i])
      train <- data[train_indices, ]
      test <- data[test_indices, ]
      form <- as.formula(paste("infected ~ strat + cd420 + z30 + time + cd40*trt - cd40 + treat*str2 - str2 + offtrt*preanti - preanti + cd820*wtkg + ",var,sep=""))
      model <- glm(form, data = train, family = "binomial")
      predictions <- predict(model, newdata = test, type = "response")
      roc <- roc(response=test$infected,predictor=predictions,levels=c("No", "Yes"),direction = "<")
      auc_scores[i] <- auc(roc)
    }, error = function(e) { 
      auc_scores <- auc_scores[-i]
    })
  }
  var_aucs$auc[var_aucs$var == var] <- mean(auc_scores)
}
# homo*trt = 0.7024102

# Try removing a variable
set.seed(22)
start_form_str <- 'infected ~ strat + cd420 + z30 + time + trt+ cd40*trt - cd40 + treat*str2 + treat - str2 + offtrt*preanti + offtrt - preanti + cd820*wtkg + cd820 + wtkg + homo*trt + homo'
vars <- unlist(strsplit(start_form_str, "[\\+\\-]"))
vars <- trimws(vars)
vars[1] <- substr(vars[1], 12, nchar(vars[1]))
vars <- vars[vars!="cd40"]
vars <- vars[vars!="preanti"]
vars <- vars[vars!="str2"]
num_vars <- length(vars)
var_aucs <- data.frame("vars" = vars)
num_folds <- 10
for (j in 1:num_vars) {
  var <- vars[j]
  print(var)
  folds <- createFolds(data$infected, k = num_folds)
  auc_scores <- numeric(num_folds)
  for (i in 1:num_folds) {
    train_indices <- unlist(folds[-i])
    test_indices <- unlist(folds[i])
    train <- data[train_indices, ]
    test <- data[test_indices, ]
    form <- as.formula(paste(start_form_str," - ",var,sep=""))
    model <- glm(form, data = train, family = "binomial")
    predictions <- predict(model, newdata = test, type = "response")
    roc <- roc(response=test$infected,predictor=predictions,levels=c("No", "Yes"),direction = "<")
    auc_scores[i] <- auc(roc)
  }
  var_aucs$auc[var_aucs$var == var] <- mean(auc_scores)
} 
# Didn't remove anything

# Add to strat + cd420 + z30 + time + cd40*trt - cd40 + treat*str2 - str2 + offtrt*preanti - preanti + cd820*wtkg + homo*trt
set.seed(23)
vars <- names(data)
vars <- vars[vars!="infected"]
allVars <- vars
for (i in 1:length(vars)){
  for (j in 1:i){
    if(vars[i]!=vars[j]) {
      allVars <- c(allVars,paste(vars[i],'*',vars[j],sep=''))
    }
  }
}
allVars <- allVars[allVars!="strat"]
allVars <- allVars[allVars!="cd420"]
allVars <- allVars[allVars!="z30"]
allVars <- allVars[allVars!="time"]
allVars <- allVars[allVars!="cd40*trt"]
allVars <- allVars[allVars!="trt"]
allVars <- allVars[allVars!="treat*str2"]
allVars <- allVars[allVars!="treat"]
allVars <- allVars[allVars!="offtrt*preanti"]
allVars <- allVars[allVars!="offtrt"]
allVars <- allVars[allVars!="cd820*wtkg"]
allVars <- allVars[allVars!="cd820"]
allVars <- allVars[allVars!="wtkg"]
allVars <- allVars[allVars!="homo*trt"]
allVars <- allVars[allVars!="homo"]
num_vars <- length(allVars)
var_aucs <- data.frame("vars" = allVars)
num_folds <- 10
for (j in 1:num_vars) {
  var <- allVars[j]
  print(paste(j,'/',num_vars,': ',var,sep=''))
  folds <- createFolds(data$infected, k = num_folds)
  auc_scores <- numeric(num_folds)
  for (i in 1:num_folds) {
    tryCatch({
      train_indices <- unlist(folds[-i])
      test_indices <- unlist(folds[i])
      train <- data[train_indices, ]
      test <- data[test_indices, ]
      form <- as.formula(paste("infected ~ strat + cd420 + z30 + time + cd40*trt - cd40 + treat*str2 - str2 + offtrt*preanti - preanti + cd820*wtkg + homo*trt + ",var,sep=""))
      model <- glm(form, data = train, family = "binomial")
      predictions <- predict(model, newdata = test, type = "response")
      roc <- roc(response=test$infected,predictor=predictions,levels=c("No", "Yes"),direction = "<")
      auc_scores[i] <- auc(roc)
    }, error = function(e) { 
      auc_scores <- auc_scores[-i]
    })
  }
  var_aucs$auc[var_aucs$var == var] <- mean(auc_scores)
}
# symptom*oprior = 0.7025499

# Try removing a variable
set.seed(24)
start_form_str <- 'infected ~ strat + cd420 + z30 + time + trt+ cd40*trt - cd40 + treat*str2 + treat - str2 + offtrt*preanti + offtrt - preanti + cd820*wtkg + cd820 + wtkg + homo*trt + homo + symptom*oprior + symptom + oprior'
vars <- unlist(strsplit(start_form_str, "[\\+\\-]"))
vars <- trimws(vars)
vars[1] <- substr(vars[1], 12, nchar(vars[1]))
vars <- vars[vars!="cd40"]
vars <- vars[vars!="preanti"]
vars <- vars[vars!="str2"]
num_vars <- length(vars)
var_aucs <- data.frame("vars" = vars)
num_folds <- 10
for (j in 1:num_vars) {
  var <- vars[j]
  print(var)
  folds <- createFolds(data$infected, k = num_folds)
  auc_scores <- numeric(num_folds)
  for (i in 1:num_folds) {
    train_indices <- unlist(folds[-i])
    test_indices <- unlist(folds[i])
    train <- data[train_indices, ]
    test <- data[test_indices, ]
    form <- as.formula(paste(start_form_str," - ",var,sep=""))
    model <- glm(form, data = train, family = "binomial")
    predictions <- predict(model, newdata = test, type = "response")
    roc <- roc(response=test$infected,predictor=predictions,levels=c("No", "Yes"),direction = "<")
    auc_scores[i] <- auc(roc)
  }
  var_aucs$auc[var_aucs$var == var] <- mean(auc_scores)
} 
# Didn't remove anything

# Add to strat + cd420 + z30 + time + cd40*trt - cd40 + treat*str2 - str2 + offtrt*preanti - preanti + cd820*wtkg + homo*trt + symptom*oprior
set.seed(25)
vars <- names(data)
vars <- vars[vars!="infected"]
allVars <- vars
for (i in 1:length(vars)){
  for (j in 1:i){
    if(vars[i]!=vars[j]) {
      allVars <- c(allVars,paste(vars[i],'*',vars[j],sep=''))
    }
  }
}
allVars <- allVars[allVars!="strat"]
allVars <- allVars[allVars!="cd420"]
allVars <- allVars[allVars!="z30"]
allVars <- allVars[allVars!="time"]
allVars <- allVars[allVars!="cd40*trt"]
allVars <- allVars[allVars!="trt"]
allVars <- allVars[allVars!="treat*str2"]
allVars <- allVars[allVars!="treat"]
allVars <- allVars[allVars!="offtrt*preanti"]
allVars <- allVars[allVars!="offtrt"]
allVars <- allVars[allVars!="cd820*wtkg"]
allVars <- allVars[allVars!="cd820"]
allVars <- allVars[allVars!="wtkg"]
allVars <- allVars[allVars!="homo*trt"]
allVars <- allVars[allVars!="homo"]
allVars <- allVars[allVars!="symptom*oprior"]
allVars <- allVars[allVars!="symptom"]
allVars <- allVars[allVars!="oprior"]
num_vars <- length(allVars)
var_aucs <- data.frame("vars" = allVars)
num_folds <- 10
for (j in 1:num_vars) {
  var <- allVars[j]
  print(paste(j,'/',num_vars,': ',var,sep=''))
  folds <- createFolds(data$infected, k = num_folds)
  auc_scores <- numeric(num_folds)
  for (i in 1:num_folds) {
    tryCatch({
      train_indices <- unlist(folds[-i])
      test_indices <- unlist(folds[i])
      train <- data[train_indices, ]
      test <- data[test_indices, ]
      form <- as.formula(paste("infected ~ strat + cd420 + z30 + time + cd40*trt - cd40 + treat*str2 - str2 + offtrt*preanti - preanti + cd820*wtkg + homo*trt + symptom*oprior + ",var,sep=""))
      model <- glm(form, data = train, family = "binomial")
      predictions <- predict(model, newdata = test, type = "response")
      roc <- roc(response=test$infected,predictor=predictions,levels=c("No", "Yes"),direction = "<")
      auc_scores[i] <- auc(roc)
    }, error = function(e) { 
      auc_scores <- auc_scores[-i]
    })
  }
  var_aucs$auc[var_aucs$var == var] <- mean(auc_scores)
}
# cd420*treat = 0.7026381

# Try removing a variable
set.seed(26)
start_form_str <- 'infected ~ strat + cd420 + z30 + time + trt+ cd40*trt - cd40 + treat*str2 + treat - str2 + offtrt*preanti + offtrt - preanti + cd820*wtkg + cd820 + wtkg + homo*trt + homo + symptom*oprior + symptom + oprior + cd420*treat'
vars <- unlist(strsplit(start_form_str, "[\\+\\-]"))
vars <- trimws(vars)
vars[1] <- substr(vars[1], 12, nchar(vars[1]))
vars <- vars[vars!="cd40"]
vars <- vars[vars!="preanti"]
vars <- vars[vars!="str2"]
num_vars <- length(vars)
var_aucs <- data.frame("vars" = vars)
num_folds <- 10
for (j in 1:num_vars) {
  var <- vars[j]
  print(var)
  folds <- createFolds(data$infected, k = num_folds)
  auc_scores <- numeric(num_folds)
  for (i in 1:num_folds) {
    train_indices <- unlist(folds[-i])
    test_indices <- unlist(folds[i])
    train <- data[train_indices, ]
    test <- data[test_indices, ]
    form <- as.formula(paste(start_form_str," - ",var,sep=""))
    model <- glm(form, data = train, family = "binomial")
    predictions <- predict(model, newdata = test, type = "response")
    roc <- roc(response=test$infected,predictor=predictions,levels=c("No", "Yes"),direction = "<")
    auc_scores[i] <- auc(roc)
  }
  var_aucs$auc[var_aucs$var == var] <- mean(auc_scores)
} 
# Didn't remove anything

# Try to add a variable
set.seed(27)
vars <- names(data)
vars <- vars[vars!="infected"]
allVars <- vars
for (i in 1:length(vars)){
  for (j in 1:i){
    if(vars[i]!=vars[j]) {
      allVars <- c(allVars,paste(vars[i],'*',vars[j],sep=''))
    }
  }
}
allVars <- allVars[allVars!="strat"]
allVars <- allVars[allVars!="cd420"]
allVars <- allVars[allVars!="z30"]
allVars <- allVars[allVars!="time"]
allVars <- allVars[allVars!="cd40*trt"]
allVars <- allVars[allVars!="trt"]
allVars <- allVars[allVars!="treat*str2"]
allVars <- allVars[allVars!="treat"]
allVars <- allVars[allVars!="offtrt*preanti"]
allVars <- allVars[allVars!="offtrt"]
allVars <- allVars[allVars!="cd820*wtkg"]
allVars <- allVars[allVars!="cd820"]
allVars <- allVars[allVars!="wtkg"]
allVars <- allVars[allVars!="homo*trt"]
allVars <- allVars[allVars!="homo"]
allVars <- allVars[allVars!="symptom*oprior"]
allVars <- allVars[allVars!="symptom"]
allVars <- allVars[allVars!="oprior"]
allVars <- allVars[allVars!="cd420*treat"]
num_vars <- length(allVars)
var_aucs <- data.frame("vars" = allVars)
num_folds <- 10
for (j in 1:num_vars) {
  var <- allVars[j]
  print(paste(j,'/',num_vars,': ',var,sep=''))
  folds <- createFolds(data$infected, k = num_folds)
  auc_scores <- numeric(num_folds)
  for (i in 1:num_folds) {
    tryCatch({
      train_indices <- unlist(folds[-i])
      test_indices <- unlist(folds[i])
      train <- data[train_indices, ]
      test <- data[test_indices, ]
      form <- as.formula(paste("infected ~ strat + cd420 + z30 + time + cd40*trt - cd40 + treat*str2 - str2 + offtrt*preanti - preanti + cd820*wtkg + homo*trt + symptom*oprior + cd420*treat + ",var,sep=""))
      model <- glm(form, data = train, family = "binomial")
      predictions <- predict(model, newdata = test, type = "response")
      roc <- roc(response=test$infected,predictor=predictions,levels=c("No", "Yes"),direction = "<")
      auc_scores[i] <- auc(roc)
    }, error = function(e) { 
      auc_scores <- auc_scores[-i]
    })
  }
  var_aucs$auc[var_aucs$var == var] <- mean(auc_scores)
}
# offtrt*gender = 0.7028418

# Try removing a variable
set.seed(28)
start_form_str <- 'infected ~ strat + cd420 + z30 + time + trt+ cd40:trt - cd40 + treat:str2 + treat - str2 + offtrt:preanti + offtrt - preanti + cd820:wtkg + cd820 + wtkg + homo:trt + homo + symptom:oprior + symptom + oprior + cd420:treat + offtrt:gender + gender'
vars <- unlist(strsplit(start_form_str, "[\\+\\-]"))
vars <- trimws(vars)
vars[1] <- substr(vars[1], 12, nchar(vars[1]))
vars <- vars[vars!="cd40"]
vars <- vars[vars!="preanti"]
vars <- vars[vars!="str2"]
num_vars <- length(vars)
var_aucs <- data.frame("vars" = vars)
num_folds <- 10
for (j in 1:num_vars) {
  var <- vars[j]
  print(var)
  folds <- createFolds(data$infected, k = num_folds)
  auc_scores <- numeric(num_folds)
  for (i in 1:num_folds) {
    train_indices <- unlist(folds[-i])
    test_indices <- unlist(folds[i])
    train <- data[train_indices, ]
    test <- data[test_indices, ]
    form <- as.formula(paste(start_form_str," - ",var,sep=""))
    model <- glm(form, data = train, family = "binomial")
    predictions <- predict(model, newdata = test, type = "response")
    roc <- roc(response=test$infected,predictor=predictions,levels=c("No", "Yes"),direction = "<")
    auc_scores[i] <- auc(roc)
  }
  var_aucs$auc[var_aucs$var == var] <- mean(auc_scores)
} 
# Didn't remove anything

model <- glm('infected ~ strat + cd420 + strat*cd420 - strat*cd420', data = data, family = "binomial")
summary(model)

# Try to add a variable
set.seed(29)
vars <- names(data)
vars <- vars[vars!="infected"]
allVars <- vars
for (i in 1:length(vars)){
  for (j in 1:i){
    if(vars[i]!=vars[j]) {
      allVars <- c(allVars,paste(vars[i],'*',vars[j],sep=''))
    }
  }
}
allVars <- allVars[allVars!="strat"]
allVars <- allVars[allVars!="cd420"]
allVars <- allVars[allVars!="z30"]
allVars <- allVars[allVars!="time"]
allVars <- allVars[allVars!="cd40*trt"]
allVars <- allVars[allVars!="trt"]
allVars <- allVars[allVars!="treat*str2"]
allVars <- allVars[allVars!="treat"]
allVars <- allVars[allVars!="offtrt*preanti"]
allVars <- allVars[allVars!="offtrt"]
allVars <- allVars[allVars!="cd820*wtkg"]
allVars <- allVars[allVars!="cd820"]
allVars <- allVars[allVars!="wtkg"]
allVars <- allVars[allVars!="homo*trt"]
allVars <- allVars[allVars!="homo"]
allVars <- allVars[allVars!="symptom*oprior"]
allVars <- allVars[allVars!="symptom"]
allVars <- allVars[allVars!="oprior"]
allVars <- allVars[allVars!="cd420*treat"]
allVars <- allVars[allVars!="offtrt*gender"]
allVars <- allVars[allVars!="gender"]
num_vars <- length(allVars)
var_aucs <- data.frame("vars" = allVars)
num_folds <- 10
for (j in 1:num_vars) {
  var <- allVars[j]
  print(paste(j,'/',num_vars,': ',var,sep=''))
  folds <- createFolds(data$infected, k = num_folds)
  auc_scores <- numeric(num_folds)
  for (i in 1:num_folds) {
    tryCatch({
      train_indices <- unlist(folds[-i])
      test_indices <- unlist(folds[i])
      train <- data[train_indices, ]
      test <- data[test_indices, ]
      form <- as.formula(paste("infected ~ strat + cd420 + z30 + time + cd40:trt + trt + treat:str2 + treat + offtrt:preanti + offtrt + cd820*wtkg + homo*trt + symptom*oprior + cd420*treat + offtrt*gender + ",var,sep=""))
      model <- glm(form, data = train, family = "binomial")
      predictions <- predict(model, newdata = test, type = "response")
      roc <- roc(response=test$infected,predictor=predictions,levels=c("No", "Yes"),direction = "<")
      auc_scores[i] <- auc(roc)
    }, error = function(e) { 
      auc_scores <- auc_scores[-i]
    })
  }
  var_aucs$auc[var_aucs$var == var] <- mean(auc_scores)
}
# Added no more variables

model <- glm('infected ~ strat + cd420 + z30 + time + cd40:trt + trt + treat:str2 + treat + offtrt:preanti + offtrt + cd820*wtkg + homo*trt + symptom*oprior + cd420*treat + offtrt*gender', data = data, family = "binomial")
summary(model)
