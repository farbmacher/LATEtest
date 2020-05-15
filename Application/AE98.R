#################################################
# Application in
# "Instrument Validity Tests with Causal Forests"
# Helmut Farbmacher et al 2020
# May 2020
#################################################

rm(list=ls())
packages <- c("LATEtest","readstata13", "dummies")
invisible(lapply(packages, library, character.only=T))

system <- Sys.info()['sysname']
if(system == 'Windows') {
  setwd("C:/Users/farbma/Documents/GitHub/LATEtest")    # my path Windows
} else if(system == 'Darwin') {
  setwd("/Users/helmut/Documents/GitHub/LATEtest")      # my path Mac
}

set.seed(10101)
setwd("./Application")
data<-read.dta13("AE98_data.dta")
nrow(data)

subsets <- 4
synthetic <- FALSE     # adds synthetic direct effects to the outcome (empirical Monte-Carlo simulation)

Xvars <- paste0(c("age","agefst","othrace","white","black","hisp","lesshs","hs","morehs"))
data$Y <- data$lincome
data$D <- data$morekids
data$Z <- data$samesex

if(synthetic==TRUE) {
  subset <- data$age<25 & data$hs==1
  table(data$hs,data$age)
  table(subset(data.frame(data$hs,data$age),subset))
  table(subset(data.frame(data$Z,data$D),subset))
  print("mean(data$Y)")
  print(mean(data$Y))
  print("sd(data$Y)")
  print(sd(data$Y))
  print("length(subset(data$Y,subset))")
  print(length(subset(data$Y,subset)))
  print("length(subset(data$Y,subset))/nrow(data)")
  print(length(subset(data$Y,subset))/nrow(data))
  print("mean(subset(data$Y,(subset)))")
  print(mean(subset(data$Y,(subset))))
  effect_size <- 1/4*sd(data$Y)
  print("effect_size")
  print(effect_size)
  data$Y <- data$Y + effect_size*subset*data$Z
  print("mean(subset(data$Y,subset))")
  print(mean(subset(data$Y,subset)))
  print("mean(data$Y)")
  print(mean(data$Y))
}

start_time = Sys.time()
test <- LATEtest(data=data, covars=Xvars, subsets=subsets, huge=TRUE, tree_fraction=0.05, alpha=0.05, minsize=800)
end_time = Sys.time()
maxtree <- eval(parse(text=paste("test$treelist$tree_",test$maxTtree$label,test$maxTtree$J,sep="")))
rpart.plot(maxtree,extra = 101,box.palette="GyRd",shadow.col="gray",nn=TRUE,roundint=FALSE)

sink(file=paste("AE98_FGK_synth_",synthetic,".txt",sep=""), split=TRUE)
print(Sys.time())
end_time - start_time
print("test$leafinfo")
test$leafinfo
print("test$results")
test$results
sink(file=NULL)
