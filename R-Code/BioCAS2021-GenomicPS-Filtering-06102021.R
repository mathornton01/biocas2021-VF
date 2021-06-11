library(YinGenomicDFTDistances);
library(dplyr);

evenlyScaleEnsemble <- function(spectraList){
  spectraLengths <- lapply(spectraList,length); 
  maxLength <- max(unlist(spectraLengths));
  maxIndex <- which(unlist(lapply(spectraList,length)) == max(unlist(lapply(spectraList,length))));
  scaledSpectra <- list(length(spectraList)); 
  scaledSpectra[[maxIndex]] <- spectraList[[maxIndex]];
  for (i in 1:length(spectraList)){
    if (i == maxIndex){
      next;
    }
    scaledSpectra[[i]] <- evenlyScaleSingle(spectraList[[i]], maxLength); 
  }
  return(scaledSpectra); 
}

sarscvaug$Sequence %>% 
  encodeGenomes() %>% 
  getPowerSpectraEnsemble() %>% 
  evenlyScaleEnsemble() -> sarscvaugPSList;

viralData <- sarscvaug;

viralData$Location <- unlist(lapply(strsplit(as.character(viralData$Header),"/"), function(x){x[2]}));
for (inc in c('tiger','env')){
  viralData[viralData$Location == inc, ]$Location <- unlist(lapply(strsplit(as.character(viralData[viralData$Location == inc,]$Header),"/"), function(x){x[3]}))
}
viralData$Region <- rep("blank",nrow(viralData)); 
if (nrow(viralData[viralData$Location %in% c("Algeria", "Egypt", "South Africa",  "DRC",  "Gambia", "Senegal"),]) > 0){
  viralData[viralData$Location %in% c("Algeria", "Egypt", "South Africa",  "DRC",  "Gambia", "Senegal"),]$Region <- "Africa"
}
if (nrow(viralData[viralData$Location %in% c("Beijing" ,"Chongqing" ,"Fujian" ,"Guangdong" ,"Hangzhou" , "Hong Kong" ,"Jiangsu" ,"Jiangxi" ,"Jingzhou" ,"Shandong" ,"Shenzhen" , "Sichuan" , "Taiwan", "Tianmen" ,"Wuhan", "Yunnan" , "Zhejiang" ,"Lishui"  ,"Japan"  ,"Guangzhou"),]) > 0){
  viralData[viralData$Location %in% c("Beijing" ,"Chongqing" ,"Fujian" ,"Guangdong" ,"Hangzhou" , "Hong Kong" ,"Jiangsu" ,"Jiangxi" ,"Jingzhou" ,"Shandong" ,"Shenzhen" , "Sichuan" , "Taiwan", "Tianmen" ,"Wuhan", "Yunnan" , "Zhejiang" ,"Lishui"  ,"Japan"  ,"Guangzhou"),]$Region <- "East Asia"
}
if (nrow(viralData[viralData$Location %in% c( "Brazil" ,"Chile" ,"Colombia" ,"Costa Rica" ,"Uruguay"),]) > 0){
  viralData[viralData$Location %in% c( "Brazil" ,"Chile" ,"Colombia" ,"Costa Rica" ,"Uruguay"),]$Region <- "South America"
}
if (nrow(viralData[viralData$Location %in% c("Austria", "Belgium", "Czech Republic", "Denmark", "Finland", "France" , "Georgia", "Germany" , "Greece" , "Spain", "Sweden" , "Hungary" , "Portugal" , "Poland", "Russia", "Romania" , "Slovakia" , "Italy"),]) > 0){
  viralData[viralData$Location %in% c("Austria", "Belgium", "Czech Republic", "Denmark", "Finland", "France" , "Georgia", "Germany" , "Greece" , "Spain", "Sweden" , "Hungary" , "Portugal" , "Poland", "Russia", "Romania" , "Slovakia" , "Italy"),]$Region <- "Europe"
}
if (nrow(viralData[viralData$Location %in% c( "Bangladesh" , "Cambodia" , "India" , "Nepal" , "Sri Lanka" , "Vietnam" , "Thailand"),]) > 0){
  viralData[viralData$Location %in% c( "Bangladesh" , "Cambodia" , "India" , "Nepal" , "Sri Lanka" , "Vietnam" , "Thailand"),]$Region <- "West Asia"
}
if (nrow(viralData[viralData$Location %in% c( "Australia", "New Zealand" , "Indonesia" , "Malaysia"),]) > 0){
  viralData[viralData$Location %in% c( "Australia", "New Zealand" , "Indonesia" , "Malaysia"),]$Region <- "Oceania"
}
if (nrow(viralData[viralData$Location %in% c("Turkey" , "Saudi Arabia" , "Kazakhstan" ,"Iran", "Israel" , "Kuwait"),]) > 0){
  viralData[viralData$Location %in% c("Turkey" , "Saudi Arabia" , "Kazakhstan" ,"Iran", "Israel" , "Kuwait"),]$Region <- "Middle East"
}
if (nrow(viralData[viralData$Location %in% c("USA", "Puerto Rico", "Guam", "Mexico"),]) > 0){
  viralData[viralData$Location %in% c("USA", "Puerto Rico", "Guam", "Mexico"),]$Region <- "North America"
}
if (nrow(viralData[viralData$Location %in% c("blank"),]) > 0){
  viralData[viralData$Location %in% c("blank"),]$Region <- viralData[viralData$Location %in% c("blank"),]$Location
}


matrix(unlist(sarscvaugPSList),
       nrow=length(sarscvaugPSList),
       byrow=TRUE) -> sarscvaugPSMatrix;

colMeans(sarscvaugPSMatrix) -> sarscvaugPSMeans;
apply(sarscvaugPSMatrix,2,var) -> sarscvaugPSVars;
sarscvaugPSMax <- apply(sarscvaugPSMatrix,2,max);

reg <- viralData$Region
regnum <- as.numeric(reg == "East Asia");
regnum[reg == "Europe"] <- 2; 
regnum[reg == "West Asia"] <- 3; 
regnum[reg == "South America"] <- 4; 
regnum[reg == "North America"] <- 5; 
regnum[reg == "Oceania"] <- 6; 
regnum[reg == "Middle East"] <- 7; 
regnum[reg == "Africa"] <- 8; 

par(mfrow=c(1,2))
image(t(sarscvaugPSMatrix[,order(sarscvaugPSVars)[1:1000]]/sarscvaugPSMax[order(sarscvaugPSVars)[1:1000]]),useRaster=TRUE,oldstyle=TRUE)
image(t(sarscvaugPSMatrix[,order(sarscvaugPSVars)[29000:30000]]/sarscvaugPSMax[order(sarscvaugPSVars)[29000:30000]]),useRaster=TRUE,oldstyle=TRUE)


library(keras)
library(tensorflow)
model <- keras_model_sequential() %>% 
  layer_conv_1d(filters=1235, 
                kernel_size=ncol(sarscvaugPSMatrix), 
                activation = "relu", 
                input_shape = c(ncol(sarscvaugPSMatrix),1)) %>% 
  layer_global_max_pooling_1d() %>% 
  layer_dense(80) %>% 
  layer_dropout(0.2) %>%   
  layer_dense(40) %>% 
  layer_dropout(0.2) %>% 
  layer_dense(20) %>% 
  layer_dropout(0.2) %>% 
  layer_activation("relu") %>% 
  layer_dense(8) %>% 
  layer_activation("softmax"); 


set.seed(0xBEEF)
group <- rep(1:5,nrow(viralData))
group <- sample(group[1:nrow(viralData)])

model %>% compile(loss="categorical_crossentropy", 
                  optimizer="adam",
                  metrics="accuracy")

shaped_data <- array_reshape(sarscvaugPSMatrix,c(dim(sarscvaugPSMatrix),1),order="C")

model %>% fit(shaped_data, to_categorical(matrix(as.numeric(as.factor(viralData$Region))-1)),
              epochs=5,
              batch_size=150) -> fitted_model

get_weights(model) ->cnvweights
cnnfilteredvars <- apply((sarscvaugPSMatrix %*% cnvweights[[1]][,,1:128]),2,var)
image((sarscvaugPSMatrix %*% cnvweights[[1]][,,1:128])[,order(cnnfilteredvars,decreasing=TRUE)[1:20]])
image(scaledPSMatrix[,order(sarscvaugPSVars,decreasing=TRUE)])


image(cnvweights[[1]][,,1:128])
image(t((sarscvaugPSMatrix %*% cnvweights[[1]][,,1:128])))
image(t(scale((sarscvaugPSMatrix %*% cnvweights[[1]][,,1:128])[order(regnum),])))
par(mfrow=c(1,2))
image(t((sarscvaugPSMatrix %*% cnvweights[[1]][,,1:128])[order(regnum),]))
image(sarscvaugPSMatrix[order(regnum),])

varfilt <- matrix(0,nrow=10,ncol=ncol(sarscvaugPSVars))
varfilt[1,which(order(sarscvaugPSVars,decreasing = TRUE) %in% 1:20000)] <- 1
varfilt[2,which(order(sarscvaugPSVars,decreasing = TRUE) %in% 1:10000)] <- 1
varfilt[3,which(order(sarscvaugPSVars,decreasing = TRUE) %in% 1:8000)] <- 1
varfilt[4,which(order(sarscvaugPSVars,decreasing = TRUE) %in% 1:5000)] <- 1
varfilt[5,which(order(sarscvaugPSVars,decreasing = TRUE) %in% 1:2500)] <- 1
varfilt[6,which(order(sarscvaugPSVars,decreasing = TRUE) %in% 1:1000)] <- 1
varfilt[7,which(order(sarscvaugPSVars,decreasing = TRUE) %in% 1:500)] <- 1
varfilt[8,which(order(sarscvaugPSVars,decreasing = TRUE) %in% 1:250)] <- 1
varfilt[9,which(order(sarscvaugPSVars,decreasing = TRUE) %in% 1:125)] <- 1
varfilt[10,which(order(sarscvaugPSVars,decreasing = TRUE) %in% 1:50)] <- 1

image(varfilt)

varMaxCoeffs <- order(sarscvaugPSVars,decreasing=TRUE)[1:1396]; 
PSPCA <- princomp(sarscvaugPSMatrix[,varMaxCoeffs])

image(scaledPSMatrix[,varMaxCoeffs])
image(PSPCA$loadings)

pdistFull <- dist(sarscvaugPSMatrix)

compmat <- sarscvaugPSMatrix[,varMaxCoeffs] %*% PSPCA$loadings

pdistComp100 <- dist(compmat[,1:100])
pdistComp200 <- dist(compmat[,1:200])
pdistComp500 <- dist(compmat[,1:500])
pdistCompFull <- dist(compmat)


knots <- seq.int(5,1235,5)
corPCA2 <- c(); 
for (i in 1235:ncol(compmat)){
  pdistComp <- dist(compmat[,1:i]);
  corPCA2 <- c(corPCA2,cor(pdistComp,pdistFull))
}

plot(corPCA,type='l')
corVarFilt <- c();
for(i in knots[199:247]){
  pdistVarFilt <- dist(sarscvaugPSMatrix[,order(sarscvaugPSVars,decreasing=TRUE)[1:i]]);
  corVarFilt <- c(corVarFilt,cor(pdistFull,pdistVarFilt))
  cat(i/1235);
}

corVarCNN <- c(); 
for(i in knots){
  pdistVarFilt <- dist(sarscvaugPSMatrix %*% cnvweights[[1]][,,1:i]);
  corVarCNN <- c(corVarCNN,cor(pdistFull,pdistVarFilt))
  cat(i/1235);
}

corKnots <- corPCA[knots]

setwd("C:/Users/s182954/Desktop/Paper-Submissions/BioCAS-2021/Variance-Filtering/Manuscript-06022021/Images/Files")
png("PSFiltMethods.PNG",height=480,width=960)
plot(knots,corKnots,type='l',ylim=c(0.4,1),
     xlab="Filtered PS Coefficient Size", 
     ylab="Correlation to Unfilitered PS", 
     main="Comparison of Three Filtering Methods \n by Correlation to Full PS Distances",
     cex.axis=2,
     cex.main=2,
     cex.lab=1.5,lwd=2)
lines(knots,corVarFilt,type='l',col='blue',lwd=2)
lines(knots,corVarCNN,type='l',col='green',lwd=2)
legend("bottomright",legend=c("Maximal Variance Principal Components","Minimal Variance Filtering","CNN Regional Differentiator Filters"), col=c("black","blue","green"),lty=1)
dev.off()

correlBioCASDF <- data.frame(knots=knots,VF = corVarFilt, CNN=corVarCNN, PCA = corKnots)
save(correlBioCASDF,file="BioCAScor.RData")
getwd()

plot(sarscvaugPSVars,ylim=c(0,500000000),cex.axis=2)
hist(sarscvaugPSVars,breaks=100000,xlim=c(0,100000000),cex.axis=2)

tedf <- ecdf(sarscvaugPSVars)
plot(seq.int(0,100000000,1000),1-tedf(seq.int(0,100000000,1000)),type='l',cex.axis=2,lwd=2)


corVFLP <- c() 
for (i in seq.int(1,ncol(sarscvaugPSMatrix),2000)){
  pdistVarFilt <- dist(sarscvaugPSMatrix[,order(sarscvaugPSVars,decreasing=TRUE)[1:i]]);
  corVFLP <- c(corVFLP,cor(pdistFull,pdistVarFilt))
  cat(i/30563);
}

plot(c(knots,seq.int(1,ncol(sarscvaugPSMatrix),2000)[2:length(corVFLP)]),c(corVarFilt,corVFLP[2:length(corVFLP)]),type='l',
     lwd=2,cex.axis=2)

biplot(PSPCA,cex.axis=2)
screeplot(PSPCA,cex.axis=2,cex.lab=2,cex=2)

Labs <- viralData$Region

# Maximal Variance Filtering (MVF) 
MVF50 <- cbind(Labs,sarscvaugPSMatrix[,order(sarscvaugPSVars,decreasing=TRUE)[1:50]])
MVF100 <- cbind(Labs,sarscvaugPSMatrix[,order(sarscvaugPSVars,decreasing=TRUE)[1:100]])
MVF250 <- cbind(Labs,sarscvaugPSMatrix[,order(sarscvaugPSVars,decreasing=TRUE)[1:250]])
MVF500 <- cbind(Labs,sarscvaugPSMatrix[,order(sarscvaugPSVars,decreasing=TRUE)[1:500]])
MVF1000 <- cbind(Labs,sarscvaugPSMatrix[,order(sarscvaugPSVars,decreasing=TRUE)[1:1000]])

# Automatic Filter Learning (AFL)
AFL50 <- cbind(Labs,(sarscvaugPSMatrix %*% cnvweights[[1]][,,1:50]))
AFL100 <- cbind(Labs,(sarscvaugPSMatrix %*% cnvweights[[1]][,,1:100]))
AFL250 <- cbind(Labs,(sarscvaugPSMatrix %*% cnvweights[[1]][,,1:250]))
AFL500 <- cbind(Labs,(sarscvaugPSMatrix %*% cnvweights[[1]][,,1:500]))
AFL1000 <- cbind(Labs,(sarscvaugPSMatrix %*% cnvweights[[1]][,,1:1000]))

# Maximal Variance Principal Components Filter (MVPCF)
MVPCF50 <- cbind(Labs,(compmat[,1:50]))
MVPCF100 <- cbind(Labs,(compmat[,1:100]))
MVPCF250 <- cbind(Labs,(compmat[,1:250]))
MVPCF500 <- cbind(Labs,(compmat[,1:500]))
MVPCF1000 <- cbind(Labs,(compmat[,1:1000]))


# Random Forest Models with five fold cross validation: 
library(randomForest)

dflist <- list(MVF50,MVF100,MVF250,MVF500,MVF1000,
               AFL50,AFL100,AFL250,AFL500,AFL1000,
               MVPCF50,MVPCF100,MVPCF250,MVPCF500,MVPCF1000);

set.seed(0xBEEF)
group <- rep(1:5,nrow(viralData))
group <- sample(group[1:nrow(viralData)])

lapply(dflist, function(x) {
  return(mean(unlist(lapply(1:5, function(y){
  rf <- randomForest(factor(Labs,,levels=unique(x[,'Labs'])) ~ .,data=x[group != y,])
  return(sum(diag(table(factor(x[group==y,1],levels=unique(x[,'Labs'])),predict(rf,newdata=x[group==y,]))))/sum(table(factor(x[group==y,1],levels=unique(x[,'Labs'])),predict(rf,newdata=x[group==y,]))))
  }))))
}) -> acclistlong

