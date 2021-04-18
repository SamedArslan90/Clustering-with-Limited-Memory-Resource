library(farff)
library(dplyr)
library(cluster)
library(clValid) 
library(funtimes )
library(aricode )
library(mcclust )
library(caret)
library(clue)
set.seed(1)
options(scipen = 999)

filelist <- list.files(".\\datasets\\artificial")
catalog <- data.frame()

# Creation of data catalog
for (name in filelist) {


data = readARFF(paste(".\\datasets\\artificial\\",name,sep = ""))


catalog <- rbind(catalog,data.frame(line = nrow(catalog)+1, name= name,rnumer = nrow(data), cnumber= ncol(data),lastcolname= tail(colnames(data),1),nclass=length(unique(data[,ncol(data)]))  ))

}





## matches the fitted classes returns real classes


classmatcher <- function(target,pred){
  
  
  confmatrix <- table(target,pred)
  
  confmatrix <- cbind(confmatrix,matrix(0,ncol=nrow(confmatrix)-ncol(confmatrix)  ,nrow = nrow(confmatrix) ))
  
  
  df.res = solve_LSAP(confmatrix, maximum = T)
  
  map <- cbind(rownames(confmatrix), colnames(confmatrix)[df.res])
  colnames(map) <- c("target","pred")
  
  res <- left_join(data.frame(pred = as.character(pred)),map,copy=T)
  
  factor(as.character(res$target))
}


## calculates nearest cluster 
findclus <- function(x) { x %>% as.numeric() %>% sweep(cent[, -1], 1, . , "-", check.margin = FALSE) %>% .^2 %>% rowSums() %>% which.min() %>% cent[.,1] }

performancelist <- data.frame()

#########################  

for(line in 1:nrow(catalog)){


data = readARFF(paste(".\\datasets\\artificial\\",catalog[line,"name"],sep = ""))

## data shuffling two times
data <- data[sample(nrow(data)),]
data <- data[sample(nrow(data)),]
data <- data[sample(nrow(data)),]
data <- data[sample(nrow(data)),]

rownames(data) <- NULL



databody <- data[,-catalog[line,"cnumber"]] 

datacluster <- factor(data[,catalog[line,"cnumber"]])


## Basic clustering

hc <- stats::hclust(dist(databody), method = "centroid")
hc.class <- cutree(hc, k = catalog[line,"nclass"])
hc.real.class <- classmatcher(datacluster,hc.class)
all_class <- union(datacluster, hc.real.class)
hc.real.class.table <- table(factor(datacluster, all_class), factor(hc.real.class, all_class))





# Clustering with SMOTE

datasmote <- data[,-catalog[line,"cnumber"]]

lastnrow <- nrow(datasmote) + 1
group <- data.frame()
grouplength <- 100
# Chosing big clusters to smote and defining little groups with a limit
minclasslimit <- round(grouplength/ifelse(catalog[line,"nclass"]<4,4,catalog[line,"nclass"]))


while (nrow(datasmote) > grouplength & nrow(datasmote)/lastnrow < 1.1 & sum(is.na(datasmote[1:grouplength,])) == 0) {
  
  group <- datasmote[1:grouplength,]
  
  
  hc <- stats::hclust(dist(group), method = "centroid")
  group$tag <- cutree(hc, k = catalog[line,"nclass"])
  
  
  smotedata <- group[group$tag %in% which(table(group$tag) > minclasslimit), ]
  
  smotedata$tag <- as.factor(smotedata$tag)
  # table(smotedata$tag)
  if (length(unique(smotedata$tag)) > 1) {
    
    tryCatch({
    aftersmote <- DMwR::SMOTE(tag ~ ., smotedata, perc.over = 25)} , error = function(e) aftersmote <- smotedata)
    colnames(aftersmote) <- colnames(group)
    
    # table(aftersmote$tag)
    smoteoutput <- rbind(aftersmote[,-length(smotedata)], group[group$tag %in% which(table(group$tag) < (minclasslimit+1)),-length(smotedata)])
    
  } else {   
    
    smoteoutput <- rbind(head(group[group$tag %in% which(table(group$tag) > minclasslimit) ,-length(smotedata)],round(grouplength/5))
                         , group[group$tag %in% which(table(group$tag) < (minclasslimit+1)) ,-length(smotedata)])
    
  }
  
  smoteoutput <- as.data.frame(na.omit(smoteoutput))
  lastnrow <- nrow(datasmote)
  
  datasmote <- data.frame(rbind(datasmote[(grouplength+1):nrow(datasmote),],smoteoutput))
  
  rownames(datasmote) <- NULL
}


############ Last clustering

hc <- stats::hclust(dist(datasmote), method = "centroid")
datasmote$tag <- cutree(hc, k = catalog[line,"nclass"])

#detect cluster centers

cent <- aggregate(. ~ tag, datasmote, mean)

#finds nearest cluster center
hc.smote.class <-  as.vector(apply(databody, 1, findclus ))
hc.smote.real.class <- classmatcher(datacluster,hc.smote.class)
all_class <- union(datacluster, hc.smote.real.class)
hc.smote.real.class.table <- table(factor(datacluster, all_class), factor(hc.smote.real.class, all_class))



results <- cbind(databody, target = datacluster, pred = hc.smote.class)




# Performances

CM.hc <- caret::confusionMatrix(hc.real.class.table)
CM.smote.hc <- caret::confusionMatrix(hc.smote.real.class.table)


#cent.real <- aggregate(databody,data.frame(datacluster), mean)
#cent.hc <- aggregate(databody,data.frame(hc.real.class), mean)
#cent.hc.smote <- aggregate(databody,data.frame(hc.smote.real.class), mean)

#mse.hc.center <- sum((cent.real[,-1]-cent.hc[,-1])^2)/(catalog[line,"nclass"]*(catalog[line,"cnumber"]-1))
#mse.hc.smote.center <- sum((cent.real[,-1]-cent.hc.smote[,-1])^2)/(catalog[line,"nclass"]*(catalog[line,"cnumber"]-1))


## metrics
distincemat <- dist(databody)


tryCatch({
performancelist <- rbind(performancelist,
data.frame(catalog[line,],
           
           #mse.hc.center = mse.hc.center,
           #mse.hc.smote.center=mse.hc.smote.center,
           
           accuracy.hc = CM.hc$overall["Accuracy"],
           accuracy.hc.smote = CM.smote.hc$overall["Accuracy"],
           
           Kappa.hc= CM.hc$overall["Kappa"],
           Kappa.hc.smote = CM.smote.hc$overall["Kappa"],
           
           sill.base = mean(silhouette(as.integer(datacluster), distincemat)[,3]),
           sill.hc = mean(silhouette(as.integer(hc.real.class), distincemat)[,3]),
           sill.hc.smote = mean(silhouette(as.integer(hc.smote.real.class), distincemat)[,3]),
           
           dunn.base =clValid::dunn(distincemat,as.integer(datacluster)),
           dunn.hc =clValid::dunn(distincemat,as.integer(hc.real.class)),
           dunn.hc.smote =clValid::dunn(distincemat,as.integer(hc.smote.real.class)),
           
           purity.hc = funtimes::purity(datacluster,hc.real.class)$pur,
           purity.hc.smote = funtimes::purity(datacluster,hc.smote.real.class)$pur,
           
           mut.inf.hc = aricode::NMI(hc.real.class,datacluster),
           mut.inf.hc.smote = aricode::NMI(hc.smote.real.class,datacluster),
           mut.inf.hc.smoteandhc = aricode::NMI(hc.smote.real.class,hc.real.class),
           
           
           var.inf.hc =mcclust::vi.dist(hc.real.class,datacluster),
           var.inf.hc.smote =mcclust::vi.dist(hc.smote.real.class,datacluster)
           
           ) )
}, error= function(e) print(paste("data has no cluster line number",line )))



if(catalog[line,"cnumber"]<4){
jpeg(paste0('plots\\',catalog[line,"line"],'.jpg'))
par(mfrow=c(3,1))

with(databody, plot(databody, col=datacluster,xlab="Base1",main=paste("Number :" ,catalog[line,"line"] ,"  Data :",as.character(catalog[line,"name"]))))
with(databody, plot(databody, col=hc.real.class,xlab="Hiearchical"))
with(databody, plot(databody, col=hc.smote.real.class,xlab="Smoted Hiearchical"))

 dev.off()
}

 saveRDS(performancelist,"performancelist.rds")

}



writexl::write_xlsx(performancelist,"performancelist.xlsx")

mean(performancelist$accuracy.hc-performancelist$accuracy.hc.smote)



## Performance Analysis

library(corrplot)

M <- cor(performancelist[,c(3:4,6:23)])
corrplot(M, method = "circle")



























