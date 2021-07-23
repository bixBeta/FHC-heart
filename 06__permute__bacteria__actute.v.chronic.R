#------------------------------------------------------------
### Acute v. Chronic
#------------------------------------------------------------
# load("all__matrices__norm__and__raw.Rdata")


acute.v.chronic.InputMatrix = cbind(acute.matrix, chronic.matrix)


library(dplyr)
library(purrr)
library(tidyr)

acute.v.chronic.randomizeLabels <- function(inputMatrix){
  # get random order of sample labels
  acute.v.chronic.permuted.labels <- sample(colnames(inputMatrix), size = ncol(inputMatrix), replace = F)
  
  # assign random labels to colnames of inputMatrix
  colnames(inputMatrix) <- acute.v.chronic.permuted.labels
  
  # subset the get the new matrix of acute
  caseM <- as.data.frame(inputMatrix) %>% select(unname(unlist(filter(acute.meta, Type != "AMCffpe") %>% select(simpleName))))
  
  # subset to get the new matrix of chronic
  ctrlM <- as.data.frame(inputMatrix) %>% select(unname(unlist(filter(chronic.meta, Type != "AMCffpe") %>% select(simpleName))))
  
  
  # RUN WILCOXON TEST
  # create an empty list
  acute.v.chronic.permuted.res <- list()
  for (i in 1:nrow(inputMatrix)) {
    acute.v.chronic.permuted.res[[i]] <- wilcox.test(x = unlist(unname(caseM[i,])), y = unlist(unname(ctrlM[i,])),
                                                     alternative = "two.sided", paired = F, conf.level = "0.95", exact = F, correct = F)
    
  }
  
  names(acute.v.chronic.permuted.res) <- rownames(inputMatrix)
  
  return(acute.v.chronic.permuted.res)
  
}

MEGA.acute.v.chronic.permuted.res = list()

for (i in 1:1000) {
  
  MEGA.acute.v.chronic.permuted.res[[i]] <- acute.v.chronic.randomizeLabels(inputMatrix = acute.v.chronic.InputMatrix)
  
}

saveRDS(MEGA.acute.v.chronic.permuted.res, "results/bacteria_MEGA.acute.v.chronic.permuted.res.Rds")


getPs.acute <- function(x){
  
  LCf = list()
  
  for (i in 1:length(MEGA.acute.v.chronic.permuted.res)) {
    
    
    LCf[[i]] <-  pluck(MEGA.acute.v.chronic.permuted.res, i, x, 'p.value')
    
  }
  
  return(unlist(LCf))
  
}


allPs.list.Acute = list()

for (i in 1:nrow(ctrl.matrix)) {
  allPs.list.Acute[[i]] <- getPs.acute(x = i)
}

names(allPs.list.Acute) <- rownames(ctrl.matrix)
acute.v.chronic.1000.table = do.call("rbind", allPs.list.Acute)

props.acute.v.chronic.list = list()
for (i in 1:nrow(acute.v.chronic.1000.table)) {
  
  props.acute.v.chronic.list[[i]] = round(as.data.frame(acute.v.chronic.1000.table) %>% slice(i),3)  %>% select_if(~any(. <= round(p.vals.acute.v.chronic, 3)[i])) %>% length() /1000
  names(props.acute.v.chronic.list)[[i]] = rownames(acute.v.chronic.1000.table)[i]
}

acute.v.chronic.FDR = as.data.frame(do.call("rbind", props.acute.v.chronic.list)) 
acute.v.chronic.FDR$p.value =  p.vals.acute
colnames(acute.v.chronic.FDR)[1] = "p.adj"
acute.v.chronic.RESULTS = cbind(acute.v.chronic.FDR, acute.v.chronic.InputMatrix)

write.csv(acute.v.chronic.RESULTS, "results/bacteria__acute.v.chronic.Results.csv")

