#------------------------------------------------------------
### chronic
#------------------------------------------------------------
# load("all__matrices__norm__and__raw.Rdata")


chronic.InputMatrix = cbind(ctrl.matrix, chronic.matrix)


library(dplyr)
library(purrr)
library(tidyr)

chronic.randomizeLabels <- function(inputMatrix){
  
  # get random order of sample labels
  chronic.permuted.labels <- sample(colnames(inputMatrix), size = ncol(inputMatrix), replace = F)
  
  # assign random labels to colnames of inputMatrix
  colnames(inputMatrix) <- chronic.permuted.labels
  
  # subset the get the new matrix of cases
  caseM <- as.data.frame(inputMatrix) %>% select(unname(unlist(filter(chronic.meta, Type != "AMCffpe") %>% select(simpleName))))
  
  # subset to get the new matrix of ctrls
  ctrlM <- as.data.frame(inputMatrix) %>% select(control.meta$simpleName)
  
  
  # RUN WILCOXON TEST
  # create an empty list
  chronic.permuted.res <- list()
  for (i in 1:nrow(inputMatrix)) {
    chronic.permuted.res[[i]] <- wilcox.test(x = unlist(unname(caseM[i,])), y = unlist(unname(ctrlM[i,])),
                                           alternative = "two.sided", paired = F, conf.level = "0.95", exact = F, correct = F)
    
  }
  
  names(chronic.permuted.res) <- rownames(inputMatrix)
  
  return(chronic.permuted.res)
  
}

MEGA.chronic.permuted.res = list()

for (i in 1:1000) {
  
  MEGA.chronic.permuted.res[[i]] <- chronic.randomizeLabels(inputMatrix = chronic.InputMatrix)
  
}

saveRDS(MEGA.chronic.permuted.res, "results/bacteria_MEGA.chronic.permuted.res.Rds")


getPs.chronic <- function(x){
  
  LCf = list()
  
  for (i in 1:length(MEGA.chronic.permuted.res)) {
    
    
    LCf[[i]] <-  pluck(MEGA.chronic.permuted.res, i, x, 'p.value')
    
  }
  
  return(unlist(LCf))
  
}


allPs.list.chronic = list()

for (i in 1:nrow(ctrl.matrix)) {
  allPs.list.chronic[[i]] <- getPs.chronic(x = i)
}

names(allPs.list.chronic) <- rownames(ctrl.matrix)
chronic.1000.table = do.call("rbind", allPs.list.chronic)

props.chronic.list = list()
for (i in 1:nrow(chronic.1000.table)) {
  
  props.chronic.list[[i]] = round(as.data.frame(chronic.1000.table) %>% slice(i),3)  %>% select_if(~any(. <= round(p.vals.chronic, 3)[i])) %>% length() /1000
  names(props.chronic.list)[[i]] = rownames(chronic.1000.table)[i]
}

chronic.FDR = as.data.frame(do.call("rbind", props.chronic.list)) 
chronic.FDR$p.value =  p.vals.chronic
colnames(chronic.FDR)[1] = "p.adj"
chronic.RESULTS = cbind(chronic.FDR, chronic.InputMatrix)

write.csv(chronic.RESULTS, "results/bacteria__chronic.v.ctrl.Results.csv")

