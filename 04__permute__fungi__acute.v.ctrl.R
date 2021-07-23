#------------------------------------------------------------
### Acute
#------------------------------------------------------------
# load("all__matrices__norm__and__raw.Rdata")

# remove outlier sample
# acute.meta = acute.meta[-14,]
# acute.matrix = acute.matrix %>% select(-s_11043640)


acute.InputMatrix = cbind(ctrl.matrix, acute.matrix)


library(dplyr)
library(purrr)
library(tidyr)

acute.randomizeLabels <- function(inputMatrix){
  
  # get random order of sample labels
  acute.permuted.labels <- sample(colnames(inputMatrix), size = ncol(inputMatrix), replace = F)
  
  # assign random labels to colnames of inputMatrix
  colnames(inputMatrix) <- acute.permuted.labels
  
  # subset the get the new matrix of cases
  caseM <- as.data.frame(inputMatrix) %>% select(unname(unlist(filter(acute.meta, Type != "AMCffpe") %>% select(simpleName))))
  
  # subset to get the new matrix of ctrls
  ctrlM <- as.data.frame(inputMatrix) %>% select(control.meta$simpleName)
  
  
  # RUN WILCOXON TEST
  # create an empty list
  acute.permuted.res <- list()
  for (i in 1:nrow(inputMatrix)) {
    acute.permuted.res[[i]] <- wilcox.test(x = unlist(unname(caseM[i,])), y = unlist(unname(ctrlM[i,])),
                                           alternative = "two.sided", paired = F, conf.level = "0.95", exact = F, correct = F)
    
  }
  
  names(acute.permuted.res) <- rownames(inputMatrix)
  
  return(acute.permuted.res)
  
}

MEGA.acute.permuted.res = list()

for (i in 1:1000) {
  
  MEGA.acute.permuted.res[[i]] <- acute.randomizeLabels(inputMatrix = acute.InputMatrix)
  
}

saveRDS(MEGA.acute.permuted.res, "results/fungi_MEGA.acute.permuted.res.Rds")


getPs.acute <- function(x){
  
  LCf = list()
  
  for (i in 1:length(MEGA.acute.permuted.res)) {
    
    
    LCf[[i]] <-  pluck(MEGA.acute.permuted.res, i, x, 'p.value')
    
  }
  
  return(unlist(LCf))
  
}


allPs.list.Acute = list()

for (i in 1:nrow(ctrl.matrix)) {
  allPs.list.Acute[[i]] <- getPs.acute(x = i)
}

names(allPs.list.Acute) <- rownames(ctrl.matrix)
acute.1000.table = do.call("rbind", allPs.list.Acute)

props.acute.list = list()
for (i in 1:nrow(acute.1000.table)) {
  
  props.acute.list[[i]] = round(as.data.frame(acute.1000.table) %>% slice(i),3)  %>% select_if(~any(. <= round(p.vals.acute, 3)[i])) %>% length() /1000
  names(props.acute.list)[[i]] = rownames(acute.1000.table)[i]
}

acute.FDR = as.data.frame(do.call("rbind", props.acute.list)) 
acute.FDR$p.value =  p.vals.acute
colnames(acute.FDR)[1] = "p.adj"
acute.RESULTS = cbind(acute.FDR, acute.InputMatrix)

write.csv(acute.RESULTS, "results/fungi__acute.v.ctrl.Results.csv")

