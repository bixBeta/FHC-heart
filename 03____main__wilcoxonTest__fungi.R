#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#eta_pavian = read.csv("global_Data/hear_meta_pavian.csv", header = T)
#rownames(meta_pavian) = meta_pavian$simpleName

norm.fungi.list = list()

for(i in 1:nrow(meta_pavian)){
  
  norm.fungi.list[[i]] <- raw.matrix.fungi[,i]/meta_pavian$Fungal.reads[i]  
  
}

names(norm.fungi.list) <- rownames(meta_pavian)
norm.fungi.counts <- do.call("rbind", norm.fungi.list)
colnames(norm.fungi.counts) <- rownames(raw.matrix.fungi)
input.matrix = t(norm.fungi.counts)
colnames(input.matrix) <- colnames(raw.matrix.fungi)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#colnames(input.matrix) <- rownames(samplemeta)

ctrl.matrix = as.data.frame(input.matrix) %>% select(control.meta$simpleName)

acute.matrix = as.data.frame(input.matrix) %>% select(acute.meta$simpleName)

chronic.matrix = as.data.frame(input.matrix) %>% select(chronic.meta$simpleName)

save(ctrl.matrix, acute.matrix, chronic.matrix, input.matrix, meta_pavian, raw.matrix.fungi, samplemeta, file = "data/02__HEART__fungi__all__matrices__norm__and__raw.Rdata")


### Remove AMC samples
acute.matrix = acute.matrix[,which(!colnames(acute.matrix) %in% unclass(unname(amc.meta %>% filter(aacc == "ACUTE") %>% select(simpleName)))[[1]])]
chronic.matrix =  chronic.matrix[,which(!colnames(chronic.matrix) %in% unclass(unname(amc.meta %>% filter(aacc == "CHRONIC") %>% select(simpleName)))[[1]])]


#------------------------------------------------------------
### Acute
#------------------------------------------------------------
res.acute = list()

for (i in 1:nrow(ctrl.matrix)) {
  res.acute[[i]] <- wilcox.test(x = unlist(unname(acute.matrix[i,])), y = unlist(unname(ctrl.matrix[i,])),
                                alternative = "two.sided", paired = F, conf.level = "0.95", exact = F, correct = F)
  
}

names(res.acute) <- rownames(ctrl.matrix)

p.vals.acute = vector()
for (i in 1:length(res.acute)) {
  p.vals.acute <- c(p.vals.acute, res.acute[[i]]$p.value)
}

#------------------------------------------------------------
### Chronic
#------------------------------------------------------------
res.chronic = list()

for (i in 1:nrow(ctrl.matrix)) {
  res.chronic[[i]] <- wilcox.test(x = unlist(unname(chronic.matrix[i,])), y = unlist(unname(ctrl.matrix[i,])),
                                  alternative = "two.sided", paired = F, conf.level = "0.95", exact = F, correct = F)
  
}

names(res.chronic) <- rownames(ctrl.matrix)

p.vals.chronic = vector()
for (i in 1:length(res.chronic)) {
  p.vals.chronic <- c(p.vals.chronic, res.chronic[[i]]$p.value)
}

#------------------------------------------------------------
### Acute v. Chronic
#------------------------------------------------------------
res.acute.v.chronic = list()

for (i in 1:nrow(ctrl.matrix)) {
  res.acute.v.chronic[[i]] <- wilcox.test(x = unlist(unname(acute.matrix[i,])), y = unlist(unname(chronic.matrix[i,])),
                                          alternative = "two.sided", paired = F, conf.level = "0.95", exact = F, correct = F)
  
}

names(res.acute.v.chronic) <- rownames(ctrl.matrix)

p.vals.acute.v.chronic = vector()
for (i in 1:length(res.acute.v.chronic)) {
  p.vals.acute.v.chronic <- c(p.vals.acute.v.chronic, res.acute.v.chronic[[i]]$p.value)
}


#------------------------------------------------------------
save(res.acute, res.chronic, res.acute.v.chronic, p.vals.acute, p.vals.chronic, p.vals.acute.v.chronic, file = "results/fungi__wilcoxon__results.Rdata")
#------------------------------------------------------------




