files = list.files("REPORT.FILES/", pattern = "\\.report$", full.names = T)
# gsub(pattern = "\\.report$", "", files)

kraken.reports = list()
for (i in 1:length(files)) {
  x <- read.delim(files[i], header=FALSE, quote="")
  colnames(x) <- c("percent.frags", "covered", "assigned", "rank", "id","name")
  kraken.reports[[i]] <- x
  names(kraken.reports)[[i]] <- gsub(pattern = "\\.report$", "", basename(files[i]))
}
names(kraken.reports)


library(dplyr)
library(purrr)


D.index = lapply(kraken.reports, FUN = function(x){
  
  x[which(x$rank == "D"),]  
  
})

all.bacteria.list = lapply(kraken.reports, FUN = function(x){
  
  x[1373507:1818835,]  
  
})

all.fungi.list = lapply(kraken.reports, FUN = function(x){
  
  x[945818:1107875,]
  
})

all.virus.list = lapply(kraken.reports, FUN = function(x){
  
  x[1829922:2044738,]
  
})

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# bacteria setup
all.bacteria.sorted.list = lapply(all.bacteria.list, function(x){x %>% arrange(name)})
all.virus.sorted.list = lapply(all.virus.list, function(x){x %>% arrange(name)})
all.fungi.sorted.list = lapply(all.fungi.list, function(x){x %>% arrange(name)})

bacteria.raw.table = do.call("cbind", all.bacteria.sorted.list)
bacteria.table = bacteria.raw.table %>% select(matches(c(".name", ".covered", "rank")))


rownames(bacteria.table) = paste0(rownames(bacteria.table), ":", bacteria.table$`Cs_10373653_CKDL200147435-1a-AK7209-AK7215_HTFC2DSXX_L2_1_val_1.report-zero-counts.conf25.filt.kraken.name`)
bacteria.rank.matrix.full = bacteria.table %>% select(matches("rank"))
bacteria.table.covs = bacteria.table %>% select(matches("covered"))
raw.matrix.bacteria = bacteria.table.covs[rowSums(bacteria.table.covs > 10)>3,]
raw.matrix.bacteria = raw.matrix.bacteria[1:nrow(raw.matrix.bacteria) -1,]

saveRDS(bacteria.table , "data/bacteria.table.RDS")
saveRDS(raw.matrix.bacteria, "data/bacteria__raw__matrix__min10__atleaset3.RDS")


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# fungi setup
fungi.raw.table = do.call("cbind", all.fungi.sorted.list)
fungi.table = fungi.raw.table %>% select(matches(c(".name", ".covered", "rank")))


rownames(fungi.table) = paste0(rownames(fungi.table), ":", fungi.table$`Cs_10373653_CKDL200147435-1a-AK7209-AK7215_HTFC2DSXX_L2_1_val_1.report-zero-counts.conf25.filt.kraken.name`)
fungi.rank.matrix.full = fungi.table %>% select(matches("rank"))
fungi.table.covs = fungi.table %>% select(matches("covered"))
raw.matrix.fungi = fungi.table.covs[rowSums(fungi.table.covs > 10)>3,]
raw.matrix.fungi = raw.matrix.fungi[1:nrow(raw.matrix.fungi) -1,]

saveRDS(fungi.table , "data/fungi.table.RDS")
saveRDS(raw.matrix.fungi, "data/fungi__raw__matrix__min10__atleaset3.RDS")



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# virus setup
virus.raw.table = do.call("cbind", all.virus.sorted.list)
virus.table = virus.raw.table %>% select(matches(c(".name", ".covered", "rank")))

rownames(virus.table) = paste0(rownames(virus.table), ":", virus.table$`Cs_10373653_CKDL200147435-1a-AK7209-AK7215_HTFC2DSXX_L2_1_val_1.report-zero-counts.conf25.filt.kraken.name`)
virus.rank.matrix.full = virus.table %>% select(matches("rank"))
virus.table.covs = virus.table %>% select(matches("covered"))
raw.matrix.virus = virus.table.covs[rowSums(virus.table.covs > 10)>3,]
# raw.matrix.virus = raw.matrix.virus[1:nrow(raw.matrix.virus) -1,]

saveRDS(virus.table , "data/virus.table.RDS")
saveRDS(raw.matrix.virus, "data/virus__raw__matrix__min10__atleaset3.RDS")




