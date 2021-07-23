meta_pavian = read.csv("global_Data/heart51___REPORT_FILES-summary-210713.csv", header = T)
colnames(raw.matrix.bacteria) == colnames(raw.matrix.fungi)

colnames(raw.matrix.fungi)

# remove unwanted samples
raw.matrix.bacteria = raw.matrix.bacteria %>% select(- "K91_1_val_1.report-zero-counts.conf25.filt.kraken.covered",
                      - "PK92_1_val_1.report-zero-counts.conf25.filt.kraken.covered",
                      - "s_10993625_CKDL200147435-1a-AK5007-AK7215_HTFC2DSXX_L2_1_val_1.report-zero-counts.conf25.filt.kraken.covered",
                      - "s_11043641_CKDL200147435-1a-AK4243-AK7216_HTFC2DSXX_L2_1_val_1.report-zero-counts.conf25.filt.kraken.covered")
raw.matrix.fungi= raw.matrix.fungi %>% select(- "K91_1_val_1.report-zero-counts.conf25.filt.kraken.covered",
                                                     - "PK92_1_val_1.report-zero-counts.conf25.filt.kraken.covered",
                                                     - "s_10993625_CKDL200147435-1a-AK5007-AK7215_HTFC2DSXX_L2_1_val_1.report-zero-counts.conf25.filt.kraken.covered",
                                                     - "s_11043641_CKDL200147435-1a-AK4243-AK7216_HTFC2DSXX_L2_1_val_1.report-zero-counts.conf25.filt.kraken.covered")
raw.matrix.virus = raw.matrix.virus %>% select(- "K91_1_val_1.report-zero-counts.conf25.filt.kraken.covered",
                                                     - "PK92_1_val_1.report-zero-counts.conf25.filt.kraken.covered",
                                                     - "s_10993625_CKDL200147435-1a-AK5007-AK7215_HTFC2DSXX_L2_1_val_1.report-zero-counts.conf25.filt.kraken.covered",
                                                     - "s_11043641_CKDL200147435-1a-AK4243-AK7216_HTFC2DSXX_L2_1_val_1.report-zero-counts.conf25.filt.kraken.covered")


colnames(raw.matrix.bacteria) = samplemeta$simpleName
colnames(raw.matrix.virus) = samplemeta$simpleName
colnames(raw.matrix.fungi) = samplemeta$simpleName


chronic.meta = filter(samplemeta, aacc== "CHRONIC") %>% select(simpleName, aacc, Type, Group.1)
acute.meta = filter(samplemeta, aacc== "ACUTE") %>% select(simpleName, aacc, Type, Group.1)
amc.meta = filter(samplemeta, Type== "AMCffpe") %>% select(simpleName, aacc, Type, Group.1)
control.meta = filter(samplemeta, Group.1 == "control") %>% select(simpleName, aacc, Type, Group.1)


chronic.meta$Type[5] = "fresh"
chronic.meta$Type[6] = "fresh"
chronic.meta$Group.1[5]="case"
chronic.meta$Group.1[6]="case"

rownames(samplemeta) = samplemeta$simpleName




