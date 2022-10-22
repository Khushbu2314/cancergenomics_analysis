library(survival)
library(survminer)
setwd("C:\\Users\\sai\\Documents\\R")
data = read.table('CGGA.mRNAseq_693_clinical.20200506.txt',
                  header = T, sep = '\t', row.names = 1)

fit =  survfit(Surv(OS, Censor..alive.0..dead.1.)  ~  Gender, 
               data = data)
fit


ggsurvplot(fit, data = data)

ggsurvplot(fit, data = data,
           surv.median.line = 'hv')

ggsurvplot(fit, data = data,
           surv.median.line = 'hv',
           pval = T)

ggsurvplot(fit, data = data,
           surv.median.line = 'hv',
           pval = T, risk.table = T)
fit =  survfit(Surv(OS, Censor..alive.0..dead.1.) ~ IDH_mutation_status, 
               data = data)

fit =  survfit(Surv(OS, Censor..alive.0..dead.1.) ~ IDH_mutation_status + Gender, 
               data = data)
