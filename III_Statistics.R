#statistics
library(tidyverse)
library(edgeR)

#set working directory to Source File: Session --> Set Working Directory --> To Source File Location

source("./Functions.R")

mixed_imp <- read.delim("./Mixed_Imputation_Table.csv", sep = ",")


fit <- diff_expr_analysis(mixed_imp)
contr <- makeContrasts(CondvsCtrl =  groupProtein.P1_dRW - groupProtein.E14,
                       levels = colnames(coef(fit))) #specify conditions to test
tmp <- contrasts.fit(fit, contr[,"CondvsCtrl"]) |>
  eBayes() 

stat_table <- topTable(tmp, sort.by = "P", n = Inf) |>
  rownames_to_column("Protein")

generate_vulcano(stat_table, 
                 is_IP = T, #put T if your dataset comes from an AP-MS, put F if your dataset is a standard proteomics
                 is_stringent = T, #put T if you want a stringent significqnce (p.adj), put F if you want relax (p.value)
                 p_adj_cutoff = 0.05, 
                 p_value_cutoff = 0.05, 
                 lfc_cutoff = 1)


#write statistical table
stat_table |>
  dplyr::select(-AveExpr, -t, -B) |>
  write.csv ("./Table_Statistics.csv", row.names = F)
 