## MISSING VALUES IMPUTATION
library(tidyverse)

#set working directory to Source File: Session --> Set Working Directory --> To Source File Location

source("./Functions.R")

LFQ_tab <- read.delim("./Filtered_proteinGroups.csv", sep = ",")

#remove dangerous columns according to qc
LFQ_tab <- LFQ_tab |>
  dplyr::select(-LFQ_P1_dRING_04)

# first count number of valid values per condition
LFQ_tab <-LFQ_tab |>
  rowwise() |> 
  mutate(count_tot = sum(!is.na(c_across(LFQ_E14_01:LFQ_P1_FL_04))),
         count_Cond1 = sum(!is.na(c_across(LFQ_E14_01:LFQ_E14_02))),
         count_Cond2 = sum(!is.na(c_across(LFQ_P1_dRING_01:LFQ_P1_dRING_03))),
         count_Cond3 = sum(!is.na(c_across(LFQ_P1_dRW_01:LFQ_P1_dRW_04))),
         count_Cond5 = sum(!is.na(c_across(LFQ_P1_FL_01:LFQ_P1_FL_04)))) |>
  ungroup()



# The first thing to do is to filter out rows withou enough number of valid replicates
# uncomment the filter that you want to use

#1 TUTTI I VALORI VALIDI IN TUTTE LE CONDIZIONI 
#table_only_valids <- 
#  table_temp |>  
#  filter(NA_count_tot == 0 ) |>
#  select(-NA_count_tot, -NA_count_Cond1, -NA_count_Cond2)

#2 50% DEI VALORI SONO NON NULLI (SUL TOTALE)
#table_50_noNA_tot <- 
#  table_temp |>  
#  filter(NA_count_tot <= 4 ) |>
#  select(-NA_count_tot, -NA_count_Cond1, -NA_count_Cond2)

#3 75% DEI VALORI SONO NON NULLI (SUL TOTALE)
#table_75_noNA_tot <- 
#  table_temp |>  
#  filter(NA_count_tot <= 2 ) |>
#  select(-NA_count_tot, -NA_count_Cond1, -NA_count_Cond2)

#4 50% DEI VALORI SONO NON NULLI (ALMENO IN UNA CONDIZIONE)
# table_50_noNA_atleast1 <- 
#   table_temp |>  
#   filter(NA_count_Cond1 <= 2 | NA_count_Cond2 <=2 ) |> #condition OR
#   select(-NA_count_tot, -NA_count_Cond1, -NA_count_Cond2)

#5 75% DEI VALORI SONO NON NULLI (ALMENO IN UNA CONDIZIONE)
filtered_table <- LFQ_tab[(LFQ_tab['count_Cond1'] >=1) |
                           (LFQ_tab['count_Cond2'] >= 2) |
                           (LFQ_tab['count_Cond3'] >= 2) |
                           (LFQ_tab['count_Cond5'] >=2),]

filtered_table <- filtered_table |> dplyr::select(-contains("count"))


## IMPUTATION
# missing values can be because they do not exists, at random, or they are too low to be detected.
# For the first and the last condition they can be simply imputed with a low number, which requires the generation of a new distribution
# For the second condition, we assume the missing value is close to the intensity of those identified, indeed their average.

## case:
# 1. sostituisco con 0 or with a fixed value
# 2. sostituisco con il valore minimo degli LFQ values della table
# 3. sostituisco con uno tra i valori che ottengo generando una distribuzione normale partendo dai valori 
#    LFQ; uso media e sd per creare gaussiana shiftata verso la parte bassa della distribuzione degli LFQ. 
#    Uso dei valori standard dati da Perseus. 30% della sd e -1.8 della media
# 4. sostituisco i missing del gruppo CTRL con la media di quella riga del CTRL e i missing del gruppo IP con
#    la media di quella riga dell'IP
# 5. MIXED IMPUTATION based on the number of missing values


#1
imputation_NA_0 <- imputation(filtered_table, "zero")
# write.csv (imputation_NA_0, "./Imputation_NA_0.csv")


#2
imputation_NA_min <- imputation(filtered_table, "min")
# write.csv (imputation_NA_min,"./Imputation_NA_min_value_table.csv")


#3
imputation_NA_mean <- imputation(filtered_table, "mean")


#3 : MIXED IMPUTATION
# Define multiple condition sets and control sets
condition_sets <- list(
  c("LFQ_P1_FL_01", "LFQ_P1_FL_02", "LFQ_P1_FL_03", "LFQ_P1_FL_04"),
  c("LFQ_P1_dRING_01", "LFQ_P1_dRING_02", "LFQ_P1_dRING_03"),
  c("LFQ_P1_dRW_01", "LFQ_P1_dRW_02", "LFQ_P1_dRW_03", "LFQ_P1_dRW_04")
  #c("LFQ_P2_dT_01", "LFQ_P2_dT_02", "LFQ_P2_dT_03", "LFQ_P2_dT_04")
)

control_sets <- list(
  c("LFQ_E14_01", "LFQ_E14_02")
)

# Define missing value thresholds for each condition and control group
missing_thresholds <- list(
  cond = c(2, 2, 1),  # Threshold for condition 1, 2
  ctrl = c(2)   # Threshold for control 1, 2
)

mixed_imp <- mixed_imputation(filtered_table, condition_sets, control_sets, missing_thresholds)



write.csv (mixed_imp,"./Mixed_Imputation_Table.csv", row.names = T)
  
