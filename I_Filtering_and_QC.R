library(tidyverse)
library(readr)
library(corrplot)
library(psych)
library(edgeR)
library(RColorBrewer)

#set working directory to Source File: Session --> Set Working Directory --> To Source File Location

source("./Functions.R")

file = "path_to_file/file.txt" #path to proteinGroups.txt

proteinGroups <- read_proteinGroups(file)
head(proteinGroups)


## LFQ table
# Modify table to have:
# - only necessary columns
# - remove common contaminants (Reverse, Only identified by site, Potential contaminants)
# - keep only proteins with unique.peptide >= 1

# check number of contaminats to remove them
sprintf("The number of contaminants are %d + %d + %d", 
      table(is.na(proteinGroups$`Only identified by site`))[1],
      table(is.na(proteinGroups$Reverse))[1],
      table(is.na(proteinGroups$`Potential contaminant`))[1])

sprintf("The number of proteins with at least on unique peptide is %d", 
        table(proteinGroups$`Unique peptides` >= 1)[1])


clean_proteinGroups <- filtering(proteinGroups)
anno_tab <- annotation_table(proteinGroups)



# Convert LFQ raw values in Log2
LFQ_tab <- clean_proteinGroups |>
  mutate(across(starts_with("LFQ"), log2)) |>
  mutate(across(starts_with ("LFQ"), ~ na_if(., -Inf)))

## QC plot
barplot(apply(LFQ_tab[-1],2,function(x) !is.na(x)), las = 2) #number of proteins per sample

LFQ_tab |>
  column_to_rownames("Protein_Gene") |>
  boxplot(outline = F, las = 2) #boxplot to see distribution is centered around the median for all samples.

LFQ_tab |>
  column_to_rownames("Protein_Gene") |>
  cor(use = "na.or.complete") |>
  corrplot(method = "square", addCoef.col = T) #correlation among replicates and across conditions

pairs.panels(LFQ_tab[,-1], 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
) #summary ploy with correlations and histograms of intensities



colCondition <- c(rep("lightblue",2), rep("violet",4),
                  rep("orange",4), rep("green",4)) #select colors for pca according to the number of replicates
plotMDS(LFQ_tab[,-1], pch = 16, col = colCondition)
#legend("topleft", lwd = 2, cex =1, col = c("lightblue", "violet"),
#       legend = c("Cond1", "Cond2")) #name of the conditions



#remove unnecessary objects prior to subsequent analysis (save space in R)
rm(clean_proteinGroups, proteinGroups, colCondition, file)


write.csv(LFQ_tab,"./Filtered_proteinGroups.csv", row.names = F)
write.csv(anno_tab,"./Annotation_proteinGroups.csv", row.names = F)
