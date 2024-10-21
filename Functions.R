read_proteinGroups <- function(file = file){
  col_names <- colnames(read_delim(file, delim = "\t", n_max = 0))
  lfq_cols <- grep("LFQ", col_names, value = TRUE)  # Get columns that start with 'LFQ'
  col_spec <- cols(.default = col_guess(),  # Default guess for other columns
                   Reverse = col_character(),
                   `Only identified by site` = col_character(),
                   `Potential contaminants` = col_character())  # Explicitly set Contiaminants column type
  
  # Add LFQ columns as col_number() dynamically
  for (col in lfq_cols) {
    col_spec[[col]] <- col_number()
  }
  
  data <- read_delim(file, 
                     delim = "\t", 
                     escape_double = FALSE, 
                     col_types = col_spec, 
                     comment = "#", 
                     trim_ws = TRUE)
  
  data <- data |>
    separate(col = `Majority protein IDs`, 
             into = "ProteinID", 
             sep = ";", remove = TRUE) |>
    separate(col = `Gene names`, 
             into = "GeneID",
             sep = ";", remove = TRUE) |>
    unite(col = "Protein_Gene", 
          ProteinID, GeneID, 
          sep = "_", remove = TRUE)
  
  return(data)
}

filtering <- function(data_table = proteinGroups){
  data_filtered <- data_table |> 
    filter(is.na(Reverse), 
           is.na(`Only identified by site`), 
           is.na(`Potential contaminant`),
           `Unique peptides` >= 1) |>
    dplyr::select(Protein_Gene, contains("LFQ ")) |>
    rename_with(~gsub("intensity ","",.)) |>
    rename_with(~gsub(" ","_",.))
  
  return(data_filtered)
}

annotation_table <- function(data_table = proteinGroups){
  data_filtered <- data_table |> 
    filter(is.na(Reverse), 
           is.na(`Only identified by site`), 
           is.na(`Potential contaminant`),
           `Unique peptides` >= 1) |>
    dplyr::select(Protein_Gene, Peptides, contains("coverage"), contains("iBAQ "))
  return(data_filtered)
}

imputation <- function(lfq_table = LFQ_tab, type_of_imputation = "mean"){
  #parameters of the distribution
  mean_table <- lfq_table |>
    dplyr::select(Protein_Gene, contains("LFQ")) |>
    as.matrix() |>
    mean(na.rm = T) 

  if (type_of_imputation == "zero") {
    imputed_tab <- lfq_table |>
      dplyr::select(Protein_Gene, contains("LFQ")) |>
      column_to_rownames("Protein_Gene") |>
      mutate(across(everything(), ~ replace_na(. , 0)))
  }
  
  else if (type_of_imputation == "min") {
    imputed_tab <- lfq_table |>
      dplyr::select(Protein_Gene, contains("LFQ")) |>
      column_to_rownames("Protein_Gene") |>
      mutate(across(everything(), ~ replace_na(. , min(lfq_table |> dplyr::select(contains("LFQ")), na.rm = T))))
  }
  
  else if (type_of_imputation == "mean") {
    imputed_tab <- lfq_table |>
      dplyr::select(Protein_Gene, contains("LFQ")) |>
      column_to_rownames("Protein_Gene") |>
      mutate(across(everything(), ~ replace_na(. , mean_table)))
  }
  
  else {
    warning("You did not insert a valid imputation option. Valid options are 'zero', 'minimum', or 'mixed'.")
    return(NULL)
  }
    
  return(imputed_tab)
}

diff_expr_analysis <- function(table){
  mat <- table |> as.matrix()
  tipo <- "Protein"
  cond <- colnames(mat) |> str_extract("(?<=_).+(?=_)")
  group <- interaction(tipo, cond)
  mm <- model.matrix(~0 + group)
  fit <- lmFit(mat, mm)
   
  return(fit)
}


vulcano_proteomics <- function(stat_table = stat_table, is_stringent = TRUE, p_adj_cutoff = 0.05, p_value_cutoff = 0.05, lfc_cutoff = 1){
  if(is_stringent){
    stat_table |>
      mutate(P = if_else(adj.P.Val < p_adj_cutoff, T, F),
             UP = if_else(logFC > lfc_cutoff, T, F),
             DOWN = if_else(logFC < -lfc_cutoff, T, F),
             Sign = case_when(
               P & UP ~ "UP",
               P & DOWN ~ "DOWN",
               TRUE ~ "nr")) |>
      # Gene = vapply(strsplit(Protein,"_"), `[`, 2, FUN.VALUE=character(1)),
      # Label = if_else(Gene %in% complex, Gene, NA)) |>
      ggplot(aes(x=logFC, y=-log10(P.Value), col=Sign)) +
      geom_point(alpha=0.5)+
      scale_color_manual(values = c("orange","grey", "purple"))+
      #geom_text(aes(label = Label), vjust = 1.5, parse = T)+
      theme_bw()
  }
  else{
    stat_table |>
      mutate(P = if_else(P.Value < p_value_cutoff, T, F),
             UP = if_else(logFC > lfc_cutoff, T, F),
             DOWN = if_else(logFC < -lfc_cutoff, T, F),
             Sign = case_when(
               P & UP ~ "UP",
               P & DOWN ~ "DOWN",
               TRUE ~ "nr")) |>
      # Gene = vapply(strsplit(Protein,"_"), `[`, 2, FUN.VALUE=character(1)),
      # Label = if_else(Gene %in% complex, Gene, NA)) |>
      ggplot(aes(x=logFC, y=-log10(P.Value), col=Sign)) +
      geom_point(alpha=0.5)+
      scale_color_manual(values = c("orange","grey", "purple"))+
      #geom_text(aes(label = Label), vjust = 1.5, parse = T)+
      theme_bw()
  }
}

vulcano_IP <- function(stat_table = stat_table, is_stringent = TRUE, p_adj_cutoff = 0.05, p_value_cutoff = 0.05, lfc_cutoff = 1){
  if(is_stringent){
    stat_table |>
      mutate(Sign = if_else(adj.P.Val < p_adj_cutoff, T, F),
             UP = if_else(logFC > lfc_cutoff, T, F),
             Sign.UP = if_else(Sign & UP, "+", "-")) |>
             # Gene = vapply(strsplit(Protein,"_"), `[`, 2, FUN.VALUE=character(1)),
             # Label = if_else(Gene %in% complex, Gene, NA)) |>
      ggplot(aes(x=logFC, y=-log10(P.Value), col=Sign.UP)) +
      geom_point(alpha=0.5)+
      scale_color_manual(values = c("grey", "blue"))+
      #geom_text(aes(label = Label), vjust = 1.5, parse = T)+
      theme_bw()
  }
  else{
    stat_table |>
      mutate(Sign = if_else(P.Value < p_value_cutoff, T, F),
             UP = if_else(logFC > lfc_cutoff, T, F),
             Sign.UP = if_else(Sign & UP, "+", "-")) |>
      # Gene = vapply(strsplit(Protein,"_"), `[`, 2, FUN.VALUE=character(1)),
      # Label = if_else(Gene %in% complex, Gene, NA)) |>
      ggplot(aes(x=logFC, y=-log10(P.Value), col=Sign.UP)) +
      geom_point(alpha=0.5)+
      scale_color_manual(values = c("grey", "blue"))+
      #geom_text(aes(label = Label), vjust = 1.5, parse = T)+
      theme_bw()
  }
}

generate_vulcano <- function(stat_table = stat_table, is_stringent = TRUE, p_adj_cutoff = 0.05, p_value_cutoff = 0.05, lfc_cutoff = 1, is_IP = TRUE){
  if(is_IP){
    print("Generating vulcano for IP...")
    vulcano_IP(stat_table, is_stringent, p_adj_cutoff, p_value_cutoff, lfc_cutoff)
  }
  else{
    print("Generating vulcano for Proteomics...")
    vulcano_proteomics(stat_table, is_stringent, p_adj_cutoff, p_value_cutoff, lfc_cutoff)
  }
}


mixed_imputation <- function(data, condition_sets, control_sets, missing_thresholds){
  #parameters of the distribution
  mean_table <- filtered_table |>
    column_to_rownames("Protein_Gene") |>
    as.matrix() |>
    mean(na.rm = T)
  
  sd_value <- filtered_table |>
    column_to_rownames("Protein_Gene") |>
    as.matrix() |>
    sd(na.rm = T)
  
  #count number of missign values per condition
  for (i in seq_along(condition_sets)) {
    cond_cols <- condition_sets[[i]]
    ctrl_cols <- control_sets[[1]]
    cond_missing_threshold <- missing_thresholds$cond[[i]]
    ctrl_missing_threshold <- missing_thresholds$ctrl[[1]]
    
    # Compute mean and count NAs for each condition and control set
    data <- data |>
      rowwise() |>
      mutate(
        !!paste0("mean_Cond_", i) := mean(c_across(all_of(cond_cols)), na.rm = TRUE),
        !!paste0("mean_CTRL_", 1) := mean(c_across(all_of(ctrl_cols)), na.rm = TRUE),
        !!paste0("NA_count_Cond_", i) := sum(is.na(c_across(all_of(cond_cols)))),
        !!paste0("NA_count_CTRL_", 1) := sum(is.na(c_across(all_of(ctrl_cols))))
      )
  }
  
  data <- data |> ungroup()
  
  #Imputation logic for each condition/control group
  for (i in seq_along(condition_sets)) {
    cond_cols <- condition_sets[[i]]
    ctrl_cols <- control_sets[[1]]
    cond_missing_threshold <- missing_thresholds$cond[[i]]
    ctrl_missing_threshold <- missing_thresholds$ctrl[[1]]
    
    # Apply imputation based on the calculated means and thresholds
    data <- data |>
      mutate(
        across(all_of(ctrl_cols), 
               ~ ifelse(is.na(.) & !!sym(paste0("NA_count_CTRL_", 1)) <= ctrl_missing_threshold, 
                        !!sym(paste0("mean_CTRL_", 1)), .)),
        across(all_of(ctrl_cols), 
               ~ ifelse(is.na(.), 
                        rnorm(sum(is.na(.)), mean = (mean_table - 1.8 * sd_value), sd = (0.3 * sd_value)), .)),
        across(all_of(cond_cols), 
               ~ ifelse(is.na(.) & !!sym(paste0("NA_count_Cond_", i)) <= cond_missing_threshold, 
                        !!sym(paste0("mean_Cond_", i)), .)),
        across(all_of(cond_cols), 
               ~ ifelse(is.na(.) & !!sym(paste0("NA_count_Cond_", i)) > cond_missing_threshold, 
                        rnorm(sum(is.na(.)), mean = (mean_table - 1.8 * sd_value), sd = (0.3 * sd_value)), .))
      )
  }
  data <- data |> column_to_rownames("Protein_Gene") |> dplyr::select(contains("LFQ"))
  
  
  return(data)
}
