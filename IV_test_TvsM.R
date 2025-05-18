library(ivreg)
library(AER)
library(dplyr)
library(openxlsx)
data <- read.csv(file ="/Users/QIUAODON/Desktop/wholeT_communication/IV_analysis/gene_df_PD1_Myeloid.csv", header = TRUE, sep = ',', stringsAsFactors = FALSE)

# Initialize an empty dataframe to store the results
results_df <- data.frame(
  gene_T = character(),
  gene_cell = character(),
  p_value = numeric(),
  r_squared = numeric(),
  stringsAsFactors = FALSE
)
cell_types = c('B', 'M', 'Epi', 'Endo', 'Fibro')

# sort the DEG pairs according to P_value < 0.01 and R squared > 0.5

for (cell_type in cell_types) {
  DEG_exp_T = select(data, ends_with('_T'))
  DEG_exp_cell = select(data, ends_with(paste0('_', cell_type)))
  
  #get through every gene pair
  for (gene_T in names(DEG_exp_T)) {
    for (gene_cell in names(DEG_exp_cell)) {
      formula <- as.formula(paste0(gene_cell, " ~ ", gene_T, " | treatment"))
      
      #IV regression
      result <- ivreg(formula, data = data)
      p_value = summary(result)$coefficients[2, 4]
      adj_r_squared <- summary(result)$adj.r.squared
      if (p_value < 0.01) {
        print(paste("gene pair:", gene_T, "vs", gene_cell))
        # IV_result_low_p[[paste("gene pair:", gene_T, "vs", gene_cell)]] <- list(result = result, p_value = p_value)
        # Append the result to the dataframe
        results_df <- rbind(results_df, data.frame(
          gene_T = gene_T,
          gene_cell = gene_cell,
          p_value = p_value,
          r_squared = adj_r_squared,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}

# Calculate q-values for the collected p-values
results_df$q_value <- p.adjust(results_df$p_value, method = "fdr")

write.xlsx(results_df, "/Users/QIUAODON/Desktop/wholeT_communication//IV_analysis/IV_regression_results_PD1vsM_GEM_including_qval.xlsx")
