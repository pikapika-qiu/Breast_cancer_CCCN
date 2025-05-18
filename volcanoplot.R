library(ggplot2)
library(ggrepel)
library(dplyr)     

# Read the data from the Excel file
data <- readxl::read_excel("/Users/QIUAODON/Desktop/Pancancer_result_data/PDCD1exp_DEG/T_PDCD1+_DEG(pre_vs_on).xlsx")

data$logP = -log10(data$pval) 
data$log2fc = -data$log2fc
data$Gene = data$...1
# Define significance thresholds
threshold_pvalue <- 0.05
threshold_logfc <- 0.5

significant_genes <- data %>% 
  filter(abs(log2fc) > threshold_logfc & pval < threshold_pvalue)
# Filter significant genes
cluster_deg_df <- data
significant_genes <- cluster_deg_df %>%
  filter(abs(log2fc) > threshold_logfc & pval < threshold_pvalue)

# List of specific genes to label
genes_to_label <- c(
  'PRDM1', 'TXNIP', 'TSC22D3', 'FKBP5', 'NFKBIA', 'ZFP36L1',  'RGS1', 'CEMIP2', 
  'IRF1', 'CXCR4', 'FASLG', 'TNFSF14', 'HMGB2', 'TNFAIP3', 'SOCS1', 'IER2', 'TAGAP', 'CTLA4',  
  'PELI1', 'KLF10', 'IL6ST', 'IL7R',  'SOCS3', 'RGS16', 'NR4A2',
  'PMAIP1', 'TNFSF8', 'IL10', 'GADD45B', 'HMOX1', 'ATF3' , 'ICAM1' , 'CEBPD', 'TRAF1',
  'IRF4', 'CSF2',  'HPGD', 'AREG', 'CCNA2', 'MMP9', 'CXCL2', 'GADD45G',
  'IER3', 'KLF2', 'CISH', 'IL2', 'CSF1', 'LTA', 'MYB', 'CXCL11',   'CCL8' 
)

# Create volcano plot
volcano_plot <- ggplot(cluster_deg_df, aes(x = log2fc, y = logP)) +
  geom_point(aes(color = log2fc > threshold_logfc & pval < threshold_pvalue), 
             color = ifelse(cluster_deg_df$log2fc > threshold_logfc & cluster_deg_df$pval < threshold_pvalue, "red", 
                            ifelse(cluster_deg_df$log2fc < -threshold_logfc & cluster_deg_df$pval < threshold_pvalue, "blue", "gray"))) +
  scale_color_manual(values = c("gray", "red", "blue")) +
  labs(title = expression("Volcano Plot of " * italic("PDCD1") * "+ T cell DEGs"),
       x = "Log2 Fold Change",
       y = "-log10 p-value") +
  theme_minimal() +
  theme(legend.position = "none") +
  xlim(-2, 2) +
  geom_text_repel(data = cluster_deg_df %>% filter(Gene %in% genes_to_label),
                  aes(label = Gene), size = 3, fontface = "italic", 
                  box.padding = 0.5, point.padding = 0.6,
                  max.overlaps = 20) +  # Control the number of annotations

  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))  # Add box lines

# Print the plot
print(volcano_plot)

#--------------


# Define dynamic criteria for labeling genes (e.g., top 10 genes by p-value)
genes_to_label <- cluster_deg_df %>%
  filter(log2fc > threshold_logfc & pval < threshold_pvalue) %>%
  top_n(-20, pval) %>%  # Select the top 10 genes with the lowest p-value
  pull(Gene)

# Create volcano plot
volcano_plot <- ggplot(cluster_deg_df, aes(x = log2fc, y = logP)) +
  geom_point(aes(color = log2fc > threshold_logfc & pval < threshold_pvalue), 
             color = ifelse(cluster_deg_df$log2fc > threshold_logfc & cluster_deg_df$pval < threshold_pvalue, "red", 
                            ifelse(cluster_deg_df$log2fc < -threshold_logfc & cluster_deg_df$pval < threshold_pvalue, "blue", "gray"))) +
  scale_color_manual(values = c("gray", "red", "blue")) +
  labs(title = "Volcano Plot of " * italic("PDCD1") * "+ T cell DEGs",
       x = "Log2 Fold Change",
       y = "-log10 p-value") +
  theme_minimal() +
  theme(legend.position = "none") +
  xlim(-2, 2) +
  geom_text_repel(data = cluster_deg_df %>% filter(Gene %in% genes_to_label),
                  aes(label = Gene), size = 5,
                  box.padding = 0.35, point.padding = 0.5,
                  max.overlaps = 20) +  # Control the number of annotations
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))  # Add box lines

# Print the plot
print(volcano_plot)


#--------CD4EX

data <- read.csv("/Users/QIUAODON/Desktop/cluster3_CD4_subtypes_communication/CD4EX_DEGs.csv")

data$logP = -log10(data$pval) 
data$log2fc = -data$log2fc
data$Gene = data$X
# Define significance thresholds
threshold_pvalue <- 0.05
threshold_logfc <- 0.5

significant_genes <- data %>% 
  filter(abs(log2fc) > threshold_logfc & pval < threshold_pvalue)
# Filter significant genes
cluster_deg_df <- data
significant_genes <- cluster_deg_df %>%
  filter(abs(log2fc) > threshold_logfc & pval < threshold_pvalue)

# List of specific genes to label
genes_to_label <- c(
  'TXNIP', 'NFKBIA', 'TSC22D3', 'ZFP36L1', 'PRDM1', 'BIRC3', 'RGS1', 'PELI1', 'ZFP36', 'BHLHE40', 'IRF1', 'FOS', 'KLF2', 'GADD45G', 'GADD45B', 'IRF4', 'TNFAIP3', 'CBLB', 'JUN', 'TNF', 'CCL2', 'HMOX1', 'ATF3', 'SOCS3', 'NR4A1'
)

# Create volcano plot
volcano_plot <- ggplot(cluster_deg_df, aes(x = log2fc, y = logP)) +
  geom_point(aes(color = log2fc > threshold_logfc & pval < threshold_pvalue), 
             color = ifelse(cluster_deg_df$log2fc > threshold_logfc & cluster_deg_df$pval < threshold_pvalue, "red", 
                            ifelse(cluster_deg_df$log2fc < -threshold_logfc & cluster_deg_df$pval < threshold_pvalue, "blue", "gray"))) +
  scale_color_manual(values = c("gray", "red", "blue")) +
  labs(title = "Volcano Plot of CD4EX T cell DEGs",
       x = "Log2 Fold Change",
       y = "-log10 p-value") +
  theme_minimal() +
  theme(legend.position = "none") +
  xlim(-2, 2) +
  geom_text_repel(data = cluster_deg_df %>% filter(Gene %in% genes_to_label),
                  aes(label = Gene), size = 3, fontface = "italic", 
                  box.padding = 0.5, point.padding = 0.6,
                  max.overlaps = 20) +  # Control the number of annotations
  
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))  # Add box lines

# Print the plot
print(volcano_plot)

#--------CD8EX

# Read the data from the Excel file
data <- readxl::read_excel("/Users/QIUAODON/Desktop/cluster4_CD8_subtypes_communication/CD8EX_DEGs.xlsx")

data$logP = -log10(data$pval) 
data$log2fc = -data$log2fc
data$Gene = data$...1
# Define significance thresholds
threshold_pvalue <- 0.05
threshold_logfc <- 0.5

significant_genes <- data %>% 
  filter(abs(log2fc) > threshold_logfc & pval < threshold_pvalue)
# Filter significant genes
cluster_deg_df <- data
significant_genes <- cluster_deg_df %>%
  filter(abs(log2fc) > threshold_logfc & pval < threshold_pvalue)

# List of specific genes to label
genes_to_label <- c(
  'TSC22D3', 'PRDM1', 'NFKBIA', 'TXNIP', 'RGS1', 'CXCR4', 'IRF1', 'SOCS1', 'FKBP5', 'TNFSF14', 'FASLG', 'KLF10', 'DUSP2', 'DUSP1', 'CBLB', 'CEMIP2', 'HMGB2', 'IL6ST', 'ICAM1', 'ATF3', 'IRF4', 'AREG',
  'PELI1', 'DUSP6', 'CEBPD', 'XAF1', 'NR4A2', 'PMAIP1', 'RGS16', 'FGR', 'CCL3L1',
  'CISH', 'CX3CR1', 'TGFBR1', 'OSM', 'SOCS3', 'CXCL9',
  'FCGR3A', 'CD83', 'SGK1', 'EBI3', 'CSF2', 'MMP2', 'IL10', 'LYZ', 'CXCL2', 'C1QB', 'C1QC', 'CCNA2', 'IFIT2'
  
)

# Create volcano plot
volcano_plot <- ggplot(cluster_deg_df, aes(x = log2fc, y = logP)) +
  geom_point(aes(color = log2fc > threshold_logfc & pval < threshold_pvalue), 
             color = ifelse(cluster_deg_df$log2fc > threshold_logfc & cluster_deg_df$pval < threshold_pvalue, "red", 
                            ifelse(cluster_deg_df$log2fc < -threshold_logfc & cluster_deg_df$pval < threshold_pvalue, "blue", "gray"))) +
  scale_color_manual(values = c("gray", "red", "blue")) +
  labs(title = "Volcano Plot of CD8EX T cell DEGs",
       x = "Log2 Fold Change",
       y = "-log10 p-value") +
  theme_minimal() +
  theme(legend.position = "none") +
  xlim(-2, 2) +
  geom_text_repel(data = cluster_deg_df %>% filter(Gene %in% genes_to_label),
                  aes(label = Gene), size = 3, fontface = "italic", 
                  box.padding = 0.5, point.padding = 0.6,
                  max.overlaps = 20) +  # Control the number of annotations
  
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))  # Add box lines

# Print the plot
print(volcano_plot)
