library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

# Read in perturbvi enrichment results from ShinyGO
enrich_perturb <- read.csv("PerturbVI/enrichment_pv.csv") %>%
  rename(FDR = `Enrichment.FDR.`,
         size = `nGenes.`,
         enrichmentRatio = `Fold.Enrichment.`,
         description = `Pathway.`) %>%
  mutate(across(where(is.character), ~ gsub("\\\\$", "", .))) %>%
  mutate(FDR = as.numeric(FDR),
         enrichmentRatio = as.numeric(enrichmentRatio),
         size = as.numeric(size),
         method = "Perturb_VI")

sig_enrich_perturb <- enrich_perturb %>%
  filter(FDR <= 0.05)

# Read in GSFA enrichment results from ShinyGO
enrich_gsfa <- read.csv("/Users/camellia/PycharmProjects/perturbvi_analysis/luhmes/luhmes_results/gsfa/gsea/enrichment_gsfa.csv") %>%
  rename(FDR = `Enrichment.FDR.`,
         size = `nGenes.`,
         enrichmentRatio = `Fold.Enrichment.`,
         description = `Pathway.`) %>%
  mutate(across(where(is.character), ~ gsub("\\\\$", "", .))) %>%
  mutate(FDR = as.numeric(FDR),
         enrichmentRatio = as.numeric(enrichmentRatio),
         size = as.numeric(size),
         method = "GSFA")

sig_enrich_gsfa <- enrich_gsfa %>%
  filter(FDR <= 0.05)

# significant enrichment results for both methods
sig_enrich = rbind(sig_enrich_perturb,sig_enrich_gsfa)

# filter on neuron-related pathways
neuron_sig_enrich <- sig_enrich %>%
  filter(grepl("neuron", description, ignore.case = TRUE))

# Get pathway names
neuron_sig_enrich <- neuron_sig_enrich %>%
  separate(description, into = c("GO_ID", "Pathway"), sep = " ", extra = "merge")
sig_enrich_perturb <- sig_enrich_perturb %>%
  separate(description, into = c("GO_ID", "Pathway"), sep = " ", extra = "merge")
neuron_pathway = neuron_sig_enrich %>% select(Pathway) %>% pull() %>% unique()

# Apply abbreviation
abbreviations <- c(
  "negative" = "neg.",
  "regulation" = "reg.",
  "positive" = "pos.",
  "development" = "dev.",
  "projection" = "proj."
)

# Function
categorize_column <- function(data_column, num_categories) {
  # Compute the breaks for the specified number of equal parts
  min_val <- min(data_column, na.rm = TRUE)
  max_val <- max(data_column, na.rm = TRUE)
  breaks <- seq(min_val, max_val, length.out = num_categories + 1) # +1 to include max value
  
  # Round up the left side of each interval to the nearest integer for labels
  rounded_breaks <- ceiling(breaks[-length(breaks)]) # Exclude the last break to avoid duplication
  
  # Convert labels to character for use as factor levels
  labels <- as.character(rounded_breaks)
  
  # Use the computed breaks and labels with the cut function
  categorized_column <- cut(data_column, breaks = breaks, labels = labels, include.lowest = TRUE)
  
  return(categorized_column)
}
#################################Figure about factor 2, 3, 5, 10######################################################
# Extract factor 2, 3, 5, 10
enrich_pv_neuron = 
  sig_enrich_perturb %>% 
  filter(Factor %in% c("2","3","5","10"))%>% 
  filter(FDR < 0.05, enrichmentRatio >= 1) %>%
  filter(Pathway %in% neuron_pathway) %>% 
  mutate(log_fdr = -log10(FDR)) %>% 
  mutate(log_fdr_cat = ifelse(log_fdr<2.5,2,
                              ifelse(log_fdr<4.5,4,
                                     ifelse(log_fdr<6.5,6,
                                            ifelse(log_fdr<8.5,8,10))))) %>% 
  mutate(log_fdr_cat = factor(log_fdr_cat))

# Apply abbreviations to the pathway names
enrich_pv_neuron$Pathway_short <- str_replace_all(enrich_pv_neuron$Pathway, abbreviations)

# short the Pathway
enrich_pv_neuron$Pathway_short <- str_wrap(enrich_pv_neuron$Pathway_short, width = 25)
enrich_pv_neuron$Pathway_short <- factor(enrich_pv_neuron$Pathway_short,
                                         levels = unique(enrich_pv_neuron$Pathway_short))
# categorize number of genes
#categorize column
enrich_pv_neuron$num_gene <- categorize_column(enrich_pv_neuron$size,num_categories = 4)

#function for gsea
plot_gsea = function(df,factor,subset){
  enrich_p = 
    df %>% 
    filter(Factor == factor,Pathway_short %in% subset) %>% 
    ggplot(aes(x = enrichmentRatio, y = reorder(Pathway_short, enrichmentRatio))) +
    geom_point(aes(color = log_fdr_cat, size = num_gene)) +
    scale_color_manual(values =  c("#a1dab4","#41b6c4","#2c7fb8","#253494","#00004d"), name = "-log10(FDR)") +
    geom_segment(aes(x = 0, xend = enrichmentRatio, y = Pathway_short, yend = Pathway_short, color = log_fdr_cat), size = 1) + 
    scale_size_manual(values = seq(2, 4, 0.5), name = "N. Genes") +
    ylab("Pathways") +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text.y = element_text(size = 14),
          axis.text = element_text(size = 14), 
          axis.title=element_text(size=14,face="bold"),
          legend.title = element_text(size=14),
          legend.text = element_text(size=14),
          legend.position = "bottom",
          legend.justification = "center",
          legend.box = "vertical",
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    guides(color = guide_legend(order = 2),
           size = guide_legend(order = 1))
  return(enrich_p)
  
}

# Extract neuron pathway names
get_neuron_pathway_subset <- function(factor, enrich_pv_neuron) {
  subset_df <- enrich_pv_neuron[enrich_pv_neuron$Factor == factor, ]
  
  neuron_pathway_subset <- as.character(subset_df$Pathway_short)
  
  # Remove "\n"
  neuron_pathway_subset <- gsub("\n", " ", neuron_pathway_subset)
  
  # Remove extra whitespace
  neuron_pathway_subset <- trimws(neuron_pathway_subset)
  
  # Dynamically name the subset
  assign(paste0("neuron_pathway_subset_", factor), neuron_pathway_subset, envir = .GlobalEnv)
  
  return(neuron_pathway_subset)
}

# clean up the pathway name
enrich_pv_neuron$Pathway_short <- gsub("\n", " ", enrich_pv_neuron$Pathway_short)
enrich_pv_neuron$Pathway_short <- trimws(enrich_pv_neuron$Pathway_short_clean)

neuron_pathway_subset_2 <- get_neuron_pathway_subset(2, enrich_pv_neuron)
neuron_pathway_subset_3 <- get_neuron_pathway_subset(3, enrich_pv_neuron)
neuron_pathway_subset_5 <- get_neuron_pathway_subset(5, enrich_pv_neuron)
neuron_pathway_subset_10 <- get_neuron_pathway_subset(10, enrich_pv_neuron)

enrich_p2 = plot_gsea(enrich_pv_neuron,"2",neuron_pathway_subset_2)+labs(title = "Factor 2")
enrich_p3 =plot_gsea(enrich_pv_neuron,"3",neuron_pathway_subset_3)+labs(title = "Factor 3")
enrich_p5 =plot_gsea(enrich_pv_neuron,"5",neuron_pathway_subset_5)+labs(title = "Factor 5")

legend = get_plot_component(enrich_p7, "guide-box",return_all = TRUE)[[3]]

