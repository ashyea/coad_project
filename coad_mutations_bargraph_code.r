#This code can be used to make a stacked bar graph comparing mutation frequencies for any studies. 

#install the mentioned packages beforehand
library(cbioportalR)
library(dplyr)
library(ggplot2)
library(tidyr)

set_cbioportal_db("https://www.cbioportal.org")

# studies for data retrieval
studies <- c(
  "coadread_tcga",
  "coadread_tcga_pan_can_atlas_2018",
  "coadread_tcga_pub"
)

#data retrieval
all_mutations <- list()

for(study in studies) {
  cat("Processing:", study, "\n")
  
  tryCatch({
    study_data <- get_genetics_by_study(study_id = study)
    mutations <- study_data$mutation #extracting mutation data
    mutations$study <- study #study identifier for mutation data
    all_mutations[[study]] <- mutations
    cat("  Retrieved", nrow(mutations), "mutations\n")
    
  }, error = function(e) {
    cat("  Error for", study, ":", e$message, "\n")
  })
  
  Sys.sleep(0.3)
}

# combining all data in one list
combined_data <- bind_rows(all_mutations)

# verification
cat("\n=== DATA SUMMARY ===\n")
cat("Studies retrieved:", paste(unique(combined_data$study), collapse = ", "), "\n")
cat("Total mutation records:", nrow(combined_data), "\n")
cat("Unique patients:", n_distinct(combined_data$patientId), "\n\n")

# table of major cancer genes
gene_table <- table(combined_data$hugoGeneSymbol)
top_genes <- names(sort(gene_table, decreasing = TRUE))[1:5] #the ranking of top genes can be adjusted
cat("Top 5 genes:", paste(top_genes, collapse = ", "), "\n") 

# calculating mutation frequencies

# one patient per altered gene: for frequency, no. of genes affected, not no. of mutations
mutation_counts <- combined_data %>%
  filter(hugoGeneSymbol %in% top_genes) %>%
  distinct(study, patientId, hugoGeneSymbol) %>%  # one count per gene in unique patient
  group_by(study, hugoGeneSymbol) %>% 
  summarise(mutated_patients = n(), .groups = "drop") 

# total patients per study
patient_counts <- combined_data %>%
  group_by(study) %>%
  summarise(total_patients = n_distinct(patientId))

# frequency table
# creating all possible study-gene combinations
all_combinations <- expand.grid(
  study = unique(combined_data$study),
  hugoGeneSymbol = top_genes,
  stringsAsFactors = FALSE
)

# merge with calculated mutation counts 	
mutation_freq <- all_combinations %>%
  left_join(mutation_counts, by = c("study", "hugoGeneSymbol")) %>%
  left_join(patient_counts, by = "study") %>%
  mutate(
    mutated_patients = ifelse(is.na(mutated_patients), 0, mutated_patients), #converting any non-tested mutations to 0 for ease of calculation
    frequency = mutated_patients / total_patients * 100
  )

# checking the frequency table
print("Mutation frequency table (first 10 rows):")
print(head(mutation_freq, 10))

# bar graph
ggplot(mutation_freq, aes(x = study, y = frequency, fill = hugoGeneSymbol)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  labs(
    title = "Mutation Frequency Across Colorectal Cancer Studies",
    subtitle = "Stacked Bar Graph Visualising the Frequency of Top 5 Mutated Genes",
    x = "Study",
    y = "Mutation Frequency (%)",
    fill = "Gene",
    caption = paste("Data from", n_distinct(combined_data$study), "studies,",
                   n_distinct(combined_data$patientId), "patients total")
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12),
    legend.position = "right",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 11, face = "bold"),
    plot.caption = element_text(size = 9, color = "gray40")
  ) +
  scale_fill_viridis_d(option = "rainbow", name = "Gene") +
  
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())


ggsave("mutation_frequency_stacked_bar.png", width = 12, height = 8, dpi = 300)
