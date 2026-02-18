{r}

### Differential Motif analysis using modisco results
##aggregating to family-level rather than individual for TF families


# Install packages if needed
packages <- c("tidyr", "rvest", "dplyr")
missing <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(missing)) install.packages(missing)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

bioc_packages <- c("universalmotif", "ComplexHeatmap")
missing_bioc <- bioc_packages[!(bioc_packages %in% installed.packages()[,"Package"])]
if(length(missing_bioc)) BiocManager::install(missing_bioc)

# Load libraries
require("universalmotif")
require("tidyr")
require("ComplexHeatmap")
require("rvest")
require("dplyr")
require("circlize")


library(dplyr)
library(tidyr)
library(ggplot2)
library(rvest)
library(universalmotif)



# Set paths
base.data.dir <- /change/path/here/
output.dir <- add/path/here
if(!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)

# Load JASPAR motif families
transfac.path <- /path/to/transfac/family/names/
transfac.motifs <- read_transfac(transfac.path)

family.table <- data.frame(
  motif.id = unlist(lapply(transfac.motifs, function(x){x@name})),
  motif.name = unlist(lapply(transfac.motifs, function(x){x@altname})),
  motif.family = gsub("tf_family\\:","", unlist(lapply(transfac.motifs, function(x){x@family})))
)

family.table$motif.family[which(family.table$motif.family == "")] <- 
  family.table$motif.name[which(family.table$motif.family == "")]

# ADD BASE ID for version matching
family.table$motif.id.base <- gsub("\\.\\d+$", "", family.table$motif.id)
family.table <- family.table %>% distinct()

cat("Loaded", length(unique(family.table$motif.family)), "TF families\n")

# Extract motif data
use_folds <- 0:4
total.peak.data <- data.frame()

for (condition in c("Control", "Treated")){
  cat("\nProcessing", condition, "...\n")
  
  for (fold in use_folds){
    motif.html <- paste0(base.data.dir, condition, "/Bcell_", condition, "_", fold, "_motifs.html")
    
    if(!file.exists(motif.html)){
      cat("  Missing motifs.html for fold", fold, "\n")
      next
    }
    
    motif.info <- as.data.frame(html_table(read_html(motif.html))[[1]])
    
    # JOIN ON BASE ID + QVAL FILTER
    motif.tally <- motif.info %>% 
      filter(!is.na(qval0), match0 != "None", qval0 <= 0.05) %>%
      mutate(
        scaled.seqlets = num_seqlets/sum(num_seqlets),
        match0.base = gsub("\\.\\d+$", "", match0)
      ) %>%
      left_join(family.table, by = join_by(match0.base == motif.id.base)) %>%
      filter(!is.na(motif.family)) %>%
      select(scaled.seqlets, match0, motif.name, motif.family, pattern, qval0)
    
    if(nrow(motif.tally) > 0) {
      motif.tally$condition <- condition
      motif.tally$fold <- fold
      total.peak.data <- rbind(total.peak.data, motif.tally)
      cat("  Fold", fold, ":", nrow(motif.tally), "motifs (qval < 0.05)\n")
    }
  }
}

cat("\nTotal instances:", nrow(total.peak.data), "\n")
cat("Unique TF families:", length(unique(total.peak.data$motif.family)), "\n")

# Aggregate by family
family.per.fold <- total.peak.data %>%
  group_by(motif.family, condition, fold) %>%
  summarise(total_scaled_seqlets = sum(scaled.seqlets), .groups = 'drop')

family.summary <- family.per.fold %>%
  group_by(motif.family, condition) %>%
  summarise(
    n_folds_detected = n(),
    mean_scaled_seqlets = mean(total_scaled_seqlets),
    sd_scaled_seqlets = sd(total_scaled_seqlets),
    .groups = 'drop'
  )

# Create comparison - CHANGED TO "Stimulated"
family.comparison <- family.summary %>%
  pivot_wider(
    id_cols = motif.family,
    names_from = condition,
    values_from = c(n_folds_detected, mean_scaled_seqlets, sd_scaled_seqlets),
    values_fill = list(n_folds_detected = 0, mean_scaled_seqlets = 0, sd_scaled_seqlets = 0)
  ) %>%
  mutate(
    robust = (n_folds_detected_Control >= 3 | n_folds_detected_Treated >= 3),
    log2FC = log2((mean_scaled_seqlets_Treated + 0.001) / 
                  (mean_scaled_seqlets_Control + 0.001)),
    category = case_when(
      log2FC > 0 ~ "Stimulated",  # change HERE
      log2FC < 0 ~ "Suppressed",  #and here
      TRUE ~ "Unchanged"
    )
  ) %>%
  filter(robust, abs(log2FC) > 0.5)

cat("\nFamilies with |log2FC| > 0.5:", nrow(family.comparison), "\n")
cat("Stimulated:", sum(family.comparison$category == "Stimulated"), "\n")
cat("Suppressed:", sum(family.comparison$category == "Suppressed"), "\n")

# Print the table
print(family.comparison %>% arrange(desc(abs(log2FC))))

# Create plot and assign to p
p <- family.comparison %>%
  arrange(log2FC) %>%
  mutate(motif.family = factor(motif.family, levels = motif.family)) %>%
  
  ggplot(aes(x = motif.family, y = log2FC, fill = category)) +
  geom_col(width = 0.75) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  coord_flip() +
  
  scale_fill_manual(
    values = c("Stimulated" = "#9B1B30",
               "Suppressed" = "skyblue"),
    name = ""
  ) +
  
  labs(
    x = "",
    y = "log2(stimulated/control)",
    title = "Differential TF Binding in..."
  ) +
  
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "top"
  )
p
# # Now save it
#ggsave(paste0(output.dir, "qfiltered_log2FC_final.png"), plot = p, width = 6, height = 6)

