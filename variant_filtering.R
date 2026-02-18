{r}

#extracting TFs with with positive pattern tag and qval < 0.05


library(dplyr)
library(rvest)


# Filter hits_unique.tsv files: qval < 0.05 + positive patterns only
# Save each fold/condition separately
# SKIPS if output file already exists


base.data.dir <- /path/to/base/dir
output.dir <- /ourput/directory/
if(!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)

use_folds <- 0:4

for (condition in c("Control", "Treated")){
  cat("\nProcessing", condition, "...\n")
  
  for (fold in use_folds){
    
    # Check if output already exists
    output.file <- paste0(output.dir, "Bcell_", condition, "_", fold, "_hits_FILTERED.tsv")
    if(file.exists(output.file)){
      cat("  Fold", fold, ": already completed, skipping\n")
      next
    }
    
    # Load motifs.html to get significant patterns
    motif.html <- paste0(base.data.dir, condition, "/Bcell_", condition, "_", fold, "_motifs.html")
    
    if(!file.exists(motif.html)){
      cat("  Missing motifs.html for fold", fold, "\n")
      next
    }
    
    motif.info <- as.data.frame(html_table(read_html(motif.html))[[1]])
    
    # Get list of significant positive patterns
    significant_patterns <- motif.info %>% 
      filter(!is.na(qval0), 
             match0 != "None", 
             qval0 < 0.05,
             grepl("^pos_", pattern)) %>%
      pull(pattern)
    
    cat("  Fold", fold, ":", length(significant_patterns), "significant positive patterns\n")
    
    if(length(significant_patterns) == 0) next
    
    # Load hits_unique file
    hits.file <- paste0(base.data.dir, condition, "/Bcell_", condition, "_", fold, "_hits_unique.tsv")
    
    if(!file.exists(hits.file)){
      cat("  Missing hits_unique.tsv for fold", fold, "\n")
      next
    }
    
    cat("  Reading hits file...\n")
    hits.data <- read.table(hits.file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    cat("  Original size:", nrow(hits.data), "rows\n")
    
    # FILTER: Keep only hits from significant positive patterns
    filtered_hits <- hits.data %>%
      filter(motif_name %in% significant_patterns)
    
    cat("  Filtered size:", nrow(filtered_hits), "rows\n")
    cat("  Reduction:", round(100 * (1 - nrow(filtered_hits)/nrow(hits.data)), 1), "%\n")
    
    # Save filtered version
    write.table(filtered_hits, output.file, sep = "\t", row.names = FALSE, quote = FALSE)
    
    # Also save as BED for UCSC upload
    bed.file <- paste0(output.dir, "Bcell_", condition, "_", fold, "_hits_FILTERED.bed")
    bed.data <- filtered_hits %>%
      select(chr, start, end, motif_name, hit_importance, strand)
    write.table(bed.data, bed.file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    cat("  Saved to", output.file, "\n")
  }
}

cat("\nDone. Filtered files saved to:", output.dir, "\n")


## also merges across 3/5 folds


library(dplyr)
library(ggplot2)
library(tidyr)


merged.dir <- /path/to/merged/directory/
plot.dir   <- /directory/for/plot/outputs/
if(!dir.exists(plot.dir)) dir.create(plot.dir, recursive = TRUE)

# LOAD DATA

ctrl.scores  <- read.table("D:/Rotation4_Gaffney/Control_averaged_annot.annotations.tsv",
                           header=TRUE, sep="\t", stringsAsFactors=FALSE)
treat.scores <- read.table("D:/Rotation4_Gaffney/Treated_averaged_annot.annotations.tsv",
                           header=TRUE, sep="\t", stringsAsFactors=FALSE)

control <- read.table(paste0(merged.dir, "Control_merged_consensus_sites.txt"),
                      header=TRUE, sep="\t", stringsAsFactors=FALSE)
treated <- read.table(paste0(merged.dir, "Treated_merged_consensus_sites.txt"),
                      header=TRUE, sep="\t", stringsAsFactors=FALSE)

cat("Control variants:", nrow(ctrl.scores), "\n")
cat("Treated variants:", nrow(treat.scores), "\n")
cat("Control TF sites:", nrow(control), "\n")
cat("Treated TF sites:", nrow(treated), "\n")


# STEP 1: FIND DIFFERENTIALLY BOUND TFs (Control vs Treated)

ctrl_tfs <- control %>%
  group_by(motif_name) %>%
  summarise(n_sites_ctrl = n(),
            mean_imp_ctrl = mean(mean_importance),
            .groups = "drop")

treat_tfs <- treated %>%
  group_by(motif_name) %>%
  summarise(n_sites_treat = n(),
            mean_imp_treat = mean(mean_importance),
            .groups = "drop")

diff_tfs <- inner_join(ctrl_tfs, treat_tfs, by="motif_name") %>%
  mutate(log2fc_sites = log2(n_sites_treat / n_sites_ctrl)) %>%
  arrange(desc(abs(log2fc_sites)))

cat("\nTop differentially bound TFs:\n")
diff_tfs %>% head(20) %>% as.data.frame() %>% print()


# STEP 2: FIND TF SITES NEAR ALL VARIANTS 

find_nearby_tf <- function(variants_df, tf_sites_df, condition_label, window = 500){
  results <- data.frame()
  for(i in 1:nrow(variants_df)){
    var_chr   <- variants_df$chr[i]
    var_pos   <- variants_df$pos[i]
    var_id    <- variants_df$variant_id[i]
    var_logfc <- variants_df$logfc.mean[i]

    nearby <- tf_sites_df %>%
      filter(chr == var_chr,
             start >= (var_pos - window),
             end   <= (var_pos + window)) %>%
      mutate(variant_id          = var_id,
             variant_pos         = var_pos,
             variant_logfc       = var_logfc,
             distance_to_variant = abs(start - var_pos),
             condition           = condition_label)

    results <- rbind(results, nearby)
  }
  return(results)
}

# Use ALL variants, not just in-peak
ctrl.variant.tf  <- find_nearby_tf(ctrl.scores,  control, "Control", window=500)
treat.variant.tf <- find_nearby_tf(treat.scores, treated, "Treated", window=500)
all.variant.tf   <- rbind(ctrl.variant.tf, treat.variant.tf)

cat("\nTotal TF-variant proximities:", nrow(all.variant.tf), "\n")
cat("Unique variants with a nearby TF:", n_distinct(all.variant.tf$variant_id), "\n")


# Join variant-TF proximities with differential TF info
linked <- all.variant.tf %>%
  left_join(diff_tfs, by="motif_name") %>%
  mutate(tf_diff     = abs(log2fc_sites),
         var_effect  = abs(variant_logfc)) %>%
  filter(!is.na(log2fc_sites))

# Summary per variant-TF pair
variant_tf_pairs <- linked %>%
  group_by(variant_id, motif_name, condition) %>%
  summarise(
    variant_pos      = first(variant_pos),
    chr              = first(chr),
    variant_logfc    = first(variant_logfc),
    tf_log2fc_sites  = first(log2fc_sites),
    mean_importance  = mean(mean_importance),
    n_sites          = n(),
    distance         = min(distance_to_variant),
    .groups = "drop"
  ) %>%
  mutate(
    var_effect = abs(variant_logfc),
    tf_effect  = abs(tf_log2fc_sites),
    combined_score = var_effect * tf_effect  # both changing = high score
  ) %>%
  arrange(desc(combined_score))

cat("\nTop variant-TF pairs (variant changing AND TF changing):\n")
variant_tf_pairs %>%
  select(variant_id, chr, variant_pos, motif_name,
         variant_logfc, tf_log2fc_sites, combined_score, distance) %>%
  head(30) %>%
  as.data.frame() %>%
  print()
