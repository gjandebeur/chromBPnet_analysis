# chromBPnet_analysis
This repository contains scripts for training ChromBPNet models on B cell ATAC-seq data (Unstimulated vs Stimulated conditions) and performing variant effect prediction at lupus-associated genetic loci.



Tested across 1, 40, and 80 individual groups

<img width="2400" height="1200" alt="chrombpnet_model_performance" src="https://github.com/user-attachments/assets/aa7b121f-91da-4c12-b38e-b0a87444c933" />


## Prerequisites
-Clone repository and activate in ChromBPnet environment

```
git clone https://github.com/gjandebeur/chromBPnet_analysis.git
cd chromBPnet_analysis
conda activate chrombpnet ### Will need to set up the conda environment first.
``` 

To run this pipeline, you will need to obtain:
  -ATAC-seq BAMs
  -reference FASTA
  -chromosome sizes
  -blacklist peaks
  -preprocessing files (folds, peak calling, etc)


### GETTING STARTED 

Step 1: Start by using ChromBPnet to bias train out the Tn5 motif insertions from its models
update bias_train script to fix the necessary directories and run

*run bias_train.sh script*

Step 2: After removing Tn5 bias, the model must be trained on data to predict chromatin accessibility from DNA sequence.

*run train.sh script*
-This produces one model per fold and per condition, so say you have 2 conditions, expect 10 models back (5 folds)

Step 3: Produce contribution bigwig score files
-compute per-base contribution score across the genome to reflect how much each nucleotide contributes to predicted accessibility and further downstream analysis.

*run Contrib_bw script*

Step 4: Using the model and predicted bigwigs, produce transcription factor pattern analysis to identify known TF motifs across the genome.

*run TF-MoDISco.sh script*

Step 5: Perform variant effect prediction using the trained models from above
-score fine-mapped variants from list by computing predicted accesibility between reference and alternate alleles.

*run ChromBPnet_variants.sh script first*
*then run variant_filtering.R*
-the R script filters and annotates scored variants, identifying major changes between conditions.


Step 6: Differential motif analysis
-Identify TF motifs with differential activity between conditions to identify changes TFs and/or variants

*run Differential_motif.R*


## Citation
I am not the original creator of ChromBPnet and only using it for my applications. If you use this pipeline, please consider citing the original authors.
```
Pampari et al. (2023). chrombpnet: bias factorized, base-resolution deep learning models of chromatin accessibility reveal cis-regulatory sequences in human diseases. bioRxiv. https://doi.org/10.1101/2023.10.19.563323
```

