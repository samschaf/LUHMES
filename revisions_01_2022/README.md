
# Alpha-synuclein induces epigenomic dysregulation of glutamate signaling and locomotor pathways

This code is for the manuscript [Alpha-synuclein induces epigenomic dysregulation of glutamate signaling and locomotor pathways.](https://www.biorxiv.org/content/10.1101/2021.06.12.448150v1)
This folder contains scripts added and/or modified during the revision stage (January 2022).

Files:

"2c-LUHMES_genomic_enrichment_background_12102021update.Rmd" was updated to add enhancers (H3K4me1 LUHMES ChIP-seq) to the background annotation assigning one gene feature per EPIC probe.

"3-LUHMES_linear_model_mC_12102021.Rmd" was updated to permute differentially methylated probes against the background with annotated enhancers (Figure 2).

"3c-LUHMES_CoMeBack_mC.Rmd" contains definition of custom DNAm CMRs, linear modeling for differential CMR methylation, permutation of gene features, and overlap with site-specific analysis (Supplementary Figure 3).

"3c-LUHMES_CoMeBack_hmC.Rmd" contains definition of custom DNAhm CMRs, linear modeling for differential CMR hydroxymethylation, permutation of gene features, and overlap with site-specific analysis (Supplementary Figure 4).

"4-LUHMES_linear_model_hmC_12142021.Rmd" was updated to permute differentially hydroxymethylated probes against the background with annotated enhancers (Figure 3).

"5-LUHMES_overlap_permute_hmC_mC.Rmd" contains code to generate the new divergent heat map in the previous Supplementary Figure 4, now numbered as Supplementary Figure 6.

"5b-LUHMES_overlap_permute_mC_hmC_mRNA.Rmd" contains permutation testing for overlap in differentially methylated probes, differentially hydroxymethylated probes, and differentially expressed genes.

"11-SMITE_WT_092020_12082021update.Rmd" contains additional SMITE models for the control vs. WT aSyn analysis with various weights applied.

"11-SMITE_A30P_092020_12102021update.Rmd" contains additional SMITE models for the control vs. A30P aSyn analysis with various weights applied.

"12-omics_overlap_12082021update.Rmd" was updated to add transcript-specific differential expression for SNCA and to update genome browser plots generated for Figure 1.

"13-LUHMES_vs_BAC_SNCA_mice.Rmd" identifies overlapping differentially methylated and hydroxymethylated genes in control vs. WT aSyn LUHMES and hippocampus of 12-month old control vs. WT aSyn BAC mice.
