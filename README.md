# Alpha-synuclein overexpression is associated with epigenomic dysregulation of glutamate signaling and locomotor pathways

This code is for the manuscript [Alpha-synuclein overexpression is associated with epigenomic dysregulation of glutamate signaling and locomotor pathways.](https://www.biorxiv.org/content/10.1101/2021.06.12.448150v1)

Files:

"1-LUHMES_clean_SS_0618020.Rmd" contains pre-processing, normalization, and hmC calculations.
"2a-LUHMES_annotate_all_transcripts.Rmd" annotates all gene transcripts to each EPIC array probe in the final dataset.
"2b-longest_transcript_per_site.Rmd" annotates the longest gene transcript mapping to each EPIC array probe in the dataset.
"2c-LUHMES_genomic_enrichment_background.Rmd" compiled results from scripts 2a and 2b to create a custom background annotation for gene onotlogy enrichment, containing one transcript annotation per CpG site (longest transcript).
"3-LUHMES_linear_model_mC_06182020.Rmd" contains differential methylation analysis for WT aSyn vs control and A30P aSyn vs control comparisons (Figure 2).
"3a-LUHMES_linear_model_mC_A30P_WT_06182020.Rmd" contains differential methylation analysis for A30P aSyn vs WT aSyn comparisons (Figure 2).
"3b-LUHMES_overlap_permute_3way.Rmd" contains permutation testing for overlap in differentially methylated probes between comparisons (Figure 2).
"4-LUHMES_linear_model_hmC_06182020.Rmd" contains differential hydroxymethylation analysis for WT aSyn vs control and A30P aSyn vs control comparisons (Figure 3).
"4a-LUHMES_linear_model_hmC_A30P_WT_06192020.Rmd" contains differential hydroxymethylation analysis for A30P aSyn vs WT aSyn comparisons (Figure 3).
"4b-LUHMES_overlap_permute_3way_hmC_06192020.Rmd" contains permutation testing for overlap in differentially hydroxymethylated probes between comparisons (Figure 3).
"5-LUHMES_overlap_permute_mC_hmC.Rmd" contains permutation testing for overlap in differentially methylated and hydroxymethylated probes (Supplementary Figure 7).
"6-ermineR_uniqueDMgenes_longest_transcript_082020.Rmd" contains gene ontology enrichment analysis, using background from script 2c (Figure 4, Supplementary Figure 10).
"7-Top_probes_tables.Rmd" was used to create tables of the top 20 differentially methylated and hydroxymethylated probes (Tables 1-2, Supplementary Tables 6-8 and 12-14).
"8-TUBA8_multi_omic_plot.Rmd" was used to generate Supplementary Figure 2.
"9-Plotting_loco_hits.Rmd" was used to plot CpGs mapping to locomotory behaviour-related genes (Figure 4).
"9a-Plotting_glutamate_hits.Rmd" was used to plot CpGs mapping to glutamate receptor signaling genes (Supplementary Figure 6).
"10-GRIK2_multi_omic_plot.Rmd" was used to generate Figure 6.
"11-SMITE_WT_092020.Rmd" contains SMITE multi-omic integration for WT aSyn vs control cells (Figures 5-6, Supplementary Figure 8, Supplementary Table 17).
"11a-SMITE_A30P_092020.Rmd" contains SMITE multi-omic integration for A30P aSyn vs control cells (Figures 5-6, Supplementary Figure 9, Supplementary Table 18).
"12-omics_overlap.Rmd" was used to identify overlapping hits from independent differential DNAm, DNAhm, and expression analysis (Supplementary Tables 15-16) and to generate Figure 1.



