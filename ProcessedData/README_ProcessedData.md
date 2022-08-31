Please download data files from Zenodo and save them in this folder.

https://doi.org/10.5281/zenodo.6993160

#######################################################
# Folder metabolomics
- metabolites_allions_combined_formulas_with_metabolite_filters.csv - all ions stably detected across tissues with different methods with annotations and filter values for annotated metabolites.
- metabolites_allions_combined_formulas_with_metabolite_filters_spatial100clusters.csv - all ions with information on which spatial cluster they belong to. 
- metabolites_allions_combined_formulas_with_metabolite_filters_spatial100clusters_HMDBsubclass.csv - all ions with spatial clusters and HMDB subclass info. 
- metabolites_allions_combined_formulas_with_metabolite_filters_spatial100clusters_with_mean.csv - all ions with spatial clusters and mean values per tissue.
- metabolites_allions_combined_norm_intensity.csv - normalized intensity per sample for all ions. 
- table_diff_abundance_metabolite_ions.csv - differential abundance of all ions with fold change, p-values (t-test) and FDR. 
- table_diff_abundance_metabolite_ions_removed2outliers.csv - differential abundance of all ions with fold change, p-values (t-test) and FDR, up to 2 outliers were removed from 5 replicates efore differential analysis in case they were not detected. 

#######################################################
# Folder sequencing
- countsMatrix_annTable_filtered.csv - annotated and filtered list of genes detected with metagenomics and metatranscriptomics.
- countsMatrixGetMM_DNA_table.txt - GetMM-normalized metagenomics data.
- countsMatrixGetMM_RNA_table.txt - GetMM-normalized metatranscriptomics data. 
- countsMatrixRAW_DNA_table.txt - raw counts for metagenomics data. 
- countsMatrixRAW_DNA_table_with_species.txt - raw counts for metagenomics data with species column.
- countsMatrixRAW_RNA_table.txt - raw counts for metatranscriptomics data.
- countsMatrixRAW_RNA_table_with_species.txt - raw counts for metatranscriptomics data with species column.
- edgeR_gene_fold_changes_and_ann.csv - fold changes, p-values and FDR for EdgeR calculated differential expression of genes.
- edgeR_gene_fold_changes_and_ann_filtered_FDRrecalculated.csv - fold changes, p-values and FDR for EdgeR calculated differential expression of genes, FDR recalculated for the entire table.
- gene_annotation_COG.tsv - gene-COG insidence table extracted from EggNOG genome annotations. 
- gene_annotation_EC.tsv - gene-EC number insidence table extracted from EggNOG genome annotations. 
- gene_annotation_EggNOG.tsv - gene-EggNOG insidence table extracted from EggNOG genome annotations. 
- gene_annotation_go.tsv - gene-GO term insidence table extracted from EggNOG genome annotations. 
- gene_annotation_kegg.tsv - gene-KEGG pathway insidence table extracted from EggNOG genome annotations. 
- geneAnnTable_full.csv - full gene annotation table. 
- merged_humann2_metaphlan_bugs_list_mouseDNA_OTUs.txt - MetaPhlan2 output - list of detected species in metagenomics samples.
- merged_humann2_metaphlan_species_filtered_mouseDNA_OTUs.txt - filtered MetaPhlan2 output (only species level, detected in >5 samples > 0.01%). 
- merged_humann2_metaphlan_species_filtered_mouseDNA_OTUs_with_FC.txt - fold changes of MetaPhlan2 species abundances. 

#######################################################
# Folder util
- hmdbPTWtables.mat - HMDB pathway annotations of metabolites. For each hmdb class, subclass and superclass, a binary table for all stably detected ions of whether they belong to thie pathway or not. 
- hmdbV4_072021_ptw.mat - HMDB pathway annotations (downloaded from HMDB) in matlab structure format. This is used to create hmdbPTWtables in script 'create_hmdb_matrices.m'.
- keggreactionsALLreactantsenzymes.mat - structure of kegg reactions with corresponding enzymes and substrate-product pairs. This file can be created with scripts 'get_kegg_reaction_list.m' to get the list of KEGG reactions, followed by 'get_kegg_reaction_reactants.m' to get ec numbers and reactant of each reaction. 

#######################################################
# Folder example_output - example output files produces by defferent analysis scripts
- cgo_clustergrams_of_model_coefficients.mat - matlab file containing clustergram of model coefficients and manually defined sub-clustergrams.
- kegg_rxn_ec_subprod.mat - matlab file containing kegg reactions with enzymes, substrates and products as a table.
- kmeans_clustering_GIT100_CTR_DC.txt - results of k-means clustering of spatial metabolite profiles for CTR DC group.
- kmeans_clustering_GIT100_CTR_GF.txt - results of k-means clustering of spatial metabolite profiles for CTR GF group.
- kmeans_clustering_GIT100_HFD_DC.txt - results of k-means clustering of spatial metabolite profiles for HFD DC group.
- kmeans_clustering_GIT100_HFD_GF.txt - results of k-means clustering of spatial metabolite profiles for HFD GF group.
- metabolites_allions_combined_formulas.csv - combined ions detected across all tissues and methods.
- metabolites_allions_combined_formulas_with_metabolite_filters.csv - combined ions with metabolite annotations and filter information. 
- metabolites_allions_combined_norm_intensity.csv - normalized intensity for combined ions. 
- model_results_SMOOTH_normbyabsmax_2LIcoefHost1LIcoefbact_allions.csv - results of the intestinal flux model parameter fitting, parameters normalized by maximum absolute value. 
- model_results_SMOOTH_normbyabsmax_ONLYMETCOEF_2LIcoefHost1LIcoefbact_allions.csv - results of the intestinal flux model parameter fitting, parameters normalized by maximum absolute value, only metabolix flux parameters (intesinal flux omitted). 
- model_results_SMOOTH_normbyabsmax_reciprocal_problem_allions.csv - results of the reciprocal problem for the intestinal flux model parameter fitting. 
- model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_allions.csv - results of the intestinal flux model parameter fitting, original parameter values. 
- ptwenr_HMDBclass_hcluster_2LIhos1LIbact_hmdbv4.csv - HMDB chemical class enrichment analysis for groups from the hierarchical clustering of model coefficients.
- ptwenr_HMDBsubclass_hcluster_2LIhos1LIbact_hmdbv4.csv - HMDB chemical subclass enrichment analysis for groups from the hierarchical clustering of model coefficients.
- ptwenr_HMDBsuperclass_hcluster_2LIhos1LIbact_hmdbv4.csv - HMDB chemical superclass enrichment analysis for groups from the hierarchical clustering of model coefficients.
- ptwenr_HMDBclass_kmeans_clustering_spatial_all_hmdbv4.csv - HMDB chemical class enrichment analysis for clusters from k-means clustering of spatial metabolite profiles. 
- ptwenr_HMDBsubclass_kmeans_clustering_spatial_all_hmdbv4.csv - HMDB chemical subclass enrichment analysis for clusters from k-means clustering of spatial metabolite profiles. 
- ptwenr_HMDBsuperclass_kmeans_clustering_spatial_all_hmdbv4.csv - HMDB chemical superclass enrichment analysis for clusters from k-means clustering of spatial metabolite profiles. 
- ptwenr_recalc_COG_DOWN_updfiltered_eggnog_per_species.csv - COG pathway enrichment for downregulated transcripts per species. 
- ptwenr_recalc_COG_DOWN_updfiltered_eggnog_total.csv - COG pathway enrichment for downregulated transcripts for all detected transcripts.
- ptwenr_recalc_COG_UP_updfiltered_eggnog_per_species.csv - COG pathway enrichment for upregulated transcripts per species.
- ptwenr_recalc_COG_UP_updfiltered_eggnog_total.csv - COG pathway enrichment for upregulated transcripts for all detected transcripts.
- ptwenr_recalc_EC_DOWN_updfiltered_eggnog_per_species.csv - EC pathway enrichment for downregulated transcripts per species.
- ptwenr_recalc_EC_DOWN_updfiltered_eggnog_total.csv - EC pathway enrichment for downregulated transcripts for all detected transcripts.
- ptwenr_recalc_EC_UP_updfiltered_eggnog_per_species.csv - EC pathway enrichment for upregulated transcripts per species.
- ptwenr_recalc_EC_UP_updfiltered_eggnog_total.csv - EC pathway enrichment for upregulated transcripts for all detected transcripts.
- ptwenr_recalc_KEGG_DOWN_updfiltered_eggnog_per_species.csv - KEGG pathway enrichment for downregulated transcripts per species.
- ptwenr_recalc_KEGG_DOWN_updfiltered_eggnog_total.csv - KEGG pathway enrichment for downregulated transcripts for all detected transcripts.
- ptwenr_recalc_KEGG_UP_updfiltered_eggnog_per_species.csv - KEGG pathway enrichment for upregulated transcripts per species.
- ptwenr_recalc_KEGG_UP_updfiltered_eggnog_total.csv - KEGG pathway enrichment for upregulated transcripts for all detected transcripts.
- ptwNchangingGenes_COG_DOWN_updfiltered_eggnog.csv - number of changing genes (downregulated) per COG pathway.
- ptwNchangingGenes_COG_UP_updfiltered_eggnog.csv - number of changing genes (upregulated) per COG pathway.
- ptwNchangingGenes_EC_DOWN_updfiltered_eggnog.csv - number of changing genes (downregulated) per EC category.
- ptwNchangingGenes_EC_UP_updfiltered_eggnog.csv - number of changing genes (upregulated) per EC category.
- ptwNchangingGenes_KEGG_DOWN_updfiltered_eggnog.csv - number of changing genes (downregulated) per KEGG pathway.
- ptwNchangingGenes_KEGG_UP_updfiltered_eggnog.csv - number of changing genes (upregulated) per KEGG pathway.
- shortest_paths_upto_length_4.csv - list of shortest enzymatic paths up to length 4 (not including 4) for potential bacterial substrates and products. 
- table_diff_abundance_metabolite_ions_removed2outliers.csv - differential abundance of metabolites (fold change, p-value, FDR), up to 2 non-detected outliers removed from 5 replicates before analysis. 
- table_hierarchical_clustering_groups.csv - table with hierarchical clustering group attributions for all analyzed metabolites. 
- table_kegg_DNA_sub_prod_products_bestcorr_EC_pos.csv - EC numbers for genes correlating best with potential bacterial products (metagenomics).
- table_kegg_DNA_sub_prod_products_bestcorr_genes_pos.csv - genes correlating best with potential bacterial products (metagenomics).
- table_kegg_DNA_sub_prod_products_bestcorrP_pos.csv - p-values for genes correlating best with potential bacterial products (metagenomics).
- table_kegg_DNA_sub_prod_products_bestcorr_pos.csv - Pearson's correlation values for genes correlating best with potential bacterial products (metagenomics).
- table_kegg_RNA_sub_prod_products_bestcorr_EC_pos.csv - EC numbers for genes correlating best with potential bacterial products (metatranscriptomics).
- table_kegg_RNA_sub_prod_products_bestcorr_genes_pos.csv - genes correlating best with potential bacterial products (metatranscriptomics).
- table_kegg_RNA_sub_prod_products_bestcorr_pos.csv - p-values for genes correlating best with potential bacterial products (metatranscriptomics).
- table_kegg_RNA_sub_prod_products_bestcorrP_pos.csv - Pearson's correlation values for genes correlating best with potential bacterial products (metatranscriptomics).
- table_kegg_rxn_ec_substrate_product.csv - kegg reactions with EC numbers and substrates and products as a table. 
- table_potential_products_ids.csv - IDs of potential bacterial products identified with intestinal flux model. 
- table_potential_substrates_ids.csv - IDs of potential bacterial substrates identified with intestinal flux model. 
- table_rxnCompounds.csv - KEGG IDs of reactants of the reactions. 
- table_shortest_path_subprod_sorted_unique.csv - table of shortest paths between potential microbial substrates and products (unique). 
- table_shortestMatrix_paths_filtered.csv - table of all shortest paths between potential microbial substrates and products. 
- table_speciesOTU_products_corr.csv - p-values values for OTUs correlating with potential bacterial products (MetaPhlan2 on metagenomics).
- table_speciesOTU_products_corrP.csv - Pearson's correlation values for OTUs correlating with potential bacterial products (MetaPhlan2 on metagenomics).
#################################################################