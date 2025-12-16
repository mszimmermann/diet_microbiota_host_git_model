% workflow merge rnaseq tables merges raw, normalized and diff abundance
% table in one with gene annotations
% call script defining file dependenciesand global variables
addpath(genpath('.\'));
add_global_and_file_dependencies

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Required files:
% 'edgeR_deseq2_gene_fold_changes_and_ann_filtered_FDRrecalculated.csv'
% 'countsMatrixGetMM_DNA_table.txt'
% 'countsMatrixGetMM_RNA_table.txt'
% 'countsMatrixRAW_DNA_table.txt'
% 'countsMatrixRAW_RNA_table.txt'
% 'geneAnnTable_filtered.csv'
% 'geneAnnTable_full.csv'
% Output: 
% Files:
% 'table_S4_merged_geneAnnTable_with_raw_norm_and_diff_counts.csv'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% also merge pathway enrichment results
% Required files:
% 'ptwNchangingGenes_KEGG_DOWN_updfiltered_eggnog_edgerdeseq_FDR005_FC05.csv'
% 'ptwNchangingGenes_KEGG_UP_updfiltered_eggnog_edgerdeseq_FDR005_FC05.csv'
% 'ptwenr_recalc_KEGG_DOWN_updfiltered_eggnog_edgerdeseq_FDR005_FC05_per_species.csv'
% 'ptwenr_recalc_KEGG_DOWN_updfiltered_eggnog_edgerdeseq_FDR005_FC05_total.csv'
% 'ptwenr_recalc_KEGG_UP_updfiltered_eggnog_edgerdeseq_FDR005_FC05_per_species.csv'
% 'ptwenr_recalc_KEGG_UP_updfiltered_eggnog_edgerdeseq_FDR005_FC05_total.csv'
% 'kegg_ptw_names_and_bact_flag.csv'
% Output: 
% Files:
% 'table_S5_merged_ptwenrTable.csv'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read gene annotation and counts tables
% get full getMM normalization file and separate by species 
% detect import options to import columns with the same type
opts = detectImportOptions([outputFolder,...
                            'geneAnnTable_full.csv']);
% manual curation of variable types
opts.VariableNamesLine=1;
opts.DataLine = 2;
opts.VariableTypes(26:end) = {'char'};
% read table
geneAnn = readtable([outputFolder,...
    'geneAnnTable_full.csv'], ...
     opts,...
     'ReadVariableNames', 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read count tables
countsMatrixGetMM_RNA_table = readtable([outputFolder,...
    'countsMatrixGetMM_RNA_table.txt']);
countsMatrixRAW_RNA_table = readtable([outputFolder,...
    'countsMatrixRAW_RNA_table.txt']);
% load gene abundance from DNA and correlate with metabolite abundances
countsMatrixGetMM_DNA_table = readtable([outputFolder,...
    'countsMatrixGetMM_DNA_table.txt']);
countsMatrixRAW_DNA_table = readtable([outputFolder,...
    'countsMatrixRAW_DNA_table.txt']);

% read updated table with edgeR and deseq2 results instead of edgeR_gene_fold_changes_and_ann.csv
edgeRTable = readtable([outputFolder, ...
                    'edgeR_deseq2_gene_fold_changes_and_ann_filtered_FDRrecalculated.csv']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% edit column names of the count tables
% metaG raw
sample_columns = cellfun(@(x) contains(x, 'MSZ'), countsMatrixRAW_DNA_table.Properties.VariableNames);
countsMatrixRAW_DNA_table.Properties.VariableNames(sample_columns) = ...
    strcat(countsMatrixRAW_DNA_table.Properties.VariableNames(sample_columns),...
          '_metaG_raw');
% make gene filter the first column
countsMatrixRAW_DNA_table = movevars(countsMatrixRAW_DNA_table,'geneFilter','Before',1);
% metaG getMM
sample_columns = cellfun(@(x) contains(x, 'MSZ'), countsMatrixGetMM_DNA_table.Properties.VariableNames);
countsMatrixGetMM_DNA_table.Properties.VariableNames(sample_columns) = ...
    strcat(countsMatrixGetMM_DNA_table.Properties.VariableNames(sample_columns),...
          '_metaG_norm_getmm');

% metaT raw
sample_columns = cellfun(@(x) contains(x, 'MSZ'), countsMatrixRAW_RNA_table.Properties.VariableNames);
countsMatrixRAW_RNA_table.Properties.VariableNames(sample_columns) = ...
    strcat(countsMatrixRAW_RNA_table.Properties.VariableNames(sample_columns),...
          '_metaT_raw');

% metaG getMM
sample_columns = cellfun(@(x) contains(x, 'MSZ'), countsMatrixGetMM_RNA_table.Properties.VariableNames);
countsMatrixGetMM_RNA_table.Properties.VariableNames(sample_columns) = ...
    strcat(countsMatrixGetMM_RNA_table.Properties.VariableNames(sample_columns),...
          '_metaT_norm_getmm');
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% in DEG table, leave only relevant columns and species_genes and indexFiltered
keep_columns = cellfun(@(x) (contains(x, 'HFDCTR') |...
                             contains(x, 'species_genes') |...
                             contains(x, 'indexFiltered')),...
                             edgeRTable.Properties.VariableNames);
edgeRTable = edgeRTable(:, keep_columns);
% rename columns to indicate method
for i=1:length(edgeRTable.Properties.VariableNames)
    edgeRTable.Properties.VariableNames{i} = strrep(edgeRTable.Properties.VariableNames{i},...
        'fcHFDCTR_DNA', 'fc_EdgeR_HFDCTR_DNA');
    edgeRTable.Properties.VariableNames{i} = strrep(edgeRTable.Properties.VariableNames{i},...
        'pHFDCTR_DNA', 'pval_EdgeR_HFDCTR_DNA');
    edgeRTable.Properties.VariableNames{i} = strrep(edgeRTable.Properties.VariableNames{i},...
        'fdrHFDCTR_DNA', 'fdr_EdgeR_HFDCTR_DNA');
    edgeRTable.Properties.VariableNames{i} = strrep(edgeRTable.Properties.VariableNames{i},...
        'fcHFDCTR_RNA', 'fc_EdgeR_HFDCTR_RNA');
    edgeRTable.Properties.VariableNames{i} = strrep(edgeRTable.Properties.VariableNames{i},...
        'pHFDCTR_RNA', 'pval_EdgeR_HFDCTR_RNA');
    edgeRTable.Properties.VariableNames{i} = strrep(edgeRTable.Properties.VariableNames{i},...
        'fdrHFDCTR_RNA', 'fdr_EdgeR_HFDCTR_RNA');
    % edit deseq columns to match the pattern of edger columns
    edgeRTable.Properties.VariableNames{i} = strrep(edgeRTable.Properties.VariableNames{i},...
        'fcDeseq2HFDCTR_DNA', 'fc_Deseq2_HFDCTR_DNA');
    edgeRTable.Properties.VariableNames{i} = strrep(edgeRTable.Properties.VariableNames{i},...
        'pDeseq2HFDCTR_DNA', 'pval_Deseq2EdgeR_HFDCTR_DNA');
    edgeRTable.Properties.VariableNames{i} = strrep(edgeRTable.Properties.VariableNames{i},...
        'fdrDeseq2HFDCTR_DNA', 'fdr_Deseq2_HFDCTR_DNA');
    edgeRTable.Properties.VariableNames{i} = strrep(edgeRTable.Properties.VariableNames{i},...
        'fcDeseq2HFDCTR_RNA', 'fc_Deseq2_HFDCTR_RNA');
    edgeRTable.Properties.VariableNames{i} = strrep(edgeRTable.Properties.VariableNames{i},...
        'pDeseq2HFDCTR_RNA', 'pval_Deseq2_HFDCTR_RNA');
    edgeRTable.Properties.VariableNames{i} = strrep(edgeRTable.Properties.VariableNames{i},...
        'fdrDeseq2HFDCTR_RNA', 'fdr_Deseq2_HFDCTR_RNA');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%merge tables
mergedTable = join(geneAnn, countsMatrixRAW_DNA_table, 'keys', 'species_genes');
mergedTable = join(mergedTable, countsMatrixGetMM_DNA_table, 'keys', 'species_genes');
mergedTable = join(mergedTable, countsMatrixRAW_RNA_table, 'keys', 'species_genes');
mergedTable = join(mergedTable, countsMatrixGetMM_RNA_table, 'keys', 'species_genes');
mergedTable = outerjoin(mergedTable, edgeRTable, 'keys', 'species_genes');

% clean up the table:
% remove species_genes from edgeRTable
mergedTable(:,ismember(mergedTable.Properties.VariableNames, 'species_genes_edgeRTable')) = [];
% get geneFilter columns
gene_filter_columns = cellfun(@(x) contains(x, 'geneFilter'), ...
                        mergedTable.Properties.VariableNames);
% test that all gene filters are the same
test_geneFilters = mergedTable{:, gene_filter_columns};
flag_mismatch = 0;
for i=2:size(test_geneFilters,2)
    if nnz(test_geneFilters(:,i)~=test_geneFilters(:,1))>0
        disp "Error: table has non-matching gene filter values"
        flag_mismatch=1;
    end
end
if flag_mismatch==0
    % remove all gene filter columns but the first one
    delete_genefilters = find(gene_filter_columns);
    delete_genefilters = delete_genefilters(2:end);
    mergedTable(:, delete_genefilters) = [];
    % rename geneFilter column
    gene_filter_columns = cellfun(@(x) contains(x, 'geneFilter'), ...
                        mergedTable.Properties.VariableNames);
    mergedTable.Properties.VariableNames{gene_filter_columns} = 'geneFilter';
end

% put indexFiltered after geneFilter columns
mergedTable = movevars(mergedTable,'indexFiltered', 'After', 'geneFilter');

% save merged table to file
writetable(mergedTable, [outputFolder,...
    'table_S4_merged_geneAnnTable_with_raw_norm_and_diff_counts.csv']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% merge pathway enrichment tables
ptwNchangingGenes_KEGG_DOWN = readtable([outputFolder,...
  'ptwNchangingGenes_KEGG_DOWN_updfiltered_eggnog_edgerdeseq_FDR005_FC05.csv']);

ptwNchangingGenes_KEGG_UP = readtable([outputFolder,...
    'ptwNchangingGenes_KEGG_UP_updfiltered_eggnog_edgerdeseq_FDR005_FC05.csv']);

ptwenr_perspecies_KEGG_DOWN = readtable([outputFolder,...
 'ptwenr_recalc_KEGG_DOWN_updfiltered_eggnog_edgerdeseq_FDR005_FC05_per_species.csv']);

ptwenr_total_KEGG_DOWN = readtable([outputFolder,...
 'ptwenr_recalc_KEGG_DOWN_updfiltered_eggnog_edgerdeseq_FDR005_FC05_total.csv']);

ptwenr_perspecies_KEGG_UP = readtable([outputFolder,...
 'ptwenr_recalc_KEGG_UP_updfiltered_eggnog_edgerdeseq_FDR005_FC05_per_species.csv']);

ptwenr_total_KEGG_UP = readtable([outputFolder,...
 'ptwenr_recalc_KEGG_UP_updfiltered_eggnog_edgerdeseq_FDR005_FC05_total.csv']);

kegg_ptw_names = readtable([rawdataFolder,...
    'kegg_ptw_names_and_bact_flag.csv']);
kegg_ptw_names{:,1} = cellfun(@(x) strrep(x, 'path:', ''), kegg_ptw_names{:,1}, 'unif', 0);
kegg_ptw_names.Properties.VariableNames = {'pathwayID', 'pathwayName', 'flag_bacterial'};
kegg_ptw_names(:,3) = []; % remove bacterial flag


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% edit column names
sample_columns = cellfun(@(x) ~contains(x, 'pathwayID'), ptwNchangingGenes_KEGG_DOWN.Properties.VariableNames);
ptwNchangingGenes_KEGG_DOWN.Properties.VariableNames(sample_columns) = ...
    strcat(ptwNchangingGenes_KEGG_DOWN.Properties.VariableNames(sample_columns),...
          '_N_DOWN');

sample_columns = cellfun(@(x) ~contains(x, 'pathwayID'), ptwNchangingGenes_KEGG_UP.Properties.VariableNames);
ptwNchangingGenes_KEGG_UP.Properties.VariableNames(sample_columns) = ...
    strcat(ptwNchangingGenes_KEGG_UP.Properties.VariableNames(sample_columns),...
          '_N_UP');

sample_columns = cellfun(@(x) ~contains(x, 'pathwayID'), ptwenr_perspecies_KEGG_DOWN.Properties.VariableNames);
ptwenr_perspecies_KEGG_DOWN.Properties.VariableNames(sample_columns) = ...
    strcat(ptwenr_perspecies_KEGG_DOWN.Properties.VariableNames(sample_columns),...
          '_DOWN');

sample_columns = cellfun(@(x) ~contains(x, 'pathwayID'), ptwenr_perspecies_KEGG_UP.Properties.VariableNames);
ptwenr_perspecies_KEGG_UP.Properties.VariableNames(sample_columns) = ...
    strcat(ptwenr_perspecies_KEGG_UP.Properties.VariableNames(sample_columns),...
          '_UP');

sample_columns = cellfun(@(x) ~contains(x, 'pathwayID'), ptwenr_total_KEGG_DOWN.Properties.VariableNames);
ptwenr_total_KEGG_DOWN.Properties.VariableNames(sample_columns) = ...
    strcat(ptwenr_total_KEGG_DOWN.Properties.VariableNames(sample_columns),...
          '_DOWN');

sample_columns = cellfun(@(x) ~contains(x, 'pathwayID'), ptwenr_total_KEGG_UP.Properties.VariableNames);
ptwenr_total_KEGG_UP.Properties.VariableNames(sample_columns) = ...
    strcat(ptwenr_total_KEGG_UP.Properties.VariableNames(sample_columns),...
          '_UP');

% merge the table
merged_ptwenrTable = join(ptwNchangingGenes_KEGG_DOWN, ptwNchangingGenes_KEGG_UP, 'keys', 'pathwayID');
merged_ptwenrTable = join(merged_ptwenrTable, ptwenr_total_KEGG_DOWN, 'keys', 'pathwayID');
merged_ptwenrTable = join(merged_ptwenrTable, ptwenr_total_KEGG_UP, 'keys', 'pathwayID');
merged_ptwenrTable = join(merged_ptwenrTable, ptwenr_perspecies_KEGG_DOWN, 'keys', 'pathwayID');
merged_ptwenrTable = join(merged_ptwenrTable, ptwenr_perspecies_KEGG_UP, 'keys', 'pathwayID');
% add pathway names

merged_ptwenrTable = outerjoin(merged_ptwenrTable, kegg_ptw_names, ...
    'Type', 'Left', 'keys', 'pathwayID', 'mergekeys', true);

% Output: 
% put indexFiltered after geneFilter columns
merged_ptwenrTable = movevars(merged_ptwenrTable,'pathwayName', 'After', 'pathwayID');

% save merged table to file
writetable(merged_ptwenrTable, [outputFolder,...
    'table_S5_merged_ptwenrTable.csv']);

