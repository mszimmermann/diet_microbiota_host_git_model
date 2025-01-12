%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% match path and gene expression data
% test which of the single enzymes are expressed

% call script defining file dependenciesand global variables
addpath(genpath('.\'));
add_global_and_file_dependencies

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Required files:
% 'edgeR_gene_fold_changes_and_ann.csv'
% 'gene_annotation_EC.tsv'
% 'shortest_paths_upto_length_4.csv'
% 'countsMatrixGetMM_RNA_table.txt'
% 'countsMatrixGetMM_DNA_table.txt'
% 'metabolites_allions_combined_norm_intensity.csv'
% 'metabolites_allions_combined_formulas_with_metabolite_filters.csv'
% 'table_diff_abundance_metabolite_ions_removed2outliers.csv'
% 'table_potential_substrates_ids.csv'
% 'table_potential_products_ids.csv'
% 'merged_humann2_metaphlan_species_filtered_mouseDNA_OTUs.txt'
% Output: 
% Files:
% table_shortest_path_subprod_sorted_unique.csv
% table_speciesOTU_products_corr.csv
% table_speciesOTU_products_corrP.csv
% table_kegg_RNA_sub_prod_products_bestcorrP_pos.csv
% table_kegg_RNA_sub_prod_products_bestcorr_pos.csv
% table_kegg_RNA_sub_prod_products_bestcorr_genes_pos.csv
% table_kegg_RNA_sub_prod_products_bestcorr_EC_pos.csv
% table_kegg_DNA_sub_prod_products_bestcorr_pos.csv
% table_kegg_DNA_sub_prod_products_bestcorrP_pos.csv
% table_kegg_DNA_sub_prod_products_bestcorr_genes_pos.csv
% table_kegg_DNA_sub_prod_products_bestcorr_EC_pos.csv
% Figures:
% fig_5c_violin_best_product_corrPOS_otuFiltered_DNA_RNA_LIchange_sub_prod_signrank.pdf
% fig_5d_clustergam_best_product_corrPOS_RNA_filteredbyRNADNA_allpos.pdf
% fig_sup_clustergam_best_product_corrPOS_DNA_filteredbyRNADNA_allpos.pdf
% fig_sup_clustergam_best_product_corrPOS_OTU_filteredbyRNADNA_allpos.pdf
% fig_5e_sup_heatmap_products_filtered_byRNAandDNA_onlyPosCorr_significant0_1_sortedX.pdf
% fig_sup_scatter_metabolites_vs_OTU_DNA_per_species.pdf
% fig_sup_scatter_metabolites_vs_best_enzymes_per_species_pos_pRNA_or_pRNA_0_1.pdf
% fig_sup_scatter_metabolites_vs_best_enzymes_DNA_per_species.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% shortestPathTable = readtable([resultsFolder, ...
%     'shortest_paths_upto_length_4.csv'], 'delim', ',');
strict_class = 0;
shortestPathTable = readtable([outputFolder, ...
    'shortest_paths_combined_IP_LI_PCC_within_high_total_strictclass'...
    num2str(strict_class) '_upto_length_4.csv'],...
    'delim', ',');

% get gene FC
geneFileName = [inputFolderSeq, ...
    'edgeR_gene_fold_changes_and_ann.csv'];
% detect import options to change some of the variable to strings
opts = detectImportOptions(geneFileName);
for i = 26:42
    opts.VariableTypes{i} = 'char';  % This will change column 15 from 'double' or whatever it is to 'char', which is what you want.
end
opts.Delimiter = ',';
geneTable = readtable(geneFileName, opts);

% get EC annotation
ecTable = readtable([inputFolderSeq, ...
    'gene_annotation_EC.tsv'], 'fileType','text', 'delim', '\t');

% remove unfiltered genes
geneTable = geneTable(geneTable.geneFilter==1,:);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load gene expression from RNA and correlate with metabolite abundances
countsMatrixGetMM_RNA_table = readtable([inputFolderSeq,...
    'countsMatrixGetMM_RNA_table.txt']);
% add filter and remove unfiltered genes
countsMatrixGetMM_RNA_table = countsMatrixGetMM_RNA_table(countsMatrixGetMM_RNA_table.geneFilter==1,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load gene abundance from DNA and correlate with metabolite abundances
countsMatrixGetMM_DNA_table = readtable([inputFolderSeq,...
    'countsMatrixGetMM_DNA_table.txt']);
% add filter and remove unfiltered genes
countsMatrixGetMM_DNA_table = countsMatrixGetMM_DNA_table(countsMatrixGetMM_DNA_table.geneFilter==1,:);

% read shortbred counts for selected proteins across all mice
% countsMatrix_shortBRED = readtable(['.\ProcessedData\shortBRED\',...
%     'merged_results_4-2-1-24.txt']);
% countsMatrix_shortBRED_DNA = readtable(['.\ProcessedData\shortBRED\',...
%     'merged_metaG_shortBRED_results.txt']);
% countsMatrix_shortBRED_RNA = readtable(['.\ProcessedData\shortBRED\',...
%     'merged_metaT_shortBRED_results.txt']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get shortbred results for selected metabolites
% folderName_shortbredDNA = 'Z:\henseler\ZH002_Revisions_Diet_Microbiome_Interactions\Data\mouse_shortBRED\rerun_060524\quantify_metaG\quantified';
% countsMatrix_shortBRED_DNA = combine_shortbred_results(folderName_shortbredDNA,...
%     ['.\ProcessedData\output\',...
%     'merged_metaG_DCCVR_shortBRED_results.txt']);
% folderName_shortbredRNA = 'Z:\bartmans\Projects\ZachRerun\quantified';
% countsMatrix_shortBRED_RNA = combine_shortbred_results(folderName_shortbredRNA,...
%     ['.\ProcessedData\shortBRED\',...
%     'merged_metaT_DCCVR_shortBRED_results.txt']);
countsMatrix_shortBRED_DNA = readtable(['.\ProcessedData\output\',...
    'merged_metaG_DCCVR_shortBRED_results.txt']);
countsMatrix_shortBRED_RNA = readtable(['.\ProcessedData\shortBRED\',...
    'merged_metaT_DCCVR_shortBRED_results.txt']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read metabolite data
metaboliteData = readtable([outputFolder,...
    'metabolites_allions_combined_norm_intensity_with_CVR.csv']);

% metaboliteData = readtable([inputFolder,...
%     'metabolites_allions_combined_norm_intensity.csv']);
    
metaboliteFilters = readtable([inputFolder,...
    'metabolites_allions_combined_formulas_with_metabolite_filters.csv']);

metaboliteDiff = readtable([resultsFolder,...
    'table_diff_abundance_metabolite_ions_removed2outliers.csv']);

% potential_substrates_ids_all_table = readtable(...
%     [resultsFolder, 'table_potential_substrates_ids.csv']);
% potential_products_ids_all_table = readtable(...
%     [resultsFolder, 'table_potential_products_ids.csv']);
sel_crit1 = 'IP';
sel_crit2 = 'LI_PCC_within_high_total';
potential_substrates_ids_all_table = readtable(...
    [outputFolder, 'table_potential_substrates_ids_combined_', ...
            sel_crit1, '_', sel_crit2,'_strictclass', num2str(strict_class), '.csv']);
potential_products_ids_all_table = readtable(...
    [outputFolder, 'table_potential_products_ids_combined_', ...
            sel_crit1, '_', sel_crit2,'_strictclass', num2str(strict_class), '.csv']);


% add direction to shortestPathTable
shortestPathTable_dir = zeros(size(shortestPathTable,1),1);
for i=1:size(shortestPathTable,1)
    if ismember(shortestPathTable.Substrate(i), potential_substrates_ids_all_table.PotentialSubstrateID) &&...
       ismember(shortestPathTable.Product(i), potential_products_ids_all_table.PotentialProductID)
        shortestPathTable_dir(i)=1;
    else
        shortestPathTable_dir(i)=-1;
    end
end
shortestPathTable.Dir = shortestPathTable_dir;
clear shortestPathTable_dir

% get metabolomics data from the large intestine
met_columns = cellfun(@(x) (contains(x,'Cec') & contains(x, '_M0')),...
    metaboliteData.Properties.VariableNames);
met_columns = metaboliteData.Properties.VariableNames(met_columns);
met_columns_mice = cellfun(@(x) x(strfind(x,'M0'):end), met_columns, 'unif',0);
% get Colon and Feces columns to add to cecum
met_Colon_columns = cellfun(@(x) (contains(x,'Colon') & contains(x, '_M0')),...
    metaboliteData.Properties.VariableNames);
met_Colon_columns = metaboliteData.Properties.VariableNames(met_Colon_columns);
met_Colon_columns_mice = cellfun(@(x) x(strfind(x,'M0'):end), met_Colon_columns, 'unif',0);

met_Feces_columns = cellfun(@(x) (contains(x,'Feces') & contains(x, '_M0')),...
    metaboliteData.Properties.VariableNames);
met_Feces_columns = metaboliteData.Properties.VariableNames(met_Feces_columns);
met_Feces_columns_mice = cellfun(@(x) x(strfind(x,'M0'):end), met_Feces_columns, 'unif',0);

% intersect mice between columns for Cecum, Colon and Feces
[~,met_Colon_idx] = intersect(met_columns_mice, met_Colon_columns_mice, 'stable');
[~,met_Feces_idx] = intersect(met_columns_mice, met_Feces_columns_mice, 'stable');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get indices of KEGG ids
kegg_ids_unique = cell(size(metaboliteFilters,1),1);
idx = 1;
for i=1:size(metaboliteFilters,1)
    currkegg = metaboliteFilters.CompoundID{i};
    currkegg = strsplit(currkegg,';');
    % leave only those which have 'cpd:'
    currkegg(cellfun(@(x) ~contains(x, 'cpd:'), currkegg))=[];
    kegg_ids_unique(idx:idx+length(currkegg)-1) = currkegg;
    idx = idx+length(currkegg);
end
kegg_ids_unique(idx:end)=[];
kegg_ids_unique = unique(kegg_ids_unique);

% get idx for each kegg ID
kegg_ids_unique_index = zeros(size(kegg_ids_unique));
for i=1:size(metaboliteFilters,1)
    if (metaboliteFilters.MetaboliteFilter(i)==1)
        currkegg = metaboliteFilters.CompoundID{i};
        currkegg = strsplit(currkegg,';');
        % leave only those which have 'cpd:'
        currkegg(cellfun(@(x) ~contains(x, 'cpd:'), currkegg))=[];
        [~, curidx] = intersect(kegg_ids_unique, currkegg);
        kegg_ids_unique_index(curidx) = i;
    end
end

% remove "cpd" from kegg_ids_unique
kegg_ids_unique = cellfun(@(x) strrep(x, 'cpd:', ''), kegg_ids_unique, 'unif', 0);

% get matrices of "substrate" and "product" intensities across mice
% note that pairs are not always oriented as sub-prod
kegg_sub_intensities = zeros(size(shortestPathTable,1), length(met_columns));
kegg_prod_intensities = zeros(size(shortestPathTable,1), length(met_columns));
kegg_sub_sum_intensities = zeros(size(shortestPathTable,1), length(met_columns));
kegg_prod_sum_intensities = zeros(size(shortestPathTable,1), length(met_columns));

for i=1:size(shortestPathTable,1)%
    % sub
    cur_sub_idx = kegg_ids_unique_index(ismember(kegg_ids_unique,...
                                        shortestPathTable.Substrate{i}));
    kegg_sub_intensities(i,:) = table2array(metaboliteData(cur_sub_idx, met_columns));
    % get intensities of cecum, colon and feces and sum them up
    cur_intensities = zeros(3, length(met_columns));
    cur_intensities(1,:) = table2array(metaboliteData(cur_sub_idx, met_columns));
    cur_intensities(2,met_Colon_idx) = table2array(metaboliteData(cur_sub_idx, met_Colon_columns));
    cur_intensities(3,met_Feces_idx) = table2array(metaboliteData(cur_sub_idx, met_Feces_columns));
    % replace noise values with zeros to not sum they up
    cur_intensities(cur_intensities==intensityNoise)=0;
    % sum up values from three tissues
    cur_intensities = sum(cur_intensities,1);
    % replace zeros back to intensitynoise, if amy
    cur_intensities(cur_intensities==0)=intensityNoise;
    kegg_sub_sum_intensities(i,:) = cur_intensities;
    
    % prod
    cur_prod_idx = kegg_ids_unique_index(ismember(kegg_ids_unique,...
                                        shortestPathTable.Product{i}));
    kegg_prod_intensities(i,:) = table2array(metaboliteData(cur_prod_idx, met_columns));
    % get intensities of cecum, colon and feces and sum them up
    cur_intensities = zeros(3, length(met_columns));
    cur_intensities(1,:) = table2array(metaboliteData(cur_prod_idx, met_columns));
    cur_intensities(2,met_Colon_idx) = table2array(metaboliteData(cur_prod_idx, met_Colon_columns));
    cur_intensities(3,met_Feces_idx) = table2array(metaboliteData(cur_prod_idx, met_Feces_columns));
    % replace noise values with zeros to not sum they up
    cur_intensities(cur_intensities==intensityNoise)=0;
    % sum up values from three tissues
    cur_intensities = sum(cur_intensities,1);
    % replace zeros back to intensitynoise, if amy
    cur_intensities(cur_intensities==0)=intensityNoise;
    kegg_prod_sum_intensities(i,:) = cur_intensities;

    if mod(i,1000)==0
        fprintf('Calculated intensities for %d metabolites\n',i);
    end
end
% get names of metabolomics samples in the same format as metagenomics
metabolomics_samples = cellfun(@(x) x(strfind(x,'_M')+1:end), met_columns, 'unif', 0);
metabolomics_samples = cellfun(@(x) strrep(x, 'M0', 'MSZ'), metabolomics_samples, 'unif', 0);
[~, metabolomics_idx, gene_idx] = intersect(metabolomics_samples, countsMatrixGetMM_RNA_table.Properties.VariableNames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get information from differential analysis on whether substrates and
% products are changing in at least one low intestine tissue in at least
% one diet between GF and colonized mice
liTissues = {'Cecum', 'Colon', 'Feces'};
dietChangeNames = {'HFD_DCGF', 'CTR_DCGF'};
fcThreshold = log2(1.5);
pfdrThreshold = 0.05;
pvalColumn = '_pFDR_'; % choose between _p_ and _pFDR_
kegg_sub_li_changes = zeros(size(shortestPathTable,1), length(liTissues)*length(dietChangeNames));
kegg_prod_li_changes = zeros(size(shortestPathTable,1), length(liTissues)*length(dietChangeNames));
% record also liver and serum changes
kegg_sub_liver_changes = zeros(size(shortestPathTable,1), length(dietChangeNames));
kegg_prod_liver_changes = zeros(size(shortestPathTable,1), length(dietChangeNames));
kegg_sub_serum_changes = zeros(size(shortestPathTable,1), length(dietChangeNames));
kegg_prod_serum_changes = zeros(size(shortestPathTable,1), length(dietChangeNames));
for i=1:size(shortestPathTable,1)
    idx=1;
    % sub
    cur_sub_idx = kegg_ids_unique_index(ismember(kegg_ids_unique,...
                                        shortestPathTable.Substrate{i}));
    % prod
    cur_prod_idx = kegg_ids_unique_index(ismember(kegg_ids_unique,...
                                        shortestPathTable.Product{i}));
    for diet_i=1:length(dietChangeNames)
        for tissue_i=1:length(liTissues)
            % sub
            kegg_sub_li_changes(i,idx) = (abs(log2(metaboliteDiff{cur_sub_idx,...
                 ismember(metaboliteDiff.Properties.VariableNames,...
                 {[liTissues{tissue_i} '_fc_', dietChangeNames{diet_i}]})})) >= fcThreshold).*...
                 (metaboliteDiff{cur_sub_idx,...
                 ismember(metaboliteDiff.Properties.VariableNames,...
                 {[liTissues{tissue_i} pvalColumn, dietChangeNames{diet_i}]})} <= pfdrThreshold);
      
            % prod
            kegg_prod_li_changes(i,idx) = (abs(log2(metaboliteDiff{cur_prod_idx,...
                 ismember(metaboliteDiff.Properties.VariableNames,...
                 {[liTissues{tissue_i} '_fc_', dietChangeNames{diet_i}]})})) >= fcThreshold).*...
                 (metaboliteDiff{cur_prod_idx,...
                 ismember(metaboliteDiff.Properties.VariableNames,...
                 {[liTissues{tissue_i} pvalColumn, dietChangeNames{diet_i}]})} <= pfdrThreshold);
            idx = idx+1;
        end
        % liver
        % sub (record sign of change)
        kegg_sub_liver_changes(i,diet_i) = (abs(log2(metaboliteDiff{cur_sub_idx,...
             ismember(metaboliteDiff.Properties.VariableNames,...
             {['Liver_fc_', dietChangeNames{diet_i}]})})) >= fcThreshold).*...
             (metaboliteDiff{cur_sub_idx,...
             ismember(metaboliteDiff.Properties.VariableNames,...
             {['Liver' pvalColumn, dietChangeNames{diet_i}]})} <= pfdrThreshold).*...
             sign(log2(metaboliteDiff{cur_sub_idx,...
             ismember(metaboliteDiff.Properties.VariableNames,...
             {['Liver_fc_', dietChangeNames{diet_i}]})}));

        % prod (record sign of change)
        kegg_prod_li_changes(i,diet_i) = (abs(log2(metaboliteDiff{cur_prod_idx,...
             ismember(metaboliteDiff.Properties.VariableNames,...
             {['Liver_fc_', dietChangeNames{diet_i}]})})) >= fcThreshold).*...
             (metaboliteDiff{cur_prod_idx,...
             ismember(metaboliteDiff.Properties.VariableNames,...
             {['Liver' pvalColumn, dietChangeNames{diet_i}]})} <= pfdrThreshold).*...
             sign(log2(metaboliteDiff{cur_prod_idx,...
             ismember(metaboliteDiff.Properties.VariableNames,...
             {['Liver_fc_', dietChangeNames{diet_i}]})}));
        % serum
        % sub (record sign of change)
        kegg_sub_serum_changes(i,diet_i) = (abs(log2(metaboliteDiff{cur_sub_idx,...
             ismember(metaboliteDiff.Properties.VariableNames,...
             {['Serum_fc_', dietChangeNames{diet_i}]})})) >= fcThreshold).*...
             (metaboliteDiff{cur_sub_idx,...
             ismember(metaboliteDiff.Properties.VariableNames,...
             {['Serum' pvalColumn, dietChangeNames{diet_i}]})} <= pfdrThreshold).*...
             sign(log2(metaboliteDiff{cur_sub_idx,...
             ismember(metaboliteDiff.Properties.VariableNames,...
             {['Serum_fc_', dietChangeNames{diet_i}]})}));

        % prod (record sign of change)
        kegg_prod_serum_changes(i,diet_i) = (abs(log2(metaboliteDiff{cur_prod_idx,...
             ismember(metaboliteDiff.Properties.VariableNames,...
             {['Serum_fc_', dietChangeNames{diet_i}]})})) >= fcThreshold).*...
             (metaboliteDiff{cur_prod_idx,...
             ismember(metaboliteDiff.Properties.VariableNames,...
             {['Serum' pvalColumn, dietChangeNames{diet_i}]})} <= pfdrThreshold).*...
             sign(log2(metaboliteDiff{cur_prod_idx,...
             ismember(metaboliteDiff.Properties.VariableNames,...
             {['Serum_fc_', dietChangeNames{diet_i}]})}));
             
    end
    
    if mod(i,1000)==0
        fprintf('Extracted fold changes for %d metabolites\n',i);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make a table of max abs enzyme fold change between the two conditions
%substrates
shortestPathTable_length = zeros(size(shortestPathTable,1),1);
for i=1:size(shortestPathTable,1)
    shortestPathTable_length(i) = length(strsplit(shortestPathTable.EC_path{i},';'));
end
% plot histogram of length distributions
figure
histogram(shortestPathTable_length)
set(gca, 'XTick', 1:max(shortestPathTable_length))
axis square
for i=1:max(shortestPathTable_length)
    text(i-0.25, nnz(shortestPathTable_length==i)+250, ...
        sprintf('N=%d', nnz(shortestPathTable_length==i)));
end
ylim([0 15000])
title('Number of substrate-product paths')
xlabel('Path length')
ylabel('Number of paths')
%print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
%    [figureFolder, 'fig_sup_histogram_path_length_max3.pdf'])
print(gcf, '-vector', '-dpdf', '-r600', '-bestfit', ...
    [figureFolder, 'fig_sup_histogram_path_length_max3_combined_', ...
            sel_crit1, '_', sel_crit2,'_strictclass', num2str(strict_class), '.pdf']);

% calculate for 1 enzyme path to begin with
select_path_idx = find(shortestPathTable_length==1);
kegg_substrate_product_1enzyme = shortestPathTable(select_path_idx,:);
[kegg_substrate_product_1enzyme, idx] = unique(kegg_substrate_product_1enzyme);
select_path_idx = select_path_idx(idx);
species_list = unique(geneTable.abbrSpecies);
geneTable_EC = geneTable.EC;
geneTable_EC_unique = cell(length(geneTable_EC),1);
idx = 1;
for i=1:length(geneTable_EC)
    curec = geneTable_EC{i};
    curec = strsplit(curec,',');
    %keep only EC that have three dots
    curec(cellfun(@(x) length(strfind(x,'.'))~=3, curec))=[];
    geneTable_EC_unique(idx:idx+length(curec)-1) = curec;
    idx = idx+length(curec);
end
geneTable_EC_unique(idx:end)=[];
geneTable_EC_unique = unique(geneTable_EC_unique);

geneTable_EC_species = cell(length(geneTable_EC_unique),...
                            length(species_list));
for i=1:length(geneTable_EC)
    curec = geneTable_EC{i};
    curec = strsplit(curec,',');
    %keep only EC that have three dots
    curec(cellfun(@(x) length(strfind(x,'.'))~=3, curec))=[];
    for j=1:length(curec)
        [~, ecidx] = intersect(geneTable_EC_unique, curec{j});
        [~, species_idx] = intersect(species_list, geneTable.abbrSpecies{i});
        geneTable_EC_species{ecidx,species_idx} = ...
            [geneTable_EC_species{ecidx,species_idx} i];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate correlations between RNA and metabolites
kegg_sub_prod_EC_sub_corr = zeros(size(kegg_substrate_product_1enzyme,1), length(species_list));
kegg_sub_prod_EC_sub_corrP = ones(size(kegg_substrate_product_1enzyme,1), length(species_list));
kegg_sub_prod_EC_sub_sum_corr = zeros(size(kegg_substrate_product_1enzyme,1), length(species_list));
kegg_sub_prod_EC_sub_sum_corrP = ones(size(kegg_substrate_product_1enzyme,1), length(species_list));
kegg_sub_prod_EC_sub_corr_id = cell(size(kegg_substrate_product_1enzyme,1), length(species_list));
kegg_sub_prod_EC_sub_corr_geneidx = cell(size(kegg_substrate_product_1enzyme,1), length(species_list));
% prepare matrices for linear model fits and r adjusted
kegg_sub_prod_EC_sub_mdlX = zeros(size(kegg_substrate_product_1enzyme,1), length(species_list));
kegg_sub_prod_EC_sub_mdlP = ones(size(kegg_substrate_product_1enzyme,1), length(species_list));
kegg_sub_prod_EC_sub_mdlRsqadj = zeros(size(kegg_substrate_product_1enzyme,1), length(species_list));

% substrates
for i=1:size(kegg_substrate_product_1enzyme,1)
    curenzymes = strsplit(kegg_substrate_product_1enzyme.EC_path{i},'|');
    if ~isempty(curenzymes)
        % intersect enzymes with the ones from the gene table
        [~, idxEC, idxGene] = intersect(curenzymes, geneTable_EC_unique, 'stable');
        if ~isempty(idxEC)
            curenzymes_geneidx = geneTable_EC_species(idxGene,:);
            % prepare tables for correlation values etc
            curenzymes_mat = zeros(length(idxEC), length(species_list));
            curenzymes_p = ones(length(idxEC), length(species_list));
            curenzymes_ec = cell(length(idxEC), length(species_list));
            % for linear model
            curenzymes_mdlX = zeros(length(idxEC), length(species_list));
            curenzymes_mdlXp = ones(length(idxEC), length(species_list));
            curenzymes_mdlRsqadj = zeros(length(idxEC), length(species_list));
            
            for k=1:size(curenzymes_geneidx,2)
                for j=1:size(curenzymes_geneidx,1)
                    if ~isempty(curenzymes_geneidx{j,k})
                        cur_expression = table2array(countsMatrixGetMM_RNA_table(curenzymes_geneidx{j,k},gene_idx))';
                        cur_expression = sum(cur_expression, 2);

                        if sum(cur_expression)>0
                            % corr with substrate
                            y = kegg_sub_sum_intensities(select_path_idx(i),metabolomics_idx)';                          

                            [curcorr, curcorrP] = corr(cur_expression,y);
                            % add to matrix
                            curenzymes_mat(j,k) = curcorr;
                            curenzymes_p(j,k) = curcorrP;
                            
                             % calculate mdl coefficients and rsquared
                            mdl = fitlm(cur_expression,y);
                            curenzymes_mdlX(j,k) = mdl.Coefficients.Estimate(2);
                            curenzymes_mdlXp(j,k) = mdl.Coefficients.pValue(2);
                            curenzymes_mdlRsqadj(j,k) = mdl.Rsquared.Adjusted;
                        end
                        curenzymes_ec{j, k} = curec;
                        
                    end
                end
                % calculate corr with sum of all enzymes
                if ~isempty(horzcat(curenzymes_geneidx{:,k}))
                    cur_expression = table2array(countsMatrixGetMM_RNA_table(horzcat(curenzymes_geneidx{:,k}),gene_idx))';
                    cur_expression = sum(cur_expression, 2);
                    
                    if sum(cur_expression)>0
                        % corr with substrate
                        y = kegg_sub_sum_intensities(select_path_idx(i),metabolomics_idx)';                          

                        [curcorr, curcorrP] = corr(cur_expression, y);
                        % add to matrix
                        kegg_sub_prod_EC_sub_sum_corr(i,k) = curcorr;
                        kegg_sub_prod_EC_sub_sum_corrP(i,k) = curcorrP;
                        
                    end
                end
            end
            % keep max absolute value per species across all relevant enzymes
            [~, idx] = max(abs(curenzymes_mat), [], 1, 'linear');
            kegg_sub_prod_EC_sub_corr(i,:) = curenzymes_mat(idx);
            kegg_sub_prod_EC_sub_corrP(i,:) = curenzymes_p(idx);
            kegg_sub_prod_EC_sub_corr_id(i,:) = curenzymes_ec(idx);
            kegg_sub_prod_EC_sub_corr_geneidx(i,:) = curenzymes_geneidx(idx);
            kegg_sub_prod_EC_sub_mdlX(i,:)  = curenzymes_mdlX(idx);
            kegg_sub_prod_EC_sub_mdlP(i,:)  = curenzymes_mdlXp(idx);
            kegg_sub_prod_EC_sub_mdlRsqadj(i,:) = curenzymes_mdlRsqadj(idx);
        end
    end
end
% products
kegg_sub_prod_EC_prod_corr = zeros(size(kegg_substrate_product_1enzyme,1), length(species_list));
kegg_sub_prod_EC_prod_corrP = ones(size(kegg_substrate_product_1enzyme,1), length(species_list));
kegg_sub_prod_EC_prod_sum_corr = zeros(size(kegg_substrate_product_1enzyme,1), length(species_list));
kegg_sub_prod_EC_prod_sum_corrP = ones(size(kegg_substrate_product_1enzyme,1), length(species_list));
kegg_sub_prod_EC_prod_corr_id = cell(size(kegg_substrate_product_1enzyme,1), length(species_list));
kegg_sub_prod_EC_prod_corr_geneidx = cell(size(kegg_substrate_product_1enzyme,1), length(species_list));
% prepare matrices for linear model fits and r adjusted
kegg_sub_prod_EC_prod_mdlX = zeros(size(kegg_substrate_product_1enzyme,1), length(species_list));
kegg_sub_prod_EC_prod_mdlP = ones(size(kegg_substrate_product_1enzyme,1), length(species_list));
kegg_sub_prod_EC_prod_mdlRsqadj = zeros(size(kegg_substrate_product_1enzyme,1), length(species_list));

for i=1:size(kegg_substrate_product_1enzyme,1)
    curenzymes = strsplit(kegg_substrate_product_1enzyme.EC_path{i},'|');
    if ~isempty(curenzymes)
        % intersect enzymes with the ones from the gene table
        [~, idxEC, idxGene] = intersect(curenzymes, geneTable_EC_unique, 'stable');
        if ~isempty(idxEC)
            curenzymes_geneidx = geneTable_EC_species(idxGene,:);
            % prepare tables for correlation values etc
            curenzymes_mat = zeros(length(idxEC), length(species_list));
            curenzymes_p = ones(length(idxEC), length(species_list));
            curenzymes_ec = cell(length(idxEC), length(species_list));
            % for linear model
            curenzymes_mdlX = zeros(length(idxEC), length(species_list));
            curenzymes_mdlXp = ones(length(idxEC), length(species_list));
            curenzymes_mdlRsqadj = zeros(length(idxEC), length(species_list));

            for k=1:size(curenzymes_geneidx,2)
                for j=1:size(curenzymes_geneidx,1)
                    if ~isempty(curenzymes_geneidx{j,k})
                        cur_expression = table2array(countsMatrixGetMM_RNA_table(curenzymes_geneidx{j,k},gene_idx))';
                        cur_expression = sum(cur_expression, 2);

                        if sum(cur_expression)>0
                            % corr with substrate
                            y = kegg_prod_sum_intensities(select_path_idx(i),metabolomics_idx)';
                            [curcorr, curcorrP] = corr(cur_expression,y);
                            % add to matrix
                            curenzymes_mat(j,k) = curcorr;
                            curenzymes_p(j,k) = curcorrP;
                            
                             % calculate mdl coefficients and rsquared
                            mdl = fitlm(cur_expression,y);
                            curenzymes_mdlX(j,k) = mdl.Coefficients.Estimate(2);
                            curenzymes_mdlXp(j,k) = mdl.Coefficients.pValue(2);
                            curenzymes_mdlRsqadj(j,k) = mdl.Rsquared.Adjusted;
                        end
                        curenzymes_ec{j, k} = curec;
                    end
                end
                % calculate corr with sum of all enzymes
                if ~isempty(horzcat(curenzymes_geneidx{:,k}))
                    cur_expression = table2array(countsMatrixGetMM_RNA_table(horzcat(curenzymes_geneidx{:,k}),gene_idx))';
                    cur_expression = sum(cur_expression, 2);

                    if sum(cur_expression)>0
                        % corr with substrate
                        [curcorr, curcorrP] = corr(cur_expression,...
                            kegg_prod_sum_intensities(select_path_idx(i),metabolomics_idx)');
                        % add to matrix
                        kegg_sub_prod_EC_prod_sum_corr(i,k) = curcorr;
                        kegg_sub_prod_EC_prod_sum_corrP(i,k) = curcorrP;
                    end
                end
            end
            % keep max absolute value per species across all relevant enzymes
            [~, idx] = max(abs(curenzymes_mat), [], 1, 'linear');
            kegg_sub_prod_EC_prod_corr(i,:) = curenzymes_mat(idx);
            kegg_sub_prod_EC_prod_corrP(i,:) = curenzymes_p(idx);
            kegg_sub_prod_EC_prod_corr_id(i,:) = curenzymes_ec(idx);
            kegg_sub_prod_EC_prod_corr_geneidx(i,:) = curenzymes_geneidx(idx);
            kegg_sub_prod_EC_prod_mdlX(i,:) = curenzymes_mdlX(idx);
            kegg_sub_prod_EC_prod_mdlP(i,:) = curenzymes_mdlXp(idx);
            kegg_sub_prod_EC_prod_mdlRsqadj(i,:) = curenzymes_mdlRsqadj(idx);
        end
    end
end

% get sub and prod LI changes in 1 enzyme paths
kegg_sub_prod_sub_LIchange = kegg_sub_li_changes(select_path_idx,:);
kegg_sub_prod_prod_LIchange = kegg_prod_li_changes(select_path_idx,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% aggregate correlation between enzymes and metabolites per metabolite, 
% keep only the strongest correlation per species
kegg_sub_prod_substrates = cell(size(kegg_substrate_product_1enzyme,1),1);
kegg_sub_prod_products = cell(size(kegg_substrate_product_1enzyme,1),1);
for i=1:size(kegg_substrate_product_1enzyme,1)
    if kegg_substrate_product_1enzyme.Dir(i)==1
        kegg_sub_prod_substrates{i} = kegg_substrate_product_1enzyme.Substrate{i};
        kegg_sub_prod_products{i} = kegg_substrate_product_1enzyme.Product{i};
    else
        kegg_sub_prod_substrates{i} = kegg_substrate_product_1enzyme.Product{i};
        kegg_sub_prod_products{i} = kegg_substrate_product_1enzyme.Substrate{i};
    end
end
% keep only unique substrates and products
kegg_sub_prod_substrates_unique = unique(kegg_sub_prod_substrates);
kegg_sub_prod_products_unique = unique(kegg_sub_prod_products);

% calculate best correlation per species across all enzymes
kegg_sub_prod_substrates_bestcorr_pos = zeros(length(kegg_sub_prod_substrates_unique),...
    size(kegg_sub_prod_EC_prod_corr,2));
kegg_sub_prod_substrates_bestcorr_neg = zeros(length(kegg_sub_prod_substrates_unique),...
    size(kegg_sub_prod_EC_prod_corr,2));
kegg_sub_prod_substrates_bestcorrP_pos = zeros(length(kegg_sub_prod_substrates_unique),...
    size(kegg_sub_prod_EC_prod_corr,2));
kegg_sub_prod_substrates_bestcorrP_neg = zeros(length(kegg_sub_prod_substrates_unique),...
    size(kegg_sub_prod_EC_prod_corr,2));
kegg_sub_prod_substrates_bestcorr_geneidx_pos = cell(length(kegg_sub_prod_substrates_unique),...
    size(kegg_sub_prod_EC_prod_corr,2));
kegg_sub_prod_substrates_bestcorr_geneidx_neg = cell(length(kegg_sub_prod_substrates_unique),...
    size(kegg_sub_prod_EC_prod_corr,2));
for i=1:length(kegg_sub_prod_substrates_unique)
    curidx = ismember(kegg_sub_prod_substrates, kegg_sub_prod_substrates_unique{i});
    curcorr = (kegg_sub_prod_EC_sub_corr.*(kegg_substrate_product_1enzyme.Dir==1)) +...
           (kegg_sub_prod_EC_prod_corr.*(kegg_substrate_product_1enzyme.Dir==-1));
    curcorr=curcorr(curidx,:);
    curcorrP = (kegg_sub_prod_EC_sub_corrP.*(kegg_substrate_product_1enzyme.Dir==1)) +...
           (kegg_sub_prod_EC_prod_corrP.*(kegg_substrate_product_1enzyme.Dir==-1));
    curcorrP=curcorrP(curidx,:);
    curgeneidx = kegg_sub_prod_EC_sub_corr_geneidx;
    
    % get info about change in LI for substrate and product
    cursubLI = sum(kegg_sub_prod_sub_LIchange(curidx,:),2)>0;
    curprodLI = sum(kegg_sub_prod_prod_LIchange(curidx,:),2)>0;
    
    %multiply corr by changes of substrate and products
    curcorr = curcorr.*repmat(cursubLI,1,size(curcorr,2)).*repmat(curprodLI,1,size(curcorr,2));
    
    for j=1:length(kegg_substrate_product_1enzyme.Dir)
        if kegg_substrate_product_1enzyme.Dir(j)==-1
            curgeneidx(j,:) = kegg_sub_prod_EC_prod_corr_geneidx(j,:);
        end
    end
    curgeneidx=curgeneidx(curidx,:);
    
    [kegg_sub_prod_substrates_bestcorr_pos(i,:),idx] = max(curcorr,[],1);
    for j=1:length(idx)
        kegg_sub_prod_substrates_bestcorrP_pos(i,j) = curcorrP(idx(j),j);
        kegg_sub_prod_substrates_bestcorr_geneidx_pos{i,j} = curgeneidx{idx(j),j};
    end
    [kegg_sub_prod_substrates_bestcorr_neg(i,:),idx] = min(curcorr,[],1);
    for j=1:length(idx)
        kegg_sub_prod_substrates_bestcorrP_neg(i,j) = curcorrP(idx(j),j);
        kegg_sub_prod_substrates_bestcorr_geneidx_neg{i,j} = curgeneidx{idx(j),j};
    end
end
% products 
kegg_sub_prod_products_bestcorr_pos = zeros(length(kegg_sub_prod_products_unique),...
    size(kegg_sub_prod_EC_prod_corr,2));
kegg_sub_prod_products_bestcorr_neg = zeros(length(kegg_sub_prod_products_unique),...
    size(kegg_sub_prod_EC_prod_corr,2));
kegg_sub_prod_products_bestcorrP_pos = zeros(length(kegg_sub_prod_products_unique),...
    size(kegg_sub_prod_EC_prod_corr,2));
kegg_sub_prod_products_bestcorrP_neg = zeros(length(kegg_sub_prod_products_unique),...
    size(kegg_sub_prod_EC_prod_corr,2));
kegg_sub_prod_products_bestcorr_geneidx_pos = cell(length(kegg_sub_prod_products_unique),...
    size(kegg_sub_prod_EC_prod_corr,2));
kegg_sub_prod_products_bestcorr_geneidx_neg = cell(length(kegg_sub_prod_products_unique),...
    size(kegg_sub_prod_EC_prod_corr,2));
for i=1:length(kegg_sub_prod_products_unique)
    curidx = ismember(kegg_sub_prod_products, kegg_sub_prod_products_unique{i});
    curcorr = (kegg_sub_prod_EC_prod_corr.*(kegg_substrate_product_1enzyme.Dir==1)) +...
           (kegg_sub_prod_EC_sub_corr.*(kegg_substrate_product_1enzyme.Dir==-1));
    curcorr=curcorr(curidx,:);
    curcorrP = (kegg_sub_prod_EC_prod_corrP.*(kegg_substrate_product_1enzyme.Dir==1)) +...
           (kegg_sub_prod_EC_sub_corrP.*(kegg_substrate_product_1enzyme.Dir==-1));
    curcorrP=curcorrP(curidx,:);
    
    % get info about change in LI for substrate and product
    cursubLI = sum(kegg_sub_prod_sub_LIchange(curidx,:),2)>0;
    curprodLI = sum(kegg_sub_prod_prod_LIchange(curidx,:),2)>0;
    
    %multiply corr by changes of substrate and products
    curcorr = curcorr.*repmat(cursubLI,1,size(curcorr,2)).*repmat(curprodLI,1,size(curcorr,2));
 
    
    curgeneidx = kegg_sub_prod_EC_prod_corr_geneidx;
    for j=1:length(kegg_substrate_product_1enzyme.Dir)
        if kegg_substrate_product_1enzyme.Dir(j)==-1
            curgeneidx(j,:) = kegg_sub_prod_EC_sub_corr_geneidx(j,:);
        end
    end
    curgeneidx=curgeneidx(curidx,:);
    
    [kegg_sub_prod_products_bestcorr_pos(i,:), idx] = max(curcorr,[],1);    
    for j=1:length(idx)
        kegg_sub_prod_products_bestcorrP_pos(i,j) = curcorrP(idx(j),j);
        kegg_sub_prod_products_bestcorr_geneidx_pos{i,j} = curgeneidx{idx(j),j};
    end
    [kegg_sub_prod_products_bestcorr_neg(i,:), idx] = min(curcorr,[],1);
    for j=1:length(idx)
        kegg_sub_prod_products_bestcorrP_neg(i,j) = curcorrP(idx(j),j);
        kegg_sub_prod_products_bestcorr_geneidx_neg{i,j} = curgeneidx{idx(j),j};
    end
end

 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % plot correlation between species enzymes and metabolites 
% % given metabolite and enzyme ID
% %figureFile = 'scatter_metabolites_vs_best_enzymes_per_species.ps';
% plot_met_id = 'C00106';% 'C00253';
% cur_met_idx = 45782
% %cur_met_idx = kegg_ids_unique_index(ismember(kegg_ids_unique,...
% %                                        plot_met_id));
% y = table2array(metaboliteData(cur_met_idx, met_columns));
% y = y(metabolomics_idx)';  
% plot_ec_id = '2.4.2.3';%'2.4.2.2';%'6.3.4.21';%'3.5.1.19';%'2.4.2.1';%
% cur_ec_idx = ismember(ecTable.PathwayID, plot_ec_id);
% cur_data_genes = ecTable.GeneIDXfiltered(cur_ec_idx);
% cur_data_genes = strsplit(cur_data_genes{1},';');
% cur_data_genes(end) = [];
% cur_data_genes = cellfun(@(x) str2num(x), cur_data_genes);
% 
% fig = figure('units','normalized','outerposition',[0 0 1 1]);
% spidx = 1;
% for k=1:length(cur_data_genes)
%     cur_expression = table2array(countsMatrixGetMM_RNA_table(...
%                 cur_data_genes(k),gene_idx))';
%     cur_expression = sum(cur_expression, 2);
%     %cur_gene_id = countsMatrixGetMM_RNA_table.species_genes{cur_data_genes(k)};
%     cur_gene_id = geneTable.gene_name{cur_data_genes(k)};
%     
%     if sum(cur_expression)>0
%         % corr with substrate
%         X = cur_expression;
%            % remove min and max values
%             %X(y==max(y))=[];
%             %y(y==max(y))=[];
%             %X(y==min(y))=[];
%             %y(y==min(y))=[];
%             [curcorr, curcorrP] = corr(X,y);
%             [curcorrS] = corr(X,y, 'type','spearman');                          
% 
% 
%             mdl = fitlm(X,y);
%             subplot(3,5,spidx);
%             scatter(X,y);
%             axis square
%             title(sprintf('%s\n%.2f %.2f %.2f %.2f',cur_gene_id,...
%                                          curcorr, curcorrP, curcorrS,...
%                                          mdl.Rsquared.Adjusted))
%    end
%    spidx = spidx+1;
% end
% suptitle(sprintf('%s %s (corr, corrP, corrS, Rsqadj)',plot_met_id, plot_ec_id));
% orient landscape 
%    print(gcf, '-painters', '-dpsc2', '-r600', '-append', '-bestfit',...
%        [figureFolder, figureFile]);
%    close(fig)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate correlations between DNA and metabolites
kegg_DNA_sub_prod_EC_sub_corr = zeros(size(kegg_substrate_product_1enzyme,1), length(species_list));
kegg_DNA_sub_prod_EC_sub_corrP = ones(size(kegg_substrate_product_1enzyme,1), length(species_list));
kegg_DNA_sub_prod_EC_sub_sum_corr = zeros(size(kegg_substrate_product_1enzyme,1), length(species_list));
kegg_DNA_sub_prod_EC_sub_sum_corrP = ones(size(kegg_substrate_product_1enzyme,1), length(species_list));
kegg_DNA_sub_prod_EC_sub_corr_id = cell(size(kegg_substrate_product_1enzyme,1), length(species_list));
kegg_DNA_sub_prod_EC_sub_corr_geneidx = cell(size(kegg_substrate_product_1enzyme,1), length(species_list));
% prepare matrices for linear model fits and r adjusted
kegg_DNA_sub_prod_EC_sub_mdlX = zeros(size(kegg_substrate_product_1enzyme,1), length(species_list));
kegg_DNA_sub_prod_EC_sub_mdlP = ones(size(kegg_substrate_product_1enzyme,1), length(species_list));
kegg_DNA_sub_prod_EC_sub_mdlRsqadj = zeros(size(kegg_substrate_product_1enzyme,1), length(species_list));

% substrates
for i=1:size(kegg_substrate_product_1enzyme,1)
    curenzymes = strsplit(kegg_substrate_product_1enzyme.EC_path{i},'|');
    if ~isempty(curenzymes)
        % intersect enzymes with the ones from the gene table
        [~, idxEC, idxGene] = intersect(curenzymes, geneTable_EC_unique, 'stable');
        if ~isempty(idxEC)
            curenzymes_geneidx = geneTable_EC_species(idxGene,:);
            % prepare tables for correlation values etc
            curenzymes_mat = zeros(length(idxEC), length(species_list));
            curenzymes_p = ones(length(idxEC), length(species_list));
            curenzymes_ec = cell(length(idxEC), length(species_list));
            % for linear model
            curenzymes_mdlX = zeros(length(idxEC), length(species_list));
            curenzymes_mdlXp = ones(length(idxEC), length(species_list));
            curenzymes_mdlRsqadj = zeros(length(idxEC), length(species_list));
            
            for k=1:size(curenzymes_geneidx,2)
                for j=1:size(curenzymes_geneidx,1)
                    if ~isempty(curenzymes_geneidx{j,k})
                        cur_expression = table2array(countsMatrixGetMM_DNA_table(curenzymes_geneidx{j,k},gene_idx))';
                        cur_expression = sum(cur_expression, 2);

                        if sum(cur_expression)>0
                            % corr with substrate
                            y = kegg_sub_sum_intensities(select_path_idx(i),metabolomics_idx)';                          

                            [curcorr, curcorrP] = corr(cur_expression,y);
                            % add to matrix
                            curenzymes_mat(j,k) = curcorr;
                            curenzymes_p(j,k) = curcorrP;
                            
                             % calculate mdl coefficients and rsquared
                            mdl = fitlm(cur_expression,y);
                            curenzymes_mdlX(j,k) = mdl.Coefficients.Estimate(2);
                            curenzymes_mdlXp(j,k) = mdl.Coefficients.pValue(2);
                            curenzymes_mdlRsqadj(j,k) = mdl.Rsquared.Adjusted;
                        end
                        curenzymes_ec{j, k} = curec;
                        
                    end
                end
                % calculate corr with sum of all enzymes
                if ~isempty(horzcat(curenzymes_geneidx{:,k}))
                    cur_expression = table2array(countsMatrixGetMM_DNA_table(horzcat(curenzymes_geneidx{:,k}),gene_idx))';
                    cur_expression = sum(cur_expression, 2);
                    
                    if sum(cur_expression)>0
                        % corr with substrate
                        y = kegg_sub_sum_intensities(select_path_idx(i),metabolomics_idx)';                          

                        [curcorr, curcorrP] = corr(cur_expression, y);
                        % add to matrix
                        kegg_DNA_sub_prod_EC_sub_sum_corr(i,k) = curcorr;
                        kegg_DNA_sub_prod_EC_sub_sum_corrP(i,k) = curcorrP;
                        
                    end
                end
            end
            % keep max absolute value per species across all relevant enzymes
            [~, idx] = max(abs(curenzymes_mat), [], 1, 'linear');
            kegg_DNA_sub_prod_EC_sub_corr(i,:) = curenzymes_mat(idx);
            kegg_DNA_sub_prod_EC_sub_corrP(i,:) = curenzymes_p(idx);
            kegg_DNA_sub_prod_EC_sub_corr_id(i,:) = curenzymes_ec(idx);
            kegg_DNA_sub_prod_EC_sub_corr_geneidx(i,:) = curenzymes_geneidx(idx);
            kegg_DNA_sub_prod_EC_sub_mdlX(i,:)  = curenzymes_mdlX(idx);
            kegg_DNA_sub_prod_EC_sub_mdlP(i,:)  = curenzymes_mdlXp(idx);
            kegg_DNA_sub_prod_EC_sub_mdlRsqadj(i,:) = curenzymes_mdlRsqadj(idx);
        end
    end
end
% products
kegg_DNA_sub_prod_EC_prod_corr = zeros(size(kegg_substrate_product_1enzyme,1), length(species_list));
kegg_DNA_sub_prod_EC_prod_corrP = ones(size(kegg_substrate_product_1enzyme,1), length(species_list));
kegg_DNA_sub_prod_EC_prod_sum_corr = zeros(size(kegg_substrate_product_1enzyme,1), length(species_list));
kegg_DNA_sub_prod_EC_prod_sum_corrP = ones(size(kegg_substrate_product_1enzyme,1), length(species_list));
kegg_DNA_sub_prod_EC_prod_corr_id = cell(size(kegg_substrate_product_1enzyme,1), length(species_list));
kegg_DNA_sub_prod_EC_prod_corr_geneidx = cell(size(kegg_substrate_product_1enzyme,1), length(species_list));
% prepare matrices for linear model fits and r adjusted
kegg_DNA_sub_prod_EC_prod_mdlX = zeros(size(kegg_substrate_product_1enzyme,1), length(species_list));
kegg_DNA_sub_prod_EC_prod_mdlP = ones(size(kegg_substrate_product_1enzyme,1), length(species_list));
kegg_DNA_sub_prod_EC_prod_mdlRsqadj = zeros(size(kegg_substrate_product_1enzyme,1), length(species_list));

for i=1:size(kegg_substrate_product_1enzyme,1)
    curenzymes = strsplit(kegg_substrate_product_1enzyme.EC_path{i},'|');
    if ~isempty(curenzymes)
        % intersect enzymes with the ones from the gene table
        [~, idxEC, idxGene] = intersect(curenzymes, geneTable_EC_unique, 'stable');
        if ~isempty(idxEC)
            curenzymes_geneidx = geneTable_EC_species(idxGene,:);
            % prepare tables for correlation values etc
            curenzymes_mat = zeros(length(idxEC), length(species_list));
            curenzymes_p = ones(length(idxEC), length(species_list));
            curenzymes_ec = cell(length(idxEC), length(species_list));
            % for linear model
            curenzymes_mdlX = zeros(length(idxEC), length(species_list));
            curenzymes_mdlXp = ones(length(idxEC), length(species_list));
            curenzymes_mdlRsqadj = zeros(length(idxEC), length(species_list));

            for k=1:size(curenzymes_geneidx,2)
                for j=1:size(curenzymes_geneidx,1)
                    if ~isempty(curenzymes_geneidx{j,k})
                        cur_expression = table2array(countsMatrixGetMM_DNA_table(curenzymes_geneidx{j,k},gene_idx))';
                        cur_expression = sum(cur_expression, 2);

                        if sum(cur_expression)>0
                            % corr with substrate
                            y = kegg_prod_sum_intensities(select_path_idx(i),metabolomics_idx)';
                            [curcorr, curcorrP] = corr(cur_expression,y);
                            % add to matrix
                            curenzymes_mat(j,k) = curcorr;
                            curenzymes_p(j,k) = curcorrP;
                            
                             % calculate mdl coefficients and rsquared
                            mdl = fitlm(cur_expression,y);
                            curenzymes_mdlX(j,k) = mdl.Coefficients.Estimate(2);
                            curenzymes_mdlXp(j,k) = mdl.Coefficients.pValue(2);
                            curenzymes_mdlRsqadj(j,k) = mdl.Rsquared.Adjusted;
                        end
                        curenzymes_ec{j, k} = curec;
                    end
                end
                % calculate corr with sum of all enzymes
                if ~isempty(horzcat(curenzymes_geneidx{:,k}))
                    cur_expression = table2array(countsMatrixGetMM_DNA_table(horzcat(curenzymes_geneidx{:,k}),gene_idx))';
                    cur_expression = sum(cur_expression, 2);

                    if sum(cur_expression)>0
                        % corr with substrate
                        [curcorr, curcorrP] = corr(cur_expression,...
                            kegg_prod_sum_intensities(select_path_idx(i),metabolomics_idx)');
                        % add to matrix
                        kegg_DNA_sub_prod_EC_prod_sum_corr(i,k) = curcorr;
                        kegg_DNA_sub_prod_EC_prod_sum_corrP(i,k) = curcorrP;
                    end
                end
            end
            % keep max absolute value per species across all relevant enzymes
            [~, idx] = max(abs(curenzymes_mat), [], 1, 'linear');
            kegg_DNA_sub_prod_EC_prod_corr(i,:) = curenzymes_mat(idx);
            kegg_DNA_sub_prod_EC_prod_corrP(i,:) = curenzymes_p(idx);
            kegg_DNA_sub_prod_EC_prod_corr_id(i,:) = curenzymes_ec(idx);
            kegg_DNA_sub_prod_EC_prod_corr_geneidx(i,:) = curenzymes_geneidx(idx);
            kegg_DNA_sub_prod_EC_prod_mdlX(i,:) = curenzymes_mdlX(idx);
            kegg_DNA_sub_prod_EC_prod_mdlP(i,:) = curenzymes_mdlXp(idx);
            kegg_DNA_sub_prod_EC_prod_mdlRsqadj(i,:) = curenzymes_mdlRsqadj(idx);
        end
    end
end

% calculate best correlation per species across all enzymes
kegg_DNA_sub_prod_substrates_bestcorr_pos = zeros(length(kegg_sub_prod_substrates_unique),...
    size(kegg_sub_prod_EC_prod_corr,2));
kegg_DNA_sub_prod_substrates_bestcorr_neg = zeros(length(kegg_sub_prod_substrates_unique),...
    size(kegg_sub_prod_EC_prod_corr,2));
kegg_DNA_sub_prod_substrates_bestcorrP_pos = zeros(length(kegg_sub_prod_substrates_unique),...
    size(kegg_sub_prod_EC_prod_corr,2));
kegg_DNA_sub_prod_substrates_bestcorrP_neg = zeros(length(kegg_sub_prod_substrates_unique),...
    size(kegg_sub_prod_EC_prod_corr,2));
kegg_DNA_sub_prod_substrates_bestcorr_geneidx_pos = cell(length(kegg_sub_prod_substrates_unique),...
    size(kegg_sub_prod_EC_prod_corr,2));
kegg_DNA_sub_prod_substrates_bestcorr_geneidx_neg = cell(length(kegg_sub_prod_substrates_unique),...
    size(kegg_sub_prod_EC_prod_corr,2));
for i=1:length(kegg_sub_prod_substrates_unique)
    curidx = ismember(kegg_sub_prod_substrates, kegg_sub_prod_substrates_unique{i});
    curcorr = (kegg_DNA_sub_prod_EC_sub_corr.*(kegg_substrate_product_1enzyme.Dir==1)) +...
           (kegg_DNA_sub_prod_EC_prod_corr.*(kegg_substrate_product_1enzyme.Dir==-1));
    curcorr=curcorr(curidx,:);
    curcorrP = (kegg_DNA_sub_prod_EC_sub_corrP.*(kegg_substrate_product_1enzyme.Dir==1)) +...
           (kegg_DNA_sub_prod_EC_prod_corrP.*(kegg_substrate_product_1enzyme.Dir==-1));
    curcorrP=curcorrP(curidx,:);
    curgeneidx = kegg_DNA_sub_prod_EC_sub_corr_geneidx;
    
    % get info about change in LI for substrate and product
    cursubLI = sum(kegg_sub_prod_sub_LIchange(curidx,:),2)>0;
    curprodLI = sum(kegg_sub_prod_prod_LIchange(curidx,:),2)>0;
    
    %multiply corr by changes of substrate and products
    curcorr = curcorr.*repmat(cursubLI,1,size(curcorr,2)).*repmat(curprodLI,1,size(curcorr,2));

    
    for j=1:length(kegg_substrate_product_1enzyme.Dir)
        if kegg_substrate_product_1enzyme.Dir(j)==-1
            curgeneidx(j,:) = kegg_DNA_sub_prod_EC_prod_corr_geneidx(j,:);
        end
    end
    curgeneidx=curgeneidx(curidx,:);
    
    [kegg_DNA_sub_prod_substrates_bestcorr_pos(i,:),idx] = max(curcorr,[],1);
    for j=1:length(idx)
        kegg_DNA_sub_prod_substrates_bestcorrP_pos(i,j) = curcorrP(idx(j),j);
        kegg_DNA_sub_prod_substrates_bestcorr_geneidx_pos{i,j} = curgeneidx{idx(j),j};
    end
    [kegg_DNA_sub_prod_substrates_bestcorr_neg(i,:),idx] = min(curcorr,[],1);
    for j=1:length(idx)
        kegg_DNA_sub_prod_substrates_bestcorrP_neg(i,j) = curcorrP(idx(j),j);
        kegg_DNA_sub_prod_substrates_bestcorr_geneidx_neg{i,j} = curgeneidx{idx(j),j};
    end
end
% products 
kegg_DNA_sub_prod_products_bestcorr_pos = zeros(length(kegg_sub_prod_products_unique),...
    size(kegg_sub_prod_EC_prod_corr,2));
kegg_DNA_sub_prod_products_bestcorr_neg = zeros(length(kegg_sub_prod_products_unique),...
    size(kegg_sub_prod_EC_prod_corr,2));
kegg_DNA_sub_prod_products_bestcorrP_pos = zeros(length(kegg_sub_prod_products_unique),...
    size(kegg_sub_prod_EC_prod_corr,2));
kegg_DNA_sub_prod_products_bestcorrP_neg = zeros(length(kegg_sub_prod_products_unique),...
    size(kegg_sub_prod_EC_prod_corr,2));
kegg_DNA_sub_prod_products_bestcorr_geneidx_pos = cell(length(kegg_sub_prod_products_unique),...
    size(kegg_sub_prod_EC_prod_corr,2));
kegg_DNA_sub_prod_products_bestcorr_geneidx_neg = cell(length(kegg_sub_prod_products_unique),...
    size(kegg_sub_prod_EC_prod_corr,2));
for i=1:length(kegg_sub_prod_products_unique)
    curidx = ismember(kegg_sub_prod_products, kegg_sub_prod_products_unique{i});
    curcorr = (kegg_DNA_sub_prod_EC_prod_corr.*(kegg_substrate_product_1enzyme.Dir==1)) +...
           (kegg_DNA_sub_prod_EC_sub_corr.*(kegg_substrate_product_1enzyme.Dir==-1));
    curcorr=curcorr(curidx,:);
    curcorrP = (kegg_DNA_sub_prod_EC_prod_corrP.*(kegg_substrate_product_1enzyme.Dir==1)) +...
           (kegg_DNA_sub_prod_EC_sub_corrP.*(kegg_substrate_product_1enzyme.Dir==-1));
    curcorrP=curcorrP(curidx,:);
    
    curgeneidx = kegg_DNA_sub_prod_EC_prod_corr_geneidx;
    
    % get info about change in LI for substrate and product
    cursubLI = sum(kegg_sub_prod_sub_LIchange(curidx,:),2)>0;
    curprodLI = sum(kegg_sub_prod_prod_LIchange(curidx,:),2)>0;
    
    %multiply corr by changes of substrate and products
    curcorr = curcorr.*repmat(cursubLI,1,size(curcorr,2)).*repmat(curprodLI,1,size(curcorr,2));

    
    for j=1:length(kegg_substrate_product_1enzyme.Dir)
        if kegg_substrate_product_1enzyme.Dir(j)==-1
            curgeneidx(j,:) = kegg_DNA_sub_prod_EC_sub_corr_geneidx(j,:);
        end
    end
    curgeneidx=curgeneidx(curidx,:);
    
    [kegg_DNA_sub_prod_products_bestcorr_pos(i,:), idx] = max(curcorr,[],1);    
    for j=1:length(idx)
        kegg_DNA_sub_prod_products_bestcorrP_pos(i,j) = curcorrP(idx(j),j);
        kegg_DNA_sub_prod_products_bestcorr_geneidx_pos{i,j} = curgeneidx{idx(j),j};
    end
    [kegg_DNA_sub_prod_products_bestcorr_neg(i,:), idx] = min(curcorr,[],1);
    for j=1:length(idx)
        kegg_DNA_sub_prod_products_bestcorrP_neg(i,j) = curcorrP(idx(j),j);
        kegg_DNA_sub_prod_products_bestcorr_geneidx_neg{i,j} = curgeneidx{idx(j),j};
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% correlate substrates and products with species abundances
speciesOTU = readtable([inputFolderSeq, ...
    'merged_humann2_metaphlan_species_filtered_mouseDNA_OTUs.txt'], 'delim', ',');
speciesOTUabbr = cell(size(speciesOTU,1),1);
for i=1:size(speciesOTU,1)
    curname = speciesOTU.Row{i};
    curname = strsplit(curname, '_');
    speciesOTUabbr{i} = [curname{1}(1) curname{2}(1:3)];
end
% change Robe to Bobe
speciesOTUabbr(ismember(speciesOTUabbr, {'Robe'})) = {'Bobe'};
% calculate correlations
speciesOTU_products_corr = zeros(length(kegg_sub_prod_products_unique),...
                                 length(speciesOTUabbr));
speciesOTU_products_corrP = ones(length(kegg_sub_prod_products_unique),...
                                 length(speciesOTUabbr));
speciesOTU_substrates_corr = zeros(length(kegg_sub_prod_substrates_unique),...
                                 length(speciesOTUabbr));
speciesOTU_substrates_corrP = ones(length(kegg_sub_prod_substrates_unique),...
                                 length(speciesOTUabbr));
                             
for i=1:length(kegg_sub_prod_products_unique)
    curmetidx = find(ismember(kegg_sub_prod_products, ...
       kegg_sub_prod_products_unique{i}));
      
   for k=1:length(speciesOTUabbr)
        cur_otu = table2array(speciesOTU(k,2:end))';
        % corr with product
        X = cur_otu;
        if shortestPathTable.Dir(select_path_idx(curmetidx(1)))==1
            y = kegg_prod_sum_intensities(select_path_idx(curmetidx(1)),...
                                      metabolomics_idx)';  
        else
            y = kegg_sub_sum_intensities(select_path_idx(curmetidx(1)),...
                                      metabolomics_idx)'; 
        end
        [curcorr, curcorrP] = corr(X,y);
        
        % get info about change in LI for substrate and product
        cursubLI = sum(kegg_sub_prod_sub_LIchange(curmetidx(1),:),2)>0;
        curprodLI = sum(kegg_sub_prod_prod_LIchange(curmetidx(1),:),2)>0;

        %multiply corr by changes of substrate and products
        curcorr = curcorr.*cursubLI.*curprodLI;
    
        speciesOTU_products_corr(i, ismember(species_list, speciesOTUabbr{k})) = curcorr;
        speciesOTU_products_corrP(i, ismember(species_list, speciesOTUabbr{k})) = curcorrP;       
   end
end
% substrates
for i=1:length(kegg_sub_prod_substrates_unique)
    curmetidx = find(ismember(kegg_sub_prod_substrates, ...
       kegg_sub_prod_substrates_unique{i}));
   for k=1:length(speciesOTUabbr)
        cur_otu = table2array(speciesOTU(k,2:end))';
        % corr with product
        X = cur_otu;
        if shortestPathTable.Dir(select_path_idx(curmetidx(1)))==1
            y = kegg_sub_sum_intensities(select_path_idx(curmetidx(1)),...
                                      metabolomics_idx)';  
        else
            y = kegg_prod_sum_intensities(select_path_idx(curmetidx(1)),...
                                      metabolomics_idx)'; 
        end
       
        [curcorr, curcorrP] = corr(X,y);
        
         % get info about change in LI for substrate and product
        cursubLI = sum(kegg_sub_prod_sub_LIchange(curmetidx(1),:),2)>0;
        curprodLI = sum(kegg_sub_prod_prod_LIchange(curmetidx(1),:),2)>0;

        %multiply corr by changes of substrate and products
        curcorr = curcorr.*cursubLI.*curprodLI;

        
        speciesOTU_substrates_corr(i, ismember(species_list, speciesOTUabbr{k})) = curcorr;
        speciesOTU_substrates_corrP(i, ismember(species_list, speciesOTUabbr{k})) = curcorrP;       
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot correlation distributions
% plot boxplot with dots
plotnames = {'OTUfiltered', 'DNA', 'RNA'};

%%% Positive correlations with products
% plot only correlations for existing genes
filter_genes = cellfun(@(x) ~isempty(x), kegg_sub_prod_products_bestcorr_geneidx_pos(:));
plotdata = [speciesOTU_products_corr(filter_genes==1)...
            kegg_DNA_sub_prod_products_bestcorr_pos(filter_genes==1)...
            kegg_sub_prod_products_bestcorr_pos(filter_genes==1)];

%%% Positive correlations with substrates
%filter_genes = cellfun(@(x) ~isempty(x), kegg_sub_prod_substrates_bestcorr_geneidx_pos(:));
% plotdata = [speciesOTU_substrates_corr(filter_genes==1)...
%             kegg_DNA_sub_prod_substrates_bestcorr_pos(filter_genes==1)...
%             kegg_sub_prod_substrates_bestcorr_pos(filter_genes==1)];
% filter_genes = cellfun(@(x) ~isempty(x), kegg_sub_prod_substrates_bestcorr_geneidx_neg(:));

%%% Negative correlations with substrates
% plotdata = [speciesOTU_substrates_corr(filter_genes==1)...
%             kegg_DNA_sub_prod_substrates_bestcorr_neg(filter_genes==1)...
%             kegg_sub_prod_substrates_bestcorr_neg(filter_genes==1)];
% remove nan rows
plotdata(isnan(sum(plotdata,2)),:)=[];
% remove zeros in both DNA and RNA
plotdata(sum(abs(plotdata(:,1:3)),2)==0,:)=[];

figure
hold on
for i=1:size(plotdata,2)
    scatter(ones(size(plotdata,1),1).*(i+(rand(size(plotdata,1),1)-0.5)/2),...
        plotdata(:,i),'k','filled')
end
%violin(plotdata,plotnames)
boxplot(plotdata,plotnames, 'Notch','on')
set(gca, 'XTick', 1:length(plotnames))
set(gca, 'XTickLabel', plotnames)

ylim([-1 1])
for i=1:size(plotdata,2)
    text(i-0.5, -0.8, sprintf('Mean %.2f', mean(plotdata(:,i))))
    text(i-0.5, -0.95, sprintf('Median %.2f', median(plotdata(:,i))))
end
idx=1;
for i=1:size(plotdata,2)-1
    for j=i+1:size(plotdata,2)
        % compare distributions with Wilcoxon signed rank test
        curp = signrank(plotdata(:,i), plotdata(:,j));
        text(idx-0.5, -0.7, sprintf('sgrk %d %d p %.2e', i,j,curp))
        idx=idx+1;
    end
end
axis square
title('Best product positive corr')
print(gcf, '-vector', '-dpdf', '-r600', '-bestfit', ...
    [figureFolder ...
    'fig_5c_boxplotnotch_best_combinedsol_strictclass'...
    num2str(strict_class), '_product_corrPOS_otuFiltered_DNA_RNA_LIchange_sub_prod_signrank.pdf'])
% print(gcf, '-vector', '-dpdf', '-r600', '-bestfit', ...
%     [figureFolder ...
%     'fig_5c_violin_best_product_corrPOS_otuFiltered_DNA_RNA_LIchange_sub_prod_signrank.pdf'])
% title('Best substrate positive corr')
% print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
%     [figureFolder ...
%      'fig_sup_violin_best_substrate_corrPOS_otuFiltered_DNA_RNA_LIchange_sub_prod_signrank.pdf'])
% title('Best substrate negative corr')
% print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
%     [figureFolder ...
%       'fig_sup_violin_best_substrate_corrNEG_otuFiltered_DNA_RNA_LIchange_sub_prod_signrank.pdf'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot heatmaps for OTU, DNA and RNA for genes passing DNA and RNA criteria

exclude_products = {'C00022', 'C00041'};%, 'C00245'};
grayredcmap = [[linspace(200, 247, 5)', linspace(200, 247, 5)', linspace(200, 247,5)'];...
               [linspace(247, 130, 6)', linspace(247, 0, 6)', linspace(247, 0, 6)']]/256;


plotdataRNA = kegg_sub_prod_products_bestcorr_pos.*...
             (kegg_sub_prod_products_bestcorrP_pos<0.1);

%DNA
plotdataDNA = kegg_DNA_sub_prod_products_bestcorr_pos.*...
             (kegg_DNA_sub_prod_products_bestcorrP_pos<0.1);
% OTU
plotdataOTU = speciesOTU_products_corr.*...
             (speciesOTU_products_corrP<0.1);
% filter by RNA only
%keep_pairs = (sum(plotdataRNA>0,2)>0);
% filter by RNA or DNA
keep_pairs = ((sum(plotdataRNA>0,2)>0) | (sum(plotdataDNA>0,2)>0) &...
    ~ismember(kegg_sub_prod_products_unique, exclude_products));

plotdataRNA = kegg_sub_prod_products_bestcorr_pos;
plotdataDNA = kegg_DNA_sub_prod_products_bestcorr_pos;
plotdataOTU = speciesOTU_products_corr;

% add nonexisting genes
plotdataRNA(cellfun(@(x) isempty(x), kegg_sub_prod_products_bestcorr_geneidx_pos))=100;
plotdataDNA(cellfun(@(x) isempty(x), kegg_sub_prod_products_bestcorr_geneidx_pos))=100;

plotdataRNA = plotdataRNA(keep_pairs,:);
plotdataDNA = plotdataDNA(keep_pairs,:);
plotdataOTU = plotdataOTU(keep_pairs,:);

plotdata_genes = kegg_sub_prod_products_bestcorr_geneidx_pos(keep_pairs,:);
plotdata_rows = kegg_sub_prod_products_unique(keep_pairs);
% keep only positive corr
plotdataDNA(plotdataDNA<0)=0;
plotdataRNA(plotdataRNA<0)=0;
plotdataOTU(plotdataOTU<0)=0;

plotdataDNA(plotdataDNA==100)=-1;
plotdataRNA(plotdataRNA==100)=-1;

cgo = clustergram(plotdataRNA, ...
    'columnLabels', species_list,...
    'rowLabels', plotdata_rows,...
    'colormap', grayredcmap,...
    'DisplayRange', 10);
species_order_cluster = cgo.ColumnLabels;%

fig = cgo.plot;
print(gcf, '-vector', '-dpdf', '-r600', '-bestfit', ...
    [figureFolder ...
    'fig_5d_clustergam_best_combined_strictclass'...
    num2str(strict_class) '_product_corrPOS_RNA_filteredbyRNADNA_allpos_withtaurine.pdf'])


cgo = clustergram(plotdataDNA, ...
    'columnLabels', species_list,...
    'rowLabels', plotdata_rows,...
    'colormap', grayredcmap);

fig = cgo.plot;
print(gcf, '-vector', '-dpdf', '-r600', '-bestfit', ...
    [figureFolder ...
    'fig_sup_clustergam_best_combined_strictclass'...
    num2str(sctrict_class) '_product_corrPOS_DNA_filteredbyRNADNA_allpos_withtaurine.pdf'])

cgo = clustergram(plotdataOTU, ...
    'columnLabels', species_list,...
    'rowLabels', plotdata_rows,...
    'colormap', grayredcmap);

fig = cgo.plot;
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
    [figureFolder ...
    'fig_sup_clustergam_best_combined_strictclass'...
    num2str(strict_class), '_product_corrPOS_OTU_filteredbyRNADNA_allpos_withtaurine.pdf'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for each selected product, plot heatmaps of corr with OTU, DNA and RNA
%figureFile = 'heatmap_products_filtered_byRNAandDNA_significant0_1.ps';
figureFile = ['fig_5e_sup_heatmap_combined_strictclass' num2str(strict_class)...
              '_products_filtered_byRNAandDNA_onlyPosCorr_significant0_1_sortedX.ps'];
[~,~,plot_order] = intersect(species_order_cluster, species_list, 'stable');
for i=1:size(plotdataRNA,1)
    fig = figure;
    heatmap(species_order_cluster,...
            {'OTU', 'DNA', 'RNA'},...
            [plotdataOTU(i,plot_order); plotdataDNA(i,plot_order); plotdataRNA(i,plot_order)])
    colormap(grayredcmap)
    clim([-1 1])
    title(plotdata_rows{i})
    orient landscape 
    print(gcf, '-vector', '-dpsc2', '-r600', '-append', '-bestfit',...
       [figureFolder, figureFile]);
    close(fig)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot correlation between species enzymes and metabolites for each
% metabolite with significant correlation
figureFile = ['fig_sup_scatter_combined_strictclass' num2str(strict_class)...
    '_metabolites_vs_best_enzymes_per_species_pos_pRNA_or_pRNA_0_1.ps'];
for i=1:length(plotdata_rows)
   curmetidx = find(ismember(kegg_sub_prod_products, plotdata_rows{i}));
   fig = figure('units','normalized','outerposition',[0 0 1 1]);
   spidx=1;   
   for k=1:size(plotdata_genes,2)
        if ~isempty(plotdata_genes{i,k})
            cur_expression = table2array(countsMatrixGetMM_RNA_table(...
                plotdata_genes{i,k},gene_idx))';
            cur_expression = sum(cur_expression, 2);

            if sum(cur_expression)>0
                % corr with substrate
                X = cur_expression;
                if shortestPathTable.Dir(select_path_idx(curmetidx(1)))==1
                    y = kegg_prod_sum_intensities(select_path_idx(curmetidx(1)),...
                                              metabolomics_idx)';  
                else
                    y = kegg_sub_sum_intensities(select_path_idx(curmetidx(1)),...
                                              metabolomics_idx)'; 
                end
                % remove min and max values
                
                [curcorr, curcorrP] = corr(X,y);
                [curcorrS] = corr(X,y, 'type','spearman');                          
                                        

                mdl = fitlm(X,y);
                subplot(3,5,spidx);
                scatter(X,y);
                axis square
                title(sprintf('%s %.2f %.2f %.2f %.2f %.2f',species_list{k},...
                                             plotdataRNA(i,k),...
                                             curcorr, curcorrP, curcorrS,...
                                             mdl.Rsquared.Adjusted))
            end
        end
        spidx = spidx+1;
   end
   sgtitle(sprintf('%s (max corr, corr, corrP, corrS, Rsqadj)',plotdata_rows{i}));
   orient landscape 
   print(gcf, '-painters', '-dpsc2', '-r600', '-append', '-bestfit',...
       [figureFolder, figureFile]);
   close(fig)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot correlation between species enzymes and metabolites for each
% metabolite with significant correlation
figureFile = ['fig_sup_scatter_combined_strictclass' num2str(strict_class),...
    '_metabolites_vs_best_enzymes_DNA_per_species.ps'];
for i=1:length(plotdata_rows)
   curmetidx = find(ismember(kegg_sub_prod_products, plotdata_rows{i}));
   fig = figure('units','normalized','outerposition',[0 0 1 1]);
   spidx=1;   
   for k=1:size(plotdata_genes,2)
        if ~isempty(plotdata_genes{i,k})
            cur_expression = table2array(countsMatrixGetMM_DNA_table(...
                plotdata_genes{i,k},gene_idx))';
            cur_expression = sum(cur_expression, 2);

            if sum(cur_expression)>0
                % corr with substrate
                X = cur_expression;
                if shortestPathTable.Dir(select_path_idx(curmetidx(1)))==1
                    y = kegg_prod_sum_intensities(select_path_idx(curmetidx(1)),...
                                              metabolomics_idx)';  
                else
                    y = kegg_sub_sum_intensities(select_path_idx(curmetidx(1)),...
                                              metabolomics_idx)'; 
                end
                [curcorr, curcorrP] = corr(X,y);
                [curcorrS] = corr(X,y, 'type','spearman');                          
                                        

                mdl = fitlm(X,y);
                subplot(3,5,spidx);
                scatter(X,y);
                axis square
                title(sprintf('%s %.2f %.2f %.2f %.2f %.2f',species_list{k},...
                                             plotdataDNA(i,k),...
                                             curcorr, curcorrP, curcorrS,...
                                             mdl.Rsquared.Adjusted))
            end
        end
        spidx = spidx+1;
   end
   sgtitle(sprintf('%s DNA (max corr, corr, corrP, corrS, Rsqadj)',plotdata_rows{i}));
   orient landscape 
   print(gcf, '-vector', '-dpsc2', '-r600', '-append', '-bestfit',...
       [figureFolder, figureFile]);
   close(fig)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% include data from CV mice as well
%figureFile = 'fig_sup_scatter_metabolites_vs_best_enzymes_DNA_per_species.ps';
%for i=1:length(plotdata_rows)
i=9; % porph
   curmetidx = find(ismember(kegg_sub_prod_products, plotdata_rows{i}));
   fig = figure('units','normalized','outerposition',[0 0 1 1]);
   spidx=1;   
   for k=1:size(plotdata_genes,2)
        if ~isempty(plotdata_genes{i,k})
            cur_expression = table2array(countsMatrixGetMM_DNA_table(...
                plotdata_genes{i,k},gene_idx))';
            cur_expression = sum(cur_expression, 2);

            if sum(cur_expression)>0
                % corr with product
                X = cur_expression;
                if shortestPathTable.Dir(select_path_idx(curmetidx(1)))==1
                    y = kegg_prod_sum_intensities(select_path_idx(curmetidx(1)),...
                                              metabolomics_idx)';  
                else
                    y = kegg_sub_sum_intensities(select_path_idx(curmetidx(1)),...
                                              metabolomics_idx)'; 
                end
                [curcorr, curcorrP] = corr(X,y);
                [curcorrS] = corr(X,y, 'type','spearman');                          
                                        

                mdl = fitlm(X,y);
                subplot(3,5,spidx);
                scatter(X,y);
                axis square
                title(sprintf('%s %.2f %.2f %.2f %.2f %.2f',species_list{k},...
                                             plotdataDNA(i,k),...
                                             curcorr, curcorrP, curcorrS,...
                                             mdl.Rsquared.Adjusted))
            end
        end
        spidx = spidx+1;
   end
   sgtitle(sprintf('%s DNA (max corr, corr, corrP, corrS, Rsqadj)',plotdata_rows{i}));

% plot shortbred results
%shortbred_genes = unique(countsMatrix_shortBRED.Family);
%shortbred_columns = cellfun(@(x) contains(x, 'hits'), countsMatrix_shortBRED.Properties.VariableNames);
%%
% shortbred_genes = unique(countsMatrix_shortBRED_DNA.Var1); 
% shortbred_columns = cellfun(@(x) contains(x, 'MSZ'), countsMatrix_shortBRED_DNA.Properties.VariableNames);
%%
shortbred_dna_flag = 1;
if shortbred_dna_flag ==1
    shortbred_genes = unique(countsMatrix_shortBRED_DNA.Family); 
    shortbred_columns = cellfun(@(x) contains(x, 'MSZ'), countsMatrix_shortBRED_DNA.Properties.VariableNames);
else
    shortbred_genes = unique(countsMatrix_shortBRED_RNA.Family); 
    shortbred_columns = cellfun(@(x) contains(x, 'MSZ'), countsMatrix_shortBRED_RNA.Properties.VariableNames);
end
% allocate row for each unique protein and calculate sum
plotdata_shortbred = zeros(length(shortbred_genes), nnz(shortbred_columns));
for i=1:length(shortbred_genes)
    %curcounts = countsMatrix_shortBRED_DNA(ismember(countsMatrix_shortBRED_DNA.Var1, shortbred_genes{i}),:);
    %curcounts = countsMatrix_shortBRED(ismember(countsMatrix_shortBRED.Family, shortbred_genes{i}),:);
    if shortbred_dna_flag==1
        curcounts = countsMatrix_shortBRED_DNA(ismember(countsMatrix_shortBRED_DNA.Family, shortbred_genes{i}),:);
    else
        curcounts = countsMatrix_shortBRED_RNA(ismember(countsMatrix_shortBRED_RNA.Family, shortbred_genes{i}),:);
    end
    plotdata_shortbred(i,:) = sum(table2array(curcounts(:,shortbred_columns)),1);
end
%remove zero rows
shortbred_genes(sum(plotdata_shortbred,2)==0)=[];
plotdata_shortbred(sum(plotdata_shortbred,2)==0,:)=[];
%shortbred_columns = countsMatrix_shortBRED.Properties.VariableNames(shortbred_columns);
if shortbred_dna_flag==1
    shortbred_columns = countsMatrix_shortBRED_DNA.Properties.VariableNames(shortbred_columns);
else
    shortbred_columns = countsMatrix_shortBRED_RNA.Properties.VariableNames(shortbred_columns);
end
shortbred_columns = cellfun(@(x) strrep(strrep(x, 'Count_', ''), '_results',''),...
                shortbred_columns, 'unif', 0);
% match shortbred columns with metabolomics columns
[~, metabolomics_shortbred_idx, shortbred_columns_idx] = intersect(metabolomics_samples, shortbred_columns);

% reformat shortbred_genes
shortbred_genes = cellfun(@(x) strrep(x, "['", ""), shortbred_genes, 'unif', 0);
shortbred_genes = cellfun(@(x) strrep(x, "']", ""), shortbred_genes, 'unif', 0);
shortbred_genes = cellfun(@(x) x{1}, shortbred_genes, 'unif', 0);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot for each metablite product shortbred results
curenzymes_locustags_combined=[];
curenzymes_locustags_total_num = 0;

for i=1:length(plotdata_rows)
%i = 7; 
    
    curmetidx = find(ismember(kegg_sub_prod_products, plotdata_rows{i}));
    
    curenzymes_geneidx = plotdata_genes(i,:);
    % remove empty values and keep only indeces of existing genes
    curenzymes_geneidx(cellfun(@(x) isempty(x), curenzymes_geneidx)) = [];
    curenzymes_geneidx = cell2mat(curenzymes_geneidx);
    curenzymes_geneinfo = geneTable(curenzymes_geneidx,:);
    curenzymes_locustags = curenzymes_geneinfo.locus_tag;
    % reformat locutags to cell array of strings without special symbols
    curenzymes_locustags = cellfun(@(x) strrep(x, "['", ""), curenzymes_locustags, 'unif', 0);
    curenzymes_locustags = cellfun(@(x) strrep(x, "']", ""), curenzymes_locustags, 'unif', 0);
    curenzymes_locustags = cellfun(@(x) x{1}, curenzymes_locustags, 'unif', 0);
    
    curenzymes_locustags_total_num = curenzymes_locustags_total_num + length(curenzymes_locustags);

    
    sprintf('%s expected enzymes: %d,  found enzymes: %d\n', ...
        plotdata_rows{i}, ...
        length(curenzymes_locustags),...
        length(intersect(shortbred_genes, curenzymes_locustags)))
    curenzymes_locustags_combined = [curenzymes_locustags_combined; ...
        [repmat(plotdata_rows(i),length(curenzymes_locustags),1), curenzymes_locustags]];
end

curenzymes_locustags_combined = table(curenzymes_locustags_combined(:,1),...
                                      curenzymes_locustags_combined(:,2),...
                                      'VariableNames', {'KEGGID_product', 'locus_tag'});
writetable(curenzymes_locustags_combined, ...
    '.\ProcessedData\shortBRED\locustags_for_plotting_products.txt');

%curenzymes_locustags_total_num =
%   363

% leave only sample names
%shortbred_columns = cellfun(@(x)x(1:strfind(x,'_')-1), shortbred_columns, 'unif',0);
min_mice_num = 5; % minimum number of mice to consider correlation

if shortbred_dna_flag==1
    figureFile = ['fig_sup_scatter_metabolites_vs_best_enzymes_shortbred_DNA_per_species_minmice'...
                   num2str(min_mice_num) '.ps'];
else
    figureFile = ['fig_sup_scatter_metabolites_vs_best_enzymes_shortbred_RNA_per_species_minmice'...
                   num2str(min_mice_num) '.ps'];
end

% plot correlation for each enzyme
% save information about metabolite-enzyme correlation for all mice, DC
% mice and CVR mice
select_path_idx_orig = select_path_idx;

%nrand_perm = 1000;
%corr_shortbred_rand = zeros(nrand_perm, 1);

%for rand_i = 1:nrand_perm
%    select_path_idx = select_path_idx_orig(randperm(length(select_path_idx_orig)));         
select_path_idx = select_path_idx_orig;

met_enz_PCCall = zeros(height(curenzymes_locustags_combined),1);
met_enz_PCCPall = zeros(height(curenzymes_locustags_combined),1);
met_enz_SCCall = zeros(height(curenzymes_locustags_combined),1);
met_enz_SCCPall = zeros(height(curenzymes_locustags_combined),1);
met_enz_RSQall = zeros(height(curenzymes_locustags_combined),1);

met_enz_PCCdc = zeros(height(curenzymes_locustags_combined),1);
met_enz_PCCPdc = zeros(height(curenzymes_locustags_combined),1);
met_enz_SCCdc = zeros(height(curenzymes_locustags_combined),1);
met_enz_SCCPdc = zeros(height(curenzymes_locustags_combined),1);
met_enz_RSQdc = zeros(height(curenzymes_locustags_combined),1);

met_enz_PCCcvr = zeros(height(curenzymes_locustags_combined),1);
met_enz_PCCPcvr = zeros(height(curenzymes_locustags_combined),1);
met_enz_SCCcvr = zeros(height(curenzymes_locustags_combined),1);
met_enz_SCCPcvr = zeros(height(curenzymes_locustags_combined),1);
met_enz_RSQcvr = zeros(height(curenzymes_locustags_combined),1);

met_enz_NMETall = zeros(height(curenzymes_locustags_combined),1);
met_enz_NENZPall = zeros(height(curenzymes_locustags_combined),1);
met_enz_NMETdc = zeros(height(curenzymes_locustags_combined),1);
met_enz_NENZdc = zeros(height(curenzymes_locustags_combined),1);
met_enz_NMETcvr = zeros(height(curenzymes_locustags_combined),1);
met_enz_NENZcvr = zeros(height(curenzymes_locustags_combined),1);

met_enz_metabolite = cell(height(curenzymes_locustags_combined),1);
met_enz_enzyme = cell(height(curenzymes_locustags_combined),1);

% save metabolite and enzyme data
met_enz_metabolite_data = zeros(height(curenzymes_locustags_combined),20);
met_enz_enzyme_data = zeros(height(curenzymes_locustags_combined),20);

met_enz_idx=1;

plot_shortbred_flag = 0;

if plot_shortbred_flag
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
end

for i=1:length(plotdata_rows)
    curmetidx = find(ismember(kegg_sub_prod_products, plotdata_rows{i}));
    
    curenzymes_geneidx = plotdata_genes(i,:);
    % remove empty values and keep only indeces of existing genes
    curenzymes_geneidx(cellfun(@(x) isempty(x), curenzymes_geneidx)) = [];
    curenzymes_geneidx = cell2mat(curenzymes_geneidx);
    curenzymes_geneinfo = geneTable(curenzymes_geneidx,:);
    curenzymes_locustags = curenzymes_geneinfo.locus_tag;
    % reformat locutags to cell array of strings without special symbols
    curenzymes_locustags = cellfun(@(x) strrep(x, "['", ""), curenzymes_locustags, 'unif', 0);
    curenzymes_locustags = cellfun(@(x) strrep(x, "']", ""), curenzymes_locustags, 'unif', 0);
    curenzymes_locustags = cellfun(@(x) x{1}, curenzymes_locustags, 'unif', 0);
   
    spidx=1;   
    for k=1:length(curenzymes_locustags)%size(plotdata_shortbred,1)%
        %cur_expression = plotdata_shortbred(k,shortbred_columns_idx)';
        shortbred_row_idx = ismember(shortbred_genes, curenzymes_locustags{k});
        cur_expression = plotdata_shortbred(shortbred_row_idx,shortbred_columns_idx)';
        
        cur_expression = sum(cur_expression, 2);
    
        if sum(cur_expression)>0

            % corr with product
            X = cur_expression;
            if shortestPathTable.Dir(select_path_idx(curmetidx(1)))==1
                y = kegg_prod_sum_intensities(select_path_idx(curmetidx(1)),...
                                          metabolomics_shortbred_idx)';  
            else
                y = kegg_sub_sum_intensities(select_path_idx(curmetidx(1)),...
                                          metabolomics_shortbred_idx)'; 
            end
          
            if (nnz(X(11:20))>=min_mice_num) &&  (nnz(y(11:20)>5000)>=min_mice_num)

                met_enz_metabolite{met_enz_idx} = plotdata_rows{i};
                met_enz_enzyme{met_enz_idx} = curenzymes_locustags{k};

                met_enz_metabolite_data(met_enz_idx, :) = X;
                met_enz_enzyme_data(met_enz_idx, :) = y; 

                % plot full dataset
                [curcorr, curcorrP] = corr(X,y);
                [curcorrS, curcorrSP] = corr(X,y, 'type','spearman'); 
                mdl = fitlm(X,y); % fit linear model
    
                if plot_shortbred_flag
                    subplot(3,5,spidx);
                    hold on
                    scatter(X,y, 'filled');
                    if shortbred_dna_flag
                        xlabel('All mice, DNA')  
                    else
                        xlabel('All mice, RNA')
                    end
    
                    axis square
                    title({shortbred_genes{shortbred_row_idx},...
                           sprintf('%.2f %.2f %.2f %.2f %.2f',curcorr, curcorrP,...
                                                 curcorrS,curcorrSP),...
                           sprintf('RSQ %.2f', mdl.Rsquared.Adjusted)}, ...
                           'Interpreter', 'none')
                end
                met_enz_PCCall(met_enz_idx) = curcorr;
                met_enz_PCCPall(met_enz_idx) = curcorrP;
                met_enz_SCCall(met_enz_idx) = curcorrS;
                met_enz_SCCPall(met_enz_idx) = curcorrSP;
                met_enz_RSQall(met_enz_idx) = mdl.Rsquared.Adjusted;
                met_enz_NMETall(met_enz_idx) = nnz(y);
                met_enz_NENZPall(met_enz_idx) = nnz(X);

                % plot only DC mice

                [curcorr, curcorrP] = corr(X(1:10),y(1:10));
                [curcorrS, curcorrSP] = corr(X(1:10),y(1:10), 'type','spearman'); 
                mdl = fitlm(X(1:10),y(1:10)); % fit linear model
    
                if plot_shortbred_flag
                    subplot(3,5,5+spidx);
                    hold on
                    scatter(X(1:10),y(1:10), 'filled');
                    if shortbred_dna_flag
                        xlabel('DC mice, DNA')  
                    else
                        xlabel('DC mice, RNA')
                    end
    
                    axis square
                    title({shortbred_genes{shortbred_row_idx},...
                           sprintf('%.2f %.2f %.2f %.2f %.2f',curcorr, curcorrP,...
                                                 curcorrS,curcorrSP),...
                           sprintf('RSQ %.2f', mdl.Rsquared.Adjusted)}, ...
                           'Interpreter', 'none')
                end

                met_enz_PCCdc(met_enz_idx) = curcorr;
                met_enz_PCCPdc(met_enz_idx) = curcorrP;
                met_enz_SCCdc(met_enz_idx) = curcorrS;
                met_enz_SCCPdc(met_enz_idx) = curcorrSP;
                met_enz_RSQdc(met_enz_idx) = mdl.Rsquared.Adjusted;
                met_enz_NMETdc(met_enz_idx) = nnz(y(1:10));
                met_enz_NENZdc(met_enz_idx) = nnz(X(1:10));

                % plot only CVR mice
                [curcorr, curcorrP] = corr(X(11:20),y(11:20));
                [curcorrS, curcorrSP] = corr(X(11:20),y(11:20), 'type','spearman'); 
                mdl = fitlm(X(11:20),y(11:20)); % fit linear model
    
                if plot_shortbred_flag
                    subplot(3,5,10+spidx);
                    hold on
                    scatter(X(11:20),y(11:20), 'filled');
                    if shortbred_dna_flag
                        xlabel('CVR mice, DNA')  
                    else
                        xlabel('CVR mice, RNA')
                    end
    
                    axis square
                    title({shortbred_genes{shortbred_row_idx},...
                           sprintf('%.2f %.2f %.2f %.2f',curcorr, curcorrP,...
                                                 curcorrS,curcorrSP),...
                           sprintf('RSQ %.2f', mdl.Rsquared.Adjusted)}, ...
                           'Interpreter', 'none')
                end
                met_enz_PCCcvr(met_enz_idx) = curcorr;
                met_enz_PCCPcvr(met_enz_idx) = curcorrP;
                met_enz_SCCcvr(met_enz_idx) = curcorrS;
                met_enz_SCCPcvr(met_enz_idx) = curcorrSP;
                met_enz_RSQcvr(met_enz_idx) = mdl.Rsquared.Adjusted;
                met_enz_NMETcvr(met_enz_idx) = nnz(y(11:20));
                met_enz_NENZcvr(met_enz_idx) = nnz(X(11:20));
              
                spidx = spidx+1;
                met_enz_idx = met_enz_idx+1;
            end
        end
        if spidx>5
            if shortbred_dna_flag
               sgtitle(sprintf('%s shortbred DNA (PCC, PCCp, SCC, SCCp, Rsqadj)',plotdata_rows{i}));
            else
               sgtitle(sprintf('%s shortbred RNA (PCC, PCCp, SCC, SCCp, Rsqadj)',plotdata_rows{i}));
            end 
  
            if plot_shortbred_flag
                print(fig, '-vector', '-dpsc2', '-r600', '-append', '-bestfit',...
                                            [figureFolder, figureFile]);
                clf(fig)
                spidx=1;
            end
        end
    end
   if shortbred_dna_flag
       sgtitle(sprintf('%s shortbred DNA (PCC, PCCp, SCC, SCCp, Rsqadj)',plotdata_rows{i}));
   else
       sgtitle(sprintf('%s shortbred RNA (PCC, PCCp, SCC, SCCp, Rsqadj)',plotdata_rows{i}));
   end    
   if plot_shortbred_flag
       orient landscape 
       print(fig, '-vector', '-dpsc2', '-r600', '-append', '-bestfit',...
           [figureFolder, figureFile]);
       clf(fig)
   end
end

met_enz_PCCall(met_enz_idx:end) = [];
met_enz_PCCPall(met_enz_idx:end) = [];
met_enz_SCCall(met_enz_idx:end) = [];
met_enz_SCCPall(met_enz_idx:end) = [];
met_enz_RSQall(met_enz_idx:end) = [];

met_enz_PCCdc(met_enz_idx:end) = [];
met_enz_PCCPdc(met_enz_idx:end) = [];
met_enz_SCCdc(met_enz_idx:end) = [];
met_enz_SCCPdc(met_enz_idx:end) = [];
met_enz_RSQdc(met_enz_idx:end) = [];

met_enz_PCCcvr(met_enz_idx:end) = [];
met_enz_PCCPcvr(met_enz_idx:end) = [];
met_enz_SCCcvr(met_enz_idx:end) = [];
met_enz_SCCPcvr(met_enz_idx:end) = [];
met_enz_RSQcvr(met_enz_idx:end) = [];

met_enz_NMETall(met_enz_idx:end) = [];
met_enz_NENZPall(met_enz_idx:end) = [];
met_enz_NMETdc(met_enz_idx:end) = [];
met_enz_NENZdc(met_enz_idx:end) = [];
met_enz_NMETcvr(met_enz_idx:end) = [];
met_enz_NENZcvr(met_enz_idx:end) = [];

met_enz_metabolite(met_enz_idx:end) = [];
met_enz_enzyme(met_enz_idx:end) = [];

met_enz_metabolite_data(met_enz_idx:end,:) = [];
met_enz_enzyme_data(met_enz_idx:end,:) = [];

%corr_shortbred_rand(rand_i) = corr(met_enz_PCCdc, met_enz_PCCcvr);
%end

figure
histogram(corr_shortbred_rand)
hold on
plot([corrC, corrC], [1, 20], 'g')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot correlations per metabolite
figure
scatter(met_enz_PCCdc, met_enz_PCCcvr, 'filled')
mdl = fitlm(met_enz_PCCdc, met_enz_PCCcvr); % fit linear model
hold on
plot([-1 1], [-1 1]*mdl.Coefficients.Estimate(2)+mdl.Coefficients.Estimate(1), 'k--')
xlim([-1 1])
ylim([-1 1])
axis square
xlabel('PCC between enzyme and product, DC')
ylabel('PCC between enzyme and product, CVR')
[corrC, corrP] = corr(met_enz_PCCdc, met_enz_PCCcvr);
[corrSC, corrSP] = corr(met_enz_PCCdc, met_enz_PCCcvr, 'type', 'spearman');

if shortbred_dna_flag   
    title({['DNA minmicenum=' num2str(min_mice_num)],...
           sprintf('PCC %.2f %.2f SCC %.2f %.2f', corrC, corrP, corrSC, corrSP)})
    figureFile = 'FigX_scatter_correnzprodDC_vs_CVR_DNA.pdf';
else
    title({['RNA minmicenum=' num2str(min_mice_num)],...
           sprintf('PCC %.2f %2f SCC %.2f %.2f', corrC, corrP, corrSC, corrSP)})
    figureFile = 'FigX_scatter_correnzprodDC_vs_CVR_RNA.pdf';
end
print(gcf, '-vector', '-dpdf', '-r600', '-bestfit',...
       [figureFolder, figureFile]);

% permute which CV enzyme is correlated to the corresponding metabolite
% (calculate correlation with a random CV enzyme)
nperm = 100;
testcorr = zeros(nperm,1);                
figure
hold on
for i=1:nperm
    met_enz_enzyme_data_shuffled = [met_enz_enzyme_data(:,1:10),...
                                met_enz_enzyme_data(randperm(size(met_enz_enzyme_data,1)),11:20)];      
    cor_dc = zeros(size(met_enz_enzyme_data_shuffled,1),1);
    cor_cvr = zeros(size(met_enz_enzyme_data_shuffled,1),1);
    
    for j=1:size(met_enz_enzyme_data_shuffled,1)
        cor_dc(j) = corr(met_enz_metabolite_data(j,1:10)', met_enz_enzyme_data_shuffled(j,1:10)');
        cor_cvr(j) = corr(met_enz_metabolite_data(j,11:20)', met_enz_enzyme_data_shuffled(j,11:20)');
    end
        
    testcorr(i) = corr(cor_dc, cor_cvr);
    scatter(cor_dc, cor_cvr, 'filled', 'MarkerFaceColor', [.5 .5 .5])
end
figure
histogram(testcorr)
hold on
plot([corrC, corrC], [1, 20], 'g')

% permute correlation coefficients (correlate random pairs of
% enzyme-metabolites that belong together)
nperm = 100;
testcorr = zeros(nperm,1);
hold on
for i=1:nperm
    met_enz_shuffled_dc = met_enz_PCCdc(randperm(length(met_enz_PCCdc)));
    met_enz_shuffled = met_enz_PCCcvr(randperm(length(met_enz_PCCcvr)));
    testcorr(i) = corr(met_enz_shuffled_dc, met_enz_shuffled);
    scatter(met_enz_shuffled_dc, met_enz_shuffled, 'filled', 'MarkerFaceColor', [.5 .5 .5])
end

figure
histogram(testcorr)
hold on
plot([corrC, corrC], [1, 20], 'g')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig = figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:length(plotdata_rows)
    curmetidx = find(ismember(met_enz_metabolite, plotdata_rows{i}));
    
    subplot(3,4,i)
    scatter(met_enz_PCCcvr, met_enz_PCCdc)
    hold on
    scatter(met_enz_PCCcvr(curmetidx), met_enz_PCCdc(curmetidx), 'r')

    title(plotdata_rows{i})
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot correlation between species enzymes and metabolites for each
% metabolite with significant correlation
figureFile = 'fig_sup_scatter_metabolites_vs_OTU_DNA_per_species.ps';
for i=1:length(plotdata_rows)
   curmetidx = find(ismember(kegg_sub_prod_products, plotdata_rows{i}));
   fig = figure('units','normalized','outerposition',[0 0 1 1]);
   spidx=1;   
   for k=1:size(plotdata_genes,2)
        X= table2array(speciesOTU(k,2:end))';
        
        if shortestPathTable.Dir(select_path_idx(curmetidx(1)))==1
            y = kegg_prod_sum_intensities(select_path_idx(curmetidx(1)),...
                                      metabolomics_idx)';  
        else
            y = kegg_sub_sum_intensities(select_path_idx(curmetidx(1)),...
                                      metabolomics_idx)'; 
        end
       
        [curcorr, curcorrP] = corr(X,y);
        [curcorrS] = corr(X,y, 'type','spearman');                          


        mdl = fitlm(X,y);
        subplot(3,5,spidx);
        scatter(X,y);
        axis square
        title(sprintf('%s %.2f %.2f %.2f %.2f',species_list{k},...
                                     curcorr, curcorrP, curcorrS,...
                                     mdl.Rsquared.Adjusted))
    
        spidx = spidx+1;
   end
   sgtitle(sprintf('%s OTU DNA (corr, corrP, corrS, Rsqadj)',plotdata_rows{i}));
   orient landscape 
   print(gcf, '-vector', '-dpsc2', '-r600', '-append', '-bestfit',...
       [figureFolder, figureFile]);
   close(fig)
end
 
% % for each pair, calculate in how many species it is present
% kegg_sub_prod_ec_count = sum(kegg_sub_prod_ECchange~=0,2);
% plotdata = kegg_sub_prod_ec_count(kegg_sub_prod_ec_count~=0);
% figure
% histogram(plotdata)
% xlabel('Number of species with enzyme')
% ylabel('Number of sub-prod pairs')
% axis square
% print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
%     'histogram_number_of_species_single_pair_enzyme.pdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate number of products correlating with enzymes
nnz(sum(kegg_sub_prod_products_bestcorrP_pos<0.1,2))
% calculate mz for each product to calculate unique ions
kegg_sub_prod_products_unique_mz = zeros(size(kegg_sub_prod_products_unique));
for i=1:length(kegg_sub_prod_products_unique)
    curidx = kegg_ids_unique_index(ismember(kegg_ids_unique, kegg_sub_prod_products_unique{i}));
    kegg_sub_prod_products_unique_mz(i) = metaboliteFilters.MZ(curidx);
end
% calculate p-values for unique ions
kegg_sub_prod_products_unique_mz_unique = unique(kegg_sub_prod_products_unique_mz);
kegg_sub_prod_products_unique_mz_unique_sigcorr = zeros(size(kegg_sub_prod_products_unique_mz_unique));
for i=1:length(kegg_sub_prod_products_unique_mz_unique)
    curidx = (kegg_sub_prod_products_unique_mz == kegg_sub_prod_products_unique_mz_unique(i));
    kegg_sub_prod_products_unique_mz_unique_sigcorr(i) = max(sum(kegg_sub_prod_products_bestcorrP_pos(curidx ,:)<0.1,2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get unique substrate mz
% calculate mz for each product to calculate unique ions
kegg_sub_prod_substrates_unique_mz = zeros(size(kegg_sub_prod_substrates_unique));
for i=1:length(kegg_sub_prod_substrates_unique)
    curidx = kegg_ids_unique_index(ismember(kegg_ids_unique, kegg_sub_prod_substrates_unique{i}));
    kegg_sub_prod_substrates_unique_mz(i) = metaboliteFilters.MZ(curidx);
end
length(unique(kegg_sub_prod_substrates_unique_mz))
% ans = 70
length(unique(kegg_sub_prod_products_unique_mz))
% ans = 66
nnz(kegg_sub_prod_products_unique_mz_unique_sigcorr)
% ans = 9
nnz(kegg_sub_prod_products_unique_mz_unique_sigcorr)/length((kegg_sub_prod_products_unique_mz_unique_sigcorr))
% ans = 0.1364

% get unique pathes and change substrates and products 
% and get MZs
[shortestPathTable_unique, testidx] = unique(shortestPathTable);
shortestPathTable_unique.Path_length = shortestPathTable_length(testidx);
shortestPathTable_unique_sub = cell(size(shortestPathTable_unique,1),1);
shortestPathTable_unique_prod = cell(size(shortestPathTable_unique,1),1);
shortestPathTable_unique_sub_mz = zeros(size(shortestPathTable_unique,1),1);
shortestPathTable_unique_prod_mz = zeros(size(shortestPathTable_unique,1),1);
shortestPathTable_unique_sub_idx = zeros(size(shortestPathTable_unique,1),1);
shortestPathTable_unique_prod_idx = zeros(size(shortestPathTable_unique,1),1);
for i=1:size(shortestPathTable_unique,1)
    if (shortestPathTable_unique.Dir(i)==1)
        cursub = shortestPathTable_unique.Substrate{i};
        curprod = shortestPathTable_unique.Product{i};
    else
        cursub = shortestPathTable_unique.Product{i};
        curprod = shortestPathTable_unique.Substrate{i};
    end
    %substrate mz
    curidx = kegg_ids_unique_index(ismember(kegg_ids_unique, cursub));
    cursub_mz = metaboliteFilters.MZ(curidx);
    shortestPathTable_unique_sub_idx(i) = curidx;
    % product mz
    curidx = kegg_ids_unique_index(ismember(kegg_ids_unique, curprod));
    curprod_mz = metaboliteFilters.MZ(curidx);
    shortestPathTable_unique_prod_idx(i) = curidx;
    
    shortestPathTable_unique_sub{i} = cursub;
    shortestPathTable_unique_prod{i} = curprod;
    shortestPathTable_unique_sub_mz(i) = cursub_mz;
    shortestPathTable_unique_prod_mz(i) = curprod_mz;
end
shortestPathTable_unique.Substrate = shortestPathTable_unique_sub;
shortestPathTable_unique.Product = shortestPathTable_unique_prod;
shortestPathTable_unique.Substrate_MZ = shortestPathTable_unique_sub_mz;
shortestPathTable_unique.Product_MZ = shortestPathTable_unique_prod_mz;
shortestPathTable_unique.Substrate_IDX = shortestPathTable_unique_sub_idx;
shortestPathTable_unique.Product_IDX = shortestPathTable_unique_prod_idx;
% write to table
writetable(shortestPathTable_unique, [outputFolder, 'table_shortest_path_subprod_sorted_unique.csv']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write correlation tables to file
% corr RNA
corrtable = array2table(kegg_sub_prod_products_bestcorr_pos,...
    'RowNames', kegg_sub_prod_products_unique,...
    'VariableNames', species_list);
writetable(corrtable, ...
    [outputFolder, 'table_kegg_RNA_sub_prod_products_bestcorr_pos.csv'],...
    'WriteRowNames',true);
% pvalues
corrtable = array2table(kegg_sub_prod_products_bestcorrP_pos,...
    'RowNames', kegg_sub_prod_products_unique,...
    'VariableNames', species_list);
writetable(corrtable, ...
    [outputFolder, 'table_kegg_RNA_sub_prod_products_bestcorrP_pos.csv'],...
    'WriteRowNames',true);
% enzymes
corrtable = cell(size(kegg_sub_prod_products_bestcorr_geneidx_pos));
corrtableEC = cell(size(kegg_sub_prod_products_bestcorr_geneidx_pos));
for i=1:size(kegg_sub_prod_products_bestcorr_geneidx_pos,1)
    for j=1:size(kegg_sub_prod_products_bestcorr_geneidx_pos,2)
        if ~isempty(kegg_sub_prod_products_bestcorr_geneidx_pos{i,j})
            corrtable{i,j} = strjoin(geneTable.species_genes(kegg_sub_prod_products_bestcorr_geneidx_pos{i,j}),"|");
            curec = geneTable.EC(kegg_sub_prod_products_bestcorr_geneidx_pos{i,j});
            curec = unique(curec);
            curec(ismember(curec, {'nan'}))=[];
            corrtableEC{i,j} = strjoin(curec, "|");
        end
    end
end
corrtable = cell2table(corrtable,...
    'RowNames', kegg_sub_prod_products_unique,...
    'VariableNames', species_list);
writetable(corrtable, ...
    [outputFolder, 'table_kegg_RNA_sub_prod_products_bestcorr_genes_pos.csv'],...
    'WriteRowNames',true);
corrtableEC = cell2table(corrtableEC,...
    'RowNames', kegg_sub_prod_products_unique,...
    'VariableNames', species_list);
writetable(corrtableEC, ...
    [outputFolder, 'table_kegg_RNA_sub_prod_products_bestcorr_EC_pos.csv'],...
    'WriteRowNames',true);

% corr DNA
corrtable = array2table(kegg_DNA_sub_prod_products_bestcorr_pos,...
    'RowNames', kegg_sub_prod_products_unique,...
    'VariableNames', species_list);
writetable(corrtable, ...
    [outputFolder, 'table_kegg_DNA_sub_prod_products_bestcorr_pos.csv'],...
    'WriteRowNames',true);
% pvalues
corrtable = array2table(kegg_DNA_sub_prod_products_bestcorrP_pos,...
    'RowNames', kegg_sub_prod_products_unique,...
    'VariableNames', species_list);
writetable(corrtable, ...
    [outputFolder, 'table_kegg_DNA_sub_prod_products_bestcorrP_pos.csv'],...
    'WriteRowNames',true);
% enzymes
corrtable = cell(size(kegg_DNA_sub_prod_products_bestcorr_geneidx_pos));
corrtableEC = cell(size(kegg_DNA_sub_prod_products_bestcorr_geneidx_pos));
for i=1:size(kegg_DNA_sub_prod_products_bestcorr_geneidx_pos,1)
    for j=1:size(kegg_DNA_sub_prod_products_bestcorr_geneidx_pos,2)
        if ~isempty(kegg_DNA_sub_prod_products_bestcorr_geneidx_pos{i,j})
            corrtable{i,j} = strjoin(geneTable.species_genes(kegg_DNA_sub_prod_products_bestcorr_geneidx_pos{i,j}),"|");
            curec = geneTable.EC(kegg_DNA_sub_prod_products_bestcorr_geneidx_pos{i,j});
            curec = unique(curec);
            curec(ismember(curec, {'nan'}))=[];
            corrtableEC{i,j} = strjoin(curec, "|");
        end
    end
end
corrtable = cell2table(corrtable,...
    'RowNames', kegg_sub_prod_products_unique,...
    'VariableNames', species_list);
writetable(corrtable, ...
    [outputFolder, 'table_kegg_DNA_sub_prod_products_bestcorr_genes_pos.csv'],...
    'WriteRowNames',true);
corrtableEC = cell2table(corrtableEC,...
    'RowNames', kegg_sub_prod_products_unique,...
    'VariableNames', species_list);
writetable(corrtableEC, ...
    [outputFolder, 'table_kegg_DNA_sub_prod_products_bestcorr_EC_pos.csv'],...
    'WriteRowNames',true);


% corr OTU
corrtable = array2table(speciesOTU_products_corr,...
    'RowNames', kegg_sub_prod_products_unique,...
    'VariableNames', species_list);
writetable(corrtable, ...
    [outputFolder, 'table_speciesOTU_products_corr.csv'],...
    'WriteRowNames',true);
% pvalues
corrtable = array2table(speciesOTU_products_corrP,...
    'RowNames', kegg_sub_prod_products_unique,...
    'VariableNames', species_list);
writetable(corrtable, ...
    [outputFolder, 'table_speciesOTU_products_corrP.csv'],...
    'WriteRowNames',true);

