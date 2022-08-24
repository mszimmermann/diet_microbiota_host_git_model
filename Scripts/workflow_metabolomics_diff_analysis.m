%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% differential analysis of metabolite intensitie sbetween the mouse gourps

% call script defining file dependenciesand global variables
addpath(genpath('.\'));
add_global_and_file_dependencies

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data requirements:
% 'metabolites_allions_combined_formulas_with_metabolite_filters.csv'
% 'metabolites_allions_combined_norm_intensity.csv'
% 'hmdbPTWtables.mat'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output: 
% Files:
% 'table_diff_abundance_metabolite_ions_removed2outliers.csv'
% or 'table_diff_abundance_metabolite_ions.csv'
% Figures:
% 'fig_sup_volcanos_combined_ann_metabolites_removed2outliers.pdf'
% or 'fig_sup_volcanos_combined_ann_metabolites.pdf'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read metabolite data
metaboliteData = readtable([inputFolder,...
    'metabolites_allions_combined_norm_intensity.csv']);
    
metaboliteFilters = readtable([inputFolder,...
    'metabolites_allions_combined_formulas_with_metabolite_filters.csv']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% separate data into a matrix of normalized intensities
% and get information on diet, tissue and mouse group
dataColumns = metaboliteData.Properties.VariableNames;
% get columns that have mouse suffix
dataColumns = dataColumns(cellfun(@(x) contains(x,'_M'), dataColumns));
combinedTissues = cell(size(dataColumns));
combinedType = cell(size(dataColumns));
combinedDiet = cell(size(dataColumns));
combinedMouse = cell(size(dataColumns));
for i=1:length(dataColumns)
    curcol = strsplit(dataColumns{i}, '_');
    combinedTissues{i} = curcol{3};
    combinedType{i} = curcol{2};
    combinedDiet{i} = curcol{1};
    combinedMouse{i} = curcol{4};
end
sampleType_unique = unique(combinedType);
sampleDiet_unique = unique(combinedDiet);
sampleTissue_unique = unique(combinedTissues);

sampleTissue_order = {'SI1','SI2','SI3','Cecum','Colon','Feces','Liver','Serum'};
% order tissues along the GIT
if length(intersect(sampleTissue_unique, sampleTissue_order))==length(sampleTissue_unique)
    sampleTissue_unique = sampleTissue_order;
else
    sprintf('Non-overlapping tissue: %s', strjoin(setdiff( sampleTissue_unique, sampleTissue_order),'; '));
end

% convert table to matrix
combinedIntensitiesNorm = table2array(metaboliteData(:,dataColumns));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% flag whether remove two outliers
remove_outliers_flag = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform differential analysis of the ions across tissues                                              
diffMatrix_all = zeros(size(combinedIntensitiesNorm,1),...
                length(sampleTissue_unique)*16);
diffMatrix_columns = cell(1,length(sampleTissue_unique)*16);           
idx=1;
for d = 1:length(sampleTissue_unique)
    % Differential analysis of ions to identify drug metabolites
    curtissueidx = ismember(combinedTissues, sampleTissue_unique{d});
    tic
    if remove_outliers_flag
        curdata = workflow_ReplaceTwoOutliers(...
                                                        combinedIntensitiesNorm(:,curtissueidx),...
                                                        combinedDiet(curtissueidx),...
                                                        combinedType(curtissueidx));
    else
        curdata = combinedIntensitiesNorm(:,curtissueidx);
    end
                                                    
    [fcMatrix_HFD_DCGF, pMatrix_HFD_DCGF,...
      fcMatrix_CTR_DCGF, pMatrix_CTR_DCGF,...
      fcMatrix_DC_HFDCTR, pMatrix_DC_HFDCTR,...
      fcMatrix_GF_HFDCTR, pMatrix_GF_HFDCTR] = workflow_DifferentialAnalysisHFD(...
                                                        curdata,...
                                                        combinedDiet(curtissueidx),...
                                                        combinedType(curtissueidx));
    % perform FDR adjustment for all ions
    pFDRMatrix_HFD_DCGF = mafdr(pMatrix_HFD_DCGF,'bhfdr',1);
    pFDRMatrix_CTR_DCGF = mafdr(pMatrix_CTR_DCGF,'bhfdr',1);
    pFDRMatrix_DC_HFDCTR = mafdr(pMatrix_DC_HFDCTR,'bhfdr',1);
    pFDRMatrix_GF_HFDCTR = mafdr(pMatrix_GF_HFDCTR, 'bhfdr',1);
    % perform FDR adjustment for filtered ions
    pFDRMatrix_HFD_DCGF_filtered = ones(size(pMatrix_HFD_DCGF));
    pFDRMatrix_CTR_DCGF_filtered = ones(size(pMatrix_CTR_DCGF));
    pFDRMatrix_DC_HFDCTR_filtered = ones(size(pMatrix_DC_HFDCTR));
    pFDRMatrix_GF_HFDCTR_filtered = ones(size(pMatrix_GF_HFDCTR));
    % adjust only filtered
    pFDRMatrix_HFD_DCGF_filtered(metaboliteFilters.MetaboliteFilter==1) =...
        mafdr(pMatrix_HFD_DCGF(metaboliteFilters.MetaboliteFilter==1),'bhfdr',1);
    pFDRMatrix_CTR_DCGF_filtered(metaboliteFilters.MetaboliteFilter==1) = ...
        mafdr(pMatrix_CTR_DCGF(metaboliteFilters.MetaboliteFilter==1),'bhfdr',1);
    pFDRMatrix_DC_HFDCTR_filtered(metaboliteFilters.MetaboliteFilter==1) = ...
        mafdr(pMatrix_DC_HFDCTR(metaboliteFilters.MetaboliteFilter==1),'bhfdr',1);
    pFDRMatrix_GF_HFDCTR_filtered(metaboliteFilters.MetaboliteFilter==1) = ...
        mafdr(pMatrix_GF_HFDCTR(metaboliteFilters.MetaboliteFilter==1), 'bhfdr',1);
    toc

    diffMatrix_all(:,idx) = fcMatrix_HFD_DCGF;
    diffMatrix_columns{idx} = [sampleTissue_unique{d}, '_fc_HFD_DCGF'];
    diffMatrix_all(:,idx+1) = pMatrix_HFD_DCGF;
    diffMatrix_columns{idx+1} = [sampleTissue_unique{d}, '_p_HFD_DCGF'];
    diffMatrix_all(:,idx+2) = pFDRMatrix_HFD_DCGF;
    diffMatrix_columns{idx+2} = [sampleTissue_unique{d}, '_pFDR_HFD_DCGF'];
    diffMatrix_all(:,idx+3) = pFDRMatrix_HFD_DCGF_filtered;
    diffMatrix_columns{idx+3} = [sampleTissue_unique{d}, '_pFDRf_HFD_DCGF'];
    diffMatrix_all(:,idx+4) = fcMatrix_CTR_DCGF;
    diffMatrix_columns{idx+4} = [sampleTissue_unique{d}, '_fc_CTR_DCGF'];
    diffMatrix_all(:,idx+5) = pMatrix_CTR_DCGF;
    diffMatrix_columns{idx+5} = [sampleTissue_unique{d}, '_p_CTR_DCGF'];
    diffMatrix_all(:,idx+6) = pFDRMatrix_CTR_DCGF;
    diffMatrix_columns{idx+6} = [sampleTissue_unique{d}, '_pFDR_CTR_DCGF'];
    diffMatrix_all(:,idx+7) = pFDRMatrix_CTR_DCGF_filtered;
    diffMatrix_columns{idx+7} = [sampleTissue_unique{d}, '_pFDRf_CTR_DCGF'];
    diffMatrix_all(:,idx+8) = fcMatrix_DC_HFDCTR;
    diffMatrix_columns{idx+8} = [sampleTissue_unique{d}, '_fc_DC_HFDCTR'];
    diffMatrix_all(:,idx+9) = pMatrix_DC_HFDCTR;
    diffMatrix_columns{idx+9} = [sampleTissue_unique{d}, '_p_DC_HFDCTR'];
    diffMatrix_all(:,idx+10) = pFDRMatrix_DC_HFDCTR;
    diffMatrix_columns{idx+10} = [sampleTissue_unique{d}, '_pFDR_DC_HFDCTR'];
    diffMatrix_all(:,idx+11) = pFDRMatrix_DC_HFDCTR_filtered;
    diffMatrix_columns{idx+11} = [sampleTissue_unique{d}, '_pFDRf_DC_HFDCTR'];
    diffMatrix_all(:,idx+12) = fcMatrix_GF_HFDCTR;
    diffMatrix_columns{idx+12} = [sampleTissue_unique{d}, '_fc_GF_HFDCTR'];
    diffMatrix_all(:,idx+13) = pMatrix_GF_HFDCTR;
    diffMatrix_columns{idx+13} = [sampleTissue_unique{d}, '_p_GF_HFDCTR'];
    diffMatrix_all(:,idx+14) = pFDRMatrix_GF_HFDCTR;
    diffMatrix_columns{idx+14} = [sampleTissue_unique{d}, '_pFDR_GF_HFDCTR'];
    diffMatrix_all(:,idx+15) = pFDRMatrix_GF_HFDCTR_filtered;
    diffMatrix_columns{idx+15} = [sampleTissue_unique{d}, '_pFDRf_GF_HFDCTR'];
    
    idx = idx+16;
    fprintf('Calculated fold changes for %s\n', sampleTissue_unique{d}); 
end

clear curtissueidx fcMatrix_HFD_DCGF pMatrix_HFD_DCGF fcMatrix_CTR_DCGF pMatrix_CTR_DCGF
clear fcMatrix_DC_HFDCTR pMatrix_DC_HFDCTR fcMatrix_GF_HFDCTR pMatrix_GF_HFDCTR


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot volcano plots for selected tissues
fig = figure('units','normalized','outerposition',[0 0 1 1]);
fcThreshold = log2(2);
fdrThreshold = 0.05;
dietcomparison = {'HFD, DC vs GF',...
                  'CTR, DC vs GF',...
                  'DC, HFD vs CTR',...
                  'GF, HFD vs CTR'};
for comparisontype = 1:length(dietcomparison)
    for i=1:length(sampleTissue_unique)
        
        switch comparisontype
            case 1
                curFC = log2(diffMatrix_all(:,ismember(diffMatrix_columns,...
                    {[sampleTissue_unique{i},'_fc_HFD_DCGF']})));
                curP = diffMatrix_all(:,ismember(diffMatrix_columns,...
                    {[sampleTissue_unique{i},'_pFDR_HFD_DCGF']}));
            case 2
                curFC = log2(diffMatrix_all(:,ismember(diffMatrix_columns,...
                    {[sampleTissue_unique{i},'_fc_CTR_DCGF']})));
                curP = diffMatrix_all(:,ismember(diffMatrix_columns,...
                    {[sampleTissue_unique{i},'_pFDR_CTR_DCGF']}));
            case 3
                 curFC = log2(diffMatrix_all(:,ismember(diffMatrix_columns,...
                    {[sampleTissue_unique{i},'_fc_DC_HFDCTR']})));
                 curP = diffMatrix_all(:,ismember(diffMatrix_columns,...
                    {[sampleTissue_unique{i},'_pFDR_DC_HFDCTR']}));
            case 4
                curFC = log2(diffMatrix_all(:,ismember(diffMatrix_columns,...
                    {[sampleTissue_unique{i},'_fc_GF_HFDCTR']})));
                curP = diffMatrix_all(:,ismember(diffMatrix_columns,...
                    {[sampleTissue_unique{i},'_pFDR_GF_HFDCTR']}));
        end
        
        % leave only annotated metabolites
        curFC = curFC(metaboliteFilters.MetaboliteFilter==1);
        curP = curP(metaboliteFilters.MetaboliteFilter==1);
               
        sp = subplot(4,length(sampleTissue_unique),...
            (comparisontype-1)*length(sampleTissue_unique)+i);
        plotVolcano(curFC, curP, [], fcThreshold, fdrThreshold, [],sp);
        if i==1
            ylabel(dietcomparison{comparisontype})
        end

        title(sampleTissue_unique{i})
        xlim([-15 15])
        ylim([0 11])
        axis square
    end
end
suptitle('Annotated filtered metabolites')
orient landscape
if remove_outliers_flag
    print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
        [figureFolder ...
        'fig_sup_volcanos_combined_ann_metabolites_removed2outliers.pdf'])
else
    print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
        [figureFolder...
        'fig_sup_volcanos_combined_ann_metabolites.pdf'])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save fold change and p-value table to file
diffMatrix_all_table = [metaboliteData(:,1:5),...
    array2table(diffMatrix_all,'VariableNames', diffMatrix_columns)];
if remove_outliers_flag
    writetable(diffMatrix_all_table, ...
        [outputFolder ,...
        'table_diff_abundance_metabolite_ions_removed2outliers.csv']);
else
    writetable(diffMatrix_all_table, ...
        [outputFolder ,...
        'table_diff_abundance_metabolite_ions.csv']);
end


% diffMatrix_all_table = readtable([resultsFolder, ...
%     'table_diff_abundance_metabolite_ions.csv']);

