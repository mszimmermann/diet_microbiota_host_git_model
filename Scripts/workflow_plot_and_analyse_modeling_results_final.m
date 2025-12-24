% Plot and analyse modeling results
% call script defining file dependenciesand global variables
addpath(genpath('.\'));
add_global_and_file_dependencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data requirements: 
% 'metabolites_allions_combined_formulas_with_metabolite_filters_spatial100clusters_with_mean.csv'
% 'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_allions.csv'
% Matlab file with manually selected clustergrams
% 'cgo_clustergrams_of_model_coefficients.mat'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output:
% Files:
% 'table_hierarchical_clustering_groups.csv'
% Figures:
% 'fig4b_clustergram_2LIhos1LIbact_model_coefs_annotated_ions_in_GIT_Rmodel_corr_0_7'
% 'fig4b_barh_2LIhos1LIbact_model_cluster_sizes_percent.pdf'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% try to cluster ions based on their GI profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read annotation from file
%annotationTableSpatialClusters = readtable([inputFolder ...
%    'metabolites_allions_combined_formulas_with_metabolite_filters_spatial100clusters_with_mean.csv']);
annotationTableSpatialClusters = readtable([outputFolder ...
    'metabolites_allions_combined_formulas_with_metabolite_filters_spatial100clusters_with_mean_with_CVR_0925.csv']);
% get clustering info
% select columns from annotationTable with cluster info
annColumns = annotationTableSpatialClusters.Properties.VariableNames;
clusterColumns = annColumns(cellfun(@(x) contains(x,'spatial_clust100'), annColumns));
[~, clusteridx] = intersect(annColumns, clusterColumns, 'stable');
spatialClusters = table2array(annotationTableSpatialClusters(:, clusteridx));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modelingResults = readtable([resultsFolder ...
%             'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_allions.csv']);
% modelingResults = readtable([outputFolder ...
%             'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_allions_with_CVR.csv']);
% x_met_smooth = modelingResults{:, width(modelingResults)-8:end};
% % get correlations calculated with reverse problem
% if isnumeric(modelingResults.ReciprocalCorr(1))
%     x_data_corr = modelingResults.ReciprocalCorr;
% else
%     % it is not numeric, probably contains NaN - convert to numeric
%     x_data_corr = cellfun(@(x) str2double(x), modelingResults.ReciprocalCorr);
% end
% coefvalues = modelingResults.Properties.VariableNames(width(modelingResults)-8:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze either DC or CVR data
%analyzeGroup = 'CVR'; 
analyzeGroup = 'DC';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = [outputFolder ...
            'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_best'...
            analyzeGroup];
corrthreshold = 0.7;
% try combining solutions
sel_crit1 = 'IP';
sel_crit2 = 'LI_PCC_within_high_total';

[met_info_combined, met_bestsol_combined] = ...
                combine_bestsols_from_file(filename, ...
                sel_crit1, sel_crit2,...
                corrthreshold);
% get model coefficients
x_met_smooth = met_bestsol_combined.x;
% get correlations calculated with reverse problem
x_data_corr = met_bestsol_combined.x_sel_CorrRev;
coefvalues = met_bestsol_combined.coefvalues;

% get all solutions from file
filename = [outputFolder ...
            'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_all'...
            analyzeGroup];
[met_info_read, met_gitfits_read] = read_allsols_from_files(filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot correlation of restored and original data
% calculate differentce in corr distrbutions
selected_mets = find((sum(spatialClusters,2)>0) &...
                     (annotationTableSpatialClusters.MetaboliteFilter==1));

filename = [figureFolder,...
    'Fig4a_histogram_MAXcorr_model_2LIcoefHost1LIcoefbact_'...
    analyzeGroup, 'gitann_all_', sel_crit1, '_', sel_crit2];
plot_gitfit_model_corr(met_gitfits_read(selected_mets), filename)

% plot correlation of restored and original data for the best solution
% calculate differentce in corr distrbutions
filename = [figureFolder,...
    'Fig4a_histogram_MAXcorr_model_2LIcoefHost1LIcoefbact_'...
    analyzeGroup, 'gitann_best_combined_', sel_crit1, '_', sel_crit2];
plot_bestfit_model_corr(filter_bestsols_by_index(met_bestsol_combined, selected_mets),...
                        met_gitfits_read(selected_mets), filename)

% save combined best solution to file
outputfilename = [outputFolder ...
            'table_model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_'...
            analyzeGroup, '_combined_',...
            sel_crit1, '_', sel_crit2];
print_combined_bestsol_to_files(met_info_combined, met_bestsol_combined,...
                                outputfilename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cluster only metabolites with corr>0.7
% uncomment if using whole metabolite table
if size(x_met_smooth,1) == size(annotationTableSpatialClusters,1)
    clustidx = find(annotationTableSpatialClusters.MetaboliteFilter==1);
else
    clustidx = 1:size(x_met_smooth,1); % for testing purposed plot all 
end
clustdata = x_met_smooth(clustidx,2:end)';
clustrows = coefvalues(2:end);
% filter by resiprocal corr
clustcorr = x_data_corr(clustidx);
clustdata = clustdata(:, clustcorr>=0.7);
clustidx = clustidx(clustcorr>=0.7);

for i=1:size(clustdata,2)
    clustdata(:,i) = clustdata(:,i)/max(abs(clustdata(:,i)));
end
clustnan = isnan(sum(clustdata,1));
clustdata(:,clustnan) = [];
clustidx(clustnan)=[];

clustdist ='cityblock';%'correlation';%'euclidean';%
cgo = clustergram(clustdata,...
            'RowLabels', clustrows,...
            'ColumnLabels', clustidx,...
            'ColumnPdist',clustdist,...
            'RowPdist', clustdist,...
            'DisplayRange', 1,...
            'colormap', redbluecmap);
% solution to turn on colormap programmatically from https://stackoverflow.com/questions/20648627/turn-on-colorbar-programmatically-in-clustergram
% combined with solution from https://de.mathworks.com/matlabcentral/answers/305274-reduce-font-size-of-column-labels-in-clustergram?#answer_533054
%cgfig = findall(0,'type','figure', 'tag', 'Clustergram'); % Figure handle
%cgax = findobj(cgfig, 'type','axes','tag','HeatMapAxes'); % main axis handle

cbButton = findall(0,'tag','HMInsertColorbar');
ccb = get(cbButton,'ClickedCallback');
set(cbButton,'State','on')
ccb{1}(cbButton,[],ccb{2})
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Manually add colorbar to clustergram and print to figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = cgo.plot;
orient landscape
% add title to colobar
C = findall(gcf,'type','ColorBar');                         
C.Label.String = ['Relative ' analyzeGroup ' model coefficient value'];

print(gcf, '-vector', '-dpdf', '-r600', '-bestfit',...
    [figureFolder,...
    'fig4b_clustergram_2LIhos1LIbact_model_coefs_annotated_ions_in_GIT_Rmodel_combined_'...
    analyzeGroup, '_', sel_crit1, '_', sel_crit2,'corr_0_7'])


      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clustergrams were manually selected from the clustergram window
% they can be loaded from 'cgo_clustergrams_of_model_coefficients.mat'
%load([resultsFolder 'cgo_clustergrams_of_model_coefficients.mat'])
clustergramFile = ['cgo_clustergrams_model_'...
                    analyzeGroup, '_combined_IP_LIPPCwhighTotalPCC.mat'];
load([outputFolder clustergramFile]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print cluster size info        
if isequal(analyzeGroup, 'DC')
    % order groups according to the full clustergram
    cluster_order = [{cellfun(@(x) (str2double(x)), cgoHColhfd.ColumnLabels)},...
            {setdiff(cellfun(@(x) (str2double(x)), cgoHighHostColon.ColumnLabels),...
                     cellfun(@(x) (str2double(x)), cgoHColhfd.ColumnLabels))},...
            {cellfun(@(x) (str2double(x)),cgoHCIhigh_BLIctrneg.ColumnLabels)},...
            {cellfun(@(x) (str2double(x)),cgoHSI_only.ColumnLabels)},...
            {setdiff(cellfun(@(x) (str2double(x)), cgoHSI_and_mix.ColumnLabels),...
                     union(cellfun(@(x) (str2double(x)), cgoHSI_only.ColumnLabels),...
                           cellfun(@(x) (str2double(x)), cgoHCIhigh_BLIctrneg.ColumnLabels)))},...
            {cellfun(@(x) (str2double(x)),cgoHLIhfd.ColumnLabels)},...
            {cellfun(@(x) (str2double(x)),cgoHLIctr.ColumnLabels)},...
            {cellfun(@(x) (str2double(x)),cgoBLIhfd.ColumnLabels)},...
            {cellfun(@(x) (str2double(x)),cgoBLIctr.ColumnLabels)},...
            {setdiff(cellfun(@(x) (str2double(x)), cgoLIhigh_all.ColumnLabels),...
                     union(cellfun(@(x) (str2double(x)), cgoHLIhfd.ColumnLabels),...
                           union(cellfun(@(x) (str2double(x)), cgoHLIctr.ColumnLabels),...
                                 union(cellfun(@(x) (str2double(x)), cgoBLIhfd.ColumnLabels),...
                                       cellfun(@(x) (str2double(x)), cgoBLIctr.ColumnLabels)))))},...
            {cellfun(@(x) (str2double(x)),cgoHColhigh_mixed.ColumnLabels)},...
            {cellfun(@(x) (str2double(x)),cgoHSIhighHLIlow.ColumnLabels)},...
            {cellfun(@(x) (str2double(x)),cgoBLIctrneg_mixed.ColumnLabels)},...
            {cellfun(@(x) (str2double(x)),cgoHLIhfdneg_mixed.ColumnLabels)},...
            {cellfun(@(x) (str2double(x)),cgoHLIctrneg_mixed.ColumnLabels)}];
    
    % get the indeces missing from subclusters
    % join all subcluster idx
    cluster_subcluster_idx = cat(2,cluster_order{:});
    % find missing ones and add them manually
    [~, idx] = setdiff(cellfun(@(x) str2double(x), cgo.ColumnLabels), cluster_subcluster_idx, 'stable');
    extra_subcluster1 = cellfun(@(x) (str2double(x)), cgo.ColumnLabels(idx(1:3)));
    extra_subcluster2 = cellfun(@(x) (str2double(x)), cgo.ColumnLabels(idx(4:25)));
    extra_subcluster3 = cellfun(@(x) (str2double(x)), cgo.ColumnLabels(idx(26:end)));
    cluster_order = [cluster_order(1:5), {extra_subcluster1}, ...
                     cluster_order(6:14), {extra_subcluster2}, ...
                     cluster_order(15:end), {extra_subcluster3}];
    
    cluster_sizes = [length(cgoHColhfd.ColumnLabels),...
            length(cgoHighHostColon.ColumnLabels)-length(cgoHColhfd.ColumnLabels),...
            length(cgoHCIhigh_BLIctrneg.ColumnLabels),...
            length(cgoHSI_only.ColumnLabels),...
            length(cgoHSI_and_mix.ColumnLabels) - length(cgoHSI_only.ColumnLabels) - length(cgoHCIhigh_BLIctrneg.ColumnLabels),...
            length(extra_subcluster1),...
            length(cgoHLIhfd.ColumnLabels),...
            length(cgoHLIctr.ColumnLabels),...
            length(cgoBLIhfd.ColumnLabels),...
            length(cgoBLIctr.ColumnLabels),...
            length(cgoLIhigh_all.ColumnLabels) - (length(cgoHLIhfd.ColumnLabels)+...
                                                  length(cgoHLIctr.ColumnLabels)+...
                                                  length(cgoBLIhfd.ColumnLabels)+...
                                                  length(cgoBLIctr.ColumnLabels)),...
            length(cgoHColhigh_mixed.ColumnLabels),...
            length(cgoHSIhighHLIlow.ColumnLabels),...
            length(cgoBLIctrneg_mixed.ColumnLabels),...
            length(cgoHLIhfdneg_mixed.ColumnLabels),...
            length(extra_subcluster2),...
            length(cgoHLIctrneg_mixed.ColumnLabels),...            
            length(extra_subcluster3),...
            ];
    
    cluster_attribution = zeros(size(annotationTableSpatialClusters,1),1);
    for i=1:length(cluster_order)
        cluster_attribution(cluster_order{i})=i;
    end
    % add all other metabolites to a cluster size plus one group
    cluster_attribution(setdiff(clustidx, find(cluster_attribution>0))) = length(cluster_order)+1; 
    hier_clust_attribution = table(annotationTableSpatialClusters.MZ,...
                                    annotationTableSpatialClusters.RT,... 
                                    annotationTableSpatialClusters.MetaboliteFilter,...
                                    cluster_attribution,...
                                    'VariableNames', {'MZ', 'RT', 'MetaboliteFilter', 'HierarchicalClustGroup'});
    writetable(hier_clust_attribution, [outputFolder, ...
        'table_met_' analyzeGroup '_hierarchical_clustering_groups.csv']);
else
    cluster_order = [{cellfun(@(x) (str2double(x)), cgoHSIonly.ColumnLabels)},...
            {setdiff(cellfun(@(x) (str2double(x)), cgoHSImixed_all.ColumnLabels),...
                     cellfun(@(x) (str2double(x)), cgoHSIonly.ColumnLabels))},...
            {cellfun(@(x) (str2double(x)),cgoHColhigh.ColumnLabels)},...
            {cellfun(@(x) (str2double(x)),cgoBLIctr.ColumnLabels)},...
            {cellfun(@(x) (str2double(x)),cgoBLIhfd.ColumnLabels)},...
            {setdiff(cellfun(@(x) (str2double(x)), cgoBLIall.ColumnLabels),...
                     union(cellfun(@(x) (str2double(x)), cgoBLIhfd.ColumnLabels),...
                           cellfun(@(x) (str2double(x)), cgoBLIctr.ColumnLabels)))},...
            {cellfun(@(x) (str2double(x)),cgoHLIhfd_low.ColumnLabels)},...
            {cellfun(@(x) (str2double(x)),cgoHLIctr_low.ColumnLabels)},...
            {setdiff(cellfun(@(x) (str2double(x)), cgoHLIlow_all.ColumnLabels),...
                     union(cellfun(@(x) (str2double(x)), cgoHLIhfd_low.ColumnLabels),...
                           cellfun(@(x) (str2double(x)), cgoHLIctr_low.ColumnLabels)))},...
            {cellfun(@(x) (str2double(x)),cgoHLIctr_high.ColumnLabels)},...
            {cellfun(@(x) (str2double(x)),cgoHLIhfd_high.ColumnLabels)},...
            {cellfun(@(x) (str2double(x)),cgoBLI_low_mixed.ColumnLabels)}];
   
    % get the indeces missing from subclusters
    % join all subcluster idx
    cluster_subcluster_idx = cat(2,cluster_order{:});
    % find missing ones and add them manually
    [~, idx] = setdiff(cellfun(@(x) str2double(x), cgo.ColumnLabels), cluster_subcluster_idx, 'stable');
    extra_subcluster1 = cellfun(@(x) (str2double(x)), cgo.ColumnLabels(idx(1:8)));
    extra_subcluster2 = cellfun(@(x) (str2double(x)), cgo.ColumnLabels(idx(9:end)));
    cluster_order = [cluster_order(1:12), {extra_subcluster1}, ...
                     cluster_order(13:end), {extra_subcluster2}];


    % calculate cluster sizes
    cluster_sizes = [length(cgoHSIonly.ColumnLabels),...
            length(cgoHSImixed_all.ColumnLabels)-length(cgoHSIonly.ColumnLabels),...
            length(cgoHColhigh.ColumnLabels),...
            length(cgoBLIctr.ColumnLabels),...
            length(cgoBLIhfd.ColumnLabels),...
            length(cgoBLIall.ColumnLabels) - length(cgoBLIhfd.ColumnLabels) - length(cgoBLIctr.ColumnLabels),...
            length(cgoHLIhfd_low.ColumnLabels),...
            length(cgoHLIctr_low.ColumnLabels),...
            length(cgoHLIlow_all.ColumnLabels)-length(cgoHLIhfd_low.ColumnLabels)-length(cgoHLIctr_low.ColumnLabels),...
            length(cgoHLIctr_high.ColumnLabels),...
            length(cgoHLIhfd_high.ColumnLabels),...
            length(extra_subcluster1),...
            length(cgoBLI_low_mixed.ColumnLabels),...
            length(extra_subcluster2)];
        
    cluster_attribution = zeros(size(annotationTableSpatialClusters,1),1);    
    for i=1:length(cluster_order)
        cluster_attribution(cluster_order{i})=i;
    end
    % add all other metabolites to a cluster size plus one group
    cluster_attribution(setdiff(clustidx, find(cluster_attribution>0))) = length(cluster_order)+1; 
    hier_clust_attribution = table(annotationTableSpatialClusters.MZ,...
                                    annotationTableSpatialClusters.RT,... 
                                    annotationTableSpatialClusters.MetaboliteFilter,...
                                    cluster_attribution,...
                                    'VariableNames', {'MZ', 'RT', 'MetaboliteFilter', 'HierarchicalClustGroup'});
    writetable(hier_clust_attribution, [outputFolder, ...
        'table_met_' analyzeGroup '_hierarchical_clustering_groups.csv']);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot cluster sizes as stacked bar plot
figure
barh([cluster_sizes;cluster_sizes], 'stacked')
xlim([0 sum(cluster_sizes)])
for i=1:length(cluster_sizes)
    text(sum(cluster_sizes(1:(i-1))), 1, num2str(cluster_sizes(i)))
    text(sum(cluster_sizes(1:(i-1))), 2, num2str(round(100*cluster_sizes(i)/length(cgo.ColumnLabels))))
end
set(gca, 'Ytick', 1:2)
set(gca, 'YtickLabel', {'Number', 'Percent'})
title(['Hierarchical clustering metabolite group sizes ' analyzeGroup])

orient landscape
print(gcf, '-vector', '-dpdf', '-r600', '-bestfit',...
    [figureFolder, ...
    ['fig4b_barh_2LIhos1LIbact_model_coefs_annotated_ions_in_GIT_Rmodel_combined_'...
    analyzeGroup, '_', sel_crit1, '_', sel_crit2,'corr_0_7.pdf']])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform pathway enrichment for clustering

load('hmdbPTWtables.mat');
for hmdbtype = 1:3
    switch hmdbtype
        case 1
            hmdbPtwIonTable = hmdbPTWclassTable;
            ptwNames = hmdbPTWclassNames;
        case 2
            hmdbPtwIonTable = hmdbPTWSubClassTable;
            ptwNames = hmdbPTWSubClassNames;
        case 3
            hmdbPtwIonTable = hmdbPTWSuperClassTable;
            ptwNames = hmdbPTWSuperClassNames;
    end
    
    % select only metabolites involved in clustering
    selected_mets = cellfun(@(x) (str2double(x)), cgo.ColumnLabels);
        
    % pathway enrichment analysis for clusters
    hclusterTables = cell(length(ptwNames), length(cluster_order)*3+1);
    ttidx = 1;
    keggCompoundsID = annotationTableSpatialClusters.CompoundID;
    keggCompoundsID = keggCompoundsID(selected_mets);
    ptwIonTable = hmdbPtwIonTable(:,selected_mets);

    for i = 1:length(cluster_order)
        tic
        labels = zeros(size(selected_mets));
        for jj=1:length(cluster_order{i})
            labels(selected_mets==cluster_order{i}(jj)) = 1;
        end
        
        [~, ~, enrichmentTable] = ptwEnrichmentMetabolitesIonsPTWTable(ptwIonTable,...
                                                keggCompoundsID, labels, ptwNames);

        %compile a joint table
        for jj=1:length(enrichmentTable)
             enrichmentTable{jj,8} = ['(',num2str(enrichmentTable{jj,4}), ',',num2str(enrichmentTable{jj,5}), ',',num2str(enrichmentTable{jj,6}), ',',num2str(enrichmentTable{jj,7}), ')' ];
        end
        enrichmentTable(:,4:7) = [];
        if ttidx == 1
            hclusterTables(:, ttidx) = enrichmentTable(:,3);
            ttidx = ttidx+1;
        end
        hclusterTables(:,ttidx:ttidx+2) = [enrichmentTable(:,1:2), enrichmentTable(:,4)];
        ttidx = ttidx+3;
    toc
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(hclusterTables{1,1})
        hclusterTables(:,1) = enrichmentTable(:,3);
    end
    switch hmdbtype
        case 1
            hclusterTables_hmdbClass = hclusterTables;
        case 2
            hclusterTables_hmdbSubClass = hclusterTables;
        case 3
            hclusterTables_hmdbSuperClass = hclusterTables;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot pathway enrichment results
% hclusterTables = hclusterTables_hmdbClass;
% hclusterTables = hclusterTables_hmdbSuperClass;
% hclusterTables = hclusterTables_hmdbSubClass;

hclusterTables = [hclusterTables_hmdbClass;...
                  hclusterTables_hmdbSuperClass;...
                  hclusterTables_hmdbSubClass];


showclust = 1;
startclust = 0;
ncluster = length(cluster_order);
numMet = zeros(size(hclusterTables,1), ncluster*showclust);
pMet = zeros(size(hclusterTables,1), ncluster*showclust); 
fdrMet = zeros(size(hclusterTables,1), ncluster*showclust);
idx = 1;
for i=1:ncluster
    numMet(:,idx) = cellfun(@(x) str2double(x(2:strfind(x,',')-1)),...
                                 hclusterTables(:,4+(i-1)*3));
    pMet(:,idx) = cell2mat(hclusterTables(:,2+(i-1)*3));
    fdrMet(:,idx) = cell2mat(hclusterTables(:,3+(i-1)*3));
    idx = idx+1;                          
end

% plot FDR values
selectpathways = sum(fdrMet<=0.1,2)>0;
displaymat = fdrMet(selectpathways,:);

% %or uncomment below to plot uncorrected p-values
%selectpathways = sum(pMet<0.01,2)>0;
%displaymat = pMet(selectpathways,:);

% display values above 0.1 as black
displaymat(displaymat>0.1)=1;

% Generate cluster names
rowLabels = cellfun(@(x) x{1}, hclusterTables(selectpathways,1), 'unif', 0);

pwenrcgo = clustergram(displaymat, ...
            'RowLabels', rowLabels,...'ColumnLabels', colLabels(1:size(displaymat,2)),...
            'ColumnPDist', 'euclidean',...
            'RowPDist', 'euclidean',...
            'Cluster', 'column',...
            'Symmetric', 0,...
            'DisplayRange', 1,...
            'Colormap', flipud(bone));
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot clustergram as imagesc with rows sorted according to clustering
% and columns proportional to cluster sizes
[~, ~, sortidx] = intersect(pwenrcgo.RowLabels, rowLabels, 'stable');
plotpwmat = displaymat(sortidx,:);
plotpwmat_rowLabels = rowLabels(sortidx);
plotpwmat_extended = zeros(size(plotpwmat,1), sum(cluster_sizes));
idx=1;
for i=1:length(cluster_sizes)
    plotpwmat_extended(:, idx:idx+cluster_sizes(i)-1) = repmat(plotpwmat(:,i),1,cluster_sizes(i));
    idx = idx+cluster_sizes(i);
end

fig = figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(plotpwmat_extended)
colormap(flipud(bone))
set(gca, 'YTick', 1:length(plotpwmat_rowLabels))
set(gca, 'YTickLabel', plotpwmat_rowLabels)
set(gca, 'XTick', cumsum(cluster_sizes))
set(gca, 'XTickLabel', cluster_sizes)
colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
orient landscape
print(gcf, '-vector', '-dpdf', '-r600', '-bestfit',...
            [figureFolder,...
            'fig4c_clustergram_' analyzeGroup,'_PTWENR_hcluster_combined_sol_hmdbCombinedClasses_2LIhos1LIbact_pFDR_le0_1.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot relative number of metabolites per group
displaymat = numMet(selectpathways,:);

displaymat = displaymat./repmat(cluster_sizes, size(displaymat,1),1);


%normalize by max in each category
for i=1:size(displaymat,1)
    displaymat(i,:) = displaymat(i,:)/max(displaymat(i,:));
end

%displaymat = pMet(selectpathways,selectclusters);
%displaymat = fdrMet(selectpathways,:);
displayvalues = '_nummet_';

%displaymat(displaymat>0.5)=1;
displaymat(displaymat>30)=30;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot clustergram as imagesc with rows sorted according to clustering
% and columns proportional to cluster sizes
[~, ~, sortidx] = intersect(pwenrcgo.RowLabels, rowLabels, 'stable');
plotpwmat = displaymat(sortidx,:);
plotpwmat_rowLabels = rowLabels(sortidx);
plotpwmat_extended = zeros(size(plotpwmat,1), sum(cluster_sizes));
% save midle of clusters for subsequent annotation of the plot
plotpwmat_extended_xcoord = zeros(1, length(cluster_sizes));
idx=1;
for i=1:length(cluster_sizes)
    plotpwmat_extended(:, idx:idx+cluster_sizes(i)-1) = repmat(plotpwmat(:,i),1,cluster_sizes(i));
    plotpwmat_extended_xcoord(i) = round(idx+cluster_sizes(i)/2);
    idx = idx+cluster_sizes(i);
end

fig = figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(plotpwmat_extended)
colormap(slanCM('greys')) % flipud(bone))
set(gca, 'YTick', 1:length(plotpwmat_rowLabels))
set(gca, 'YTickLabel', plotpwmat_rowLabels)
set(gca, 'XTick', cumsum(cluster_sizes))
set(gca, 'XTickLabel', cluster_sizes)
colorbar

%%%%%%%%%%%%%%%%%%%%%%
% Annotate imagesc with FDR values
%resort FDR values based on clustergram ordering
clusterFDR = fdrMet(selectpathways,:);
clusterFDR = clusterFDR(sortidx,:);
% round to two decimals
%clusterFDR = round(clusterFDR*100)/100;

for i=1:size(clusterFDR,1)
    for j=1:size(clusterFDR,2)
    
        curcolor = 'green';
        if clusterFDR(i,j)<=0.1
            curcolor = 'white';
        
            text(plotpwmat_extended_xcoord(j)-20,i+0.2, '*',...
                      'Color', curcolor, 'FontSize', 20);
        end
    end
end

% add title to colobar
C = findall(gcf,'type','ColorBar');                         
C.Label.String = 'Relative number of metabolites';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
orient landscape
print(gcf, '-vector', '-dpdf', '-r600', '-bestfit',...
            [figureFolder,...
            'fig4c_clustergram_' analyzeGroup, '_PTWENR_hcluster_combined_sol_',...
            displayvalues,'_hmdbCombinedClasses_2LIhos1LIbact_pFDR_le0_1.pdf'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save pathway enrichment to file
hclusterTables_colLabels = [ {'HMDBclass'},...
    strcat( arrayfun(@(x) ['Cluster ',num2str(x)], repelem(1:length(cluster_order),3), 'unif', 0),...
            {' '},...
            repmat([{'pvalue'},{'FDR'},{'stats'}],1,length(cluster_order)))];
        
for hmdbtype = 1:3
    switch hmdbtype
        case 1
            hclusterTables = hclusterTables_hmdbClass;
            curvarName = 'HMDBclass';
        case 2
            hclusterTables = hclusterTables_hmdbSubClass;
            curvarName = 'HMDBsubclass';
        case 3
            hclusterTables = hclusterTables_hmdbSuperClass;
            curvarName = 'HMDBsuperclass';
    end
    
    hclusterTables= cell2table(hclusterTables, 'VariableNames',...
                                hclusterTables_colLabels);
    writetable(hclusterTables,...
               [outputFolder...
               'ptwenr_', curvarName, '_GIT_Rmodel_combined_'...
                analyzeGroup, '_', sel_crit1, '_', sel_crit2,...
                'corr_0_7_hcluster_2LIhos1LIbact_hmdbv4.csv']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % calculate changing metabolites in the system
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % read metabolite data
% metaboliteDiff = readtable([resultsFolder,...
%     'table_diff_abundance_metabolite_ions_removed2outliers.csv']);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % get information from differential analysis on whether metabolites are 
% % differentially abundant in either serum or liver in at least
% % one diet between GF and colonized mice
% dietChangeNames = {'HFD_DCGF', 'CTR_DCGF'};
% fcThreshold = log2(1.5);
% pfdrThreshold = 0.05;
% pvalColumn = '_pFDR_'; % choose between _p_ and _pFDR_
% % record also liver and serum changes
% met_liver_changes = zeros(size(metaboliteDiff,1), length(dietChangeNames));
% met_serum_changes = zeros(size(metaboliteDiff,1), length(dietChangeNames));
% for diet_i=1:length(dietChangeNames)
%         % liver
%         % (record sign of change)
%         met_liver_changes(:,diet_i) = (abs(log2(metaboliteDiff{:,...
%              ismember(metaboliteDiff.Properties.VariableNames,...
%              {['Liver_fc_', dietChangeNames{diet_i}]})})) >= fcThreshold).*...
%              (metaboliteDiff{:,...
%              ismember(metaboliteDiff.Properties.VariableNames,...
%              {['Liver' pvalColumn, dietChangeNames{diet_i}]})} <= pfdrThreshold).*...
%              sign(log2(metaboliteDiff{:,...
%              ismember(metaboliteDiff.Properties.VariableNames,...
%              {['Liver_fc_', dietChangeNames{diet_i}]})}));
%         % serum
%         % (record sign of change)
%         met_serum_changes(:,diet_i) = (abs(log2(metaboliteDiff{:,...
%              ismember(metaboliteDiff.Properties.VariableNames,...
%              {['Serum_fc_', dietChangeNames{diet_i}]})})) >= fcThreshold).*...
%              (metaboliteDiff{:,...
%              ismember(metaboliteDiff.Properties.VariableNames,...
%              {['Serum' pvalColumn, dietChangeNames{diet_i}]})} <= pfdrThreshold).*...
%              sign(log2(metaboliteDiff{:,...
%              ismember(metaboliteDiff.Properties.VariableNames,...
%              {['Serum_fc_', dietChangeNames{diet_i}]})}));
% end
% % check whether metabolites are present in liver and serum clusters
% met_liver_cluster = sum(annotationTableSpatialClusters{:,cellfun(@(x) contains(x, 'spatial_Liver'),...
%     annotationTableSpatialClusters.Properties.VariableNames)},2);
% met_serum_cluster = sum(annotationTableSpatialClusters{:,cellfun(@(x) contains(x, 'spatial_Serum'),...
%     annotationTableSpatialClusters.Properties.VariableNames)},2);
% % total annotated liver metabolites
% nnz(met_liver_cluster & annotationTableSpatialClusters.MetaboliteFilter)
% %ans = 1733
% 
% % changing between mouse groups
% nnz(met_liver_cluster & annotationTableSpatialClusters.MetaboliteFilter &...
%     (sum(abs(met_liver_changes),2)>0))
% %ans = 160
% 
% %fraction of changing liver metabolites between mouse groups
% nnz(met_liver_cluster & annotationTableSpatialClusters.MetaboliteFilter &...
%     (sum(abs(met_liver_changes),2)>0)) / ...
%     nnz(met_liver_cluster & annotationTableSpatialClusters.MetaboliteFilter)
% % ans = 0.0923
% 
% 
% % total annotates serum metabolites
% % nnz(met_serum_cluster & annotationTableSpatialClusters.MetaboliteFilter)
% % ans = 1241
% 
% % changing between mouse groups
% nnz(met_serum_cluster & annotationTableSpatialClusters.MetaboliteFilter &...
%     (sum(abs(met_serum_changes),2)>0))
% % ans = 294
% 
% %fraction of changing serum metabolites between mouse groups
% nnz(met_serum_cluster & annotationTableSpatialClusters.MetaboliteFilter &...
%     (sum(abs(met_serum_changes),2)>0)) / ...
%     nnz(met_serum_cluster & annotationTableSpatialClusters.MetaboliteFilter)
% % ans = 0.2369
% 
% % serum and liver together (fraction)
% nnz((met_liver_cluster & annotationTableSpatialClusters.MetaboliteFilter &...
%     (sum(abs(met_liver_changes),2)>0)) |...
%     (met_serum_cluster & annotationTableSpatialClusters.MetaboliteFilter &...
%     (sum(abs(met_serum_changes),2)>0))) / ...
%     nnz((met_liver_cluster | met_serum_cluster) & annotationTableSpatialClusters.MetaboliteFilter)
% % ans = 0.2032
% 
% % check for clusters
% hier_clust_attribution = readtable([resultsFolder, 'table_hierarchical_clustering_groups.csv']);
% % number of bacterial products
% nnz((hier_clust_attribution.HierarchicalClustGroup==10) | (hier_clust_attribution.HierarchicalClustGroup==11))
% % ans = 244
% 
% % of those, changing in system
% nnz(((hier_clust_attribution.HierarchicalClustGroup==10) |...
%     (hier_clust_attribution.HierarchicalClustGroup==11)) &...
%     ( (sum(abs(met_liver_changes),2)>0) |...
%       (sum(abs(met_serum_changes),2)>0) ) & ...
%       annotationTableSpatialClusters.MetaboliteFilter )
% % ans = 20
% 
% % of those, higher in system in colonized mice
% nnz(( (hier_clust_attribution.HierarchicalClustGroup==10) |...
%       (hier_clust_attribution.HierarchicalClustGroup==11) ) &...
%     ( ((met_liver_changes(:,1)>0) | (met_liver_changes(:,2)>0) |...
%        (met_serum_changes(:,1)>0) | (met_serum_changes(:,2)>0) ) ) )
% % ans = 14
%    
% % check with substrates and products   
% x_met_smooth_norm = x_met_smooth(:,2:end);
% for i=1:size(x_met_smooth_norm,1)
%     x_met_smooth_norm(i,:) = x_met_smooth_norm(i,:)/max(abs(x_met_smooth_norm(i,:)));
% end
% 
% % products
% nnz(((x_met_smooth_norm(:,end)>=0.5)|x_met_smooth_norm(:,end-1)>=0.5) & (x_data_corr>=0.7) &...
%     (modelingResults.MetaboliteFilter==1))
% % ans = 569
% nnz(((x_met_smooth_norm(:,end)>=0.5)|x_met_smooth_norm(:,end-1)>=0.5) &...
%      (x_data_corr>=0.7) &...
%      (modelingResults.MetaboliteFilter==1) &...
%      ((met_liver_changes(:,1)>0) | (met_liver_changes(:,2)>0) |...
%       (met_serum_changes(:,1)>0) | (met_serum_changes(:,2)>0) ) )
% % ans = 27
% % any changes
% nnz(((x_met_smooth_norm(:,end)>=0.5)|x_met_smooth_norm(:,end-1)>=0.5) &...
%      (x_data_corr>=0.7) &...
%      (modelingResults.MetaboliteFilter==1) &...
%      ( (sum(abs(met_liver_changes),2)>0) |...
%        (sum(abs(met_serum_changes),2)>0)) )
% % ans = 45
% % substrates
% nnz(((x_met_smooth_norm(:,end)<=-0.5)|x_met_smooth_norm(:,end-1)<=-0.5) &...
%     (x_data_corr>=0.7) & (modelingResults.MetaboliteFilter==1))
% % ans = 547
% nnz(((x_met_smooth_norm(:,end)<=-0.5)|x_met_smooth_norm(:,end-1)<=-0.5) &...
%      (x_data_corr>=0.7) &...
%      (modelingResults.MetaboliteFilter==1) &...
%      ((met_liver_changes(:,1)<0) | (met_liver_changes(:,2)<0) |...
%       (met_serum_changes(:,1)<0) | (met_serum_changes(:,2)<0) ))
% % ans = 21
% % any change
% nnz(((x_met_smooth_norm(:,end)<=-0.5)|x_met_smooth_norm(:,end-1)<=-0.5) &...
%      (x_data_corr>=0.7) &...
%      (modelingResults.MetaboliteFilter==1) &...
%      ( (sum(abs(met_liver_changes),2)>0) |...
%        (sum(abs(met_serum_changes),2)>0)) )
% % ans = 60
