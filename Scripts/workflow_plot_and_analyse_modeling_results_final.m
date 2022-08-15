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
annotationTableSpatialClusters = readtable([inputFolder ...
    'metabolites_allions_combined_formulas_with_metabolite_filters_spatial100clusters_with_mean.csv']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modelingResults = readtable([outputFolder ...
            'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_allions.csv']);
x_met_smooth = modelingResults{:, width(modelingResults)-8:end};
% get correlations calculated with reverse problem
if isnumeric(modelingResults.ReciprocalCorr(1))
    x_data_corr = modelingResults.ReciprocalCorr;
else
    % it is not numeric, probably contains NaN - convert to numeric
    x_data_corr = cellfun(@(x) str2double(x), modelingResults.ReciprocalCorr);
end
coefvalues = modelingResults.Properties.VariableNames(width(modelingResults)-8:end);

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
% at this point colorbar and print to figure have to be done manually in
% the clustergram window
% print to figure
orient landscape
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
    [figureFolder,...
    'fig4b_clustergram_2LIhos1LIbact_model_coefs_annotated_ions_in_GIT_Rmodel_corr_0_7'])

      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clustergrams were manually selected from the clustergram window
% they can be loaded from 'cgo_clustergrams_of_model_coefficients.mat'
load([outputFolder 'cgo_clustergrams_of_model_coefficients.mat'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print cluster size info
cluster_order = [{cellfun(@(x) (str2double(x)), cgoHColhfd.ColumnLabels)},...
        {setdiff(cellfun(@(x) (str2double(x)), cgoHColboth.ColumnLabels),...
        cellfun(@(x) (str2double(x)), cgoHColhfd.ColumnLabels))},...
        {cellfun(@(x) (str2double(x)),cgoHSIboth_neg.ColumnLabels)},...
        {cellfun(@(x) (str2double(x)),cgoBLIctr_neg.ColumnLabels)},...
        {cellfun(@(x) (str2double(x)),cgoHLIhfd_neg.ColumnLabels)},...
        {cellfun(@(x) (str2double(x)),cgoHLIctr_neg.ColumnLabels)},...
        {cellfun(@(x) (str2double(x)),cgoHCol_both_slight_pos.ColumnLabels)},...
        {cellfun(@(x) (str2double(x)),cgoHLIhfd_pos.ColumnLabels)},...
        {cellfun(@(x) (str2double(x)),cgoHLIctr_pos.ColumnLabels)},...
        {cellfun(@(x) (str2double(x)),cgoBLIhfd_pos.ColumnLabels)},...
        {cellfun(@(x) (str2double(x)),cgoBLIctr_pos.ColumnLabels)},...
        {cellfun(@(x) (str2double(x)),cgo_mixed.ColumnLabels)}];
    
cluster_sizes = [length(cgoHColhfd.ColumnLabels),...
        length(cgoHColboth.ColumnLabels)-length(cgoHColhfd.ColumnLabels),...
        length(cgoHSIboth_neg.ColumnLabels),...
        length(cgoBLIctr_neg.ColumnLabels),...
        length(cgoHLIhfd_neg.ColumnLabels),...
        length(cgoHLIctr_neg.ColumnLabels),...
        length(cgoHCol_both_slight_pos.ColumnLabels),...
        length(cgoHLIhfd_pos.ColumnLabels),...
        length(cgoHLIctr_pos.ColumnLabels),...
        length(cgoBLIhfd_pos.ColumnLabels),...
        length(cgoBLIctr_pos.ColumnLabels),...
        length(cgo_mixed.ColumnLabels)];
        
cluster_attribution = zeros(size(annotationTableSpatialClusters,1),1);
for i=1:length(cluster_order)
    cluster_attribution(cluster_order{i})=i;
end
% add all other metabolites to a thirteens group
cluster_attribution(setdiff(clustidx, find(cluster_attribution>0))) = length(cluster_order)+1; 
hier_clust_attribution = table(annotationTableSpatialClusters.MZ,...
                                annotationTableSpatialClusters.RT,... 
                                annotationTableSpatialClusters.MetaboliteFilter,...
                                cluster_attribution,...
                                'VariableNames', {'MZ', 'RT', 'MetaboliteFilter', 'HierarchicalClustGroup'});
writetable(hier_clust_attribution, [outputFolder, ...
    'table_hierarchical_clustering_groups.csv']);
                            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot cluster sizes as stacked bar plot
figure
barh([cluster_sizes;cluster_sizes], 'stacked')
xlim([0 sum(cluster_sizes)])
for i=1:length(cluster_sizes)
    text(sum(cluster_sizes(1:(i-1))), 1, num2str(cluster_sizes(i)))
    text(sum(cluster_sizes(1:(i-1))), 2, num2str(round(100*cluster_sizes(i)/length(cgo.ColumnLabels))))
end
orient landscape
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
    [figureFolder, ...
    'fig4b_barh_2LIhos1LIbact_model_cluster_sizes_percent.pdf'])

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
hclusterTables = hclusterTables_hmdbClass;
hclusterTables = hclusterTables_hmdbSuperClass;
hclusterTables = hclusterTables_hmdbSubClass;

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

%selectpathways = sum(fdrMet<=0.1,2)>0;
selectpathways = sum(pMet<0.01,2)>0;

% Generate cluster names
displaymat = pMet(selectpathways,:);
%displaymat = fdrMet(selectpathways,:);

displaymat(displaymat>0.1)=1;

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
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
            [figureFolder,...
            'clustergram_PTWENR_hcluster_hmdbCombinedClasses_2LIhos1LIbact_pFDR_le0_1.pdf'])

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
               'ptwenr_', curvarName, '_hcluster_2LIhos1LIbact_hmdbv4.csv']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(2,3,1)
i=1;
j=4;
%scatter(clustdata(end-1,:), clustdata(end,:))
scatter(clustdata(i,:), clustdata(j,:),...
    'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
xlabel(clustrows(i))
ylabel(clustrows(j))
[corrc, corrp] = corr(clustdata(i,:)', clustdata(j,:)', 'rows', 'complete');
title(sprintf('PCC=%.2f, p=%.2f', corrc, corrp))

axis square

subplot(2,3,2)
%scatter(clustdata(2,:), clustdata(5,:))
i=2;
j=5;
%scatter(clustdata(end-1,:), clustdata(end,:))
scatter(clustdata(i,:), clustdata(j,:),...
    'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
xlabel(clustrows(i))
ylabel(clustrows(j))
[corrc, corrp] = corr(clustdata(i,:)', clustdata(j,:)', 'rows', 'complete');
title(sprintf('PCC=%.2f, p=%.2f', corrc, corrp))

axis square

subplot(2,3,3)
%scatter(clustdata(3,:), clustdata(6,:))
i=3;
j=6;
%scatter(clustdata(end-1,:), clustdata(end,:))
scatter(clustdata(i,:), clustdata(j,:),...
    'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
xlabel(clustrows(i))
ylabel(clustrows(j))
[corrc, corrp] = corr(clustdata(i,:)', clustdata(j,:)', 'rows', 'complete');
title(sprintf('PCC=%.2f, p=%.2f', corrc, corrp))

axis square

subplot(2,3,4)
%scatter(clustdata(2,:), clustdata(3,:))
i=2;
j=3;
%scatter(clustdata(end-1,:), clustdata(end,:))
scatter(clustdata(i,:), clustdata(j,:),...
    'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
xlabel(clustrows(i))
ylabel(clustrows(j))
[corrc, corrp] = corr(clustdata(i,:)', clustdata(j,:)', 'rows', 'complete');
title(sprintf('PCC=%.2f, p=%.2f', corrc, corrp))

axis square

subplot(2,3,5)
%scatter(clustdata(5,:), clustdata(6,:))
i=5;
j=6;
%scatter(clustdata(end-1,:), clustdata(end,:))
scatter(clustdata(i,:), clustdata(j,:),...
    'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
xlabel(clustrows(i))
ylabel(clustrows(j))
[corrc, corrp] = corr(clustdata(i,:)', clustdata(j,:)', 'rows', 'complete');
title(sprintf('PCC=%.2f, p=%.2f', corrc, corrp))
axis square

subplot(2,3,6)
i=size(clustdata,1)-1;
j=size(clustdata,1);
%scatter(clustdata(end-1,:), clustdata(end,:))
scatter(clustdata(i,:), clustdata(j,:),...
    'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
xlabel(clustrows(i))
ylabel(clustrows(j))
[corrc, corrp] = corr(clustdata(i,:)', clustdata(j,:)', 'rows', 'complete');
title(sprintf('PCC=%.2f, p=%.2f', corrc, corrp))
axis square

orient landscape
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
    [figureFolder...
    'scatter_model2LIhost1LIbact_coefs_circles_transparent.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot bacterial coefficients separately with marking of potential
% substrates and products
figure
i=size(clustdata,1)-1;
j=size(clustdata,1);
%scatter(clustdata(i,:), clustdata(j,:), 'k')
potential_products = (clustdata(i,:)>=0.5) | (clustdata(j,:)>=0.5);
potential_substrates = (clustdata(i,:)<=-0.5) | (clustdata(j,:)<=-0.5);
potential_nothing = ~(potential_substrates | potential_products);
scatter(clustdata(i,potential_nothing), clustdata(j,potential_nothing),...
    'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
hold on
scatter(clustdata(i,potential_products), clustdata(j,potential_products),...
    'MarkerFaceColor','r','MarkerEdgeColor','r',...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
scatter(clustdata(i,potential_substrates), clustdata(j,potential_substrates),...
    'MarkerFaceColor','b','MarkerEdgeColor','b',...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)

xlabel(clustrows(i))
ylabel(clustrows(j))
[corrc, corrp] = corr(clustdata(i,:)', clustdata(j,:)', 'rows', 'complete');
title({sprintf('PCC=%.2f, p=%.2f', corrc, corrp),...
       sprintf('SubN=%d, ProdN=%d',nnz(potential_substrates),nnz(potential_products))});
   
axis square

orient landscape
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
    [figureFolder...
    'scatter_model_coefs_annotated_ions_in_GIT_Rmodel_corr_0_7'])  
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate mass differences between substrates and products
% not such a good idea as there are too many cases where the mass shift
% fits some pattern, but without structure there are too many candidates
% met_pos_hfd = clustidx(clustdata(ismember(clustrows, 'B1LIhfd'),:)>0.5);
% met_neg_hfd = clustidx(clustdata(ismember(clustrows, 'B1LIhfd'),:)<-0.5);
% 
% single_mzdelta = zeros(length(met_neg_hfd), length(met_pos_hfd));
% for i=1:length(met_neg_hfd)
%     for j=1:length(met_pos_hfd)
%         single_mzdelta(i,j) = annotationTableSpatialClusters.MZ(met_neg_hfd(i))-...
%             annotationTableSpatialClusters.MZ(met_pos_hfd(j));
%     end
% end

cmpd_interest =  'C00931';%porph' 'C00430'; %5-aminolevulinate 'C00025'; % glu'C00327'; %citrulline    

cmpd_interest_idx = find(cellfun(@(x) contains(x,cmpd_interest), annotationTableSpatialClusters.CompoundID) &...
                    (annotationTableSpatialClusters.MetaboliteFilter==1));

nnz(clustidx(potential_substrates)==cmpd_interest_idx)
nnz(clustidx(potential_products)==cmpd_interest_idx)
clustdata(end-1:end, clustidx==cmpd_interest_idx)

% i = find(met_neg_hfd==cmpd_interest_idx);
% test = sort(single_mzdelta(i,:))';
% %j = find(abs(single_mzdelta(i,:)+18.010)<=0.002)
% j = find(met_pos_hfd==cmpd_interest_idx)
% single_mzdelta(i,j)
% 
% annotationTableSpatialClusters.CompoundName(met_pos_hfd(j))
% annotationTableSpatialClusters.CompoundID(met_pos_hfd(j))
% annotationTableSpatialClusters.CompoundID(met_neg_hfd(i))
cmpd_interest = 'C00931';%porphobilonogen
cmpd_interest = 'C03207';%Dehydroprogesterone
cmpd_interest = 'C00025';%Glutamate
cmpd_interest = 'HMDB00552'; %fatty acid esters
cmpd_interest = 'C15636'; %ethylchenodeoxycholic acid
cmpd_interest = 'C01595'; %linoleate
cmpd_interest = 'C16901'; %latrunculin
cmpd_interest = 'C14983'; %quinone
cmpd_interest = 'C06492'; %Compound II(R/S)
cmpd_interest = 'C13824'; %R-L3

cmpd_interest = 'C00065'; %serine
cmpd_interest = 'C00078'; %tryptophane
cmpd_interest = 'C00258'; %glycerate
cmpd_interest = 'C00214'; %thymidine
cmpd_interest = 'C00178'; %thymine
cmpd_interest = 'C09675'; %pentadecaheptane
cmpd_interest = 'C00299'; %uridine
cmpd_interest = 'C00106'; %uracil
cmpd_interest = 'C01909'; %dethiobiotin
cmpd_interest = 'C00120'; %biotin
cmpd_interest = 'C21357'; %dethiobiotin
cmpd_interest = 'C00041'; %alanine
cmpd_interest = 'C00123'; %leucine
cmpd_interest = 'C05589'; %normethanephrine
cmpd_interest = 'C00355'; %L-dopa
cmpd_interest = 'C00547'; %noradrenaline
cmpd_interest = 'C00082'; %Tyrosine

cmpd_interest_idx = cellfun(@(x) contains(x,cmpd_interest), annotationTable.CompoundID) &...
                    (annotationTable.MetaboliteFilter==1);
                
  
pos_neg_corr_all = zeros(size(single_mzdelta));
pos_neg_corr_li = zeros(size(single_mzdelta));
pos_neg_corr_limax = zeros(size(single_mzdelta));

for i=1:length(met_neg_hfd)
    for j=1:length(met_pos_hfd)
        dataOrigNeg = reshape(kmean_vector_joint_orig(:,selected_mets==met_neg_hfd(i)),4,[]);
        dataOrigPos = reshape(kmean_vector_joint_orig(:,selected_mets==met_pos_hfd(j)),4,[]);

        pos_neg_corr_all(i,j) = corr(dataOrigNeg(:), dataOrigPos(:));
        pos_neg_corr_li(i,j) = corr(reshape(dataOrigNeg(:,4:6),[],1), reshape(dataOrigPos(:,4:6),[],1));
        corrdiag = diag(corr(dataOrigNeg(:,4:6), dataOrigPos(:,4:6)));
        [~, maxidx] = max(-(corrdiag));
        pos_neg_corr_limax(i,j) = corrdiag(maxidx);
     
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot all annotated metabolites across tissues with model coefficients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data and restored data from file
% save model results to file - reciprocal data restoration
modelData = readtable([outputFolder...
    'model_results_SMOOTH_normbyabsmax_reciprocal_problem_allions.csv']);
modelData_data = modelData(:, 9:end);
modelData_orig = modelData_data{:, cellfun(@(x) ~(contains(x, 'Recip') |...
                                             contains(x, 'Random') ),...
                                             modelData_data.Properties.VariableNames)};
modelData_orig_cols = modelData_data.Properties.VariableNames(cellfun(@(x) ~(contains(x, 'Recip') |...
                                             contains(x, 'Random') ),...
                                             modelData_data.Properties.VariableNames));
modelData_recip = modelData_data{:, cellfun(@(x) contains(x, 'Recip'),...
                                             modelData_data.Properties.VariableNames)};

spatialClusters_names = annotationTableSpatialClusters.Properties.VariableNames(...
    cellfun(@(x) contains(x, 'spatial_clust100'), annotationTableSpatialClusters.Properties.VariableNames));
spatialClusters = annotationTableSpatialClusters{:, spatialClusters_names};                                   
spatialClusters_names = cellfun(@(x) strrep(x, 'spatial_clust100_', ''), spatialClusters_names, 'unif', 0);

%fileNameprofiles = 'profiles_annotated_mets_modelSMOOTH_2LIcoefHost_1LIbact_R_corr07_part11.ps';
%fileNameprofiles = 'profiles_propionate_modelSMOOTH_2LIcoefHost_1LIbact.ps';
%targetMZ = 74.037;
%fileNameprofiles = 'profiles_riboflavin_modelSMOOTH_2LIcoefHost_1LIbact.ps';
%targetMZ = 376.139;
%fileNameprofiles = 'profiles_biotin_modelSMOOTH_2LIcoefHost_1LIbact.ps';
%targetMZ = 244.088;

%fileNameprofiles = 'profiles_niacin_modelSMOOTH_2LIcoefHost_1LIbact.ps';
%targetMZ = 123.032;

%fileNameprofiles = 'profiles_thymidine_modelSMOOTH_2LIcoefHost_1LIbact.ps';
%targetMZ = 242.090;

%fileNameprofiles = 'profiles_thymine_modelSMOOTH_2LIcoefHost_1LIbact.ps';
%targetMZ = 126.042;

%fileNameprofiles = 'profiles_cytidine_modelSMOOTH_2LIcoefHost_1LIbact.ps';
%targetMZ = 243.085;

% fileNameprofiles = 'profiles_cytosine_modelSMOOTH_2LIcoefHost_1LIbact.ps';
%targetMZ = 111.043;

%fileNameprofiles = 'profiles_uracil_modelSMOOTH_2LIcoefHost_1LIbact.ps'; 	
%targetMZ = 112.027;

%fileNameprofiles = 'profiles_uridine_modelSMOOTH_2LIcoefHost_1LIbact.ps'; 	
%targetMZ = 244.069;

% fileNameprofiles = 'profiles_propionylPhosphate_modelSMOOTH_2LIcoefHost_1LIbact.ps'; 	
% targetMZ = 154.003;

fileNameprofiles = 'profiles_figureselected_modelSMOOTH_2LIcoefHost_1LIbact.ps'; 	
targetMZ = 147.053;%74.037;%287.210;%

mycolors = [0 115 178;... %dark blue
            204 227 240;...%light blue
            211 96 39;... %dark orange
            246 223 212]/256;%light orange
git_labels = {'Du', 'Je', 'Il', 'Cec', 'Col', 'Fec'};


%compoundsInterest = find((annotationTableSpatialClusters.MetaboliteFilter>0) &...
%                (x_data_corr>=0.7));
%compoundsInterest = find((abs(metaboliteFilters.MZ-targetMZ)<=0.001));

compoundsInterest = find((annotationTableSpatialClusters.MetaboliteFilter>0) &...
               (x_data_corr>=0.7) &...
               (abs(metaboliteFilters.MZ-targetMZ)<=0.001));
           
curData_cols = reshape(modelData_orig_cols, 6, 4);
fig = figure('units','normalized','outerposition',[0 0 1 1]);

for cpdix=1:length(compoundsInterest) %

    
    testidx = compoundsInterest(cpdix);

    testann = annotationTableSpatialClusters.CompoundName(testidx,:);
    testannID = annotationTableSpatialClusters.CompoundID(testidx,:);
    testmz = annotationTableSpatialClusters.MZ(testidx,:);
    testrt = annotationTableSpatialClusters.RT(testidx,:);
    testmethod = annotationTableSpatialClusters.Method(testidx,:);
    testmode = annotationTableSpatialClusters.Mode(testidx,:);
 
    spx=1;
    spy=3;
    spidx = 1;
    coloridx = 1;
    idx=1;
    curmat = zeros(4,6);
    cur_data = reshape(modelData_orig(testidx,:), 4, 6)';
    cur_rdata = reshape(modelData_recip(testidx,:), 6, 4);
    
    %normalize rdata to max
    cur_rdata = (cur_rdata-0.8);
    cur_rdata = cur_rdata/max(max(cur_rdata));
    
    legend_entries = cell(4,1);
    for i = 1:size(cur_data,2)

        subplot(spx,spy,idx)
        hold on
        h(i) = plot(cur_data(:,i),...
             'LineWidth', 2,...
             'Color', mycolors(i,:));
        ylabel('Original normbymax')

        subplot(spx,spy,idx+1);
        hold on
        plot((cur_rdata(:,i)),...
                 'LineWidth', 2,...
                 'Color', mycolors(i,:));
        ylabel('Restored normbymax')

    end
    title(sprintf('%s %s %s %s: %d %d %d %d PCC=%.2f',...
                        spatialClusters_names{:},...
                        spatialClusters(testidx,:),...
                        x_data_corr(testidx)))

    legend(h, curData_cols(1,:))%, 'Location', 'bestoutside');
    
    subplot(spx,spy,idx+2);
    
    curcoefs = x_met_smooth(testidx, 2:end);
    barh(curcoefs./max(abs(curcoefs)))
    set(gca, 'YTick', 1:length(curcoefs));
    set(gca, 'YTickLabel', coefvalues(2:end));
    set(gca, 'YDir','reverse')
    ylim([0.5 length(curcoefs)+0.5])
    xlim([-1 1]);
    axis square
    
    for spi = 1:(spx*spy)-1
        subplot(spx,spy,spi)
        set(gca, 'XTick', 1:6)
        xlim([1 6])
        ylim([0 1])
        set(gca, 'XTick', 1:length(git_labels))
        set(gca, 'XTickLabel', git_labels)
        
        axis square
    end
    spt = suptitle({sprintf('MZ=%.3f',testmz(idx,1)),...
                                        testannID{1},...
                                        testann{1}});
    set(spt,'FontSize',8,'FontWeight','normal')
    orient landscape
    %export_fig output.pdf -pdf –append
    print(gcf, '-painters', '-dpsc2', '-r600', '-append', '-bestfit',...
            [figureFolder, fileNameprofiles])
   % print('-dpsc2', '-painters', '-noui', '-append', fileNameprofiles);
    clf('reset')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze bacterial coefficients
plot(x_met_mean(9,:))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate changing metabolites in the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read metabolite data
metaboliteFolder = '.\metabolomics\ProcessedData\';
metaboliteDiff = readtable([metaboliteFolder,...
    'table_diff_abundance_metabolite_ions_removed2outliers.csv']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get information from differential analysis on whether metabolites are 
% differentially abundant in either serum or liver in at least
% one diet between GF and colonized mice
dietChangeNames = {'HFD_DCGF', 'CTR_DCGF'};
fcThreshold = log2(1.5);
pfdrThreshold = 0.05;
pvalColumn = '_pFDR_'; % choose between _p_ and _pFDR_
% record also liver and serum changes
met_liver_changes = zeros(size(metaboliteDiff,1), length(dietChangeNames));
met_serum_changes = zeros(size(metaboliteDiff,1), length(dietChangeNames));
for diet_i=1:length(dietChangeNames)
        % liver
        % (record sign of change)
        met_liver_changes(:,diet_i) = (abs(log2(metaboliteDiff{:,...
             ismember(metaboliteDiff.Properties.VariableNames,...
             {['Liver_fc_', dietChangeNames{diet_i}]})})) >= fcThreshold).*...
             (metaboliteDiff{:,...
             ismember(metaboliteDiff.Properties.VariableNames,...
             {['Liver' pvalColumn, dietChangeNames{diet_i}]})} <= pfdrThreshold).*...
             sign(log2(metaboliteDiff{:,...
             ismember(metaboliteDiff.Properties.VariableNames,...
             {['Liver_fc_', dietChangeNames{diet_i}]})}));
        % serum
        % (record sign of change)
        met_serum_changes(:,diet_i) = (abs(log2(metaboliteDiff{:,...
             ismember(metaboliteDiff.Properties.VariableNames,...
             {['Serum_fc_', dietChangeNames{diet_i}]})})) >= fcThreshold).*...
             (metaboliteDiff{:,...
             ismember(metaboliteDiff.Properties.VariableNames,...
             {['Serum' pvalColumn, dietChangeNames{diet_i}]})} <= pfdrThreshold).*...
             sign(log2(metaboliteDiff{:,...
             ismember(metaboliteDiff.Properties.VariableNames,...
             {['Serum_fc_', dietChangeNames{diet_i}]})}));
end
% check whether metabolites are present in liver and serum clusters
met_liver_cluster = sum(annotationTableSpatialClusters{:,cellfun(@(x) contains(x, 'spatial_Liver'),...
    annotationTableSpatialClusters.Properties.VariableNames)},2);
met_serum_cluster = sum(annotationTableSpatialClusters{:,cellfun(@(x) contains(x, 'spatial_Serum'),...
    annotationTableSpatialClusters.Properties.VariableNames)},2);
% total annotated liver metabolites
nnz(met_liver_cluster & annotationTableSpatialClusters.MetaboliteFilter)
% changing between mouse groups
nnz(met_liver_cluster & annotationTableSpatialClusters.MetaboliteFilter &...
    (sum(abs(met_liver_changes),2)>0))
%fractions
nnz(met_liver_cluster & annotationTableSpatialClusters.MetaboliteFilter &...
    (sum(abs(met_liver_changes),2)>0)) / ...
    nnz(met_liver_cluster & annotationTableSpatialClusters.MetaboliteFilter)

% total annotates serum metabolites
nnz(met_serum_cluster & annotationTableSpatialClusters.MetaboliteFilter)
% changing between mouse groups
nnz(met_serum_cluster & annotationTableSpatialClusters.MetaboliteFilter &...
    (sum(abs(met_serum_changes),2)>0))
%fractions
nnz(met_serum_cluster & annotationTableSpatialClusters.MetaboliteFilter &...
    (sum(abs(met_serum_changes),2)>0)) / ...
    nnz(met_serum_cluster & annotationTableSpatialClusters.MetaboliteFilter)

% serum and liver together
nnz((met_liver_cluster & annotationTableSpatialClusters.MetaboliteFilter &...
    (sum(abs(met_liver_changes),2)>0)) |...
    (met_serum_cluster & annotationTableSpatialClusters.MetaboliteFilter &...
    (sum(abs(met_serum_changes),2)>0))) / ...
    nnz((met_liver_cluster | met_serum_cluster) & annotationTableSpatialClusters.MetaboliteFilter)

% check for clusters
hier_clust_attribution = readtable([outputFolder, 'table_hierarchical_clustering_groups.csv']);
% number of bacterial products
nnz((hier_clust_attribution.HierarchicalClustGroup==10) | (hier_clust_attribution.HierarchicalClustGroup==11))

% changing in system
nnz((hier_clust_attribution.HierarchicalClustGroup==10) | (hier_clust_attribution.HierarchicalClustGroup==11) &...
    ( ((sum(abs(met_liver_changes),2)>0) | (sum(abs(met_serum_changes),2)>0)) & ...
      annotationTableSpatialClusters.MetaboliteFilter ) )
% higher in system
nnz(( (hier_clust_attribution.HierarchicalClustGroup==10) |...
      (hier_clust_attribution.HierarchicalClustGroup==11) ) &...
    ( ((met_liver_changes(:,1)>0) | (met_liver_changes(:,2)>0) |...
       (met_serum_changes(:,1)>0) | (met_serum_changes(:,2)>0) ) ) )
   
% check with substrates and products   
x_met_smooth_norm = x_met_smooth(:,2:end);
for i=1:size(x_met_smooth_norm,1)
    x_met_smooth_norm(i,:) = x_met_smooth_norm(i,:)/max(abs(x_met_smooth_norm(i,:)));
end

% products
nnz(((x_met_smooth_norm(:,end)>=0.5)|x_met_smooth_norm(:,end-1)>=0.5) & (x_data_corr>=0.7) & (metaboliteFilters.MetaboliteFilter==1))
% ans = 569
nnz(((x_met_smooth_norm(:,end)>=0.5)|x_met_smooth_norm(:,end-1)>=0.5) &...
     (x_data_corr>=0.7) &...
     (metaboliteFilters.MetaboliteFilter==1) &...
     ((met_liver_changes(:,1)>0) | (met_liver_changes(:,2)>0) |...
      (met_serum_changes(:,1)>0) | (met_serum_changes(:,2)>0) ) )
% ans = 27
% any changes
nnz(((x_met_smooth_norm(:,end)>=0.5)|x_met_smooth_norm(:,end-1)>=0.5) &...
     (x_data_corr>=0.7) &...
     (metaboliteFilters.MetaboliteFilter==1) &...
     ( (sum(abs(met_liver_changes),2)>0) |...
       (sum(abs(met_serum_changes),2)>0)) )
% ans = 45
% substrates
nnz(((x_met_smooth_norm(:,end)<=-0.5)|x_met_smooth_norm(:,end-1)<=-0.5) & (x_data_corr>=0.7) & (metaboliteFilters.MetaboliteFilter==1))
% ans = 547
nnz(((x_met_smooth_norm(:,end)<=-0.5)|x_met_smooth_norm(:,end-1)<=-0.5) &...
     (x_data_corr>=0.7) &...
     (metaboliteFilters.MetaboliteFilter==1) &...
     ((met_liver_changes(:,1)<0) | (met_liver_changes(:,2)<0) |...
      (met_serum_changes(:,1)<0) | (met_serum_changes(:,2)<0) ))
% 21
% any change
nnz(((x_met_smooth_norm(:,end)<=-0.5)|x_met_smooth_norm(:,end-1)<=-0.5) &...
     (x_data_corr>=0.7) &...
     (metaboliteFilters.MetaboliteFilter==1) &...
     ( (sum(abs(met_liver_changes),2)>0) |...
       (sum(abs(met_serum_changes),2)>0)) )
% 60
