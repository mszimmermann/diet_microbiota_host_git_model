%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get substrate-product pathes for potential substrates and products
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% call script defining file dependenciesand global variables
addpath(genpath('.\'));
add_global_and_file_dependencies

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Required files:
% 'model_results_SMOOTH_normbyabsmax_ONLYMETCOEF_2LIcoefHost1LIcoefbact_allions.csv'
% matrices in the form 'keggCompoundRXNadjacency_shortest_path_1_2.mat' in  'KEGGreaction_path\' folder
% Output: 
% Files:
% 'table_potential_substrates_ids.csv'
% 'table_potential_products_ids.csv'
% 'table_shortestMatrix_paths_filtered.csv'
% Figures:
% 'fig_5a_histogram_bacterial_coefs_filtered_corr_0_7_ybreak.pdf'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model_coefs_annotated = readtable([resultsFolder ...
    'model_results_SMOOTH_normbyabsmax_ONLYMETCOEF_2LIcoefHost1LIcoefbact_allions.csv'],...
    'TreatAsEmpty', '#NUM!');

endcol = length(model_coefs_annotated.Properties.VariableNames);
startcol = endcol-8;
for i=startcol:endcol
    testval = model_coefs_annotated(:,model_coefs_annotated.Properties.VariableNames{i});
    if ~isnumeric(testval{1,1})
        model_coefs_annotated(:,model_coefs_annotated.Properties.VariableNames{i}) = ...
            array2table(cellfun(@(x) str2double(x), testval{:,1}));
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot coefficients of bacterial metabolism
model_coefs_annotated_filtered = model_coefs_annotated(...
                     (model_coefs_annotated.MetaboliteFilter==1) & ...
                     (model_coefs_annotated.ReciprocalCorr>=0.7),:);
figure
subplot(2,1,1)
histogram(model_coefs_annotated_filtered.B1LIctr,50)
hold on
histogram(model_coefs_annotated_filtered.B1LIhfd,50)
%xlabel('Normalized bacterial coefficient')
ylabel('Number of metabolites')
legend({'CTR', 'HFD'}, 'location', 'EastOutside')
%axis square
plot([-0.5, -0.5], [0, 1200], 'k--')
plot([0.5, 0.5], [0, 1200], 'k--')
ylim([800, 1100])
% add text with numbers
text(0.6, 1000, sprintf('N = %d', nnz( (model_coefs_annotated_filtered.B1LIctr>=0.5) |...
                                        (model_coefs_annotated_filtered.B1LIhfd>=0.5) )));
text(-0.9, 1000, sprintf('N = %d', nnz( (model_coefs_annotated_filtered.B1LIctr<=-0.5) |...
                                        (model_coefs_annotated_filtered.B1LIhfd<=-0.5) )));
                                    
subplot(2,1,2)
histogram(model_coefs_annotated_filtered.B1LIctr,50)
hold on
histogram(model_coefs_annotated_filtered.B1LIhfd,50)
xlabel('Normalized bacterial coefficient')
ylabel('Number of metabolites')
legend({'CTR', 'HFD'}, 'location', 'EastOutside')
%axis square
plot([-0.5, -0.5], [0, 1200], 'k--')
plot([0.5, 0.5], [0, 1200], 'k--')
ylim([0, 200])
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
    [figureFolder, ...
    'fig_5a_histogram_bacterial_coefs_filtered_corr_0_7_ybreak.pdf']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate among potential substrates and products, how many of them
% significantly differ in serum and  liver
potential_substrate_idx = (model_coefs_annotated.MetaboliteFilter==1) & ...
                          (model_coefs_annotated.ReciprocalCorr>=0.7) & ...
                          ((model_coefs_annotated.B1LIctr<=-0.5) |...
                           (model_coefs_annotated.B1LIhfd<=-0.5));
potential_product_idx = (model_coefs_annotated.MetaboliteFilter==1) & ...
                          (model_coefs_annotated.ReciprocalCorr>=0.7) & ...
                          ((model_coefs_annotated.B1LIctr>=0.5) |...
                           (model_coefs_annotated.B1LIhfd>=0.5));

% % check how many metabolites change in the system
% metaboliteDiff = readtable([inputFolder,...
%     'table_diff_abundance_metabolite_ions.csv']);
% pfdrThreshold = 0.1;
% % substrates
% liver_substrate_significance  = ((metaboliteDiff{potential_substrate_idx,...
%              ismember(metaboliteDiff.Properties.VariableNames,...
%              {'Liver_pFDR_HFD_DCGF'})} <= pfdrThreshold) |...
%     (metaboliteDiff{potential_substrate_idx,...
%              ismember(metaboliteDiff.Properties.VariableNames,...
%              {'Liver_pFDR_CTR_DCGF'})} <= pfdrThreshold));
% 
% serum_substrate_significance = ((metaboliteDiff{potential_substrate_idx,...
%              ismember(metaboliteDiff.Properties.VariableNames,...
%              {'Serum_pFDR_HFD_DCGF'})} <= pfdrThreshold) |...
%     (metaboliteDiff{potential_substrate_idx,...
%              ismember(metaboliteDiff.Properties.VariableNames,...
%              {'Serum_pFDR_CTR_DCGF'})} <= pfdrThreshold));
% 
%          
% % products
% liver_product_significance = ((metaboliteDiff{potential_product_idx,...
%              ismember(metaboliteDiff.Properties.VariableNames,...
%              {'Liver_pFDR_HFD_DCGF'})} <= pfdrThreshold) |...
%     (metaboliteDiff{potential_product_idx,...
%              ismember(metaboliteDiff.Properties.VariableNames,...
%              {'Liver_pFDR_CTR_DCGF'})} <= pfdrThreshold));
% 
% serum_product_significance = ((metaboliteDiff{potential_product_idx,...
%              ismember(metaboliteDiff.Properties.VariableNames,...
%              {'Serum_pFDR_HFD_DCGF'})} <= pfdrThreshold) |...
%     (metaboliteDiff{potential_product_idx,...
%              ismember(metaboliteDiff.Properties.VariableNames,...
%              {'Serum_pFDR_CTR_DCGF'})} <= pfdrThreshold));
% nnz(liver_substrate_significance | serum_substrate_significance)
% nnz(liver_substrate_significance)
% nnz(serum_substrate_significance)
% nnz(liver_product_significance)
% nnz(serum_product_significance)
% nnz(liver_product_significance | serum_product_significance)        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get potential substrates and products
potential_products = (model_coefs_annotated.MetaboliteFilter==1) & ...
                     (model_coefs_annotated.ReciprocalCorr>=0.7) & ...
                     ((model_coefs_annotated.B1LIctr>=0.5) | (model_coefs_annotated.B1LIhfd>=0.5));
potential_substrates = (model_coefs_annotated.MetaboliteFilter==1) & ...
                     (model_coefs_annotated.ReciprocalCorr>=0.7) & ...
                     ((model_coefs_annotated.B1LIctr<=-0.5) | (model_coefs_annotated.B1LIhfd<=-0.5));
potential_products_ids = model_coefs_annotated.CompoundID(potential_products);
potential_substrates_ids = model_coefs_annotated.CompoundID(potential_substrates);

potential_products_ids_all = cell(length(potential_products_ids)*20,1);
potential_substrates_ids_all = cell(length(potential_substrates_ids)*20,1);
% get all product ids
idx=1;
for i=1:length(potential_products_ids)
    curid = strsplit(potential_products_ids{i},';');
    potential_products_ids_all(idx:idx+length(curid)-1) = curid;
    idx = idx+length(curid);
end
potential_products_ids_all(idx:end) = [];

% get all substrate ids
idx=1;
for i=1:length(potential_substrates_ids)
    curid = strsplit(potential_substrates_ids{i},';');
    potential_substrates_ids_all(idx:idx+length(curid)-1) = curid;
    idx = idx+length(curid);
end
potential_substrates_ids_all(idx:end) = [];

% combine kegg substrate and produc ids
potential_sp_keggids = [potential_products_ids_all; potential_substrates_ids_all];
potential_sp_keggids(cellfun(@(x) ~contains(x, 'cpd:'), potential_sp_keggids))=[];
potential_sp_keggids = unique(potential_sp_keggids);
potential_sp_keggids = cellfun(@(x) strrep(x, 'cpd:', ''), potential_sp_keggids, 'unif', 0);
potential_substrates_ids_all = cellfun(@(x) strrep(x, 'cpd:', ''), potential_substrates_ids_all, 'unif', 0);
potential_products_ids_all = cellfun(@(x) strrep(x, 'cpd:', ''), potential_products_ids_all, 'unif', 0);

% save potential substrate and product IDs to file
temptable = potential_substrates_ids_all;
temptable(cellfun(@(x) isempty(x), temptable))=[];
temptable = unique(temptable);
potential_substrates_ids_all_table = table(temptable,...
    'VariableNames', {'PotentialSubstrateID'});
temptable = potential_products_ids_all;
temptable(cellfun(@(x) isempty(x), temptable))=[];
temptable = unique(temptable);
potential_products_ids_all_table = table(temptable,...
    'VariableNames', {'PotentialProductID'});
clear temptable
writetable(potential_substrates_ids_all_table, ...
    [outputFolder, 'table_potential_substrates_ids.csv']);
writetable(potential_products_ids_all_table, ...
    [outputFolder, 'table_potential_products_ids.csv']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load one of the KEGG adjacency matrices to get IDs
curname = 'keggCompoundRXNadjacency_shortest_path_1_2.mat';
curmat = load([rawdataFolder, 'KEGGreaction_path\',...
                    curname]);
rxnCompounds = curmat.rxnCompounds;

% try to load only pathes for one cluster
[shortestMatrix_sources,~,matrixIDs] = intersect(potential_sp_keggids,rxnCompounds);
shortestMatrix_cluster = cell(length(matrixIDs), length(rxnCompounds));
source_nodes = intersect(potential_sp_keggids,...
                    rxnCompounds);

maxid = max(matrixIDs);
for j=1:2:maxid
    curname = ['keggCompoundRXNadjacency_shortest_path_', ...
                num2str(j), '_',...
                num2str(j+1), '.mat'];
    curmat = load(['C:\Users\mazimmer\Documents\Documents\MATLAB\GutModels\AGORAreaction_path\',...
                    curname]);
    curids = matrixIDs(matrixIDs>(j+1));
    for i=1:length(curids)
        shortestMatrix_cluster(matrixIDs==curids(i),j:j+1) = ...
            curmat.shortestMatrix(:,curids(i));
    end
    % get all the distances between smaller indeces and current index
    curids = matrixIDs((matrixIDs==j) | (matrixIDs==(j+1)));
    for i=1:length(curids)
        if mod(curids(i),2)==0
            shortestMatrix_cluster(matrixIDs==curids(i),(curids(i)+1):end) = ...
                curmat.shortestMatrix(2,(curids(i)+1):end);
        else
            shortestMatrix_cluster(matrixIDs==curids(i),(curids(i)+1):end) = ...
                curmat.shortestMatrix(1,(curids(i)+1):end);
        end
    end
    if mod(j+1,100)==0
        fprintf('Done with matrix %d of %d\n', j, maxid)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% leave only targets that belong to the clusters
shortestMatrix_cluster = shortestMatrix_cluster(:, matrixIDs);

% calculate length of shortest path
shortesMatrix_distance = zeros(size(shortestMatrix_cluster));
for i=1:size(shortestMatrix_cluster)
    for j=1:size(shortestMatrix_cluster)
        curpaths = shortestMatrix_cluster{i,j};
        if ~isempty(curpaths)
            curdist = cellfun(@(x) length(x), curpaths);
            shortesMatrix_distance(i,j) = min(curdist);
        end
    end
end

% plot distribution of distances between cluster members
figure
histogram(shortesMatrix_distance(:))
xlabel('Length of path between two metabolites')
ylabel('Number of metabolite pairs')
title('Shortest distances between all cluster metabolites')
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
    [figureFolder 'fig_sup_histogram_cluster_met_pairs_distances_model2LIhost1LUbact_subprod_together.pdf'])

% save subset of pathes that belong to the cluster to file
% save('shortestMatrix_cluster_model2LIhost1LUbact_subprod.mat',...
%     'shortestMatrix_cluster', 'shortestMatrix_sources','rxnCompounds');

% filter shortest matrix to only connect substrates and products, and not
% substrates and substrates or products and products
shortestMatrix_distance_filtered = shortesMatrix_distance;

for i=1:size(shortesMatrix_distance)
    for j=i:size(shortesMatrix_distance)
        if ~(( ismember(shortestMatrix_sources{i}, potential_substrates_ids_all) &&...
               ismember(shortestMatrix_sources{j}, potential_products_ids_all) ) ||...
             ( ismember(shortestMatrix_sources{i}, potential_products_ids_all) &&...
               ismember(shortestMatrix_sources{j}, potential_substrates_ids_all)) )
                shortestMatrix_distance_filtered(i,j) = 0;
                shortestMatrix_distance_filtered(j,i) = 0;
        end
    end
end

%remove paths between non substrates and products
shortestMatrix_cluster_filtered = shortestMatrix_cluster;
for i=1:size(shortestMatrix_distance_filtered)
    for j=i:size(shortestMatrix_distance_filtered)
        if shortestMatrix_distance_filtered(i,j)==0
            shortestMatrix_cluster_filtered{i,j} = [];
            shortestMatrix_cluster_filtered{j,i} = [];
        end
    end
end     

% plot distribution of distances between cluster members
figure
plotdistances = triu(shortestMatrix_distance_filtered);
plotdistances  = plotdistances(:);
plotdistances(plotdistances==0)=[];
histogram(plotdistances)
xlabel('Length of path between two metabolites')
ylabel('Number of metabolite pairs')
title({'Shortest distances between', 'potential substrates and products'})
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
    [figureFolder,...
    'fig_sup_histogram_cluster_met_pairs_distances_model2LIhost1LUbact_subprod_filtered.pdf'])

% save subset of pathes that belong to the cluster to file
% save('shortestMatrix_cluster_model2LIhost1LUbact_subprod_filtered.mat',...
%     'shortestMatrix_cluster_filtered', 'shortestMatrix_sources','rxnCompounds',...
%     'shortestMatrix_distance_filtered',...
%     'potential_products_ids_all', 'potential_substrates_ids_all');
% 

% save path analysis results to file
% shortest distance
table_shortest_distance = array2table(shortestMatrix_distance_filtered,...
    'VariableNames', shortestMatrix_sources);
table_shortest_distance.Source = shortestMatrix_sources;
table_shortest_distance = table_shortest_distance(:,['Source'; shortestMatrix_sources]);
writetable(table_shortest_distance, [outputFolder, filesep, 'table_shortestMatrix_distance_filtered.csv']);

% rxnCompounds
table_rxnCompounds = cell2table(rxnCompounds,...
    'VariableNames', {'KEGGID'});
writetable(table_rxnCompounds, [outputFolder, filesep, 'table_rxnCompounds.csv']);

% all pathes
fid = fopen([outputFolder, filesep, 'table_shortestMatrix_paths_filtered.csv'],'w');
fprintf(fid, 'row,column,path\n');
for i=1:size(shortestMatrix_cluster_filtered,1)
    for j=1:size(shortestMatrix_cluster_filtered,2)
        if ~isempty(shortestMatrix_cluster_filtered{i,j})
            for k=1:length(shortestMatrix_cluster_filtered{i,j})
                curpath = shortestMatrix_cluster_filtered{i,j}{k};
                fprintf(fid, '%d,%d', i, j);
                for l=1:length(curpath)
                    fprintf(fid, ',%d', curpath(l));
                end
                fprintf(fid, '\n');
            end
        end
    end
end
    