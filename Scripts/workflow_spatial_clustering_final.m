%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cluster ions based on their GI profiles

% call script defining file dependenciesand global variables
addpath(genpath('.\'));
add_global_and_file_dependencies

% add_global_and_file_dependencies contains variable
% testing_mode_flag (=1)
% that determines that scripts runs only for annotated metabolites and not the
% full dataset (lines 72-75)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data requirements:
% 'metabolites_allions_combined_formulas_with_metabolite_filters.csv'
% 'metabolites_allions_combined_norm_intensity.csv'
% 'hmdbPTWtables.mat'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output: 
% Files:
% 'metabolites_allions_combined_formulas_with_metabolite_filters_spatial100clusters.csv'
% 'metabolites_allions_combined_formulas_with_metabolite_filters_spatial100clusters.csv'
% 'metabolites_allions_combined_formulas_with_metabolite_filters_spatial100clusters_HMDBsubclass.csv'
% 'metabolites_allions_combined_formulas_with_metabolite_filters_spatial100clusters_with_mean.csv'
% 'ptwenr_HMDBsuperclass_kmeans_clustering_spatial_all_hmdbv4.csv'
% Figures:
% 'fig_sup_kmeans_silhouette_evaluation.pdf'
% 'fig2d_kmeans_profiled_6_Mouse_Diet_GI_profiles_updated_filter_100repeat.pdf'
% 'fig2e_clustergram_spatial_clustering_hmdb_v4_SubClass';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read annotation from file
annotationTable = readtable([outputFolder,...inputFolder,...
    'metabolites_allions_combined_formulas_with_metabolite_filters_updated_filtering_0925.csv']);   
%readtable([inputFolder,...
%    'metabolites_allions_combined_formulas_with_metabolite_filters.csv']);
% read metabolite normalized data from file
metaboliteData = readtable([outputFolder,...inputFolder,...
    'metabolites_allions_combined_norm_intensity_with_CVR_300825.csv']);
%readtable([outputFolder,...
%    'metabolites_allions_combined_norm_intensity_with_CVR.csv']);
%metaboliteData = readtable([inputFolder,...
    %'metabolites_allions_combined_norm_intensity.csv']);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% get mz and RT of ions as numeric vector
joinedMzRT = table2array(metaboliteData(:,{'MZ','RT'}));
% updated filter of annotated metabolites
% use all ions, not only annotated
if testing_mode_flag
    binaryFilter = annotationTable.MetaboliteFilter;
else
    % all ions
    binaryFilter = ones(size(annotationTable,1),1);
end

% calculate mean and std of intensities in each tissue and condition
meanMatrix = zeros(size(combinedIntensitiesNorm,1),...
                   length(sampleType_unique)*length(sampleDiet_unique)*length(sampleTissue_unique));
stdMatrix = zeros(size(combinedIntensitiesNorm,1),...
                  length(sampleType_unique)*length(sampleDiet_unique)*length(sampleTissue_unique));
repMatrix = zeros(size(combinedIntensitiesNorm,1),...
                  length(sampleType_unique)*length(sampleDiet_unique)*length(sampleTissue_unique));              
fcidx = 1;
for d = 1:length(sampleTissue_unique)
    curtissueIDX = ismember(combinedTissues, sampleTissue_unique{d});
    for pert = 1:length(sampleType_unique)
        curpertIDX = ismember(combinedType, sampleType_unique{pert});
        for i = 1:length(sampleDiet_unique)
            curdietIDX = ismember(combinedDiet, sampleDiet_unique{i});
            % calculate mean and std per metabolite removing unmeasured but
            % leaving at least 3 replicates
            for met_i = 1:size(combinedIntensitiesNorm,1)
                curdata = combinedIntensitiesNorm(met_i,curtissueIDX &...
                                                                 curpertIDX & curdietIDX);
                % save number of nonzero replicates
                repMatrix(met_i, fcidx) = nnz(curdata>intensityNoise);
                if repMatrix(met_i, fcidx)>2
                    curdata = curdata(curdata>intensityNoise);
                else
                    curdata = sort(curdata, 'descend');
                    curdata = curdata(1:3);
                end
                meanMatrix(met_i,fcidx) = mean(curdata);
                stdMatrix(met_i,fcidx) = std(curdata);
            end
            fcidx = fcidx+1;
        end
    end
end
clear curcol curdata curdietIDX curpertIDX curtissueIDX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create vector of mean conditions
meanConditions = cell(size(meanMatrix,2),1);
idx = 1;
for d = 1:length(sampleTissue_unique)
    for pert = 1:length(sampleType_unique)
        for diet_i = 1:length(sampleDiet_unique)
            meanConditions{idx} = [sampleTissue_unique{d} ' ' ...
                                   sampleType_unique{pert} ' ' ...
                                   sampleDiet_unique{diet_i}];
            idx = idx+1;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select conditions for clustering
conditionsArray = cell(2, length(sampleDiet_unique)*length(sampleType_unique));
%conditionsArray = {'HFD', 'CTR', 'HFD', 'CTR';
%                   'GF', 'GF', 'DC', 'DC'};
i=1;
for type_i = 1:length(sampleType_unique)
    for diet_i = 1:length(sampleDiet_unique)
        conditionsArray{1,i} = sampleDiet_unique{diet_i};
        conditionsArray{2,i} = sampleType_unique{type_i};
        i = i+1;
    end
end

% find optimal value of K
numiter = 20;
kmeansEva_cell = cell(size(conditionsArray,2),numiter);
idx=1;
fprintf('Selecting optimal k for different conditions\n');
for cond_i = 3:size(conditionsArray,2)
    for iter_i = 1:numiter
        selectDiet = conditionsArray{1,cond_i};
        selectMouse = conditionsArray{2,cond_i};
        kmeanMatrix = meanMatrix(:,cellfun(@(x) contains(x,selectMouse) & ...
                                   contains(x,selectDiet),meanConditions));
        kmeanConditions = meanConditions(cellfun(@(x) contains(x,selectMouse) & ...
                                    contains(x,selectDiet),meanConditions));
        kmeanMatrix = kmeanMatrix(:,1:6);
        kmeanConditions = kmeanConditions(1:6);
    
        selected_mets = find(sum(kmeanMatrix~=intensityNoise,2)~=0 &...
                             binaryFilter);
    
        kmeans_ions = joinedMzRT(selected_mets);
        kmeans_means = zscore(kmeanMatrix(selected_mets,:),[],2);
        eva = evalclusters(kmeans_means,'kmeans','silhouette', 'KList', 1:30);
        kmeansEva_cell{idx, iter_i} = eva;
        fprintf('Selected optimal k for %s %s (%d of %d) iteration %d of %d\n', selectDiet, selectMouse, cond_i, size(conditionsArray,2), iter_i, numiter);
    end
    idx = idx+1;
end

% delete empty rows if any
empty_row =find(cellfun(@(x) isempty(x), kmeansEva_cell(:,1)),1);
kmeansEva_cell(empty_row:end,:) = [];

% make a table with all criteria and write to file
colLabels = strcat(conditionsArray(1,:),'-', ...
    conditionsArray(2,:));
if length(colLabels)>size(kmeansEva_cell,1)
    colLabels = colLabels(3:end);
end
criterion_mat = zeros(numiter*size(kmeansEva_cell,1), length(kmeansEva_cell{1,1}.CriterionValues));
criterion_cols = arrayfun(@(x) num2str(x), kmeansEva_cell{1,1}.InspectedK, 'unif', 0);
criterion_rows = cell(numiter*size(kmeansEva_cell,1),1);
idx=1;
for i=1:size(kmeansEva_cell,1)
    for iter_i=1:numiter
        criterion_mat(idx,:) = kmeansEva_cell{i, iter_i}.CriterionValues;
        criterion_rows{idx} = strcat(colLabels{i}, '_iter', num2str(iter_i));
        idx = idx+1;
    end
end
criterion_mat = array2table(criterion_mat,...
    'RowNames',criterion_rows, 'VariableNames', criterion_cols);
if nnz(binaryFilter==0)==0
    writetable(criterion_mat, [outputFolder, 'kmeans_silhouette_allmets_20iter_092025.csv'],...
    'WriteRowNames', 1);
else
     writetable(criterion_mat, [outputFolder, 'kmeans_silhouette_annotmets_20iter_092025.csv'],...
    'WriteRowNames', 1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mycolors = [ 0    0.4492    0.6953;...
%     0.8242    0.3750    0.1523;...
%     0.7969    0.8867    0.9375;...
%     0.9609    0.8711    0.8281];
mycolors = [ 0    0.4492    0.6953;...
    0.8242    0.3750    0.1523];
mylinestyles = {':', '-', '--'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot eval silhouette
conditionsArray_plot = conditionsArray(:,3:end);
figure
hold on
for i=1:size(kmeansEva_cell,1)
    if isequal(conditionsArray_plot{1,i}, sampleDiet_unique{1})
        curcolor = mycolors(1,:);
    else
        curcolor = mycolors(2,:);
    end
    if isequal(conditionsArray_plot{2,i}, sampleType_unique{1})
        continue; % skip CVR
        curstyle = mylinestyles{1};
    elseif isequal(conditionsArray_plot{2,i}, sampleType_unique{2})
        curstyle = mylinestyles{2};
    else
        curstyle = mylinestyles{3};
    end

    criterion_mat = zeros(numiter, length(kmeansEva_cell{i,1}.CriterionValues));
    for iter_i=1:numiter
        criterion_mat(iter_i,:) = kmeansEva_cell{i, iter_i}.CriterionValues;
    end
    criterion_mean = mean(criterion_mat);
    criterion_std = std(criterion_mat);

    % plot(kmeansEva_cell{i}.InspectedK,...
    %      criterion_mean,...
    %      'Color', curcolor, 'LineStyle', curstyle,  'LineWidth', 2)
    errorbar(kmeansEva_cell{i,1}.InspectedK,...
         criterion_mean,...
         criterion_std,...
         'Color', curcolor, 'LineStyle', curstyle,  'LineWidth', 2)


    % plot(kmeansEva_cell{i}.InspectedK,...
    %      kmeansEva_cell{i}.CriterionValues,...
    %      'Color', curcolor, 'LineStyle', curstyle,  'LineWidth', 2)
end
xlim([1 30])
set(gca, 'XTick', 1:30)
set(gca, 'XTickLabel', 1:30)
xlabel('Number of clusters')
ylabel(kmeansEva_cell{1}.CriterionName)
title('Kmeans cluster evaluation')

%legend(colLabels([2 3 5 6]), 'location', 'southeast')
legend(colLabels, 'location', 'southeast')
if nnz(binaryFilter==0)==0
    print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
            [figureFolder,...
            'fig_sup_kmeans_silhouette_evaluation_numiter20_allions_noCVR_0925.pdf'])
else
    print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
            [figureFolder,...
            'fig_sup_kmeans_silhouette_evaluation_numiter20_annotions_noCVR_0925.pdf'])
end    
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test cluster robustness
nrepeat = 100;
ncluster = 6;
kmeans_clusters_repeat_cell = cell(length(sampleDiet_unique)*...
                                   length(sampleType_unique),1);
cluster_occurance = zeros(ncluster,length(kmeans_clusters_repeat_cell));
idx = 1;
for diet_i = 1:length(sampleDiet_unique)
    for type_i = 1:length(sampleType_unique)
        selectDiet = sampleDiet_unique{diet_i};
        selectMouse = sampleType_unique{type_i};
        kmeanMatrix = meanMatrix(:,cellfun(@(x) contains(x,selectMouse) & contains(x,selectDiet),meanConditions));
        kmeanConditions = meanConditions(cellfun(@(x) contains(x,selectMouse) & contains(x,selectDiet),meanConditions));
        giConditions = cellfun(@(x) ~(contains(x,'Serum') | contains(x,'Liver')),...
                            kmeanConditions);
        kmeanMatrix = kmeanMatrix(:,giConditions);
        kmeanConditions = kmeanConditions(giConditions);
        % print tissue order to doublecheck
        fprintf('Conditions: %s\n', strjoin(kmeanConditions,';'));
        
        selected_mets = find(sum(kmeanMatrix~=intensityNoise,2)~=0 &...
                             binaryFilter);


        kmeans_ions = joinedMzRT(selected_mets);
        kmeans_means = zscore(kmeanMatrix(selected_mets,:),[],2);
        % prepare matrix to store repeat cluster results
        kmeans_clusters_repeat = zeros(size(kmeans_means,1),nrepeat);
        maxpeakrepeat = zeros(nrepeat, 6);
        for j=1:nrepeat
            kmeans_clusters = kmeans(kmeans_means, ncluster);
            % sort clusters by peak time point
            [~, maxpeak] = arrayfun(@(x) max(mean(kmeans_means(kmeans_clusters==x,:))), 1:ncluster);
            maxpeakrepeat(j,:) = maxpeak;
            kmeans_clusters_repeat(:,j) = kmeans_clusters;
        end
        % save clusters in the order of peak appearance
        kmeans_clusters_repeat_cell{idx} = zeros(size(kmeans_clusters_repeat));
        for j=1:nrepeat
            maxpeak = maxpeakrepeat(j,:);
            [~, plotptw_rows] = sort(maxpeak);
            for i=1:ncluster
                kmeans_clusters_repeat_cell{idx}(kmeans_clusters_repeat(:,j)==plotptw_rows(i),j) = i;
            end
        end 
        idx = idx+1;
    end
end
        
% determine clustering for each ion based on majority repeats
kmeans_clusters_spatial_conditions = zeros(length(binaryFilter), length(kmeans_clusters_repeat_cell));
idx=1;
for diet_i = 1:length(sampleDiet_unique)
    for type_i = 1:length(sampleType_unique)
        selectDiet = sampleDiet_unique{diet_i};
        selectMouse = sampleType_unique{type_i};
        % get the matrix to get the filtering
        kmeanMatrix = meanMatrix(:,cellfun(@(x) contains(x,selectMouse) & contains(x,selectDiet),meanConditions));
        kmeanConditions = meanConditions(cellfun(@(x) contains(x,selectMouse) & contains(x,selectDiet),meanConditions));
        giConditions = cellfun(@(x) ~(contains(x,'Serum') | contains(x,'Liver')),...
                            kmeanConditions);
        kmeanMatrix = kmeanMatrix(:,giConditions);
        kmeanConditions = kmeanConditions(giConditions);
        % print tissue order to doublecheck interpretation
        fprintf('Conditions: %s\n', strjoin(kmeanConditions,';'));
     
        
        selected_mets = find(sum(kmeanMatrix~=intensityNoise,2)~=0 &...
                             binaryFilter);
        % for each met, mode selects the most frequent value
        final_clusters = mode(kmeans_clusters_repeat_cell{idx},2);
        kmeans_clusters_spatial_conditions(selected_mets,idx) = final_clusters;
        
        % save repeated clustering results to file
        kmeans_clusters_spatial_conditions_table = zeros(length(binaryFilter),...
                                    size(kmeans_clusters_repeat_cell{idx},2));
        kmeans_clusters_spatial_conditions_table(selected_mets,:) = kmeans_clusters_repeat_cell{idx};
        kmeans_clusters_spatial_conditions_table = array2table(kmeans_clusters_spatial_conditions_table);
        if nnz(binaryFilter==0)==0
            writetable(kmeans_clusters_spatial_conditions_table,...
                [outputFolder, 'kmeans_clustering_allions_GIT100_',...
                selectDiet, '_', selectMouse, '.txt']);
        else
            writetable(kmeans_clusters_spatial_conditions_table,...
                [outputFolder, 'kmeans_clustering_annotions_GIT100_',...
                selectDiet, '_', selectMouse, '.txt']);
        end

        idx = idx+1;
    end
end

% add clustering to annotation table and save the results
annotationTable = [ annotationTable,...
                    array2table(kmeans_clusters_spatial_conditions,...
                                 'VariableNames',...
                                strcat('spatial_clust100_', conditionsArray(1,:), '_', ...
                                                         conditionsArray(2,:)))];
% add index to track alternative annotations
if nnz(binaryFilter==0)==0
     write(annotationTable,...
        [outputFolder,...
    'metabolites_allions_combined_formulas_with_metabolite_filters_spatial100clusters_0925.csv']);
else
         write(annotationTable,...
        [outputFolder,...
    'metabolites_annotions_combined_formulas_with_metabolite_filters_spatial100clusters_0925.csv']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot cluster centers for each condition
if nnz(binaryFilter==0)>0 % plot only for filtered metabolites
    idx=1;
    for diet_i = 1:length(sampleDiet_unique)
        for type_i = 1:length(sampleType_unique)
            selectDiet = sampleDiet_unique{diet_i};
            selectMouse = sampleType_unique{type_i};
            kmeanMatrix = meanMatrix(:,cellfun(@(x) contains(x,selectMouse) & contains(x,selectDiet),meanConditions));
            kmeanConditions = meanConditions(cellfun(@(x) contains(x,selectMouse) & contains(x,selectDiet),meanConditions));
            giConditions = cellfun(@(x) ~(contains(x,'Serum') | contains(x,'Liver')),...
                                kmeanConditions);
            kmeanMatrix = kmeanMatrix(:,giConditions);
            kmeanConditions = kmeanConditions(giConditions);
            % print tissue order to doublecheck interpretation
            fprintf('Conditions: %s\n', strjoin(kmeanConditions,';'));
            
            kmeans_means = zscore(kmeanMatrix,[],2);
           
            kmeans_clusters = kmeans_clusters_spatial_conditions(:,idx);
            figure
            for i=1:ncluster
                subplot(2,3,i)
                curmeans = kmeans_means(kmeans_clusters==i,:);
                hold on
                h = fill([1,3.5,3.5,1],[-3,-3,3,3],'yellow', 'LineStyle','none');
                h.FaceAlpha=0.3;
                h = fill([3.5,6,6,3.5],[-3,-3,3,3],'green', 'LineStyle','none');
                h.FaceAlpha=0.3;
    
                for j=1:size(curmeans,1)
                    plot((1:size(curmeans,2)), curmeans(j,:), '-','Color', [.5 .5 .5]);
                end
                plot((1:size(curmeans,2)),mean(curmeans), 'o-', 'Color', mycolors(1,:), 'LineWidth', 2);
                xlim([1 size(curmeans,2)])
                ylim([-3 3])
                title(sprintf('Cluster %d',i))
    
                set(gca, 'XTick', 1:size(curmeans,2))
                set(gca, 'XTickLabel', kmeanConditions)
            end
            sgtitle(['Metabolite profiles of ' selectMouse ' group under ' selectDiet ' diet'])
            orient landscape
           print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
               [figureFolder,...
               'fig2d_kmeans_profiled_6_' selectMouse '_' selectDiet '_GI_profiles_updated_filter_100repeat.pdf'])
            idx = idx+1;
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find serum and liver-specific metabolites
kmeans_clustersLiver = zeros(length(joinedMzRT),...
    length(sampleDiet_unique)*length(sampleType_unique));
kmeans_clustersSerum = zeros(length(joinedMzRT),...
    length(sampleDiet_unique)*length(sampleType_unique));

idx = 1;
for diet_i = 1:length(sampleDiet_unique)
    for type_i = 1:length(sampleType_unique)
        selectDiet = sampleDiet_unique{diet_i};
        selectMouse = sampleType_unique{type_i};
        kmeanMatrix = meanMatrix(:,cellfun(@(x) contains(x,selectMouse) & contains(x,selectDiet),meanConditions));
        kmeanConditions = meanConditions(cellfun(@(x) contains(x,selectMouse) & contains(x,selectDiet),meanConditions));
        kmeanidxLiver = cellfun(@(x) contains(x, 'Liver'),...
                                kmeanConditions);
        kmeanidxSerum = cellfun(@(x) contains(x,'Serum'),...
                                kmeanConditions);
       
        kmeanMatrixGI = kmeanMatrix(:,~(kmeanidxLiver | kmeanidxSerum));
        selected_metsGI = (sum(kmeanMatrixGI~=intensityNoise,2)~=0 &...
                             binaryFilter);
        selected_metsLiver = (sum(kmeanMatrix(:,kmeanidxLiver)~=intensityNoise,2)~=0 &...
                             binaryFilter);
        selected_metsSerum = (sum(kmeanMatrix(:,kmeanidxSerum)~=intensityNoise,2)~=0 &...
                             binaryFilter);
        % first indicate ions mesured in serum or liver
        kmeans_clustersLiver(selected_metsLiver, idx) = 1;
        kmeans_clustersSerum(selected_metsSerum, idx) = 1;
        % next indicate ions mesured in serum and liver
        kmeans_clustersLiver(selected_metsLiver & selected_metsSerum, idx) = 2;
        kmeans_clustersSerum(selected_metsLiver & selected_metsSerum, idx) = 2;
        % finally indicate ions mesured in serum or liver and GIT
        kmeans_clustersLiver(selected_metsLiver & selected_metsGI, idx) = 6;
        kmeans_clustersSerum(selected_metsSerum & selected_metsGI, idx) = 6;
        
        idx = idx+1;
    end
end

% add serum and liver clustering to annotation table and save the results
annotationTable = [ annotationTable,...
                    array2table(kmeans_clustersLiver,...
                                 'VariableNames',...
                                strcat('spatial_Liver_', conditionsArray(1,:), '_', ...
                                                         conditionsArray(2,:)))];
annotationTable = [ annotationTable,...
                    array2table(kmeans_clustersSerum,...
                                 'VariableNames',...
                                strcat('spatial_Serum_', conditionsArray(1,:), '_', ...
                                                         conditionsArray(2,:)))];
% add index to track alternative annotations
if nnz(binaryFilter==0)==0
    write(annotationTable,...
    [outputFolder,...
    'metabolites_allions_combined_formulas_with_metabolite_filters_spatial100clusters_0925.csv']);
else
     write(annotationTable,...
    [outputFolder,...
    'metabolites_annotions_combined_formulas_with_metabolite_filters_spatial100clusters_0925.csv']);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate number of filtered metabolites
annotationTable = readtable([outputFolder,...
    'metabolites_allions_combined_formulas_with_metabolite_filters_spatial100clusters_0925.csv']);

% select columns from annotationTable with cluster info
annColumns = annotationTable.Properties.VariableNames';
%nnColumnsTEST = annotationTableTEST.Properties.VariableNames';
clusterColumns = annColumns(cellfun(@(x) contains(x,'spatial_'), annColumns));
for i=1:length(clusterColumns)
    curcol = strsplit(clusterColumns{i}, '_');
    combinedTissues{i} = curcol{2};
    combinedType{i} = curcol{4};
    combinedDiet{i} = curcol{3};
end
sampleType_unique = unique(combinedType, 'stable');
sampleDiet_unique = unique(combinedDiet, 'stable');
sampleTissue_unique = unique(combinedTissues, 'stable');
sampleTissue_unique = {'Duo','Jej','Ile','Cec','Col','Fec'};

[~, clusteridx] = intersect(annColumns, clusterColumns);

test = table2array(annotationTable(:, clusteridx));
% GIT ions
fprintf('Ions with GIT clustering N=%d\n', nnz(sum(test(:,1:4),2)))
% GIT metabolites
fprintf('Metabolites with GIT clustering N=%d\n', nnz((sum(test(:,1:4),2)>0) &...
    (annotationTable.MetaboliteFilter>0)))
% ans Ions with GIT clustering N=14156
%     Metabolites with GIT clustering N=1859

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform pathway enrichment for clustering

%load('hmdbPTWtables.mat');
load('hmdbPTWmimedbPTWtables.mat');
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
    % prepare enrichment table
    kmeansTables_current = cell(length(clusterColumns),1);

    idx = 1;
    for col_i = 1:length(clusterColumns)
        selected_mets = find((annotationTable.(clusterColumns{col_i})~=0) &...
            (annotationTable.MetaboliteFilter~=0));
        
        kmeans_clusters = annotationTable.(clusterColumns{col_i});
        kmeans_clusters = kmeans_clusters(selected_mets);

        % perform enrichment for each cluster
        plotptw_rows = unique(kmeans_clusters);
        
        % pathway enrichment analysis for clusters
        kmeansTables = cell(length(ptwNames), length(plotptw_rows)*3+1);
        ttidx = 1;
        keggCompoundsID = annotationTable.CompoundID;
        keggCompoundsID = keggCompoundsID(selected_mets);
        ptwIonTable = hmdbPtwIonTable(:,selected_mets);
        
        for i = 1:length(plotptw_rows)
            tic
            labels = kmeans_clusters==plotptw_rows(i);

            [~, ~, enrichmentTable] = ptwEnrichmentMetabolitesIonsPTWTable(ptwIonTable,...
                                                    keggCompoundsID, labels, ptwNames);

            %compile a joint table
            for jj=1:length(enrichmentTable)
                 enrichmentTable{jj,8} = ['(',num2str(enrichmentTable{jj,4}), ',',num2str(enrichmentTable{jj,5}), ',',num2str(enrichmentTable{jj,6}), ',',num2str(enrichmentTable{jj,7}), ')' ];
            end
            enrichmentTable(:,4:7) = [];
            if ttidx == 1
                kmeansTables(:, ttidx) = enrichmentTable(:,3);
                ttidx = ttidx+1;
            end
            kmeansTables(:,ttidx:ttidx+2) = [enrichmentTable(:,1:2), enrichmentTable(:,4)];
            ttidx = ttidx+3;
        toc
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isempty(kmeansTables{1,1})
            kmeansTables(:,1) = enrichmentTable(:,3);
        end
        % save enrichment tables
        kmeansTables_current{idx} = kmeansTables;
        idx = idx+1;
    end
    switch hmdbtype
        case 1
            kmeansTables_hmdbClass = kmeansTables_current;
        case 2
            kmeansTables_hmdbSubClass = kmeansTables_current;
        case 3
            kmeansTables_hmdbSuperClass = kmeansTables_current;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add hmdb subclass to annotation
annotationHMDBsubclass = cell(size(annotationTable,1),1);
annotationHMDBsuperclass = cell(size(annotationTable,1),1);
annotationHMDBclass = cell(size(annotationTable,1),1);
for i=1:size(annotationTable,1)
    if nnz(hmdbPTWclassTable(:,i))
        % class
        curHMDB = hmdbPTWclassNames(hmdbPTWclassTable(:,i)~=0);
        annotationHMDBclass{i} =  strjoin(curHMDB, ';');
        % subclass
        curHMDB = hmdbPTWSubClassNames(hmdbPTWSubClassTable(:,i)~=0);
        annotationHMDBsubclass{i} =  strjoin(curHMDB, ';');
        % superclass
        curHMDB = hmdbPTWSuperClassNames(hmdbPTWSuperClassTable(:,i)~=0);
        annotationHMDBsuperclass{i} =  strjoin(curHMDB, ';');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add HMDB annotation to the table and save to file
annotationTable.HMDBSuperClass = annotationHMDBsuperclass;
annotationTable.HMDBclass = annotationHMDBclass;
annotationTable.HMDBSubClass = annotationHMDBsubclass;
writetable(annotationTable,...
    [outputFolder,...
    'metabolites_allions_combined_formulas_with_metabolite_filters_spatial100clusters_HMDBsubclass_0925.csv']);


annotationTableSpatialClusters = readtable(...
    [outputFolder,...
    'metabolites_allions_combined_formulas_with_metabolite_filters_spatial100clusters_HMDBsubclass_0925.csv']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot pathway enrichment results
plot_hmdb = 3;%1-Class 2-Superclass 3-SubClass
plotclusterFlag = 0;%0-all, 1-serum and liver, 6-only GI
plotCVRflag = 0; 

switch plot_hmdb
    case 1
        kmeansTables_cell = kmeansTables_hmdbClass;
        printFileName = 'clustergram_spatial_clustering_hmdb_v4_Class';
    case 2
        kmeansTables_cell = kmeansTables_hmdbSuperClass;
         printFileName = 'clustergram_spatial_clustering_hmdb_v4_SuperClass';
    case 3
        kmeansTables_cell = kmeansTables_hmdbSubClass;
         printFileName = 'clustergram_spatial_clustering_hmdb_v4_SubClass';
end

showclust = 4;
ncluster = 12;
numMet = zeros(size(kmeansTables_cell{1},1), ncluster*showclust);
groupSize = zeros(size(kmeansTables_cell{1},1), ncluster*showclust);
pMet = zeros(size(kmeansTables_cell{1},1), ncluster*showclust); 
fdrMet = zeros(size(kmeansTables_cell{1},1), ncluster*showclust);
idx = 1;
for j=1:length(kmeansTables_cell)
    ncluster = (size(kmeansTables_cell{j},2)-1)/3;
    for i=1:ncluster
        pMet(:,idx) = cell2mat(kmeansTables_cell{j}(:,2+(i-1)*3));
        fdrMet(:,idx) = cell2mat(kmeansTables_cell{j}(:,3+(i-1)*3));
        for k=1:size(kmeansTables_cell{j},1)
            curstat = kmeansTables_cell{j}{k,4+(i-1)*3};
            curstat = strsplit(curstat(2:end-1), ',');
            numMet(k,idx) = str2double(curstat{1});
            groupSize(k,idx) = str2double(curstat{2});
        end
        idx = idx+1;
    end                          
end


% GI clusters
ncluster_total = 6;
% Generate cluster names
[bb,aa] = ndgrid(sampleType_unique, sampleDiet_unique); 
colLabels = strcat(aa(:),'-', bb(:));
[bb,aa] = ndgrid(colLabels, arrayfun(@(x) num2str(x), 1:ncluster_total, 'unif', 0)); 
bb = bb';
aa = aa';
colLabels = strcat(aa(:),'-', bb(:));
clusterColLabels = colLabels;
% Liver clusters
ncluster_total = 3;
% Generate cluster names
[bb,aa] = ndgrid(sampleType_unique, sampleDiet_unique); 
colLabels = strcat(aa(:),'-', bb(:));
[bb,aa] = ndgrid(colLabels, arrayfun(@(x) num2str(x), 1:ncluster_total, 'unif', 0)); 
bb = bb';
aa = aa';
colLabels = strcat(aa(:),'-liver-', bb(:));
clusterColLabels = [clusterColLabels; colLabels];
% serum clusters
ncluster_total = 3;
% Generate cluster names
[bb,aa] = ndgrid(sampleType_unique, sampleDiet_unique); 
colLabels = strcat(aa(:),'-', bb(:));
[bb,aa] = ndgrid(colLabels, arrayfun(@(x) num2str(x), 1:ncluster_total, 'unif', 0)); 
bb = bb';
aa = aa';
colLabels = strcat(aa(:),'-serum-', bb(:));
clusterColLabels = [clusterColLabels; colLabels];


% leave only GF and DC GI tract
curselect = 1:length(clusterColLabels);%[7:18 25:36];%[1:12 19:30];%

numMet = numMet(:, curselect);
pMet = pMet(:, curselect); 
fdrMet = fdrMet(:,curselect);
clusterColLabels = clusterColLabels(curselect);


selectclusters = 1:length(clusterColLabels);

switch plotclusterFlag
    case 0
        % plot all clusters
        %selectclusters = 1:48;
        plotted_clusters = '_all_clusters_';
    case 1
        selectclusters = cellfun(@(x) contains(x, '-serum-') | ...
                            contains(x, '-liver-'),  reshape(clusterColLabels,1,[]));
        plotted_clusters = '_Liver_Serum_clusters_';
    case 6
        selectclusters = ~cellfun(@(x) contains(x, '-serum-') | ...
                            contains(x, '-liver-'),  reshape(clusterColLabels,1,[]));
        plotted_clusters = '_GIT_clusters_';
end
if ~plotCVRflag
    selectclusters = selectclusters & ~cellfun(@(x) contains(x, 'CVR'), reshape(clusterColLabels,1,[]));
end
% if both 2-serum and 2-liver cluster names exist, leave only one of them
% as the yare identical criteria
if sum(cellfun(@(x) contains(x, '2-serum-'),  clusterColLabels))>0 &&...
   sum(cellfun(@(x) contains(x, '2-liver-'),  clusterColLabels))>0 
    selectclusters = selectclusters & ~cellfun(@(x) contains(x, '2-serum-'), ...
        reshape(clusterColLabels,1,[]));
end


selectpathways = sum(fdrMet(:, selectclusters)<0.1,2)>0;
%displaymat = numMet(selectpathways,selectclusters);
%normalize number of metabolites per cluster size
displaymat = numMet(selectpathways,selectclusters)./...
             groupSize(selectpathways,selectclusters);


%displaymat = zscore(displaymat, dim, 2);
%normalize by max in each category
for i=1:size(displaymat,1)
    displaymat(i,:) = displaymat(i,:)/max(displaymat(i,:));
end

%displaymat = pMet(selectpathways,selectclusters);
%displaymat = fdrMet(selectpathways,:);
%displayvalues = '_nummet_';
displayvalues = '_fractmetperclustsize_';

%displaymat(displaymat>0.5)=1;
displaymat(displaymat>30)=30;
displaymat(isnan(displaymat))=0; % put 0 instead of NaN

if iscell(kmeansTables_cell{1}{1,1})
    rowLabels = cellfun(@(x) x{1}, kmeansTables_cell{1}(selectpathways,1), 'unif', 0);
else
    rowLabels = kmeansTables_cell{1}(selectpathways,1);
end

cgo = clustergram(displaymat, ...
            'RowLabels', rowLabels,...
            'ColumnLabels', clusterColLabels(selectclusters),...
            'ColumnPDist', 'correlation',...'euclidean',...
            'RowPDist', 'correlation',... 'euclidean',...
            'Symmetric', 0,...
            'DisplayRange', 3,...
            'Colormap', slanCM('greys'));%flipud(bone));


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

% Annotate clustergram with FDR values
%resort FDR values based on clustergram ordering
[~, ~, orderrow] = intersect(cgo.RowLabels, rowLabels, 'stable');
[~, ~, ordercol] = intersect(cgo.ColumnLabels, clusterColLabels(selectclusters), 'stable');
clusterFDR = fdrMet(selectpathways,selectclusters);
clusterFDR = clusterFDR(orderrow, ordercol); 
% round to two decimals
clusterFDR = round(clusterFDR*100)/100;

% adding pause to get the correct figure positions from here https://de.mathworks.com/matlabcentral/answers/454417-figure-position-property-returning-incorrect-values
pause(0.01); % correct the drawing process---
clustcoord = fig.Position;

% offsets in the next lines are pure heuristics to plot FDR values on top of the heatmap, 
% and might need to be adapted manually in trial and error if displayed wrongly
clustercoordstep_x = clustcoord(3)/(size(clusterFDR,1)-0.7);
clustercoordstep_y = clustcoord(4)/(size(clusterFDR,2)+1);
% add offset to get to the heatmap position
clustcoord(1) = clustcoord(1)-0.01;
clustcoord(2) = clustcoord(2)-0.04;

for i=1:size(clusterFDR,1)
    for j=1:size(clusterFDR,2)
    
        curcolor = 'green';
        if clusterFDR(i,j)<=0.1
            curcolor = 'white';
        
            annotation(gcf, 'textbox', [clustcoord(1)+(j-1)*clustercoordstep_y,...
                                        clustcoord(2)+(i-1)*clustercoordstep_x,...
                                        0.1 0.1],...
                        'String', {'*'},...
                        'LineStyle', 'none', 'Color', curcolor, 'FontSize', 20)
        end
    end
end

% add title to colobar
C = findall(gcf,'type','ColorBar');                         
C.Label.String = 'Relative number of metabolites';


if plotclusterFlag == 6
    print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
        [figureFolder, 'fig2e_', printFileName, plotted_clusters, displayvalues, '_0925.pdf'])
    saveas(fig, [figureFolder, 'fig2e_', printFileName, plotted_clusters, displayvalues, '_0925.png']);
else
    print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
        [figureFolder, 'fig_sup_', printFileName, plotted_clusters, displayvalues, '_0925.pdf'])
     saveas(fig, [figureFolder, 'fig_sup_', printFileName, plotted_clusters, displayvalues, '_0925.png'])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save hmdb pw enrichment to file
% create names for joint table columns
[bb,aa] = ndgrid({'pval', 'FDR', 'STATS'},clusterColLabels); 
colLabels = strcat(aa(:),'-', bb(:));
colLabels = cellfun(@(x) strrep(x, '-', '_'), colLabels, 'unif', 0)';

for hmdbtype = 1:3
    switch hmdbtype
        case 1
            kmeansTables_cell = kmeansTables_hmdbClass;
            curvarName = 'HMDBclass';
        case 2
            kmeansTables_cell = kmeansTables_hmdbSubClass;
            curvarName = 'HMDBsubclass';
        case 3
            kmeansTables_cell = kmeansTables_hmdbSuperClass;
            curvarName = 'HMDBsuperclass';
    end
    % combine all tables together
    kmeansTables = kmeansTables_cell{1};
    for i=2:length(kmeansTables_cell)
        kmeansTables = [kmeansTables kmeansTables_cell{i}(:,2:end)];
    end
    kmeansTables_table = cell2table(kmeansTables, 'VariableNames',...
                                [curvarName colLabels]);
    writetable(kmeansTables_table,...
               [outputFolder,...
               'ptwenr_', curvarName, '_kmeans_clustering_spatial_all_hmdbv4_0925.csv']);
end
clear kmeansTables kmeansTables_cell kmeansTables_current kmeansTables_table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read hmdb pw enrichment to file
% create names for joint table columns
for hmdbtype = 1:3
    switch hmdbtype
        case 1
            curvarName = 'HMDBclass';
        case 2
            curvarName = 'HMDBsubclass';
        case 3
            curvarName = 'HMDBsuperclass';
    end
    kmeansTables_table = readtable( [outputFolder,...
               'ptwenr_', curvarName, '_kmeans_clustering_spatial_all_hmdbv4_0925.csv']);
 
    switch hmdbtype
        case 1
            kmeansTables_hmdbClass = cell(length(clusterColumns),1);
            idx=1;
            for i=2:3:length(kmeansTables_table.Properties.VariableNames)-1
                kmeansTables_hmdbClass{idx} = [table2cell(kmeansTables_table(:,1)),...
                    table2cell(kmeansTables_table(:,i:i+2))];
                idx=idx+1;
            end
        case 2
            kmeansTables_hmdbSubClass = cell(length(clusterColumns),1);
            idx=1;
            for i=2:3:length(kmeansTables_table.Properties.VariableNames)-1
                kmeansTables_hmdbSubClass{idx} = [table2cell(kmeansTables_table(:,1)),...
                    table2cell(kmeansTables_table(:,i:i+2))];
                idx=idx+1;
            end
        case 3
            kmeansTables_hmdbSuperClass = cell(length(clusterColumns),1);
            idx=1;
            for i=2:3:length(kmeansTables_table.Properties.VariableNames)-1
                kmeansTables_hmdbSuperClass{idx} = [table2cell(kmeansTables_table(:,1)),...
                    table2cell(kmeansTables_table(:,i:i+2))];
                idx=idx+1;
            end
    end
end
clear kmeansTables_table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save mean values per tissue to file
annotationTable_withmean = [annotationTable,...
    array2table(meanMatrix, 'VariableNames', ...
    cellfun(@(x) strcat('mean_', strrep(x,' ', '_')), meanConditions, 'unif', 0))];
% write(annotationTable_withmean,...
%     [outputFolder,...
%     'metabolites_allions_combined_formulas_with_metabolite_filters_spatial100clusters_with_mean.csv']);
write(annotationTable_withmean,...
    [outputFolder,...
    'metabolites_allions_combined_formulas_with_metabolite_filters_spatial100clusters_with_mean_with_CVR_0925.csv']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
