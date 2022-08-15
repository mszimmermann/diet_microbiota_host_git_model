% addpath('C:\Users\mazimmer\Documents\Documents\MATLAB\GutModels\Analysis_final')
% 
% inputFolder = '.\metabolomics\ProcessedData\';
% outputFolder = '.\ProccesedData\';
% figureFolder = '.\Figures_final\';

% call script defining file dependenciesand global variables
addpath(genpath('.\'));
add_global_and_file_dependencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data requirements:
% 'metabolites_allions_combined_formulas_with_metabolite_filters_spatial100clusters_with_mean.csv'
% 'metabolites_allions_combined_norm_intensity.csv'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output: 
% Files:
% 'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_allions.csv'
% 'model_results_SMOOTH_normbyabsmax_2LIcoefHost1LIcoefbact_allions.csv'
% 'model_results_SMOOTH_normbyabsmax_ONLYMETCOEF_2LIcoefHost1LIcoefbact_allions.csv'
% 'model_results_SMOOTH_normbyabsmax_reciprocal_problem_allions.csv'
% Figures:
% 'Fig4a_histogram_corr_model_2LIcoefHost1LIcoefbact_reversedata_annotated'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% try to cluster ions based on their GI profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read annotation from file
annotationTableSpatialClusters = readtable([inputFolder ...
    'metabolites_allions_combined_formulas_with_metabolite_filters_spatial100clusters_with_mean.csv']);
% read metabolite normalized data from file
metaboliteData = readtable([inputFolder ...
    'metabolites_allions_combined_norm_intensity.csv']);
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
% get clustering info
% select columns from annotationTable with cluster info
annColumns = annotationTableSpatialClusters.Properties.VariableNames;
clusterColumns = annColumns(cellfun(@(x) contains(x,'spatial_clust100'), annColumns));

[~, clusteridx] = intersect(annColumns, clusterColumns, 'stable');

spatialClusters = table2array(annotationTableSpatialClusters(:, clusteridx));

meanConditions = annColumns(cellfun(@(x) contains(x,'mean_'), annColumns));
[~, meanidx] = intersect(annColumns, meanConditions, 'stable');

meanMatrix = table2array(annotationTableSpatialClusters(:, meanidx));

mycolors = [0 115 178;... %dark blue
            204 227 240;...%light blue
            211 96 39;... %dark orange
            246 223 212]/256;%light orange
[bb,aa] = ndgrid(sampleType_unique, sampleDiet_unique); 
condLabels = strcat(aa(:),'-', bb(:));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% selected_mets are metabolites detected along the GI tract

% define which metabolites to model
%%%% all metabolites
% selected_mets = 1:size(meanMatrix,1);
%%%% metabolites detected in the GIT and annotated
selected_mets = find((sum(spatialClusters,2)>0) &...
                     (annotationTableSpatialClusters.MetaboliteFilter==1));
%%%% testing mode: hundred of annotated metabolites
selected_mets = selected_mets(1:100);

% define tissue order
sampleTissue_unique = {'SI1', 'SI2', 'SI3', 'Cecum', 'Colon', 'Feces'};
nrand=100;
% call calculateAmatrix_final with zero matrices of correct size
% to get coefvalues output, which is used to allocate variables further
[~, coefvalues] = calculateAmatrix_final(zeros(length(condLabels), length(sampleTissue_unique)));

% allocate variables for model prediction results
x_met_mean = zeros(length(coefvalues),length(selected_mets));
x_met_std = zeros(length(coefvalues),length(selected_mets));
x_met_smooth = zeros(length(coefvalues),length(selected_mets));

% allocate variables for correlation calculated from values restored with
% reverse problem with original and randomly shuffled parameters
x_data_corr = zeros(length(selected_mets),1);
x_resid = zeros(length(selected_mets),1);
x_data_corr_shuffled = zeros(length(selected_mets),1);

% calculate reverse A matrix with zero matrix of the right size to allocate
% variables according to Ra size
[Ra] = calculateRAmatrix_final(zeros(size(coefvalues)));
x_rdata = zeros(size(Ra,2),length(selected_mets));
x_rdata_shuffled = zeros(size(Ra,2),length(selected_mets));

minval = 0.01;

for met_i = 1:length(selected_mets)
   
    if sum(spatialClusters(met_i,:))>0
        cmpd_interest_idx = selected_mets(met_i);
        % calculate mean profiles            
        idx=1;
        kmeanMatrix_joint = zeros(4,6);
        for diet_i = 1:length(sampleDiet_unique)
            for type_i = 1:length(sampleType_unique)
                selectDiet = sampleDiet_unique{diet_i};
                selectMouse = sampleType_unique{type_i};
                kmeanMatrix = meanMatrix(:,cellfun(@(x) contains(x,selectMouse) & contains(x,selectDiet),meanConditions));
                kmeanMatrix = kmeanMatrix(cmpd_interest_idx,1:6);
                kmeanMatrix_joint(idx,:) = kmeanMatrix;
                idx = idx+1;
            end
        end
        % normalize by max intensity
        kmeanMatrix_joint = kmeanMatrix_joint./max(max(kmeanMatrix_joint));
        % replace small values with noise
        kmeanMatrix_joint_orig = kmeanMatrix_joint;
        kmeanMatrix_joint(kmeanMatrix_joint<minval)=...
        kmeanMatrix_joint(kmeanMatrix_joint<minval)+rand(nnz(kmeanMatrix_joint<minval),1)*minval;

        % calculate fluxes for average profile
        Aorig = calculateAmatrix_final(kmeanMatrix_joint_orig);
        [A,coefvalues] = calculateAmatrix_final(kmeanMatrix_joint);
        A(Aorig(:,1)==0,:)=[];
        b = zeros(size(A,1),1);
        options = optimoptions(@lsqlin,'Display', 'off');

        xlowerlim = zeros(size(A,2),1);
        xlowerlim(2:end)=-Inf;

        %%%%%%%%%%%%%%%%%%%%%%%%%
        %set bact to zero
        xupperlim = Inf*ones(size(A,2),1);
        %%%%%%%%%%%%%%%%%%%%%%%%%

        [x, xres] = lsqlin(A,b,[],[],[],[],xlowerlim, xupperlim, [], options);
        x_met_smooth(:, met_i) = x;
        x_resid(met_i) = xres;


    [Ra,rb] = calculateRAmatrix_final(x*1000);
    dataR = lsqlin(Ra,rb,[],[],[],[],zeros(1,size(Ra,2)), [], [], options);
    % save reconstructed data
    x_rdata(:,met_i) = dataR;

    %rsshape to matrix form
    dataR = reshape(dataR,[],4)';

    x_data_corr(met_i) = corr(kmeanMatrix_joint_orig(:), dataR(:));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate reverse problem for shuffled x
    x_shuffled = x(randperm(length(x)));
    [Ra_shuffled,rb_shuffled] = calculateRAmatrix_final(x_shuffled*1000);
    dataR_shuffled = lsqlin(Ra_shuffled,rb_shuffled,[],[],[],[],...
                            zeros(1,size(Ra_shuffled,2)), [], [], options);
    % save reconstructed data
    x_rdata_shuffled(:,met_i) = dataR_shuffled;

    %rsshape to matrix form
    dataR_shuffled = reshape(dataR_shuffled,[],4)';

    x_data_corr_shuffled(met_i) = corr(kmeanMatrix_joint_orig(:), dataR_shuffled(:));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % calculate joint matrix per mouse
        kmeanMatrix_joint_mouse = cell(5,1);
        for mouse_i = 1:5
            kmeanMatrix_joint = zeros(4,6);
            kmeanMatrix_joint_samples = cell(4,6);
            idx=1;
            for diet_i = 1:length(sampleDiet_unique)
                for type_i = 1:length(sampleType_unique)
                    for tissue_i = 1:length(sampleTissue_unique)
                        selectDiet = sampleDiet_unique{diet_i};
                        selectType = sampleType_unique{type_i};
                        selectTissue = sampleTissue_unique{tissue_i};
                        kmeanMatrix = combinedIntensitiesNorm(cmpd_interest_idx,...
                                                cellfun(@(x) contains(x,selectType),combinedType) & ...
                                                cellfun(@(x) contains(x,selectDiet),combinedDiet) & ...
                                                cellfun(@(x) contains(x,selectTissue),combinedTissues));
                        if length(kmeanMatrix)<mouse_i
                            kmeanMatrix_joint(idx,tissue_i) = mean(kmeanMatrix);
                        else
                            kmeanMatrix_joint(idx,tissue_i) = kmeanMatrix(mouse_i);
                        end
                        kmeanMatrix_joint_samples{idx,tissue_i} = ['P_' lower(strrep(selectTissue,'DC','feces')) ...
                                                                   '_' lower(strrep(selectType,'DC','cv'))...
                                                                   '_' lower(selectDiet)];
                    end
                    idx = idx+1;
                end
            end

            kmeanMatrix_joint = kmeanMatrix_joint./max(max(kmeanMatrix_joint));
            kmeanMatrix_joint_mouse{mouse_i} = kmeanMatrix_joint;
        end

        x_rand = zeros(length(coefvalues), nrand);
        for i=1:nrand
            curmice = randi(5,1,4);
            kmeanMatrix_joint = zeros(4,6);
            for j=1:length(curmice)
                kmeanMatrix_joint(j,:) = kmeanMatrix_joint_mouse{curmice(j)}(j,:);
            end
            % replace small values with noise
            kmeanMatrix_joint_orig = kmeanMatrix_joint;
            kmeanMatrix_joint(kmeanMatrix_joint<minval)=...
                kmeanMatrix_joint(kmeanMatrix_joint<minval)+rand(nnz(kmeanMatrix_joint<minval),1)*minval;


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Aorig = calculateAmatrix_final(kmeanMatrix_joint_orig);
            [A,coefvalues] = calculateAmatrix_final(kmeanMatrix_joint);
            A(Aorig(:,1)==0,:)=[];
            b = zeros(size(A,1),1);
            options = optimoptions(@lsqlin,'Display', 'off');

            xlowerlim = zeros(size(A,2),1);
            xlowerlim(2:end)=-Inf;

            %%%%%%%%%%%%%%%%%%%%%%%%%
            %set bact to zero
            xupperlim = Inf*ones(size(A,2),1);
            %%%%%%%%%%%%%%%%%%%%%%%%%

            if size(A,1)>0
                x = lsqlin(A,b,[],[],[],[],xlowerlim, xupperlim, [], options);
                x_rand(:, i) = x;
            end

        end
        %  remove all zero solutions
        x_rand(:,sum(x_rand==0,1)==size(x,1))=[];
        if size(x_rand,1)>0
            x_met_mean(:,met_i) = mean(x_rand,2);
            x_met_std(:,met_i) = std(x_rand,[],2);
        end
    end
    % display status every 100 metabolites
    if mod(met_i,100)==0
        fprintf('Calculated parameters for %d metabolites\n',i);
    end
end

% get the original metabolomics data as vectors to calculate correlations
% between original and restored measurements
kmean_vector_joint_orig = zeros(size(Ra,2),length(selected_mets));
kmean_vector_joint = zeros(size(Ra,2),length(selected_mets));
for met_i = 1:length(selected_mets)
   
    cmpd_interest_idx = selected_mets(met_i);
    % calculate mean profiles            
    idx=1;
    kmeanMatrix_joint = zeros(4,6);
    for diet_i = 1:length(sampleDiet_unique)
        for type_i = 1:length(sampleType_unique)
            selectDiet = sampleDiet_unique{diet_i};
            selectMouse = sampleType_unique{type_i};
            kmeanMatrix = meanMatrix(:,cellfun(@(x) contains(x,selectMouse) & contains(x,selectDiet),meanConditions));
            % get serum and liver data
            kmeanMatrix = kmeanMatrix(cmpd_interest_idx,1:6);
            kmeanMatrix_joint(idx,:) = kmeanMatrix;
            idx = idx+1;
        end
    end
    % normalize by max intensity
    kmeanMatrix_joint = kmeanMatrix_joint./max(max(kmeanMatrix_joint));
    % replace small values with noise
    kmeanMatrix_joint_orig = kmeanMatrix_joint;
    kmeanMatrix_joint(kmeanMatrix_joint<minval)=...
    kmeanMatrix_joint(kmeanMatrix_joint<minval)+rand(nnz(kmeanMatrix_joint<minval),1)*minval;
    
    kmean_vector_joint_orig(:, met_i) = kmeanMatrix_joint_orig(:);
    kmean_vector_joint(:, met_i) = kmeanMatrix_joint(:);
end

% calculate Spearman correlations
x_corr_Spearman = zeros(size(x_data_corr));
corr_Mouse = zeros(length(x_data_corr),4);
corr_Mouse_Spearman = zeros(length(x_data_corr),4);

for i = 1:length(x_data_corr)
    dataR = reshape(x_rdata(:,i),[],4)';
    dataOrig = reshape(kmean_vector_joint_orig(:,i),4,[]);

    x_corr_Spearman(i) = corr(dataOrig(:), dataR(:), 'type', 'Spearman');
    corr_Mouse(i,:) = diag(corr(dataOrig', dataR'));
    corr_Mouse_Spearman(i,:) = diag(corr(dataOrig', dataR', 'type', 'Spearman'));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save model results to file - raw coefficients
% add gut filter to file - flag indicating whether metabolites were
% detected in the GIT
gut_filter = (sum(spatialClusters,2)>0);

fid = fopen([outputFolder ...
            'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_allions.csv'], 'w');
fprintf(fid, 'MZ\tRT\tCompoundID\tCompoundName\tMetaboliteFilter\tSumGITclusters\tReciprocalCorr');
for i=1:length(coefvalues)
    fprintf(fid, '\t%s', coefvalues{i});
end
fprintf(fid, '\n');
for i=1:size(x_met_smooth,2)
    fprintf(fid, '%.3f\t%3f', annotationTableSpatialClusters.MZ(i),...
        annotationTableSpatialClusters.RT(i));
    fprintf(fid, '\t%s\t%s\t%d', annotationTableSpatialClusters.CompoundID{i},...
        annotationTableSpatialClusters.CompoundName{i},...
        annotationTableSpatialClusters.MetaboliteFilter(i));
    fprintf(fid, '\t%d', gut_filter(i));
    fprintf(fid, '\t%.3f', x_data_corr(i));
    for j=1:size(x_met_smooth,1)
        fprintf(fid, '\t%e', x_met_smooth(j,i));
    end
    fprintf(fid, '\n');
end
fclose(fid);

% save model results to file - normalized by max coefficient
fid = fopen([outputFolder...
    'model_results_SMOOTH_normbyabsmax_2LIcoefHost1LIcoefbact_allions.csv'], 'w');
fprintf(fid, 'MZ\tRT\tCompoundID\tCompoundName\tMetaboliteFilter\tSumGITclusters\tReciprocalCorr');
for i=1:length(coefvalues)
    fprintf(fid, '\t%s', coefvalues{i});
end
fprintf(fid, '\n');
for i=1:size(x_met_smooth,2)
    fprintf(fid, '%.3f\t%3f', annotationTableSpatialClusters.MZ(i),...
        annotationTableSpatialClusters.RT(i));
    fprintf(fid, '\t%s\t%s\t%d', annotationTableSpatialClusters.CompoundID{i},...
        annotationTableSpatialClusters.CompoundName{i},...
        annotationTableSpatialClusters.MetaboliteFilter(i));
    fprintf(fid, '\t%d', gut_filter(i));
    fprintf(fid, '\t%.3f', x_data_corr(i));
    for j=1:size(x_met_smooth,1)
        fprintf(fid, '\t%.3f', x_met_smooth(j,i)./max(abs(x_met_smooth(:,i))));
    end
    fprintf(fid, '\n');
end
fclose(fid);

% save model results to file - metabolism_coefficients_only
fid = fopen([outputFolder...
    'model_results_SMOOTH_normbyabsmax_ONLYMETCOEF_2LIcoefHost1LIcoefbact_allions.csv'], 'w');
fprintf(fid, 'MZ\tRT\tCompoundID\tCompoundName\tMetaboliteFilter\tSumGITclusters\tReciprocalCorr');
for i=2:length(coefvalues)
    fprintf(fid, '\t%s', coefvalues{i});
end
fprintf(fid, '\n');
for i=1:size(x_met_smooth,2)
    fprintf(fid, '%.3f\t%3f', annotationTableSpatialClusters.MZ(i),...
        annotationTableSpatialClusters.RT(i));
    fprintf(fid, '\t%s\t%s\t%d', annotationTableSpatialClusters.CompoundID{i},...
        annotationTableSpatialClusters.CompoundName{i},...
        annotationTableSpatialClusters.MetaboliteFilter(i));
    fprintf(fid, '\t%d', gut_filter(i));
    fprintf(fid, '\t%.3f', x_data_corr(i));
    for j=2:size(x_met_smooth,1)
        fprintf(fid, '\t%.3f', x_met_smooth(j,i)./max(abs(x_met_smooth(2:end,i))));
    end
    fprintf(fid, '\n');
end
fclose(fid);
    
% save model results to file - reciprocal data restoration
fid = fopen([outputFolder...
    'model_results_SMOOTH_normbyabsmax_reciprocal_problem_allions.csv'], 'w');
columnNames = cell(size(x_rdata,1),1);
idx=1;
for diet_i = 1:length(sampleDiet_unique)
    for type_i = 1:length(sampleType_unique)
        for tiss_i = 1:length(sampleTissue_order)-2
                selectDiet = sampleDiet_unique{diet_i};
                selectMouse = sampleType_unique{type_i};
                selectTissue = sampleTissue_order{tiss_i};
                columnNames{idx} = [selectDiet '_' selectMouse '_' selectTissue];
                idx = idx+1;
        end
    end
end
         
fprintf(fid, 'MZ\tRT\tCompoundID\tCompoundName\tMetaboliteFilter\tSumGITclusters\tReciprocalCorr\tRandomCorr');
for i=1:length(columnNames)
    fprintf(fid, '\t%s', columnNames{i});
end
for i=1:length(columnNames)
    fprintf(fid, '\tRecip_%s', columnNames{i});
end
for i=1:length(columnNames)
    fprintf(fid, '\tRandom_%s', columnNames{i});
end
fprintf(fid, '\n');
for i=1:size(x_rdata,2)
    fprintf(fid, '%.3f\t%3f', annotationTableSpatialClusters.MZ(i),...
        annotationTableSpatialClusters.RT(i));
    fprintf(fid, '\t%s\t%s\t%d', annotationTableSpatialClusters.CompoundID{i},...
        annotationTableSpatialClusters.CompoundName{i},...
        annotationTableSpatialClusters.MetaboliteFilter(i));
    fprintf(fid, '\t%d', gut_filter(i));
    fprintf(fid, '\t%.3f', x_data_corr(i));
    fprintf(fid, '\t%.3f', x_data_corr_shuffled(i));
    for j=1:size(kmean_vector_joint_orig(:,i),1)
        fprintf(fid, '\t%.3f', kmean_vector_joint_orig(j,i));
    end
    for j=1:size(x_rdata(:,i),1)
        fprintf(fid, '\t%.3f', x_rdata(j,i));
    end
    for j=1:size(x_rdata_shuffled(:,i),1)
        fprintf(fid, '\t%.3f', x_rdata_shuffled(j,i));
    end
    fprintf(fid, '\n');
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot correlation of restored and original data
% %metabolites detected in the GIT and annotated
% met_filter = ((annotationTableSpatialClusters.MetaboliteFilter==1) &...
%               (sum(spatialClusters,2)>0));
% %testing mode: all calculated
met_filter = ones(size(x_data_corr,1),1);
% calculate differentce in corr distrbutions
p_corr_diff = ranksum(x_data_corr_shuffled(met_filter==1),...
    x_data_corr(met_filter==1));
figure
% perason corr all 
h = histogram(x_data_corr_shuffled(met_filter==1),100);
hold on
histogram(x_data_corr(met_filter==1),100)
xlim([-1 1])
axis square
xlabel('Pearson correlation between metabolomics data and model estimate')
ylabel('Number of ions')
title('All data together')
orient landscape
legend({'Random coefficients', 'Model coefficients'},...
        'Location', 'NorthWest')
    
% compare distributions of correlations
pval = signrank(x_data_corr_shuffled(met_filter==1),...
                x_data_corr(met_filter==1));
% print Wilcoxon signed rank test o-value on the plot
text(-0.9, 0.7*max(h.Values), sprintf('signrank p = %.2e', pval))
% save figure to file
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
    [figureFolder,...
    'Fig4a_histogram_corr_model_2LIcoefHost1LIcoefbact_reversedata_annotated'])


