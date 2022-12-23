% Fit steady state fluxes to GIT metabolites

% call script defining file dependenciesand global variables
addpath(genpath('.\'));
add_global_and_file_dependencies
% add_global_and_file_dependencies contains variable
% testing_mode_flag (=1)
% that determines that scripts runs only for 100 metabolites and not the
% full dataset (lines 100-105)

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
sampleType_unique = {'CV', 'GF'};
sampleDiet_unique = {'Chow1', 'Chow2'};

meanConditions = [cellfun(@(x) strrep(x, 'Chow', 'Chow1'), meanData_columns, 'unif', 0),...
                  cellfun(@(x) strrep(x, 'Chow', 'Chow2'), meanData_columns, 'unif', 0)];
              
meanMatrix = vertcat(meanData_cell{:});
% duplicate since we have only one diet
meanMatrix = [meanMatrix meanMatrix];
meanMatrix_mets = horzcat(meanMets_cell{:})';

mycolors = [0 115 178;... %dark blue
            204 227 240;...%light blue
            211 96 39;... %dark orange
            246 223 212]/256;%light orange
[bb,aa] = ndgrid(sampleType_unique, sampleDiet_unique); 
condLabels = strcat(aa(:),'-', bb(:));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% selected_mets are metabolites detected along the GI tract
selected_mets = 1:length(meanMatrix_mets);
% define tissue order
sampleTissue_unique = {'SI', 'SII', 'SIII', 'Cecum', 'Colon', 'Feces'};
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

fprintf('Starting parameter estimation for %d metabolites\n',length(selected_mets));

for met_i = 1:length(selected_mets)
   
        cmpd_interest_idx = selected_mets(met_i);
        % calculate mean profiles            
        idx=1;
        kmeanMatrix_joint = zeros(4,6);
        for diet_i = 1:length(sampleDiet_unique)
            for type_i = 1:length(sampleType_unique)
                selectDiet = sampleDiet_unique{diet_i};
                selectMouse = sampleType_unique{type_i};
                kmeanMatrix = meanMatrix(:,cellfun(@(x) contains(x,selectMouse) & ...
                                        contains(x,selectDiet),meanConditions));
                kmeanMatrix = kmeanMatrix(cmpd_interest_idx,1:6);
                kmeanMatrix_joint(idx,:) = kmeanMatrix;
                idx = idx+1;
            end
        end
        % normalize by max intensity
        kmeanMatrix_joint = kmeanMatrix_joint./max(max(kmeanMatrix_joint));
        kmeanMatrix_joint(isnan(kmeanMatrix_joint))=0;
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
        % check if A is empty an if yes move to the next metabolite
        if isempty(A)
            continue;
        end

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
            'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_allions_with_CVR.csv'], 'w');
             %'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_allions.csv'], 'w');
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
    'model_results_SMOOTH_normbyabsmax_2LIcoefHost1LIcoefbact_allions_with_CVR.csv'], 'w');
    %'model_results_SMOOTH_normbyabsmax_2LIcoefHost1LIcoefbact_allions.csv'], 'w');
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
    'model_results_SMOOTH_normbyabsmax_ONLYMETCOEF_2LIcoefHost1LIcoefbact_allions_with_CVR.csv'], 'w');
    %'model_results_SMOOTH_normbyabsmax_ONLYMETCOEF_2LIcoefHost1LIcoefbact_allions.csv'], 'w');
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
    'model_results_SMOOTH_normbyabsmax_reciprocal_problem_allions_with_CVR.csv'], 'w');
    %'model_results_SMOOTH_normbyabsmax_reciprocal_problem_allions.csv'], 'w');
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
% %testing mode: all calculated
if testing_mode_flag
    met_filter = ones(size(x_data_corr,1),1);
else
    % %metabolites detected in the GIT and annotated
    met_filter = ((annotationTableSpatialClusters.MetaboliteFilter==1) &...
               (sum(spatialClusters,2)>0));
end

% calculate differentce in corr distrbutions
p_corr_diff = ranksum(x_data_corr_shuffled(met_filter==1),...
    x_data_corr(met_filter==1));
figure
% pearson corr all 
h = histogram(x_data_corr_shuffled(met_filter==1),100);
hold on
histogram(x_data_corr(met_filter==1),100)
xlim([-1 1])
axis square
xlabel('Pearson correlation between metabolomics data and model estimate')
ylabel('Number of ions')
title('All data together')
orient landscape
% print line for PCC=0.7
plot([0.7, 0.7], [0, max(h.Values)], 'k--')

legend({'Random coefficients', 'Model coefficients', 'PCC=0.7'},...
        'Location', 'NorthWest')
    
% compare distributions of correlations
pval = signrank(x_data_corr_shuffled(met_filter==1),...
                x_data_corr(met_filter==1));
% print Wilcoxon signed rank test o-value on the plot
text(-0.9, 0.7*max(h.Values), sprintf('signrank p = %.2e', pval))
% save figure to file
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
    [figureFolder,...
    'Fig4a_histogram_corr_model_2LIcoefHost1LIcoefbact_reversedata_annotated_CVR'])
       %'Fig4a_histogram_corr_model_2LIcoefHost1LIcoefbact_reversedata_annotated'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot data and restored intensities and model coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define colora and GIT section names for plotting
mycolors = [0 115 178;... %dark blue
            204 227 240;...%light blue
            211 96 39;... %dark orange
            246 223 212]/256;%light orange
git_labels = {'Du', 'Je', 'Il', 'Cec', 'Col', 'Fec'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting file name
fileNameprofiles = 'fig_profiles_figureselected_modelSMOOTH_2LIcoefHost_1LIbact_drugs.ps';
          
curData_cols = reshape(meanConditions, 6, 4);
compoundsInterest = 1:length(meanMatrix_mets);

fig = figure('units','normalized','outerposition',[0 0 1 1]);

fprintf('Plotting profiles for %d metabolites\n',length(compoundsInterest));
for cpdix=1:length(compoundsInterest)
    
    testidx = compoundsInterest(cpdix);

    spx=1;
    spy=3;
    spidx = 1;
    coloridx = 1;
    idx=1;
    curmat = zeros(4,6);
    %cur_data = reshape(kmean_vector_joint_orig(:,testidx)', 4, 6)';
    cur_data = reshape(meanMatrix(testidx,:), 6, 4);
    cur_rdata = reshape(x_rdata(:,testidx)', 6, 4);
    
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
    title(sprintf('%s PCC=%.2f',...
                  meanMatrix_mets{testidx},...
                        x_data_corr(testidx)),...
                        'Interpreter', 'none');

    legend(h, curData_cols(1,:),...
             'Interpreter', 'none');%, 'Location', 'bestoutside');
    
    subplot(spx,spy,idx+2);
    
    curcoefs = x_met_smooth(2:end, testidx);
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
        %ylim([0 1])
        set(gca, 'XTick', 1:length(git_labels))
        set(gca, 'XTickLabel', git_labels)
        
        axis square
    end
%     spt = suptitle({sprintf('MZ=%.3f',testmz(idx,1)),...
%                                         testannID{1},...
%                                         testann{1}});
%     set(spt,'FontSize',8,'FontWeight','normal')
    orient landscape
    %print to figure
    print(gcf, '-painters', '-dpsc2', '-r600', '-append', '-bestfit',...
            fileNameprofiles);%[figureFolder,...
            % fileNameprofiles])
    clf('reset')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%