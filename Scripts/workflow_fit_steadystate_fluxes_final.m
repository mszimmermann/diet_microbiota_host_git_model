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
% read annotation from file
% annotationTableSpatialClusters = readtable([inputFolder ...
%     'metabolites_allions_combined_formulas_with_metabolite_filters_spatial100clusters_with_mean.csv']);
annotationTableSpatialClusters = readtable([outputFolder ...
    'metabolites_allions_combined_formulas_with_metabolite_filters_spatial100clusters_with_mean_with_CVR.csv']);
% read metabolite normalized data from file
% metaboliteData = readtable([inputFolder ...
%     'metabolites_allions_combined_norm_intensity.csv']);
metaboliteData = readtable([outputFolder ...
    'metabolites_allions_combined_norm_intensity_with_CVR.csv']);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% selected_mets are metabolites detected along the GI tract
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the mean normalized value across replicates

% define tissue order
sampleTissue_unique = {'SI1', 'SI2', 'SI3', 'Cecum', 'Colon', 'Feces'};
git_labels = {'Du', 'Je', 'Il', 'Cec', 'Col', 'Fec'};


% define which metabolites to model
%%%% all metabolites
selected_mets = 1:size(meanMatrix,1);
%%%% metabolites detected in the GIT and annotated
%selected_mets = find((sum(spatialClusters,2)>0) &...
%                     (annotationTableSpatialClusters.MetaboliteFilter==1));
% 
met_gitfits = cell(length(selected_mets),1);
met_bestsols = cell(length(selected_mets),1);

met_gitfitsCVR = cell(length(selected_mets),1);
met_bestsolsCVR = cell(length(selected_mets),1);

met_gitfitsDC = cell(length(selected_mets),1);
met_bestsolsDC = cell(length(selected_mets),1);

diag_plot_flag = 0; % diagnostic plotting flag

for flagCVR = 0:1
    sampleType_unique = unique(combinedType);
    if flagCVR
        % remove DC and leave only CVR
        sampleType_unique(ismember(sampleType_unique, {'DC'}))=[];
    
        % figure name to store profiles and model coefficients
        fileNameFigure = [figureFolder...
             'fig_CVRann_modelSMOOTH_2LIcoefHost_1LIbact.ps'];
         % file name to store solutions
        filenameAllModel = [outputFolder ...
                'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_allCVR'];
        filenameBestModel = [outputFolder ...
                'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_bestCVR'];
    else
        % remove CVR and leave only DC
        sampleType_unique(ismember(sampleType_unique, {'CVR'}))=[];
    
        % figure name to store profiles and model coefficients
        fileNameFigure = [figureFolder...
             'fig_DCann_modelSMOOTH_2LIcoefHost_1LIbact.ps'];
        % file name to store solutions
        filenameAllModel = [outputFolder ...
                'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_allDC'];
        filenameBestModel = [outputFolder ...
                'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_bestDC'];
    
    end
    
    if diag_plot_flag
        fig = figure('units','normalized','outerposition',[0 0 1 1]);
    end
    
    gut_filter = zeros(length(selected_mets),1);
    for met_i = 1:length(selected_mets)
        cmpd_interest_idx = selected_mets(met_i);
        if sum(spatialClusters(cmpd_interest_idx,:))>0
            % this metabolite is detected across the gut
            gut_filter(met_i)=1;
            % calculate mean profiles            
            idx=1;
            kmeanMatrix_joint = zeros(4,6);
            kmeanMatrix_joint_names = cell(4,6);
            for diet_i = 1:length(sampleDiet_unique)
                for type_i = 1:length(sampleType_unique)
                    selectDiet = sampleDiet_unique{diet_i};
                    selectMouse = sampleType_unique{type_i};
                    kmeanMatrix = meanMatrix(:,cellfun(@(x) contains(x,selectMouse) & contains(x,selectDiet),meanConditions));
                    kmeanMatrix = kmeanMatrix(cmpd_interest_idx,1:6);
                    kmeanMatrix_joint(idx,:) = kmeanMatrix;
                    kmeanMatrix_joint_names(idx,:) = strcat(selectDiet,...
                                                            '_',...
                                                            selectMouse,...
                                                            '_',...
                                                            git_labels);
                    idx = idx+1;
                end
            end
            % normalize by max intensity
            kmeanMatrix_joint = kmeanMatrix_joint./max(max(kmeanMatrix_joint));
            kmeanMatrix_joint(isnan(kmeanMatrix_joint))=0;
    
            if nnz(kmeanMatrix_joint)
                %[gitfit] = fitGITmodel(kmeanMatrix_joint, ncond, shuffle_flag)
                [gitfit] = fitGITmodel(kmeanMatrix_joint, kmeanMatrix_joint_names, 4, 1);
                % get best solution
                [bestsol] = select_gitfit_sol(gitfit);
                
                % save gitfit and best solution for the current metabolite
                met_gitfits{met_i} = gitfit;
                met_bestsols{met_i} = bestsol;
                
                if diag_plot_flag
                    gitfit_diag_plot(gitfit,...
                        [num2str(annotationTableSpatialClusters.MZ(cmpd_interest_idx)),'_',...
                         num2str(annotationTableSpatialClusters.RT(cmpd_interest_idx)),'_',...
                         annotationTableSpatialClusters.CompoundName{cmpd_interest_idx}],...
                         fig)
    
                    orient landscape
                    %print to figure
                    print(fig, '-vector', '-dpsc2', '-r600', '-append', '-bestfit',...
                            fileNameFigure);
                    clf('reset')
                end    
            end
        end
        % display status every 100 metabolites
        if mod(met_i,100)==0
            fprintf('Calculated parameters for %d metabolites\n',met_i);
        end
    end
    % print solutions to file
    % print best solutions to files
    % create met_info object needed for the printing function
    met_info = annotationTableSpatialClusters(selected_mets,:);
    met_info.gut_filter = gut_filter;
    % print solutions to files
    print_bestsol_to_files(met_info, met_bestsols, filenameBestModel);

    % print solutions to files
    print_allsols_to_files(met_info, met_gitfits, filenameAllModel);

    if flagCVR
        met_gitfitsCVR = met_gitfits;
        met_bestsolsCVR = met_bestsols;
    else
        met_gitfitsDC = met_gitfits;
        met_bestsolsDC = met_bestsols;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print best solutions to files
% create met_info object needed for the printing function
met_info = annotationTableSpatialClusters(selected_mets,:);
% print solutions to files
print_bestsol_to_files(met_info, met_bestsols, filenameModel);

filename = [outputFolder ...
            'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_allDC'];
% create met_info object needed for the printing function
met_info = annotationTableSpatialClusters(selected_mets,:);
% print solutions to files
print_bestsol_to_files(met_info, met_bestsols_dc, filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print all solutions to files
filename = [outputFolder ...
            'model_results_ALL_SMOOTH_raw_2LIcoefHost1LIcoefbact_annDC'];
% print solutions to files
print_allsols_to_files(met_info, met_gitfits_dc, filename);
% CVR
filename = [outputFolder ...
            'model_results_ALL_SMOOTH_raw_2LIcoefHost1LIcoefbact_annCVR'];
% print solutions to files
print_allsols_to_files(met_info, met_gitfits, filename);

% test reading from file
[met_info_read, met_gitfits_read] = read_allsols_from_files(filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fix_mz_rt_issue = 0;
if fix_mz_rt_issue==1
    % add mz and rt to the tables which were missing due to a mistake in the
    % print function
    modelFileFolder = '.\ProcessedData\output\model_results_no_mzrt\';
    modelFiles = dir(modelFileFolder);
    modelFileNames = cell(size(modelFiles,1),1);
    for i=1:length(modelFiles)
        modelFileNames{i} = modelFiles(i).name;
    end
    % leave only model files
    modelFileNames = modelFileNames(cellfun(@(x) contains(x,'model_results'),modelFileNames));
    
    for i=1:length(modelFileNames)
        currentData = readtable([modelFileFolder modelFileNames{i}]);
        %check if MZ column is indeed empty
        if sum(currentData.MZ)==0
            %compare compound IDs between currentData and annotation
            % check if each ID is the same
            if (sum(arrayfun(@(x) isequal(currentData.CompoundID{x}, ...
                    annotationTableSpatialClusters.CompoundID{x}), ...
                    1:size(currentData,1))) == size(currentData,1))
                currentData.MZ = annotationTableSpatialClusters.MZ;
                currentData.RT = annotationTableSpatialClusters.RT;
                writetable(currentData, [outputFolder, modelFileNames{i}])
                fprintf('Updated %s\n', modelFileNames{i});
            else
                fprintf('warning: MZ and RT could not be matched, file %s not updated', modelFileNames{i});
            end
        end
    end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % define which metabolites to model
% %%%% all metabolites
% selected_mets = 1:size(meanMatrix,1);
% %%%% metabolites detected in the GIT and annotated
% % selected_mets = find((sum(spatialClusters,2)>0) &...
% %                      (annotationTableSpatialClusters.MetaboliteFilter==1));
% % %%%% testing mode: hundred of annotated metabolites
% % if testing_mode_flag
% %     selected_mets = selected_mets(1:100);
% % end
% 
% % if select_by_ionmz
% %     targetMZ = [147.053; 74.037; 287.210;...]; %glutamate, propionate, l-octanoylcarnitine 
% %         499.297; 125.015; 131.058; 226.095]; %taurodeoxycholate, taurine, 5-aminolevulinate, porphobilonogen
% %     % find annotated compound with this MZ
% %     % for which modelling results correlate >0.7 with original
% %     selected_mets = find((annotationTableSpatialClusters.MetaboliteFilter>0) &...
% %                          (arrayfun(@(x) sum(abs(x-targetMZ)<=0.001),...
% %                                   annotationTableSpatialClusters.MZ)>0));
% % end
% 
% % define tissue order
% sampleTissue_unique = {'SI1', 'SI2', 'SI3', 'Cecum', 'Colon', 'Feces'};
% nrand=100;
% % call calculateAmatrix_final with zero matrices of correct size
% % to get coefvalues output, which is used to allocate variables further
% [~, coefvalues] = calculateAmatrix_final(zeros(length(condLabels), length(sampleTissue_unique)));
% 
% % allocate variables for model prediction results
% x_met_mean = zeros(length(coefvalues),length(selected_mets));
% x_met_std = zeros(length(coefvalues),length(selected_mets));
% x_met_smooth = zeros(length(coefvalues),length(selected_mets));
% 
% % allocate variables for correlation calculated from values restored with
% % reverse problem with original and randomly shuffled parameters
% x_data_corr = zeros(length(selected_mets),1);
% x_resid = zeros(length(selected_mets),1);
% x_data_corr_shuffled = zeros(length(selected_mets),1);
% 
% % calculate reverse A matrix with zero matrix of the right size to allocate
% % variables according to Ra size
% [Ra] = calculateRAmatrix_final(zeros(size(coefvalues)));
% x_rdata = zeros(size(Ra,2),length(selected_mets));
% x_rdata_shuffled = zeros(size(Ra,2),length(selected_mets));
% 
% minval = 0.01;
% 
% fprintf('Starting parameter estimation for %d metabolites\n',length(selected_mets));
% 
% % remove DC and leave only CVR
% sampleType_unique(ismember(sampleType_unique, {'DC'}))=[];
% 
% for met_i = 1:length(selected_mets)
% 
%     if sum(spatialClusters(met_i,:))>0
%         cmpd_interest_idx = selected_mets(met_i);
%         % calculate mean profiles            
%         idx=1;
%         kmeanMatrix_joint = zeros(4,6);
%         for diet_i = 1:length(sampleDiet_unique)
%             for type_i = 1:length(sampleType_unique)
%                 selectDiet = sampleDiet_unique{diet_i};
%                 selectMouse = sampleType_unique{type_i};
%                 kmeanMatrix = meanMatrix(:,cellfun(@(x) contains(x,selectMouse) & contains(x,selectDiet),meanConditions));
%                 kmeanMatrix = kmeanMatrix(cmpd_interest_idx,1:6);
%                 kmeanMatrix_joint(idx,:) = kmeanMatrix;
%                 idx = idx+1;
%             end
%         end
%         % normalize by max intensity
%         kmeanMatrix_joint = kmeanMatrix_joint./max(max(kmeanMatrix_joint));
%         % replace small values with noise
%         kmeanMatrix_joint_orig = kmeanMatrix_joint;
%         kmeanMatrix_joint(kmeanMatrix_joint<minval)=...
%         kmeanMatrix_joint(kmeanMatrix_joint<minval)+rand(nnz(kmeanMatrix_joint<minval),1)*minval;
% 
%         % calculate fluxes for average profile
%         Aorig = calculateAmatrix_final(kmeanMatrix_joint_orig);
%         [A,coefvalues] = calculateAmatrix_final(kmeanMatrix_joint);
%         A(Aorig(:,1)==0,:)=[];
%         b = zeros(size(A,1),1);
%         options = optimoptions(@lsqlin,'Display', 'off');
% 
%         xlowerlim = zeros(size(A,2),1);
%         xlowerlim(2:end)=-Inf;
% 
%         %%%%%%%%%%%%%%%%%%%%%%%%%
%         %set bact to zero
%         xupperlim = Inf*ones(size(A,2),1);
%         %%%%%%%%%%%%%%%%%%%%%%%%%
%         % check if A is empty an if yes move to the next metabolite
%         if isempty(A)
%             continue;
%         end
% 
%         [x, xres] = lsqlin(A,b,[],[],[],[],xlowerlim, xupperlim, [], options);
%         x_met_smooth(:, met_i) = x;
%         x_resid(met_i) = xres;
% 
% 
%     [Ra,rb] = calculateRAmatrix_final(x*1000);
%     dataR = lsqlin(Ra,rb,[],[],[],[],zeros(1,size(Ra,2)), [], [], options);
%     % save reconstructed data
%     x_rdata(:,met_i) = dataR;
% 
%     %rsshape to matrix form
%     dataR = reshape(dataR,[],4)';
% 
%     x_data_corr(met_i) = corr(kmeanMatrix_joint_orig(:), dataR(:));
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % calculate reverse problem for shuffled x
%     x_shuffled = x(randperm(length(x)));
%     [Ra_shuffled,rb_shuffled] = calculateRAmatrix_final(x_shuffled*1000);
%     dataR_shuffled = lsqlin(Ra_shuffled,rb_shuffled,[],[],[],[],...
%                             zeros(1,size(Ra_shuffled,2)), [], [], options);
%     % save reconstructed data
%     x_rdata_shuffled(:,met_i) = dataR_shuffled;
% 
%     %rsshape to matrix form
%     dataR_shuffled = reshape(dataR_shuffled,[],4)';
% 
%     x_data_corr_shuffled(met_i) = corr(kmeanMatrix_joint_orig(:), dataR_shuffled(:));
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % calculate joint matrix per mouse
%         kmeanMatrix_joint_mouse = cell(5,1);
%         for mouse_i = 1:5
%             kmeanMatrix_joint = zeros(4,6);
%             kmeanMatrix_joint_samples = cell(4,6);
%             idx=1;
%             for diet_i = 1:length(sampleDiet_unique)
%                 for type_i = 1:length(sampleType_unique)
%                     for tissue_i = 1:length(sampleTissue_unique)
%                         selectDiet = sampleDiet_unique{diet_i};
%                         selectType = sampleType_unique{type_i};
%                         selectTissue = sampleTissue_unique{tissue_i};
%                         kmeanMatrix = combinedIntensitiesNorm(cmpd_interest_idx,...
%                                                 cellfun(@(x) contains(x,selectType),combinedType) & ...
%                                                 cellfun(@(x) contains(x,selectDiet),combinedDiet) & ...
%                                                 cellfun(@(x) contains(x,selectTissue),combinedTissues));
%                         if length(kmeanMatrix)<mouse_i
%                             kmeanMatrix_joint(idx,tissue_i) = mean(kmeanMatrix);
%                         else
%                             kmeanMatrix_joint(idx,tissue_i) = kmeanMatrix(mouse_i);
%                         end
%                         kmeanMatrix_joint_samples{idx,tissue_i} = ['P_' lower(strrep(selectTissue,'DC','feces')) ...
%                                                                    '_' lower(strrep(selectType,'DC','cv'))...
%                                                                    '_' lower(selectDiet)];
%                     end
%                     idx = idx+1;
%                 end
%             end
% 
%             kmeanMatrix_joint = kmeanMatrix_joint./max(max(kmeanMatrix_joint));
%             kmeanMatrix_joint_mouse{mouse_i} = kmeanMatrix_joint;
%         end
% 
%         x_rand = zeros(length(coefvalues), nrand);
%         for i=1:nrand
%             curmice = randi(5,1,4);
%             kmeanMatrix_joint = zeros(4,6);
%             for j=1:length(curmice)
%                 kmeanMatrix_joint(j,:) = kmeanMatrix_joint_mouse{curmice(j)}(j,:);
%             end
%             % replace small values with noise
%             kmeanMatrix_joint_orig = kmeanMatrix_joint;
%             kmeanMatrix_joint(kmeanMatrix_joint<minval)=...
%                 kmeanMatrix_joint(kmeanMatrix_joint<minval)+rand(nnz(kmeanMatrix_joint<minval),1)*minval;
% 
% 
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Aorig = calculateAmatrix_final(kmeanMatrix_joint_orig);
%             [A,coefvalues] = calculateAmatrix_final(kmeanMatrix_joint);
%             A(Aorig(:,1)==0,:)=[];
%             b = zeros(size(A,1),1);
%             options = optimoptions(@lsqlin,'Display', 'off');
% 
%             xlowerlim = zeros(size(A,2),1);
%             xlowerlim(2:end)=-Inf;
% 
%             %%%%%%%%%%%%%%%%%%%%%%%%%
%             %set bact to zero
%             xupperlim = Inf*ones(size(A,2),1);
%             %%%%%%%%%%%%%%%%%%%%%%%%%
% 
%             if size(A,1)>0
%                 x = lsqlin(A,b,[],[],[],[],xlowerlim, xupperlim, [], options);
%                 x_rand(:, i) = x;
%             end
% 
%         end
%         %  remove all zero solutions
%         x_rand(:,sum(x_rand==0,1)==size(x,1))=[];
%         if size(x_rand,1)>0
%             x_met_mean(:,met_i) = mean(x_rand,2);
%             x_met_std(:,met_i) = std(x_rand,[],2);
%         end
%     end
%     % display status every 100 metabolites
%     if mod(met_i,100)==0
%         fprintf('Calculated parameters for %d metabolites\n',met_i);
%     end
% end
% 
% % get the original metabolomics data as vectors to calculate correlations
% % between original and restored measurements
% kmean_vector_joint_orig = zeros(size(Ra,2),length(selected_mets));
% kmean_vector_joint = zeros(size(Ra,2),length(selected_mets));
% for met_i = 1:length(selected_mets)
% 
%     cmpd_interest_idx = selected_mets(met_i);
%     % calculate mean profiles            
%     idx=1;
%     kmeanMatrix_joint = zeros(4,6);
%     for diet_i = 1:length(sampleDiet_unique)
%         for type_i = 1:length(sampleType_unique)
%             selectDiet = sampleDiet_unique{diet_i};
%             selectMouse = sampleType_unique{type_i};
%             kmeanMatrix = meanMatrix(:,cellfun(@(x) contains(x,selectMouse) & contains(x,selectDiet),meanConditions));
%             % get serum and liver data
%             kmeanMatrix = kmeanMatrix(cmpd_interest_idx,1:6);
%             kmeanMatrix_joint(idx,:) = kmeanMatrix;
%             idx = idx+1;
%         end
%     end
%     % normalize by max intensity
%     kmeanMatrix_joint = kmeanMatrix_joint./max(max(kmeanMatrix_joint));
%     % replace small values with noise
%     kmeanMatrix_joint_orig = kmeanMatrix_joint;
%     kmeanMatrix_joint(kmeanMatrix_joint<minval)=...
%     kmeanMatrix_joint(kmeanMatrix_joint<minval)+rand(nnz(kmeanMatrix_joint<minval),1)*minval;
% 
%     kmean_vector_joint_orig(:, met_i) = kmeanMatrix_joint_orig(:);
%     kmean_vector_joint(:, met_i) = kmeanMatrix_joint(:);
% end
% 
% % calculate Spearman correlations
% x_corr_Spearman = zeros(size(x_data_corr));
% corr_Mouse = zeros(length(x_data_corr),4);
% corr_Mouse_Spearman = zeros(length(x_data_corr),4);
% 
% for i = 1:length(x_data_corr)
%     dataR = reshape(x_rdata(:,i),[],4)';
%     dataOrig = reshape(kmean_vector_joint_orig(:,i),4,[]);
% 
%     x_corr_Spearman(i) = corr(dataOrig(:), dataR(:), 'type', 'Spearman');
%     corr_Mouse(i,:) = diag(corr(dataOrig', dataR'));
%     corr_Mouse_Spearman(i,:) = diag(corr(dataOrig', dataR', 'type', 'Spearman'));
% 
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % save model results to file - raw coefficients
% % add gut filter to file - flag indicating whether metabolites were
% % detected in the GIT
% gut_filter = (sum(spatialClusters,2)>0);
% 
% fid = fopen([outputFolder ...
%             'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_allions_with_CVR.csv'], 'w');
%              %'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_allions.csv'], 'w');
% fprintf(fid, 'MZ\tRT\tCompoundID\tCompoundName\tMetaboliteFilter\tSumGITclusters\tReciprocalCorr');
% for i=1:length(coefvalues)
%     fprintf(fid, '\t%s', coefvalues{i});
% end
% fprintf(fid, '\n');
% for i=1:size(x_met_smooth,2)
%     fprintf(fid, '%.3f\t%3f', annotationTableSpatialClusters.MZ(i),...
%         annotationTableSpatialClusters.RT(i));
%     fprintf(fid, '\t%s\t%s\t%d', annotationTableSpatialClusters.CompoundID{i},...
%         annotationTableSpatialClusters.CompoundName{i},...
%         annotationTableSpatialClusters.MetaboliteFilter(i));
%     fprintf(fid, '\t%d', gut_filter(i));
%     fprintf(fid, '\t%.3f', x_data_corr(i));
%     for j=1:size(x_met_smooth,1)
%         fprintf(fid, '\t%e', x_met_smooth(j,i));
%     end
%     fprintf(fid, '\n');
% end
% fclose(fid);
% 
% % save model results to file - normalized by max coefficient
% fid = fopen([outputFolder...
%     'model_results_SMOOTH_normbyabsmax_2LIcoefHost1LIcoefbact_allions_with_CVR.csv'], 'w');
%     %'model_results_SMOOTH_normbyabsmax_2LIcoefHost1LIcoefbact_allions.csv'], 'w');
% fprintf(fid, 'MZ\tRT\tCompoundID\tCompoundName\tMetaboliteFilter\tSumGITclusters\tReciprocalCorr');
% for i=1:length(coefvalues)
%     fprintf(fid, '\t%s', coefvalues{i});
% end
% fprintf(fid, '\n');
% for i=1:size(x_met_smooth,2)
%     fprintf(fid, '%.3f\t%3f', annotationTableSpatialClusters.MZ(i),...
%         annotationTableSpatialClusters.RT(i));
%     fprintf(fid, '\t%s\t%s\t%d', annotationTableSpatialClusters.CompoundID{i},...
%         annotationTableSpatialClusters.CompoundName{i},...
%         annotationTableSpatialClusters.MetaboliteFilter(i));
%     fprintf(fid, '\t%d', gut_filter(i));
%     fprintf(fid, '\t%.3f', x_data_corr(i));
%     for j=1:size(x_met_smooth,1)
%         fprintf(fid, '\t%.3f', x_met_smooth(j,i)./max(abs(x_met_smooth(:,i))));
%     end
%     fprintf(fid, '\n');
% end
% fclose(fid);
% 
% % save model results to file - metabolism_coefficients_only
% fid = fopen([outputFolder...
%     'model_results_SMOOTH_normbyabsmax_ONLYMETCOEF_2LIcoefHost1LIcoefbact_allions_with_CVR.csv'], 'w');
%     %'model_results_SMOOTH_normbyabsmax_ONLYMETCOEF_2LIcoefHost1LIcoefbact_allions.csv'], 'w');
% fprintf(fid, 'MZ\tRT\tCompoundID\tCompoundName\tMetaboliteFilter\tSumGITclusters\tReciprocalCorr');
% for i=2:length(coefvalues)
%     fprintf(fid, '\t%s', coefvalues{i});
% end
% fprintf(fid, '\n');
% for i=1:size(x_met_smooth,2)
%     fprintf(fid, '%.3f\t%3f', annotationTableSpatialClusters.MZ(i),...
%         annotationTableSpatialClusters.RT(i));
%     fprintf(fid, '\t%s\t%s\t%d', annotationTableSpatialClusters.CompoundID{i},...
%         annotationTableSpatialClusters.CompoundName{i},...
%         annotationTableSpatialClusters.MetaboliteFilter(i));
%     fprintf(fid, '\t%d', gut_filter(i));
%     fprintf(fid, '\t%.3f', x_data_corr(i));
%     for j=2:size(x_met_smooth,1)
%         fprintf(fid, '\t%.3f', x_met_smooth(j,i)./max(abs(x_met_smooth(2:end,i))));
%     end
%     fprintf(fid, '\n');
% end
% fclose(fid);
% 
% % save model results to file - reciprocal data restoration
% fid = fopen([outputFolder...
%     'model_results_SMOOTH_normbyabsmax_reciprocal_problem_allions_with_CVR.csv'], 'w');
%     %'model_results_SMOOTH_normbyabsmax_reciprocal_problem_allions.csv'], 'w');
% columnNames = cell(size(x_rdata,1),1);
% idx=1;
% for diet_i = 1:length(sampleDiet_unique)
%     for type_i = 1:length(sampleType_unique)
%         for tiss_i = 1:length(sampleTissue_order)-2
%                 selectDiet = sampleDiet_unique{diet_i};
%                 selectMouse = sampleType_unique{type_i};
%                 selectTissue = sampleTissue_order{tiss_i};
%                 columnNames{idx} = [selectDiet '_' selectMouse '_' selectTissue];
%                 idx = idx+1;
%         end
%     end
% end
% 
% fprintf(fid, 'MZ\tRT\tCompoundID\tCompoundName\tMetaboliteFilter\tSumGITclusters\tReciprocalCorr\tRandomCorr');
% for i=1:length(columnNames)
%     fprintf(fid, '\t%s', columnNames{i});
% end
% for i=1:length(columnNames)
%     fprintf(fid, '\tRecip_%s', columnNames{i});
% end
% for i=1:length(columnNames)
%     fprintf(fid, '\tRandom_%s', columnNames{i});
% end
% fprintf(fid, '\n');
% for i=1:size(x_rdata,2)
%     fprintf(fid, '%.3f\t%3f', annotationTableSpatialClusters.MZ(i),...
%         annotationTableSpatialClusters.RT(i));
%     fprintf(fid, '\t%s\t%s\t%d', annotationTableSpatialClusters.CompoundID{i},...
%         annotationTableSpatialClusters.CompoundName{i},...
%         annotationTableSpatialClusters.MetaboliteFilter(i));
%     fprintf(fid, '\t%d', gut_filter(i));
%     fprintf(fid, '\t%.3f', x_data_corr(i));
%     fprintf(fid, '\t%.3f', x_data_corr_shuffled(i));
%     for j=1:size(kmean_vector_joint_orig(:,i),1)
%         fprintf(fid, '\t%.3f', kmean_vector_joint_orig(j,i));
%     end
%     for j=1:size(x_rdata(:,i),1)
%         fprintf(fid, '\t%.3f', x_rdata(j,i));
%     end
%     for j=1:size(x_rdata_shuffled(:,i),1)
%         fprintf(fid, '\t%.3f', x_rdata_shuffled(j,i));
%     end
%     fprintf(fid, '\n');
% end
% fclose(fid);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % plot correlation of restored and original data
% % %testing mode: all calculated
% if testing_mode_flag
%     met_filter = ones(size(x_data_corr,1),1);
% else
%     % %metabolites detected in the GIT and annotated
%     met_filter = ((annotationTableSpatialClusters.MetaboliteFilter==1) &...
%                (sum(spatialClusters,2)>0));
% end
% 
% % calculate differentce in corr distrbutions
% p_corr_diff = ranksum(x_data_corr_shuffled(met_filter==1),...
%     x_data_corr(met_filter==1));
% figure
% % pearson corr all 
% h = histogram(x_data_corr_shuffled(met_filter==1),100);
% hold on
% histogram(x_data_corr(met_filter==1),100)
% xlim([-1 1])
% axis square
% xlabel('Pearson correlation between metabolomics data and model estimate')
% ylabel('Number of ions')
% title('All data together')
% orient landscape
% % print line for PCC=0.7
% plot([0.7, 0.7], [0, max(h.Values)], 'k--')
% 
% legend({'Random coefficients', 'Model coefficients', 'PCC=0.7'},...
%         'Location', 'NorthWest')
% 
% % compare distributions of correlations
% pval = signrank(x_data_corr_shuffled(met_filter==1),...
%                 x_data_corr(met_filter==1));
% % print Wilcoxon signed rank test o-value on the plot
% text(-0.9, 0.7*max(h.Values), sprintf('signrank p = %.2e', pval))
% % save figure to file
% print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
%     [figureFolder,...
%     'Fig4a_histogram_corr_model_2LIcoefHost1LIcoefbact_reversedata_annotated_CVR'])
%        %'Fig4a_histogram_corr_model_2LIcoefHost1LIcoefbact_reversedata_annotated'])


