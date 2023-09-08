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

% mycolors = [0 115 178;... %dark blue
%             204 227 240;...%light blue
%             211 96 39;... %dark orange
%             246 223 212]/256;%light orange
% [bb,aa] = ndgrid(sampleType_unique, sampleDiet_unique); 
% condLabels = strcat(aa(:),'-', bb(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% selected_mets are metabolites detected along the GI tract
selected_mets = 1:length(meanMatrix_mets);
% define tissue order
sampleTissue_unique = {'SI', 'SII', 'SIII', 'Cecum', 'Colon', 'Feces'};
%nrand=100;
% call calculateAmatrix_final with zero matrices of correct size
% to get coefvalues output, which is used to allocate variables further
%[~, coefvalues] = calculateAmatrix_final(zeros(length(condLabels), length(sampleTissue_unique)));

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
% % save Rsquared of total
% % calculate corr and Rsquared of SI and LI separately
% x_data_Rsq = zeros(length(selected_mets),1);
% x_data_Rsq_shuffled = zeros(length(selected_mets),1);
% x_data_Rsq_SI = zeros(length(selected_mets),1);
% x_data_Rsq_SI_shuffled = zeros(length(selected_mets),1);
% x_data_Rsq_LI = zeros(length(selected_mets),1);
% x_data_Rsq_LI_shuffled = zeros(length(selected_mets),1);
% x_data_corr_SI = zeros(length(selected_mets),1);
% x_data_corr_SI_shuffled = zeros(length(selected_mets),1);
% x_data_corr_LI = zeros(length(selected_mets),1);
% x_data_corr_LI_shuffled = zeros(length(selected_mets),1);

% calculate reverse A matrix with zero matrix of the right size to allocate
% variables according to Ra size
% [Ra] = calculateRAmatrix_final(zeros(size(coefvalues)));
% x_rdata = zeros(size(Ra,2),length(selected_mets));
% x_rdata_shuffled = zeros(size(Ra,2),length(selected_mets));

%minval = 0.01;

fprintf('Starting parameter estimation for %d metabolites\n',length(selected_mets));

volumeMatrix_CVR = repmat([0.3 0.3 0.3 0.3 0.3 0.3;...
                           0.3 0.3 0.3 3 3 3], 2, 1);
volumeMatrix_GF = repmat([0.3 0.3 0.3 3 3 3], 4, 1);


% use volume flag
use_volume_flag = 1;
if use_volume_flag
%     fileNameFigure = [figureFolder...
%         'fig_2023_drugdiag_volume_profiles_modelSMOOTH_2LIcoefHost_1LIbact.ps'];
    fileNameFigure = [figureFolder...
        'fig_2023_drugdiag_volume_4profiles_plus_mut_modelSMOOTH_2LIcoefHost_1LIbact.ps'];
else
%     fileNameFigure = [figureFolder...
%         'fig_2023_drugdiag_conc_profiles_modelSMOOTH_2LIcoefHost_1LIbact.ps'];

    fileNameFigure = [figureFolder...
        'fig_2023_drugdiag_conc_4profiles_plus_mut_modelSMOOTH_2LIcoefHost_1LIbact.ps'];
end

diag_plot_flag = 0; % diagnostic plotting flag
if diag_plot_flag
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
end

met_gitfits = cell(length(selected_mets),1);
met_bestsols = cell(length(selected_mets),1);

for met_i = 1:length(selected_mets)
   
        cmpd_interest_idx = selected_mets(met_i);
        % set volume to CV or GF/WT
        if contains(meanMatrix_mets{cmpd_interest_idx}, '_CV')
            volumeMatrix = volumeMatrix_CVR;
        else
            volumeMatrix = volumeMatrix_GF;
        end
        % calculate mean profiles            
        idx=1;
        kmeanMatrix_joint = zeros(4,6);
        kmeanMatrix_joint_names = cell(4,6);
        for diet_i = 1:length(sampleDiet_unique)
            for type_i = 1:length(sampleType_unique)
                selectDiet = sampleDiet_unique{diet_i};
                selectMouse = sampleType_unique{type_i};
                kmeanMatrix = meanMatrix(:,cellfun(@(x) contains(x,selectMouse) & ...
                                        contains(x,selectDiet),meanConditions));
                kmeanMatrix = kmeanMatrix(cmpd_interest_idx,1:6);
                kmeanMatrix_joint(idx,:) = kmeanMatrix;
                kmeanMatrix_joint_names(idx,:) = strcat(selectDiet,...
                                                        '_',...
                                                        selectMouse,...
                                                        '_',...
                                                        sampleTissue_unique);
                idx = idx+1;
            end
        end
        
        if use_volume_flag
            % multiply by volume to get amounts
            kmeanMatrix_joint = kmeanMatrix_joint.*volumeMatrix; 
        end
        % normalize by max intensity
        kmeanMatrix_joint = kmeanMatrix_joint./max(max(kmeanMatrix_joint));
        kmeanMatrix_joint(isnan(kmeanMatrix_joint))=0;

       
        %[gitfit] = fitGITmodel(kmeanMatrix_joint, ncond, shuffle_flag)
        [gitfit] = fitGITmodel(kmeanMatrix_joint, kmeanMatrix_joint_names, 2, 1);
        % get best solution
        [bestsol] = select_gitfit_sol(gitfit);
        
        % save gitfit and best solution for the current metabolite
        met_gitfits{met_i} = gitfit;
        met_bestsols{met_i} = bestsol;
        
        if diag_plot_flag
            gitfit_diag_plot(gitfit,meanMatrix_mets{met_i},fig)

            orient landscape
            %print to figure
            print(fig, '-painters', '-dpsc2', '-r600', '-append', '-bestfit',...
                    fileNameFigure);
            clf('reset')
        end    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot correlation of restored and original data
% calculate differentce in corr distrbutions
filename = [figureFolder,...
    'Fig4a_histogram_MAXcorr_model_2LIcoefHost1LIcoefbact_drugs'];
plot_gitfit_model_corr(met_gitfits, filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select LI within best total solution
x_selected = zeros(size(met_bestsols{1}.x,1), length(met_bestsols));
x_data_corr = zeros(1, length(met_bestsols));
x_data_corr_SI = zeros(1, length(met_bestsols));
x_data_corr_LI = zeros(1, length(met_bestsols));
x_data_corr_mean = zeros(1, length(met_bestsols));

for i=1:length(met_bestsols)
    picksol = ismember(met_bestsols{i}.selection_criterion,...
                                {'LI PCC within high total'});
    if size(met_bestsols{i}.x,2)==length(picksol)
        x_selected(:,i) = met_bestsols{i}.x(:,picksol);
        x_data_corr(i) = met_bestsols{i}.x_sel_CorrRev(picksol);
        x_data_corr_SI(i) = met_bestsols{i}.x_sel_CorrRevSI(picksol);
        x_data_corr_LI(i) = met_bestsols{i}.x_sel_CorrRevLI(picksol);
        x_data_corr_mean(i) = mean([x_data_corr(i) x_data_corr_SI(i) x_data_corr_LI(i)]);
    else
        % there was no total PCC above thresold to select LI PCC,
        % get the best possible total PCC
        picksol = ismember(met_bestsols{i}.selection_criterion,...
                                {'total PCC'});
        picksol = picksol(1:size(met_bestsols{i}.x,2));
        x_selected(:,i) = met_bestsols{i}.x(:,picksol);
        x_data_corr(i) = met_bestsols{i}.x_sel_CorrRev(picksol);
        x_data_corr_SI(i) = met_bestsols{i}.x_sel_CorrRevSI(picksol);
        x_data_corr_LI(i) = met_bestsols{i}.x_sel_CorrRevLI(picksol);
        x_data_corr_mean(i) = mean([x_data_corr(i) x_data_corr_SI(i) x_data_corr_LI(i)]);
    end        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print best solutions to files
filename = [outputFolder ...
            'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_drugs'];
% create met_info object needed for the printing function
met_info.CompoundID = meanMatrix_mets;
met_info.CompoundName = meanMatrix_mets;          
% print solutions to files
print_bestsol_to_files(met_info, met_bestsols, filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% heatmap coefficients and metabolites
plotorder = [1:3:24 2:3:24 3:3:24 25:length(meanMatrix_mets)];
plotdata = x_data_corr(plotorder);
% set values below threshold to 0
% plotdata = x_data_corr(plotorder).*...
%             ( (x_data_corr(plotorder)>=0.7) &...
%               ((x_data_corr_SI(plotorder)>=0.7) |...
%                (x_data_corr_LI(plotorder)>=0.7)));
plotdata_names = meanMatrix_mets(plotorder);
% reshape in matrix form
plotdata = reshape(plotdata,4,[])';
plotdata_names = reshape(plotdata_names,4,[])';

% obsolete reshaping when clonazepam samples contain one more timepoint
% plotdata = [[reshape(plotdata(1:32),4,[])' zeros(8,1)];...
%             reshape(plotdata(33:end),5,[])'];
% plotdata_names = [[reshape(plotdata_names(1:32),4,[])' repmat({'NaN'},8,1)];...
%                    reshape(plotdata_names(33:end,:),5,[])'];
plotdata_names_compounds = cellfun(@(x) strrep(x,'_','-'),plotdata_names(:,1),'unif',0);
plotdata_times = {'3', '5', '7', '9'};%, '12'};

figure
heatmap(plotdata,...
        'YDisplayLabels', plotdata_names(:,1),...
        'XDisplayLabels', plotdata_times,...
        'CellLabelFormat', '%0.2g')     
caxis([0 1])
colormap(flipud(parula))
title('Max total PCC')  
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
     [figureFolder,...
     'FigSX_heatmap_maxTotalPCC_volume_drug_data'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot four creteria for each metabolite: 
% total corr, LI corr, SI corr and mean corr
plotdataSI = x_data_corr_SI(plotorder);
plotdataSI = reshape(plotdataSI,4,[])';
% plotdataSI = [[reshape(plotdataSI(1:32),4,[])' zeros(8,1)];...
%                reshape(plotdataSI(33:end),5,[])'];
plotdataLI = x_data_corr_LI(plotorder);
plotdataLI = reshape(plotdataLI,4,[])';
% plotdataLI = [[reshape(plotdataLI(1:32),4,[])' zeros(8,1)];...
%                reshape(plotdataLI(33:end),5,[])'];
plotdataMEAN = x_data_corr_mean(plotorder);
plotdataMEAN = reshape(plotdataMEAN,4,[])';
% plotdataMEAN = [[reshape(plotdataMEAN(1:32),4,[])' zeros(8,1)];...
%                reshape(plotdataMEAN(33:end),5,[])'];
       
AB = reshape([plotdata;plotdataMEAN], size(plotdata,1), []);
CD = reshape([plotdataLI; plotdataSI], size(plotdataLI,1), []);
[m,n] = size(AB);
ABCD = zeros(2*m,n);
ABCD(1:2:end,:) = AB;
ABCD(2:2:end,:) = CD;
figure
heatmap(ABCD,...
        'YDisplayLabels', reshape(repmat(plotdata_names_compounds',2,1),[],1),...
        'XDisplayLabels', reshape(repmat(plotdata_times,2,1),[],1),...
        'CellLabelFormat', '%0.2g');
caxis([-1 1])
caxis([0 1])
colormap(flipud(parula))
title({'total PCC | mean PCC', 'LI PCC | SI PCC'})
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
     [figureFolder,...
     'FigSX_heatmap_4criteria_volume_drug_data'])
  
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
fileNameprofiles = 'fig_2023_moreGOF_volume_profiles_figureselected_modelSMOOTH_2LIcoefHost_1LIbact_drugs.ps';
          
curData_cols = reshape(meanConditions, 6, 4);
compoundsInterest = 1:length(meanMatrix_mets);

fig = figure('units','normalized','outerposition',[0 0 1 1]);

fprintf('Plotting profiles for %d metabolites\n',length(compoundsInterest));
for cpdix=1:length(compoundsInterest)
    
    testidx = compoundsInterest(cpdix);

    % set volume to CV or GF/WT
%         if contains(meanMatrix_mets{testidx}, '_CV')
%             volumeMatrix = volumeMatrix_CVR;
%         else
%             volumeMatrix = volumeMatrix_GF;
%         end
        
    spx=1;
    spy=3;
    spidx = 1;
    coloridx = 1;
    idx=1;
    curmat = zeros(4,6);
    cur_data = reshape(kmean_vector_joint_orig(:,testidx)', 4, 6);
    %cur_data=cur_data.*volumeMatrix;
    cur_data = cur_data';
    %cur_data = reshape(meanMatrix(testidx,:), 6, 4);
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
        ylabel('Original normbymax x volume')

        subplot(spx,spy,idx+1);
        hold on
        plot((cur_rdata(:,i)),...
                 'LineWidth', 2,...
                 'Color', mycolors(i,:));
        ylabel('Restored normbymax')

    end
    title({sprintf('%s', meanMatrix_mets{testidx}),...
           sprintf('PCC=%.2f PCC_SI=%.2f PCC_LI=%.2f',...
                   x_data_corr(testidx),x_data_corr_SI(testidx),x_data_corr_LI(testidx)),...
           sprintf('RSQ=%.2f RSQ_SI=%.2f RSQ_LI=%.2f',...
                   x_data_Rsq(testidx),x_data_Rsq_SI(testidx),x_data_Rsq_LI(testidx))},...
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