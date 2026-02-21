% Analyze and model public data from Meier et al 2023 to compare model
% coefficients for overlapping metabolites

% call script defining file dependenciesand global variables

addpath(genpath('.\'));
add_global_and_file_dependencies

% load data from Meier et al and reformat for the model

fileName = '.\InputData\public_data\2023_Meier_NatMet_42255_2023_802_MOESM3_ESM.xlsx';
 
 meier_data = readtable(fileName,...
                        'Sheet', 'intensities');
 %meier_data = readtable(fileName, 'Sheet', 'concentration');
 % select only colonic measurements (c-colonic, m-mucosal)
 select_content = ismember(meier_data.habitat, 'c'); 
 meier_matrix = meier_data{select_content, 5:end}';
 meier_mets = meier_data.Properties.VariableNames(5:end)';
 meier_condition = meier_data.Var1(select_content);
 meier_location = meier_data.site(select_content);
 meier_replicate = meier_data.replicate(select_content);
 
 meier_condition_unique = unique(meier_condition);
 meier_location_unique = unique(meier_location);
 select_locations = [3 6 9 12 14 15];
 
 % make a boxplot of all data
 fig = figure('units','normalized','outerposition',[0 0 1 1]);
 boxplot(log10(meier_matrix))
 set(gca, 'XTick', 1:length(meier_location))
 set(gca, 'XTickLabel', meier_location)
 
 % perform quantile normalization
 meier_matrix_norm =  meier_matrix;
 % replace 0 with nan
 meier_isnan = (meier_matrix_norm==0);
 meier_matrix_norm = quantilenorm(meier_matrix_norm);
 meier_matrix_norm(meier_isnan) = nan;
  
 fig = figure('units','normalized','outerposition',[0 0 1 1]);
 boxplot(log10(meier_matrix_norm))
 set(gca, 'XTick', 1:length(meier_location))
 set(gca, 'XTickLabel', meier_location)

% define colors for different mouse groups and tissue names
mycolors = [0 115 178; %dark blue for both 204 227 240;...%light blue
            0 115 178]/256;%dark blue
mylinestyle = {'--', '-'};
git_labels = {'Du', 'Je', 'Il', 'Cec', 'Col', 'Fec'};

meier_figure_file_name = [figureFolder...
    'figSX_Meier_data_plots_smooth_by_sum.ps'];

% smooth first metabolite to create placeholder matrices
[smoothloc, smoothcond, smoothrep, smoothdata] = ...
        smooth_logitudinal_data(meier_location,...
        meier_condition,...
        meier_replicate,...
        meier_matrix_norm(1, :));
meier_data_smooth = zeros(length(meier_mets), length(smoothloc));   
 % plot all datapoints separately
printflag=0;
if printflag
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
end
for met_i=1:length(meier_mets)
     % smooth data for the current metabolite
     [smoothloc, smoothcond, smoothrep, smoothdata] = ...
        smooth_logitudinal_data(meier_location,...
        meier_condition,...
        meier_replicate,...
        meier_matrix_norm(met_i, :));
     smoothdata_normbymax = smoothdata./max(max(smoothdata));
     
     if printflag
         clf
         hold on
         for cond_i = 1:length(meier_condition_unique)
             subplot(2,2,1)
             % original data
             hold on
             curdata = meier_matrix(met_i, ...
                 ismember(meier_condition, meier_condition_unique{cond_i}));
             curloc = meier_location(ismember(meier_condition, meier_condition_unique{cond_i}));
             lg1(cond_i) = scatter(curloc,...
                 curdata, 'filled');
             title('Raw intensities')

             subplot(2,2,2)
             % quantile normalized
             hold on
             curdata = meier_matrix_norm(met_i, ...
                 ismember(meier_condition, meier_condition_unique{cond_i}));
             curloc = meier_location(ismember(meier_condition, meier_condition_unique{cond_i}));
             lg2(cond_i) = scatter(curloc,...
                 curdata, 'filled');
             title('Quantile-normalized')

             subplot(2,2,3)
             % selected sections, normalized
             hold on
             curdata = smoothdata(ismember(smoothcond, meier_condition_unique{cond_i}));
             curloc = smoothloc(ismember(smoothcond, meier_condition_unique{cond_i}));
             plot(curloc, curdata,...
                 'o', 'MarkerSize', 5, 'MarkerEdgeColor', mycolors(cond_i,:),...
                                     'MarkerFaceColor', mycolors(cond_i,:));

             kmeanMatrix = zeros(length(unique(curloc)),1);
             kstdMatrix = zeros(length(unique(curloc)),1);
             for tissue_i=1:length(unique(curloc))
                kmeanMatrix(tissue_i) = nanmean(curdata(curloc==tissue_i));
                kstdMatrix(tissue_i) = nanstd(curdata(curloc==tissue_i));
             end
             lg3(cond_i) = plot(1:length(unique(curloc)), kmeanMatrix,...
                   mylinestyle{cond_i},...
                  'LineWidth', 2, 'Color', mycolors(cond_i,:));
              errorbar(1:length(unique(curloc)), kmeanMatrix, kstdMatrix,...
                  '.','LineWidth', 2, 'Color', mycolors(cond_i,:)) 
              axis square
              set(gca, 'XTick', 1:length(git_labels))
              set(gca, 'XTickLabel', git_labels)
              xlim([1,length(git_labels)])  
              title('Smoothed')
             % normalize by max
             subplot(2,2,4)
             hold on
             % selected sections, normalized
             curdata = smoothdata_normbymax(ismember(smoothcond,...
                 meier_condition_unique{cond_i}));
             plot(curloc,curdata,...
                 'o', 'MarkerSize', 5, 'MarkerEdgeColor', mycolors(cond_i,:),...
                                     'MarkerFaceColor', mycolors(cond_i,:));

             kmeanMatrix = zeros(length(unique(curloc)),1);
             kstdMatrix = zeros(length(unique(curloc)),1);
             for tissue_i=1:length(unique(curloc))
                kmeanMatrix(tissue_i) = nanmean(curdata(curloc==tissue_i));
                kstdMatrix(tissue_i) = nanstd(curdata(curloc==tissue_i));
             end
             lg4(cond_i) = plot(1:length(unique(smoothloc)), kmeanMatrix,...
                   mylinestyle{cond_i},...
                  'LineWidth', 2, 'Color', mycolors(cond_i,:));
             errorbar(1:length(unique(smoothloc)), kmeanMatrix, kstdMatrix,...
                  '.','LineWidth', 2, 'Color', mycolors(cond_i,:)) 
            axis square
            set(gca, 'XTick', 1:length(git_labels))
            set(gca, 'XTickLabel', git_labels)
            xlim([1,length(git_labels)])
            title('Normalized to max')
         end
         sgtitle(meier_mets{met_i});
         legend(lg1,meier_condition_unique, 'location', 'eastoutside');
         legend(lg2,meier_condition_unique, 'location', 'eastoutside');
         legend(lg3,meier_condition_unique, 'location', 'eastoutside');
         legend(lg4,meier_condition_unique, 'location', 'eastoutside');
         print(gcf, '-vector', '-dpsc2', '-r600', '-append', '-bestfit',...
             meier_figure_file_name);
     end
     % save smoothed data
     meier_data_smooth(met_i,:) = smoothdata;
 end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% apply modelling framework to Meier data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sampleType_unique = {'SPF', 'GF'};
sampleDiet_unique = {'Chow1', 'Chow2'};

% calculate mean across replicates
smoothloc_unique = unique(smoothloc);
meier_smooth_mean = zeros(size(meier_data_smooth,1),...
    length(smoothloc_unique)*length(meier_condition_unique));
meier_smooth_mean_condition = cell(length(smoothloc_unique)*length(meier_condition_unique),1);
idx=1;
for cond_i = 1:length(meier_condition_unique)
    for loc_i = 1:length(smoothloc_unique)
        meier_smooth_mean(:,idx) = nanmean(meier_data_smooth(:,...
            ismember(smoothcond, meier_condition_unique{cond_i}) &...
            (smoothloc == smoothloc_unique(loc_i))),2);
        meier_smooth_mean_condition{idx} = [meier_condition_unique{cond_i},...
            '_', git_labels{loc_i}];
    idx = idx+1;
    end
end

meanConditions = [strcat(meier_smooth_mean_condition, '_Chow1');
                  strcat(meier_smooth_mean_condition, '_Chow2')];              
% duplicate since we have only one diet
meanMatrix = [meier_smooth_mean meier_smooth_mean];
meanMatrix_mets = meier_mets;

% prepare labels for plotting
[bb,aa] = ndgrid(sampleType_unique, sampleDiet_unique); 
condLabels = strcat(aa(:),'-', bb(:));

% make a table of input data to the modelling framework
meier_smooth_mean_table = array2table(meier_smooth_mean, ...
    'VariableNames', meier_smooth_mean_condition);
meier_smooth_mean_table = addvars(meier_smooth_mean_table,...
    meier_mets, 'Before',meier_smooth_mean_table.Properties.VariableNames{1});

writetable(meier_smooth_mean_table,...
    [outputFolder,...
    'tableS13_Meier_input_data.csv']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileNameFigure = [figureFolder...
    'figSX_meierInt_diag_modelSMOOTH_2LIcoefHost_1LIbact.ps'];

diag_plot_flag = 0; % diagnostic plotting flag
if diag_plot_flag
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
end


% % % selected_mets are metabolites detected along the GI tract
selected_mets = 1:length(meier_mets)-1; % last metabolite is not measured

met_gitfits = cell(length(selected_mets),1);
met_bestsols = cell(length(selected_mets),1);

for met_i = 1:length(selected_mets)
   
        cmpd_interest_idx = selected_mets(met_i);
%         % set volume to CV or GF/WT
%         if contains(meanMatrix_mets{cmpd_interest_idx}, '_CV')
%             volumeMatrix = volumeMatrix_CVR;
%         else
%             volumeMatrix = volumeMatrix_GF;
%         end
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
                                                        git_labels);
                idx = idx+1;
            end
        end
        
%         if use_volume_flag
%             % multiply by volume to get amounts
%             kmeanMatrix_joint = kmeanMatrix_joint.*volumeMatrix; 
%         end
        % normalize by max intensity
        kmeanMatrix_joint = kmeanMatrix_joint./max(max(kmeanMatrix_joint));
        kmeanMatrix_joint(isnan(kmeanMatrix_joint))=0;

        if nnz(kmeanMatrix_joint)
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
                print(fig, '-vector', '-dpsc2', '-r600', '-append', '-bestfit',...
                        fileNameFigure);
                clf('reset')
            end    
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print best solutions to files
filename = [outputFolder ...
            'table_model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_Meier'];
% create met_info object needed for the printing function
met_info = struct;
met_info.CompoundID = meier_mets;
met_info.CompoundName = meier_mets;          
% print solutions to files
print_bestsol_to_files(met_info, met_bestsols, filename);
% test reading from file
%[met_info_bestsols_read, met_bestsols_IP_read] = read_bestsol_from_file(filename, 'IP');
%[met_info_bestsols_read, met_bestsols_LIwtPCC_read] = read_bestsol_from_file(filename, 'LI_PCC_within_high_total');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print all solutions to files
filename = [outputFolder ...
            'table_model_results_ALL_SMOOTH_raw_2LIcoefHost1LIcoefbact_Meier'];
% print solutions to files
print_allsols_to_files(met_info, met_gitfits, filename);
% test reading from file
%[met_info_read, met_gitfits_read] = read_allsols_from_files(filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% combine IP with LI PCC within high total PCC solutions
corrthreshold = 0.7; % threshold above which to combine solutions
inputfilename = [outputFolder ...
            'table_model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_Meier'];
sel_crit1 = 'IP';
sel_crit2 = 'LI_PCC_within_high_total';
[met_info_combined, met_bestsol_combined] = ...
                combine_bestsols_from_file(inputfilename, sel_crit1, sel_crit2,...
                                           corrthreshold);

% save combined best solution to file
outputfilename = [outputFolder ...
            'table_model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_Meier_combined_',...
            sel_crit1, '_', sel_crit2];
print_combined_bestsol_to_files(met_info_combined, met_bestsol_combined,...
                                outputfilename);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot correlation of restored and original data
% calculate differentce in corr distrbutions
filename = [figureFolder,...
    'Fig4a_histogram_MAXcorr_model_2LIcoefHost1LIcoefbact_meier_all_', sel_crit1, '_', sel_crit2];
plot_gitfit_model_corr(met_gitfits, filename)

% plot correlation of restored and original data for the best solution
% calculate differentce in corr distrbutions
filename = [figureFolder,...
    'Fig4a_histogram_MAXcorr_model_2LIcoefHost1LIcoefbact_meier_best_combined_', sel_crit1, '_', sel_crit2];
plot_bestfit_model_corr(met_bestsol_combined, met_gitfits, filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot data and restored intensities and model coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting file name
fileNameprofiles = [figureFolder, 'fig_profiles_combined_', sel_crit1, '_', sel_crit2,...
                        '_modelSMOOTH_2LIcoefHost_1LIbact_MeierDataNORMsmoothbysum.ps'];
     
gitfit_bestsol_combined_plot(met_info_combined, met_bestsol_combined,...
    fileNameprofiles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cluster only metabolites with corr>0.7
% uncomment if using whole metabolite table
clustidx = 1:size(met_bestsol_combined.x,1); % for testing purposed plot all 
clustdata = met_bestsol_combined.x(clustidx,2:end)';
clustrows = met_bestsol_combined.coefvalues(2:end);
clustcols = met_info_combined.CompoundName; 
% filter by resiprocal corr
clustcorr = met_bestsol_combined.x_sel_CorrRev(clustidx);
clustcorrLI = met_bestsol_combined.x_sel_CorrRevLI(clustidx);
clustdata = clustdata(:,clustcorr>=0.7);
clustidx = clustidx(clustcorr>=0.7);
clustcorrLI = clustcorrLI(clustcorr>=0.7);
clustcols = clustcols(clustcorr>=0.7);

for i=1:size(clustdata,2)
    clustdata(:,i) = clustdata(:,i)/max(abs(clustdata(:,i)));
end
clustnan = isnan(sum(clustdata,1));
clustdata(:,clustnan) = [];
clustidx(clustnan)=[];
% for high bacterial coefs, add additional filter for LI corr
remove_sols = (abs(clustdata(4,:))>=0.5)' & (clustcorrLI<0.7);
clustdata(:,remove_sols) = [];
clustidx(remove_sols)=[];
clustcols(remove_sols)=[];

% flip clustdata to have SI coefs on top
clustdata = flipud(clustdata);
clustrows = fliplr(clustrows);

clustdist ='cityblock';%'correlation';%'euclidean';%
cgo = clustergram(clustdata,...
            'RowLabels', clustrows,...
            'ColumnLabels', clustcols,...
            'ColumnPdist',clustdist,...
            'RowPdist', clustdist,...
            'DisplayRange', 1,...
            'colormap', redbluecmap,...
            'Cluster', 'row'); % cluster only columns and not coefficients

% add colorbar
cbButton = findall(0,'tag','HMInsertColorbar');
ccb = get(cbButton,'ClickedCallback');
set(cbButton,'State','on')
ccb{1}(cbButton,[],ccb{2})
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Manually add colorbar to clustergram and print to figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = cgo.plot;
orient landscape

print(gcf, '-vector', '-dpdf', '-r600', '-bestfit',...
         [figureFolder,...
    'fig_clustergram_2LIhos1LIbact_model_coefs_MeierData_combined',...
    sel_crit1, '_', sel_crit2, 'strictbactclass_0_7_smoothbysum'])
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
