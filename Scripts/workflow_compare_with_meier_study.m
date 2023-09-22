% Plot metabolites common with Meier study 
% across tissues with model coefficients

% call script defining file dependenciesand global variables
addpath(genpath('.\'));
add_global_and_file_dependencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data requirements: 
% 'metabolites_allions_combined_formulas_with_metabolite_filters_spatial100clusters_with_mean.csv'
% 'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_allions.csv'
% 'model_results_SMOOTH_normbyabsmax_reciprocal_problem_allions.csv'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output:
% Figures:
% 'fig3ade_profiles_figureselected_modelSMOOTH_2LIcoefHost_1LIbact.ps'; 	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read annotation from file with Meier IDs
annotationTableMeier = readtable([inputFolder ...
    'metabolites_allions_combined_formulas_with_metabolite_filters_meier.csv']);

%select_mets = cellfun(@(x) ~isempty(x), annotationTableMeier.Meier_compound);
select_mets = annotationTableMeier.MetaboliteFilter==1;

met_info.MZ =annotationTableMeier.MZ(select_mets);
met_info.RT =annotationTableMeier.RT(select_mets);
met_info.CompoundID = annotationTableMeier.CompoundID(select_mets);
met_info.CompoundName = annotationTableMeier.Meier_compound(select_mets);
% remove spaces and capitalize first letters
for i=1:length(met_info.CompoundName)
    str=lower(met_info.CompoundName{i});
    idx=regexp([' ' str],'(?<=\s+)\S','start')-1;
    str(idx)=upper(str(idx));
    str = strrep(str, ' ','');
    str = strrep(str, '-','_');
    met_info.CompoundName{i} = str;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modelingResults = readtable([resultsFolder ...
            'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_allions.csv']);
x_met_smooth = modelingResults{:, width(modelingResults)-8:end};
coefvalues = modelingResults.Properties.VariableNames(width(modelingResults)-8:end);

met_bestsol_DC.coefvalues = coefvalues;
met_bestsol_DC.x = x_met_smooth(select_mets,:);
%met_bestsol_DC.ReciprocalCorr = modelingResults.ReciprocalCorr(select_mets);
met_bestsol_DC.x_sel_CorrRev = cellfun(@(x) str2double(x),...
    modelingResults.ReciprocalCorr(select_mets));
met_bestsol_DC.x_sel_CorrRevSI = modelingResults.ReciprocalCorr(select_mets);
met_bestsol_DC.x_sel_CorrRevLI = modelingResults.ReciprocalCorr(select_mets);
met_bestsol_DC.modelname = 'Old_IP_DC';
met_bestsol_DC.selection_criterion = 'IP';


modelingResults = readtable([outputFolder...
            'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_allions_with_CVR.csv']);
x_met_smooth_CVR = modelingResults{:, width(modelingResults)-8:end};
coefvalues_CVR = modelingResults.Properties.VariableNames(width(modelingResults)-8:end);

met_bestsol_CVR.coefvalues = coefvalues_CVR;
met_bestsol_CVR.x = x_met_smooth_CVR(select_mets,:);
met_bestsol_CVR.x_sel_CorrRev = modelingResults.ReciprocalCorr(select_mets);
met_bestsol_CVR.x_sel_CorrRevSI = modelingResults.ReciprocalCorr(select_mets);
met_bestsol_CVR.x_sel_CorrRevLI = modelingResults.ReciprocalCorr(select_mets);
met_bestsol_CVR.modelname = 'Old_IP_CVR';
met_bestsol_CVR.selection_criterion = 'IP';

% load data and restored data from file
% save model results to file - reciprocal data restoration
modelData = readtable([resultsFolder...
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
% get correlations calculated with reverse problem
if isnumeric(modelData.ReciprocalCorr(1))
    x_data_corr = modelData.ReciprocalCorr;
else
    % it is not numeric, probably contains NaN - convert to numeric
    x_data_corr = cellfun(@(x) str2double(x), modelData.ReciprocalCorr);
end

met_bestsol_DC.kmeanMatrix_joint_names = modelData_data.Properties.VariableNames(1:24);
met_bestsol_DC.kmeanMatrix_joint_orig = modelData_data{select_mets, 1:24};
met_bestsol_DC.x_sel_dataR = modelData_data{select_mets, 25:48};
%%%%%%%%%%%%%%%%%%
% get CVR modelling results
modelData_CVR = readtable([outputFolder...
    'model_results_SMOOTH_normbyabsmax_reciprocal_problem_allions_with_CVR.csv']);
modelData_data_CVR = modelData_CVR(:, 9:end);
modelData_orig_CVR = modelData_data_CVR{:, cellfun(@(x) ~(contains(x, 'Recip') |...
                                             contains(x, 'Random') ),...
                                             modelData_data_CVR.Properties.VariableNames)};
modelData_orig_cols_CVR = modelData_data_CVR.Properties.VariableNames(cellfun(@(x) ~(contains(x, 'Recip') |...
                                             contains(x, 'Random') ),...
                                             modelData_data_CVR.Properties.VariableNames));
modelData_recip_CVR = modelData_data_CVR{:, cellfun(@(x) contains(x, 'Recip'),...
                                             modelData_data_CVR.Properties.VariableNames)};
% get correlations calculated with reverse problem
if isnumeric(modelData_CVR.ReciprocalCorr(1))
    x_data_corr_CVR = modelData_CVR.ReciprocalCorr;
else
    % it is not numeric, probably contains NaN - convert to numeric
    x_data_corr_CVR = cellfun(@(x) str2double(x), modelData_CVR.ReciprocalCorr);
end

met_bestsol_CVR.kmeanMatrix_joint_names = modelData_data_CVR.Properties.VariableNames(1:24);
met_bestsol_CVR.kmeanMatrix_joint_orig = modelData_data_CVR{select_mets, 1:24};
met_bestsol_CVR.x_sel_dataR = modelData_data_CVR{select_mets, 25:48};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read meier modelling results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = '.\ProcessedData\model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_Meier';
sel_crit = 'LI PCC within high total';%'IP';%
[met_info_meier, met_bestsol_meier] = read_bestsol_from_file(filename,...
                                                sel_crit);
%edit meier names to remove x in front of metabolites
for i=1:length(met_info_meier.CompoundName)
    curname = met_info_meier.CompoundName{i};
    if curname(1)=='x'
        if isstrprop(curname(2),'digit')
            curname = curname(2:end);
            met_info_meier.CompoundName{i} = curname;
        end
    end
end
met_bestsol_meier.modelname = 'Meier_LIPCCwT';%'Meier_IP';

compare_two_model_results(met_info, met_bestsol_DC,...
                          met_info_meier, met_bestsol_meier,...
                          figureFolder);

[confmat, classlabels] = compare_two_model_results(met_info, met_bestsol_DC,...
                          met_info, met_bestsol_CVR,...
                          figureFolder);

compare_two_model_reciprocal_data(met_info, met_bestsol_DC,...
                          met_info_meier, met_bestsol_meier,...
                          figureFolder);
                                                         
x_met_smooth_meier = met_bestsol_meier.sol_coefs;
coefvalues_meier = met_bestsol_meier.coefvalues;

% % load data and restored data from file
% % save model results to file - reciprocal data restoration
% modelData_meier = readtable([outputFolder...
%     'model_results_SMOOTH_normbyabsmax_reciprocal_problem_MeierData_smoothbysum.csv']);
% % modelData = readtable([outputFolder...
% %     'model_results_SMOOTH_normbyabsmax_reciprocal_problem_allions_with_CVR.csv']);
% modelData_data_meier = modelData_meier(:, 3:end);
% modelData_meier.CompoundName = cellfun(@(x) strrep(x, 'x1','1'), modelData_meier.CompoundName, 'unif',0);
% modelData_meier.CompoundName = cellfun(@(x) strrep(x, 'x2','2'), modelData_meier.CompoundName, 'unif',0);
% modelData_meier.CompoundName = cellfun(@(x) strrep(x, 'x3','3'), modelData_meier.CompoundName, 'unif',0);
% modelData_meier.CompoundName = cellfun(@(x) strrep(x, 'x4','4'), modelData_meier.CompoundName, 'unif',0);
% modelData_meier.CompoundName = cellfun(@(x) strrep(x, 'x5','5'), modelData_meier.CompoundName, 'unif',0);
% 
% 
% modelData_orig_meier = modelData_data_meier{:, cellfun(@(x) ~(contains(x, 'Recip') |...
%                                              contains(x, 'Random') ),...
%                                              modelData_data_meier.Properties.VariableNames)};
% modelData_orig_cols_meier = modelData_data_meier.Properties.VariableNames(cellfun(@(x) ~(contains(x, 'Recip') |...
%                                              contains(x, 'Random') ),...
%                                              modelData_data_meier.Properties.VariableNames));
% modelData_recip_meier = modelData_data_meier{:, cellfun(@(x) contains(x, 'Recip'),...
%                                              modelData_data_meier.Properties.VariableNames)};
% % get correlations calculated with reverse problem
% if isnumeric(modelData_meier.ReciprocalCorr(1))
%     x_data_corr_meier = modelData_meier.ReciprocalCorr;
% else
%     % it is not numeric, probably contains NaN - convert to numeric
%     x_data_corr_meier = cellfun(@(x) str2double(x), modelData_meier.ReciprocalCorr);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define colora and GIT section names for plotting
mycolors = [0 115 178;... %dark blue
            204 227 240;...%light blue
            211 96 39;... %dark orange
            246 223 212]/256;%light orange
git_labels = {'Du', 'Je', 'Il', 'Cec', 'Col', 'Fec'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting file name
%fileNameprofiles = 'fig3ade_profiles_figureselected_modelSMOOTH_2LIcoefHost_1LIbact_selectedCVR.ps';
%fileNameprofiles = 'fig_profiles_MeierComparison_smoothbysum_modelSMOOTH_2LIcoefHost_1LIbact_nonunique.ps';
%fileNameprofiles = 'fig_profiles_MeierComparison_smoothbysum_modelSMOOTH_2LIcoefHost_1LIbact.ps';
%fileNameprofiles = 'fig_profiles_MeierComparison_smoothbysum_modelSMOOTH_2LIcoefHost_1LIbact_corr07.ps';
%fileNameprofiles = 'fig_profiles_MeierComparison_smoothbysum_modelSMOOTH_2LIcoefHost_1LIbact_corr07_CVR.ps';
%fileNameprofiles = 'fig_profiles_MeierComparison_smoothbysum_modelSMOOTH_2LIcoefHost_1LIbact_DC_CVR.ps';
fileNameprofiles = 'fig_profiles_MeierComparison_smoothbysum_modelSMOOTH_2LIcoefHost_1LIbact_DC_CVR_nonunique.ps';

% select ions to plot by ion MZ
% targetMZ = [147.053; 74.037; 287.210;...]; %glutamate, propionate, l-octanoylcarnitine 
%             499.297; 125.015; 131.058; 226.095;...]; %taurodeoxycholate, taurine, 5-aminolevulinate, porphobilonogen
%             181.074; 468.272; 483.301; 245.163; 576.512; 430.345;... %tyrosine, 3,17-androstnediol glucuronide, taurolitocholate, isovalerylcarnitine, cohibin, 4a-carboxy-5a-cholesta-8-en-3b-ol
%             386.355; 99.068;... % 5alpha-cholestan-3-one, hydroxy-methylbutanitrile
%             119.058]; %threonine
% % find annotated compound with this MZ
% % for which modelling results correlate >0.7 with original
% compoundsInterest = find((annotationTableSpatialClusters.MetaboliteFilter>0) &...
%                (x_data_corr>=0.7) &...
%                (arrayfun(@(x) sum(abs(x-targetMZ)<=0.001),...
%                               annotationTableSpatialClusters.MZ)>0));
compoundsInterest = find(cellfun(@(x) ~isempty(x),...
    annotationTableSpatialClusters.Meier_compound) &...
    (annotationTableSpatialClusters.MetaboliteFilter==1));
% compoundsInterest = find(cellfun(@(x) ~isempty(x),...
%     annotationTableSpatialClusters.Meier_compound));

curData_cols = reshape(modelData_orig_cols, 6, 4);

fprintf('Plotting profiles for %d metabolites\n',length(compoundsInterest));
plotflag = 1;

if plotflag
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
end
% save all coefs for original and Meier data and correlation coefs
compare_coefs_orig = zeros(size(x_met_smooth,2)-1, length(compoundsInterest));
compare_PCC_orig = zeros(1, length(compoundsInterest));
% also add profiles
compare_profiles_orig = zeros(size(modelData_orig, 2), length(compoundsInterest));
% CVR
compare_coefs_orig_CVR = zeros(size(x_met_smooth,2)-1, length(compoundsInterest));
compare_PCC_orig_CVR = zeros(1, length(compoundsInterest));
% also add profiles
compare_profiles_CVR = zeros(size(modelData_orig, 2), length(compoundsInterest));
% Meier
compare_coefs_meier = zeros(size(x_met_smooth,2)-1, length(compoundsInterest));
compare_PCC_meier = zeros(1, length(compoundsInterest));
compare_cpd_names = cell(1, length(compoundsInterest));
% also add profiles
compare_profiles_meier = zeros(size(modelData_orig, 2), length(compoundsInterest));

for cpdix=1:length(compoundsInterest)
    
    testidx = compoundsInterest(cpdix);

    testann = annotationTableSpatialClusters.CompoundName(testidx,:);
    testannID = annotationTableSpatialClusters.CompoundID(testidx,:);
    testmz = annotationTableSpatialClusters.MZ(testidx,:);
    testrt = annotationTableSpatialClusters.RT(testidx,:);
    testmethod = annotationTableSpatialClusters.Method(testidx,:);
    testmode = annotationTableSpatialClusters.Mode(testidx,:);
    testfilter = annotationTableSpatialClusters.MetaboliteFilter(testidx,:);
    
    % gete Meier idx
    meiercompound = annotationTableSpatialClusters.Meier_compound(testidx);
    meieridx = find(ismember(cellfun(@(x) lower(x),modelData_meier.CompoundName, 'unif',0),...
        lower(strrep( strrep(meiercompound, '-','_'),' ',''))));

 
    spx=3;
    spy=3;
    spidx = 1;
    %curmat = zeros(4,6);
    % community data
    cur_data = reshape(modelData_orig(testidx,:), 4, 6)';
    cur_rdata = reshape(modelData_recip(testidx,:), 6, 4);
    %normalize rdata to max
    cur_rdata = (cur_rdata-0.8);
    cur_rdata = cur_rdata/max(max(cur_rdata));
    
    % CVR data
    cur_data_CVR = reshape(modelData_orig_CVR(testidx,:), 4, 6)';
    cur_rdata_CVR = reshape(modelData_recip_CVR(testidx,:), 6, 4);
    %normalize rdata to max
    cur_rdata_CVR = (cur_rdata_CVR-0.8);
    cur_rdata_CVR = cur_rdata_CVR/max(max(cur_rdata_CVR));
  
    curcoefs = x_met_smooth(testidx, 2:end);
    curcoefs_CVR = x_met_smooth_CVR(testidx, 2:end);
  
    
    if (((x_data_corr(testidx)>=0.7) || (x_data_corr_CVR(testidx)>=0.7)) &&...
            (x_data_corr_meier(meieridx)>=0.7))
        plotflag=1;
    else
        plotflag=0;
    end

    % do not plot
    plotflag=0;
    
    if plotflag
        legend_entries = cell(4,1);
        for i = 1:size(cur_data,2)

            subplot(spx,spy,1)
            hold on
            h(i) = plot(cur_data(:,i),...
                 'LineWidth', 2,...
                 'Color', mycolors(i,:));
            ylabel('Original normbymax')

            subplot(spx,spy,2);
            hold on
            plot((cur_rdata(:,i)),...
                     'LineWidth', 2,...
                     'Color', mycolors(i,:));
            ylabel('Restored normbymax')

        end
        title(sprintf('PCC=%.2f',...
                            x_data_corr(testidx)))

        legend(h, curData_cols(1,:))%, 'Location', 'bestoutside');

        subplot(spx,spy,3);
   
        barh(curcoefs./max(abs(curcoefs)))
        set(gca, 'YTick', 1:length(curcoefs));
        set(gca, 'YTickLabel', coefvalues(2:end));
        set(gca, 'YDir','reverse')
        ylim([0.5 length(curcoefs)+0.5])
        xlim([-1 1]);
        axis square
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plot CVR
        legend_entries = cell(4,1);
        for i = 1:size(cur_data_CVR,2)

            subplot(spx,spy,4)
            hold on
            h(i) = plot(cur_data_CVR(:,i),...
                 'LineWidth', 2,...
                 'Color', mycolors(i,:));
            ylabel('Original normbymax')

            subplot(spx,spy,5);
            hold on
            plot((cur_rdata_CVR(:,i)),...
                     'LineWidth', 2,...
                     'Color', mycolors(i,:));
            ylabel('Restored normbymax')

        end
        title(sprintf('PCC=%.2f',...
                            x_data_corr_CVR(testidx)))

        legend(h, curData_cols(1,:))%, 'Location', 'bestoutside');

        subplot(spx,spy,6);
   
        barh(curcoefs_CVR./max(abs(curcoefs_CVR)))
        set(gca, 'YTick', 1:length(curcoefs_CVR));
        set(gca, 'YTickLabel', coefvalues_CVR(2:end));
        set(gca, 'YDir','reverse')
        ylim([0.5 length(curcoefs_CVR)+0.5])
        xlim([-1 1]);
        axis square
    end
    
    % save current coefs and corr coef
    compare_coefs_orig(:,cpdix) = curcoefs./max(abs(curcoefs));
    compare_PCC_orig(cpdix) = x_data_corr(testidx);
    compare_profiles_orig(:, cpdix) = modelData_orig(testidx,:);
    % save current coefs and corr coef for CVR
    compare_coefs_orig_CVR(:,cpdix) = curcoefs_CVR./max(abs(curcoefs_CVR));
    compare_PCC_orig_CVR(cpdix) = x_data_corr_CVR(testidx);
    compare_profiles_CVR(:, cpdix) = modelData_orig_CVR(testidx,:);
   
    % plot Meier models
  
    curmat = zeros(4,6);
    cur_data = reshape(modelData_orig_meier(meieridx,:), 4, 6)';
    cur_rdata = reshape(modelData_recip_meier(meieridx,:), 6, 4);
    
    %normalize rdata to max
    cur_rdata = (cur_rdata-0.8);
    cur_rdata = cur_rdata/max(max(cur_rdata));
    curcoefs = x_met_smooth_meier(meieridx, 2:end);
  
    if plotflag
        legend_entries = cell(4,1);
        for i = 1:size(cur_data,2)

            subplot(spx,spy,7)
            hold on
            h(i) = plot(cur_data(:,i),...
                 'LineWidth', 2,...
                 'Color', mycolors(i,:));
            ylabel('Original normbymax')

            subplot(spx,spy,8);
            hold on
            plot((cur_rdata(:,i)),...
                     'LineWidth', 2,...
                     'Color', mycolors(i,:));
            ylabel('Restored normbymax')

        end
        title(sprintf('PCC=%.2f',...
                            x_data_corr_meier(meieridx)))

        legend(h, curData_cols(1,:))%, 'Location', 'bestoutside');

        subplot(spx,spy,9);
        barh(curcoefs./max(abs(curcoefs)))
        set(gca, 'YTick', 1:length(curcoefs));
        set(gca, 'YTickLabel', coefvalues(2:end));
        set(gca, 'YDir','reverse')
        ylim([0.5 length(curcoefs)+0.5])
        xlim([-1 1]);
        axis square

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for spi = 1:(spx*spy)
            if mod(spi, 3)~=0
                subplot(spx,spy,spi)
                set(gca, 'XTick', 1:6)
                xlim([1 6])
                ylim([0 1])
                set(gca, 'XTick', 1:length(git_labels))
                set(gca, 'XTickLabel', git_labels)

                axis square
            end
        end
        spt = suptitle({meiercompound{1},...
                        sprintf('MZ=%.3f RT=%.2f %s %d filter=%d',testmz(1),...
                        testrt(1), testmethod{1},testmode(1),testfilter(1)),...
                                            testannID{1},...
                                            testann{1}});
        set(spt,'FontSize',8,'FontWeight','normal')
        orient landscape
        %print to figure
        print(gcf, '-painters', '-dpsc2', '-r600', '-append', '-bestfit',...
                [figureFolder,...
                 fileNameprofiles])
        clf('reset')
    end
    % save meier coefs
    compare_coefs_meier(:,cpdix) = curcoefs./max(abs(curcoefs));
    compare_PCC_meier(cpdix) = x_data_corr_meier(meieridx);
    compare_profiles_meier(:, cpdix) = modelData_orig_meier(meieridx,:);
    compare_cpd_names{cpdix} = meiercompound{1};

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare Meier and original and CVR coefs and profiles
coef_corr_orig_meier_coefs = zeros(size(compare_coefs_orig,2),1);
coef_corr_orig_meier_profiles = zeros(size(compare_coefs_orig,2),1);
coef_corr_CVR_meier_coefs = zeros(size(compare_coefs_orig,2),1);
coef_corr_CVR_meier_profiles = zeros(size(compare_coefs_orig,2),1);
coef_corr_orig_CVR_coefs = zeros(size(compare_coefs_orig,2),1);
coef_corr_orig_CVR_profiles = zeros(size(compare_coefs_orig,2),1);
%coef_corrS = zeros(size(compare_coefs_orig,2),1);
for i=1:size(compare_coefs_orig,2)
    coef_corr_orig_meier_coefs(i) = corr(compare_coefs_orig(:,i),...
                                         compare_coefs_meier(:,i));
    coef_corr_orig_meier_profiles(i) = corr(compare_profiles_orig(:,i),...
                                         compare_profiles_meier(:,i));
    %coef_corrS(i) = corr(compare_coefs_orig(:,i), compare_coefs_meier(:,i),...
    %    'type', 'Spearman');
    coef_corr_CVR_meier_coefs(i) = corr(compare_coefs_orig_CVR(:,i),...
                                         compare_coefs_meier(:,i));
    coef_corr_CVR_meier_profiles(i) = corr(compare_profiles_CVR(:,i),...
                                         compare_profiles_meier(:,i));

    coef_corr_orig_CVR_coefs(i) = corr(compare_coefs_orig(:,i),...
                                         compare_coefs_orig_CVR(:,i));
    coef_corr_orig_CVR_profiles(i) = corr(compare_profiles_orig(:,i),...
                                         compare_profiles_CVR(:,i));

end

figure
histogram(coef_corr,20)
hold on
select_corr = (compare_PCC_meier>0.7) & (compare_PCC_orig>0.7);
histogram(coef_corr(select_corr),20)

figure
scatter(compare_coefs_orig(end,:), compare_coefs_meier(end,:))

figure
scatter(compare_coefs_orig(end,select_corr), compare_coefs_meier(end,select_corr))

figure
scatter(compare_coefs_orig(end,:), compare_coefs_orig_CVR(end,:))

figure
scatter(compare_coefs_orig_CVR(end,:), compare_coefs_meier(end,:))

nnz(((compare_coefs_orig(end,select_corr)<0) | ...
     (compare_coefs_orig(end-1,select_corr)<0)) & ...
    (compare_coefs_meier(end,select_corr)<0) )

nnz(((compare_coefs_orig(end,select_corr)>0) | ...
     (compare_coefs_orig(end-1,select_corr)>0)) & ...
    (compare_coefs_meier(end,select_corr)>0) )

threshold = 0.25;
nnz((((compare_coefs_orig(end,select_corr)<-threshold) | ...
     (compare_coefs_orig(end-1,select_corr)<-threshold)) & ...
    (compare_coefs_meier(end,select_corr)<-threshold) ) | ...
    (((compare_coefs_orig(end,select_corr)>threshold) | ...
     (compare_coefs_orig(end-1,select_corr)>threshold)) & ...
    (compare_coefs_meier(end,select_corr)>threshold) ))

nnz((((compare_coefs_orig(end,select_corr)<-threshold) | ...
     (compare_coefs_orig(end-1,select_corr)<-threshold)) | ...
    (((compare_coefs_orig(end,select_corr)>threshold) | ...
     (compare_coefs_orig(end-1,select_corr)>threshold)) )))

nnz((compare_coefs_meier(end,select_corr)<-threshold) | ...
    (compare_coefs_meier(end,select_corr)>threshold) )

nnz((((compare_coefs_orig(1,select_corr)<-threshold) | ...
     (compare_coefs_orig(4,select_corr)<-threshold)) & ...
    (compare_coefs_meier(1,select_corr)<-threshold) ) | ...
    (((compare_coefs_orig(1,select_corr)>threshold) | ...
     (compare_coefs_orig(4,select_corr)>threshold)) & ...
    (compare_coefs_meier(1,select_corr)>threshold) ))

nnz((((compare_coefs_orig(1,select_corr)<-threshold) | ...
     (compare_coefs_orig(4,select_corr)<-threshold)) | ...
    (((compare_coefs_orig(1,select_corr)>threshold) | ...
     (compare_coefs_orig(4,select_corr)>threshold)) )))

select_DC = compare_coefs_orig(:, select_corr);
select_CVR = compare_coefs_orig_CVR(:, select_corr);
select_meier = compare_coefs_meier(:, select_corr);

dist_DC_CVR = zeros(nnz(select_corr),1); 
dist_DC_meier = zeros(nnz(select_corr),1); 
dist_CVR_meier = zeros(nnz(select_corr),1); 
for i=1:nnz(select_corr)
    G1 = select_DC(:,i);
    G2 = select_CVR(:,i);
    dist_DC_CVR(i) = sqrt(sum((G1 - G2) .^ 2));
    G1 = select_DC(:,i);
    G2 = select_meier(:,i);
    dist_DC_meier(i) = sqrt(sum((G1 - G2) .^ 2));
    G1 = select_CVR(:,i);
    G2 = select_meier(:,i);
    dist_CVR_meier(i) = sqrt(sum((G1 - G2) .^ 2));
end

figure
hold on
nbin=10;
histogram(dist_DC_CVR, nbin)
histogram(dist_DC_meier, nbin)
histogram(dist_CVR_meier, nbin)
    
    
    
    