% Plot all annotated metabolites across tissues with model coefficients

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
% 'fig3aef_profiles_figureselected_modelSMOOTH_2LIcoefHost_1LIbact.ps'; 	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read annotation from file
% annotationTableSpatialClusters = readtable([inputFolder ...
%     'metabolites_allions_combined_formulas_with_metabolite_filters_spatial100clusters_with_mean.csv']);
annotationTableSpatialClusters = readtable([outputFolder ...
    'metabolites_allions_combined_formulas_with_metabolite_filters_spatial100clusters_with_mean_with_CVR_0925.csv']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modelingResults = readtable([resultsFolder ...
%             'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_allions.csv']);
%modelingResults = readtable([outputFolder...
%            'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_allions_with_CVR.csv']);
modelingResults = readtable([outputFolder...
    'table_model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_DC_combined_IP_LI_PCC_within_high_total.csv']);
x_met_smooth = modelingResults{:, width(modelingResults)-8:end};
coefvalues = modelingResults.Properties.VariableNames(width(modelingResults)-8:end);

% load data and restored data from file
% save model results to file - reciprocal data restoration
% modelData = readtable([resultsFolder...
%     'model_results_SMOOTH_normbyabsmax_reciprocal_problem_allions.csv']);
%modelData = readtable([outputFolder...
%    'model_results_SMOOTH_normbyabsmax_reciprocal_problem_allions_with_CVR.csv']);
modelData = readtable([outputFolder...
    'table_model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_DC_combined_IP_LI_PCC_within_high_total_reciprocal.csv']);
modelData_data = modelData(:, 12:end);
modelData_orig = modelData_data{:, cellfun(@(x) ~(contains(x, 'Recip') |...
                                             contains(x, 'Random') ),...
                                             modelData_data.Properties.VariableNames)};
modelData_orig_cols = modelData_data.Properties.VariableNames(cellfun(@(x) ~(contains(x, 'Recip') |...
                                             contains(x, 'Random') ),...
                                             modelData_data.Properties.VariableNames));
modelData_recip = modelData_data{:, cellfun(@(x) contains(x, 'Recip'),...
                                             modelData_data.Properties.VariableNames)};
modelData_recip_cols = modelData_data.Properties.VariableNames(cellfun(@(x) contains(x, 'Recip') ,...
                                             modelData_data.Properties.VariableNames));
% get correlations calculated with reverse problem
if isnumeric(modelData.ReciprocalCorr(1))
    x_data_corr = modelData.ReciprocalCorr;
else
    % it is not numeric, probably contains NaN - convert to numeric
    x_data_corr = cellfun(@(x) str2double(x), modelData.ReciprocalCorr);
end

% get names of spatial clusters                                         
spatialClusters_names = annotationTableSpatialClusters.Properties.VariableNames(...
    cellfun(@(x) contains(x, 'spatial_clust100'), annotationTableSpatialClusters.Properties.VariableNames));
spatialClusters = annotationTableSpatialClusters{:, spatialClusters_names};                                   
spatialClusters_names = cellfun(@(x) strrep(x, 'spatial_clust100_', ''), spatialClusters_names, 'unif', 0);

% remove CVR clusters
spatialClusters(:, cellfun(@(x) contains(x, 'CVR'), spatialClusters_names)) = [];
spatialClusters_names(cellfun(@(x) contains(x, 'CVR'), spatialClusters_names)) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define colors and GIT section names for plotting
mycolors = [0 115 178;... %dark blue
            0 115 178;... %dark blue 204 227 240;...%light blue
            211 96 39;... %dark orange
            211 96 39]/256; %dark orange % 246 223 212light orange
mylinestyles = {'-', '--', '-', '--'};

git_labels = {'Du', 'Je', 'Il', 'Cec', 'Col', 'Fec'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting file name
%fileNameprofiles = 'fig3ade_profiles_figureselected_modelSMOOTH_2LIcoefHost_1LIbact.ps';
%fileNameprofiles = 'fig3ade_profiles_figureselected_modelSMOOTH_2LIcoefHost_1LIbact_selectedCVR.ps';
fileNameprofiles = 'fig3aef_profiles_figureselected_modelSMOOTH_2LIcoefHost_1LIbact_DC_combined_IP_LI_PCC_within_high_total.ps';

% select ions to plot by ion MZ
targetMZ = [147.053; 74.037; 287.210;...]; %glutamate, propionate, l-octanoylcarnitine 
            499.297; 125.015; 131.058; 226.095;...]; %taurodeoxycholate, taurine, 5-aminolevulinate, porphobilonogen
            181.074; 468.272; 483.301; 245.163; 576.512; 430.345;... %tyrosine, 3,17-androstnediol glucuronide, taurolitocholate, isovalerylcarnitine, cohibin, 4a-carboxy-5a-cholesta-8-en-3b-ol
            384.339; 386.355; 99.068;... % Cholestenone;, 5alpha-cholestan-3-one, hydroxy-methylbutanitrile
            119.058; 138.043]; %threonine, urocanate  
% find annotated compound with this MZ
% for which modelling results correlate >0.7 with original
compoundsInterest = find((annotationTableSpatialClusters.MetaboliteFilter>0) &...
               (x_data_corr>=0.7) &...
               (arrayfun(@(x) sum(abs(x-targetMZ)<=0.001),...
                              annotationTableSpatialClusters.MZ)>0));

% %plot all annotated metabolites with high qialuty modelling results
% compoundsInterest = find((annotationTableSpatialClusters.MetaboliteFilter>0) &...
%                          (x_data_corr>=0.7));

           
curData_cols = reshape(modelData_orig_cols, 4, 6)';
fig = figure('units','normalized','outerposition',[0 0 1 1]);

fprintf('Plotting profiles for %d metabolites\n',length(compoundsInterest));
for cpdix=1:length(compoundsInterest)
    
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
    cur_rdata = reshape(modelData_recip(testidx,:), 4, 6)';
    
    %normalize rdata to max
    cur_rdata = (cur_rdata-0.8);
    cur_rdata = cur_rdata/max(max(cur_rdata));
    
    legend_entries = cell(4,1);
    for i = 1:size(cur_data,2)

        subplot(spx,spy,idx)
        hold on
        h(i) = plot(cur_data(:,i),...
             'LineWidth', 2,...
             'Color', mycolors(i,:),...
             'LineStyle', mylinestyles{i});
        ylabel('Original normbymax')

        subplot(spx,spy,idx+1);
        hold on
        plot((cur_rdata(:,i)),...
                 'LineWidth', 2,...
                 'Color', mycolors(i,:),...
                 'LineStyle', mylinestyles{i});
        ylabel('Restored normbymax')

    end
    title(sprintf('%s %s %s %s: %d %d %d %d PCC=%.2f',...
                        spatialClusters_names{:},...
                        spatialClusters(testidx,:),...
                        x_data_corr(testidx)),...
                        'Interpreter', 'none')

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
    spt = sgtitle({sprintf('MZ=%.3f',testmz(idx,1)),...
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
