%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot metabolite profiles across tissues

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output: 
% Figures:
% 'Fig4a_histogram_corr_model_2LIcoefHost1LIcoefbact_reversedata_annotated'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read metabolite data
metaboliteFolder = '.\metabolomics\ProcessedData\';
metaboliteData = readtable([inputFolder,...
    'metabolites_allions_combined_norm_intensity.csv']);
    
metaboliteFilters = readtable([inputFolder,...
    'metabolites_allions_combined_formulas_with_metabolite_filters.csv']);

% define colors for different nouse groups and tissue names
mycolors = [0 115 178;... %dark blue
            211 96 39;... %dark orange
            204 227 240;...%light blue
            246 223 212]/256;%light orange
git_labels = {'Du', 'Je', 'Il', 'Cec', 'Col', 'Fec'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot all annotated metabolites across tissues
% fileNameprofiles = 'profiles_annotated_mets_across_tissues_with_datapoints_updated_filter_std_selected_removeoutliers.ps';
% compoundsInterest = find((metaboliteFilters.MetaboliteFilter>0));

%%%
% fileNameprofiles = 'profiles_annotated_mets_across_tissues_with_datapoints_hier_clust_bact_and_system.ps';
% % from workflow_plot_and_analyze_modeling_results
% compoundsInterest = find( ( (hier_clust_attribution.HierarchicalClustGroup==10) |...
%                             (hier_clust_attribution.HierarchicalClustGroup==11) ) &...
%     (  (met_liver_changes(:,1)>0) | (met_liver_changes(:,2)>0) |...
%        (met_serum_changes(:,1)>0) | (met_serum_changes(:,2)>0) ));


% plot compounds with specified masses
% plotting file name
fileNameprofiles = 'fig2d_5fg_profiles_mets_across_tissues_with_datapoints_hier_clust_bact_and_system.ps'; 	
% select ions to plot by ion MZ
targetMZ = [147.053; 74.037; 287.210;... %glutamate, propionate, l-octanoylcarnitine 
            499.297; 125.015; 131.058; 226.095;... %aurodeoxycholate, taurine, 5-aminolevulinate, porphobilinogen
            181.074; 468.272; 483.301; 245.163; 576.512; 430.345;... %tyrosine, 3,17-androstnediol glucuronide, taurolitocholate, isovalerylcarnitine, cohibin, 4a-carboxy-5a-cholesta-8-en-3b-ol
            386.355; 99.068]; % 5alpha-cholestan-3-one, hydroxy-methylbutanitrile
% find annotated compound with this MZ
% for which modelling results correlate >0.7 with original
compoundsInterest = find((metaboliteFilters.MetaboliteFilter>0) &...
               (arrayfun(@(x) sum(abs(x-targetMZ)<=0.001),...
                              metaboliteFilters.MZ)>0));

% prepare legends for the plots
randparam = 0.25;
legend_entries = cell(length(sampleType_unique)*length(sampleDiet_unique),1);
idx=1;
for type_i = 1:length(sampleType_unique)
    for diet_i = 1:length(sampleDiet_unique)
        legend_entries{idx} = [sampleType_unique{type_i}, ' ', sampleDiet_unique{diet_i}];
        idx=idx+1;
    end
end

fig = figure('units','normalized','outerposition',[0 0 1 1]);
% flag whether two remove 1-2 ourliers out of 5
remove_outliers_flag = 1;

fprintf('Plotting profiles for %d metabolites\n',length(compoundsInterest));
for cpdix=1:length(compoundsInterest)
  
    testidx = compoundsInterest(cpdix);

    testann = metaboliteFilters.CompoundID(testidx);
    testmz = metaboliteFilters.MZ(testidx);
    testrt = metaboliteFilters.RT(testidx);
    testmethod = metaboliteFilters.Method(testidx);
    testmode = metaboliteFilters.Mode(testidx);
    testname = metaboliteFilters.CompoundName(testidx);

    spx=2;
    spy=2;
   
    coloridx = 1;
    idx=1;
    curmat = zeros(4,6);
    for type_i = 1:length(sampleType_unique)
        for diet_i = 1:length(sampleDiet_unique)
            kmeanMatrix = zeros(6,1);
            kstdMatrix = zeros(6,1);
            for tissue_i = 1:6
                selectCols = ismember(combinedTissues, sampleTissue_order{tissue_i}) &...
                             ismember(combinedDiet, sampleDiet_unique{diet_i}) &...
                             ismember(combinedType, sampleType_unique{type_i});
                curMatrix = combinedIntensitiesNorm(testidx,selectCols);
                         
                % remove two 5000
                if remove_outliers_flag
                    indx5000 = find(curMatrix==5000);
                    curMatrix(indx5000(1:min(2,length(indx5000))))=[];
                end
                
                kmeanMatrix(tissue_i) = nanmean(curMatrix);
                kstdMatrix(tissue_i) = nanstd(curMatrix);

                subplot(spx,spy,1)
            
                hold on
                plot(tissue_i*ones(size(curMatrix)) + randparam*(rand(size(curMatrix))-0.5),...
                     curMatrix,...
                                 'o', 'MarkerSize', 5, 'MarkerEdgeColor', mycolors(coloridx,:),...
                                 'MarkerFaceColor', mycolors(coloridx,:));
            end
            plot(1:6, kmeanMatrix, 'LineWidth', 2, 'Color', mycolors(coloridx,:));
            errorbar(1:6, kmeanMatrix, kstdMatrix, '.','LineWidth', 2, 'Color', mycolors(coloridx,:))
            set(gca, 'XTick', 1:length(git_labels))
            set(gca, 'XTickLabel', git_labels)
            ylabel('Quantile normalized')
            title('GIT')
            
            subplot(spx,spy,2)
            hold on
            h(coloridx) = plot(1:6, kmeanMatrix, 'LineWidth', 2, 'Color', mycolors(coloridx,:));
            
            % liver
            selectCols = ismember(combinedTissues, sampleTissue_order{7}) &...
                             ismember(combinedDiet, sampleDiet_unique{diet_i}) &...
                             ismember(combinedType, sampleType_unique{type_i});
            curMatrix = combinedIntensitiesNorm(testidx,selectCols);
                
            if remove_outliers_flag
                    indx5000 = find(curMatrix==5000);
                    curMatrix(indx5000(1:min(2,length(indx5000))))=[];
            end
                
            subplot(spx,spy,3);
            hold on
            plot((coloridx+1)*ones(size(curMatrix))+ randparam*(rand(size(curMatrix))-0.5),...
                curMatrix,...
                                 'o', 'MarkerSize', 5, 'MarkerEdgeColor', mycolors(coloridx,:),...
                                 'MarkerFaceColor', mycolors(coloridx,:));
            plot(coloridx+1, mean(curMatrix),...
                 'LineWidth', 2,...
                 'Color', mycolors(coloridx,:));
            errorbar(coloridx+1, nanmean(curMatrix), nanstd(curMatrix), '.',...
                     'LineWidth', 2,...
                     'Color', mycolors(coloridx,:))
            ylabel('Quantile normalized')
            title('Liver')
            
            % serum
            selectCols = ismember(combinedTissues, sampleTissue_order{8}) &...
                             ismember(combinedDiet, sampleDiet_unique{diet_i}) &...
                             ismember(combinedType, sampleType_unique{type_i});
            curMatrix = combinedIntensitiesNorm(testidx,selectCols);
           
            if remove_outliers_flag
                    indx5000 = find(curMatrix==5000);
                    curMatrix(indx5000(1:min(2,length(indx5000))))=[];
            end
                
            subplot(spx,spy,4);
            hold on
            plot((coloridx+1)*ones(size(curMatrix))+ randparam*(rand(size(curMatrix))-0.5),...
                curMatrix,...
                                 'o', 'MarkerSize', 5, 'MarkerEdgeColor', mycolors(coloridx,:),...
                                 'MarkerFaceColor', mycolors(coloridx,:));
            plot(coloridx+1, mean(curMatrix),...
                 'LineWidth', 2,...
                 'Color', mycolors(coloridx,:));
            errorbar(coloridx+1, nanmean(curMatrix), nanstd(curMatrix), '.',...
                     'LineWidth', 2,...
                     'Color', mycolors(coloridx,:))
            ylabel('Quantile normalized')
            title('Serum')
            
            %disp(kmeanMatrix(testidx(idx),:))
            coloridx = coloridx+1;
        end
    end

    subplot(spx, spy,2)
    legend(h, legend_entries)
    
    for spi = 1:spx*spy
        cursp = subplot(spx,spy,spi);
        set(cursp, 'XTick', 1:6)
        xlim([1 6])
        cur_ylim = get(cursp, 'ylim');
        if cur_ylim(2)==5001
            cur_ylim(2)=10000;
        end
        ylim([0 cur_ylim(2)])
        axis square
    end
    spt = suptitle(sprintf('MZ=%.3f, RT=%.2f %s %d\n%s\n%s',...
                                        testmz,testrt,testmethod{1}, testmode,...
                                        testann{1},...
                                        testname{1}));
    set(spt,'FontSize',8,'FontWeight','normal')
    orient landscape
    print(gcf, '-painters', '-dpsc2', '-r600', '-append', '-bestfit',...
        [figureFolder fileNameprofiles]);
    clf('reset')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
