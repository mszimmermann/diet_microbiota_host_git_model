%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze species composition by metaphlan2

% call script defining file dependenciesand global variables
addpath(genpath('.\'));
add_global_and_file_dependencies

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Required files:
% 'merged_humann2_metaphlan_species_filtered_mouseDNA_OTUs.txt'
% 'mouse_info.csv'
% Output:
% Files:
% 'merged_humann2_metaphlan_species_filtered_mouseDNA_OTUs.txt'
% 'merged_humann2_metaphlan_species_filtered_mouseDNA_OTUs_with_FC.txt'
% Figures:
% 'fig_1c_barplot_OTU_per_mouse.pdf'
% 'fig_1d_clustergram_DC_DNA_humann2_metaphlan_selected_filterd.pdf'
% 'barplot_DNA_otu_PHYLUM_CV_ranksum_pvalue_FDR.pdf'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%read_otu_table;
filename = [inputFolderSeq,...
    'merged_humann2_metaphlan_bugs_list_mouseDNA_OTUs.txt'];

mergeddietmetagenomeabundancetable = readtable(filename);
% remove first row containing metaphlan headers
mergeddietmetagenomeabundancetable(1,:)=[];
tablecols = mergeddietmetagenomeabundancetable.Properties.VariableNames;
dataidx = cellfun(@(x) contains(x, 'MSZ'), tablecols);
mousecols = tablecols(dataidx);
mousecols = cellfun(@(x) x(1:strfind(x, '_')-1), mousecols, 'unif', 0);
mousecols = cellfun(@(x) strrep(x, 'MSZ', 'M0'), mousecols, 'unif', 0);

mouseInfo = readtable([rawdataFolder 'mouse_info.csv'], 'delim', ',');

mousecolsdiet = mousecols;
for i=1:length(mousecols)
    mousecolsdiet(i) = mouseInfo.Diet(ismember(mouseInfo.Mouse_number, mousecols{i}));
end

if isnumeric(mergeddietmetagenomeabundancetable{:,dataidx})
    clustmat = table2array(mergeddietmetagenomeabundancetable(:,dataidx));
else

    clustmat = cellfun(@(x) str2double(x), table2cell(mergeddietmetagenomeabundancetable(:,dataidx)));
end
clustCols = mousecols;
clustRows = cellstr(mergeddietmetagenomeabundancetable.ID);

contains_used = ones(size(clustRows));

selectRows = sum(clustmat,2) ~=0 & ...
             cellfun(@(x) contains(x, '|s'), clustRows) & ...
             cellfun(@(x) ~contains(x, '|t'), clustRows) &...
             cellfun(@(x) ~contains(x, 'virus'), clustRows) &...
             cellfun(@(x) ~contains(x, 'unclassified'), clustRows);


clustRowNames = clustRows(selectRows);
clustRowNames = cellfun(@(x) x(strfind(x, '|s')+4:end), clustRowNames,'unif',0);
clustRowAbbr = cell(size(clustRowNames));
for i=1:length(clustRowNames)
    curabbr = strsplit(clustRowNames{i}, '_');
    clustRowAbbr{i} = [curabbr{1}(1), curabbr{2}(1:3)];
end
clustdist = 'correlation';%'cityblock';%'euclidean';%'spearman';%'hamming';%

clustmat_mean = [mean(clustmat(:,1:5),2), mean(clustmat(:,6:10),2)];

% filter based on values (>=5 samples)
clustmat = clustmat(selectRows,:);

removeRows = sum(clustmat==0, 2)>5;
clustmat = clustmat(~removeRows,:);
clustRowNames = clustRowNames(~removeRows);
clustRowAbbr = clustRowAbbr(~removeRows);

selectRows = cellfun(@(x) ~contains(x, 'lactis'), clustRowNames) &...
             cellfun(@(x) ~contains(x, 'epidermidis'), clustRowNames);

fprintf('Removed rows: \n%s\n%s\n', clustRowNames{~selectRows})
fprintf('Removed rows max: \n%.4f\n%.4f\n', max(clustmat(~selectRows,:),[],2))
fprintf('Removed rows mean: \n%.4f\n%.4f\n', mean(clustmat(~selectRows,:),2))
fprintf('Removed rows std: \n%.4f\n%.4f\n', std(clustmat(~selectRows,:),[],2))

clustmat = clustmat(selectRows,:);
clustRowNames = clustRowNames(selectRows);
clustRowAbbr = clustRowAbbr(selectRows);

% normalize with z-score
clustergrammat = zscore(clustmat,0,2);

set(groot, 'DefaultTextInterpreter', 'none')
set(groot, 'DefaultAxesTickLabelInterpreter', 'none')

cgo=clustergram(clustergrammat,...'DisplayRange', truncVal,...
             'RowLabels', clustRowNames,...
             'ColumnLabels', clustCols,...
             'Symmetric', 1,...
             'ColumnPdist', clustdist, 'RowPdist', clustdist,...
             'Colormap', slanCM('vik'));
% using new divergent colormap slanCM('vik') from 
% Zhaoxu Liu / slandarer (2024). 200 colormap (https://www.mathworks.com/matlabcentral/fileexchange/120088-200-colormap), MATLAB Central File Exchange. Accessed 21. August 2024. 

% solution to turn on colormap programmatically from https://stackoverflow.com/questions/20648627/turn-on-colorbar-programmatically-in-clustergram
cbButton = findall(gcf,'tag','HMInsertColorbar');
ccb = get(cbButton,'ClickedCallback');
set(cbButton,'State','on')
ccb{1}(cbButton,[],ccb{2})
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot figure from clustergram
fig = cgo.plot;
C = findall(gcf,'type','ColorBar');                         
C.Label.String = 'Relative species abundance';
    
set(0,'ShowHiddenHandles','on')
% Get all handles from root
allhnds = get(0,'Children');
% Find hearmap axis and change the font size
h = findall(allhnds, 'Tag', 'HeatMapAxes');
set(h, 'FontSize', 6)

% print to figure
orient landscape

print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
    [figureFolder, ...
    'fig_1d_clustergram_DC_DNA_humann2_metaphlan_selected_filterd_vikcmap_distanceCorr.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
otuMAT = clustmat;
otuCols = clustCols;
otuRows = clustRowAbbr;
otuRowNames = clustRowNames;
otuTable = array2table(clustmat, 'RowNames', otuRowNames, 'VariableNames', otuCols);
% write filtered OTU table to file
writetable(otuTable, ...
    [outputFolder,...
    'merged_humann2_metaphlan_species_filtered_mouseDNA_OTUs.txt'],...
    'WriteRowNames', 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot OTUs as bar plot 
sumRows = sum(clustmat~=0,2);
sumRows(~selectRows)=0;
sumRows_sort = sort(sumRows, 'descend');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
ba = bar((fliplr(clustmat))','stacked', 'FaceColor','flat');
mycolors = [70 170 150;...
            0.8*[70 170 150];...
            222 45 38;...
            0.9*[222 45 38];...
            0.8*[222 45 38];...
            0.7*[222 45 38];...
            0.6*[222 45 38];...
            0.5*[222 45 38];...
            8 81 156;...
            0.9*[8 81 156];...
            0.8*[8 81 156];...
            0.7*[8 81 156];...
            0.6*[8 81 156];...
            205 145 60]/255;
for i=1:size(mycolors,1)
    ba(i).CData = mycolors(i,:);
end
legend(clustRowAbbr, 'Location', 'EastOutside')
set(gca, 'XTick', 1:10)
set(gca, 'XTickLabel', fliplr(clustCols))
xlim([0.5 10.5])
orient landscape
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
           [figureFolder, ...
           'fig_1c_barplot_OTU_per_mouse.pdf'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate p-values of the OTUs across groups
% select only species names from rows
meanTable = [mean(otuMAT(:,ismember(mousecolsdiet, {'HFD'})),2) ...
             mean(otuMAT(:,ismember(mousecolsdiet, {'CTRL'})),2)];
stdTable = [std(otuMAT(:,ismember(mousecolsdiet, {'HFD'})),[],2) ...
            std(otuMAT(:,ismember(mousecolsdiet, {'CTRL'})),[],2)];
pTable = zeros(size(meanTable,1),1);
meanFC = zeros(size(meanTable,1),1);
stdFC = zeros(size(meanTable,1),1);
for i=1:size(meanTable,1)
    pTable(i) = ranksum(otuMAT(i,ismember(mousecolsdiet, {'HFD'})),...
                        otuMAT(i,ismember(mousecolsdiet, {'CTRL'})));
    meanFC(i) = mean(otuMAT(i,ismember(mousecolsdiet, {'HFD'})))/...
                mean(otuMAT(i,ismember(mousecolsdiet, {'CTRL'})));
    stdFC(i) = meanFC(i)*sqrt((std(otuMAT(i,ismember(mousecolsdiet, {'HFD'})))/...
                                mean(otuMAT(i,ismember(mousecolsdiet, {'HFD'}))))^2 +...
                              (std(otuMAT(i,ismember(mousecolsdiet, {'CTRL'})))/...
                               mean(otuMAT(i,ismember(mousecolsdiet, {'CTRL'}))))^2);
    
    stdFC(i) = stdFC(i)/(meanFC(i)*log(2));
    meanFC(i) = log2(meanFC(i));
end
fdrTable = mafdr(pTable, 'bhfdr', 1);
% bar plot species abundances
figure;
bar(meanTable(:, [2 1]))
hold on
errorbar([0.85:1:length(meanTable(:,2))], meanTable(:,2), stdTable(:,2), 'k.')
errorbar([1.15:1:length(meanTable(:,1))+1], meanTable(:,1), stdTable(:,1), 'k.')

for i=1:length(fdrTable)
    text(i-0.25, max(meanTable(i,:))+1, sprintf('%.2f',fdrTable(i)))
end
set(gca, 'XTick', 1:length(clustRowAbbr))
set(gca, 'XTickLabel', clustRowAbbr)

legend({'CTR', 'HFD'})
ylim([0 35])
ylabel('Relative species abundance, %')

print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
    [figureFolder, ...
    'barplot_DNA_otu_PHYLUM_CV_ranksum_pvalue_FDR.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add differential abundance analysis as well
otuTable.mean_HFD = meanTable(:,1);
otuTable.std_HFD = stdTable(:,1);
otuTable.mean_HCD = meanTable(:,2);
otuTable.std_HCD = stdTable(:,2);
otuTable.FC_HFD_CTR = meanFC;
otuTable.stdFC_HFD_CTR = stdFC;
otuTable.pval_ranksum_HFD_CTR = pTable;
otuTable.FDR_HFD_CTR = fdrTable;

% write filtered OTU table to file
writetable(otuTable, [outputFolder, ...
    'merged_humann2_metaphlan_species_filtered_mouseDNA_OTUs_with_FC.txt'],...
    'WriteRowNames', 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bar plot species abundances
% figure;
% barh(meanFC)
% hold on
% errorbar(meanFC, [1:1:length(stdFC)], stdFC, 'k.', 'horizontal')
% for i=1:length(fdrTable)
%     text(max(meanFC(i))+1, i-0.25, sprintf('%.2f',fdrTable(i)))
% end
% set(gca, 'YTick', 1:length(clustRowAbbr))
% set(gca, 'YTickLabel', clustRowAbbr)
% xlim([-3 3])
% xlabel('Fold change to CTR, log2')
% orient landscape
% 
% print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
%     [figureFolder, ...
%     'fig_sup_barplot_fc_DNA_OTUPHYLUM_CV_HFD_to_CTRL_log2_ranksum_pvalueFDR.pdf'])

