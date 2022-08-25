%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot metaG vs metaT without ribosomal RNA

% call script defining file dependenciesand global variables
addpath(genpath('.\'));
add_global_and_file_dependencies

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Required files:
% 'geneAnnTable_full.csv'
% 'countsMatrixRAW_DNA_table_with_species.txt'
% 'countsMatrixRAW_RNA_table_with_species.txt'
% 'merged_humann2_metaphlan_species_filtered_mouseDNA_OTUs.txt'
% Output: 
% Figures:
% 'fig_sup_scatter_relab_metaphlan_comparison_all_DNA_RNA.pdf'
% 'fig_sup_scatter_relab_metaphlan_comparison_all_DNA_RNA_noribo.pdf'
% 'fig_1b_histfit_metaG_metaT_genome_coverage_rnaRESEQ.pdf'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read gene annotation file
annTable = readtable([inputFolderSeq, 'geneAnnTable_full.csv']);

countsMatrixDNAtable = readtable([inputFolderSeq,...
                'countsMatrixRAW_DNA_table_with_species.txt']);
countsMatrixRNAtable = readtable([inputFolderSeq,...
                'countsMatrixRAW_RNA_table_with_species.txt']);

filenameOTU = [inputFolderSeq,...
        'merged_humann2_metaphlan_species_filtered_mouseDNA_OTUs.txt'];
mergeddietmetagenomeabundancetable = readtable(filenameOTU);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare raw reads to the whole transcriptome mapping
% get raw counts of species per method

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate back abundance per species
% get species names
species_names_unique = unique(countsMatrixDNAtable.species_name);
% format species names with strain name and underscore
species_names_DNARNA = species_names_unique;
for i=1:length(species_names_unique)
    curspecies = strsplit(species_names_unique{i},' ');
    species_names_DNARNA{i} = strjoin(curspecies(1:2),'_');
end

% get data columns
mouse_columns = cellfun(@(x) contains(x, 'MSZ'), countsMatrixDNAtable.Properties.VariableNames);
mouse_columns_names = countsMatrixDNAtable.Properties.VariableNames(mouse_columns);
% reformat to the same sample names in M0XX format
mouse_columns_names = cellfun(@(x) strrep(x, 'MSZ', 'M0'), mouse_columns_names, 'unif', 0);

speciesAbundance_DNA = zeros(length(species_names_unique),nnz(mouse_columns));
% calculate total counts
totalCounts = sum(table2array(countsMatrixDNAtable(:, mouse_columns)));

for species_i = 1:length(species_names_unique)
    curmatrix = table2array(countsMatrixDNAtable(ismember(countsMatrixDNAtable.species_name, species_names_unique{species_i}),...
        mouse_columns));
    speciesAbundance_DNA(species_i,:) = sum(curmatrix)./totalCounts;
end

speciesAbundance_RNA = zeros(length(species_names_unique),nnz(mouse_columns));
% calculate total counts
totalCounts = sum(table2array(countsMatrixRNAtable(:, mouse_columns)));

for species_i = 1:length(species_names_unique)
    curmatrix = table2array(countsMatrixRNAtable(ismember(countsMatrixRNAtable.species_name, species_names_unique{species_i}),...
        mouse_columns));
    speciesAbundance_RNA(species_i,:) = sum(curmatrix)./totalCounts;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find intesection between metaG and metaT and metaphlan OTU table
metaphlanOTU = table2array(mergeddietmetagenomeabundancetable(:,2:end));
metaphlanOTUnames = mergeddietmetagenomeabundancetable.Row;
metaphlanOTUsamples = mergeddietmetagenomeabundancetable.Properties.VariableNames(2:end);

% sort all tables like metaphlan OTU (by phylum)
[mergedNames, row1, row2] = intersect(metaphlanOTUnames, species_names_DNARNA, 'stable');
metaphlanOTU = metaphlanOTU(row1,:);
speciesAbundance_DNA = speciesAbundance_DNA(row2,:);
speciesAbundance_RNA = speciesAbundance_RNA(row2,:);
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot metaphlan vs metaG vs metaT vs OTU
fig = figure('units','normalized','outerposition',[0 0 1 1]);
for idx = 1:3
    switch idx
        case 1
            method1 = 'Metaphlan OTU';
            method2 = 'DNA';
            method1table = metaphlanOTU;
            method1samples = metaphlanOTUsamples;
            method2table = speciesAbundance_DNA;
            method2samples = mouse_columns_names;
        case 2
            method1 = 'Metaphlan OTU';
            method2 = 'RNA';
            method1table = metaphlanOTU;
            method1samples = metaphlanOTUsamples;
            method2table = speciesAbundance_RNA;
            method2samples = mouse_columns_names;
        case 3
            method1 = 'DNA';
            method2 = 'RNA';
            method1table = speciesAbundance_DNA;
            method1samples = mouse_columns_names;
            method2table = speciesAbundance_RNA;
            method2samples = mouse_columns_names;
    end
    
    [~, col1, col2] = intersect(method1samples, method2samples);
    
    % regular scale
    subplot(2,4,idx)
    hold on
    for i=1:length(col1)
          scatter((method1table(:, col1(i))),...
                  (method2table(:, col2(i))),...
                  15,...
                  species_colors, 'filled');

    end
    xlabel(method1)
    ylabel(method2)
    datacorr = ([reshape(method1table,[],1),...
                      reshape(method2table,[],1)]);
    datacorr(isinf(datacorr))=nan;
    corrPC = corr(datacorr,...
                    'rows', 'complete');
    corrSC = corr(datacorr,...
                   'rows', 'complete', 'type', 'Spearman');
    title(sprintf('PCC=%.2f SCC=%.2f',corrPC(1,2), corrSC(1,2)))
    axis square
    
    % log scale
    subplot(2,4,4+idx)
    hold on
    for i=1:length(col1)
          scatter(log10(method1table(:, col1(i))),...
                  log10(method2table(:, col2(i))),...
                  15,...
                  species_colors, 'filled');

    end
    xlabel([method1, ', log10'])
    ylabel([method2, ', log10'])
    datacorr = log10([reshape(method1table,[],1),...
                      reshape(method2table,[],1)]);
    datacorr(isinf(datacorr))=nan;
    corrPC = corr(datacorr,...
                    'rows', 'complete');
    corrSC = corr(datacorr,...
                   'rows', 'complete', 'type', 'Spearman');
    title(sprintf('PCC=%.2f SCC=%.2f',corrPC(1,2), corrSC(1,2)))
    axis square
            
end
% make legend subplot
subplot(2,4,idx+1)
hold on
for species_i=1:length(row1)
    scatter(log10(method1table(species_i, col1(i))),...
                  log10(method2table(species_i, col2(i))),...
                  15,...
                  species_colors(species_i,:), 'filled');
end
lgd = legend(cellfun(@(x) strrep(x, '_', ' '), mergedNames, 'unif', 0));
%lgd.NumColumns = 2;
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
        [figureFolder, ...
        'fig_sup_scatter_relab_metaphlan_comparison_all_DNA_RNA.pdf'])
%close(fig)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate metaG vs metaT without ribosomal RNA

speciesAbundance_DNA_noribo = zeros(length(species_names_unique),nnz(mouse_columns));
% calculate total counts (no ribosomal)
curgenes = annTable.product;
ribosomalRNAflag = cellfun(@(x) contains(x, 'ribosomal RNA'), curgenes);
totalCounts_noribo = sum(table2array(countsMatrixDNAtable(~ribosomalRNAflag, mouse_columns)));

for species_i = 1:length(species_names_unique)
    curmatrix = table2array(countsMatrixDNAtable(ismember(countsMatrixDNAtable.species_name, species_names_unique{species_i}),...
        mouse_columns));
    curgenes = annTable(ismember(countsMatrixDNAtable.species_name, species_names_unique{species_i}),:);
    curgenes = curgenes.product;
    ribosomalRNAflag = cellfun(@(x) contains(x, 'ribosomal RNA'), curgenes);
    % remove ribosomal RNA
    speciesAbundance_DNA_noribo(species_i,:) = sum(curmatrix(~ribosomalRNAflag,:))./totalCounts_noribo;
end

speciesAbundance_RNA_noribo = zeros(length(species_names_unique),nnz(mouse_columns));
% calculate total counts
curgenes = annTable.product;
ribosomalRNAflag = cellfun(@(x) contains(x, 'ribosomal RNA'), curgenes);
totalCounts_noribo = sum(table2array(countsMatrixRNAtable(~ribosomalRNAflag, mouse_columns)));

for species_i = 1:length(species_names_unique)
    curmatrix = table2array(countsMatrixRNAtable(ismember(countsMatrixRNAtable.species_name, species_names_unique{species_i}),...
        mouse_columns));
    curgenes = annTable(ismember(countsMatrixRNAtable.species_name, species_names_unique{species_i}),:);
    curgenes = curgenes.product;
    ribosomalRNAflag = cellfun(@(x) contains(x, 'ribosomal RNA'), curgenes);
    % remove ribosomal RNA
    speciesAbundance_RNA_noribo(species_i,:) = sum(curmatrix(~ribosomalRNAflag,:))./totalCounts_noribo;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot scatter of species abundances without ribosomal RNA

[mergedNames, row1, row2] = intersect(metaphlanOTUnames, species_names_DNARNA, 'stable');
metaphlanOTU = metaphlanOTU(row1,:);
speciesAbundance_DNA_noribo = speciesAbundance_DNA_noribo(row2,:);
speciesAbundance_RNA_noribo = speciesAbundance_RNA_noribo(row2,:);
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot metaphlan vs metaG vs metaT vs OTU
fig = figure('units','normalized','outerposition',[0 0 1 1]);
for idx = 1:3
    switch idx
        case 1
            method1 = 'Metaphlan OTU';
            method2 = 'DNA';
            method1table = metaphlanOTU;
            method1samples = metaphlanOTUsamples;
            method2table = speciesAbundance_DNA_noribo;
            method2samples = mouse_columns_names;
        case 2
            method1 = 'Metaphlan OTU';
            method2 = 'RNA';
            method1table = metaphlanOTU;
            method1samples = metaphlanOTUsamples;
            method2table = speciesAbundance_RNA_noribo;
            method2samples = mouse_columns_names;
        case 3
            method1 = 'DNA';
            method2 = 'RNA';
            method1table = speciesAbundance_DNA_noribo;
            method1samples = mouse_columns_names;
            method2table = speciesAbundance_RNA_noribo;
            method2samples = mouse_columns_names;
    end
    
    [~, col1, col2] = intersect(method1samples, method2samples);
    
    % regular scale
    subplot(2,4,idx)
    hold on
    for i=1:length(col1)
          scatter((method1table(:, col1(i))),...
                  (method2table(:, col2(i))),...
                  15,...
                  species_colors, 'filled');

    end
    xlabel(method1)
    ylabel(method2)
    datacorr = ([reshape(method1table,[],1),...
                      reshape(method2table,[],1)]);
    datacorr(isinf(datacorr))=nan;
    corrPC = corr(datacorr,...
                    'rows', 'complete');
    corrSC = corr(datacorr,...
                   'rows', 'complete', 'type', 'Spearman');
    title(sprintf('PCC=%.2f SCC=%.2f',corrPC(1,2), corrSC(1,2)))
    axis square
    
    % log scale
    subplot(2,4,4+idx)
    hold on
    for i=1:length(col1)
          scatter(log10(method1table(:, col1(i))),...
                  log10(method2table(:, col2(i))),...
                  15,...
                  species_colors, 'filled');

    end
    xlabel([method1, ', log10'])
    ylabel([method2, ', log10'])
    datacorr = log10([reshape(method1table,[],1),...
                      reshape(method2table,[],1)]);
    datacorr(isinf(datacorr))=nan;
    corrPC = corr(datacorr,...
                    'rows', 'complete');
    corrSC = corr(datacorr,...
                   'rows', 'complete', 'type', 'Spearman');
    title(sprintf('PCC=%.2f SCC=%.2f',corrPC(1,2), corrSC(1,2)))
    axis square
            
end
% make legend subplot
subplot(2,4,idx+1)
hold on
for species_i=1:length(row1)
    scatter(log10(method1table(species_i, col1(i))),...
                  log10(method2table(species_i, col2(i))),...
                  15,...
                  species_colors(species_i,:), 'filled');
end
lgd = legend(cellfun(@(x) strrep(x, '_', ' '), mergedNames, 'unif', 0));
%lgd.NumColumns = 2;
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
        [figureFolder, ...
        'fig_sup_scatter_relab_metaphlan_comparison_all_DNA_RNA_noribo.pdf'])
close(fig)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENOME COVERAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate genome coverage per species
% metaG coverage
speciesmetaGcoverage = zeros(length(species_names_unique),10);
speciesmetaGcoverage_fraction = zeros(length(species_names_unique),10);
for species_i = 1:length(species_names_unique)
    curmatrix = table2array(countsMatrixDNAtable(ismember(countsMatrixDNAtable.species_name, species_names_unique{species_i}),...
        mouse_columns));
    
    speciesmetaGcoverage(species_i,:) = ...
        sum(curmatrix>0);
    speciesmetaGcoverage_fraction(species_i,:) = ...
        sum(curmatrix>0)/...
        size(curmatrix,1);
end

% metaT coverage
speciesmetaTcoverage = zeros(length(species_names_unique),10);
speciesmetaTcoverage_fraction = zeros(length(species_names_unique),10);
for species_i = 1:length(species_names_unique)
    curmatrix = table2array(countsMatrixRNAtable(ismember(countsMatrixRNAtable.species_name, species_names_unique{species_i}),...
        mouse_columns));
    
    speciesmetaTcoverage(species_i,:) = ...
        sum(curmatrix>0);
    speciesmetaTcoverage_fraction(species_i,:) = ...
        sum(curmatrix>0)/...
        size(curmatrix,1);
end

% plot coverage per sample
mycolors = [0 115 178;... %dark blue
            211 96 39;... %dark orange
            204 227 240;...%light blue
            246 223 212]/256;%light orange
[~, idxsort] = sort(mean(speciesmetaTcoverage_fraction,2));
sfig = figure('units','normalized','innerposition',[0.1 0.1 0.8 0.5]);
subplot(1,2,1)
plotdata = speciesmetaTcoverage(idxsort,:);
scatter(plotdata(:),...
        repmat([1:size(speciesmetaTcoverage_fraction,1)]',...
        size(speciesmetaTcoverage_fraction,2),1),...
        50,...
        repmat([repmat(mycolors(2,:),5,1);...
                repmat(mycolors(1,:),5,1)],...
                size(speciesmetaTcoverage_fraction,1),1),'filled' )
axis square
xlabel('Number of genes')
set(gca, 'YTick', 1:length(idxsort))
set(gca, 'YTickLabel', species_names_unique(idxsort))
ylim([0 length(idxsort)+1])
title('Number of genes detected')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,2)
plotdata = speciesmetaTcoverage_fraction(idxsort,:);
scatter(plotdata(:),...
        5*repmat([1:size(speciesmetaTcoverage_fraction,1)]',...
        size(speciesmetaTcoverage_fraction,2),1),...
        50,...
        repmat([repmat(mycolors(2,:),5,1);...
                repmat(mycolors(1,:),5,1)],...
                size(speciesmetaTcoverage_fraction,1),1),'filled' )
axis square
xlabel('Number of genes')
set(gca, 'YTick', 5*(1:length(idxsort)))
set(gca, 'YTickLabel', species_names_unique(idxsort))
ylim([0 5*(length(idxsort)+1)])
title('Fraction of genes detected')
suptitle('Genome coverage by metaG and metaT')
orient landscape
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
            [figureFolder, ...
            'fig_sup_scatter_metaG_metaT_genome_coverage_rnaRESEQ_persample.pdf'])
%close(sfig);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plt distributions
for method_type = 1:2
    if method_type == 1
        plotdata = speciesmetaGcoverage_fraction(idxsort,:);
    else
        plotdata = speciesmetaTcoverage_fraction(idxsort,:);
    end
    figure
    plotcurve = cell(size(plotdata,1),1);
    for i=1:size(plotdata,1)
        curfit = plotdata(i,:);
        curfit(curfit==max(curfit))=[];
        %curfit(curfit==min(curfit))=[];
        h = histfit(curfit,10, 'kernel');
        curline = h(2);
        plotcurve{i} = [curline.XData; curline.YData];
    end
    if method_type==1
        plotcurveDNA = plotcurve;
    else
        plotcurveRNA = plotcurve;
    end
end

fig = figure('units','normalized','innerposition',[0.1 0.1 0.8 0.5]);
hold on
for i=1:size(plotdata,1)
    curcurve = plotcurveDNA{i};
    fill(curcurve(1,:),...
         i*5 + 4*[curcurve(2,:)/max(curcurve(2,:))],...
         'k')
    curcurve = plotcurveRNA{i};
    fill(curcurve(1,:),...
         i*5 + 4*[curcurve(2,:)/max(curcurve(2,:))],...
         [.8 .8 .8]) 
end
xlim([0 1])
ylim([0 5*length(idxsort)+5])
axis square
xlabel('Fraction of genes')
set(gca, 'YTick', 5:5:5*length(idxsort))
set(gca, 'YTickLabel', species_names_unique(idxsort))
legend({'DNA', 'RNA'},'location','northwest')
suptitle('Genome coverage by metaG and metaT')
orient landscape
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
            [figureFolder,...
            'fig_1b_histfit_metaG_metaT_genome_coverage_rnaRESEQ.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate coverage and intensity per pathway or gene type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%pwName = 'ribosomal protein';
%pwName = 'ribosomal RNA';
% pwName = 'ribosomal RNA';
% speciesmetaTcoveragePW_num = zeros(length(species_names_unique),10);
% speciesmetaTcoveragePW_intensity = zeros(length(species_names_unique),10);
% speciesmetaTcoveragePW_intensityratio = zeros(length(species_names_unique),10);
% for species_i = 1:length(species_names_unique)
%     curmatrix = table2array(countsMatrixDNAtable(ismember(countsMatrixDNAtable.species_name, species_names_unique{species_i}),...
%         mouse_columns));
%     curgenes = annTable(ismember(countsMatrixDNAtable.species_name, species_names_unique{species_i}),:);
% %     pwidx = find(cellfun(@(x) ~contains(x,pwName),...
% %                              curgenes.product));
%    
%     pwidx = find(cellfun(@(x) contains(x,pwName),...
%                              curgenes.product));
%                          
%     speciesmetaTcoveragePW_num(species_i,:) = ...
%         nansum(curmatrix(pwidx,:)>0);
%     speciesmetaTcoveragePW_intensity(species_i,:) = ...
%         nansum(curmatrix(pwidx,:));
%     speciesmetaTcoveragePW_intensityratio(species_i,:) = ...
%         nansum(curmatrix(pwidx,:))./...
%         nansum(curmatrix);
% end
% fig = figure('units','normalized','innerposition',[0.1 0.1 0.8 0.5]);
% subplot(1,3,1)
% plotdata = speciesmetaTcoveragePW_intensityratio(idxsort,:);
% scatter(plotdata(:),...
%         repmat([1:size(speciesmetaTcoverage_fraction,1)]',...
%         size(speciesmetaTcoverage_fraction,2),1),...
%         50,...
%         repmat([repmat(mycolors(2,:),5,1);...
%                 repmat(mycolors(1,:),5,1)],...
%                 size(speciesmetaTcoverage_fraction,1),1),'filled' )
% axis square
% xlabel(sprintf('Fraction of %s', pwName))
% set(gca, 'YTick', 1:length(idxsort))
% set(gca, 'YTickLabel', species_names_unique(idxsort))
% ylim([0 length(idxsort)+1])
% 
% subplot(1,3,2)
% plotdata = log10(speciesmetaTcoveragePW_intensity(idxsort,:));
% scatter(plotdata(:),...
%         repmat([1:size(speciesmetaTcoverage_fraction,1)]',...
%         size(speciesmetaTcoverage_fraction,2),1),...
%         50,...
%         repmat([repmat(mycolors(2,:),5,1);...
%                 repmat(mycolors(1,:),5,1)],...
%                 size(speciesmetaTcoverage_fraction,1),1),'filled' )
% axis square
% xlabel(sprintf('Summed intensity, log10, %s', pwName))
% set(gca, 'YTick', 1:length(idxsort))
% set(gca, 'YTickLabel', species_names_unique(idxsort))
% ylim([0 length(idxsort)+1])
% 
% figure
% plotdata =  speciesmetaTcoveragePW_intensityratio(idxsort,:);
% plotcurve = cell(size(plotdata,1),1);
% for i=1:size(plotdata,1)
%     curfit = plotdata(i,:);
%     curfit(curfit==max(curfit))=[];
%     %curfit(curfit==min(curfit))=[];
%     h = histfit(curfit,10, 'kernel');
%     curline = h(2);
%     plotcurve{i} = [curline.XData; curline.YData];
% end
% subplot(1,3,3)
% hold on
% for i=1:size(plotdata,1)
%     curcurve = plotcurve{i};
%     fill(curcurve(1,:),...
%          i*5 + 4*[curcurve(2,:)/max(curcurve(2,:))],...
%          [.8 .8 .8]) 
% end
% xlim([0 1])
% ylim([0 5*length(idxsort)+5])
% axis square
% xlabel('Fraction of intensity')
% set(gca, 'YTick', 5:5:5*length(idxsort))
% set(gca, 'YTickLabel', species_names_unique(idxsort))
% suptitle(sprintf('NOT %s coverage by metaT (CTR:blue, HFD:orange)',...
%             pwName))
% orient landscape
%print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
%            'scatter_and_dist_metaG_genome_coverage_NONribosomal_rnaRESEQ.pdf')

         
         
         

