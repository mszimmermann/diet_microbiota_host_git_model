%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot PCA of metagenomics and metatranscriptomics

% call script defining file dependenciesand global variables
addpath(genpath('.\'));
add_global_and_file_dependencies

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Required files:
% 'countsMatrix_annTable_filtered.csv'
% 'countsMatrixGetMM_DNA_table.txt'
% 'countsMatrixGetMM_RNA_table.txt'
% Output: 
% Figures:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% perform PCA on filtered metatranscriptomics 

annTablefiltered = readtable([inputFolderSeq,...
                        'countsMatrix_annTable_filtered.csv']);
countsMatrixDNAtable = readtable([inputFolderSeq,...
                        'countsMatrixGetMM_DNA_table.txt']);
countsMatrixRNAtable = readtable([inputFolderSeq,...
                        'countsMatrixGetMM_RNA_table.txt']);

countsMatrixSamples = countsMatrixDNAtable.Properties.VariableNames(...
    cellfun(@(x) contains(x, 'MSZ'), countsMatrixDNAtable.Properties.VariableNames));
countsMatrixSamplesIDs = cellfun(@(x) strrep(x, 'MSZ','M0'), countsMatrixSamples, 'unif', 0);

mouseInfo = readtable([rawdataFolder 'mouse_info.csv'], 'delim', ',');

mycolors = [0 115 178;... %dark blue
            211 96 39;... %dark orange
            204 227 240;...%light blue
            246 223 212]/256;%light orange
mymarkers = {'p','o', 's','d', 'v','>', '<', '^'};
markersize = 50;
tissue_i=2;

pca_diet = cellfun(@(x) mouseInfo.Diet(ismember(mouseInfo.Mouse_number, x)),...
    countsMatrixSamplesIDs);
        
pca_diet_unique = unique(pca_diet);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% metagenomics
pcamatrix = cell2mat(table2cell(countsMatrixDNAtable(:,countsMatrixSamples)));
pcamatrix = pcamatrix(countsMatrixDNAtable.geneFilter==1,:);

pcamatrix = log10(pcamatrix);
pcamatrix(isinf(pcamatrix))=nan;

% Keep rows complete because there are enough full rows
[coeff,~,~,~,explained] = pca(pcamatrix, 'Rows','complete');
figure;
hold on
sch = [];
for diet_i = 1:length(pca_diet_unique)
    cursamples = ismember(pca_diet, pca_diet_unique{diet_i});
    scatter(coeff(cursamples,1),coeff(cursamples,2),...
                            markersize,...
                            mymarkers{tissue_i},...
                            'MarkerFaceColor', mycolors(diet_i,:),...
                            'MarkerEdgeColor', 'k');

end
axis square
xlabel(sprintf('PC1 explained %.2f %%', explained(1)))
ylabel(sprintf('PC2 explained %.2f %%', explained(2)))
legend(pca_diet_unique)
title('All filtered genes, getMM normalization, log10, DNA')
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
            [figureFolder ...
            'fig_sup_pca_DNA_all_raw_combined_getMMnorm_log10.pdf'])
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% metatranscriptomics
pcamatrix = cell2mat(table2cell(countsMatrixRNAtable(:,countsMatrixSamples)));
pcamatrix = pcamatrix(countsMatrixRNAtable.geneFilter==1,:);

pcamatrix = log10(pcamatrix);
pcamatrix(isinf(pcamatrix))=nan;

[coeff,~,~,~,explained] = pca(pcamatrix, 'Rows','complete');
figure;
hold on
sch = [];
for diet_i = 1:length(pca_diet_unique)
    cursamples = ismember(pca_diet, pca_diet_unique{diet_i});
    scatter(coeff(cursamples,1),coeff(cursamples,2),...
                            markersize,...
                            mymarkers{tissue_i},...
                            'MarkerFaceColor', mycolors(diet_i,:),...
                            'MarkerEdgeColor', 'k');

end
axis square
xlabel(sprintf('PC1 explained %.2f %%', explained(1)))
ylabel(sprintf('PC2 explained %.2f %%', explained(2)))
legend(pca_diet_unique)
title('All filtered genes, getmm norm log10, RNA')
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
            [figureFolder ...
            'fig_sup_pca_RNA__combined_getmm_norm_log10.pdf'])
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot number of reads per species
% get chromosome id to species apping from the annotation table
species_genes_genomes = [annTablefiltered.species_genes, annTablefiltered.genome_name];
species_genes_genomes(:,1) = cellfun(@(x) x(1:strfind(x,'.')-1), species_genes_genomes(:,1), 'unif', 0);
%remove nan
species_genes_genomes(ismember(species_genes_genomes(:,2),{'nan'}),:)=[];
species_genes_genomes(:,3) = strcat(species_genes_genomes(:,1), '_', species_genes_genomes(:,2));
[~, idx] = unique(species_genes_genomes(:,3));
species_genes_genomes = species_genes_genomes(idx,:);

% % get raw reads
% countsMatrixDNAtable = readtable([inputFolderSeq,...
%                                 'countsMatrixRAW_DNA_table.txt']);
% countsMatrixRNAtable = readtable([inputFolderSeq,...
%                                 'countsMatrixRAW_RNA_table.txt']);
% 
% % get rows of sample data
% countsMatrixSamples = countsMatrixDNAtable.Properties.VariableNames(...
%     cellfun(@(x) contains(x, 'MSZ'), countsMatrixDNAtable.Properties.VariableNames));
% 
% % add species names to the table
% countsMatrixDNAtable.species_chr = cellfun(@(x) x(1:strfind(x,'.')-1), countsMatrixDNAtable.species_genes,'unif',0);
% countsMatrixRNAtable.species_chr = cellfun(@(x) x(1:strfind(x,'.')-1), countsMatrixRNAtable.species_genes,'unif',0);
% 
% % get species genomes names from the dictionary
% countsMatrixDNAtable.species_name = cellfun(@(x) x(1:strfind(x,'.')-1), countsMatrixDNAtable.species_genes,'unif',0);
% countsMatrixRNAtable.species_name = cellfun(@(x) x(1:strfind(x,'.')-1), countsMatrixRNAtable.species_genes,'unif',0);
% 
% for i=1:length(countsMatrixDNAtable.species_name)
%     curnameDNA = species_genes_genomes(ismember(species_genes_genomes(:,1),...
%         countsMatrixDNAtable.species_chr(i)), 2);
%     curnameRNA = species_genes_genomes(ismember(species_genes_genomes(:,1),...
%         countsMatrixRNAtable.species_chr(i)), 2);
%     if ~isempty(curnameDNA)
%         countsMatrixDNAtable.species_name(i) = curnameDNA;
%     else
%         countsMatrixDNAtable.species_name(i) = {'nan'};
%     end
%     if ~isempty(curnameRNA)
%         countsMatrixRNAtable.species_name(i) = curnameRNA;
%     else
%         countsMatrixRNAtable.species_name(i) = {'nan'};
%     end
%         
% end
% % fill nans with species names if previous and next genomes are the same
% nanpositions = find(ismember(countsMatrixDNAtable.species_name, {'nan'}));
% % get positions of previous and next non-nan
% nanpositions = [nanpositions zeros(length(nanpositions),2)];
% for i=1:length(nanpositions)
%     nanpositions(i,2) = nanpositions(i,1)-1;
%     while sum(nanpositions(:,1)==nanpositions(i,2))>0
%         nanpositions(i,2) = nanpositions(i,2)-1;
%     end
%     nanpositions(i,3) = nanpositions(i,1)+1;
%     % add more if next one is also nan
%     while sum(nanpositions(:,1)==nanpositions(i,3))>0
%         nanpositions(i,3) = nanpositions(i,3)+1;
%     end
% end
% % fill nan positions
% for i=1:length(nanpositions)
%     if nanpositions(i,3)<=height(countsMatrixDNAtable)
%         if isequal(countsMatrixDNAtable.species_name{nanpositions(i,2)},...
%                 countsMatrixDNAtable.species_name{nanpositions(i,3)})
%             countsMatrixDNAtable.species_name{nanpositions(i,1)} = ...
%                 countsMatrixDNAtable.species_name{nanpositions(i,2)};
%         else
%             if isequal(countsMatrixDNAtable.species_chr{nanpositions(i,1)}(1:5),...
%                     countsMatrixDNAtable.species_chr{nanpositions(i,2)}(1:5))
%                 countsMatrixDNAtable.species_name{nanpositions(i,1)} = ...
%                 countsMatrixDNAtable.species_name{nanpositions(i,2)};
%             end
%         end
%         
%     else
%         countsMatrixDNAtable.species_name{nanpositions(i,1)} = ...
%             countsMatrixDNAtable.species_name{nanpositions(i,2)};
%     end
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % fill RNA species
% % fill nans with species names if previous and next genomes are the same
% nanpositions = find(ismember(countsMatrixRNAtable.species_name, {'nan'}));
% % get positions of previous and next non-nan
% nanpositions = [nanpositions zeros(length(nanpositions),2)];
% for i=1:length(nanpositions)
%     nanpositions(i,2) = nanpositions(i,1)-1;
%     while sum(nanpositions(:,1)==nanpositions(i,2))>0
%         nanpositions(i,2) = nanpositions(i,2)-1;
%     end
%     nanpositions(i,3) = nanpositions(i,1)+1;
%     % add more if next one is also nan
%     while sum(nanpositions(:,1)==nanpositions(i,3))>0
%         nanpositions(i,3) = nanpositions(i,3)+1;
%     end
% end
% % fill nan positions
% for i=1:length(nanpositions)
%     if nanpositions(i,3)<=height(countsMatrixRNAtable)
%         if isequal(countsMatrixRNAtable.species_name{nanpositions(i,2)},...
%                 countsMatrixRNAtable.species_name{nanpositions(i,3)})
%             countsMatrixRNAtable.species_name{nanpositions(i,1)} = ...
%                 countsMatrixRNAtable.species_name{nanpositions(i,2)};
%         else
%             if isequal(countsMatrixRNAtable.species_chr{nanpositions(i,1)}(1:5),...
%                     countsMatrixRNAtable.species_chr{nanpositions(i,2)}(1:5))
%                 countsMatrixRNAtable.species_name{nanpositions(i,1)} = ...
%                 countsMatrixRNAtable.species_name{nanpositions(i,2)};
%             end
%         end
%         
%     else
%         countsMatrixRNAtable.species_name{nanpositions(i,1)} = ...
%             countsMatrixRNAtable.species_name{nanpositions(i,2)};
%     end
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % save count tables with species names to file
% writetable(countsMatrixDNAtable, ...
%     [outputFolder, 'countsMatrixRAW_DNA_table_with_species.txt']);
% writetable(countsMatrixRNAtable, ...
%     [outputFolder, 'countsMatrixRAW_RNA_table_with_species.txt']);

% get raw reads with species info
countsMatrixDNAtable = readtable([inputFolderSeq,...
                       'countsMatrixRAW_DNA_table_with_species.txt']);
countsMatrixRNAtable = readtable([inputFolderSeq,...
                       'countsMatrixRAW_RNA_table_with_species.txt']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
species_unique = unique(countsMatrixDNAtable.species_name);

reads_per_speciesDNA = zeros(length(species_unique),length(countsMatrixSamples));
reads_per_speciesRNA = zeros(length(species_unique),length(countsMatrixSamples));

for i=1:length(species_unique)
    % DNA
    curspeciesids = ismember(countsMatrixDNAtable.species_name, species_unique{i});
    
    pcamatrix = cell2mat(table2cell(countsMatrixDNAtable(curspeciesids,countsMatrixSamples)));

    reads_per_speciesDNA(i,:) = sum(pcamatrix);
    
    % RNA
    curspeciesids = ismember(countsMatrixRNAtable.species_name, species_unique{i});
    
    pcamatrix = cell2mat(table2cell(countsMatrixRNAtable(curspeciesids,countsMatrixSamples)));

    reads_per_speciesRNA(i,:) = sum(pcamatrix);
end

figure
barh(([sum(reads_per_speciesDNA)',...
    sum(reads_per_speciesRNA)']))
axis square
legend({'DNA', 'RNA'})
set(gca, 'YTick', 1:10)
set(gca, 'YTickLabel', countsMatrixSamplesIDs)
ylim([0 11])
xlabel('Number of total reads')
title('metaG and metaT reads per sample')
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
            [figureFolder 'fig_sup_barh_metaG_metaT_reads_per_sample.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot reads per species
reads_per_speciesDNA = reads_per_speciesDNA';

fig = figure('units','normalized','outerposition',[0 0 1 1]);
boxplot(reads_per_speciesDNA); 
hold on;
x=repmat(1:size(reads_per_speciesDNA,2),...
        size(reads_per_speciesDNA,1),1);
scatter(x(:),reads_per_speciesDNA(:),'filled',...
    'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);

set(gca, 'XTick', 1:length(species_unique))
set(gca, 'XTickLabel', species_unique)
xticklabel_rotate([], 90)

title('Number of total DNA reads')
% print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
%             [figureFolder 'boxplot_metaG_reads_per_species.pdf'])

% plot reads per species metaT
reads_per_speciesRNA = reads_per_speciesRNA';

fig = figure('units','normalized','outerposition',[0 0 1 1]);
boxplot(reads_per_speciesRNA); 
hold on;
x=repmat(1:size(reads_per_speciesRNA,2),...
        size(reads_per_speciesRNA,1),1);
scatter(x(:),reads_per_speciesRNA(:),'filled',...
    'MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.15);

set(gca, 'XTick', 1:length(species_unique))
set(gca, 'XTickLabel', species_unique)
xticklabel_rotate([], 90)

title('Number of total RNA reads')
% print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
%             [figureFolder 'boxplot_metaT_reads_per_species.pdf'])
