% Fit steady state fluxes to GIT metabolites

% call script defining file dependenciesand global variables
addpath(genpath('.\'));
add_global_and_file_dependencies

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data requirements:
% 'mouse_info.csv'
% 'metabolites_allions_combined_formulas_with_metabolite_filters.csv'
% 'metabolites_allions_combined_norm_intensity.csv'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output: 
% Figures:
% 'fig2a_barh_annotated_ions_per_tissue'
% 'fig2b_pca_metabolomics_all_annotated_ions_updated_metfilter_lognorm_Impute5000'
% 'fig2b_sup_pca_metabolomics_per_tissue_annotated_ions_updated_metfilter_lognorm_Impute5000.ps'
% 'fig_sup_heatmap_fraction_of_annotated_ion_overlap_per_tissue_updated_metfilter'
% 'fig_sup_heatmap_number_of_annotated_ion_overlap_per_tissue_updated_metfilter'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mouseInfo = readtable([rawdataFolder 'mouse_info.csv'], 'delim', ',');
annotationTable = readtable([inputFolder ...
        'metabolites_allions_combined_formulas_with_metabolite_filters.csv']);
combinedIntensitiesData = readtable([inputFolder ...
        'metabolites_allions_combined_norm_intensity.csv']);


% get the filter vector
binaryFilter = annotationTable.MetaboliteFilter;

% get intensities from the table as a separae matrix
sample_columns = cellfun(@(x) contains(x,'_M'), combinedIntensitiesData.Properties.VariableNames);
combinedIntensitiesNorm = cell2mat(table2cell(combinedIntensitiesData(:,sample_columns)));
% get information on sample type, tissue and condition
sampleNames = combinedIntensitiesData.Properties.VariableNames(sample_columns);
sampleNames_parsed = cellfun(@(x) strsplit(x, '_'), sampleNames, 'unif', 0);
combinedTissues = cellfun(@(x) x{3}, sampleNames_parsed, 'unif', 0);
combinedType = cellfun(@(x) x{2}, sampleNames_parsed, 'unif', 0);
combinedDiet = cellfun(@(x) x{1}, sampleNames_parsed, 'unif', 0);

% get info on mouse gender
mouseGender = mouseInfo.Gender;
combinedGender = repmat(mouseGender, round(length(combinedTissues)/length(mouseGender)),1);

% calculate number of ions per sample
combined_ions_per_sample = sum(combinedIntensitiesNorm>5000,1);

% calculate number of ions per condition 
pca_type = combinedType;
pca_tissue = combinedTissues;
pca_diet = combinedDiet;
% calculate unique values
pca_diet_unique = unique(pca_diet);
pca_type_unique = unique(pca_type);
pca_tissue_unique = unique(pca_tissue);
% sort according to GI propagation
pca_tissue_unique = pca_tissue_unique([5 6 7 1 2 3 4 8]);
% prepare data
combined_cond_data = zeros(size(combinedIntensitiesNorm,1),...
    length(pca_diet_unique)*length(pca_type_unique)*length(pca_tissue_unique));
combined_cond_diet = cell(length(pca_diet_unique)*length(pca_type_unique)*length(pca_tissue_unique),1);
combined_cond_type = cell(length(pca_diet_unique)*length(pca_type_unique)*length(pca_tissue_unique),1);
combined_cond_tissue = cell(length(pca_diet_unique)*length(pca_type_unique)*length(pca_tissue_unique),1);

% calculate mean abundance of each ion per tissue and condition
% only calculate mean of >=3 samples
idx = 1;
for diet_i = 1:length(pca_diet_unique)
    for type_i = 1:length(pca_type_unique)
        for tissue_i = 1:length(pca_tissue_unique)
            cursamples = ismember(pca_diet, pca_diet_unique{diet_i}) &...
                         ismember(pca_type, pca_type_unique{type_i}) &...
                         ismember(pca_tissue, pca_tissue_unique{tissue_i});

            curdata = combinedIntensitiesNorm(:, cursamples);
            % replace 5000 to nan
            curdata(curdata<=5000) = nan;
            combined_cond_data(:, idx) = nanmean(curdata,2);
            % replace to nan if in less than 3 samples
            replace_to_nan = sum(isnan(curdata),2)>2;
            combined_cond_data(replace_to_nan, idx) = nan;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            combined_cond_diet{idx} = pca_diet_unique{diet_i};
            combined_cond_type{idx} = pca_type_unique{type_i};
            combined_cond_tissue{idx} = pca_tissue_unique{tissue_i};
            idx = idx+1;
        end
    end
end

% calculate number of ions per sample
combined_ions_per_sample = sum(~isnan(combined_cond_data),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate shared ions between conditions
combined_cond_data_all = combined_cond_data;
% leave only annotated ions
combined_cond_data = combined_cond_data(binaryFilter==1,:);
combined_cond_data_ann = annotationTable(binaryFilter==1,:);

combined_ions_per_tissue = zeros(length(pca_tissue_unique),5);
for tissue_i = 1:length(pca_tissue_unique)
    cursamples = ismember(combined_cond_tissue, pca_tissue_unique{tissue_i});
    curdata_binary = ~isnan(combined_cond_data(:,cursamples));
    curdata_diet = combined_cond_diet(cursamples);
    curdata_type = combined_cond_type(cursamples);
    
    % total ions per tissue
    combined_ions_per_tissue(tissue_i,1) = ...
                sum(sum(curdata_binary,2) > 0);
    % shared ions across conditions
    combined_ions_per_tissue(tissue_i,2) = ...
        sum((sum(curdata_binary(:,ismember(curdata_type,pca_type_unique{2})),2) > 0) &...
            (sum(curdata_binary(:,ismember(curdata_type,pca_type_unique{1})),2) > 0 ));
    % type 1
    combined_ions_per_tissue(tissue_i,3) = ...
        sum((sum(curdata_binary(:,ismember(curdata_type,pca_type_unique{2})),2) == 0) &...
            (sum(curdata_binary(:,ismember(curdata_type,pca_type_unique{1})),2) <= 2) &...
            (sum(curdata_binary(:,ismember(curdata_type,pca_type_unique{1})),2) > 0 ));
    % type 2
    combined_ions_per_tissue(tissue_i,4) = ...
        sum((sum(curdata_binary(:,ismember(curdata_type,pca_type_unique{1})),2) == 0) &...
            (sum(curdata_binary(:,ismember(curdata_type,pca_type_unique{2})),2) <= 2) &...
            (sum(curdata_binary(:,ismember(curdata_type,pca_type_unique{2})),2) > 0 ));
end
combined_ions_per_tissue(:,5) = combined_ions_per_tissue(:,1) -...
                                 sum(combined_ions_per_tissue(:,2:end),2);

figure
barh(combined_ions_per_tissue(:,2:end), 'stacked')
set(gca, 'YTick', 1:length(pca_tissue_unique))
set(gca, 'YTickLabel', (pca_tissue_unique))
legend([{'Shared'}; pca_type_unique'],...
       'location', 'EastOutside')
axis square
axis ij
title('Number of annotated ions per tissue, n>=3')
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
        [figureFolder...
        'fig2a_barh_annotated_ions_per_tissue'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot number of ions shared across tissues
tissue_ion_overlap = zeros(length(pca_tissue_unique));
tissue_ion_overlap_fraction = zeros(length(pca_tissue_unique));
for tissue_i = 1:length(pca_tissue_unique)
    for tissue2 = tissue_i:length(pca_tissue_unique)
        cursamples = ismember(combined_cond_tissue, pca_tissue_unique{tissue_i});
        curdata_tissue1 = sum(sum(~isnan(combined_cond_data(:,cursamples)),2)>0);
        cursamples2 = ismember(combined_cond_tissue, pca_tissue_unique{tissue2});
        curdata_overlap = sum((sum(~isnan(combined_cond_data(:,cursamples)),2)>0) &...
                              (sum(~isnan(combined_cond_data(:,cursamples2)),2)>0));
        tissue_ion_overlap(tissue_i, tissue2) = curdata_overlap;
        tissue_ion_overlap_fraction(tissue_i, tissue2) = curdata_overlap/curdata_tissue1;
    end
end
figure
heatmap(pca_tissue_unique, pca_tissue_unique,tissue_ion_overlap_fraction)
colormap(flipud(gray))
caxis([0 1])

title('Fraction of overlapping annotated ions')
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
       [figureFolder...
       'fig_sup_heatmap_fraction_of_annotated_ion_overlap_per_tissue_updated_metfilter'])
% plot number of ions
figure
heatmap(pca_tissue_unique, pca_tissue_unique,tissue_ion_overlap)
colormap(flipud(gray))
caxis([0 2000])
title('Number of overlapping annotated ions')
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
        [figureFolder...
        'fig_sup_heatmap_number_of_annotated_ion_overlap_per_tissue_updated_metfilter'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform PCA analysis
pca_gender = combinedGender;

% perform PCA on annotated ions only 
pcamatrix = combinedIntensitiesNorm;
pcamatrix = pcamatrix(binaryFilter==1,:);

% transform intensities to log10 scale
pcamatrix = log10(pcamatrix);
% remove all empty rows, if any
pcamatrix(sum(pcamatrix==log10(5000),2) == size(pcamatrix,2),:)=[];
 
%Find the principal components
[coeff,score,latent,tsquared,explained] = pca(pcamatrix);

% define colors for different conditions
mycolors = [0 115 178;... %dark blue
            204 227 240;...%light blue
            211 96 39;... %dark orange
            246 223 212]/256;%light orange
pca_tissue_unique = unique(pca_tissue)';
pca_type_unique = unique(pca_type);
pca_diet_unique = unique(pca_diet);
pca_tissue_unique = [pca_tissue_unique(cellfun(@(x) contains(lower(x), 'serum'),pca_tissue_unique));...
                     pca_tissue_unique(cellfun(@(x) contains(lower(x), 'liver'),pca_tissue_unique));...
                     sort(pca_tissue_unique(cellfun(@(x) contains(lower(x), 'si'),pca_tissue_unique)));...
                     sort(pca_tissue_unique(cellfun(@(x) contains(lower(x), 'c'),pca_tissue_unique)))];
% define markers for different tissues
mymarkers = {'p','s','d', 'v','<','o', '>','^'};
markersize = 100;

figure;
hold on
sch = [];
coloridx=1;
for diet_i = 1:length(pca_diet_unique)
    for type_i = 1:length(pca_type_unique)
        for tissue_i = 1:length(pca_tissue_unique)
            cursamples = ismember(pca_diet, pca_diet_unique{diet_i}) &...
                         ismember(pca_type, pca_type_unique{type_i}) &...
                         ismember(pca_tissue, pca_tissue_unique{tissue_i});
            switch type_i
                case 3
                    scatter(coeff(cursamples,1),coeff(cursamples,2),...
                            markersize,...
                            mymarkers{tissue_i},...
                            'MarkerFaceColor', mycolors(coloridx,:),...
                            'MarkerEdgeColor', 'none');%'k');
                case 1
                    sch(tissue_i) = scatter(coeff(cursamples,1),coeff(cursamples,2),...
                            markersize,...
                            mymarkers{tissue_i},...
                            'MarkerFaceColor', mycolors(coloridx,:),...
                            'MarkerEdgeColor', 'none');%mycolors(coloridx,:));
                case 2
                    scatter(coeff(cursamples,1),coeff(cursamples,2),...
                            markersize,...
                            mymarkers{tissue_i},...
                            'MarkerFaceColor', mycolors(coloridx,:),...[1 1 1],...
                            'MarkerEdgeColor', 'none');%mycolors(coloridx,:));
            end
        end
        coloridx = coloridx+1;
    end
end
% plot std ellipses around tissues
tisue_color = jet;
tisue_color = flipud(tisue_color(1:round(length(tisue_color)/length(pca_tissue_unique)):end,:));
for tissue_i = 1:length(pca_tissue_unique)
    cursamples = ismember(pca_tissue, pca_tissue_unique{tissue_i});
    tc(tissue_i) = plot_error_ellipse([coeff(cursamples,1), coeff(cursamples,2)],...
                            2.4477, tisue_color(tissue_i,:));%, 2.6)
end
        
xlabel(sprintf('1st PC %.2f', (explained(1))));
ylabel(sprintf('2st PC %.2f', (explained(2))));
ylim([-0.15 0.2])
xlim([0.05 0.12])
axis square
legend([sch tc], [pca_tissue_unique;pca_tissue_unique], 'Location', 'EastOutside')
title('PCA of metabolomics samples, annotated ions, log10intensities, Impute5000')
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
         [figureFolder...
         'fig2b_pca_metabolomics_all_annotated_ions_updated_metfilter_lognorm_Impute5000'])

   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% perfrom PCA per tissue
cursamples = ismember(pca_tissue, pca_tissue_unique{tissue_i});
coeff_tissues1 = nan(nnz(cursamples), length(pca_tissue_unique));
coeff_tissues2 = nan(nnz(cursamples), length(pca_tissue_unique));
explained_tissues = zeros(2, length(pca_tissue_unique));

figure;
for tissue_i = 1:length(pca_tissue_unique)
    cursamples = ismember(pca_tissue, pca_tissue_unique{tissue_i});
    
    pcamatrix_tissue = pcamatrix(:, cursamples);
    % remove rows with more than half nan values
    pcamatrix_tissue(sum(pcamatrix_tissue==log10(5000),2) == size(pcamatrix_tissue,2),:)=[];
     
    pca_diet_tissue = pca_diet(cursamples);
    pca_type_tissue = pca_type(cursamples);
    
    [coeff,~,~,~,explained] = pca(pcamatrix_tissue);
    
    % save coeff for further analysis
    coeff_tissues1(1:size(coeff,1), tissue_i) = coeff(:,1);
    coeff_tissues2(1:size(coeff,1), tissue_i) = coeff(:,2);
    explained_tissues(:, tissue_i) = explained(1:2);
    
    hold on
    sch = [];
    sch_legend = {}; 
    schidx=1;
    coloridx=1;
    for diet_i = 1:length(pca_diet_unique)
        for type_i = 1:length(pca_type_unique)
            cursamples = ismember(pca_diet_tissue, pca_diet_unique{diet_i}) &...
                         ismember(pca_type_tissue, pca_type_unique{type_i});
            switch type_i
                case 3
                    scatter(coeff(cursamples,1),coeff(cursamples,2),...
                            markersize,...
                            mymarkers{tissue_i},...
                            'MarkerFaceColor', mycolors(coloridx,:),...
                            'MarkerEdgeColor', 'none');%'k');
                case 1
                    sch(schidx) = scatter(coeff(cursamples,1),coeff(cursamples,2),...
                            markersize,...
                            mymarkers{tissue_i},...
                            'MarkerFaceColor', mycolors(coloridx,:),...
                            'MarkerEdgeColor', 'none');%mycolors(diet_i,:));
                case 2
                    sch(schidx) = scatter(coeff(cursamples,1),coeff(cursamples,2),...
                            markersize,...
                            mymarkers{tissue_i},...
                            'MarkerFaceColor', mycolors(coloridx,:),...%[1 1 1],...
                            'MarkerEdgeColor', 'none');%mycolors(diet_i,:));
            end
            plot_error_ellipse([coeff(cursamples,1), coeff(cursamples,2)],...
                                2.4477, 'k');%, 2.6)

            sch_legend = [sch_legend;
                          strcat(pca_diet_unique(diet_i), ...
                                        pca_type_unique(type_i))];
            schidx = schidx + 1;
            coloridx = coloridx+1;
        end
    end
    xlabel(sprintf('1st PC %.2f', (explained(1))));
    ylabel(sprintf('2st PC %.2f', (explained(2))));
    title(pca_tissue_unique(tissue_i))
    
    scalecoef = 100*((-log10(max(coeff(:,1)) - min(coeff(:,1))))>1) + ...
                10*((-log10(max(coeff(:,1)) - min(coeff(:,1))))<1);    
    xlim([floor(min(coeff(:,1))*scalecoef)/scalecoef ...
          ceil(max(coeff(:,1))*scalecoef)/scalecoef])

    scalecoef = 100*((-log10(max(coeff(:,2)) - min(coeff(:,2))))>1) + ...
                10*((-log10(max(coeff(:,2)) - min(coeff(:,2))))<1);    
    ylim([floor(min(coeff(:,2))*scalecoef)/scalecoef ...
          ceil(max(coeff(:,2))*scalecoef)/scalecoef])
    
    axis square
    legend(sch, sch_legend, 'Location', 'EastOutside')
    print(gcf, '-painters', '-dpsc', '-r600', '-bestfit', ...
              '-append',...
              [figureFolder...
              'fig2b_sup_pca_metabolomics_per_tissue_annotated_ions_updated_metfilter_lognorm_Impute5000'])
    clf
end
    