%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot diet formulations

% call script defining file dependenciesand global variables
addpath(genpath('.\'));
add_global_and_file_dependencies

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output: 
% Figures:
% 'fig_sup_diet_formulation_comparison'
% 'fig_sup_diet_formulation_comparison_zoomed.pdf'
% 'fig_sup_diet_nutrientinfo_comparison'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dietFormulation_HFD = {'Casein' 265.0;...
    'L-Cystine' 4.0;...
    'Maltodextrin' 160.0;...
    'Sucrose' 90.0;...
    'Lard' 310.0;...
    'Soybean Oil' 30.0;...
    'Cellulose' 65.5;...
    'Mineral Mix AIN-93G-MX (94046)' 48.0;...
    'Calcium Phosphate, dibasic' 3.4;...
    'Vitamin Mix, AIN-93-VX (94047)' 21.0;...
    'Choline Bitartrate' 3.0;...
    'Blue Food Color' 0.1};

dietFormulation_CTR = {'Casein' 140.0;...
    'L-Cystine' 1.8;...
    'Corn Starch' 460.427;...
    'Maltodextrin' 155.0;...
    'Sucrose' 100.0;...
    'Soybean Oil' 40.0;...
    'Cellulose' 50.0;...
    'Mineral Mix, AIN-93M-MX (94049)' 35.0;...
    'Vitamin Mix, AIN-93-VX (94047)' 15.0;...
    'Thiamin (81%)' 0.013;...
    'Vitamin K1, phylloquinone' 0.002;...
    'Choline Bitartrate' 2.75;...
    'TBHQ, antioxidant' 0.008};

% combine the two in one plot
dietFormulation_combined = unique([dietFormulation_HFD(:,1); dietFormulation_CTR(:,1)]);
dietFormulation_combined_numbers = zeros(length(dietFormulation_combined),2);
for i=1:length(dietFormulation_combined)
    idx = find(ismember(dietFormulation_CTR(:,1), dietFormulation_combined{i}));
    if ~isempty(idx)
        dietFormulation_combined_numbers(i,1) = dietFormulation_CTR{idx,2};
    end
    idx = find(ismember(dietFormulation_HFD(:,1), dietFormulation_combined{i}));
    if ~isempty(idx)
        dietFormulation_combined_numbers(i,2) = dietFormulation_HFD{idx,2};
    end
end
dietlabels = {'HCD','HFD'};
% sort by max
[~, sortidx] = sort(max(dietFormulation_combined_numbers,[],2));
figure
barh((dietFormulation_combined_numbers(sortidx,:)))
set(gca, 'YTick', 1:length(dietFormulation_combined));
set(gca, 'YTickLabel', dietFormulation_combined(sortidx));
legend(dietlabels)
title('HFD and HCD diet comparison')
xlabel('g/Kg')
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
    [figureFolder 'fig_sup_diet_formulation_comparison'])

% plot with splitting low and high concentrations
fig = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,1,1)
barh(dietFormulation_combined_numbers(sortidx,:))
set(gca, 'YTick', 1:length(dietFormulation_combined));
set(gca, 'YTickLabel', dietFormulation_combined(sortidx));
legend(dietlabels, 'Location', 'SouthEast')
title('HFD and HCD diet comparison')
xlabel('g/Kg')
xlim([0 10])
ylim([0 7.5])
axis square

subplot(2,1,2)
barh(dietFormulation_combined_numbers(sortidx,:))
set(gca, 'YTick', 1:length(dietFormulation_combined));
set(gca, 'YTickLabel', dietFormulation_combined(sortidx));
legend(dietlabels, 'Location', 'SouthEast')
title('HFD and HCD diet comparison')
xlabel('g/Kg')
xlim([0 500])
ylim([7.5 17.5])
axis square
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
    [figureFolder, 'fig_sup_diet_formulation_comparison_zoomed.pdf'])


dietNutrientInfo = {'Protein'; 'Carbohydrate'; 'Fat'};
dietNutrientPercentByWeight = [12.4 23.5;...
                               68.4 27.3;...
                               4.1 34.3];
dietNutrientPercentCaloriesFrom = [13.7 18.3;...
                                   75.9 21.4;...
                                   10.3 60.3];
figure
subplot(1,2,1)
barh(dietNutrientPercentByWeight)
set(gca, 'YTick', 1:length(dietNutrientInfo));
set(gca, 'YTickLabel', dietNutrientInfo);
legend(dietlabels, 'Location', 'SouthEast')
title('Nutrient percent by weight')
xlabel('%')
axis square

subplot(1,2,2)
barh(dietNutrientPercentCaloriesFrom)
set(gca, 'YTick', 1:length(dietNutrientInfo));
set(gca, 'YTickLabel', dietNutrientInfo);
legend(dietlabels, 'Location', 'SouthEast')
title('Nutrient percent calories')
xlabel('%')
suptitle('HFD and CTR diet comparison')
axis square
orient landscape
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
    [figureFolder 'fig_sup_diet_nutrientinfo_comparison'])
   

    