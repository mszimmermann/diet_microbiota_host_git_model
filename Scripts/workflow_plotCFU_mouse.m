%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot CFU per mouse

% call script defining file dependenciesand global variables
addpath(genpath('.\'));
add_global_and_file_dependencies

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: 
% mouse_CFU_DC.csv
% mouse_fecal_weights.csv
% Output: 
% Figures:
% 'fig_sup_scatter_mean_mouseCFU_diet.pdf'
% 'fig_sup_scatter_cfu_per_g_feces_per_mouse.pdf'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate CFU/mg fecal mas 
liquidV = 4; %uL of plated cellular volume
initialVolume = 1000; %total uL

DCcfu_table = readtable([rawdataFolder,...
                         'mouse_CFU_DC.csv']);
mouseFecalWeights = readtable([rawdataFolder,...
    'mouse_fecal_weights.csv']);
DCcfu = DCcfu_table{3:end, 2:end};
DC_mouse = DCcfu_table{3:end,1};
DC_dilution = DCcfu_table{1,2:end};

% calculate CFUs
CURcfu = DCcfu;
CURcfu_mouse = DC_mouse;
CURcfu_dilution = DC_dilution;

CURcfu_calculated = CURcfu*initialVolume/liquidV;
for i=1:size(CURcfu_calculated,2)
    CURcfu_calculated(:,i) = CURcfu_calculated(:,i).*prod(CURcfu_dilution(1:i));
end
for i=1:size(CURcfu_calculated,1)
    CURcfu_calculated(i,:) = CURcfu_calculated(i,:)*1000./...
                             mouseFecalWeights{ismember(mouseFecalWeights{:,1},...
                             CURcfu_mouse(i)),2};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate mean CFU per mouse
curmouseNumbers = unique(CURcfu_mouse);

curmouseMeanCFU = zeros(size(curmouseNumbers));
curmouseSTDCFU = zeros(size(curmouseNumbers));
for i=1:length(curmouseNumbers)
    curCFU = CURcfu_calculated(ismember(CURcfu_mouse,curmouseNumbers(i)),:);
    curCFU = curCFU(:);
    % remove nonmeasured zeros
    curCFU(curCFU==0) = [];
    curCFU = log10(curCFU);
    curmouseMeanCFU(i) = mean(curCFU(:), "omitnan");
    curmouseSTDCFU(i) = std(curCFU(:),"omitnan");
end
% plot CFUs per group
figure
scatter(rand(1,5), curmouseMeanCFU(6:end), 'filled')
hold on
plot([0 1], [mean(curmouseMeanCFU(6:end)) mean(curmouseMeanCFU(6:end))], 'k')
scatter(2+rand(1,5), curmouseMeanCFU(1:5), 'filled')
plot([2 3], [mean(curmouseMeanCFU(1:5)) mean(curmouseMeanCFU(1:5))], 'k')

ylim([5 15])
xlim([-1 4])
axis square
set(gca, 'XTick', [0.5 2.5])
set(gca, 'XTickLabel', {'CTR', 'HFD'})
[~, pvalue] = ttest2(curmouseMeanCFU(1:5),curmouseMeanCFU(6:end));
text(1.5,13, sprintf('p=%.3f',pvalue))
ylabel('CFU/g feces, log10')
title('Bacterial load after 4 weeks of diet')
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
    [figureFolder 'fig_sup_scatter_mean_mouseCFU_diet.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT cfu PER MOUSE
figure;
curmouseLabels = curmouseNumbers;
% plot first control and thenb HFD
curmouseLabels = curmouseLabels([6:10 1:5]);
for i=1:length(curmouseNumbers)
    curCFU = CURcfu_calculated(ismember(CURcfu_mouse,curmouseLabels(i)),:);
    curCFU = curCFU(:);
    curCFU(curCFU==0) = [];
    scatter(i*10+5*rand(size(curCFU)), log10(curCFU), 'filled')
    hold on
    plot([i*10 i*10+5], [mean(log10(curCFU),"omitnan") ...
                         mean(log10(curCFU), "omitnan")], 'k')
end

set(gca, 'XTick', 12.5:10:110)
set(gca, 'XTickLabel',  curmouseLabels)
ylabel('CFU / g feces, log10')
xlim([5 110])
ylim([5 15])
title('Mouse CFU/g feces after 4 weeks on HFD or control diet')
axis square
orient landscape
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
    [figureFolder 'fig_sup_scatter_cfu_per_g_feces_per_mouse.pdf'])