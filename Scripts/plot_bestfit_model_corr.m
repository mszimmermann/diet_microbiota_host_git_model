function plot_bestfit_model_corr(met_bestfit, met_gitfits, filename)

% plot correlation of restored and original data
% and save plot to file filename
% between best solution and random solutions

% calculate difference in corr distrbutions
% get best correlation for each solution
x_data_corr = met_bestfit.x_sel_CorrRev;
% get best correlation for each random solution
x_data_corr_shuffled_best = cellfun(@(y)max(y.testCorrRev_shuffled),met_gitfits);

% get random random solution
x_data_corr_shuffled = cellfun(@(y) y.testCorrRev_shuffled(randi(length(y.testCorrRev_shuffled))),met_gitfits);

% compare corr distributions
p_corr_diff = ranksum(x_data_corr_shuffled, x_data_corr);
figure
nbins=10;
% pearson corr all 
histogram(x_data_corr_shuffled, nbins);
hold on
h = histogram(x_data_corr, nbins);
xlim([-1 1])
axis square
xlabel('MAX PCC between metabolomics data and model estimate')
ylabel('Number of ions')
orient landscape
% print line for PCC=0.7
plot([0.7, 0.7], [0, max(h.Values)], 'k--', 'LineWidth',2)

legend({'Random coefficients', 'Model coefficients', 'PCC=0.7'},...
        'Location', 'NorthWest')
    
% compare distributions of correlations
pval = signrank(x_data_corr_shuffled,...
                x_data_corr);
% print Wilcoxon signed rank test o-value on the plot
text(-0.9, 0.7*max(h.Values), sprintf('signrank p = %.2e', pval))
text(-0.9, 0.6*max(h.Values), sprintf('ranksum p = %.2e', p_corr_diff))
text(-0.9, 0.5*max(h.Values), sprintf('median(model) = %.2f', median(x_data_corr)))
text(-0.9, 0.4*max(h.Values), sprintf('median(random) = %.2f', median(x_data_corr_shuffled)))

% save figure to file
print(gcf, '-vector', '-dpdf', '-r600', '-bestfit', ...
      [filename, '_random_random_sol'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot best random solutions
% compare corr distributions
p_corr_diff = ranksum(x_data_corr_shuffled_best, x_data_corr);
figure
nbins=10;
% pearson corr all 
histogram(x_data_corr_shuffled_best, nbins);
hold on
h = histogram(x_data_corr, nbins);
xlim([-1 1])
axis square
xlabel('MAX PCC between metabolomics data and model estimate')
ylabel('Number of ions')
orient landscape
% print line for PCC=0.7
plot([0.7, 0.7], [0, max(h.Values)], 'k--', 'LineWidth',2)

legend({'Best random coefficients', 'Model coefficients', 'PCC=0.7'},...
        'Location', 'NorthWest')
    
% compare distributions of correlations
pval = signrank(x_data_corr_shuffled_best,...
                x_data_corr);
% print Wilcoxon signed rank test o-value on the plot
text(-0.9, 0.7*max(h.Values), sprintf('signrank p = %.2e', pval))
text(-0.9, 0.6*max(h.Values), sprintf('ranksum p = %.2e', p_corr_diff))
text(-0.9, 0.5*max(h.Values), sprintf('median(model) = %.2f', median(x_data_corr)))
text(-0.9, 0.4*max(h.Values), sprintf('median(random) = %.2f', median(x_data_corr_shuffled_best)))

% save figure to file
print(gcf, '-vector', '-dpdf', '-r600', '-bestfit', ...
      [filename, '_best_random_sol'])

