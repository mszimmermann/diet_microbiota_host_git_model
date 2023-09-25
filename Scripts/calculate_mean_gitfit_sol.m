function [x_meansol] = calculate_mean_gitfit_sol(gitfit)
% calculate mean solution among the ones with correlation
% above corr_threshold
% report mean and std of the normalized solution
% and the number of supporting solutions
corr_threshold = 0.7; % for subsetting

% mean total correlations
testmean_tot = gitfit.testx(:, gitfit.testCorrRev>=corr_threshold);
if ~isempty(testmean_tot)
    %nomalize by max
    testmean_tot = testmean_tot(2:end,:)./abs(max(testmean_tot(2:end,:)));
    x_mean_tot = mean(testmean_tot,2);
    x_std_tot = std(testmean_tot,[],2);
    x_sup_tot = size(testmean_tot,2);
else
    x_mean_tot = zeros(size(gitfit.testx,1),1);
    x_std_tot = zeros(size(gitfit.testx,1),1);
    x_sup_tot = 0;
end

% max LI correlations within high total correlations
subset_total = (gitfit.testCorrRev>=corr_threshold &...
                    gitfit.testCorrRevLI>=corr_threshold);
testmean_LI_wtot = gitfit.testx(:, subset_total);
if ~isempty(testmean_LI_wtot)
    %nomalize by max
    testmean_LI_wtot = testmean_LI_wtot(2:end,:)./abs(max(testmean_LI_wtot(2:end,:)));
    x_mean_LI_wtot = mean(testmean_LI_wtot,2);
    x_std_LI_wtot = std(testmean_LI_wtot,[],2);
    x_sup_LI_wtot = size(testmean_LI_wtot,2);
else
    x_mean_LI_wtot = zeros(size(gitfit.testx,1),1);
    x_std_LI_wtot = zeros(size(gitfit.testx,1),1);
    x_sup_LI_wtot = 0;
end

% merge selected solutions
x_mean = [x_mean_tot x_mean_LI_wtot];
x_std = [x_std_tot x_std_LI_wtot];
x_support = [x_sup_tot x_sup_LI_wtot];
selection_criterion = {'mean total PCC', 'mean LI PCC within high total'};

x_meansol.x_mean = x_mean;
x_meansol.x_std = x_std;
x_meansol.x_support = x_support;
x_meansol.selection_criterion = selection_criterion;
x_meansol.corr_threshold = corr_threshold;
% add original metabolite data to the solution
x_meansol.kmeanMatrix_joint_orig = gitfit.kmeanMatrix_joint_orig;
x_meansol.kmeanMatrix_joint_names = gitfit.kmeanMatrix_joint_names;
% add coefficient names to the object
x_meansol.coefvalues = gitfit.coefvalues;