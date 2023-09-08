function [x_selected] = select_gitfit_sol(gitfit)
% select gitfit solution
corr_threshold = 0.7; % for subsetting

% max total correlations
[testmax_tot, testargmax_tot] = max(gitfit.testCorrRev);
x_tot = gitfit.testx(:,testargmax_tot);

% max SI correlations
[testmax_SI, testargmax_SI] = max(gitfit.testCorrRevSI);
x_SI = gitfit.testx(:,testargmax_SI);

% max LI correlations
[testmax_LI, testargmax_LI] = max(gitfit.testCorrRevLI);
x_LI = gitfit.testx(:,testargmax_LI);

% max sum correlations
[testmax_sum, testargmax_sum] = max(gitfit.testCorrRev + ...
                            gitfit.testCorrRevSI + ...
                            gitfit.testCorrRevLI);
x_sum = gitfit.testx(:,testargmax_sum);

% max SI correlations within high total correlations
subset_total = find(gitfit.testCorrRev>=corr_threshold);
[testmax_SI_wtot, testargmax_SI_wtot] = max(gitfit.testCorrRevSI(subset_total));
testargmax_SI_wtot = subset_total(testargmax_SI_wtot);
x_SI_wtot = gitfit.testx(:,testargmax_SI_wtot);

% max LI correlations within high total correlations
[testmax_LI_wtot, testargmax_LI_wtot] = max(gitfit.testCorrRevLI(subset_total));
testargmax_LI_wtot = subset_total(testargmax_LI_wtot);
x_LI_wtot = gitfit.testx(:,testargmax_LI_wtot);

% merge selected solutions
x_sel = [gitfit.x_ip x_tot x_SI x_LI x_sum x_SI_wtot x_LI_wtot];
selection_value = [gitfit.x_ip_CorrRev testmax_tot testmax_SI testmax_LI testmax_sum,...
                    testmax_SI_wtot testmax_LI_wtot];
selection_criterion = {'IP' 'total PCC', 'SI PCC', 'LI PCC', 'sum PCC',...
                       'SI PCC within high total',...
                       'LI PCC within high total'};
selection_arg = [1 testargmax_tot testargmax_SI testargmax_LI testargmax_sum,...
                testargmax_SI_wtot testargmax_LI_wtot];
x_sel_CorrRev = [gitfit.x_ip_CorrRev gitfit.testCorrRev(selection_arg(2:end))];
x_sel_CorrRevSI = [gitfit.x_ip_CorrRevSI gitfit.testCorrRevSI(selection_arg(2:end))];
x_sel_CorrRevLI = [gitfit.x_ip_CorrRevLI gitfit.testCorrRevLI(selection_arg(2:end))];
% add resiprocal data for each selected solution
x_sel_dataR = [gitfit.dataR_ip(:) gitfit.testdataR(:,selection_arg(2:end))];


x_selected.x = x_sel;
x_selected.selection_value = selection_value;
x_selected.selection_criterion = selection_criterion;
x_selected.selection_arg = selection_arg;
x_selected.x_sel_CorrRev = x_sel_CorrRev;
x_selected.x_sel_CorrRevSI = x_sel_CorrRevSI;
x_selected.x_sel_CorrRevLI = x_sel_CorrRevLI;
x_selected.x_sel_dataR = x_sel_dataR;    
% add coefficient names to the object
x_selected.coefvalues = gitfit.coefvalues;
              