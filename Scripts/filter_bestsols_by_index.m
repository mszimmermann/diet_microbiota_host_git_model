function met_bestsols = filter_bestsols_by_index(met_bestsols, filter_index)

met_bestsols.x_sel_CorrRev = met_bestsols.x_sel_CorrRev(filter_index,:);
met_bestsols.x_sel_CorrRevSI = met_bestsols.x_sel_CorrRevSI(filter_index,:);
met_bestsols.x_sel_CorrRevLI = met_bestsols.x_sel_CorrRevLI(filter_index,:);
met_bestsols.x_sel_CorrRevMean = met_bestsols.x_sel_CorrRevMean(filter_index,:);
met_bestsols.x = met_bestsols.x(filter_index,:);

met_bestsols.kmeanMatrix_joint_orig = met_bestsols.kmeanMatrix_joint_orig(filter_index,:);
met_bestsols.x_sel_dataR = met_bestsols.x_sel_dataR(filter_index,:);
met_bestsols.selection_criteria = met_bestsols.selection_criteria(filter_index,:);


