% Compare modelling results for metabolites between different studies

% call script defining file dependencies and global variables
addpath(genpath('.\'));
add_global_and_file_dependencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data requirements: 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output:
% Figures:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read annotation from file with Meier IDs
annotationTableMeier = readtable([inputFolder ...
    'metabolites_allions_combined_formulas_with_metabolite_filters_meier.csv']);

select_mets = cellfun(@(x) ~isempty(x), annotationTableMeier.Meier_compound);
%select_mets = annotationTableMeier.MetaboliteFilter==1;
% git_mets = (annotationTableSpatialClusters.spatial_clust100_CTR_DC + ...
%             annotationTableSpatialClusters.spatial_clust100_CTR_GF + ...
%             annotationTableSpatialClusters.spatial_clust100_HFD_DC + ...
%             annotationTableSpatialClusters.spatial_clust100_HFD_GF)>0;
% select_mets = (annotationTableSpatialClusters.MetaboliteFilter==1) &...
%                (git_mets>0);

met_info=[];
met_info.MZ =annotationTableMeier.MZ(select_mets);
met_info.RT =annotationTableMeier.RT(select_mets);
met_info.CompoundID = annotationTableMeier.CompoundID(select_mets);
met_info.CompoundName = annotationTableMeier.Meier_compound(select_mets);
% remove spaces and capitalize first letters
for i=1:length(met_info.CompoundName)
    str=lower(met_info.CompoundName{i});
    idx=regexp([' ' str],'(?<=\s+)\S','start')-1;
    str(idx)=upper(str(idx));
    str = strrep(str, ' ','');
    str = strrep(str, '-','_');
    met_info.CompoundName{i} = str;
end



% read CVR results
filename = [outputFolder ...
            'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_allCVR'];
sel_crit = 'LI PCC within high total';%'total PCC';%'IP';%
[met_info_cvr, met_bestsol_cvr] = read_bestsol_from_file(filename,...
                                                sel_crit);
met_bestsol_cvr.modelname = ['CVR_' sel_crit];

% read CVR results
filename = [outputFolder ...
            'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_allDC'];
sel_crit = 'total PCC';%'IP';%'LI PCC within high total';%'IP';%
[met_info_cvr, met_bestsol_cvr] = read_bestsol_from_file(filename,...
                                                sel_crit);
met_bestsol_cvr.modelname = ['DC_' sel_crit];

% keep only selected metabolites
met_bestsol_cvr.x_sel_CorrRev(~select_mets)=[];
met_bestsol_cvr.x_sel_CorrRevLI(~select_mets)=[];
met_bestsol_cvr.x_sel_CorrRevSI(~select_mets)=[];
met_bestsol_cvr.x_sel_CorrRevMean(~select_mets)=[];
met_bestsol_cvr.x(~select_mets,:)=[];
met_bestsol_cvr.kmeanMatrix_joint_orig(~select_mets,:)=[];
met_bestsol_cvr.x_sel_dataR(~select_mets,:)=[];
% also in met_info
met_info_cvr.MZ(~select_mets)=[];
met_info_cvr.RT(~select_mets)=[];
met_info_cvr.CompoundID(~select_mets)=[];
met_info_cvr.CompoundName(~select_mets)=[];
met_info_cvr.MetaboliteFilter(~select_mets)=[];
met_info_cvr.SumGITclusters(~select_mets)=[];

[confmat, classlabels] = compare_two_model_results(met_info_dc, met_bestsol_dc,...
                          met_info_cvr, met_bestsol_cvr,...
                          figureFolder);


                      
        x_met_smooth = modelingResults{:, width(modelingResults)-8:end};
coefvalues = modelingResults.Properties.VariableNames(width(modelingResults)-8:end);

met_bestsol_DC.coefvalues = coefvalues;
met_bestsol_DC.x = x_met_smooth(select_mets,:);
%met_bestsol_DC.ReciprocalCorr = modelingResults.ReciprocalCorr(select_mets);
met_bestsol_DC.x_sel_CorrRev = cellfun(@(x) str2double(x),...
    modelingResults.ReciprocalCorr(select_mets));
met_bestsol_DC.x_sel_CorrRevSI = modelingResults.ReciprocalCorr(select_mets);
met_bestsol_DC.x_sel_CorrRevLI = modelingResults.ReciprocalCorr(select_mets);
met_bestsol_DC.modelname = 'Old_IP_DC';
met_bestsol_DC.selection_criterion = 'IP';


modelingResults = readtable([outputFolder...
            'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_allions_with_CVR.csv']);
x_met_smooth_CVR = modelingResults{:, width(modelingResults)-8:end};
coefvalues_CVR = modelingResults.Properties.VariableNames(width(modelingResults)-8:end);

met_bestsol_CVR.coefvalues = coefvalues_CVR;
met_bestsol_CVR.x = x_met_smooth_CVR(select_mets,:);
met_bestsol_CVR.x_sel_CorrRev = modelingResults.ReciprocalCorr(select_mets);
met_bestsol_CVR.x_sel_CorrRevSI = modelingResults.ReciprocalCorr(select_mets);
met_bestsol_CVR.x_sel_CorrRevLI = modelingResults.ReciprocalCorr(select_mets);
met_bestsol_CVR.modelname = 'Old_IP_CVR';
met_bestsol_CVR.selection_criterion = 'IP';

% load data and restored data from file
% save model results to file - reciprocal data restoration
modelData = readtable([resultsFolder...
    'model_results_SMOOTH_normbyabsmax_reciprocal_problem_allions.csv']);
modelData_data = modelData(:, 9:end);
modelData_orig = modelData_data{:, cellfun(@(x) ~(contains(x, 'Recip') |...
                                             contains(x, 'Random') ),...
                                             modelData_data.Properties.VariableNames)};
modelData_orig_cols = modelData_data.Properties.VariableNames(cellfun(@(x) ~(contains(x, 'Recip') |...
                                             contains(x, 'Random') ),...
                                             modelData_data.Properties.VariableNames));
modelData_recip = modelData_data{:, cellfun(@(x) contains(x, 'Recip'),...
                                             modelData_data.Properties.VariableNames)};
% get correlations calculated with reverse problem
if isnumeric(modelData.ReciprocalCorr(1))
    x_data_corr = modelData.ReciprocalCorr;
else
    % it is not numeric, probably contains NaN - convert to numeric
    x_data_corr = cellfun(@(x) str2double(x), modelData.ReciprocalCorr);
end

met_bestsol_DC.kmeanMatrix_joint_names = modelData_data.Properties.VariableNames(1:24);
met_bestsol_DC.kmeanMatrix_joint_orig = modelData_data{select_mets, 1:24};
met_bestsol_DC.x_sel_dataR = modelData_data{select_mets, 25:48};
%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare Meier and DC/CVR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read DC results
filename = [outputFolder ...
            'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_allCVR'];
sel_crit = 'IP';%'total PCC';%'LI PCC within high total';%
[met_info_dc, met_bestsol_dc] = read_bestsol_from_file(filename,...
                                                sel_crit);
% combine IP and LI within PCC criteria
sel_crit1 = 'IP';
sel_crit2 = 'LI PCC within high total';%
[met_info_dc_combined, met_bestsol_dc_combined] = ...
                combine_bestsols_from_file(filename, sel_crit1, sel_crit2);
met_bestsol_dc_combined.modelname = 'CVR_IP_LIPCCwT';%

% set compound name to Meier
met_info_dc.CompoundName = annotationTableMeier.Meier_compound;
met_bestsol_dc.modelname = ['CVR_' sel_crit];
% keep only selected metabolites
met_bestsol_dc.x_sel_CorrRev(~select_mets)=[];
met_bestsol_dc.x_sel_CorrRevLI(~select_mets)=[];
met_bestsol_dc.x_sel_CorrRevSI(~select_mets)=[];
met_bestsol_dc.x_sel_CorrRevMean(~select_mets)=[];
met_bestsol_dc.x(~select_mets,:)=[];
met_bestsol_dc.kmeanMatrix_joint_orig(~select_mets,:)=[];
met_bestsol_dc.x_sel_dataR(~select_mets,:)=[];
% also in met_info
met_info_dc.MZ(~select_mets)=[];
met_info_dc.RT(~select_mets)=[];
met_info_dc.CompoundID(~select_mets)=[];
met_info_dc.CompoundName(~select_mets)=[];
met_info_dc.MetaboliteFilter(~select_mets)=[];
met_info_dc.SumGITclusters(~select_mets)=[];
% remove spaces and capitalize first letters
for i=1:length(met_info_dc.CompoundName)
    str=lower(met_info_dc.CompoundName{i});
    idx=regexp([' ' str],'(?<=\s+)\S','start')-1;
    str(idx)=upper(str(idx));
    str = strrep(str, ' ','');
    str = strrep(str, '-','_');
    met_info_dc.CompoundName{i} = str;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read meier modelling results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = '.\ProcessedData\model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_Meier';
sel_crit = 'IP';%'LI PCC within high total';%
[met_info_meier, met_bestsol_meier] = read_bestsol_from_file(filename,...
                                                sel_crit);
                                            
% combine IP and LI within PCC criteria
sel_crit1 = 'IP';
sel_crit2 = 'LI PCC within high total';%
[met_info_meier_combined, met_bestsol_meier_combined] = ...
                combine_bestsols_from_file(filename, sel_crit1, sel_crit2);
met_bestsol_meier_combined.modelname = 'Meier_IP_LIPCCwT';%
            
%edit meier names to remove x in front of metabolites
for i=1:length(met_info_meier.CompoundName)
    curname = met_info_meier.CompoundName{i};
    if curname(1)=='x'
        if isstrprop(curname(2),'digit')
            curname = curname(2:end);
            met_info_meier.CompoundName{i} = curname;
        end
    end
end
met_bestsol_meier.modelname = 'Meier_IP';%'Meier_LIPCCwT';%

[confmat, classlabels] = compare_two_model_results(met_info_dc, met_bestsol_dc,...
                          met_info_meier, met_bestsol_meier,...
                          figureFolder,0);

[confmat, classlabels] = compare_two_model_results(met_info_dc, met_bestsol_dc_combined,...
                          met_info_meier, met_bestsol_meier_combined,...
                          figureFolder,1);

