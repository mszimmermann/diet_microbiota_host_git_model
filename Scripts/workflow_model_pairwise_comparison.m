% Compare modelling results for metabolites between different studies

% call script defining file dependencies and global variables
addpath(genpath('.\'));
add_global_and_file_dependencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data requirements: 
% metabolite annotation with mapping to metabolites from the Meier et al 2023 study
%'metabolites_allions_combined_formulas_with_metabolite_filters_meier.csv
% modelling outputs
% 'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_bestDC'
% 'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_bestCVR'
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read CVR modelling results
filename = [outputFolder ...
            'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_bestCVR'];
corrthreshold = 0.7;
% try combining solutions
sel_crit1 = 'IP';
sel_crit2 = 'LI_PCC_within_high_total';

[met_info_combined_CVR, met_bestsol_combined_CVR] = ...
                combine_bestsols_from_file(filename, ...
                sel_crit1, sel_crit2,...
                corrthreshold);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read DC modelling results
filename = [outputFolder ...
            'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_bestDC'];
corrthreshold = 0.7;
% try combining solutions
sel_crit1 = 'IP';
sel_crit2 = 'LI_PCC_within_high_total';

[met_info_combined_DC, met_bestsol_combined_DC] = ...
                combine_bestsols_from_file(filename, ...
                sel_crit1, sel_crit2,...
                corrthreshold);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check why GIT clusers and met filter ==1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% compare all ions
select_mets = 1:size(annotationTableMeier,1);

% keep only selected metabolites
met_bestsol_combined_DC_filtered = filter_bestsols_by_index(met_bestsol_combined_DC, select_mets);
met_bestsol_combined_CVR_filtered = filter_bestsols_by_index(met_bestsol_combined_CVR, select_mets);

% also in met_info
met_info_combined_DC_filtered = filter_metinfo_by_index(met_info_combined_DC, select_mets);
met_info_combined_CVR_filtered = filter_metinfo_by_index(met_info_combined_CVR, select_mets);

% add model names for plotting and saving to file
met_bestsol_combined_DC_filtered.modelname = 'bestsol_combined_DC_allions';
met_bestsol_combined_CVR_filtered.modelname = 'bestsol_combined_CVR_allions';

[confmatDCCVRallions, classlabelsDCCVRallions] = compare_two_model_results(...
                          met_info_combined_DC_filtered,...
                          met_bestsol_combined_DC_filtered,...
                          met_info_combined_CVR_filtered,...
                          met_bestsol_combined_CVR_filtered,...
                          figureFolder, 1,... %flag strictclass
                          1); %compare by mzrt and not compound name

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare only annotated ions
select_mets = annotationTableMeier.MetaboliteFilter==1;

% keep only selected metabolites
met_bestsol_combined_DC_filtered = filter_bestsols_by_index(met_bestsol_combined_DC, select_mets);
met_bestsol_combined_CVR_filtered = filter_bestsols_by_index(met_bestsol_combined_CVR, select_mets);

% also in met_info
met_info_combined_DC_filtered = filter_metinfo_by_index(met_info_combined_DC, select_mets);
met_info_combined_CVR_filtered = filter_metinfo_by_index(met_info_combined_CVR, select_mets);

% add model names for plotting and saving to file
met_bestsol_combined_DC_filtered.modelname = 'bestsol_combined_DC_annmets';
met_bestsol_combined_CVR_filtered.modelname = 'bestsol_combined_CVR_annmets';

[confmatDCCVRannmets, classlabelsDCCVRannmets] = compare_two_model_results(...
                          met_info_combined_DC_filtered,...
                          met_bestsol_combined_DC_filtered,...
                          met_info_combined_CVR_filtered,...
                          met_bestsol_combined_CVR_filtered,...
                          figureFolder, 1,... %flag strictclass
                          0); %compare by mzrt and not compound name

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read meier modelling results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% combine IP with LI PCC within high total PCC solutions
inputfilename = [outputFolder ...
            'table_model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_Meier'];
[met_info_meier, met_bestsol_meier] = ...
                combine_bestsols_from_file(inputfilename, sel_crit1, sel_crit2,...
                                           corrthreshold);

% select metabolites that map to Meier et al metabolites in the annotation table
select_mets = cellfun(@(x) ~isempty(x), annotationTableMeier.Meier_compound);

met_info_intersect_meier=[];
met_info_intersect_meier.MZ =annotationTableMeier.MZ(select_mets);
met_info_intersect_meier.RT =annotationTableMeier.RT(select_mets);
met_info_intersect_meier.CompoundID = annotationTableMeier.CompoundID(select_mets);
met_info_intersect_meier.CompoundName = annotationTableMeier.Meier_compound(select_mets);
% remove spaces and capitalize first letters
for i=1:length(met_info_intersect_meier.CompoundName)
    str=lower(met_info_intersect_meier.CompoundName{i});
    idx=regexp([' ' str],'(?<=\s+)\S','start')-1;
    str(idx)=upper(str(idx));
    str = strrep(str, ' ','');
    str = strrep(str, '-','_');
    met_info_intersect_meier.CompoundName{i} = str;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
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
met_bestsol_meier.modelname = 'Meier_LIPCCwIP';

% leave only those Meier metabolites that are overlapping with this study
[~, ~, selectidx] = intersect(cellfun(@(x) lower(x), met_info_intersect_meier.CompoundName, 'unif',0),...
    cellfun(@(x) lower(x), met_info_meier.CompoundName, 'unif', 0),'stable');

met_bestsol_meier = filter_bestsols_by_index(met_bestsol_meier, selectidx);
met_info_meier = filter_metinfo_by_index(met_info_meier, selectidx);

% keep only selected metabolites
met_bestsol_combined_DC_filtered = filter_bestsols_by_index(met_bestsol_combined_DC, select_mets);
met_bestsol_combined_CVR_filtered = filter_bestsols_by_index(met_bestsol_combined_CVR, select_mets);

% also in met_info
met_info_combined_DC_filtered = filter_metinfo_by_index(met_info_combined_DC, select_mets);
met_info_combined_CVR_filtered = filter_metinfo_by_index(met_info_combined_CVR, select_mets);

% add model names for plotting and saving to file
met_bestsol_combined_DC_filtered.modelname = 'bestsol_combined_DC_meiermets';
met_bestsol_combined_CVR_filtered.modelname = 'bestsol_combined_CVR_meiermets';


[confmatDCMeier, classlabelsDCMeier] = compare_two_model_results(...
    met_info_intersect_meier, met_bestsol_combined_DC_filtered,...
    met_info_intersect_meier, met_bestsol_meier,...
    figureFolder,1,0);

[confmatCVRMeier, classlabelsCVRMeier] = compare_two_model_results(...
    met_info_intersect_meier, met_bestsol_combined_CVR_filtered,...
    met_info_intersect_meier, met_bestsol_meier,...
    figureFolder,1,0);

% [confmatDCMeier, classlabelsDCMeier] = compare_two_model_results(...
%     met_info_intersect_meier, met_bestsol_combined_DC_filtered,...
%     met_info_intersect_meier, met_bestsol_meier,...
%     figureFolder,0,0);
% 
% [confmatCVRMeier, classlabelsCVRMeier] = compare_two_model_results(...
%     met_info_intersect_meier, met_bestsol_combined_CVR_filtered,...
%     met_info_intersect_meier, met_bestsol_meier,...
%     figureFolder,0,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate agreement between diets as accuracy / F1 score
                      
classification_report_DC_CV_allions = create_classification_report(...
    confmatDCCVRallions{1,1}, classlabelsDCCVRallions);

classification_report_DC_CV_annmets = create_classification_report(...
    confmatDCCVRannmets{1,1}, classlabelsDCCVRallions);

classification_report_DC_Meier_HCD = create_classification_report(...
    confmatDCMeier{1,1}, classlabelsDCMeier);

classification_report_DC_Meier_HFD = create_classification_report(...
    confmatDCMeier{2,1}, classlabelsDCMeier);

classification_report_CV_Meier = create_classification_report(...
    confmatDCCVRannmets{1,1}, classlabelsDCCVRallions);
