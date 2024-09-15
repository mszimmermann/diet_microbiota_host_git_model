function print_combined_bestsol_to_files(met_info, met_bestsol_combined, filename)
% met_info contains information about metabolites
% met_bestsol_combined has best combined solutions
% filename is the output file name

% check whether met_info contains MZ and RT fiels
% if not (e.g. drug info) set them to 0
if isfield(met_info, 'MZ') == 0
    met_info.MZ = zeros(size(met_bestsol_combined.x,1),1);
end
if isfield(met_info, 'RT') == 0
    met_info.RT = zeros(size(met_bestsol_combined.x,1),1);
end

% check whether met_info contains MetaboliteFilter - flag indicating
% annotation filtering
if isfield(met_info, 'MetaboliteFilter') == 0
    met_info.MetaboliteFilter = ones(size(met_bestsol_combined.x,1),1);
end

% check whether met_info contains gut_filter - flag indicating
% whether metabolites were detected in the GIT
if isfield(met_info, 'gut_filter') == 0
    met_info.gut_filter = ones(size(met_bestsol_combined.x,1),1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save model results to file - raw coefficients

% save selected solution to a separate file
curfilename = [filename, '.csv'];
fid = fopen(curfilename, 'w');
fprintf(fid, ['MZ\tRT\tCompoundID\tCompoundName\tMetaboliteFilter\tSumGITclusters\t'...
              'SolutionType\tReciprocalCorr\tReciprocalCorrSI\tReciprocalCorrLI\tReciprocalCorrMean']);
for i=1:length(met_bestsol_combined.coefvalues)
    fprintf(fid, '\t%s', met_bestsol_combined.coefvalues{i});
end
fprintf(fid, '\n');
for i=1:length(met_bestsol_combined.selection_criteria)
    fprintf(fid, '%.3f\t%3f',   met_info.MZ(i),...
                                met_info.RT(i));
    fprintf(fid, '\t%s\t%s\t%d',met_info.CompoundID{i},...
                                met_info.CompoundName{i},...
                                met_info.MetaboliteFilter(i));
    fprintf(fid, '\t%d', met_info.gut_filter(i));
    
    %print solution type
    fprintf(fid, '\t%s', met_bestsol_combined.selection_criteria{i});
    
    fprintf(fid, '\t%.3f', met_bestsol_combined.x_sel_CorrRev(i));
    fprintf(fid, '\t%.3f', met_bestsol_combined.x_sel_CorrRevSI(i));
    fprintf(fid, '\t%.3f', met_bestsol_combined.x_sel_CorrRevLI(i));
    fprintf(fid, '\t%.3f', met_bestsol_combined.x_sel_CorrRevMean(i));
    for j=1:size(met_bestsol_combined.x,2)
        fprintf(fid, '\t%e', met_bestsol_combined.x(i,j));
    end
    
    fprintf(fid, '\n');
end
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% save model results to file - reciprocal data restoration
% save each type of selected solution to a separate file
curfilename = [filename...
               '_reciprocal.csv'];
fid = fopen(curfilename, 'w');

% column names for the data
columnNames = met_bestsol_combined.kmeanMatrix_joint_names(:);

fprintf(fid, ['MZ\tRT\tCompoundID\tCompoundName\tMetaboliteFilter\tSumGITclusters\t'...
              'SolutionType\tReciprocalCorr\tReciprocalCorrSI\tReciprocalCorrLI\tReciprocalCorrMean']);
for i=1:length(columnNames)
    fprintf(fid, '\t%s', columnNames{i});
    
end
for i=1:length(columnNames)
    fprintf(fid, '\tRecip_%s', columnNames{i});
end
fprintf(fid, '\n');
for i=1:length(met_bestsol_combined.selection_criteria)
    fprintf(fid, '%.3f\t%3f',   met_info.MZ(i),...
                                met_info.RT(i));
    fprintf(fid, '\t%s\t%s\t%d',met_info.CompoundID{i},...
                                met_info.CompoundName{i},...
                                met_info.MetaboliteFilter(i));
    fprintf(fid, '\t%d', met_info.gut_filter(i));
    
    fprintf(fid, '\t%s', met_bestsol_combined.selection_criteria{i});

    fprintf(fid, '\t%.3f', met_bestsol_combined.x_sel_CorrRev(i));
    fprintf(fid, '\t%.3f', met_bestsol_combined.x_sel_CorrRevSI(i));
    fprintf(fid, '\t%.3f', met_bestsol_combined.x_sel_CorrRevLI(i));
    fprintf(fid, '\t%.3f', met_bestsol_combined.x_sel_CorrRevMean(i));
            
    for j=1:length(met_bestsol_combined.kmeanMatrix_joint_names)
        fprintf(fid, '\t%.3f', met_bestsol_combined.kmeanMatrix_joint_orig(i,j));
    end
    for j=1:length(met_bestsol_combined.kmeanMatrix_joint_names)
        fprintf(fid, '\t%.3f', met_bestsol_combined.x_sel_dataR(i,j));
    end
        
    fprintf(fid, '\n');
end
fclose(fid);
