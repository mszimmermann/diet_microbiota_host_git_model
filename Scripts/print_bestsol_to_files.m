function print_bestsol_to_files(met_info, met_bestsols, filename)
% met_info contains information about metabolites
% met_bestsols has best solutions
% filename is the output file name

% check whether met_info contains MZ and RT fiels
% if not (e.g. drug info) set them to 0
if isfield(met_info, 'MZ') == 0
    met_info.MZ = zeros(size(met_bestsols,1),1);
end
if isfield(met_info, 'RT') == 0
    met_info.RT = zeros(size(met_bestsols,1),1);
end

% check whether met_info contains MetaboliteFilter - flag indicating
% annotation filtering
if isfield(met_info, 'MetaboliteFilter') == 0
    met_info.MetaboliteFilter = ones(size(met_bestsols,1),1);
end

% check whether met_info contains gut_filter - flag indicating
% whether metabolites were detected in the GIT
if isfield(met_info, 'gut_filter') == 0
    met_info.gut_filter = ones(size(met_bestsols,1),1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save model results to file - raw coefficients

% save each type of selected solution to a separate file
for sel_crit = 1:length(met_bestsols{1}.selection_criterion)
    curfilename = [filename...
                   strrep(met_bestsols{1}.selection_criterion{sel_crit},' ','_'),...
                   '.csv'];
    fid = fopen(curfilename, 'w');
    fprintf(fid, ['MZ\tRT\tCompoundID\tCompoundName\tMetaboliteFilter\tSumGITclusters\t'...
                  'ReciprocalCorr\tReciprocalCorrSI\tReciprocalCorrLI\tReciprocalCorrMean\t']);
    for i=1:length(met_bestsols{1}.coefvalues)
        fprintf(fid, '\t%s', met_bestsols{1}.coefvalues{i});
    end
    fprintf(fid, '\n');
    for i=1:length(met_bestsols)
        fprintf(fid, '%.3f\t%3f',   met_info.MZ(i),...
                                    met_info.RT(i));
        fprintf(fid, '\t%s\t%s\t%d',met_info.CompoundID{i},...
                                    met_info.CompoundName{i},...
                                    met_info.MetaboliteFilter(i));
        fprintf(fid, '\t%d', met_info.gut_filter(i));
        
        if sel_crit <= length(met_bestsols{i}.x_sel_CorrRev)
            fprintf(fid, '\t%.3f', met_bestsols{i}.x_sel_CorrRev(sel_crit));
            fprintf(fid, '\t%.3f', met_bestsols{i}.x_sel_CorrRevSI(sel_crit));
            fprintf(fid, '\t%.3f', met_bestsols{i}.x_sel_CorrRevLI(sel_crit));
            fprintf(fid, '\t%.3f', mean([met_bestsols{i}.x_sel_CorrRev(sel_crit),...
                                         met_bestsols{i}.x_sel_CorrRevSI(sel_crit),...
                                         met_bestsols{i}.x_sel_CorrRevLI(sel_crit)]));
            for j=1:size(met_bestsols{i}.x,1)
                fprintf(fid, '\t%e', met_bestsols{i}.x(j,sel_crit));
            end
        else
            % best solution of this type was not possible to select for
            % this metabolite
            fprintf(fid, '\t0\t0\t0\t0');
            for j=1:size(met_bestsols{i}.x,1)
                fprintf(fid, '\t0');
            end
        end
            
        fprintf(fid, '\n');
    end
    fclose(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% save model results to file - reciprocal data restoration
% save each type of selected solution to a separate file
for sel_crit = 1:length(met_bestsols{1}.selection_criterion)
    curfilename = [filename...
                   '_reciprocal_'... 
                   strrep(met_bestsols{1}.selection_criterion{sel_crit},' ','_'),...
                   '.csv'];
    fid = fopen(curfilename, 'w');

    % column names for the data
    columnNames = met_bestsols{1}.kmeanMatrix_joint_names;
    
    fprintf(fid, ['MZ\tRT\tCompoundID\tCompoundName\tMetaboliteFilter\tSumGITclusters\t'...
                  'ReciprocalCorr\tReciprocalCorrSI\tReciprocalCorrLI\tReciprocalCorrMean\t']);
    for i=1:length(columnNames)
        fprintf(fid, '\t%s', columnNames{i});
    end
    for i=1:length(columnNames)
        fprintf(fid, '\tRecip_%s', columnNames{i});
    end
    fprintf(fid, '\n');
    for i=1:length(met_bestsols)
        fprintf(fid, '%.3f\t%3f',   met_info.MZ(i),...
                                    met_info.RT(i));
        fprintf(fid, '\t%s\t%s\t%d',met_info.CompoundID{i},...
                                    met_info.CompoundName{i},...
                                    met_info.MetaboliteFilter(i));
        fprintf(fid, '\t%d', met_info.gut_filter(i));
        
        if sel_crit <= length(met_bestsols{i}.x_sel_CorrRev)
        
            fprintf(fid, '\t%.3f', met_bestsols{i}.x_sel_CorrRev(sel_crit));
            fprintf(fid, '\t%.3f', met_bestsols{i}.x_sel_CorrRevSI(sel_crit));
            fprintf(fid, '\t%.3f', met_bestsols{i}.x_sel_CorrRevLI(sel_crit));
            fprintf(fid, '\t%.3f', mean([met_bestsols{i}.x_sel_CorrRev(sel_crit),...
                                         met_bestsols{i}.x_sel_CorrRevSI(sel_crit),...
                                         met_bestsols{i}.x_sel_CorrRevLI(sel_crit)]));
            kmean_vector_joint_orig = met_bestsols{i}.kmeanMatrix_joint_orig(:);
            for j=1:length(kmean_vector_joint_orig)
                fprintf(fid, '\t%.3f', kmean_vector_joint_orig(j));
            end
            for j=1:size(met_bestsols{i}.x_sel_dataR,1)
                fprintf(fid, '\t%.3f', met_bestsols{i}.x_sel_dataR(j));
            end
            fprintf(fid, '\n');
        else
            % best solution of this type was not possible to select for
            % this metabolite
            fprintf(fid, '\t0\t0\t0\t0');
            kmean_vector_joint_orig = met_bestsols{i}.kmeanMatrix_joint_orig(:);
            for j=1:length(kmean_vector_joint_orig)
                fprintf(fid, '\t%.3f', kmean_vector_joint_orig(j));
            end
            for j=1:length(kmean_vector_joint_orig)
                fprintf(fid, '\t0');
            end
        end
    end
    fclose(fid);
end