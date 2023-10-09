function [met_info_combined, met_bestsol_combined] = ...
                combine_bestsols_from_file(filename, sel_crit1, sel_crit2)
% corr threshold for reliable solutions
corrthreshold = 0.7;
% combine solutions from criteria sel_crit1 and sel_crit2
[met_info1, met_bestsol1] = read_bestsol_from_file(filename,...
                                                sel_crit1);
                                            
[met_info2, met_bestsol2] = read_bestsol_from_file(filename,...
                                                sel_crit2);

met_info_combined = met_info1;
met_bestsol_combined = met_bestsol1;
met_bestsol_combined_criterion = cell(length(met_bestsol1.x_sel_CorrRev),1);

% check that met_info structures are the same
if size(met_info1.MZ)== size(met_info2.MZ)
    if (met_info1.MZ ~= met_info2.MZ)
        fprintf(['Warning: metabolite lists do not match between files %s %s\n%s\n'...
                 'Loading solution from file %s\n'],...
                 filename, sel_crit1, sel_crit2, sel_crit1);
         return
    end
else
    fprintf(['Warning: metabolite lists do not match between files %s %s\n%s\n'...
                 'Loading solution from file %s\n'],...
                 filename, sel_crit1, sel_crit2, sel_crit1);
    return
end    
for i=1:length(met_bestsol1.x_sel_CorrRev)
    % if both types of best solutions exist, take the first type as
    % long as it passess corr thresholds, otherwise take the second
    % type
    if (met_bestsol1.x_sel_CorrRev(i)>=corrthreshold) &&...
       (met_bestsol1.x_sel_CorrRevLI(i)>=corrthreshold)
        %keep first solution as it passes the threshold
        met_bestsol_combined_criterion{i} = sel_crit1; 
    else
        % second solution type passess both thresholds - take
        % this one
        if (met_bestsol2.x_sel_CorrRev(i)>=corrthreshold) &&...
           (met_bestsol2.x_sel_CorrRevLI(i)>=corrthreshold)
            met_bestsol_combined_criterion{i} = sel_crit2; 
        else
            % first solution passes total corr threshold but
            % not LI corr threshold
            if (met_bestsol1.x_sel_CorrRev(i)>=corrthreshold) 
                met_bestsol_combined_criterion{i} = sel_crit1; 
            else
                % second solution passes the threshold
                if (met_bestsol2.x_sel_CorrRev(i)>=corrthreshold) 
                    met_bestsol_combined_criterion{i} = sel_crit2;
                else
                    % second solution is not better - keep the
                    % first one
                    met_bestsol_combined_criterion{i} = sel_crit1;
                end
            end
        end
    end

    % check if first solution type shuld be replaced and replace if yes
    if isequal(met_bestsol_combined_criterion{i}, sel_crit2)
        % replace all items with the ones from the second silution type
        met_bestsol_combined.x_sel_CorrRev(i) = met_bestsol2.x_sel_CorrRev(i);
        met_bestsol_combined.x_sel_CorrRevSI(i) = met_bestsol2.x_sel_CorrRevSI(i);
        met_bestsol_combined.x_sel_CorrRevLI(i) = met_bestsol2.x_sel_CorrRevLI(i);
        met_bestsol_combined.x_sel_CorrRevMean(i) = met_bestsol2.x_sel_CorrRevMean(i);
        met_bestsol_combined.x(i,:) = met_bestsol2.x(i,:);
        met_bestsol_combined.x_sel_dataR(i,:) = met_bestsol2.x_sel_dataR(i,:);
    end
end

met_bestsol_combined.selection_criterion = {sel_crit1, sel_crit2};
% add information about criterion for each metabolite
met_bestsol_combined.selection_criteria = met_bestsol_combined_criterion;
