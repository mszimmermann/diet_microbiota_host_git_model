function [met_info, met_bestsol] = read_bestsol_from_file(filename, sel_crit)
% sel_crit can take the following arguments:
% 'IP' 'total_PCC' 'SI_PCC' 'LI_PCC' 'sum_PCC'
% 'SI_PCC_within_high_total' 'LI_PCC_within_high_total'

% read specific solution type
curfilename = [filename...
               strrep(sel_crit,' ','_'),...
               '.csv'];
soldata = readtable(curfilename,...
    'ReadVariableNames',1, 'HeaderLines',0);
% get metabolite info
met_info.MZ = soldata.MZ;
met_info.RT = soldata.RT;
met_info.CompoundID = soldata.CompoundID;
met_info.CompoundName = soldata.CompoundName;
met_info.MetaboliteFilter = soldata.MetaboliteFilter;
met_info.SumGITclusters = soldata.SumGITclusters;

% get solution correlation info and solution as matrix
columnNames = soldata.Properties.VariableNames;
columnNames_reciprocal_idx = find(cellfun(@(x) contains(lower(x), ...
                    'reciprocal'), columnNames));
% save the filed names in the same way as in bestsol
met_bestsol.x_sel_CorrRev = soldata.ReciprocalCorr;
met_bestsol.x_sel_CorrRevSI = soldata.ReciprocalCorrSI;
met_bestsol.x_sel_CorrRevLI = soldata.ReciprocalCorrLI;
met_bestsol.x_sel_CorrRevMean = soldata.ReciprocalCorrMean;
% get the x values that come after correlation columns
met_bestsol.x = soldata{:, columnNames_reciprocal_idx(end)+1:end};
met_bestsol.coefvalues = columnNames(columnNames_reciprocal_idx(end)+1:end);
% save criterion and file
met_bestsol.selection_criterion = sel_crit;
met_bestsol.fileName = curfilename;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read reciprocal data
curfilenameR = [filename...
               '_reciprocal_'... 
               strrep(sel_crit,' ','_'),...
               '.csv'];
soldataR = readtable(curfilenameR,...
    'delim', '\t', 'ReadVariableNames',1, 'HeaderLines',0);
if length(soldataR.Properties.VariableNames)==1
    soldataR = readtable(curfilenameR,...
    'delim', ',', 'ReadVariableNames',1, 'HeaderLines',0);
end
% check that the metabolites are the same in soldata and soldataR
if ~isequal(soldata.CompoundID, soldataR.CompoundID)
    fprintf(['Warning: metabolite IDs do not match between:\n%s\nand\n%s\n'...
             'omit loading reciprocal data\n'], curfilename, curfilenameR)
else
    % load original data
    % get colunms with original data
    columnNamesR = soldataR.Properties.VariableNames;
    columnNamesR_reciprocal_idx = find(cellfun(@(x) contains(lower(x), ...
                    'reciprocal'), columnNamesR));
    columnNamesR_reciprocal_data_idx = find(cellfun(@(x) contains(lower(x), ...
                    'recip_'), columnNamesR));
                
    met_bestsol.kmeanMatrix_joint_names = columnNamesR(columnNamesR_reciprocal_idx(end)+1:...
                        columnNamesR_reciprocal_data_idx(1)-1);
    met_bestsol.kmeanMatrix_joint_orig = soldataR{:,columnNamesR_reciprocal_idx(end)+1:...
                        columnNamesR_reciprocal_data_idx(1)-1};

    % check if the correlations are the same in the reciprocal file
    % if not, this might be a different solution
    if max(abs(soldata.ReciprocalCorr-soldataR.ReciprocalCorr))>eps
        fprintf(['Warning: reciprocal correlations do not match between:\n%s\nand\n%s\n'...
          'omit loading reciprocal data\n'], curfilename, curfilenameR)
    else
        % add reciprocal data to met_bestsol
        met_bestsol.x_sel_dataR = soldataR{:,columnNamesR_reciprocal_data_idx};
    end
end  
