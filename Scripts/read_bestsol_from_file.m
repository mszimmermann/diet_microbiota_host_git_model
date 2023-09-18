function [met_info, met_bestsol] = read_bestsol_from_file(filename, sel_crit)

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
for i=1:length(columnNames_reciprocal_idx)
    met_bestsol.(columnNames{columnNames_reciprocal_idx(i)}) = ...
        soldata.(columnNames{columnNames_reciprocal_idx(i)});
end
met_bestsol.sol_coefs = soldata{:, columnNames_reciprocal_idx(end)+1:end};
met_bestsol.coefvalues = columnNames(columnNames_reciprocal_idx(end)+1:end);
% save criterion and file
met_bestsol.sel_crit = sel_crit;
met_bestsol.fileName = curfilename;

    
