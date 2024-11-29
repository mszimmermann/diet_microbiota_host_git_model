function merged_table = combine_shortbred_results(folderName, filePathName)

% combine shortbred results from multiple folders and save to one file if
% filePathName provided

shortbredFolders = dir(folderName);
shortbredFolderNames = cell(size(shortbredFolders,1),1);
for i=1:length(shortbredFolders)
    shortbredFolderNames{i} = shortbredFolders(i).name;
end
%only keep folder names that have _results in the name
shortbredFolderNames = shortbredFolderNames(...
                                    cellfun(@(x) contains(x,'_results'),...
                                                 shortbredFolderNames))';
% shortbred file names are inside the folders with shortBRED_results ending
shortbredFileNames = cellfun(@(x) strrep(x, '_results', '_shortBRED_results.txt'),...
                                    shortbredFolderNames, 'unif', 0);

id_column = 'Family';
select_column = 'Count';
% read and join all results files
merged_table = [];
for i=1:length(shortbredFileNames)
    curtable = readtable([folderName, filesep, shortbredFolderNames{i},...
        filesep, shortbredFileNames{i}]);
    % select and rename count column
    curtable = curtable(:,{id_column, select_column});
    curtable.Properties.VariableNames = cellfun(@(x) strrep(x, select_column,...
        [select_column, '_', shortbredFolderNames{i}]),...
        curtable.Properties.VariableNames, 'unif', 0);
    if isempty(merged_table)
        merged_table = curtable;
    else
        merged_table = outerjoin(merged_table, curtable, 'MergeKeys', 1);
    end
end

% save to file if filepathname is provided
if nargin>1
    writetable(merged_table, filePathName);
end

    