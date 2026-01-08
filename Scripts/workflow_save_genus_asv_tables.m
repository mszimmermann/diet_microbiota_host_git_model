%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save OTU/ASV tables from 16S sequencing to supplementary tables

% call script defining file dependenciesand global variables
addpath(genpath('.\'));
add_global_and_file_dependencies

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Required files:
% 'genus_RA_table.txt'
% 'mouse_metadata.txt'
% '221216_miSeq_MZK_Goodman.txt'
% Output: 
% Files:
% table_S3_genus_RA_feces_over_time.csv
% table_S3_genus_RA_across_tissues_final.csv

genus_RA_table = readtable([inputFolderSeq,...
    '\DADA2\genus_RA_table.txt']);

metadata_table = readtable([inputFolderSeq,...
    '\DADA2\mouse_metadata.txt'], 'ReadVariableNames',false);

sample_table = readtable([inputFolderSeq,...
    '\DADA2\221216_miSeq_MZK_Goodman.txt']);

% create sample names according to mouse and time or tissue
sample_name_merged = sample_table.Animal;
for i=1:length(sample_name_merged)
    sample_name_merged{i} = strcat(sample_name_merged{i}, '_', sample_table.Timepoint_Tissue{i});
end
% add merged sample name to the table
sample_table = addvars(sample_table, sample_name_merged);

% change sample names in genus table to sample_name_merged
new_sample_names = genus_RA_table.Properties.VariableNames;
% add underscore to match the sample table format
new_sample_names = cellfun(@(x) strrep(x, 'MG', 'MG_'), new_sample_names, 'unif', 0);
for i=1:length(new_sample_names)
    if sum(ismember(sample_table.SampleID, new_sample_names(i)))>0
        new_sample_names(i) = sample_table.sample_name_merged(ismember(sample_table.SampleID, new_sample_names{i}));
    end
end
% check if sample corresponds to a timepoint
sample_type_time = cellfun(@(x) contains(x, '_T'), new_sample_names);

% separate genus RA tables into time course and tissue tables
% time course table
genus_RA_table_timecourse = genus_RA_table(:,sample_type_time==1);
genus_RA_table_timecourse.Properties.VariableNames = new_sample_names(sample_type_time==1);
% sort columns by alphabet
[~, sortidx] = sort(genus_RA_table_timecourse.Properties.VariableNames);
genus_RA_table_timecourse = genus_RA_table_timecourse(:, sortidx);
genus_RA_table_timecourse = addvars(genus_RA_table_timecourse,...
    genus_RA_table.Var1, 'NewVariableNames', 'Genus',...
    'before', genus_RA_table_timecourse.Properties.VariableNames{1});

writetable(genus_RA_table_timecourse ,...
    [outputFolder,...
    'table_S3_genus_RA_feces_over_time.csv']);

% tissue table
genus_RA_table_tissues = genus_RA_table(:,sample_type_time==0);
genus_RA_table_tissues.Properties.VariableNames = new_sample_names(sample_type_time==0);
genus_RA_table_tissues.Properties.VariableNames{1} = 'Genus';
% sort columns by alphabet
[~, sortidx] = sort(genus_RA_table_tissues.Properties.VariableNames);
genus_RA_table_tissues = genus_RA_table_tissues(:, sortidx);

writetable(genus_RA_table_tissues ,...
    [outputFolder,...
    'table_S3_genus_RA_across_tissues_final.csv']);