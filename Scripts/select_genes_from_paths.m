%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% match path and gene expression data
% test which of the single enzymes are expressed

% call script defining file dependenciesand global variables
addpath(genpath('.\'));
add_global_and_file_dependencies

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shortestPathTable = readtable([resultsFolder, ...
    'shortest_paths_upto_length_4.csv'], 'delim', ',');
% get gene FC
geneTable = readtable([inputFolderSeq, ...
    'edgeR_gene_fold_changes_and_ann.csv'], 'delim', ',');
% get EC annotation
ecTable = readtable([inputFolderSeq, ...
    'gene_annotation_EC.tsv'], 'fileType','text', 'delim', '\t');

% get list of genes of interest for checking in CVR expression
% select all genes for which there is a path of length 1
shortestPathTable_length = zeros(size(shortestPathTable,1),1);
for i=1:size(shortestPathTable,1)
    shortestPathTable_length(i) = length(strsplit(shortestPathTable.EC_path{i},';'));
end

shortestPathTable_length1 = shortestPathTable(shortestPathTable_length==1,:);
target_ec_ids = shortestPathTable_length1.EC_path;
target_ec_ids = cellfun(@(x) strsplit(x, '|')', target_ec_ids, 'unif', 0);
target_ec_ids = cat(1,target_ec_ids{:});
target_ec_ids = unique(target_ec_ids);

target_gene_ids = zeros(size(geneTable,1),1);
for i=1:length(target_ec_ids)
    ec_lookup = ecTable(ismember(ecTable.PathwayID, target_ec_ids{i}),:);
    if ~isempty(ec_lookup)
        gene_ids = strsplit(ec_lookup.GeneIDXfiltered{1},';');
        gene_ids(cellfun(@(x) isempty(x), gene_ids))=[];
        gene_ids = cellfun(@(x) str2num(x), gene_ids); 
        
        target_gene_ids(gene_ids)=1;
    end
end

target_gene_table = geneTable(target_gene_ids==1,:);    

% write target table to file
writetable(target_gene_table, 'ProcessedData\target_gene_table_path1.csv');
    