%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract genes that correlate with metabolites and save them to a separate file

% call script defining file dependenciesand global variables
addpath(genpath('.\'));
add_global_and_file_dependencies

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Required files:
% Output: 
% Files:
% Figures:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get gene ID table
met_gene_table = readtable([outputFolder, ...
    'table_DNA_kegg_RNA_sub_prod_products_bestcorr_genes_pos.csv'], 'delim', ',');
% get EC table
met_ec_table = readtable([outputFolder, ...
    'table_DNA_kegg_RNA_sub_prod_products_bestcorr_EC_pos.csv'], 'delim', ',');
% get gene annotation
% detect import options to change some of the variable to strings
%annfileName = [inputFolderSeq, 'countsMatrix_annTable_filtered.csv'];
annfileName =  [outputFolder 'geneAnnTable_filtered.csv'];
opts = detectImportOptions(annfileName);
if (ismember(opts.VariableNames, 'genome_id'))
    opts.VariableTypes{ismember(opts.VariableNames, 'genome_id')} = 'char';
end
for i = find(ismember(opts.VariableNames, 'seed_eggNOG_ortholog')):length(opts.VariableTypes)
    if (~isequal(opts.VariableNames{i}, 'seed_ortholog_evalue')) &&...
       (~isequal(opts.VariableNames{i}, 'seed_ortholog_score'))
        opts.VariableTypes{i} = 'char';  % This will change column 15 from 'double' or whatever it is to 'char', which is what you want.
    end
end
opts.Delimiter = ',';
geneAnnotationTable = readtable(annfileName, opts);

% get EC annotation
ecTable = readtable([inputFolderSeq, ...
    'gene_annotation_EC.tsv'], 'fileType','text', 'delim', '\t');

% get all species genes from met_gene_table
gene_list = met_gene_table{:,2:end};
gene_list = unique(gene_list(:));
gene_list(cellfun(@(x) isempty(x), gene_list)) = [];
% for genes with multiple names split
gene_list = cellfun(@(x) strsplit(x,'|'), gene_list, 'unif',0);
% concatenate into one gene list
gene_list = gene_list';
gene_list = unique(cat(2,gene_list{:}))';

gene_list_ann = geneAnnotationTable(...
    ismember(geneAnnotationTable.species_genes, gene_list),:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get all species genes from met_EC_table
ec_list = met_ec_table{:,2:end};
ec_list = unique(ec_list(:));
ec_list(cellfun(@(x) isempty(x), ec_list)) = [];
% for EC with multiple names split
ec_list = cellfun(@(x) strsplit(x,','), ec_list, 'unif',0);
% concatenate into one EC list
ec_list = ec_list';
ec_list = unique(cat(2,ec_list{:}))';

selected_genes = ecTable(ismember(ecTable.PathwayID,ec_list),:);
selected_genes_idx = cat(2,selected_genes.GeneIDXfiltered{:});
selected_genes_idx = strsplit(selected_genes_idx,';');
selected_genes_idx(cellfun(@(x) isempty(x), selected_genes_idx))=[];
selected_genes_idx = cellfun(@(x) str2num(x), selected_genes_idx);
selected_genes_idx = unique(selected_genes_idx');

selected_ec_ann = geneAnnotationTable(selected_genes_idx,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add genes for biotin enzyme
selected_ec = {'2.8.1.6', '2.4.2.1', '2.4.2.2', '2.4.2.3', '2.4.2.4'};
selected_genes = ecTable(ismember(ecTable.PathwayID,selected_ec),:);
selected_genes_idx = cat(2,selected_genes.GeneIDXfiltered{:});
selected_genes_idx = strsplit(selected_genes_idx,';');
selected_genes_idx(cellfun(@(x) isempty(x), selected_genes_idx))=[];
selected_genes_idx = cellfun(@(x) str2num(x), selected_genes_idx);

selected_genes_ann = geneAnnotationTable(selected_genes_idx,:);

%gene_list_ann = [gene_list_ann; selected_genes_ann];
gene_list_ann = [selected_ec_ann; selected_genes_ann];
 
%writetable(gene_list_ann, [resultsFolder, 'selected_genes_gitann.csv']);
writetable(gene_list_ann, [resultsFolder, 'selected_genes_from_ec_ann.csv']);