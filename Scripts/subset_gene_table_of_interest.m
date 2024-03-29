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

met_gene_table = readtable([outputFolder, ...
    'table_DNA_kegg_RNA_sub_prod_products_bestcorr_genes_pos.csv'], 'delim', ',');
% get gene annotation
geneAnnotationTable = readtable([inputFolderSeq, ...
    'countsMatrix_annTable_filtered.csv'], 'delim', ',');
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
    ismember(geneAnnotationTable.species_gene, gene_list),:);
% add genes for biotin enzyme
selected_ec = {'2.8.1.6', '2.4.2.1', '2.4.2.2', '2.4.2.3', '2.4.2.4'};
selected_genes = ecTable(ismember(ecTable.PathwayID,selected_ec),:);
selected_genes_idx = cat(2,selected_genes.GeneIDXfiltered{:});
selected_genes_idx = strsplit(selected_genes_idx,';');
selected_genes_idx(cellfun(@(x) isempty(x), selected_genes_idx))=[];
selected_genes_idx = cellfun(@(x) str2num(x), selected_genes_idx);

selected_genes_ann = geneAnnotationTable(selected_genes_idx,:);
% remove rows without EC
selected_genes_ann(ismember(selected_genes_ann.EC, {'nan'}),:)=[];

gene_list_ann = [gene_list_ann; selected_genes_ann];

writetable(gene_list_ann, [resultsFolder, 'selected_genes_ann.csv']);