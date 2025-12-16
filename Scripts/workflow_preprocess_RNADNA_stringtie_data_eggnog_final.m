%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preprocess raw sequencing reads, edgeR and GetMM results and 
% format them into fold change and read tables
% extract kegg, COG, eggNOG, GO terms from annotation and create 
% gene-category binary tables

% call script defining file dependenciesand global variables
addpath(genpath('.\'));
add_global_and_file_dependencies

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Required files:
% Raw counts, edgeR and GetMM-processed files in the folders
% '.\InputData\sequencing_data\ballGown_RNA\eggNOGann\'
% and '.\InputData\sequencing_data\ballGown_DNA\eggNOGann\'
% .\InputData\genomes14_eggnog\genome_list.csv contains list of genome names

% Output: 
% Files:
% 'countsMatrixGetMM_DNA_table.txt'
% 'countsMatrixGetMM_RNA_table.txt'
% 'countsMatrixRAW_DNA_table.txt'
% 'countsMatrixRAW_RNA_table.txt'
% 'geneAnnTable_filtered.csv'
% 'geneAnnTable_full.csv'
% 'edgeR_gene_fold_changes_and_ann.csv'
% 'gene_annotation_kegg.tsv'
% 'gene_annotation_go.tsv'
% 'gene_annotation_EC.tsv'
% 'gene_annotation_EggNOG.tsv'
% 'gene_annotation_COG.tsv'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RNA results

% read genome abbreviations and names
genome_abbr_names = readtable(['.\InputData\genomes14_eggnog\'...
                                'genome_list.csv'], 'delim', ',');
% combine name and strain
genome_abbr_names.("genome_name") = strcat(genome_abbr_names.Species_name,...
                                    {' '}, genome_abbr_names.Strain);

% read eggnog annotations                   
annAbbr = 'EGGNOGann';

stringtieFiles = dir(stringtieFolderRNA);
stringtieFileNames = cell(size(stringtieFiles ,1),1);
for i=1:length(stringtieFiles)
    stringtieFileNames{i} = stringtieFiles(i).name;
end
%extract file names of raw stringtie files
stringtieFileNamesRNA_raw = stringtieFileNames(...
                                    cellfun(@(x) contains(x,'stringtie') & ... 
                                                 contains(x,annAbbr) & ...
                                                 contains(x,'ucount') & ...
                                                 contains(x,'RNA'),...
                                                           stringtieFileNames))';
% get the merged raw reads file name
stringtieFileNamesRNA_raw = stringtieFileNamesRNA_raw(...
                                    cellfun(@(x) contains(x(1:length('stringtie')),'stringtie') &...
                                                 contains(x,'RNA100_t_e'),...
                                    stringtieFileNamesRNA_raw));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read differential expression analysis results
stringtieFilesDEG = dir(degFolderRNA);
stringtieFileNamesDEG = cell(size(stringtieFilesDEG,1),1);
for i=1:length(stringtieFilesDEG)
    stringtieFileNamesDEG{i} = stringtieFilesDEG(i).name;
end

%extract file names of edgeR results
stringtieFileNamesRNA_edgeR2results = stringtieFileNamesDEG(...
                                    cellfun(@(x) contains(x,'edgeR_LRT_results') & ... 
                                                 contains(x,annAbbr) & ...
                                                 contains(x,'ucount') & ...
                                                 contains(x,'RNA'),...
                                                           stringtieFileNamesDEG))';
%extract file names of deseq2 results
stringtieFileNamesRNA_deseq2results = stringtieFileNamesDEG(...
                                    cellfun(@(x) contains(x,'deseq_results') & ... 
                                                 contains(x,annAbbr) & ...
                                                 contains(x,'ucount') & ...
                                                 contains(x,'RNA'),...
                                                           stringtieFileNamesDEG))';

%extract file names of GetMM normalized counts
% get only the merged getmm file
getmmFileNamesRNA_norm = stringtieFileNamesDEG(cellfun(@(x) contains(x,'getmm_normalized') & ...
                                                           contains(x,annAbbr) & ...
                                                           contains(x,'ucount') & ...
                                                           contains(x,'RNA100_t_e'),...
                                                           stringtieFileNamesDEG))';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get DNA results
stringtieFiles = dir(stringtieFolderDNA);
stringtieFileNames = cell(size(stringtieFolderDNA,1),1);
for i=1:length(stringtieFiles)
    stringtieFileNames{i} = stringtieFiles(i).name;
end                                                       
%extract file names of raw stringtie files
stringtieFileNamesDNA_raw = stringtieFileNames(...
                                    cellfun(@(x) contains(x,'stringtie') & ... 
                                                 contains(x,annAbbr) & ...
                                                 contains(x,'ucount') & ...
                                                 contains(x,'DNA'),...
                                                           stringtieFileNames))';
% keep only the merged DNA file
stringtieFileNamesDNA_raw = stringtieFileNamesDNA_raw(...
                                    cellfun(@(x) contains(x(1:length('stringtie')),'stringtie') &...
                                                 contains(x,'DNA_t_e'),...
                                    stringtieFileNamesDNA_raw));
%extract file names of edgeR results
stringtieFileNamesDNA_edgeR2results = stringtieFileNames(...
                                    cellfun(@(x) contains(x,'edgeR_LRT_results') & ... 
                                                 contains(x,annAbbr) & ...
                                                 contains(x,'ucount') & ...
                                                 contains(x,'DNA'),...
                                                           stringtieFileNames))';
%extract file names of deseq2 results
stringtieFileNamesDNA_deseq2results = stringtieFileNames(...
                                    cellfun(@(x) contains(x,'deseq_results') & ... 
                                                 contains(x,annAbbr) & ...
                                                 contains(x,'ucount') & ...
                                                 contains(x,'DNA'),...
                                                           stringtieFileNames))';

%extract file names of GetMM normalized counts
getmmFileNamesDNA_norm = stringtieFileNames(cellfun(@(x) contains(x,'getmm_normalized') & ...
                                                           contains(x,annAbbr) & ...
                                                           contains(x,'ucount') & ...
                                                           contains(x,'DNA_t_e'),...
                                                           stringtieFileNames))';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                     
% extract species abbreviations from DNA file names
% (is there is no DNA coverage, ignore those species)
edgeRFileNames_abbr = cell(size(stringtieFileNamesDNA_edgeR2results));
for i=1:length(stringtieFileNamesDNA_edgeR2results)
    curabbr = strsplit(stringtieFileNamesDNA_edgeR2results{i}, '_');
    edgeRFileNames_abbr{i} = curabbr{12};
end    
edgeRFileNames_abbr(ismember(edgeRFileNames_abbr, {'t'})) = [];
% remove Bvulgatus
edgeRFileNames_abbr(ismember(edgeRFileNames_abbr, 'Bvul'))=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load differential analysis matrices for edgeR
edgeResultsDNA = cell(length(edgeRFileNames_abbr),1);
edgeResultsRNA = cell(length(edgeRFileNames_abbr),1);

% detect import options to import columns with the same type
opts = detectImportOptions([stringtieFolderDNA...
                            stringtieFileNamesDNA_edgeR2results{1}]);
% manual curation of variable types
opts.VariableNamesLine = 1;
opts.VariableTypes(16:end) = {'char'};
opts.VariableTypes([30 31]) = {'double'};

for method_i = 1:2
    for i=1:length(edgeRFileNames_abbr)
        % read edgeR results
        switch method_i
            case 1
                curfile = cellfun(@(x) contains(x, edgeRFileNames_abbr{i}),...
                                                stringtieFileNamesDNA_edgeR2results);
                curtable = readtable([stringtieFolderDNA...
                                      stringtieFileNamesDNA_edgeR2results{curfile}],...
                                      opts,...
                                      'ReadVariableNames', 1);
            case 2
                curfile = cellfun(@(x) contains(x, edgeRFileNames_abbr{i}),...
                                                stringtieFileNamesRNA_edgeR2results);
                curtable = readtable([degFolderRNA...
                                      stringtieFileNamesRNA_edgeR2results{curfile}],...
                                      opts,...
                                      'ReadVariableNames', 1);
        end
        % shift column names one to the right
        curtable.Properties.VariableNames = ['species_genes',...
                        curtable.Properties.VariableNames(1:end-1)];
        switch method_i
            case 1
                edgeResultsDNA{i} = curtable;
            case 2
                edgeResultsRNA{i} = curtable;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load differential analysis matrices for deseq2
deseqResultsDNA = cell(length(edgeRFileNames_abbr),1);
deseqResultsRNA = cell(length(edgeRFileNames_abbr),1);

% detect import options to import columns with the same type
opts = detectImportOptions([stringtieFolderDNA...
                            stringtieFileNamesDNA_deseq2results{1}]);
% manual curation of variable types
opts.VariableNamesLine = 1;
opts.VariableTypes(16:end) = {'char'};
opts.VariableTypes([30 31]) = {'double'};

for method_i = 1:2
    for i=1:length(edgeRFileNames_abbr)
        % read edgeR results
        switch method_i
            case 1
                curfile = cellfun(@(x) contains(x, edgeRFileNames_abbr{i}),...
                                                stringtieFileNamesDNA_deseq2results);
                curtable = readtable([stringtieFolderDNA...
                                      stringtieFileNamesDNA_deseq2results{curfile}],...
                                      opts,...
                                      'ReadVariableNames', 1);
            case 2
                curfile = cellfun(@(x) contains(x, edgeRFileNames_abbr{i}),...
                                                stringtieFileNamesRNA_deseq2results);
                curtable = readtable([degFolderRNA...
                                      stringtieFileNamesRNA_deseq2results{curfile}],...
                                      opts,...
                                      'ReadVariableNames', 1);
        end
        % shift column names one to the right
        curtable.Properties.VariableNames = ['species_genes',...
                        curtable.Properties.VariableNames(1:end-1)];
        switch method_i
            case 1
                deseqResultsDNA{i} = curtable;
            case 2
                deseqResultsRNA{i} = curtable;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get full getMM normalization file and separate by species 
% detect import options to import columns with the same type
opts = detectImportOptions([stringtieFolderDNA...
                          getmmFileNamesDNA_norm{1}],...
                          'Delimiter', '\t');
% manual curation of variable types
opts.VariableNamesLine=1;
opts.DataLine = 2;
opts.VariableTypes(2:11) = repmat({'double'},1,10);
opts.VariableTypes(15:end) = {'char'};
opts.VariableTypes([17, 18, 20, 35, 36]) = {'double'};
% read DNA getMM
curtable = readtable([stringtieFolderDNA...
                          getmmFileNamesDNA_norm{1}],...
                          opts,...
                          'ReadVariableNames', 1);
% shift column names one to the right
curtable.Properties.VariableNames = ['species_genes',...
                curtable.Properties.VariableNames(1:end-1)];
% get annotation from single tables
for i=1:length(edgeResultsRNA)
    curcols = intersect(edgeResultsRNA{i}.Properties.VariableNames, curtable.Properties.VariableNames);
    [~, ~, currows] = intersect(edgeResultsRNA{i}{:,1}, curtable{:,1}, 'stable');
    curtable(currows, curcols) = edgeResultsRNA{i}(:,curcols);
end
getmmNormDNA_full = curtable;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read RNA getMM
curtable = readtable([degFolderRNA...
                          getmmFileNamesRNA_norm{1}],...
                          opts,...
                          'ReadVariableNames', 1);
% shift column names one to the right
curtable.Properties.VariableNames = ['species_genes',...
                curtable.Properties.VariableNames(1:end-1)];
% get annotation from single tables
for i=1:length(edgeResultsRNA)
    curcols = intersect(edgeResultsRNA{i}.Properties.VariableNames, curtable.Properties.VariableNames);
    [~, ~, currows] = intersect(edgeResultsRNA{i}{:,1}, curtable{:,1}, 'stable');
    curtable(currows, curcols) = edgeResultsRNA{i}(:,curcols);
end
getmmNormRNA_full = curtable;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read DNA stringtie raw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get full raw stringtie file and separate by species 
% detect import options to import columns with the same type
opts = detectImportOptions([stringtieFolderDNA...
                           stringtieFileNamesDNA_raw{1}],...
                          'Delimiter', ',');
% manual curation of variable types
opts.VariableNamesLine=1;
opts.DataLine = 2;
opts.VariableTypes(1:3) = {'char'};
opts.VariableTypes(11:20) = repmat({'double'},1,10);
opts.VariableTypes(21:end) = {'char'};
opts.VariableTypes([34,35]) = {'double'};

curtable = readtable([stringtieFolderDNA...
                          stringtieFileNamesDNA_raw{1}],...
                          opts,...
                          'ReadVariableNames', 1);
curtable(:,1:2)=[];
curtable(:,1) = strcat(curtable.chr, '_',...
                                curtable.gene_id);
curtable.Properties.VariableNames{1} = 'species_genes';

% get annotation from single tables
for i=1:length(edgeResultsRNA)
    curcols = intersect(edgeResultsRNA{i}.Properties.VariableNames, curtable.Properties.VariableNames);
    [~, test, currows] = intersect(edgeResultsRNA{i}{:,1}, curtable{:,1}, 'stable');
    curtable(currows, curcols) = edgeResultsRNA{i}(:,curcols);
end
stringtieRawDNA_full = curtable;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read raw RNA stringtie
curtable = readtable([stringtieFolderRNA...
                          stringtieFileNamesRNA_raw{1}],...
                          opts,...
                          'ReadVariableNames', 1);
curtable(:,1:2)=[];
curtable(:,1) = strcat(curtable.chr, '_',...
                                curtable.gene_id);
curtable.Properties.VariableNames{1} = 'species_genes';
% add genome abbreviation and name columns
curtable.("genome_abbr") = cell(height(curtable),1);
curtable.("genome_name") = cell(height(curtable),1);

% get annotation from single tables
for i=1:length(edgeResultsRNA)
    curcols = intersect(edgeResultsRNA{i}.Properties.VariableNames, curtable.Properties.VariableNames);
    [~, test, currows] = intersect(edgeResultsRNA{i}{:,1}, curtable{:,1}, 'stable');
    curtable(currows, curcols) = edgeResultsRNA{i}(:,curcols);
    % add genome abbrevuation as in the file name
    % get genome name from the genome abbreviation table
    curname = genome_abbr_names.genome_name(ismember(genome_abbr_names.Species_abbr,...
                                            edgeRFileNames_abbr(i)));
    curtable(currows, "genome_abbr") = repmat(edgeRFileNames_abbr(i),length(currows),1);
    curtable(currows, "genome_name") = repmat(curname,length(currows),1);
end
stringtieRawRNA_full = curtable;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% combine all fold changes, p-values, pfdr and annotations in one matrix
% calculate total size
totalsize = 0;
for i=1:length(edgeResultsDNA)
    totalsize = totalsize+size(edgeResultsDNA{i},1);
end
fcMatrixHFDCTR_DNA = zeros(totalsize,1);
pMatrixHFDCTR_DNA = zeros(totalsize,1);
padjMatrixHFDCTR_DNA= zeros(totalsize,1);
fcMatrixHFDCTR_RNA = zeros(totalsize,1);
pMatrixHFDCTR_RNA = zeros(totalsize,1);
padjMatrixHFDCTR_RNA = zeros(totalsize,1);
% create matrices for deseq2 results 
fcMatrixDeseq2HFDCTR_DNA = zeros(totalsize,1);
pMatrixDeseq2HFDCTR_DNA = zeros(totalsize,1);
padjMatrixDeseq2HFDCTR_DNA= zeros(totalsize,1);
fcMatrixDeseq2HFDCTR_RNA = zeros(totalsize,1);
pMatrixDeseq2HFDCTR_RNA = zeros(totalsize,1);
padjMatrixDeseq2HFDCTR_RNA = zeros(totalsize,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get gene names and info
geneNamesHFDCTR = cell(totalsize,1);
abbrMatrix = cell(totalsize,1);
idx = 1;
for i=1:length(edgeRFileNames_abbr)
    % check gene ids and remove nans
    curidx = 1:length(edgeResultsDNA{i}.species_genes);
    curgeneidx = edgeResultsDNA{i}.genes;
    curidx(cellfun(@(x) isequal(x,'nan'), curgeneidx))=[];
    curgeneidx(cellfun(@(x) isequal(x,'nan'), curgeneidx))=[];
    cursize = length(curidx);
    % get DNA edgeR data
    fcMatrixHFDCTR_DNA(idx:idx+cursize-1) = edgeResultsDNA{i}.logFC(curidx);
    pMatrixHFDCTR_DNA(idx:idx+cursize-1) = edgeResultsDNA{i}.PValue(curidx);
    padjMatrixHFDCTR_DNA(idx:idx+cursize-1) = edgeResultsDNA{i}.FDR(curidx);
    % get RNA edgeR data
    % find corresponding gene ids
    [~, ~, rnaidx] = intersect(curgeneidx, edgeResultsRNA{i}.species_genes,'stable');
    
    fcMatrixHFDCTR_RNA(idx:idx+cursize-1) = edgeResultsRNA{i}.logFC(rnaidx);
    pMatrixHFDCTR_RNA(idx:idx+cursize-1) = edgeResultsRNA{i}.PValue(rnaidx);
    padjMatrixHFDCTR_RNA(idx:idx+cursize-1) = edgeResultsRNA{i}.FDR(rnaidx);

    % get DNA deseq2 data
    [~, ~, deseqidx] = intersect(curgeneidx, deseqResultsDNA{i}.species_genes,'stable');
    fcMatrixDeseq2HFDCTR_DNA(idx:idx+cursize-1) = deseqResultsDNA{i}.log2FoldChange(deseqidx);
    pMatrixDeseq2HFDCTR_DNA(idx:idx+cursize-1) = deseqResultsDNA{i}.pvalue(deseqidx);
    padjMatrixDeseq2HFDCTR_DNA(idx:idx+cursize-1) = deseqResultsDNA{i}.padj(deseqidx);
    % get RNA deseq data
    % find corresponding gene ids
    [~, ~, rnaidx] = intersect(curgeneidx, deseqResultsRNA{i}.species_genes,'stable');
    
    fcMatrixDeseq2HFDCTR_RNA(idx:idx+cursize-1) = deseqResultsRNA{i}.log2FoldChange(rnaidx);
    pMatrixDeseq2HFDCTR_RNA(idx:idx+cursize-1) = deseqResultsRNA{i}.pvalue(rnaidx);
    padjMatrixDeseq2HFDCTR_RNA(idx:idx+cursize-1) = deseqResultsRNA{i}.padj(rnaidx);

    % gene names and species abbreviations
    abbrMatrix(idx:idx+cursize-1) = repmat(edgeRFileNames_abbr(i),cursize,1);
    geneNamesHFDCTR(idx:idx+cursize-1) = curgeneidx;
   
    idx = idx+cursize;
end

% intersect full raw and normalized tables and combined FC table
% merge counts and results by gene id
[~,~,mergeidx] = intersect(geneNamesHFDCTR, ...
                           getmmNormDNA_full.species_genes,...
                           'stable');
countcols = contains(getmmNormDNA_full.Properties.VariableNames, 'MSZ');                       
countsMatrixGetMM_DNA = ...
            table2array(getmmNormDNA_full(mergeidx, countcols));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,~,mergeidx] = intersect(geneNamesHFDCTR, ...
                           getmmNormRNA_full.species_genes,...
                           'stable');
countcols = contains(getmmNormRNA_full.Properties.VariableNames, 'MSZ');                       
countsMatrixGetMM_RNA = ...
            table2array(getmmNormRNA_full(mergeidx, countcols));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,~,mergeidx] = intersect(geneNamesHFDCTR, ...
                           stringtieRawDNA_full.species_genes,...
                           'stable');
countcols = contains(stringtieRawDNA_full.Properties.VariableNames, 'MSZ');                       
countsMatrixRAW_DNA = ...
            table2array(stringtieRawDNA_full(mergeidx, countcols));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,~,mergeidx] = intersect(geneNamesHFDCTR, ...
                           stringtieRawRNA_full.species_genes,...
                           'stable');
countcols = contains(stringtieRawRNA_full.Properties.VariableNames, 'MSZ');                       
countsMatrixRAW_RNA = ...
            table2array(stringtieRawRNA_full(mergeidx, countcols));
annTable = stringtieRawRNA_full(mergeidx, ~countcols);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filter out rare transcipts                 
ndetected = sum(countsMatrixRAW_RNA~=0,2);
ndetectedThreshold = 5;
% filter DNA and RNA counts
filter_index = (ndetected>=ndetectedThreshold);

countsMatrixGetMM_DNA_filtered = countsMatrixGetMM_DNA(filter_index,:);
countsMatrixGetMM_RNA_filtered = countsMatrixGetMM_RNA(filter_index,:);
countsMatrixRAW_DNA_filtered = countsMatrixRAW_DNA(filter_index,:);
countsMatrixRAW_RNA_filtered = countsMatrixRAW_RNA(filter_index,:);

fcMatrixHFDCTR_DNA_filtered = fcMatrixHFDCTR_DNA(filter_index,:);
pMatrixHFDCTR_DNA_filtered = pMatrixHFDCTR_DNA(filter_index,:);
padjMatrixHFDCTR_DNA_filtered= padjMatrixHFDCTR_DNA(filter_index,:);
fcMatrixHFDCTR_RNA_filtered = fcMatrixHFDCTR_RNA(filter_index,:);
pMatrixHFDCTR_RNA_filtered = pMatrixHFDCTR_RNA(filter_index,:);
padjMatrixHFDCTR_RNA_filtered = padjMatrixHFDCTR_RNA(filter_index,:);

fcMatrixDeseq2HFDCTR_DNA_filtered = fcMatrixDeseq2HFDCTR_DNA(filter_index,:);
pMatrixDeseq2HFDCTR_DNA_filtered = pMatrixDeseq2HFDCTR_DNA(filter_index,:);
padjMatrixDeseq2HFDCTR_DNA_filtered= padjMatrixDeseq2HFDCTR_DNA(filter_index,:);
fcMatrixDeseq2HFDCTR_RNA_filtered = fcMatrixDeseq2HFDCTR_RNA(filter_index,:);
pMatrixDeseq2HFDCTR_RNA_filtered = pMatrixDeseq2HFDCTR_RNA(filter_index,:);
padjMatrixDeseq2HFDCTR_RNA_filtered = padjMatrixDeseq2HFDCTR_RNA(filter_index,:);

geneNamesHFDCTR_filtered = geneNamesHFDCTR(filter_index);
abbrMatrix_filtered = abbrMatrix(filter_index);
annTable_filtered = annTable(filter_index,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract pathway annotation from genes
keggPathways = cell(200,1);
keggPathways_geneidx = cell(200,1);
idx=1;
for i=1:size(annTable_filtered,1)
    %curpathway = annTable_filtered.pathway{i};
    %curpathway = strsplit(curpathway,';');
    curpathway = annTable_filtered.KEGG_Pathway{i}; %eggNOG
    curpathway = strsplit(curpathway,',');
    % remove starting with ko
    if length(curpathway)>1
        curpathway(cellfun(@(x) isequal(x(1:2),'ko'), curpathway))=[];
    end
    % remove duplicates if they exist
    curpathway = unique(curpathway);
    % check whether pathway is already in the array
    for j=1:length(curpathway)
        [~,checkidx] = intersect(keggPathways(1:idx-1),curpathway{j});
        if isempty(checkidx)
            % add new patway to the list and corresponding gene idx
            keggPathways{idx} = curpathway{j};
            keggPathways_geneidx{idx} = i;
            idx = idx+1;
        else
            % add gene idx to existing pathway
            keggPathways_geneidx{checkidx} = [keggPathways_geneidx{checkidx};...
                                              i];
        end
    end
end
keggPathways(idx:end) = [];
keggPathways_geneidx(idx:end) = [];
% remove empty pathways
keggPathways_geneidx(cellfun(@(x) ~contains(x,'map'),keggPathways))=[];
keggPathways(cellfun(@(x) ~contains(x,'map'),keggPathways))=[];
[keggPathways, sortidx] = sort(keggPathways);
keggPathways_geneidx = keggPathways_geneidx(sortidx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract EC annotation from genes
ecPathways = cell(2000,1);
ecPathways_geneidx = cell(2000,1);
idx=1;
for i=1:size(annTable_filtered,1)
    %curpathway = annTable_filtered.ec{i};
    %curpathway = strsplit(curpathway,';');
    curpathway = annTable_filtered.EC{i}; %eggNOG
    curpathway = strsplit(curpathway,',');
    % remove duplicates if they exist
    curpathway = unique(curpathway);
    % check whether pathway is already in the array
    for j=1:length(curpathway)
        [~,checkidx] = intersect(ecPathways(1:idx-1),curpathway{j});
        if isempty(checkidx)
            % add new patway to the list and corresponding gene idx
            ecPathways{idx} = curpathway{j};
            ecPathways_geneidx{idx} = i;
            idx = idx+1;
        else
            % add gene idx to existing pathway
            ecPathways_geneidx{checkidx} = [ecPathways_geneidx{checkidx};...
                                              i];
        end
    end
end
ecPathways(idx:end) = [];
ecPathways_geneidx(idx:end) = [];
% remove empty pathways
ecPathways_geneidx(cellfun(@(x) ~contains(x,'.'),ecPathways))=[];
ecPathways(cellfun(@(x) ~contains(x,'.'),ecPathways))=[];
ecPathways_geneidx(cellfun(@(x) length(x)<4,ecPathways))=[];
ecPathways(cellfun(@(x) length(x)<4,ecPathways))=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract GO annotation from genes
goPathways = cell(2000,1);
goPathways_geneidx = cell(2000,1);
idx=1;
for i=1:size(annTable_filtered,1)
    curpathway = annTable_filtered.GOs{i};
    curpathway = strsplit(curpathway,',');
    % remove duplicates if they exist
    curpathway = unique(curpathway);
    % check whether pathway is already in the array
    for j=1:length(curpathway)
        [~,checkidx] = intersect(goPathways(1:idx-1),curpathway{j});
        if isempty(checkidx)
            % add new patway to the list and corresponding gene idx
            goPathways{idx} = curpathway{j};
            goPathways_geneidx{idx} = i;
            idx = idx+1;
        else
            % add gene idx to existing pathway
            goPathways_geneidx{checkidx} = [goPathways_geneidx{checkidx};...
                                              i];
        end
    end
end
goPathways(idx:end) = [];
goPathways_geneidx(idx:end) = [];
% remove empty pathways
goPathways_geneidx(cellfun(@(x) ~contains(x,'GO'),goPathways))=[];
goPathways(cellfun(@(x) ~contains(x,'GO'),goPathways))=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract egggNOG annotation from genes
eggnogPathways = cell(10000,1);
eggnogPathways_geneidx = cell(10000,1);
idx=1;
for i=1:size(annTable_filtered,1)
    curpathway = annTable_filtered.matchingOGs{i};
    curpathway = strsplit(curpathway,',');
    % remove duplicates if they exist
    curpathway = unique(curpathway);
    % check whether pathway is already in the array
    for j=1:length(curpathway)
        [~,checkidx] = intersect(eggnogPathways(1:idx-1),curpathway{j});
        if isempty(checkidx)
            % add new patway to the list and corresponding gene idx
            eggnogPathways{idx} = curpathway{j};
            eggnogPathways_geneidx{idx} = i;
            idx = idx+1;
        else
            % add gene idx to existing pathway
            eggnogPathways_geneidx{checkidx} = [eggnogPathways_geneidx{checkidx};...
                                              i];
        end
    end
end
eggnogPathways(idx:end) = [];
eggnogPathways_geneidx(idx:end) = [];
% remove empty pathways
eggnogPathways_geneidx(cellfun(@(x) ~contains(x,'@'),eggnogPathways))=[];
eggnogPathways(cellfun(@(x) ~contains(x,'@'),eggnogPathways))=[];

% calculate how many species/parents are in eggNOG
eggnogPathways_parents = unique(cellfun(@(x) x(strfind(x,'@'):end), eggnogPathways, 'unif',0));

% get only COG pathways
cogPathways_geneidx = eggnogPathways_geneidx(cellfun(@(x) contains(x,'COG') & isequal(x(end-1:end),'@2'),eggnogPathways));
cogPathways = eggnogPathways(cellfun(@(x) contains(x,'COG') & isequal(x(end-1:end),'@2'),eggnogPathways));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save fold change tables to file
annTable_FCinfo = annTable;
% add information on species abbreviations
annTable_FCinfo.abbrSpecies = abbrMatrix;
% add filter
annTable_FCinfo.geneFilter = filter_index;
% add index after filtering for pathway info
index_after_filter = zeros(size(filter_index));
idx=1;
for i=1:length(filter_index)
    if filter_index(i)==1
        index_after_filter(i) = idx;
        idx = idx+1;
    end
end
annTable_FCinfo.indexFiltered = index_after_filter;
% add fold change, p-values and fdr
%dna
annTable_FCinfo.fcHFDCTR_DNA = fcMatrixHFDCTR_DNA;
annTable_FCinfo.pHFDCTR_DNA = pMatrixHFDCTR_DNA;
annTable_FCinfo.fdrHFDCTR_DNA = padjMatrixHFDCTR_DNA;
%rna
annTable_FCinfo.fcHFDCTR_RNA = fcMatrixHFDCTR_RNA;
annTable_FCinfo.pHFDCTR_RNA = pMatrixHFDCTR_RNA;
annTable_FCinfo.fdrHFDCTR_RNA = padjMatrixHFDCTR_RNA;
% add fold change, p-values and fdr from deseq2
%dna
annTable_FCinfo.fcDeseq2HFDCTR_DNA = fcMatrixDeseq2HFDCTR_DNA;
annTable_FCinfo.pDeseq2HFDCTR_DNA = pMatrixDeseq2HFDCTR_DNA;
annTable_FCinfo.fdrDeseq2HFDCTR_DNA = padjMatrixDeseq2HFDCTR_DNA;
%rna
annTable_FCinfo.fcDeseq2HFDCTR_RNA = fcMatrixDeseq2HFDCTR_RNA;
annTable_FCinfo.pDeseq2HFDCTR_RNA = pMatrixDeseq2HFDCTR_RNA;
annTable_FCinfo.fdrDeseq2HFDCTR_RNA = padjMatrixDeseq2HFDCTR_RNA;
% save table to file
writetable(annTable_FCinfo, [outputFolder 'table_S4_edgeR_deseq2_gene_fold_changes_and_ann.csv']);
%writetable(annTable_FCinfo, [outputFolder 'edgeR_deseq2_gene_fold_changes_and_ann_20250713_newV.csv']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write pathway-gene relationships to file
% kegg
fid = fopen([outputFolder 'gene_annotation_kegg.tsv'], 'w');
fprintf(fid, 'PathwayID\tGeneIDXfiltered\n');
for i=1:length(keggPathways)
    fprintf(fid, '%s\t', keggPathways{i});
    fprintf(fid, '%d;', keggPathways_geneidx{i}(1:end));
    fprintf(fid, '\n');
end
fclose(fid);
% go
fid = fopen([outputFolder 'gene_annotation_go.tsv'], 'w');
fprintf(fid, 'PathwayID\tGeneIDXfiltered\n');
for i=1:length(goPathways)
    fprintf(fid, '%s\t', goPathways{i});
    fprintf(fid, '%d;', goPathways_geneidx{i}(1:end));
    fprintf(fid, '\n');
end
fclose(fid);
% EC
fid = fopen([outputFolder 'gene_annotation_EC.tsv'], 'w');
fprintf(fid, 'PathwayID\tGeneIDXfiltered\n');
for i=1:length(ecPathways)
    fprintf(fid, '%s\t', ecPathways{i});
    fprintf(fid, '%d;', ecPathways_geneidx{i}(1:end));
    fprintf(fid, '\n');
end
fclose(fid);
% EggNOG
fid = fopen([outputFolder 'gene_annotation_EggNOG.tsv'], 'w');
fprintf(fid, 'PathwayID\tGeneIDXfiltered\n');
for i=1:length(eggnogPathways)
    fprintf(fid, '%s\t', eggnogPathways{i});
    fprintf(fid, '%d;', eggnogPathways_geneidx{i}(1:end));
    fprintf(fid, '\n');
end
fclose(fid);
% COG
fid = fopen([outputFolder 'gene_annotation_COG.tsv'], 'w');
fprintf(fid, 'PathwayID\tGeneIDXfiltered\n');
for i=1:length(cogPathways)
    fprintf(fid, '%s\t', cogPathways{i});
    fprintf(fid, '%d;', cogPathways_geneidx{i}(1:end));
    fprintf(fid, '\n');
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save gene counts matrices to files
postfix = '';
% get column names
countcols = contains(getmmNormDNA_full.Properties.VariableNames, 'MSZ');                       
countcols = getmmNormDNA_full.Properties.VariableNames(countcols);
countcols = cellfun(@(x) x(1:strfind(x,'_')-1), countcols, 'unif', 0);
% getMM DNA
countsMatrixGetMM_DNA_table = array2table(countsMatrixGetMM_DNA,...
                                    'VariableNames',countcols);
countsMatrixGetMM_DNA_table.species_genes = annTable.species_genes;
countsMatrixGetMM_DNA_table.geneFilter = annTable_FCinfo.geneFilter;
writetable(countsMatrixGetMM_DNA_table, [outputFolder 'countsMatrixGetMM_DNA_table' postfix '.txt']);
% getMM RNA
countsMatrixGetMM_RNA_table = array2table(countsMatrixGetMM_RNA,...
                                    'VariableNames',countcols);
countsMatrixGetMM_RNA_table.species_genes = annTable.species_genes;
countsMatrixGetMM_RNA_table.geneFilter = annTable_FCinfo.geneFilter;
writetable(countsMatrixGetMM_RNA_table, [outputFolder 'countsMatrixGetMM_RNA_table' postfix '.txt']);
% raw DNA
countsMatrixRAW_DNA_table = array2table(countsMatrixRAW_DNA,...
                                    'VariableNames',countcols);
countsMatrixRAW_DNA_table.species_genes = annTable.species_genes;
countsMatrixRAW_DNA_table.geneFilter = annTable_FCinfo.geneFilter;
writetable(countsMatrixRAW_DNA_table, [outputFolder 'countsMatrixRAW_DNA_table' postfix '.txt']);
% raw RNA
countsMatrixRAW_RNA_table = array2table(countsMatrixRAW_RNA,...
                                    'VariableNames',countcols);
countsMatrixRAW_RNA_table.species_genes = annTable.species_genes;
countsMatrixRAW_RNA_table.geneFilter = annTable_FCinfo.geneFilter;
writetable(countsMatrixRAW_RNA_table, [outputFolder 'countsMatrixRAW_RNA_table' postfix '.txt']);

% write annotations to file
writetable(annTable_filtered, [outputFolder 'geneAnnTable_filtered' postfix '.csv']);
writetable(annTable, [outputFolder 'geneAnnTable_full' postfix '.csv']);

