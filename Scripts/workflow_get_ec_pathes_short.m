%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get EC pathes from substrate-product pathes

% call script defining file dependenciesand global variables
addpath(genpath('.\'));
add_global_and_file_dependencies

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Required files:
% 'table_kegg_rxn_ec_substrate_product.csv'
% 'table_shortestMatrix_paths_filtered.csv'
% Output: 
% Files:
% 'shortest_paths_upto_length_4.csv'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kegg_rxn_ec_subprod = readtable([resultsFolder,...
    'table_kegg_rxn_ec_substrate_product.csv']);
kegg_rxn_ec_subprod = table2cell(kegg_rxn_ec_subprod);

% [fid, errmsg] = fopen([resultsFolder,...
%     'table_shortestMatrix_paths_filtered.csv'], 'r');
% tline = fgetl(fid);
% maxlen = 0;
% while ischar(tline)
%     linelen = length(strsplit(tline, ','));
%     if linelen>maxlen
%         maxlen=linelen;
%     end
%     tline = fgetl(fid);
% end
% fclose(fid);

strict_class = 0; %1
% [fid, errmsg] = fopen([resultsFolder,...
%     'table_shortestMatrix_paths_filtered.csv'], 'r');
[fid, errmsg] = fopen([outputFolder, filesep, ...
    'table_shortestMatrix_paths_filtered_combined_IP_LI_PCC_within_high_total_strictclass'...
    num2str(strict_class) '.csv'],...
    'r');
% scan all the length including zeros/empty values as each line has
% variable length corresponding to the path length
shortestMatrix_cluster_paths = textscan( fid,'%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d' ...
            ,   'CollectOutput' ,   true    ...
            ,   'Delimiter'     , ',' ...
            ,'EmptyValue',0,'HeaderLines',1);
[~] = fclose( fid );
shortestMatrix_cluster_paths = shortestMatrix_cluster_paths{1};
shortestMatrix_cluster_paths_row_column = shortestMatrix_cluster_paths(:,1:2);
shortestMatrix_cluster_paths = shortestMatrix_cluster_paths(:,3:end);



rnxCompounds = readtable([resultsFolder, 'table_RXNCompounds.csv']);
rnxCompounds = table2cell(rnxCompounds);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate enzymatic length of paths
shortestMatrix_cluster_paths_length = sum(shortestMatrix_cluster_paths>0,2);

% keep only paths with length 1-4
keep_paths = (shortestMatrix_cluster_paths_length>0) &...
             (shortestMatrix_cluster_paths_length<5);


shortestMatrix_cluster_paths = shortestMatrix_cluster_paths(keep_paths,:);
shortestMatrix_cluster_paths_row_column = shortestMatrix_cluster_paths_row_column(keep_paths,:);
shortestMatrix_cluster_paths_length = shortestMatrix_cluster_paths_length(keep_paths);


% get ec corresponding to each metabolite pair
shortestMatrix_enzymes_filtered = cell(size(shortestMatrix_cluster_paths_length));
shortestMatrix_compounds_filtered = cell(size(shortestMatrix_cluster_paths_length));
shortestMatrix_substrates_filtered = cell(size(shortestMatrix_cluster_paths_length));
shortestMatrix_products_filtered = cell(size(shortestMatrix_cluster_paths_length));
for i=1:size(shortestMatrix_cluster_paths,1)
    curpaths = shortestMatrix_cluster_paths(i,:);
    % remove zeros
    curpaths(curpaths==0)=[];
    curpaths_compounds = rnxCompounds(curpaths);
    curpaths_enzymes = cell(length(curpaths_compounds)-1,1);
   
    for j=1:length(curpaths_compounds)-1
        cursub = curpaths_compounds{j};
        curprod = curpaths_compounds{j+1};
        cur_ec_idx = find((ismember(kegg_rxn_ec_subprod(:,3), cursub) & ...
                            ismember(kegg_rxn_ec_subprod(:,4), curprod)) |...
                            (ismember(kegg_rxn_ec_subprod(:,3), curprod) & ...
                            ismember(kegg_rxn_ec_subprod(:,4), cursub)) );
        if isempty(cur_ec_idx)
            curpaths_enzymes{j}='NoEC';
        else
            curpaths_enzymes{j} = strjoin(unique(kegg_rxn_ec_subprod(cur_ec_idx,2)),'|');
        end
    end
    shortestMatrix_enzymes_filtered{i} = strjoin(curpaths_enzymes, ';');
    shortestMatrix_compounds_filtered{i} = strjoin(curpaths_compounds, ';');
    shortestMatrix_substrates_filtered{i} = curpaths_compounds{1};
    shortestMatrix_products_filtered{i} = curpaths_compounds{end};
   
    if mod(i,1000)==0
        fprintf('Calculated enzymes for %d of %d substrates\n',i, size(shortestMatrix_cluster_paths,1));
    end
end

shortesMatrix_table = table(shortestMatrix_substrates_filtered,...
                             shortestMatrix_products_filtered,...
                             shortestMatrix_compounds_filtered,...
                             shortestMatrix_enzymes_filtered,...
                             'VariableNames',...
                             {'Substrate','Product','Compound_path', 'EC_path'});
%writetable(shortesMatrix_table, [outputFolder, 'shortest_paths_upto_length_4.csv']);                         
writetable(shortesMatrix_table, [outputFolder, 'shortest_paths_combined_IP_LI_PCC_within_high_total_strictclass'...
    num2str(strict_class) '_upto_length_4.csv']);                         