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

shortestPathTable_unique = readtable([outputFolder,...
    'table_shortest_path_subprod_sorted_unique.csv']);

% keep only path of length 1
shortestPathTable_unique = shortestPathTable_unique(...
    shortestPathTable_unique.Path_length==1,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read metabolite data
   
metaboliteFilters = readtable([inputFolder,...
    'metabolites_allions_combined_formulas_with_metabolite_filters.csv']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get indeces of metabolites of interest
select_mets_ID = [shortestPathTable_unique.Substrate; shortestPathTable_unique.Product];
select_mets_IDX = [shortestPathTable_unique.Substrate_IDX;...
                   shortestPathTable_unique.Product_IDX];
% remove redundant metabolites
[select_mets_ID, uniqueidx] = unique(select_mets_ID);
select_mets_IDX = select_mets_IDX(uniqueidx);

select_mets_IDX_unique = unique(select_mets_IDX);
select_metaboliteTable = metaboliteFilters(select_mets_IDX_unique,:);

writetable(select_metaboliteTable, ...
    [outputFolder, 'table_metabolites_substrate_products_one_enzyme_selected.csv'],...
    'WriteRowNames',true);
