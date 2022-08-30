%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load kegg reaction info (ec and substrates-products)

% call script defining file dependenciesand global variables
addpath(genpath('.\'));
add_global_and_file_dependencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data requirements: 
% 'keggreactionsALLreactantsenzymes.mat'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output:
% Files:
%
% Figures:
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load([processeddataFolder, filesep, 'util', filesep,...
    'keggreactionsALLreactantsenzymes.mat'])

% for each kegg reactions, make a table of reaction ID, ecID, substrate and
% product
kegg_rxn_ec_subprod = cell(length(keggreactionsALL)*5, 4);
kegg_idx = 1;
for i=1:length(keggreactionsALL)
    cur_ec = keggreactionsALL{i}.ecnumbers;
    cur_subprod = keggreactionsALL{i}.sub_prod;
    if (numel(cur_ec)>0) && (numel(cur_subprod)>0)
        [~,idx]=unique(cell2mat(cur_subprod),'rows');
        cur_subprod =  cur_subprod(idx,:);
        cur_ec = unique(cur_ec);
        for j=1:length(cur_ec)
            for k=1:size(cur_subprod,1)
                kegg_rxn_ec_subprod{kegg_idx,1} = keggrxn{i};
                kegg_rxn_ec_subprod{kegg_idx,2} = cur_ec{j};
                kegg_rxn_ec_subprod(kegg_idx,3:4) = cur_subprod(k,:);
                kegg_idx = kegg_idx+1;
            end
        end
    end
end
kegg_rxn_ec_subprod(kegg_idx:end,:)=[];        

save([outputFolder, 'kegg_rxn_ec_subprod.mat'],...
    'kegg_rxn_ec_subprod');

% save as table
table_kegg_rxn_ec_subprod = cell2table(kegg_rxn_ec_subprod,...
    'VariableNames', {'KEGGrxn', 'KEGGec', 'KEGGsubstrate', 'KEGGproduct'});
writetable(table_kegg_rxn_ec_subprod,...
    [outputFolder, 'table_kegg_rxn_ec_substrate_product.csv']);
