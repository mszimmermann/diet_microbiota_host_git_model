function [ptwIonTable, ptwNames] = create_pathway_ion_table(ptw, mets)
% given ptw - pathway structure, make binary incidense matrices for
% metabolites in different pathways

% create pathway-metabolite table
pathway_mets = cell(200000,1);
ptwNames = cell(length(ptw),1);
ptwTable = zeros(length(ptw), length(pathway_mets));
idx=1;
for iPw = 1:length(ptw)
    % if the pathway has a field 'geneName'
    ptwNames{iPw} = ptw(iPw).ptw.name;
    curcompounds = ptw(iPw).ptw.cmpdID;
    
    if idx==1
        pathway_mets(idx:idx+length(curcompounds)-1) = curcompounds;
        ptwTable(iPw, idx:idx+length(curcompounds)-1) = 1;
        idx = idx + length(curcompounds);
        
    else
        [~,existing_compounds] = intersect(pathway_mets(1:idx-1),...
            curcompounds);
        ptwTable(iPw, existing_compounds) = 1;
        
        % add new
        newcompounds = setdiff(curcompounds, pathway_mets(1:idx-1));
        pathway_mets(idx:idx+length(newcompounds)-1) = newcompounds;
        ptwTable(iPw, idx:idx+length(newcompounds)-1) = 1;
        idx = idx + length(newcompounds);
    end
end
ptwTable(:, idx:end) = [];
pathway_mets(idx:end) = [];

ptwIonTable = zeros(length(ptwNames), length(mets));
% split the genes with several names into separate entries
for i=1:length(mets)
    curmets = strsplit(mets{i}, ';');
    [~, ~, hmdbmetspos] = intersect(curmets, pathway_mets);
    ptwIonTable(:, i) = sum(ptwTable(:, hmdbmetspos), 2);
end

ptwIonTable(ptwIonTable>1)=1;

% remove pathways without ions
nIons = sum(ptwIonTable,2);
ptwNames(nIons==0)=[];
ptwIonTable(nIons==0,:)=[];
