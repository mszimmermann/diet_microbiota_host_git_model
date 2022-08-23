function [enrichmentScores, groupLabels, enrichmentTable] = ptwEnrichmentMetabolitesIonsPTWTable(ptw, mets, groups, ptwNames)
% ptw is a matrix of pathways vs metabolites
% ptwNames is pathway names
groupLabels = unique(groups);
nGroups = length(groupLabels);

nDetectedGenes = length(mets);

enrichmentScores = zeros(size(ptw,1), nGroups);
enrichmentTable = cell(size(ptw,1), 7);
% calculate the enrichment score for each pathway for each group
for iPw = 1:size(ptw,1)
    % if the pathway has annotated ions
    if nnz(ptw(iPw,:))
                   
        % calculate how many unique ions are in the pathway
        nPwDetectedGenes = nnz(ptw(iPw,:)); % this is K - number of pathway genes
       
        bestScore = 1; % probability to draw this number of genes if H0 is true
        % this is k - number of pathway genes which are in the group
        genesFound  = nnz(groups(ptw(iPw,:)==1)==1);
        nGroup = nnz(groups==1); % number of ions in the group
        % calculate P(x<=k) - hypergeometric distribution f(k,N,K,n)
        % where N - total number of genes, n - size of the group

        score = sum( hygepdf(genesFound:nPwDetectedGenes, nDetectedGenes, nPwDetectedGenes, nGroup) );
        if score < bestScore
            bestScore = score;
        end
        enrichmentScores(iPw, :) = bestScore;
        enrichmentTable{iPw, 1} = bestScore;
        enrichmentTable{iPw, 2} = bestScore;
        enrichmentTable{iPw, 3} = ptwNames(iPw);
        enrichmentTable{iPw, 4} = genesFound;
        enrichmentTable{iPw, 5} = nGroup;
        enrichmentTable{iPw, 6} = nPwDetectedGenes;
        enrichmentTable{iPw, 7} = nDetectedGenes;
    else
        enrichmentScores(iPw,:) = 1;

        enrichmentTable{iPw, 1} = 1;
        enrichmentTable{iPw, 2} = 1;
        enrichmentTable{iPw, 3} = ptwNames(iPw);
        enrichmentTable{iPw, 4} = 0;
        enrichmentTable{iPw, 5} = 0;
        enrichmentTable{iPw, 6} = 0;
        enrichmentTable{iPw, 7} = nDetectedGenes;
    end
end

% multiple hypothesis adjustment
% padj = mafdr(cell2mat(enrichmentTable(:,1)), 'BHFDR', 1);
% for i=1:length(padj)
%     enrichmentTable{i,2} = padj(i);
% end
%multiple hypothesis adjustment -my BH
padj = my_bhfdr(cell2mat(enrichmentTable(:,1)));
for i=1:length(padj)
    enrichmentTable{i,2} = padj(i);
end