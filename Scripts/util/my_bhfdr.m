function fdr = my_bhfdr(p)
% Adopted from mafdr function in bioinformatics toolbox
p = reshape(p, [], 1); % ake sure p is a column
m = numel(p);

% get ranks counting the max rank for tied values 
p_ranks = ceil(tiedrank(p));
p_ranks = reshape(p_ranks, [], 1); % make sure it's a column

% based on this code example:
% https://www.r-bloggers.com/2023/07/the-benjamini-hochberg-procedure-fdr-and-p-value-adjusted-explained/
% and comparison with mafdr function results

% for each tied rank, make sure the highest possible rank is assigned
p_ranks_unique = unique(p_ranks);
if length(p_ranks_unique) < length(p_ranks)
    for i=length(p_ranks_unique)-1:-1:1
        % for all items with the current rank, set it to next biggest rank minus 1
        % only for cases where multiple p-values with the same value exist
        if nnz(p_ranks==p_ranks_unique(i))>1
            p_ranks(p_ranks==p_ranks_unique(i)) = p_ranks_unique(i+1)-1;
        end
    end
end

fdr_ord =  p(:) .* (m./(p_ranks));

% Running min in reverse order (in-place)
%bioinfoprivate.cumminmaxmex(fdr_ord,'min','reverse');
% check if the first FDR is already 1, then set everything to 1
fdr = ones(size(p));
for i=1:length(p)
    % get the rank
    p_rank = p_ranks(i);
    % get all the adjusted values for p-values that have greater or equal 
    % to this rank and get the min value
    fdr(i) = min(1, min(fdr_ord(p_ranks>=p_rank)));
end

% testing
% pvalues = [0.01,0.001, 0.05, 0.20, 0.15, 0.15]
% expected FDR = [0.030 0.006 0.100 0.200 0.180 0.180]
