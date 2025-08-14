function fdr = my_bhfdr(p)
% Adopted from mafdr function in bioinformatics toolbox
p = reshape(p, [], 1); % ake sure p is a column
m = numel(p);

% get ranks counting the max rank for tied values 
p_ranks = ceil(tiedrank(p));
p_ranks = reshape(p_ranks, [], 1); % make sure it's a column

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
