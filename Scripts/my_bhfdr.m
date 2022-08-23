function fdr = my_bhfdr(p)
% Adopted from mafdr function in bioinformatics toolbox
m = numel(p);
[p_ord, idx] = sort(p);
fdr_ord =  p_ord(:) .* (m./(1:m))';
% Running min in reverse order (in-place)
%bioinfoprivate.cumminmaxmex(fdr_ord,'min','reverse');
fdr(idx) = fdr_ord;
fdr(fdr>1) = 1;
