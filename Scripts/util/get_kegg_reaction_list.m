function [keggreactions] = get_kegg_reaction_list()

%get list of enzymes
keggurl = 'http://rest.kegg.jp/list/reaction';
try
    ret = webread(keggurl); %retrieve reaction list
    ret = strsplit(ret, '\n');
    ret(end) = [];

    keggreactions = cell(length(ret), 2);
    for i=1:length(ret)
        currxn = strsplit(ret{i}, '\t');
        keggreactions{i,1} = currxn{1};
        keggreactions{i,2} = currxn{2};
    end
catch
     disp('KEGG connection error\n')
end