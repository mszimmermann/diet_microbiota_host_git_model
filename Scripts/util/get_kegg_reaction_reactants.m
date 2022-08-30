function [keggreactions] = get_kegg_reaction_reactants(keggrxn)

%get list of enzymes
keggurl = 'http://rest.kegg.jp/get/';
keggreactions = cell(size(keggrxn));
for i=1:length(keggrxn)
    try
        ret = webread([keggurl, keggrxn{i}]); %retrieve all the pathways in the org
        ret = strsplit(ret, '\n');
        ret(end) = [];

        equation = ret(cellfun(@(x) contains(x, 'EQUATION'),ret));
        if ~isempty(equation)
            equation = equation{1};
            compoundidx = strfind(equation, 'C');
            equationsign = strfind(equation, '>');
            substrateidx = compoundidx(compoundidx<equationsign);
            productidx = compoundidx(compoundidx>equationsign);
            substrates = arrayfun(@(x) equation(x:(x+5)), substrateidx, 'unif', 0);
            products = arrayfun(@(x) equation(x:(x+5)), productidx, 'unif', 0);
            sub_prod_pairs = cell(length(substrates)*length(products),2);
            idx = 1;
            % check whether there is rpair
            rpairidx = find(cellfun(@(x) contains(x, 'RCLASS'),ret));
            if ~isempty(rpairidx)
                for j=rpairidx:length(ret)
                    if contains(ret{j}, 'RCLASS') || (ret{j}(1)==' ')
                        rpair = ret{j};
                        rpair = strrep(rpair, 'RCLASS', '');
                        rpair = strrep(rpair, 'RC', '');
                        compoundidx = strfind(rpair, 'C');
                        sub_prod = [arrayfun(@(x) rpair(x:(x+5)), compoundidx(1:2:end), 'unif', 0)'...
                                    arrayfun(@(x) rpair(x:(x+5)), compoundidx(2:2:end), 'unif', 0)'];
                        % check whether substrates and products from
                        % different pairs are on the same side of equation
                        for spi = 1:size(sub_prod,1)
                            if ~ismember(sub_prod{spi,1}, substrates)
                                temp = sub_prod{spi,1};
                                sub_prod{spi,1} = sub_prod{spi,2};
                                sub_prod{spi,2} = temp;
                            end
                        end     
                        sub_prod_pairs(idx:idx+size(sub_prod,1)-1,:) = sub_prod;
                        idx = idx+size(sub_prod,1);
                    else
                        break;
                    end
                end
                % remove sub prod from substrates and products
%                 substrates = setdiff(substrates, sub_prod(:));
%                 products = setdiff(products, sub_prod(:));
            end
            % only use RCLASS reactants
%             for j=1:length(substrates)
%                 for k=1:length(products)
%                     sub_prod_pairs{idx,1} = substrates{j};
%                     sub_prod_pairs{idx,2} = products{k};
%                     idx = idx+1;
%                 end
%             end
            sub_prod_pairs(idx:end,:) = [];
            % get enzymes
            ecnumbers = '';
            ecidx = find(cellfun(@(x) contains(x, 'ENZYME'),ret));
            if ~isempty(ecidx)
                for j=ecidx:length(ret)
                    if contains(ret{j}, 'ENZYME') || (ret{j}(1)==' ')
                        rpair = ret{j};
                        rpair = strrep(rpair, 'ENZYME', '');
                        rpair  = strtrim(rpair);
                        ecnumbers = [ecnumbers strsplit(rpair)];
                    else
                        break;
                    end
                end
            end
            rxninfo.ecnumbers = ecnumbers;
            rxninfo.sub_prod = sub_prod_pairs;
            keggreactions{i} = rxninfo;
        end
    catch
         disp(['KEGG connection error ' keggrxn{i}])
    end
end