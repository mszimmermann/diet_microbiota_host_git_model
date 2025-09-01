function [exactMass, ncarbon] = calculateExactMass(formulas, h_isotope)
%calculate exact masses given molecular formula
CHNOPS = zeros(length(formulas), 6);
atomMasses = [  12.000000;...
                1.007825;...
                14.003074;...
                15.994915;...
                30.973763;...
                31.972072;...
                ];
atomNames = {'C' 'H' 'N' 'O' 'P' 'S'};   

formula_passcheck = ones(length(formulas),1);
for i=1:length(formulas)
    %formula_isalphanum = isstrprop(formulas{i}, 'alphanum');
    formula_isnumeric = (formulas{i}>='0')&(formulas{i}<='9');
    formula_is_allowed_atoms = zeros(1,length(formulas{i}));
    for j=1:length(atomNames)
        formula_is_allowed_atoms = formula_is_allowed_atoms |...
            ((formulas{i}-atomNames{j})==0);
    end
    formula_check = formula_isnumeric | formula_is_allowed_atoms;
    % if there are unallowed characters, continue without processing the
    % formula
    if sum(formula_check==0)>0
       % fprintf('Can not calculate mass for %s\n', formulas{i});
        formula_passcheck(i) = 0;
        continue
    end
    %fprintf('Calculating mass for %s\n', formulas{i});
        
    curatomNamePos = cellfun(@(x) find((formulas{i}-x)==0), atomNames, 'UniformOutput', 0);
    curatomNames = find(cellfun(@(x) ~isempty(x), curatomNamePos));
    curatomNamePos = sort(cell2mat(curatomNamePos(curatomNames)));
    % get the position of each atom in the CHNOPS matrix
    curatomNames = arrayfun(@(x) find(ismember(atomNames, formulas{i}(x))),curatomNamePos);
    
    curatomQuantities = zeros(size(CHNOPS,2),1);
    for j=1:length(curatomNamePos)
        if j==length(curatomNamePos)
            curquant = formulas{i}(curatomNamePos(j)+1:end);
            if curatomNamePos(j) == length(formulas{i})
                curquant = '1';
            end
        else
            if curatomNamePos(j)+1 == curatomNamePos(j+1)
                curquant = '1';
            else
                curquant = formulas{i}(curatomNamePos(j)+1:curatomNamePos(j+1)-1);
            end
        end
        % add instead of setting in case formulas contain multiple
        % instances of the same atom
        curatomQuantities(curatomNames(j)) = curatomQuantities(curatomNames(j))  + str2double(curquant);
    end
    CHNOPS(i,:) = curatomQuantities;
end
%subtract hydrogen isotope
CHNOPS(:,2) = CHNOPS(:,2) + h_isotope;
exactMass = CHNOPS*atomMasses;
ncarbon = CHNOPS(:,1);
        