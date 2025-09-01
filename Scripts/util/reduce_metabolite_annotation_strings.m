function [CompoundID_updated, metabolites_unique] = reduce_metabolite_annotation_strings(CompoundID)
% reduce metabolite annotation:
% for items present in multiple types of annotatin strings, keep the
% shorter one

CompoundID_updated = CompoundID;

maxiter = 10; %maximal number of iterations

while maxiter>0

    % make a list of unique metabolites
    metabolites_unique = unique(CompoundID_updated);
    metabolites_unique(cellfun(@(x) isempty(x), metabolites_unique))=[];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fore each metabolite_unique, check if it is a substring of any other
    metabolites_unique_issubstr = zeros(length(metabolites_unique), 1);
    for i=1:length(metabolites_unique)
        metabolites_unique_issubstr(i) = sum(cellfun(@(x) contains(x, metabolites_unique(i)),...
            metabolites_unique));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fore each metabolite_unique, which is a substring, keep the shorter string and remove it from the longer string
    metabolites_unique_tosearch = metabolites_unique(metabolites_unique_issubstr>1);
    %%%%%%%%%%%%% exit loop if there are no substrings
    if isempty(metabolites_unique_tosearch)
        break;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    metabolites_unique_toupdate = cell(length(metabolites_unique_tosearch)*5, 3);
    idx=1;
    for i=1:length(metabolites_unique_tosearch)
        containing_mets = metabolites_unique(cellfun(@(x) (contains(x, metabolites_unique_tosearch(i)) &...
                                        (length(x) > length(metabolites_unique_tosearch{i}))),...
                                        metabolites_unique));
        containing_mets_toupdate = cellfun(@(x) strrep(x, metabolites_unique_tosearch{i}, ''),...
                                        containing_mets, 'unif', 0);
        metabolites_unique_toupdate(idx:(idx+length(containing_mets)-1),:) = ...
            [repmat(metabolites_unique_tosearch(i), length(containing_mets),1) containing_mets containing_mets_toupdate];
        idx = idx+length(containing_mets_toupdate);
    end
    metabolites_unique_toupdate(idx:end,:)=[];
    
    % sort metabolites_unique_toupdate in the order of length of the sequence
    % that should be replaced
    metabolites_unique_toupdate_strlen = cellfun(@(x) length(x), metabolites_unique_toupdate(:,2));
    [~, sortidx] = sort(metabolites_unique_toupdate_strlen, 'descend');
    metabolites_unique_toupdate = metabolites_unique_toupdate(sortidx,:);
    % replace annotation of compounds
    for i=1:size(metabolites_unique_toupdate,1)
        CompoundID_updated(cellfun(@(x) ismember(x, metabolites_unique_toupdate(i,2)),...
                    CompoundID_updated)) = metabolites_unique_toupdate(i,3);
    end
    % check if unique metabolites improved in the next loop iteration
    sprintf('Countdown from %d: reduced %d substrings\n', maxiter, length(metabolites_unique_tosearch))
    maxiter = maxiter-1; % to avoid endless loops
end
if ~isempty(metabolites_unique_tosearch)
    sprintf('Used 10 iterations but not all annotation strings are reduced, %d are left\n', length(metabolites_unique_tosearch))
else
    sprintf('Reduced substrings counting down to %d iterations\n', maxiter)
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check if there are compounds that are in multiple annotations
% check of there are cases when the same ID is part of two or more unique
% metabolite annotation strings (which are not substrings of each other as
% a whole)
maxiter = 10; %maximal number of iterations

while maxiter>0 % one iteration should be sufficient but it is important to check
    metabolites_unique = unique(CompoundID_updated);
    metabolites_IDs = cell(length(metabolites_unique)*5,1);
    idx = 1;
    % make a matrix indicating in which metabolites_unique each ID is included
    metabolites_IDs_annotation_mat = zeros(length(metabolites_unique)*5,1);
    % record how many metabolites are present in one annotation string
    metabolites_IDs_annotation_nmets = zeros(length(metabolites_unique),1);
    for i=1:length(metabolites_unique)
        curmets = strsplit(metabolites_unique{i}, ';');
        % remove empty
        curmets(cellfun(@(x) isempty(strtrim(x)), curmets)) = [];
        metabolites_IDs(idx:idx+length(curmets)-1) = curmets;
        metabolites_IDs_annotation_mat(idx:idx+length(curmets)-1) = i;
        metabolites_IDs_annotation_nmets(i) = length(curmets);
        idx = idx+length(curmets);
    end
    metabolites_IDs(idx:end) = [];
    metabolites_IDs_annotation_mat(idx:end) = [];
    [metabolites_IDs_unique, ~, idxAll] = unique(metabolites_IDs, 'stable');
    % remove empty
    empty_idx = find(cellfun(@(x) isempty(strtrim(x)), metabolites_IDs_unique));
    if empty_idx
        sprintf("Warning: empty metabolites not removed")
    end
    % for each metabolite ID check if it is in multiple annotation strings
    metabolites_IDs_unique_num_annotations = zeros(size(metabolites_IDs_unique));
    metabolites_IDs_unique_annotations_idx = cell(size(metabolites_IDs_unique));
    for i=1:length(metabolites_IDs_unique_num_annotations)
        metabolites_IDs_unique_num_annotations(i) = ...
            length(unique(metabolites_IDs_annotation_mat(idxAll==i)));
        metabolites_IDs_unique_annotations_idx{i} = unique(metabolites_IDs_annotation_mat(idxAll==i));
    end
    % for each metabolite, keep only one annotation
    % filter unique metabolites plus metabolite IDs that belong to several
    % annotation strings
    metabolites_IDs_multiple = metabolites_IDs_unique(metabolites_IDs_unique_num_annotations>1);
    %%%%%%%%%%%%% exit loop if there are no multiple annotations for one compound
    if isempty(metabolites_IDs_multiple)
        break;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    metabolites_IDs_multiple_annotations_idx = metabolites_IDs_unique_annotations_idx(metabolites_IDs_unique_num_annotations>1);
    % keep annotation in the shorter annotation string
    % go through unique combinations of annotation idx
    metabolite_names_toreplace = cell(length(metabolites_IDs_multiple),2);
    curidx=1;
    for i=1:length(metabolites_IDs_multiple)
        curannidx = metabolites_IDs_multiple_annotations_idx{i};
        curannlength = metabolites_IDs_annotation_nmets(curannidx);
        keepannidx = curannidx(curannlength==min(curannlength));
        % check if there are multiple, keep one
        if length(keepannidx)>1
            keepannidx = keepannidx(1);
        end
        % set the other indeces to remove the current ID
        removeannidx = setdiff(curannidx, keepannidx);
        % iterate over and replace all in mdified compound ID
        sprintf('Cleaning metabolite %s keep %d\n', metabolites_IDs_multiple{i}, keepannidx)
        sprintf('Keeping %s\n',metabolites_unique{keepannidx})
        for j=1:length(removeannidx)
            curstring_toreplace = metabolites_unique{removeannidx(j)};
            curmet_toreplace = metabolites_IDs_multiple{i};
            %check if this is already in the list
            if ~isempty(metabolite_names_toreplace{1,1})
                curidx = find(cellfun(@(x) isempty(x), metabolite_names_toreplace(:,1)),1);
                curidx = find(ismember(metabolite_names_toreplace(1:curidx-1,1), curstring_toreplace));
                if ~isempty(curidx)
                    % update already updated string by another metabolite
                    curstring_toreplace = metabolite_names_toreplace{curidx,2};
                else
                    curidx = find(cellfun(@(x) isempty(x), metabolite_names_toreplace(:,1)),1);
                    metabolite_names_toreplace{curidx,1} = curstring_toreplace;
                end
            else
                curidx = find(cellfun(@(x) isempty(x), metabolite_names_toreplace(:,1)),1);
                metabolite_names_toreplace{curidx,1} = curstring_toreplace;
            end
            curstring_new = strrep(curstring_toreplace, [curmet_toreplace ';'], '');
            metabolite_names_toreplace{curidx,2} = curstring_new;
            sprintf('Replaceing %d %s\n', removeannidx(j), metabolites_unique{removeannidx(j)})
            sprintf('Old %s\n', curstring_toreplace)
            sprintf('New %s\n', curstring_new)
        end
    end
    % remove empty slots
    curidx = find(cellfun(@(x) isempty(x), metabolite_names_toreplace(:,1)),1);
    metabolite_names_toreplace(curidx:end,:) = [];
    for i=1:size(metabolite_names_toreplace,1)
       % sprintf('Found idx %d\n', nnz(cellfun(@(x) ismember(x, metabolite_names_toreplace(i,1)),...
       %                     CompoundID_updated)))
       % sprintf('Replacing %s\n',metabolite_names_toreplace{i,1})
       % sprintf('Replacing by %s\n',metabolite_names_toreplace{i,2})
        CompoundID_updated(cellfun(@(x) ismember(x, metabolite_names_toreplace(i,1)),...
                            CompoundID_updated)) = metabolite_names_toreplace(i,2);
    end
    % check if unique metabolites improved in the next loop iteration
    sprintf('Countdown from %d: replaced %d metabolites in multiple strings\n', maxiter, length(metabolites_IDs_multiple))
    maxiter = maxiter-1; % to avoid endless loops
end
if ~isempty(metabolites_unique_tosearch)
    sprintf('Used 10 iterations but not all annotation strings are reduced, %d are left\n', length(metabolites_IDs_multiple))
else
    sprintf('Reduced substrings counting down to %d iterations\n', maxiter)
end
% make a final list of unique metabolites
metabolites_unique = unique(CompoundID_updated);
metabolites_unique(cellfun(@(x) isempty(x), metabolites_unique))=[];
