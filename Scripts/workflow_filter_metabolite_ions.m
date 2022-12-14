%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze annotation table to score and filter multi-annotated compounds

% call script defining file dependenciesand global variables
addpath(genpath('.\'));
add_global_and_file_dependencies

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data requirements:
% 'kegg_agora_09_2020_hmdb_06_2021_table.csv'
% 'metabolites_allions_combined_formulas.csv'
% 'metabolites_allions_combined_norm_intensity.csv'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output: 
% Files:
% 'metabolites_allions_combined_formulas_with_metabolite_filters.csv'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load annotation from joint KEGG/AGORA/HMDB file
annKEGGAGORA = readtable([rawdataFolder 'kegg_agora_09_2020_hmdb_06_2021_table.csv']);

annotationTable = readtable([ resultsFolder...
                              'metabolites_allions_combined_formulas.csv']);
combinedIntensitiesTable = readtable([resultsFolder...
                              'metabolites_allions_combined_norm_intensity.csv']);
intensity_columns = cellfun(@(x) contains(x, '_M'), combinedIntensitiesTable.Properties.VariableNames);
combinedIntensitiesNorm = table2cell(combinedIntensitiesTable(:,intensity_columns));
combinedIntensitiesNorm  = cell2mat(combinedIntensitiesNorm );


Hmass = 1.007825;
% make a list of unique metabolites
metabolites_unique = unique(annotationTable.CompoundID);
metabolites_unique(cellfun(@(x) isempty(x), metabolites_unique))=[];
% for each compound, count how many ions fit
metabolites_unique_counts = zeros(length(metabolites_unique),1);
for i=1:length(metabolites_unique)
    metabolites_unique_counts(i) = nnz(ismember(annotationTable.CompoundID, ...
                                       metabolites_unique{i}));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find injection peaks
nbin = 50;
injection_peak_RT = zeros(3,1);
injection_peak_method = {'HILIC', 'C18', 'C08'};
for i=1:length(injection_peak_method)
    [hNum, hRT] = histcounts(annotationTable.RT(ismember(annotationTable.Method,...
                                              injection_peak_method{i})),...
                                              nbin);
    % find peaks
    hRT = hRT(1:end-1) + (hRT(2:end)-hRT(1:end-1))/2;
    [pks,locs] = findpeaks(hNum);
    % round injection peak
    injection_peak_RT(i) = round(hRT(locs(1)+1)*10)/10;
end

annotationFilter = zeros(size(annotationTable,1), 10);
annotationFilter_norank = zeros(size(annotationTable,1), 10);
corrThreshold = 0.75;
RTmaxthreshold = 1.5;
CisotopeDifferenceThreshold = 0.002;
CisotopeMassThreshold = 0.03;
% for each ion, calculate filters
for i=1:length(metabolites_unique)%1373%205%720%1160%514:534%1:100
    curidx = find(ismember(annotationTable.CompoundID, ...
                                      metabolites_unique{i}));
    curmets = annotationTable(curidx, :);
    %%%%%%%% Filter 1: NOT in injection peak
    annotationFilter(curidx,1) = curmets.RT > cellfun(@(x) injection_peak_RT(ismember(injection_peak_method,x)), curmets.Method);
    %%%%%%%% Filter 2: correlation with other methods
    nmet = length(curidx);
    if size(curmets,1)>1
        corrmat = zeros(nmet);
        for j=1:nmet
            for k=j:nmet
                nonzero = (combinedIntensitiesNorm(curidx(j),:)>5000) &...
                          (combinedIntensitiesNorm(curidx(k),:)>5000);                            
                if nnz(nonzero)
                    corrmat(j,k) = corr(reshape(log10(combinedIntensitiesNorm(curidx(k),nonzero)),[],1),...
                                        reshape(log10(combinedIntensitiesNorm(curidx(j),nonzero)),[],1));
                end
            end
        end
        sumcorr = triu(corrmat,1);
        sumcorr = sum(sumcorr>corrThreshold) + sum(sumcorr'>corrThreshold);
        annotationFilter(curidx,2) = sumcorr;    
    end
    %%%%%%%% Filter 3: correlation with opposite column injection peak
    for j=1:size(curmets,1)
        if curmets.Method{j}(1)=='H'
            specific_corr = (ismember(cellfun(@(x) x(1), curmets.Method), 'C')) &...
                            (curmets.RT<cellfun(@(x) injection_peak_RT(ismember(injection_peak_method,x)), curmets.Method));
        else
            specific_corr = (ismember(cellfun(@(x) x(1), curmets.Method), 'H')) &...
                            (curmets.RT<cellfun(@(x) injection_peak_RT(ismember(injection_peak_method,x)), curmets.Method));
        end
        if nnz(specific_corr)
            specific_corr = corrmat(specific_corr,specific_corr); 
            sumcorr = triu(specific_corr,1);
            annotationFilter(curidx(j),3) = sum(sum(sumcorr>corrThreshold));
        end
    end
    %%%%%%%% Filter 4, 5, 6: ranking of max, mean and median intensities
    [~, sortedidx] = sort(curmets.maxIntensity, 'ascend');
    [~,~,annotationFilter(curidx,4)] = intersect(1:length(sortedidx),...
                                                 sortedidx);
    [~, sortedidx] = sort(curmets.meanintensity, 'ascend');
    [~,~,annotationFilter(curidx,5)] = intersect(1:length(sortedidx),...
                                                 sortedidx);
                                   
    [~, sortedidx] = sort(curmets.medianIntensity,'ascend');
    [~,~,annotationFilter(curidx,6)] = intersect(1:length(sortedidx),...
                                                 sortedidx);
    %%%%%%%% Filter 4, 5, 6: record max, mean and median intensities without ranking 
    annotationFilter_norank(curidx,4) = curmets.maxIntensity;
    annotationFilter_norank(curidx,5) = curmets.meanintensity;
    annotationFilter_norank(curidx,6) = curmets.medianIntensity;
    
    %%%%%%%% Filter 7: correlation with opposite ionization mode
    for j=1:size(curmets,1)
        opposite_mode_idx = ismember(curmets.Method, curmets.Method{j}) &...
                            (curmets.Mode == -curmets.Mode(j));
        [minRTdiff,opposite_mode_idx] = min(abs(curmets.RT(opposite_mode_idx)-...
                                        curmets.RT(j)));
                                            
        if minRTdiff < RTmaxthreshold
            annotationFilter(curidx(j),7) = (corrmat(j, opposite_mode_idx)>corrThreshold) |...
                                            (corrmat(opposite_mode_idx,j)>corrThreshold);
        end
    end
    %%%%%%%% Filter 8: ranking of number of samples
    [~, sortedidx] = sort(curmets.NumDetectedSamples,'ascend');
    [~,~,annotationFilter(curidx,8)] = intersect(1:length(sortedidx),...
                                                 sortedidx);
    %%%%%%%% Filter 8: number of samples without ranking 
    annotationFilter_norank(curidx,8) = curmets.NumDetectedSamples;
   
    %%%%%%%% Filter 9: number of carbon isotopes
    % calculate spectral ratios
    curformula = curmets.CompoundFormula{1};
    if curformula(1)=='C'
        curcarbon = '';
        j=2;
        while ~isnan(str2double(curformula(j)))
            curcarbon = [curcarbon curformula(j)];
            j=j+1;
        end
        curcarbon = str2double(curcarbon);
        if isnan(curcarbon)
            curcarbon = 1;
        end
        if curcarbon>1
            theoreticalC = [((1-0.01109)^(curcarbon));...
                            ((0.01109)*(1-0.01109)^(curcarbon-1))*nchoosek(curcarbon,1);...
                            ((0.01109^2)*(1-0.01109)^(curcarbon-2))*nchoosek(curcarbon,2);...
                            ];
        else
            theoreticalC = [((1-0.01109)^(curcarbon));...
                            ((0.01109)*(1-0.01109)^(curcarbon-1))*nchoosek(curcarbon,1);...
                            0;...
                            ];
        end
        theoreticalC = theoreticalC/theoreticalC(1);
        theoreticalCmass = [0 1.003355 2*1.003355];

        experimentalC = zeros(3,size(curmets,1));
        experimentalCmass = zeros(3,size(curmets,1));
        for j=1:size(curmets,1)
            curspectrum = curmets.Spectrum{j};
            curspectrum = strsplit(curspectrum, {'(', ','});
            curspectrum = curspectrum(2:end);
            curspectrum = cellfun(@(x) strrep(x, ')', ''), curspectrum, 'unif', 0);
            curspectrum = cellfun(@(x) str2double(x), curspectrum);

            curexpCmass = curspectrum(1:2:end)-curmets.MZ(j)-curmets.Mode(j)*Hmass;
            [~, curexpCmassIDX] = arrayfun(@(x) min(abs(curexpCmass-x)), theoreticalCmass);

            curexpCmass = curexpCmass(curexpCmassIDX);
            curexpC = curspectrum(2:2:end);
            curexpC = curexpC(curexpCmassIDX);
            curexpC = curexpC/curexpC(1);
            % check that the candidate is close to the theoretical
            correctmass = abs(curexpCmass-theoreticalCmass)<0.5;
            experimentalC(correctmass,j) = curexpC(correctmass);
            experimentalCmass(correctmass,j) = curexpCmass(correctmass);
        end   
        annotationFilter(curidx,9) = sum(abs(bsxfun(@minus,experimentalCmass,theoreticalCmass'))<CisotopeDifferenceThreshold)-1;
        annotationFilter(curidx,10) = sum(abs(bsxfun(@minus,experimentalC,theoreticalC))<CisotopeMassThreshold)-1;
    end
    annotationFilter_norank(curidx,[1:3 7 9 10]) =  annotationFilter(curidx,[1:3 7 9 10]);
    
    annotationFilter(curidx,:) = bsxfun(@rdivide,annotationFilter(curidx,:), max(annotationFilter(curidx,:)));
end
annotationFilter(isnan(annotationFilter))=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for each metabolite, leave only one annotation depending on ion ranking
metabolites_IDs = cell(length(metabolites_unique)*5,1);
idx = 1;
for i=1:length(metabolites_unique)
    curmets = strsplit(metabolites_unique{i}, ';');
    metabolites_IDs(idx:idx+length(curmets)-1) = curmets;
    idx = idx+length(curmets);
end
metabolites_IDs(idx:end) = [];
metabolites_IDs_unique = unique(metabolites_IDs);
metabolites_IDs_unique(cellfun(@(x) isempty(strtrim(x)), metabolites_IDs_unique))=[];

% for each metabolite, keep only one annotation
filteredCompoundID = cell(size(annotationTable,1),1);
filteredCompoundName = cell(size(annotationTable,1),1);
filteredCompoundFormula = cell(size(annotationTable,1),1);
filteredCompoundMZdelta = cell(size(annotationTable,1),1);
filteredCompoundIDX = cell(size(annotationTable,1),1);
filteredCompoundFilterValue = cell(size(annotationTable,1),1);
filteredCompoundAlternativeIDX = cell(size(annotationTable,1),1);
filteredCompoundAlternativeFilters = cell(size(annotationTable,1),1);
for i=1:size(annotationTable,1)
    filteredCompoundID{i} = ' ';
    filteredCompoundName{i} = ' ';
    filteredCompoundFormula{i} = ' ';
    filteredCompoundMZdelta{i} = ' ';
    filteredCompoundIDX{i} = ' ';
    filteredCompoundFilterValue{i} = ' ';
    filteredCompoundAlternativeIDX{i} = ' ';
    filteredCompoundAlternativeFilters{i} = ' ';
end
% filter metabolites
for i=1:length(metabolites_IDs_unique)
    curionidx = find(cellfun(@(x) contains(x, metabolites_IDs_unique{i}),...
                                         annotationTable.CompoundID));
    curfilter = annotationFilter_norank(curionidx,:);
    % transform max, mean and median intensity to log scale before ranking
    % to avoid very large differences between ions
    curfilter(:, 4:6) = log10(curfilter(:,4:6));
    curfilter = bsxfun(@rdivide,curfilter, max(curfilter));
    % remove nans resulted from division by 0
    curfilter(isnan(curfilter)) = 0;
    % give more weight to the number of sample filter
    curfilter(:, 8) = 2*curfilter(:, 8);
    % calculate final sum of filters
    curfilter_sum = sum(curfilter,2);
    
    filteredIDX = (curfilter_sum == max(curfilter_sum));
    % if there is more than one filtered IDX, select one with more samples
    if nnz(filteredIDX) > 1
        filteredIDX = ((curfilter_sum == max(curfilter_sum)) &...
                                 curfilter(:, 8) == max(curfilter(:, 8)));
    end
    % get current compound info from corresponding annotation row
    curAnnIDX = ismember(strsplit(annotationTable.CompoundID{curionidx(filteredIDX)},';',...
                         'CollapseDelimiters',0),...
                         metabolites_IDs_unique{i});
    curAnnName = strsplit(annotationTable.CompoundName{curionidx(filteredIDX)},';',...
                         'CollapseDelimiters',0);
    curAnnFormula = strsplit(annotationTable.CompoundFormula{curionidx(filteredIDX)},';',...
                         'CollapseDelimiters',0);
    curAnnMZdelta = strsplit(annotationTable.MZdelta{curionidx(filteredIDX)},';',...
                         'CollapseDelimiters',0);
    curAnnCompoundIDX = strsplit(annotationTable.IDX{curionidx(filteredIDX)},';',...
                         'CollapseDelimiters',0);                
    % add annotation to filtered annotation vectors
    filteredCompoundID{curionidx(filteredIDX)} = ...
        strjoin([filteredCompoundID(curionidx(filteredIDX)),...
                metabolites_IDs_unique(i)], ';');
    filteredCompoundName{curionidx(filteredIDX)} = ...
        strjoin([filteredCompoundName(curionidx(filteredIDX)) ,...
                curAnnName(curAnnIDX)],';');
    filteredCompoundFormula{curionidx(filteredIDX)} = ...
        strjoin([filteredCompoundFormula(curionidx(filteredIDX)),...
                curAnnFormula(curAnnIDX)],';');
    filteredCompoundMZdelta{curionidx(filteredIDX)} = ...
        strjoin([filteredCompoundMZdelta(curionidx(filteredIDX)),...
                curAnnMZdelta(curAnnIDX)],';');
    filteredCompoundIDX{curionidx(filteredIDX)} = ...
        strjoin([filteredCompoundIDX(curionidx(filteredIDX)),...
                curAnnCompoundIDX(curAnnIDX)],';');
    % filter values
    filteredCompoundFilterValue{curionidx(filteredIDX)} = ...
        strjoin([filteredCompoundFilterValue(curionidx(filteredIDX)),...
                num2str(curfilter_sum(filteredIDX))], ';');
    filteredCompoundAlternativeIDX{curionidx(filteredIDX)} = ...
        strjoin([filteredCompoundAlternativeIDX(curionidx(filteredIDX)),...
        ['(', strjoin(arrayfun(@(x) num2str(x), curionidx(filteredIDX==0), 'unif', 0),';'),')']], ';');
    filteredCompoundAlternativeFilters{curionidx(filteredIDX)} = ...
        strjoin([filteredCompoundAlternativeFilters(curionidx(filteredIDX)),...
        ['(', strjoin(arrayfun(@(x) num2str(x), curfilter_sum(filteredIDX==0), 'unif', 0),';'),')']], ';');
    if mod(i,100)==0
        fprintf('Done with %d of %d metabolites \n', i, length(metabolites_IDs_unique))
    end
end                                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check whether any multiple annotations or empty annotations are left
% make a list of unique metabolites
metabolites_unique = unique(filteredCompoundID);
metabolites_unique(cellfun(@(x) isempty(x), metabolites_unique))=[];
metabolites_unique(cellfun(@(x) isequal(x,' '), metabolites_unique))=[];
% for each compound, count how many ions fit
metabolites_unique_counts = zeros(length(metabolites_unique),1);
for i=1:length(metabolites_unique)
    metabolites_unique_counts(i) = nnz(ismember(filteredCompoundID, ...
                                       metabolites_unique{i}));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for each metabolite, leave only one annotation depending on ion ranking
metabolites_IDs = cell(length(metabolites_unique)*5,1);
idx = 1;
for i=1:length(metabolites_unique)
    curmets = unique(strsplit(metabolites_unique{i}, ';'));
    metabolites_IDs(idx:idx+length(curmets)-1) = curmets;
    idx = idx+length(curmets);
end
metabolites_IDs(idx:end) = [];
metabolites_IDs_unique = unique(metabolites_IDs);
metabolites_IDs_unique(cellfun(@(x) isempty(strtrim(x)), metabolites_IDs_unique))=[];


metabolites_IDs_unique_count = cellfun(@(x) nnz(ismember(metabolites_IDs,x)),...
                                            metabolites_IDs_unique);

% check mz differences between ions that share some metabolites, but not all
similar_idx = find(metabolites_IDs_unique_count>1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate metabolite filter as non-empty ions
metaboliteFilter = cellfun(@(x) ~isempty(strtrim(x)), filteredCompoundID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get name from annotation table instead
for i=1:length(filteredCompoundID)
    curid = filteredCompoundID{i};
    if ~isempty(curid)
        curid = strsplit(curid,';');
        curid(cellfun(@(x) isempty(strtrim(x)), curid))=[];
        curnames = annKEGGAGORA.CompoundNAME(cellfun(@(x) ...
                                    find(ismember(annKEGGAGORA.KEGGID,x)),curid));
        curnames = cellfun(@(x) ['(' x ')'], curnames, 'unif',0);                        
        filteredCompoundName{i} = strjoin(curnames, ';');
    end
end
% remove first ; from filtered values
filteredCompoundID = cellfun(@(x) strrep(x,' ;',''), filteredCompoundID, 'unif',0);
filteredCompoundFormula = cellfun(@(x) strrep(x,' ;',''), filteredCompoundFormula, 'unif',0);
filteredCompoundMZdelta = cellfun(@(x) strrep(x,' ;',''), filteredCompoundMZdelta, 'unif',0);
filteredCompoundIDX = cellfun(@(x) strrep(x,' ;',''), filteredCompoundIDX, 'unif',0);
filteredCompoundFilterValue = cellfun(@(x) strrep(x,' ;',''), filteredCompoundFilterValue, 'unif',0);
filteredCompoundAlternativeFilters = cellfun(@(x) strrep(x,' ;',''), filteredCompoundAlternativeFilters, 'unif',0);
filteredCompoundAlternativeIDX = cellfun(@(x) strrep(x,' ;',''), filteredCompoundAlternativeIDX, 'unif',0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% manual correction of selected compounds
% propionate
targetMZ = 74.037; % propionate
targetRT = 0.431;
targetMethod = 'C18';
% get indeces of these masses
selectedIDX = find(abs(annotationTable.MZ-targetMZ)<=massThreshold);
selectedRT = annotationTable.RT(selectedIDX);
selectedMethod = annotationTable.Method(selectedIDX);
selectedFilter = metaboliteFilter(selectedIDX);
% get the index of manually selected compound
manualIDX = selectedIDX((abs(selectedRT-targetRT)<0.1) &...
                        ismember(selectedMethod, targetMethod));
filteredIDX = selectedIDX(selectedFilter==1);
if filteredIDX ~= manualIDX
    % swap annotation of current index and manual index
    % filter
    temp = metaboliteFilter(filteredIDX);
    metaboliteFilter(filteredIDX) = metaboliteFilter(manualIDX);
    metaboliteFilter(manualIDX) = temp;
    % compoundID
    temp = filteredCompoundID(filteredIDX);
    filteredCompoundID(filteredIDX) = filteredCompoundID(manualIDX);
    filteredCompoundID(manualIDX) = temp;
    % compoundName    
    temp = filteredCompoundName(filteredIDX);
    filteredCompoundName(filteredIDX) = filteredCompoundName(manualIDX);
    filteredCompoundName(manualIDX) = temp;
    % compoundFormula    
    temp = filteredCompoundFormula(filteredIDX);
    filteredCompoundFormula(filteredIDX) = filteredCompoundFormula(manualIDX);
    filteredCompoundFormula(manualIDX) = temp;
    % compoundMZdelta    
    temp = filteredCompoundMZdelta(filteredIDX);
    filteredCompoundMZdelta(filteredIDX) = filteredCompoundMZdelta(manualIDX);
    filteredCompoundMZdelta(manualIDX) = temp;
    % compoundIDX    
    temp = filteredCompoundIDX(filteredIDX);
    filteredCompoundIDX(filteredIDX) = filteredCompoundIDX(manualIDX);
    filteredCompoundIDX(manualIDX) = temp;
    % compoundFilterValue    
    temp = filteredCompoundFilterValue(filteredIDX);
    filteredCompoundFilterValue(filteredIDX) = filteredCompoundFilterValue(manualIDX);
    filteredCompoundFilterValue(manualIDX) = temp;
    % compoundAlternativeFilters    
    temp = filteredCompoundAlternativeFilters(filteredIDX);
    filteredCompoundAlternativeFilters(filteredIDX) = filteredCompoundAlternativeFilters(manualIDX);
    filteredCompoundAlternativeFilters(manualIDX) = temp;
    % compoundAlternativeIDX    
    temp = filteredCompoundAlternativeIDX(filteredIDX);
    filteredCompoundAlternativeIDX(filteredIDX) = filteredCompoundAlternativeIDX(manualIDX);
    filteredCompoundAlternativeIDX(manualIDX) = temp;
end                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add metabolite annotation filter to the columns
annotationTable = [ annotationTable,...
                    array2table([filteredCompoundID,...
                                 filteredCompoundName,...
                                 filteredCompoundFormula,...
                                 filteredCompoundMZdelta,...
                                 filteredCompoundIDX,...
                                 filteredCompoundFilterValue,...
                                 filteredCompoundAlternativeFilters,...
                                 filteredCompoundAlternativeIDX],...
                                 'VariableNames',...
                                            {'FilteredCompoundID',...
                                             'FilteredCompoundName',...
                                             'FilteredCompoundFormula',...
                                             'FilteredCompoundMZdelta',...
                                             'FilteredCompoundIDX',...
                                             'FilteredCompoundFilterValue',...
                                             'FilteredCompoundAlternativeFilters',...
                                             'FilteredCompoundAlternativeIDX'})];   
annotationTable = [ annotationTable,...
                    array2table(metaboliteFilter,...
                                 'VariableNames',{'MetaboliteFilter'})];
% add index to track alternative annotations
annotationTable.Index = [1:length(filteredCompoundFilterValue)]';
% write table to file
write(annotationTable, ...
    [outputFolder, 'metabolites_allions_combined_formulas_with_metabolite_filters.csv']);

% remove ; in the beginning of metabolites_unique
metabolites_unique_clean = metabolites_unique;
for i=1:length(metabolites_unique_clean)
    if isequal(metabolites_unique_clean{i}(1:2), ' ;')
        metabolites_unique_clean{i} = metabolites_unique_clean{i}(3:end);
    end
end    
