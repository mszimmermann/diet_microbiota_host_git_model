%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% merge untargeted metabolomics data from different tissues together

% call script defining file dependenciesand global variables
addpath(genpath('.\'));
add_global_and_file_dependencies

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data requirements:
% 'metabolites_allions_combined_formulas_with_metabolite_filters.csv'
% 'metabolites_allions_combined_norm_intensity.csv'
% 'hmdbPTWtables.mat'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output: 
% Files:
% 'table_diff_abundance_metabolite_ions_removed2outliers.csv'
% or 'table_diff_abundance_metabolite_ions.csv'
% Figures:
% 'fig_sup_volcanos_combined_ann_metabolites_removed2outliers.pdf'
% or 'fig_sup_volcanos_combined_ann_metabolites.pdf'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load annotation from joint KEGG/AGORA/HMDB file
annKEGGAGORA = readtable([rawdataFolder 'kegg_agora_09_2020_hmdb_06_2021_table.csv']);

% Import experimental data
filenameP1 = {[rawdataFolder 'metabolomics_data\C18_Neg_Alignments\C18-Neg_'],...
              [rawdataFolder 'metabolomics_data\C18_Pos_Alignments\C18-Pos_'],...
              [rawdataFolder 'metabolomics_data\C08_Neg_Alignments\C08_Neg_'],...
              [rawdataFolder 'metabolomics_data\C08_Pos_Alignments\C08_Pos_'],...
              [rawdataFolder 'metabolomics_data\HILIC_Neg_Alignments\HILIC-Neg_'],...
              [rawdataFolder 'metabolomics_data\HILIC_Pos_Alignments\HILIC-Pos_'],...
              };
filenameP3 = [repmat({'_DIL05.txt'},1,6), {'.txt'}, {'_DIL02.txt'}];
fileTissues = {'SI1' 'SI2' 'SI3' 'Cecum' 'Colon' 'DC' 'Serum' 'Liver'};

tissueWeights = readtable([rawdataFolder 'metabolomics_data\tissue_weights.txt'],...
                            'delim', '\t');
tissueWeights = table2cell(tissueWeights);

curCompounds_cell = cell(length(filenameP1)*length(fileTissues),1);
curCompoundsSpectrum_cell = cell(length(filenameP1)*length(fileTissues),1);
curIntensitiesRaw_cell = cell(length(filenameP1)*length(fileTissues),1);
curSampleNames_cell = cell(length(filenameP1)*length(fileTissues),1);

idx = 1;
% read compounds and filter before merging
for fi = 1:length(filenameP1)
    for d = 1:length(fileTissues)
        datafile = [filenameP1{fi} fileTissues{d} filenameP3{d}];
        if ~exist(datafile, 'file')
            % maybe it is the serum file with DIL02 extension
            datafile = strrep(datafile, '.txt', '_DIL02.txt');
            if ~exist(datafile, 'file')
                fprintf('File %s does not exist!\n', datafile);
                idx = idx+1;
                continue;
            end
        end
     % Import data of one experimental condition
        TempStruct = ReadMixedTxtUntargeted(datafile, '\t');
        TempStruct.SampleNames = TempStruct.SampleNames'; 

        sampleType = cell(size(TempStruct.SampleNames));
        sampleDiet = cell(size(TempStruct.SampleNames));
        % get sample type and diet from sample name
        for i = 1:length(sampleType)
            curname = strsplit(TempStruct.SampleNames{i}, '_');
            sampleType{i} = curname{4};
            sampleDiet{i} = curname{2};
        end

        curIntensitiesRaw =  TempStruct.IntensitiesRaw;
        curCompounds =  TempStruct.Compounds;
        curCompoundsSpectrum =  TempStruct.CompositeSpectrum;

        %check if intensities are likely on logarithmic scale
        if nnz(curIntensitiesRaw==0)>0 && nnz(curIntensitiesRaw>100000)<1000
            curIntensitiesRaw = 10.^(curIntensitiesRaw);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        % get sample weights
        [~, curWeights] = intersect(tissueWeights(:,1),...
                              cellfun(@(x) strrep(x,'(raw)',''), TempStruct.SampleNames, 'unif', 0),...
                              'stable');
        if isempty(curWeights)
            [~, curWeights] = intersect(tissueWeights(:,1),...
                              cellfun(@(x) strrep(x,'-Neg-DIL05(raw)',''), TempStruct.SampleNames, 'unif', 0),...
                              'stable');
        end
        if isempty(curWeights)
            [~, curWeights] = intersect(tissueWeights(:,1),...
                              cellfun(@(x) strrep(x,'-Pos(raw)',''), TempStruct.SampleNames, 'unif', 0),...
                              'stable');
        end
        if isempty(curWeights)
            [~, curWeights] = intersect(tissueWeights(:,1),...
                              cellfun(@(x) strrep(x,'-Pos-DIL05(raw)',''), TempStruct.SampleNames, 'unif', 0),...
                              'stable');
        end
        if isempty(curWeights)
            [~, curWeights] = intersect(tissueWeights(:,1),...
                              cellfun(@(x) strrep(x,'-HILIC-Neg-DIL05(raw)',''), TempStruct.SampleNames, 'unif', 0),...
                              'stable');
        end
        if isempty(curWeights)
            [~, curWeights] = intersect(tissueWeights(:,1),...
                              cellfun(@(x) strrep(x,'-HILIC-Neg(raw)',''), TempStruct.SampleNames, 'unif', 0),...
                              'stable');
        end
        if isempty(curWeights)
            [~, curWeights] = intersect(tissueWeights(:,1),...
                              cellfun(@(x) strrep(x,'-HILIC-Pos-DIL05_1(raw)',''), TempStruct.SampleNames, 'unif', 0),...
                              'stable');
        end
        if isempty(curWeights)
            [~, curWeights] = intersect(tissueWeights(:,1),...
                              cellfun(@(x) strrep(x,'-HILIC-Pos(raw)',''), TempStruct.SampleNames, 'unif', 0),...
                              'stable');
        end
        if isempty(curWeights)
            [~, curWeights] = intersect(tissueWeights(:,1),...
                              cellfun(@(x) strrep(x,'-C08-Neg-DIL05(raw)',''), TempStruct.SampleNames, 'unif', 0),...
                              'stable');
        end
        if isempty(curWeights)
            [~, curWeights] = intersect(tissueWeights(:,1),...
                              cellfun(@(x) strrep(x,'-C08-Neg(raw)',''), TempStruct.SampleNames, 'unif', 0),...
                              'stable');
        end
        if isempty(curWeights)
            [~, curWeights] = intersect(tissueWeights(:,1),...
                              cellfun(@(x) strrep(x,'-C08-Pos-DIL05(raw)',''), TempStruct.SampleNames, 'unif', 0),...
                              'stable');
        end
        if isempty(curWeights)
            [~, curWeights] = intersect(tissueWeights(:,1),...
                              cellfun(@(x) strrep(x,'-C08-Pos(raw)',''), TempStruct.SampleNames, 'unif', 0),...
                              'stable');
        end
        if isempty(curWeights)
           fprintf('Could not find weights for iteration %d\n', idx);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % find outliers due to bad injections, if sum of all compounds is
        % lower than the mean of all 48 data sets minus three times their std
        % threshold value
        Threshold = mean(sum(curIntensitiesRaw)-(3*(std(sum(curIntensitiesRaw)))));
        curExcludedDatasets = (sum(curIntensitiesRaw) < Threshold) |...
                              cell2mat(tissueWeights(curWeights,2))'<10; % WEIGHT THRESHOLD
        curIntensitiesNorm = curIntensitiesRaw;
        curIntensitiesNorm(:,curExcludedDatasets) = [];
        curWeights(curExcludedDatasets) = [];
        curSampleNames = TempStruct.SampleNames;
        curSampleNames(curExcludedDatasets) = [];
        
        % exclude ions that were unique to the excluded dataset in all
        % parameters
        curIdxExcluded = find(sum(curIntensitiesNorm,2) ==0);
        curIntensitiesNorm(curIdxExcluded,:) =  [];
        curCompounds(curIdxExcluded,:) =  [];
        curCompoundsSpectrum(curIdxExcluded,:) =  [];

        % Only leave compounds that are present in >=half of samples!
        numabsentThreshold = nnz(~ismember(sampleType, 'CVR'))/2;
        numcompounds = sum(curIntensitiesNorm(:, ~ismember(sampleType, 'CVR'))==1,2);

        curCompounds = curCompounds(numcompounds<=numabsentThreshold);
        curCompoundsSpectrum = curCompoundsSpectrum(numcompounds<=numabsentThreshold);
        curIntensitiesNorm = curIntensitiesNorm(numcompounds<=numabsentThreshold,:);

        
        % normalize intensity in each sample by sampe weight
        % set ones to zeros 
        curIntensitiesNorm(curIntensitiesNorm==1)=0;
        for i=1:length(curWeights)
            curIntensitiesNorm(:,i) = curIntensitiesNorm(:,i)*1000/tissueWeights{curWeights(i),2};
        end
        
        curCompounds_cell{idx} = curCompounds;
        curCompoundsSpectrum_cell{idx} = curCompoundsSpectrum;
        curIntensitiesRaw_cell{idx} = curIntensitiesNorm;
        curSampleNames_cell{idx} = curSampleNames;
        fprintf('Filtered compounds for %s %s\n', filenameP1{fi}, fileTissues{d});
        idx = idx+1;
    end
end
clear datafile curIntensitiesRaw curCompounds curSamples curCompoundsSpectrum
clear Threshold curExcludedDatasets curIntensitiesNorm curIdxExcluded curname
clear totTime totDataPools totDataMat numcompounds TempStruct
% remove HILIC DC because there is no samples
empty_data = cellfun(@(x) isempty(x), curCompounds_cell);
curCompounds_cell(empty_data) = [];
curCompoundsSpectrum_cell(empty_data) = [];
curIntensitiesRaw_cell(empty_data) = [];
curSampleNames_cell(empty_data) = [];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Plot box plot of metabolite intensities across samples
% col=@(x)reshape(x,numel(x),1);
% boxplot2=@(C,varargin)boxplot(log10(cell2mat(cellfun(col,col(C),'uni',0))),...
%     cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:});
% boxplot3=@(C,varargin)boxplot(log10(cell2mat(cellfun(col,col(C),'uni',0))),varargin{1});
% test = curIntensitiesRaw_cell;
% for i=1:length(test)
%     test{i}(test{i}==0) = nan;
%     test{i}(~isnan(test{i})) = i;
% end
% test = cell2mat(cellfun(col,col(test),'uni',0));
% boxplot3(curIntensitiesRaw_cell, test);

% idx = 1;
% % read compounds and filter before merging
% fileNames = cell(size(curIntensitiesRaw_cell));
% for fi = 1:length(filenameP1)
%     for d = 1:length(fileTissues)
%         datafile = strsplit(filenameP1{fi}, '\');
%         fileNames{idx} = [datafile{end} fileTissues{d}];
%         idx = idx+1;
%     end
% end
% set(gca, 'XTick', 1:length(fileNames));
% set(gca, 'XTickLabel', (fileNames));
% xticklabel_rotate([], 90)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get mass spec method from sample names
curMethod = cell(size(curSampleNames_cell));
for i=1:length(curSampleNames_cell)
    curMethod{i} = curSampleNames_cell{i}{1};
end
curMode = -ones(size(curMethod));
curMode(cellfun(@(x) contains(x, 'Pos'), curMethod))=1;
curMethod(cellfun(@(x) contains(x, 'HILIC'), curMethod))={'HILIC'};
curMethod(cellfun(@(x) contains(x, 'C08'), curMethod))={'C08'};
curMethod(cellfun(@(x) ~(contains(x, 'HILIC') |...
                         contains(x, 'C08')), curMethod))={'C18'};
curMethod_unique = unique(curMethod);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a matrix with metabolites across tissues
% first merge metabolites across tissues
combinedDataNorm_cell = cell(length(curMethod_unique),1);
combinedDataTissues_cell = cell(length(curMethod_unique),1);
combinedDataDiet_cell = cell(length(curMethod_unique),1);
combinedDataMouse_cell = cell(length(curMethod_unique),1);
combinedDataionMode_cell = cell(length(curMethod_unique),1);
combinedIonMZ_cell = cell(length(curMethod_unique),1);
combinedIonRT_cell = cell(length(curMethod_unique),1);
combinedIonSpectrum_cell = cell(length(curMethod_unique),1);
for method_i = 1:length(curMethod_unique)
    changingMets = cell(100000,1);
    changingMetsSpectrum = cell(100000,1);
    changingMetsIonMode = zeros(100000,1); %ionization mode for comparison of spectra
    curmethodidx = find(ismember(curMethod, curMethod_unique{method_i}));
    
    idx = 1;
    for i=1:length(curmethodidx)
        curMets = curCompounds_cell{curmethodidx(i)};
        curMetsSpectrum = curCompoundsSpectrum_cell{curmethodidx(i)};

        changingMets(idx:idx+length(curMets)-1) = curMets;
        changingMetsSpectrum(idx:idx+length(curMets)-1) = curMetsSpectrum;
        changingMetsIonMode(idx:idx+length(curMets)-1) = curMode(curmethodidx(i))*ones(size(curMets));
        idx = idx+length(curMets);
    end
    changingMets(idx:end) = [];
    changingMetsSpectrum(idx:end) = [];
    changingMetsIonMode(idx:end) = [];
    clear curMets curMetsSpectrum curChange metFC metPadj
    [changingMets,idx] = unique(changingMets);
    changingMetsSpectrum = changingMetsSpectrum(idx);
    changingMetsIonMode = changingMetsIonMode(idx);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % merge metabolites that are close in mass and RT
    [changingMets_merged,...
              changingMets_merged_idx,...
              changingMets_merged_idx_unique,...
              changingMets_merged_spectrum,...
              changingMets_merged_mass,...
              changingMets_merged_RT,...
              changingMets_merged_mass_delta,...
              changingMets_merged_RT_delta,...
              changingMets_merged_number,...
              changingMets_merged_mode] = merge_changing_metabolites(changingMets,...
                                                                       changingMetsSpectrum,...
                                                                       changingMetsIonMode);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get idx conversion from each ion file to the merged
    curCompoundsMergedIDXconversion = zeros(length(changingMets_merged),...
                                            length(curmethodidx));
    for d = 1:length(curmethodidx)
        curMZ = cellfun(@(x) str2double(x(1:strfind(x, '@')-1)), curCompounds_cell{curmethodidx(d)});
        curRT = cellfun(@(x) str2double(x(strfind(x, '@')+1:end)), curCompounds_cell{curmethodidx(d)});

        tic
        for i=1:size(changingMets_merged_mass,1)

            curMetMass = changingMets_merged_mass(i);
            curMetRT = changingMets_merged_RT(i);
            % add ion mode into merging
            
            if changingMets_merged_mode(i) == curMode(curmethodidx(d))

                metIdx = find( (abs(curMZ-curMetMass)<=max(massThreshold,...
                                                           changingMets_merged_mass_delta(i))) &...
                                (abs(curRT-curMetRT)<=max(RTthreshold,...
                                                          changingMets_merged_RT_delta(i))));

                if length(metIdx)>1
                    rtDiff = abs(curRT(metIdx)-curMetRT);
                    metIdx = metIdx(rtDiff == min(rtDiff));
                    if length(metIdx)>1
                        mzDiff = abs(curMZ(metIdx)-curMetMass);
                        metIdx = metIdx(mzDiff == min(mzDiff));
                    end
                end
                if ~isempty(metIdx)    
                    curCompoundsMergedIDXconversion(i,d) = metIdx;
                end
            end
        end
        toc
        fprintf('Merged metabolites for tissue number %d\n', d);
    end        

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate intensity matrix for all metabolites
    combinedIntensity = zeros(size(curCompoundsMergedIDXconversion,1),...
                                30*length(curmethodidx));
    combinedTissues = cell(30*length(curmethodidx),1);
    combinedDiet = cell(30*length(curmethodidx),1);
    combinedType = cell(30*length(curmethodidx),1);
    combinedMouse = cell(30*length(curmethodidx),1);
    combinedionMode = zeros(30*length(curmethodidx),1);
    idx = 1;
    for d = 1:length(curmethodidx)
        curIntensitiesRaw =  curIntensitiesRaw_cell{curmethodidx(d)};

        sampleType = cell(size(curSampleNames_cell{curmethodidx(d)}));
        sampleTissue = cell(size(curSampleNames_cell{curmethodidx(d)}));
        sampleDiet = cell(size(curSampleNames_cell{curmethodidx(d)}));
        sampleMouse = cell(size(curSampleNames_cell{curmethodidx(d)}));
        % get sample type and diet from sample name
        for i = 1:length(sampleType)
            curname = strsplit(curSampleNames_cell{curmethodidx(d)}{i}, '_');
            sampleType{i} = curname{4};
            sampleTissue{i} = curname{3};
            sampleDiet{i} = curname{2};
            sampleMouse{i} = curname{1};
        end

        combinedIntensity(curCompoundsMergedIDXconversion(:,d)~=0, idx:idx+size(sampleType)-1) = ...
             curIntensitiesRaw(curCompoundsMergedIDXconversion(curCompoundsMergedIDXconversion(:,d)~=0,d),:);

        combinedTissues(idx:idx+size(sampleType)-1) = sampleTissue;
        combinedDiet(idx:idx+size(sampleType)-1) = sampleDiet;
        combinedType(idx:idx+size(sampleType)-1) = sampleType;
        combinedMouse(idx:idx+size(sampleType)-1) = sampleMouse;
        combinedionMode(idx:idx+size(sampleType)-1) = curMode(curmethodidx(d))*ones(size(sampleType));
        idx=idx+length(sampleType);
        fprintf('Merged selected metabolites for %s \n', sampleTissue{1});

    end
    combinedIntensity(:, idx:end) = [];
    combinedTissues(idx:end) = [];
    combinedDiet(idx:end) = [];
    combinedType(idx:end) = [];
    combinedMouse(idx:end) = [];
    combinedionMode(idx:end) = [];
    clear datafile curIntensitiesRaw curCompounds curSamples curCompoundsSpectrum
    clear Threshold curExcludedDatasets curIntensitiesNorm curIdxExcluded 
    clear TempResultsMetabolites TempStruct

    % calculate mean intensity across ions for POS and NEG samples
    if length(unique(combinedionMode))>1
        meanPosNeg = [mean(combinedIntensity(:,combinedionMode==-1),2)...
                  mean(combinedIntensity(:, combinedionMode==1),2)];      
        % for each ion, leave intensities of the method with higher intensity
        currentIonMode = -1 + 2*(meanPosNeg(:,2) > meanPosNeg(:,1));
        for i=1:size(meanPosNeg,1)
            if meanPosNeg(i,2) > meanPosNeg(i,1)
                for j=1:length(fileTissues)
                      if nnz(ismember(combinedTissues, fileTissues{j}) & ...
                                           combinedionMode==1)
                          combinedIntensity(i,ismember(combinedTissues, fileTissues{j}) &...
                                        combinedionMode==-1) = ...
                          combinedIntensity(i, ismember(combinedTissues, fileTissues{j}) & ...
                                               combinedionMode==1);
                      end
                end
            end
        end
        % remove second half
        combinedTissues(combinedionMode==1)=[];
        combinedDiet(combinedionMode==1)=[];
        combinedType(combinedionMode==1)=[];
        combinedMouse(combinedionMode==1)=[];
        combinedIntensity(:, combinedionMode==1)=[];
    else
        currentIonMode = unique(combinedionMode)*ones(size(combinedIntensity,1),1);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perform quantile normalization of the joint compounds
    combinedIntensitiesNorm = combinedIntensity;
    combinedIntensitiesNorm(:, ismember(combinedTissues, 'Liver') |...
                               ismember(combinedTissues, 'Serum')) = [];
    combinedIntensitiesNorm(combinedIntensitiesNorm==1) = NaN;
    combinedIntensitiesNorm(combinedIntensitiesNorm==5000) = NaN;
    combinedIntensitiesNorm(combinedIntensitiesNorm==0) = NaN;
    combinedIntensitiesNorm = quantilenorm(combinedIntensitiesNorm);
    combinedIntensitiesNormGI = combinedIntensitiesNorm;
    % normalize liver separately
    combinedIntensitiesNorm = combinedIntensity;
    combinedIntensitiesNorm(:, ~ismember(combinedTissues, 'Liver')) = [];
    combinedIntensitiesNorm(combinedIntensitiesNorm==1) = NaN;
    combinedIntensitiesNorm(combinedIntensitiesNorm==5000) = NaN;
    combinedIntensitiesNorm(combinedIntensitiesNorm==0) = NaN;
    combinedIntensitiesNorm = quantilenorm(combinedIntensitiesNorm);
    combinedIntensitiesNormLiver = combinedIntensitiesNorm;
    % normalize serum separately
    combinedIntensitiesNorm = combinedIntensity;
    combinedIntensitiesNorm(:, ~ismember(combinedTissues, 'Serum')) = [];
    combinedIntensitiesNorm(combinedIntensitiesNorm==1) = NaN;
    combinedIntensitiesNorm(combinedIntensitiesNorm==5000) = NaN;
    combinedIntensitiesNorm(combinedIntensitiesNorm==0) = NaN;
    combinedIntensitiesNorm = quantilenorm(combinedIntensitiesNorm);
    combinedIntensitiesNormSerum = combinedIntensitiesNorm;
    % combine all tissues
    combinedIntensitiesNorm = combinedIntensity;
    combinedIntensitiesNorm(:, ~(ismember(combinedTissues, 'Liver') |...
                               ismember(combinedTissues, 'Serum'))) = combinedIntensitiesNormGI;
    combinedIntensitiesNorm(:, ismember(combinedTissues, 'Liver')) = combinedIntensitiesNormLiver;
    combinedIntensitiesNorm(:, ismember(combinedTissues, 'Serum')) = combinedIntensitiesNormSerum;

    combinedIntensitiesNorm(isnan(combinedIntensitiesNorm)) = intensityNoise; %noise level
    
    combinedDataNorm_cell{method_i} = combinedIntensitiesNorm;
    combinedDataTissues_cell{method_i} = combinedTissues;
    combinedDataDiet_cell{method_i} = combinedDiet;
    combinedDataMouse_cell{method_i} = combinedType;
    combinedDataionMode_cell{method_i} = currentIonMode;
    combinedIonMZ_cell{method_i} = changingMets_merged_mass;
    combinedIonRT_cell{method_i} = changingMets_merged_RT;
    combinedIonSpectrum_cell{method_i} = changingMets_merged_spectrum;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% merge three methods
changingMets_merged_RT  = [combinedIonRT_cell{1};...
                           combinedIonRT_cell{2};...
                           combinedIonRT_cell{3};...
                           ];
changingMets_merged_mass = [combinedIonMZ_cell{1};...
                          combinedIonMZ_cell{2};...
                          combinedIonMZ_cell{3};...
                           ];
changingMets_merged_Spectrum  = [combinedIonSpectrum_cell{1};...
                           combinedIonSpectrum_cell{2};...
                           combinedIonSpectrum_cell{3};...
                           ];
combinedIntensitiesNorm = [combinedDataNorm_cell{1};...
                           combinedDataNorm_cell{2};...
                           combinedDataNorm_cell{3};...
                           ]; 
changingMets_merged_method = [repmat(curMethod_unique(1), size(combinedIonRT_cell{1}));...
                              repmat(curMethod_unique(2), size(combinedIonRT_cell{2}));...
                              repmat(curMethod_unique(3), size(combinedIonRT_cell{3}))...
                              ];
changingMets_merged_mode = [combinedDataionMode_cell{1};...
                           combinedDataionMode_cell{2};...
                           combinedDataionMode_cell{3};...
                           ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%annotate samples based on MZ
joinedMzRT = [changingMets_merged_mass changingMets_merged_RT];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[joinedAnn] = annotateIonsByMZ(joinedMzRT, annKEGGAGORA.EXACT_MASS,...
                                            annKEGGAGORA.KEGGID,...
                                            annKEGGAGORA.CompoundNAME,...
                                            massThreshold,[],[],...
                                            annKEGGAGORA.FORMULA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear negMzRT posMzRT negAnnotation posAnnotation posData negData negSamples posSamples  
% create list of KEGG IDs for pathway projector
keggCompoundsID = {joinedAnn(:).annID}';
allKEGGID = [];
for i=1:length(keggCompoundsID)
    allKEGGID = [allKEGGID; keggCompoundsID{i}];
end
allKEGGID = unique(allKEGGID);
% prepare for kegg mapping
allKEGGID = cellfun(@(x) strcat(x, ' red'), allKEGGID, 'unif', 0);
% leave only KEGG ids
allKEGGID(cellfun(@(x) ~contains(x,'cpd:'), allKEGGID)) = [];
% replace cpd:
allKEGGID = cellfun(@(x) strrep(x,'cpd:', ''), allKEGGID, 'unif',0);
% number of annotated ions
nnz(cellfun(@(x) ~isempty(x), keggCompoundsID))
% number of uniquely annotated ions 
nnz(arrayfun(@(x) length(x.annIDX), joinedAnn) == 1)
clear allKEGGID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save annotation to file
fid = fopen([outputFolder 'metabolites_allions_combined_formulas.csv'], 'w');
fprintf(fid, 'MZ\tRT\tMethod\tMode\tSpectrum\tCompoundID\tCompoundName\tCompoundFormula\tMZdelta\tIDX\tmaxIntensity\tmeanintensity\tmedianIntensity\tNumDetectedSamples\n');
onlyAnnotated = 0;%1;
for i=1:length(joinedAnn)
    if onlyAnnotated == 1
        if ~isempty(joinedAnn(i).annID)
            fprintf(fid, '%.3f\t',changingMets_merged_mass(i));
            fprintf(fid, '%.3f\t',changingMets_merged_RT(i));
            fprintf(fid, '%s\t',changingMets_merged_method{i});
            fprintf(fid, '%d\t',changingMets_merged_mode(i));
            fprintf(fid, '%s\t',changingMets_merged_Spectrum{i});
           
            for j=1:length(joinedAnn(i).annID)
                fprintf(fid, '%s;',joinedAnn(i).annID{j});
            end
            fprintf(fid, '\t');
            for j=1:length(joinedAnn(i).annNames)
                fprintf(fid, '%s;',joinedAnn(i).annNames{j});
            end
            fprintf(fid, '\t');
            for j=1:length(joinedAnn(i).annFormula)
                fprintf(fid, '%s;',joinedAnn(i).annFormula{j});
            end
            fprintf(fid, '\t');
            for j=1:length(joinedAnn(i).annMzdelta)
                fprintf(fid, '%3f;',joinedAnn(i).annMzdelta(j));
            end
            fprintf(fid, '\t');
            for j=1:length(joinedAnn(i).annIDX)
                fprintf(fid, '%d;',joinedAnn(i).annIDX(j));
            end
            fprintf(fid, '\t');
            curIntensity = combinedIntensitiesNorm(i,:);
            curIntensity(curIntensity==intensityNoise)=NaN;
            fprintf(fid, '%.3f\t',nanmax(curIntensity));
            fprintf(fid, '%.3f\t',nanmean(curIntensity));
            fprintf(fid, '%.3f\t', nanmedian(curIntensity));
            fprintf(fid, '%d', nnz(~isnan(curIntensity)));
            fprintf(fid, '\n');
        end
    else
         if ~isempty(joinedAnn(i).annID)
            fprintf(fid, '%.3f\t',changingMets_merged_mass(i));
            fprintf(fid, '%.3f\t',changingMets_merged_RT(i));
            fprintf(fid, '%s\t',changingMets_merged_method{i});
            fprintf(fid, '%d\t',changingMets_merged_mode(i));
            fprintf(fid, '%s\t',changingMets_merged_Spectrum{i});
           
            for j=1:length(joinedAnn(i).annID)
                fprintf(fid, '%s;',joinedAnn(i).annID{j});
            end
            fprintf(fid, '\t');
            for j=1:length(joinedAnn(i).annNames)
                fprintf(fid, '%s;',joinedAnn(i).annNames{j});
            end
            fprintf(fid, '\t');
            for j=1:length(joinedAnn(i).annFormula)
                fprintf(fid, '%s;',joinedAnn(i).annFormula{j});
            end
            fprintf(fid, '\t');
            for j=1:length(joinedAnn(i).annMzdelta)
                fprintf(fid, '%3f;',joinedAnn(i).annMzdelta(j));
            end
            fprintf(fid, '\t');
            for j=1:length(joinedAnn(i).annIDX)
                fprintf(fid, '%d;',joinedAnn(i).annIDX(j));
            end
            fprintf(fid, '\t');
            curIntensity = combinedIntensitiesNorm(i,:);
            curIntensity(curIntensity==intensityNoise)=NaN;
            fprintf(fid, '%.3f\t',nanmax(curIntensity));
            fprintf(fid, '%.3f\t',nanmean(curIntensity));
            fprintf(fid, '%.3f\t', nanmedian(curIntensity));
            fprintf(fid, '%d', nnz(~isnan(curIntensity)));
            fprintf(fid, '\n');
         else
            fprintf(fid, '%.3f\t',changingMets_merged_mass(i));
            fprintf(fid, '%.3f\t',changingMets_merged_RT(i));
            fprintf(fid, '%s\t',changingMets_merged_method{i});
            fprintf(fid, '%d\t',changingMets_merged_mode(i));
            fprintf(fid, '%s\t',changingMets_merged_Spectrum{i});
           
            fprintf(fid, '\t\t\t\t\t');
            curIntensity = combinedIntensitiesNorm(i,:);
            curIntensity(curIntensity==intensityNoise)=NaN;
            fprintf(fid, '%.3f\t',nanmax(curIntensity));
            fprintf(fid, '%.3f\t',nanmean(curIntensity));
            fprintf(fid, '%.3f\t', nanmedian(curIntensity));
            fprintf(fid, '%d', nnz(~isnan(curIntensity)));
            fprintf(fid, '\n');
         end
    end
end
fclose(fid);
            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save normalized data to file

% rename distal colon (DC) to feces
combinedTissues(cellfun(@(x) isequal(x, 'DC'), combinedTissues)) = {'Feces'};
%  remove CVR
removeCVR = cellfun(@(x) contains(x,'CVR'), combinedType);
combinedIntensitiesNorm(:, removeCVR)=[];
combinedDiet(removeCVR)=[];
combinedTissues(removeCVR)=[];
combinedType(removeCVR)=[];
combinedMouse(removeCVR) = [];

fid = fopen([outputFolder 'metabolites_allions_combined_norm_intensity.csv'], 'w');
fprintf(fid, 'MZ\tRT\tMethod\tMode\tSpectrum');
for i=1:length(combinedMouse)
    fprintf(fid, '\t%s_%s_%s_%s', combinedDiet{i}, ...
        combinedType{i}, combinedTissues{i}, combinedMouse{i});
end
fprintf(fid, '\n');
for i=1:size(combinedIntensitiesNorm,1)
    fprintf(fid, '%.3f\t',changingMets_merged_mass(i));
    fprintf(fid, '%.3f\t',changingMets_merged_RT(i));
    fprintf(fid, '%s\t',changingMets_merged_method{i});
    fprintf(fid, '%d\t',changingMets_merged_mode(i));
    fprintf(fid, '%s',changingMets_merged_Spectrum{i});
%     fprintf(fid, '%.3f\t',metaboliteData.MZ(i));
%     fprintf(fid, '%.3f\t',metaboliteData.RT(i));
%     fprintf(fid, '%s\t',metaboliteData.Method{i});
%     fprintf(fid, '%d\t',metaboliteData.Mode(i));
%     fprintf(fid, '%s',metaboliteData.Spectrum{i});
    
    for j=1:size(combinedIntensitiesNorm,2)
        fprintf(fid, '\t%f',combinedIntensitiesNorm(i,j));
    end
    fprintf(fid, '\n');
end
fclose(fid);
