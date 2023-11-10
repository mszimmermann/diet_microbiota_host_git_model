%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Pipeline for the analysis of community profiles (untargeted metabolites)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('C:\Users\mazimmer\Documents\Documents\MATLAB\KEGGquery\')
addpath('C:\Users\mazimmer\Documents\Documents\MATLAB\diffAnalysis\\')
addpath(genpath('C:\Users\mazimmer\Documents\Documents\MATLAB\GutModels\metabolomics'))

% read annotation file
annKEGGAGORA = readtable('C:\Users\mazimmer\Documents\Documents\MATLAB\GutModels\compound_masses_merged_KEGG_AGORA_HMDB_082019.csv');

fileFolder = 'C:\Users\mazimmer\Documents\Documents\MATLAB\GutModels\metabolomics\CommunityData\';
folderFiles = dir(fileFolder);
fileNames = cell(size(folderFiles));
for i=1:length(fileNames)
    fileNames{i} = folderFiles(i).name;
end
% leave only MZcomm files
fileNames(cellfun(@(x) ~contains(x, 'MZcomm'), fileNames))=[];
% remove mouse batch
fileNames(cellfun(@(x) contains(x, 'M001_'), fileNames))=[];

massThreshold = 0.001;
RTthreshold = 0.15;
intensityNoise = 5000;

% prepare cells for batches
dataMatrix_batch = cell(size(fileNames));
dataTime_batch = cell(size(fileNames));
dataSamples_batch = cell(size(fileNames));
compound_batch = cell(size(fileNames));
compoundSpectrum_batch = cell(size(fileNames));
        
for fi = 1:length(fileNames)
    curTable = readtable([fileFolder fileNames{fi}]);
    
%         % correct by internal standard
%         % do not use all the STD!!!!
%         IS = {'IS_YOH', 'IS_CAF', 'IS_SUL'};%, 'IS_IPR'}; % internal standard with which to correctintensities
%         for i = 1:length(TimePoints)
%             IdxIS = ismember((cellfun(@(x) x(1:6), RawData.(TimePoints{i}).Compounds, 'UniformOutput', false)), IS);
%             TempMatrix = RawData.(TimePoints{i}).IntensitiesRaw;
%             TempIS = TempMatrix(:,IdxIS); %intensities of internal std of this time point
%             MeanIS = mean(TempIS);
%             CorrectionMatrix = ones(size(TempIS));
%             for p = 1:length(MeanIS)
%                 CorrectionMatrix(:,p)= TempIS(:,p)./MeanIS(p); %calcuate fold change from the mean
%             end
%             Correction = mean(CorrectionMatrix,2);
%             for p = 1: length(Correction)
%                 TempMatrix(p,:)= TempMatrix(p,:)./Correction(p);
%                 %TempMatrix(p,:)= TempMatrix(p,:)./1;
%             end
%             % eliminate samples with bad injections
%             TempInjectionSum = sum(TempMatrix,2);
%             TempMean = mean(TempInjectionSum);
%             TempSTD = std(TempInjectionSum);
%             Outlayers = find(TempInjectionSum < (TempMean/2));
%             TempMatrix(Outlayers, :)=NaN;
% 
%             %RawData.(PlateCount).(PoolCount).IntensitiesCorrected = struct;
%             RawData.(TimePoints{i}).IntensitiesCorrected = TempMatrix;
%             %reverse normalization
%             %RawData.(TimePoints{i}).IntensitiesCorrected = RawData.(TimePoints{i}).IntensitiesRaw;
%         end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create matrix of community and time versus compounds
    % select compounds
    
    dataColumns = contains(curTable.Properties.VariableNames, 'MZ012');
    dataMatrix = table2array(curTable(:, dataColumns));
    dataSamples = curTable.Properties.VariableNames(dataColumns)';
    dataCondition = cell(size(dataSamples));
    dataTime = zeros(size(dataSamples));
    for j=1:length(dataSamples)
        curSample = strsplit(dataSamples{j}, '_');
        dataCondition{j} = curSample{3};
        dataTime(j) = str2double(curSample{4}(2:3));
    end

    dataMatrix_batch{fi} = dataMatrix;
    dataTime_batch{fi} = dataTime;
    dataSamples_batch{fi} = dataCondition;
    compound_batch{fi} = curTable.Compound;
    compoundSpectrum_batch{fi} = curTable.MS1CompositeSpectrum;
end
   
clear dataMatrix dataTime dataCondition dataSamples
clear dataColumns folderFiles
clear curSample curTable

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% merge metabolites across batches
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a matrix with metabolites across tissues
% first merge metabolites across tissues
changingMets = cell(100000,1);
changingMetsSpectrum = cell(100000,1);

idx = 1;
for i=1:length(compound_batch)
    curMets = compound_batch{i};
    curMetsSpectrum = compoundSpectrum_batch{i};

    changingMets(idx:idx+length(curMets)-1) = curMets;
    changingMetsSpectrum(idx:idx+length(curMets)-1) = curMetsSpectrum;
    idx = idx+length(curMets);
end
changingMets(idx:end) = [];
changingMetsSpectrum(idx:end) = [];
clear curMets curMetsSpectrum

% remove identical compouds and sort by name
[changingMets,idx] = unique(changingMets);
changingMetsSpectrum = changingMetsSpectrum(idx);
%changingMetsIonMode = changingMetsIonMode(idx);
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
                                                                   ones(size(changingMets)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get idx conversion from each ion file to the merged
curCompoundsMergedIDXconversion = zeros(length(changingMets_merged),...
                                        length(compound_batch));
for d = 1:length(compound_batch)
    curMZ = cellfun(@(x) str2double(x(1:strfind(x, '@')-1)), compound_batch{d});
    curRT = cellfun(@(x) str2double(x(strfind(x, '@')+1:end)), compound_batch{d});

    tic
    for i=1:size(changingMets_merged_mass,1)

        curMetMass = changingMets_merged_mass(i);
        curMetRT = changingMets_merged_RT(i);

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
    toc
    fprintf('Merged metabolites for batch number %d\n', d);
end        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate intensity matrix for all metabolites
batchSize = 300;
combinedIntensity = zeros(size(curCompoundsMergedIDXconversion,1),...
                            batchSize * length(compound_batch));
combinedSamples = cell(batchSize*length(compound_batch),1);
combinedTime = zeros(batchSize*length(compound_batch),1);
idx = 1;
for d = 1:length(compound_batch)
    curIntensitiesRaw =  dataMatrix_batch{d};

    sampleCondition = dataSamples_batch{d};
    sampleTime = dataTime_batch{d};
    
    combinedIntensity(curCompoundsMergedIDXconversion(:,d)~=0,...
                      idx:idx+size(sampleTime)-1) = ...
         curIntensitiesRaw(curCompoundsMergedIDXconversion(curCompoundsMergedIDXconversion(:,d)~=0,d),:);

    combinedSamples(idx:idx+size(sampleTime)-1) = sampleCondition;
    combinedTime(idx:idx+size(sampleTime)-1) = sampleTime;
    idx=idx+length(sampleTime);
    fprintf('Merged selected metabolites for batch %d \n', d);
end
combinedIntensity(:, idx:end) = [];
combinedSamples(idx:end) = [];
combinedTime(idx:end) = [];
clear datafile curIntensitiesRaw sampleCondition sampleTime

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform quantile normalization of the joint compounds
% combinedIntensitiesNorm(combinedIntensitiesNorm==1) = NaN;
% combinedIntensitiesNorm(combinedIntensitiesNorm==intensityNoise) = NaN;
% combinedIntensitiesNorm(combinedIntensitiesNorm==0) = NaN;
% combinedIntensitiesNorm = quantilenorm(combinedIntensitiesNorm);
% combinedIntensitiesNorm(isnan(combinedIntensitiesNorm)) = intensityNoise;
combinedIntensity(combinedIntensity==1) = NaN;
combinedIntensity(combinedIntensity==intensityNoise) = NaN;
combinedIntensity(combinedIntensity==0) = NaN;

combinedIntensity(isnan(combinedIntensity)) = intensityNoise;
combinedIntensity(combinedIntensity==intensityNoise) = NaN;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% try to normalize by IS 
IS = {194.080 'IS_CAF';...
      253.052  'IS_SUL';...390.171 'IS_YOH';...
      280.109 'IS_IPR'};%, 'IS_IPR'}; % internal standard with which to correctintensities
[~,ISidx] = cellfun(@(x) min(abs(changingMets_merged_mass-x)),...
                     IS(:,1));

TempIS = combinedIntensity(ISidx,:); %intensities of internal std of this time point
MeanIS = mean(TempIS,2);
CorrectionMatrix = ones(size(TempIS))';
for p = 1:length(MeanIS)
    CorrectionMatrix(:,p)= TempIS(p,:)./MeanIS(p); %calcuate fold change from the mean
end
Correction = mean(CorrectionMatrix,2);
combinedIntensity_norm = combinedIntensity;
for p = 1: length(Correction)
    combinedIntensity_norm(:,p) = combinedIntensity_norm(:,p) / Correction(p);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%annotate samples based on MZ
negMzRT = [changingMets_merged_mass changingMets_merged_RT];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[negAnnotation] = annotateIonsByMZ(negMzRT, annKEGGAGORA.EXACT_MASS,...
                                            annKEGGAGORA.KEGGID,...
                                            annKEGGAGORA.CompoundNAME,...
                                            massThreshold,[],[],...
                                            annKEGGAGORA.FORMULA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
joinedAnn = negAnnotation;
joinedMzRT = negMzRT;
%joinedSamples = SamplesNEG;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mycolormap = jet;
%permjet = fliplr(randperm(size(mycolormap,1)));
mycolormap = mycolormap(permjet,:);

% extract unique time and samples
combinedCondition = cellfun(@(x) x(1:end-1), combinedSamples, 'unif', 0);
combinedCondition_unique = unique(combinedCondition);
combinedTime_unique = unique(combinedTime);
% remove bad controls and mouse communities
% dataSamplesType_unique(ismember(dataSamplesType_unique,...
%                        {'CTR04', 'CTR06' 'CTR02' 'CTR08' ...
%                         'M001' 'M002' 'M003' 'M004'})) = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for each compound, plot time profile of each community over time
myfontSize = 8;
mycolors = [repmat([0 0 0], nnz(cellfun(@(x) contains(x, 'CTR'), combinedCondition_unique)),1);...
            mycolormap];
fileName = ['profiles_normalizedByIS_annotated_metabolite_profiles_at_least_50_per_community.ps'];

    
%compoundNames = {'102.03@1.5627'};
compoundMass = '214.132@1.9';% dethiobiotin'226.0956';%porphobilinogen' '244.088';%biotin '465.3';%glycocholate '244.08';%biotin '219.110';%pantothenate '187.0636';%indoleacryllate '111.079';%histamine '99.068';%hydroxymethylbutanetitrile '414.16';%diltiazem '376.1';%'146.0';%''150.05';%340.203';%'321.14';%'372.15';%'414.16';%'342.2';%'340.20';%'155.0';%'121.0';%'180.06';%'111.04';%'112.02';%'244.06';%'149.0';% '244.08';%'214.13';%'120.04';%'252.0';%'147.05';%glu %'867.';%succ-coA '75.03';%glycine '131.05';%aminolevulinate '226.0';%porphobilinogen
compoundNames = changingMets_merged(contains(changingMets_merged, compoundMass));

%compoundNames = changingMets_merged;%

mycolors = mycolormap;
lwidth = 1;

fig = figure('units','normalized','outerposition',[0 0 1 1]);
spX = 5;
spY = 8;
nonnanThreshold = 50;
for i = 1:length(compoundNames)
    cmpd_i = find(cellfun(@(x) isequal(x, compoundNames{i}),changingMets_merged));
   
    % IF ANNOTATED, PLOT
    if (~isempty(joinedAnn(cmpd_i).annID))% && ...
       %(nnz(~isnan(combinedIntensity_norm(cmpd_i,:)))>nonnanThreshold)
        figure;%clf
        hold on
        %figure
    spidx = 1;
    meanPerTime_comm = zeros(length(combinedTime_unique),length(combinedCondition_unique));
    stdPerTime_comm = zeros(length(combinedTime_unique),length(combinedCondition_unique));
    for sample_i=[2 9:length(combinedCondition_unique)]
        h=[];
        hidx = 1;
       % subplot(spX,spY,spidx)
        hold on 
    
        % save FC at last time point for correlation
        curData = combinedIntensity_norm(cmpd_i,...
                                    ismember(combinedCondition,...
                                    combinedCondition_unique(sample_i)));
        curTime = combinedTime(ismember(combinedCondition,...
                                             combinedCondition_unique(sample_i)));
        %calculate mean value per time point and plot as horizontal lines
        curTime_unique = unique(curTime);
        meanPerTime = zeros(length(curTime_unique),1);
        stdPerTime = zeros(length(curTime_unique),1);
        meanFCPerTime = zeros(length(curTime_unique),1);
        stdFCPerTime = zeros(length(curTime_unique),1);
        for j=1:length(curTime_unique)
            meanPerTime(j) = nanmean(curData(curTime==curTime_unique(j)));
            stdPerTime(j) = nanstd(curData(curTime==curTime_unique(j)));
            meanFCPerTime(j) = nanmean(log2(curData(curTime==curTime_unique(j)) ./...
                                    curData(curTime==0)));
            stdFCPerTime(j) = nanstd(log2(curData(curTime==curTime_unique(j)) ./...
                                    curData(curTime==0)));
        end
        % plot the single datapoints and means and fit
            %h(sample_i) = plot(curTime_unique, meanPerTime, 'Color', mycolors(sample_i,:));
            %errorbar(curTime_unique, meanPerTime, stdPerTime, '.', 'Color', mycolors(sample_i,:))
            
%             h(sample_i) = plot(curTime_unique, meanFCPerTime,...
%                                 'Color', mycolors(sample_i,:),...
%                                 'LineWidth', lwidth);
%             errorbar(curTime_unique, meanFCPerTime, stdFCPerTime, '.',...
%                     'Color', mycolors(sample_i,:),...
%                     'LineWidth', lwidth)
            %meanPerTime(isnan(meanPerTime))=0;
            
%             h(hidx) = plot(curTime_unique, meanPerTime,...
%                            'Color', mycolors(sample_i,:),...
%                            'LineWidth', lwidth);
%              hold on
%             hidx = hidx+1;
%             
%             errorbar(curTime_unique, meanPerTime, stdPerTime, '.',...
%                     'Color', mycolors(sample_i,:),...
%                     'LineWidth', lwidth)
%           xlim([0 max(curTime_unique)])

            h(hidx) = bar(sample_i, meanPerTime(2),...
                           'FaceColor', mycolors(sample_i,:));
            hold on
            hidx = hidx+1;
            
            errorbar(sample_i, meanPerTime(2), stdPerTime(2), '.',...
                    'Color', 'k',...mycolors(sample_i,:),...
                    'LineWidth', lwidth)
            xlim([0 length(combinedCondition_unique)])
            axis square
            spidx = spidx+1;
            meanPerTime_comm(:,sample_i) = meanPerTime;
            stdPerTime_comm(:,sample_i) = stdPerTime;
    end
    %suptitle(changingMets_merged(cmpd_i));
    title([changingMets_merged(cmpd_i),...
           strjoin(joinedAnn(cmpd_i).annID),...
           strjoin(joinedAnn(cmpd_i).annNames)]);
    orient landscape
%     print(gcf, '-painters', '-dpsc', '-r600', '-append', '-bestfit', ...
%         fileName)
    end
    %spidx = spidx+1;
%     if spidx > 4
%         orient landscape
%         legend(h, [drugName; compoundNames(mycompounds)]);
% %         print(gcf, '-painters', '-dpsc', '-r600', '-append', '-bestfit', ...
% %             fileName)
% %         clf
%         figure
%         spidx = 1;
%     end
end
% orient landscape
% print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
%               'profile_porphobilinogen_communities.pdf')
%             'profile_dethiobiotin_communities.pdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sort communities by time 1h and plot
sample_i=[2 9:length(combinedCondition_unique)];
plotdata = meanPerTime_comm(2,sample_i);
plotstd = stdPerTime_comm(2,sample_i);
[plotdata, sortidx] = sort(plotdata);
plotstd = plotstd(sortidx);

figure
hidx = 1;
for sample_i=1:length(plotdata)
    h(hidx) = bar(sample_i, plotdata(sample_i),...
                           'FaceColor', mycolors(sample_i+8,:));
    hold on
    hidx = hidx+1;

    errorbar(sample_i, plotdata(sample_i), plotstd(sample_i), '.',...
            'Color', 'k',...mycolors(sample_i,:),...
            'LineWidth', lwidth)
end
xlim([0 nnz(~isnan(plotdata))+1])
axis square
xlabel('Communities')
ylabel('Porphobilinogen after 1h of incubation')
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
    'barplot_communities_porphobilinogen_1h_incubation.pdf')
            

% correlate with all OTUs
metOTUcorr = zeros(size(MVOTU,1),1);
metOTUcorrP = zeros(size(MVOTU,1),1);
metTimepointIDX = 2;
for i=1:size(MVOTUdata,1)
    [metOTUcorr(i),...
     metOTUcorrP(i)] = corr(MVOTUdata(i,5:end)',...
                        meanPerTime_comm(metTimepointIDX,9:end)',...
                        'rows','complete');
end
    




% fig = figure('units','normalized','outerposition',[0 0 1 1]);
% spX = 2;
% spY = 2;
% spidx = 1;
% 
% for sample_i=1:length(dataSamplesType_unique)
%     h=[];
%     hidx = 1;
%     subplot(spX,spY,spidx)
%     hold on 
%     for fi=1:length(fileNames)
%         if fi==1
%             dataMatrix = dataMatrixDRUG;
%             dataTime = dataTimeDRUG;
%             dataSamples = dataSamplesDRUG;
%             compoundNames = compoundNamesDRUG;
%             mycompounds = find(contains(upper(compoundNames), drugName));
%             yyaxis right
%         elseif fi==2
%             dataMatrix = dataMatrixMET;
%             dataTime = dataTimeMET;
%             dataSamples = dataSamplesMET;
%             compoundNames = compoundNamesMET;
%             mymasses = zeros(size(compoundNames));
%             for j=1:length(drugMetMasses)
%                 mymasses = mymasses + contains(compoundNames, drugMetMasses{j});
%             end
%             mycompounds = find(contains(upper(compoundNames), drugName));% & mymasses>0);
%             yyaxis left
%         end
%         
%         for i = 1:length(mycompounds)
%             cmpd_i = mycompounds(i);
%             
%             % save FC at last time point for correlation
%             curData = dataMatrix(ismember(dataSamplesType, dataSamplesType_unique(sample_i)),...
%                                  cmpd_i);
%             curTime = dataTime(ismember(dataSamplesType, dataSamplesType_unique(sample_i)));
% 
%             %calculate mean value per time point and plot as horizontal lines
%             curTime_unique = unique(curTime);
%             meanPerTime = zeros(length(curTime_unique),1);
%             stdPerTime = zeros(length(curTime_unique),1);
%             meanFCPerTime = zeros(length(curTime_unique),1);
%             stdFCPerTime = zeros(length(curTime_unique),1);
%             for j=1:length(curTime_unique)
%                 meanPerTime(j) = mean(curData(curTime==curTime_unique(j)));
%                 stdPerTime(j) = std(curData(curTime==curTime_unique(j)));
%                 meanFCPerTime(j) = mean(log2(curData(curTime==curTime_unique(j)) ./...
%                                         curData(curTime==0)));
%                 stdFCPerTime(j) = std(log2(curData(curTime==curTime_unique(j)) ./...
%                                         curData(curTime==0)));
%             end
%             % plot the single datapoints and means and fit
%             %h(sample_i) = plot(curTime_unique, meanPerTime, 'Color', mycolors(sample_i,:));
%             %errorbar(curTime_unique, meanPerTime, stdPerTime, '.', 'Color', mycolors(sample_i,:))
%             
% %             h(sample_i) = plot(curTime_unique, meanFCPerTime,...
% %                                 'Color', mycolors(sample_i,:),...
% %                                 'LineWidth', lwidth);
% %             errorbar(curTime_unique, meanFCPerTime, stdFCPerTime, '.',...
% %                     'Color', mycolors(sample_i,:),...
% %                     'LineWidth', lwidth)
%             h(hidx) = plot(curTime_unique, meanPerTime,...
%                                 mylines{min(hidx, length(mylines))}, 'Color', mycolors(sample_i,:),...
%                                 'LineWidth', lwidth);
%             hidx = hidx+1;
%             errorbar(curTime_unique, meanPerTime, stdPerTime, '.',...
%                     'Color', mycolors(sample_i,:),...
%                     'LineWidth', lwidth)
%             xlim([0 max(curTime_unique)])
%             axis square
%         end
%     end
%     spidx = spidx+1;
%     if spidx > 4
%         orient landscape
%         legend(h, [drugName; compoundNames(mycompounds)]);
% %         print(gcf, '-painters', '-dpsc', '-r600', '-append', '-bestfit', ...
% %             fileName)
% %         clf
%         figure
%         spidx = 1;
%     end
% end
% legend(h, [drugName; compoundNames(mycompounds)]);
% orient landscape
% print(gcf, '-painters', '-dpsc', '-r600', '-append', '-bestfit', ...
%             fileName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% correlate with shortBRED
shortBREDtable = readtable('C:\Users\mazimmer\Documents\Documents\MATLAB\GutModels\metabolomics\CommunityData\shortbred_dietenzyme_results__merged.csv');
shortBREDcounts = table2array(shortBREDtable(:,contains(shortBREDtable.Properties.VariableNames, '_Count')));
shortBREDprots = shortBREDtable{:,1};
% remove zero proteins
shortBREDprots(sum(shortBREDcounts,2)==0)=[];
shortBREDcounts(sum(shortBREDcounts,2)==0,:)=[];
% correlate shortBRED with metabolites
fig = figure('units','normalized','outerposition',[0 0 1 1]);
spX = 3;
spY = 4;
spidx = 1;
metTimepointIDX=2;
for i=1:size(shortBREDprots,1)
    subplot(spX, spY, spidx)
    hold on
    for j=1:length(shortBREDcounts(i,:))
        plot(shortBREDcounts(i,j), meanPerTime_comm(metTimepointIDX,8+j),...
             'o', 'MarkerEdgeColor', mycolors(8+j,:),...
             'MarkerFaceColor', mycolors(8+j,:));
    end
    axis square
    [pcc, pccP] = corr(shortBREDcounts(i,:)', meanPerTime_comm(metTimepointIDX,9:end)', 'rows','complete');
    title({shortBREDprots{i}(1:min(25, length(shortBREDprots{i}))),...
           sprintf('%.2f %.2f', pcc, pccP)});
    spidx = spidx+1;
end
% orient landscape
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
            'shortbred_corr_porphobilinogen_communities.pdf')       
%           'shortred_corr_dethiobiotin_communities.pdf')
figure
hold on
i=[2,4,12];
plotenzyme = sum(shortBREDcounts(i,:),1);
plot(log2(plotenzyme), meanPerTime_comm(metTimepointIDX,9:end), '.');
lsline
for j=1:length(plotenzyme)
    plot(log2(plotenzyme(j)), meanPerTime_comm(metTimepointIDX,8+j),...
         'o', 'MarkerEdgeColor', mycolors(8+j,:),...
         'MarkerFaceColor', mycolors(8+j,:));
end
axis square
[pcc, pccP] = corr(plotenzyme', meanPerTime_comm(metTimepointIDX,9:end)',...
                    'rows','complete', 'type', 'spearman');
logexp = log2(plotenzyme)';
logexp(isinf(logexp)) = nan;
[pcc, pccP] = corr(logexp, (meanPerTime_comm(metTimepointIDX,9:end))',...
                    'rows','complete');
title({shortBREDprots{i}(1:min(25, length(shortBREDprots{i}))),...
       sprintf('%.2f %.2f', pcc, pccP)});
       
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
     'shortbred_corr_porphobilinogen_communities_singlegene.pdf')       
    % 'shortbred_corr_dethiobiotin_communities_singlegene.pdf')
    
figure
hold on
hidx = 1;
i=12;
logexp = log2(shortBREDcounts(i,:))';
logexp(isinf(logexp)) = 0;

plot(logexp, meanPerTime_comm(metTimepointIDX,9:end), '.');
lsline

sample_i=[2 9:length(combinedCondition_unique)];
for j=1:length(sample_i)
    if sortidx(j)~=1
        plot(logexp(sortidx(j)-1), plotdata(j),...
             'o', 'MarkerEdgeColor', mycolors(8+j,:),...
             'MarkerFaceColor', mycolors(8+j,:));
    end
end
axis square
[pcc, pccP] = corr(shortBREDcounts(i,:)', meanPerTime_comm(metTimepointIDX,9:end)',...
                    'rows','complete', 'type', 'spearman');
logexp = log2(shortBREDcounts(i,:))';
logexp(isinf(logexp)) = 0;
[pcc, pccP] = corr(logexp, (meanPerTime_comm(metTimepointIDX,9:end))',...
                    'rows','complete');
title({shortBREDprots{i}(1:min(25, length(shortBREDprots{i}))),...
       sprintf('%.2f %.2f', pcc, pccP)});
xlim([0 10])
ylim([0 700000])
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
     'scatter_porphobilinogen_synthase_shortbred_communities_log.pdf')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot corr histogram and feature
figure
histogram(metOTUcorr, 100)
hold on 
plot([pcc pcc], [0 20])
axis square
xlabel('Correlation between feature and metabolite')
ylabel('Number of features')
xlim([-0.8 0.8])
% percentile of pcc
%pccPerc = nnz(metOTUcorr<pcc)/length(metOTUcorr)*100;
pccPerc = nnz(metOTUcorr>pcc)/length(metOTUcorr)*100;
text(pcc, 20, sprintf('%d%%', round(pccPerc)))
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
         'histogram_corr_porphobilinogen_and_gene.pdf')
     %'histogram_corr_dethiobiotin_and_gene.pdf')
       
% plot cor rbetween identified gene and features
corridx = find(metOTUcorr<pcc);
fig = figure('units','normalized','outerposition',[0 0 1 1]);
spX = 3;
spY = 4;
spidx = 1;
idx = 12;
for i=1:length(corridx)
    subplot(spX, spY, spidx)
    hold on
    for j=1:length(shortBREDcounts(idx,:))
        plot(shortBREDcounts(idx,j),...
             MVOTUdata(corridx(i),4+j),...
             'o', 'MarkerEdgeColor', mycolors(8+j,:),...
             'MarkerFaceColor', mycolors(8+j,:));
    end
    axis square
    [pcc, pccP] = corr(shortBREDcounts(idx,:)',...
                       MVOTUdata(corridx(i),5:end)',...
                       'rows','complete');
    title({shortBREDprots{idx}(1:min(25, length(shortBREDprots{idx}))),...
           sprintf('%.2f %.2f', pcc, pccP)});
    ylabel(corridx(i))
    spidx = spidx+1;
end
suptitle('Correlation OTU vs porphobilinogen synthase')
orient landscape
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit', ...
        'scatter_corr_OTU_porphobilinogen_gene.pdf')

figure
scatter(shortBREDcounts(1,:)+...
        shortBREDcounts(4,:)+...
        shortBREDcounts(12,:), meanPerTime_comm(metTimepointIDX,9:end));
[pcc, pccP] = corr((shortBREDcounts(1,:)+...
        shortBREDcounts(4,:)+...
        shortBREDcounts(12,:))',...
        meanPerTime_comm(metTimepointIDX,9:end)',...
        'rows','complete');
       
figure
scatter(shortBREDcounts(10,:)+...
        ...shortBREDcounts(3,:)+...
        shortBREDcounts(9,:)+...
        shortBREDcounts(8,:)+...
        shortBREDcounts(5,:)+...
        shortBREDcounts(7,:),-meanPerTime_comm(metTimepointIDX,9:end));

figure
scatter(log2(shortBREDcounts(10,:)+...
        ...shortBREDcounts(3,:)+...
        shortBREDcounts(9,:)+...
        shortBREDcounts(8,:)+...
        shortBREDcounts(5,:)+...
        shortBREDcounts(7,:)),-meanPerTime_comm(metTimepointIDX,9:end));

[pcc, pccP] = corr((shortBREDcounts(10,:)+...
        ...shortBREDcounts(3,:)+...
        shortBREDcounts(9,:)+...
        shortBREDcounts(8,:)+...
        shortBREDcounts(5,:)+...
        shortBREDcounts(7,:))',...
        meanPerTime_comm(metTimepointIDX,9:end)',...
        'rows','complete');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make a table with metabolite abundances for selected metabolites
%compoundNames = {'102.03@1.5627'};
compoundMass = '214.132@1.9';% dethiobiotin'226.0956';%porphobilinogen' '244.088';%biotin '465.3';%glycocholate '244.08';%biotin '219.110';%pantothenate '187.0636';%indoleacryllate '111.079';%histamine '99.068';%hydroxymethylbutanetitrile '414.16';%diltiazem '376.1';%'146.0';%''150.05';%340.203';%'321.14';%'372.15';%'414.16';%'342.2';%'340.20';%'155.0';%'121.0';%'180.06';%'111.04';%'112.02';%'244.06';%'149.0';% '244.08';%'214.13';%'120.04';%'252.0';%'147.05';%glu %'867.';%succ-coA '75.03';%glycine '131.05';%aminolevulinate '226.0';%porphobilinogen
compoundNames = changingMets_merged(contains(changingMets_merged, compoundMass));

select_metaboliteTable = readtable([outputFolder,...
    'table_metabolites_substrate_products_one_enzyme_selected.csv']);

candidate_compound_idx = zeros(size(select_metaboliteTable,1)*3,1);
selected_compound_idx = zeros(size(select_metaboliteTable,1)*3,1);
mz_threshold = 0.002;
idx=1;
for i=1:size(select_metaboliteTable,1)
    cur_mass = select_metaboliteTable.MZ(i);
    find_compound = find(abs(changingMets_merged_mass-cur_mass)<mz_threshold);
    nfound = length(find_compound); 
    if nfound
        candidate_compound_idx(idx:idx+nfound-1) = find_compound;
        selected_compound_idx(idx:idx+nfound-1) = repmat(i, nfound,1);
        idx = idx+nfound;
    end
end
candidate_compound_idx(idx:end) = [];
selected_compound_idx(idx:end) = [];

test = select_metaboliteTable(selected_compound_idx,:);
test2 = joinedMzRT(candidate_compound_idx,:);

combinedIntensity_norm28 = combinedIntensity_norm;
combinedSamples28 = combinedSamples;
combinedTime28 = combinedTime;
% remove 25 and 30
toremove = (contains(combinedSamples28, 'MV25') |...
            contains(combinedSamples28, 'MV30'));
combinedIntensity_norm28(:, toremove)=[];
combinedSamples28(toremove)=[];
combinedTime28(toremove)=[];
% rename communities from 26 to 29
combinedSamples28 = cellfun(@(x) strrep(x, 'MV26', 'MV25'), combinedSamples28, 'unif', 0);
combinedSamples28 = cellfun(@(x) strrep(x, 'MV27', 'MV26'), combinedSamples28, 'unif', 0);
combinedSamples28 = cellfun(@(x) strrep(x, 'MV28', 'MV27'), combinedSamples28, 'unif', 0);
combinedSamples28 = cellfun(@(x) strrep(x, 'MV29', 'MV28'), combinedSamples28, 'unif', 0);


combinedIntensity_norm_table = array2table(combinedIntensity_norm28);
combinedIntensity_norm_table.Properties.VariableNames = strcat(combinedSamples28, '_', ...
    arrayfun(@(x) num2str(x), combinedTime28, 'unif', 0))';
        
candidate_compound_intensity = combinedIntensity_norm_table(candidate_compound_idx,:);
candidate_compound_MZRT = array2table(joinedMzRT(candidate_compound_idx,:));
candidate_compound_MZRT.Properties.VariableNames = {'Matched_MZ', 'Matched_RT'};
candidate_compound_ann = joinedAnn(:,candidate_compound_idx);

selected_compoundAnn = select_metaboliteTable(selected_compound_idx,[1 2 6 7]);

candidateTable = [selected_compoundAnn candidate_compound_MZRT candidate_compound_intensity];

writetable(candidateTable, ...
    [outputFolder, 'table_candidate_sub_prod_pair_metabolites_in_communities.csv'],...
    'WriteRowNames',true);

