%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [meanMatrix, meanMatrix_mets, meanData_columns,meanMatrix_STD] = ...
    prepare_drug_data_for_modeling(outputFolder, figureFolder,...
                                   printflag,use_volume_flag)
% printflag indicates whether to print data to figure

% prepare data from BRV paper for modeling
% table can be downloaded from https://www.science.org/doi/full/10.1126/science.aat9931
dataFolder = '.\InputData\public_data\';
dataFilenames = 'aat9931_tables_s1-26.xlsx';
dataSheetnames = {'Table S4';'Table S17'};

meanData_cell = cell(size(dataSheetnames));
meanDataSTD_cell = cell(size(dataSheetnames));
meanMets_cell = cell(size(dataSheetnames));
Mean_tissues = {'Mean_SI','Mean_SII','Mean_SIII',...
                 'Mean_Cecum', 'Mean_Colon', 'Mean_Feces'};
STD_tissues = {'STD_SI','STD_SII','STD_SIII',...
                 'STD_Cecum', 'STD_Colon', 'STD_Feces'};
% CV stands for any kind of colonized mice in this case
meanData_columns = [strcat('CV_Chow_', Mean_tissues)...
                    strcat('GF_Chow_', Mean_tissues)];

for ds_i = 1:length(dataSheetnames)
    dataSheetname = dataSheetnames{ds_i};
              
    % read drug data
    opts = detectImportOptions([dataFolder dataFilenames],'Sheet', dataSheetname, 'NumHeaderLines',0);
    % read all columns as characters and reformat later on
    opts = setvartype(opts, 'char');
    t = readtable([dataFolder dataFilenames], opts);%'Sheet', dataSheetname, 'NumHeaderLines',0, opts);
    % get column names which are in the third row
    t_colnames = t{3,:};
    % remove spaces in column names
    t_colnames = cellfun(@(x) strtrim(x), t_colnames, 'unif', 0);
    t_colnames = cellfun(@(x) strrep(x,' ', '_'), t_colnames, 'unif', 0);

    if isequal(dataSheetname, 'Table S4')
        % remove first five lines
        t = t(5:end,:);
    end
    if isequal(dataSheetname, 'Table S17')
        t = t(8:end,:);
        % rename Abbreviation and SampleTime
        t_colnames{ismember(t_colnames, 'Abbreviation')} = 'Measured_metabolite';
        t_colnames{ismember(t_colnames, 'SampleTime')} = 'Sample_Time';
        t_colnames = cellfun(@(x) strrep(x, '_treated_', '_'), t_colnames, 'unif', 0);
        t_colnames = cellfun(@(x) strrep(x, '_SI1', '_SI'), t_colnames, 'unif', 0);
        t_colnames = cellfun(@(x) strrep(x, '_SI2', '_SII'), t_colnames, 'unif', 0);
        t_colnames = cellfun(@(x) strrep(x, '_SI3', '_SIII'), t_colnames, 'unif', 0);
    end

    % set correct variable names that were in row 3 in the original table
    t.Properties.VariableNames = t_colnames;
    % convert time from text to numbers
    t.Sample_Time = cellfun(@(x) str2double(x), t.Sample_Time);

    if  isequal(dataSheetname, 'Table S4')
        Measured_mets = strcat(t.Administered_Drug, '_',t.Measured_metabolite);
    else
        Measured_mets = t.Measured_metabolite;
    end
    Measured_mets_unique = unique(Measured_mets);
    % remove Thymine as it is measured only in liver
    Measured_mets_unique(cellfun(@(x) contains(x, 'Thymine'), Measured_mets_unique))=[];
    % remove 7-acetamido-CLZ
    Measured_mets_unique(cellfun(@(x) contains(x, '7-acetamido-CLZ'), Measured_mets_unique))=[];
    % remove phenolic and glucuronyl compounds
    Measured_mets_unique(cellfun(@(x) contains(x, 'Glucuronyl'), Measured_mets_unique))=[];
    Measured_mets_unique(cellfun(@(x) contains(x, 'Phenolic'), Measured_mets_unique))=[];
    % get unique time points
    t_time_unique = unique(t.Sample_Time);
    % remove time points between 0 and 3
    t_time_unique(t_time_unique>0 & t_time_unique<3)=[];
    t_sample_type_unique = reshape(unique(t.Sample_Type), 2,[]);

    % for BRV data, add d4554 vs GF condition
    if ismember('d4554', t_sample_type_unique) && ismember('GF', t_sample_type_unique)
        addcol = size(t_sample_type_unique,2);
        t_sample_type_unique{1,addcol+1} = 'd4554';
        t_sample_type_unique{2,addcol+1} = 'GF';
    end

    meanTable_data = zeros(length(Measured_mets_unique) * length(t_time_unique)*...
                           size(t_sample_type_unique,2),...
                           length(Mean_tissues)*2);
    meanTable_dataSTD = zeros(length(Measured_mets_unique) * length(t_time_unique)*...
                           size(t_sample_type_unique,2),...
                           length(Mean_tissues)*2);
    meanTable_mets = cell(length(Measured_mets_unique) * length(t_time_unique)*...
                           size(t_sample_type_unique,2));                   
    idx=1;
    for met_i=1:length(Measured_mets_unique)
        for time_i=1:length(t_time_unique)
            for type_i = 1:size(t_sample_type_unique,2)
                curdataCV = t{ismember(Measured_mets, Measured_mets_unique{met_i}) &...
                    (t.Sample_Time == t_time_unique(time_i)) &...
                    ismember(t.Sample_Type, t_sample_type_unique{1,type_i}),...
                    Mean_tissues};
                curdataCV = cellfun(@(x) str2double(x), curdataCV);
                curdataCV_STD = t{ismember(Measured_mets, Measured_mets_unique{met_i}) &...
                    (t.Sample_Time == t_time_unique(time_i)) &...
                    ismember(t.Sample_Type, t_sample_type_unique{1,type_i}),...
                    STD_tissues};
                curdataCV_STD = cellfun(@(x) str2double(x), curdataCV_STD);
                
                curdataGF = t{ismember(Measured_mets, Measured_mets_unique{met_i}) &...
                    (t.Sample_Time == t_time_unique(time_i)) &...
                    ismember(t.Sample_Type, t_sample_type_unique{2,type_i}),...
                    Mean_tissues};
                curdataGF = cellfun(@(x) str2double(x), curdataGF);
                curdataGF_STD = t{ismember(Measured_mets, Measured_mets_unique{met_i}) &...
                    (t.Sample_Time == t_time_unique(time_i)) &...
                    ismember(t.Sample_Type, t_sample_type_unique{2,type_i}),...
                    STD_tissues};
                curdataGF_STD = cellfun(@(x) str2double(x), curdataGF_STD);

                if ~isempty(curdataCV)
                    meanTable_data(idx,:) = [curdataCV curdataGF];
                    meanTable_dataSTD(idx,:) = [curdataCV_STD curdataGF_STD];
                    meanTable_mets{idx} = strcat(Measured_mets_unique{met_i}, '_',...
                        num2str(t_time_unique(time_i)),'_',...
                        t_sample_type_unique{1,type_i});
                    idx = idx+1;
                end
            end
        end
    end
    meanTable_data(idx:end,:) = [];
    meanTable_dataSTD(idx:end,:) = [];
    meanTable_mets(idx:end) = [];
    
    % remove T12 metabolites
    t0t12metabolites = cellfun(@(x) contains(x, '_12'), meanTable_mets);
    meanTable_data(t0t12metabolites,:) = [];
    meanTable_dataSTD(t0t12metabolites,:) = [];
    
    meanTable_mets(t0t12metabolites) = [];
    
    % save dataset to cell
    meanData_cell{ds_i} = meanTable_data;
    meanDataSTD_cell{ds_i} = meanTable_dataSTD;
    meanMets_cell{ds_i} = meanTable_mets;
end

% join both datasets in one matrix
meanMatrix = vertcat(meanData_cell{:});
meanMatrix_STD = vertcat(meanDataSTD_cell{:});
meanMatrix_mets = horzcat(meanMets_cell{:})';

if use_volume_flag
    applyVolumes = repmat([0.3 0.3 0.3 3 3 3 0.3 0.3 0.3 3 3 3], ...
                            length(meanMatrix_mets), 1);
    % correct volumes for CV mice
    for i=1:length(meanMatrix_mets)
        if contains(meanMatrix_mets{i}, '_CV')
            applyVolumes(i,1:6) = [0.3 0.3 0.3 0.3 0.3 0.3];
        end
    end
    meanMatrix = meanMatrix .* applyVolumes;
    meanMatrix_STD = meanMatrix_STD .* applyVolumes;
end

if (printflag==1)
    % print table to file
    meanData_table = array2table([meanMatrix meanMatrix_STD],...
        'RowNames', meanMatrix_mets,...
        'VariableNames', [strcat(meanData_columns, '_MEAN') strcat(meanData_columns, '_STD')]);
    writetable(meanData_table, ...
        [outputFolder, 'mean_drug_data_for_modelling'],...
        "WriteRowNames",true);
    % print figures to file
    timepoints = zeros(size(meanMatrix_mets));
    metnames = cell(size(meanMatrix_mets));

    for i=1:length(meanMatrix_mets)
        cursplit = strsplit(meanMatrix_mets{i}, '_');
        timepoints(i) = str2double(cursplit{end-1});
        % combine metabolite name with condition
        metnames{i} = strjoin(cursplit([1:length(cursplit)-2, length(cursplit)]),'-');
    end

    % get tissue names
    tissues = cellfun(@(x) strsplit(x, '_'), meanData_columns, 'unif', 0);
    tissues = cellfun(@(x) x{end}, tissues, 'unif', 0);
    % plot kinetic and tissue-distributed data for each metabolite
    metnames_unique = unique(metnames);
    columnsGF = cellfun(@(x) contains(x, 'GF'), meanData_columns);
    
    for i=1:length(metnames_unique)
        curdata = meanMatrix(ismember(metnames, metnames_unique{i}),:);
        curdataSTD = meanMatrix_STD(ismember(metnames, metnames_unique{i}),:);
        curtime = timepoints(ismember(metnames, metnames_unique{i}));
        
        fig = figure('units','normalized','outerposition',[0 0 1 1]);
       
        % plot original kinetic profiles
        for j=1:nnz(columnsGF)
            subplot(2, nnz(columnsGF), j)
            plot(curtime, curdata(:,j), 'k');
            hold on
            plot(curtime, curdata(:, nnz(columnsGF)+j), 'k--')
            errorbar(curtime, curdata(:,j), curdataSTD(:,j), 'k.')
            errorbar(curtime, curdata(:, nnz(columnsGF)+j), curdataSTD(:, nnz(columnsGF)+j), 'k.')
            title(tissues{j})
            xticks(curtime)
            xlabel('Time, h')
            axis square
            ylim([0 (max(max(curdata(:,[j nnz(columnsGF)+j]))) +...
                          max(max(curdataSTD(:,[j nnz(columnsGF)+j]))) )])
        end
        legend({'Colonized', 'GF'})
        % plot metabolites per time point across tissues
        for j=1:length(curtime)
            subplot(2, nnz(columnsGF), nnz(columnsGF)+j)
            plot(1:nnz(columnsGF), curdata(j,~columnsGF), 'k')
            hold on
            plot(1:nnz(columnsGF), curdata(j, columnsGF), 'k--')
            errorbar(1:nnz(columnsGF), curdata(j,~columnsGF),...
                curdataSTD(j, ~columnsGF), '.k')
            errorbar(1:nnz(columnsGF), curdata(j, columnsGF),...
                curdataSTD(j, columnsGF), 'k.')
            title(sprintf('Time %d h', curtime(j)))
            ylim([0 max(max(curdata(j,:))) + max(max(curdataSTD(j,:)))])
            xticks(1:nnz(columnsGF))
            xticklabels(tissues(~columnsGF))
            xlim([1, nnz(columnsGF)])
            axis square
        end
        legend({'Colonized', 'GF'})
        sgtitle(metnames_unique{i})

        % print to file
        print(gcf, '-vector', '-dpdf', '-r600', '-bestfit',...
            [figureFolder, 'figSX_plots_drug_data_for_modelling_',...
            metnames_unique{i}, 'volumeflag_', num2str(use_volume_flag)]);
 
    end
end

% remove t- metabolite sas it does not make sense to model them
 t0t12metabolites = cellfun(@(x) contains(x, '_0_'), meanMatrix_mets);
 meanMatrix(t0t12metabolites,:) = [];
 meanMatrix_STD(t0t12metabolites,:) = [];
 meanMatrix_mets(t0t12metabolites) = [];

% cell(length(modelfilenames),1);
% t_mets = cell(length(modelfilenames),1);
% P_tissues = {'P_si1','P_si2','P_si3',...
%              'P_cecum', 'P_colon', 'P_feces'};
% M_tissues = {'M_cecum', 'M_colon', 'M_feces'};
% t_conditions = [strcat('GF_', ...
%     cellfun(@(x) strrep(x, 'P_', 'Chow1_'), P_tissues, 'unif', 0)),...
%     strcat('BAC_', ...
%     cellfun(@(x) strrep(x, 'P_', 'Chow1_'), P_tissues, 'unif', 0))];
% 
%     
% for files_i = 1:length(datafilenames)
%     datafilename1 = [dataFolder datafilenames1{files_i}];
%     datafilename2 = [dataFolder datafilenames2{files_i}];
%     
%     % build model according to the definition in the table 
%     [modelGutUniversal] = create_model_from_file(modelfilename);
%     % load experimental data
%     [t_GF,metNamesMap_GF, gd_GF, useForFitting_GF] = load_data_from_file(datafilename1);
% 
%     [t_MB,metNamesMap_MB, gd_MB, useForFitting_MB] = load_data_from_file(datafilename2);
%              
%     [~, ~, idxP_GF] = intersect(P_tissues,...
%         t_GF.Properties.VariableNames, 'stable');
%     [~, ~, idxM_GF] = intersect(M_tissues,...
%         t_GF.Properties.VariableNames, 'stable');
%     [~, ~, idxP_MB] = intersect(P_tissues,...
%         t_MB.Properties.VariableNames, 'stable');
%     [~, ~, idxM_MB] = intersect(M_tissues,...
%         t_MB.Properties.VariableNames, 'stable');
%   
%     time3 = find(t_GF.Time==3);
%    
%   
%     figure
%     spidx = 1;
%     for i=time3:size(t_GF,1)
%         subplot(2,3,spidx)
%         plot(t_GF{i,idxP_GF})
%         hold on
%         plot(t_MB{i,idxP_MB})
%         title(sprintf('Time = %d h', t_GF{i,1}))
%         set(gca,'XTick', 1:6)
%         set(gca,'XTickLabel', t_GF.Properties.VariableNames(idxP_GF))
%         spidx = spidx+1;
%         legend({'GF', 'MB'})
%     end
%     suptitle('Parent drug')
%         
%     figure
%     spidx = 1;
%     for i=time3:size(t_GF,1)
%         subplot(2,3,spidx)
%         plot([0 0 0 t_GF{i,idxM_GF}])
%         hold on
%         plot([0 0 0 t_MB{i,idxM_MB}])
%         title(sprintf('Time = %d h', t_GF{i,1}))
%         set(gca,'XTick', 1:6)
%         set(gca,'XTickLabel', t_GF.Properties.VariableNames(idxP_GF))
%         spidx = spidx+1;
%         legend({'GF', 'MB'})
%     end
%     suptitle('Drug metabolite')
% 
%     t{files_i} = [[t_GF{time3:end,idxP_GF};...
%                      [zeros(size(t_GF,1)-time3+1,3)...
%                       t_GF{time3:end,idxM_GF}]],...   
%                   [t_MB{time3:end,idxP_MB};...
%                      [zeros(size(t_MB,1)-time3+1,3)...
%                      t_MB{time3:end,idxM_MB}]]];   
% 
%     
%     t_mets{files_i} = [strcat(strrep(...
%                            strrep(datafilenames1{files_i}, '.csv',''),...
%                            'example_data_',''),'_',...
%                            arrayfun(@(x) num2str(x), t_GF.Time(time3:end),...
%                             'unif', 0));...
%                        strcat('met',strrep(...
%                            strrep(datafilenames1{files_i}, '.csv',''),...
%                            'example_data_',''),'_',...
%                            arrayfun(@(x) num2str(x), t_GF.Time(time3:end),...
%                             'unif', 0))];     
% end