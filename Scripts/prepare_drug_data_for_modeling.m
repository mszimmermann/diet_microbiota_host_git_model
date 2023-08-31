%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prepare data from BRV paper for modeling
dataFolder = 'C:\Users\mazimmer\Documents\Desktop\MikeBrivudine\REVISIONS\';
dataFilenames = 'Brivudine_Supp_Tables_V26.xlsx';
dataSheetnames = {'Table S4';'Table S17'};

meanData_cell = cell(size(dataSheetnames));
meanMets_cell = cell(size(dataSheetnames));
Mean_tissues = {'Mean_SI','Mean_SII','Mean_SIII',...
                 'Mean_Cecum', 'Mean_Colon', 'Mean_Feces'};
meanData_columns = [strcat('CV_Chow_', Mean_tissues)...
                    strcat('GF_Chow_', Mean_tissues)];

for ds_i = 1:length(dataSheetnames)
    dataSheetname = dataSheetnames{ds_i};
              
    % read drug data
    t = readtable([dataFolder dataFilenames], 'Sheet', dataSheetname);
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
    % remove Thymin as it is measured only in liver
    Measured_mets_unique(cellfun(@(x) contains(x, 'Thymine'), Measured_mets_unique))=[];
    t_time_unique = unique(t.Sample_Time);
    % remove time points between 0 and 3
    t_time_unique(t_time_unique>0 & t_time_unique<3)=[];
    t_sample_type_unique = reshape(unique(t.Sample_Type), 2,[]);


    meanTable_data = zeros(length(Measured_mets_unique) * length(t_time_unique)*...
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
                curdataGF = t{ismember(Measured_mets, Measured_mets_unique{met_i}) &...
                    (t.Sample_Time == t_time_unique(time_i)) &...
                    ismember(t.Sample_Type, t_sample_type_unique{2,type_i}),...
                    Mean_tissues};
                curdataGF = cellfun(@(x) str2double(x), curdataGF);
                if ~isempty(curdataCV)
                    meanTable_data(idx,:) = [curdataCV curdataGF];
                    meanTable_mets{idx} = strcat(Measured_mets_unique{met_i}, '_',...
                        num2str(t_time_unique(time_i)),'_',...
                        t_sample_type_unique{1,type_i});
                    idx = idx+1;
                end
            end
        end
    end
    meanTable_data(idx:end,:) = [];
    meanTable_mets(idx:end) = [];
    
    meanData_cell{ds_i} = meanTable_data;
    meanMets_cell{ds_i} = meanTable_mets;
end

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