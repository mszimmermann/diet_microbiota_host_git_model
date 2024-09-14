function [solution_types, met_names_predicted, ...
    classification_reports, confusion_matrices] = ...
    assess_model_for_drugs_per_drug(met_names_full, met_bestsols, flag_strictclass, ...
                                    figureFolder, outputFolder)
% assess the model fitting to drugs data 
% based on information of bacterial drug metabolism
% figureFolder is the folder to store plots

% use classification criteria to 
% substrate (-1)
% product (1)
% neither (0)
% of the model coefficients

% flag_strictclass determines if calling positive or negative class has to pass an
% additional requirement of LI PCC passing the threshold as well

% corr threshold for reliable solutions
corrthreshold = 0.7;
classthreshold = 0.7;


solution_types = met_bestsols{1}.selection_criterion;
% make all combinations with IP
classification_reports = cell(length(solution_types)*2-1,1);
confusion_matrices = cell(length(solution_types)*2-1,1);
met_names_predicted = cell(length(solution_types)*2-1,1);

met_all_correlations = zeros(length(met_names_full), 2*length(solution_types)-1);
met_all_correlationsLI = zeros(length(met_names_full), 2*length(solution_types)-1);
met_all_classes = zeros(length(met_names_full), 2*length(solution_types)-1);

% set expected class labels in case some solutions are missing one
classlabels_expected = [-1; 0; 1];

% try combining solutions
combine_solutions = {'IP', 'LI PCC within high total'};
%combine_solutions = {'IP', 'total PCC'};
% calculate accuracy per drug
drug_types = {'BRV', 'SRV', 'CLZ'};
for sol_type = 1:2*length(solution_types)-1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get all metabolite names for the current solution type
    met_names = met_names_full;
    % set classes for drugs and drug metabolites
    met_classes = zeros(size(met_names));
    % BRV is a substrate unless colonization is d4554
    met_classes(contains(met_names, {'BRV_BRV'}) &... 
                ~contains(met_names,{'d4554'})) = -1;
    % SRV is a substrate
    met_classes(contains(met_names, {'SRV_SRV'}))= -1;
    % BRV_BVU and SRV_BVU is a product unless colonization is d4554
    met_classes(contains(met_names, {'BRV_BVU', 'SRV_BVU'}) &... 
                ~contains(met_names,{'d4554'})) = 1;
    % CLZ and 3-OH-CLZ is a substrate
    met_classes(contains(met_names, {'CLZ', '3-OH-CLZ'}) &... 
                ~contains(met_names,{'NH2', 'NH2OH'})) = -1;
    % NH2-CLZ and NH2OH-CLZ is a product
    met_classes(contains(met_names,{'NH2', 'NH2OH'})) = 1;

    % save classes for all input mets
    met_names_full_classes = met_classes;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get classification by the model
    if sol_type<=length(solution_types)
        picksol = find(ismember(met_bestsols{1}.selection_criterion,...
                          solution_types(sol_type)));
        model_classes = zeros(size(met_classes));
        model_corr = zeros(size(met_classes));
        model_corrLI = zeros(size(met_classes));
        for i=1:length(met_bestsols)
            if size(met_bestsols{i}.x,2)>=picksol
                bestsol = met_bestsols{i}.x(:, picksol);
                model_corr(i) = met_bestsols{i}.selection_value(picksol);
                model_corrLI(i) = met_bestsols{i}.x_sel_CorrRevLI(picksol);
            else
                bestsol = zeros(size(met_bestsols{i}.x,1));
                model_corr(i) = 0;
                model_corrLI(i) = 0;
            end
            %normalize without the first coefficient (f)
            bestsol = bestsol(2:end)/max(abs(bestsol(2:end)));
            model_classes(i) = bestsol(end);
        end
    else
        picksol1 = find(ismember(met_bestsols{1}.selection_criterion,...
                          combine_solutions(1)));
        %picksol2 = find(ismember(met_bestsols{1}.selection_criterion,...
        %                  combine_solutions(2)));  
        picksol2 = sol_type - length(solution_types)+1;
        model_classes = zeros(size(met_classes));
        model_corr = zeros(size(met_classes));
        model_corrLI = zeros(size(met_classes));
        
        for i=1:length(met_bestsols)
            % combine two solutions
            [bestsol, picksol] = combine_bestsols(met_bestsols{i},...
                                                  picksol1, picksol2,...
                                                  corrthreshold);
            if picksol==0
                model_corr(i) = 0;
                model_corrLI(i) = 0;
            else
                model_corr(i) = met_bestsols{i}.x_sel_CorrRev(picksol);
                model_corrLI(i) = met_bestsols{i}.x_sel_CorrRevLI(picksol);
            end
            %normalize without the first coefficient (f)
            bestsol = bestsol(2:end)/max(abs(bestsol(2:end)));
            model_classes(i) = bestsol(end);
        end
    end

    % set model classifiers based on threshold
    model_classes(model_classes<=-classthreshold)=-1;
    model_classes(model_classes>=classthreshold)=1;
    model_classes((model_classes>-classthreshold) &...
                  (model_classes<classthreshold))=0;

    % save correlations and predictions for all metabolites
    met_all_correlations(:, sol_type) = model_corr;
    met_all_correlationsLI(:, sol_type) = model_corrLI;
    met_all_classes(:, sol_type) = model_classes;
    
    % keep only positions that pass the corr threshold
    model_classes(model_corr<corrthreshold)=[];
    met_classes(model_corr<corrthreshold)=[];
    model_corrLI(model_corr<corrthreshold)=[];
    met_names(model_corr<corrthreshold)=[];
    model_corr(model_corr<corrthreshold)=[];
    
    
    if flag_strictclass
        % keep only positions that pass the LI corr threshold
        % for nozero classes
        unreliable_classes = (model_classes~=0) & (model_corrLI<corrthreshold);
        model_classes(unreliable_classes)=[];
        met_classes(unreliable_classes)=[];
        model_corrLI(unreliable_classes)=[];
        model_corr(unreliable_classes)=[];
        met_names(unreliable_classes)=[];
    end
    % get the confusion matrix  
    %Brivudine (BRV)
    select_drug = cellfun(@(x) contains(x,'BRV'), met_names);
    [confmat_BRV, classlabels_BRV] = confusionmat(met_classes(select_drug),...
                                          model_classes(select_drug));
    % check if all expected classes are present
    [~, idx_exp, idx_cur] = intersect(classlabels_expected, classlabels_BRV, 'stable');
    if ~isequal(idx_exp', 1:length(classlabels_expected))
        confmat_expected = zeros(length(classlabels_expected));
        confmat_expected(idx_exp, idx_exp) = confmat_BRV(idx_cur, idx_cur);
        confmat_BRV = confmat_expected;
    end
    % Sorivudine SRV
    select_drug = cellfun(@(x) contains(x,'SRV'), met_names);
    [confmat_SRV, classlabels_SRV] = confusionmat(met_classes(select_drug),...
                                          model_classes(select_drug));
    % check if all expected classes are present
    [~, idx_exp, idx_cur] = intersect(classlabels_expected, classlabels_SRV, 'stable');
    if ~isequal(idx_exp', 1:length(classlabels_expected))
        confmat_expected = zeros(length(classlabels_expected));
        confmat_expected(idx_exp, idx_exp) = confmat_SRV(idx_cur, idx_cur);
        confmat_SRV = confmat_expected;
    end
    % Clonazepan CLZ
    select_drug = cellfun(@(x) contains(x,'CLZ'), met_names);
    [confmat_CZL, classlabels_CZL] = confusionmat(met_classes(select_drug),...
                                                  model_classes(select_drug));
    % check if all expected classes are present
    [~, idx_exp, idx_cur] = intersect(classlabels_expected, classlabels_CZL, 'stable');
    if ~isequal(idx_exp', 1:length(classlabels_expected))
        confmat_expected = zeros(length(classlabels_expected));
        confmat_expected(idx_exp, idx_exp) = confmat_CZL(idx_cur, idx_cur);
        confmat_CZL = confmat_expected;
    end
    % create a matrix with
    % accuracy, precision, recall, specificity, F1, support 
    report_labels = {'accuracy', 'precision', 'recall',...
                     'specificity', 'F1', 'support'}; 
    
    classification_report_drug = cell(3,1);
    for drug_type = 1:length(drug_types)
        if drug_types{drug_type} == 'BRV'
            confmat = confmat_BRV;
        elseif drug_types{drug_type} == 'SRV'
            confmat = confmat_SRV;
        else
            confmat = confmat_CZL;
        end
        classification_report = zeros(length(classlabels_expected),length(report_labels));
        for i = 1:length(classlabels_expected)
            TP = confmat(i,i);
            FP = sum(confmat(:, i), 1) - TP;
            FN = sum(confmat(i, :), 2) - TP;
            TN = sum(confmat(:)) - TP - FP - FN;
    
            Accuracy = (TP+TN)./(TP+FP+TN+FN);
            classification_report(i,1) = Accuracy;
    
            PPV = TP./ (TP + FP); % tp / predicted positive PRECISION
            if isnan(PPV)
                PPV = 0;
            end
            classification_report(i,2) = PPV;
    
            TPR = TP./(TP + FN);%tp/actual positive  RECALL SENSITIVITY
            if isnan(TPR)
                TPR = 0;
            end
            classification_report(i,3) = TPR;
    
            TNR = TN./ (TN+FP); %tn/ actual negative  SPECIFICITY
            if isnan(TNR)
                TNR = 0;
            end
            classification_report(i,4) = TNR;
    
        %     FPR = FP./ (TN+FP);
        %     if isnan(FPR)
        %         FPR = 0;
        %     end
            FScore = (2*(PPV * TPR)) / (PPV+TPR);
            if isnan(FScore)
                FScore = 0;
            end
            classification_report(i,5) = FScore;
    
            % number of class instances
            classification_report(i,6) = sum(confmat(i,:));
    
        end
        % set rows to zero if support is 0 (no elements of tis class)
        classification_support = repmat(classification_report(:,end),1,...
                                    size(classification_report,2));
        classification_report_drug{drug_type} = classification_report.*...
                                            (classification_support>0);
    end
    % save classificatin report and confusion matrices
    classification_reports{sol_type} = classification_report_drug;
    confusion_matrices{sol_type} = {confmat_BRV, confmat_SRV, confmat_CZL};
    % save names of metabolites that were confidently predicted
    met_names_predicted{sol_type} = met_names;
end

% add combined solution to solution types
%solution_types{end+1} = sprintf('Combined %s + %s', combine_solutions{1},...
%                                                combine_solutions{2});
for i=2:length(met_bestsols{1}.selection_criterion)
    solution_types{end+1} = sprintf('Combined %s + %s', combine_solutions{1},...
                                    met_bestsols{1}.selection_criterion{i});
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot confusion matrices for different solutions

spx = 3;%length(drug_types);
spy = 5;%length(solution_types);


for drug_type = 1:length(drug_types)    
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    
    for sol_type=1:length(solution_types)
        
        %subplot(spx, spy, (drug_type-1)*length(solution_types)+sol_type)
        subplot(spx, spy, sol_type)
        confmat = confusion_matrices{sol_type}{drug_type};
        classrep = classification_reports{sol_type}{drug_type};
        F1average = 0;
        for class_i=1:size(classrep,1)
            F1average = F1average + classrep(class_i, 5)*classrep(class_i,6);
        end
        F1average = F1average/sum(classrep(:,6));
        heatmap(confmat, ...
            'YDisplayLabels', classlabels_expected,...
            'XDisplayLabels', classlabels_expected)
        xlabel('Predicted class')
        ylabel('True class')
        title({solution_types{sol_type}, ...
            sprintf('%s acc=%.3f F1av=%.3f', drug_types{drug_type},...
                                sum(diag(confmat))/sum(sum(confmat)),...
                                F1average)})
        
    end
    sgtitle(sprintf('Confusion matrices for %s', drug_types{drug_type}))
    orient landscape
    
    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);
    
    print(fig, '-painters', '-dpdf', '-r600', '-bestfit',...
        [figureFolder, 'figSX_confusion_matrices_classification_drugs_per_drug_',...
        drug_types{drug_type},...
        '_class', strrep(num2str(classthreshold),'.','_'),...
        '_corr', strrep(num2str(corrthreshold),'.','_'),...
        '_strictclass', num2str(flag_strictclass)])

end

%plot confusion matrices for all drugs together
spx = 3;%length(drug_types);
spy = 5;%length(solution_types);

fig = figure('units','normalized','outerposition',[0 0 1 1]);
    
for sol_type=1:length(solution_types)
    
    subplot(spx, spy, sol_type)
    confmat = zeros(size(confusion_matrices{sol_type}{1}));
    for drug_type = 1:length(drug_types)    
        confmat = confmat + confusion_matrices{sol_type}{drug_type};
    end
    classrep = classification_reports{sol_type}{drug_type};
    F1average = 0;
    for class_i=1:size(classrep,1)
        TP = confmat(class_i,class_i);
        FP = sum(confmat(:, class_i), 1) - TP;
        FN = sum(confmat(class_i, :), 2) - TP;
        
        PPV = TP./ (TP + FP); % tp / predicted positive PRECISION
        if isnan(PPV)
            PPV = 0;
        end
      
        TPR = TP./(TP + FN);%tp/actual positive  RECALL SENSITIVITY
        if isnan(TPR)
            TPR = 0;
        end
     
        FScore = (2*(PPV * TPR)) / (PPV+TPR);
        F1average = F1average + FScore*sum(confmat(class_i,:));
    end
    F1average = F1average/sum(sum(confmat));
    heatmap(confmat, ...
        'YDisplayLabels', classlabels_expected,...
        'XDisplayLabels', classlabels_expected)
    xlabel('Predicted class')
    ylabel('True class')
    title({solution_types{sol_type}, ...
        sprintf('All drugs acc=%.3f F1av=%.3f', sum(diag(confmat))/sum(sum(confmat)),...
                            F1average)})
    
end

orient landscape

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);

print(fig, '-painters', '-dpdf', '-r600', '-bestfit',...
    [figureFolder, 'figSX_confusion_matrices_classification_all_combos_drugs_all_drugs_',...
    '_class', strrep(num2str(classthreshold),'.','_'),...
    '_corr', strrep(num2str(corrthreshold),'.','_'),...
    '_strictclass', num2str(flag_strictclass)])


    
% sgtitle('Confusion matrices')
% orient landscape
% 
% set(gcf,'Units','inches');
% screenposition = get(gcf,'Position');
% set(gcf,...
%     'PaperPosition',[0 0 screenposition(3:4)],...
%     'PaperSize',[screenposition(3:4)]);
% 
% print(fig, '-painters', '-dpdf', '-r600', '-bestfit',...
%     [figureFolder, 'figSX_confusion_matrices_classification_drugs_per_drug',...
%     '_class', strrep(num2str(classthreshold),'.','_'),...
%     '_corr', strrep(num2str(corrthreshold),'.','_'),...
%     '_strictclass', num2str(flag_strictclass)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot scores vs fraction of metabolites for which classification was
% possible
support_idx = 6;
F1score_idx = 5;
marker_shape = {'o', '+', '*', 'x', 's', 'd', '^', 'v'};
marker_color = {'b','k','r'};
marker_color_num = [0 0 1;...
                    0 0 0;...
                    1 0 0];
fig = figure('units','normalized','outerposition',[0 0 1 1]);

for drug_type = 1:length(drug_types)
    for sol_type = 1:size(classification_reports, 1)
        curmets_total = contains(met_names_full, drug_types{drug_type});
        curclassrep = classification_reports{sol_type}{drug_type};
      
        curscore_weighted_avg = 0;
      
        for cur_label = 1:length(classlabels_expected)
            nmets_total = nnz(met_names_full_classes(curmets_total)==classlabels_expected(cur_label));
            cursupport = curclassrep(cur_label, support_idx);
            curscore = curclassrep(cur_label, F1score_idx);
            
            % calculate weighted average score
            curscore_weighted_avg = curscore_weighted_avg + curscore*cursupport;

            % plot per drug
            subplot(2, length(drug_types), drug_type)
            hold on
            
            if sol_type <= length(met_bestsols{1}.selection_criterion)
                scatter(cursupport/nmets_total, curscore,...
                        20,...
                       marker_color{cur_label},...
                       marker_shape{sol_type})
             else
             scatter(cursupport/nmets_total, curscore,...
                       60,...
                       marker_color{cur_label},...
                       marker_shape{sol_type-length(met_bestsols{1}.selection_criterion)+1})
             
             xlim([0 1])
             ylim([0 1])
             axis square 
             xlabel(sprintf('Fraction metabolite covered (out of n=%d)', nnz(curmets_total)))
             ylabel('Average weighted F1 score')
             title({drug_types{drug_type}, sprintf('b(%d) k(%d) r(%d), size: 20(single) 50(combined)',classlabels_expected)})
            end
        end
        % plot class-weighted F1 score
        cursupport = sum(curclassrep(:, support_idx))/nnz(curmets_total);
        curscore_weighted_avg = curscore_weighted_avg/sum(curclassrep(:,support_idx));    
        % plot per drug
        subplot(2, length(drug_types), length(drug_types) + drug_type)
        hold on
        
        %fprintf('%s %s %.2f\n', drug_types{drug_type}, solution_types{sol_type}, curscore_weighted_avg)
        if sol_type <= length(met_bestsols{1}.selection_criterion)
            scatter(cursupport, curscore_weighted_avg,...
               'k',...
               marker_shape{sol_type})
        else
             scatter(cursupport, curscore_weighted_avg,...
               'r',...
               marker_shape{sol_type-length(met_bestsols{1}.selection_criterion)+1})
        end
       
    end
    xlabel(sprintf('Fraction metabolite covered (out of n=%d)', nnz(curmets_total)))
    ylabel('Average weighted F1 score')
    title(drug_types{drug_type})
    xlim([0 1])
    ylim([0 1])
    axis square
%    xlabel(report_labels{plotcrit(1,1)})    ylabel(report_labels{plotcrit(1,2)})
end
legend(solution_types, 'Location','best')

orient landscape
    
print(fig, '-painters', '-dpdf', '-r600', '-bestfit',...
    [figureFolder, 'figSX_F1score_vs_support_per_drug_',...
    '_class', strrep(num2str(classthreshold),'.','_'),...
    '_corr', strrep(num2str(corrthreshold),'.','_'),...
    '_strictclass', num2str(flag_strictclass)])

% save to file predictions, scores and confusion matrices
% save metabolite class predictions
curtable = array2table(met_all_classes,...
                      'RowNames', met_names_full,...
                      'VariableNames', solution_types);
writetable(curtable, ...
    [outputFolder, 'table_drug_metabolite_classes',...
    '_class', strrep(num2str(classthreshold),'.','_'),...
    '_corr', strrep(num2str(corrthreshold),'.','_'),...
    '_strictclass', num2str(flag_strictclass)],...
    'WriteRowNames',true,...
    'FileType', 'text');
% save metabolite reciprocal correlations
curtable = array2table(met_all_correlations,...
                      'RowNames', met_names_full,...
                      'VariableNames', solution_types);
writetable(curtable, ...
    [outputFolder, 'table_drug_metabolite_reciprocal_correlations',...
    '_class', strrep(num2str(classthreshold),'.','_'),...
    '_corr', strrep(num2str(corrthreshold),'.','_'),...
    '_strictclass', num2str(flag_strictclass)],...
    'WriteRowNames',true,...
    'FileType', 'text');
curtable = array2table(met_all_correlationsLI,...
                      'RowNames', met_names_full,...
                      'VariableNames', solution_types);
writetable(curtable, ...
    [outputFolder, 'table_drug_metabolite_reciprocal_correlationsLI',...
    '_class', strrep(num2str(classthreshold),'.','_'),...
    '_corr', strrep(num2str(corrthreshold),'.','_'),...
    '_strictclass', num2str(flag_strictclass)],...
    'WriteRowNames',true,...
    'FileType', 'text');
% save classification reports
fileName = [outputFolder, 'table_drug_classification_reports',...
    '_class', strrep(num2str(classthreshold),'.','_'),...
    '_corr', strrep(num2str(corrthreshold),'.','_'),...
    '_strictclass', num2str(flag_strictclass),'.csv'];
fid = fopen(fileName, 'w');
colnames = [{'Solution', 'Drug', 'Class'}, report_labels];
for i=1:length(colnames)
    if i<length(colnames)
        fprintf(fid, '%s,', colnames{i});
    else
        fprintf(fid, '%s\n', colnames{i});
    end
end
for i=1:length(classification_reports)
    for j=1:length(classification_reports{i})
        currep = classification_reports{i};
        currep = currep{j};
        for k=1:size(currep,1)
            fprintf(fid, '%s,%s,%d', solution_types{i},...
                drug_types{j}, classlabels_expected(k));
            for m=1:size(currep,2)
                fprintf(fid, ',%.3f', currep(k,m));
            end
            fprintf(fid, '\n');
        end
    end
end
fclose(fid);

% save confusion matrices
fileName = [outputFolder, 'table_drug_confusion_matrices',...
    '_class', strrep(num2str(classthreshold),'.','_'),...
    '_corr', strrep(num2str(corrthreshold),'.','_'),...
    '_strictclass', num2str(flag_strictclass), '.csv'];
fid = fopen(fileName, 'w');
colnames = [{'Solution', 'Drug', 'Class'}, ...
    arrayfun(@(x) ['Predicted ' num2str(x)], classlabels_expected, 'UniformOutput', false)'];
for i=1:length(colnames)
    if i<length(colnames)
        fprintf(fid, '%s,', colnames{i});
    else
        fprintf(fid, '%s\n', colnames{i});
    end
end
for i=1:length(confusion_matrices)
    for j=1:length(confusion_matrices{i})
        curmat = confusion_matrices{i};
        curmat = curmat{j};
        for k=1:size(curmat,1)
            fprintf(fid, '%s,%s,%d', solution_types{i},...
                drug_types{j}, classlabels_expected(k));
            for m=1:size(curmat,2)
                fprintf(fid, ',%.3f', curmat(k,m));
            end
            fprintf(fid, '\n');
        end
    end
end
fclose(fid);

% print metabolite classes to file
met_classes_table = table(met_names_full_classes,...
    'RowNames',met_names_full,...
    'VariableNames', {'Met_class'});
writetable(met_classes_table, [outputFolder, 'drug_metabolite_class_table',...
    '_class', strrep(num2str(classthreshold),'.','_'),...
    '_corr', strrep(num2str(corrthreshold),'.','_'),...
    '_strictclass', num2str(flag_strictclass)],...
    "WriteRowNames",true);
