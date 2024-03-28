function assess_model_for_drugs(met_names, met_bestsols, figureFolder)
% assess the model fitting to drugs data 
% based on information of bacterial drug metabolism
% figureFolder is the folder to store plots

% use classification criteria to 
% substrate (-1)
% product (1)
% neither (0)
% of the model coefficients

% corr threshold for reliable solutions
corrthreshold = 0.7;
classthreshold = 0.7;
% make a flag for positive or negative class additional requirement of 
% LI PCC passing the threshold as well
flag_strictclass = 1;

solution_types = met_bestsols{1}.selection_criterion;
classification_reports = cell(length(solution_types)+1,1);
confusion_matrices = cell(length(solution_types)+1,1);

% set expected class labels in case some solutions are missing one
classlabels_expected = [-1; 0; 1];

% try combining solutions
combine_solutions = {'IP', 'LI PCC within high total'};
%combine_solutions = {'IP', 'total PCC'};

for sol_type = 1:length(solution_types)+1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        picksol2 = find(ismember(met_bestsols{1}.selection_criterion,...
                          combine_solutions(2)));              
        model_classes = zeros(size(met_classes));
        model_corr = zeros(size(met_classes));
        model_corrLI = zeros(size(met_classes));
        
        for i=1:length(met_bestsols)
            % if both types of best solutions exist, take the first type as
            % long as it passess corr thresholds, otherwise take the second
            % type
            if size(met_bestsols{i}.x,2)>=max(picksol1,picksol2)
                if (met_bestsols{i}.x_sel_CorrRev(picksol1)>=corrthreshold) &&...
                   (met_bestsols{i}.x_sel_CorrRevLI(picksol1)>=corrthreshold)
                    picksol = picksol1;
                else
                    % second solution type passess both threshlds - take
                    % this one
                    if (met_bestsols{i}.x_sel_CorrRev(picksol2)>=corrthreshold) &&...
                       (met_bestsols{i}.x_sel_CorrRevLI(picksol2)>=corrthreshold)
                        picksol = picksol2;
                    else
                        % first solution passes total corr threshold but
                        % not LI corr threshold
                        if (met_bestsols{i}.x_sel_CorrRev(picksol1)>=corrthreshold) 
                            picksol = picksol1;
                        else
                            % second solution passes the threshold
                            if (met_bestsols{i}.x_sel_CorrRev(picksol2)>=corrthreshold) 
                                picksol = picksol2;
                            else
                                % second solution is not better - keep the
                                % first one
                                picksol = picksol1;
                            end
                        end
                    end
                end
                 
                bestsol = met_bestsols{i}.x(:, picksol);
                model_corr(i) = met_bestsols{i}.selection_value(picksol);
                model_corrLI(i) = met_bestsols{i}.x_sel_CorrRevLI(picksol);
            else
                % only one type sof solution exist for this metabolite,
                % pick this one
                if size(met_bestsols{i}.x,2)>=min(picksol1,picksol2)
                    bestsol = met_bestsols{i}.x(:, min(picksol1,picksol2));
                    model_corr(i) = met_bestsols{i}.selection_value(min(picksol1,picksol2));
                    model_corrLI(i) = met_bestsols{i}.x_sel_CorrRevLI(min(picksol1,picksol2));
                else
                    % both types do not exist - set to 0
                    bestsol = zeros(size(met_bestsols{i}.x,1));
                    model_corr(i) = 0;
                    model_corrLI(i) = 0;
                end
            end
            %normalize without the first coefficient (f)
            bestsol = bestsol(2:end)/max(abs(bestsol(2:end)));
            model_classes(i) = bestsol(end);
        end
    end
    
    % keep only positions that pass the corr threshold
    model_classes(model_corr<corrthreshold)=[];
    met_classes(model_corr<corrthreshold)=[];
    model_corrLI(model_corr<corrthreshold)=[];
    met_names(model_corr<corrthreshold)=[];
    model_corr(model_corr<corrthreshold)=[];
    
    % set model classifiers based on threshold
    model_classes(model_classes<=-classthreshold)=-1;
    model_classes(model_classes>=classthreshold)=1;
    model_classes((model_classes>-classthreshold) &...
                  (model_classes<classthreshold))=0;

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
    [confmat, classlabels] = confusionmat(met_classes,model_classes);
    % check if all expected classes are present
    [~, idx_exp, idx_cur] = intersect(classlabels_expected, classlabels, 'stable');
    if ~isequal(idx_exp', 1:length(classlabels_expected))
        confmat_expected = zeros(length(classlabels_expected));
        confmat_expected(idx_exp, idx_exp) = confmat(idx_cur, idx_cur);
        confmat = confmat_expected;
    end
    % create a matrix with
    % accuracy, precision, recall, specificity, F1, support 
    report_labels = {'accuracy', 'precision', 'recall',...
                     'specificity', 'F1', 'support'}; 
    classification_report = zeros(length(classlabels),length(report_labels));
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
    classification_reports{sol_type} = classification_report.*...
                                            (classification_support>0);
    confusion_matrices{sol_type} = confmat;
end

% add combined solution to solution types
solution_types{end+1} = sprintf('Combined %s + %s', combine_solutions{1},...
                                                combine_solutions{2});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot classification report for different solutions
ST = reshape(repmat(solution_types, length(classlabels_expected),1),[],1); 
C = repmat(arrayfun(@(x) num2str(x), classlabels_expected, 'unif', 0),...
            length(solution_types),1);
combined_report = cell2mat(vertcat(classification_reports));
figure
h = heatmap(combined_report, ...
    'YDisplayLabels', strcat(ST, ': class ', C),...
    'XDisplayLabels', report_labels);
h.CellLabelFormat = '%.2g';
caxis([0 1]) 
title('Classification report')

print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
    [figureFolder, 'figSX_classification_report_drugs',...
    '_class', strrep(num2str(classthreshold),'.','_'),...
    '_corr', strrep(num2str(corrthreshold),'.','_'),...
    '_strictclass', num2str(flag_strictclass)])
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot confusion matrices for different solutions
fig = figure('units','normalized','outerposition',[0 0 1 1]);

spx = 3;
spy = 3;

for sol_type=1:length(solution_types)
    subplot(spx, spy, sol_type)
    confmat = confusion_matrices{sol_type};
    heatmap(confmat, ...
        'YDisplayLabels', classlabels,...
        'XDisplayLabels', classlabels)
    xlabel('Predicted class')
    ylabel('True class')
    title(solution_types{sol_type})
end
% plot one more subplot as the last one appears to be of different size
subplot(spx, spy, sol_type+1)
confmat = confusion_matrices{sol_type};
heatmap(confmat, ...
    'YDisplayLabels', classlabels,...
    'XDisplayLabels', classlabels)
xlabel('Predicted class')
ylabel('True class')
title(solution_types{sol_type})
    
sgtitle('Confusion matrices')
orient landscape

print(fig, '-painters', '-dpdf', '-r600', '-bestfit',...
    [figureFolder, 'figSX_confusion_matrices_classification_drugs',...
    '_class', strrep(num2str(classthreshold),'.','_'),...
    '_corr', strrep(num2str(corrthreshold),'.','_'),...
    '_strictclass', num2str(flag_strictclass)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot all solutions in precision recall and specificity sensitivity
% metrics
marker_shape = {'o', '+', '*', 'x', 's', 'd', '^', 'v'};
marker_color = {'b','k','r'};
marker_color_num = [0 0 1;...
                    0 0 0;...
                    1 0 0];
plotcrit = [3 2;
            4 3;
            1 1
            5 5];

fig = figure('units','normalized','outerposition',[0 0 1 1]);

spx=1;
spy=2;
spidx=1;

subplot(spx, spy, spidx)
spidx=spidx+1;
hold on
for class_i = 1:length(classlabels)
    gscatter(combined_report(class_i:3:end,plotcrit(1,1)),...
             combined_report(class_i:3:end,plotcrit(1,2)),...
        solution_types',...
        repmat(marker_color{class_i},1,length(marker_shape)),...
        strjoin(marker_shape,''))
end
legend('Location','best')
xlim([0 1])
ylim([0 1])
axis square
xlabel(report_labels{plotcrit(1,1)})
ylabel(report_labels{plotcrit(1,2)})

subplot(spx, spy, spidx)
spidx=spidx+1;
hold on
for class_i = 1:length(classlabels)
    gscatter(1-combined_report(class_i:3:end,plotcrit(2,1)),...
             combined_report(class_i:3:end,plotcrit(2,2)),...
        solution_types',...
        repmat(marker_color{class_i},1,length(marker_shape)),...
        strjoin(marker_shape,''),...
        [],[], 'off') % no legend as it is the same as previous subplot
end
xlim([0 1])
ylim([0 1])
axis square
xlabel(report_labels{plotcrit(2,1)})
ylabel(report_labels{plotcrit(2,2)})

orient landscape

print(fig, '-painters', '-dpdf', '-r600', '-bestfit',...
    [figureFolder, 'figSX_selection_criteria_classification_drugs',...
    '_class', strrep(num2str(classthreshold),'.','_'),...
    '_corr', strrep(num2str(corrthreshold),'.','_'),...
    '_strictclass', num2str(flag_strictclass)])

fig = figure('units','normalized','outerposition',[0 0 1 1]);

spx=2;
spy=2;
spidx=1;

subplot(spx, spy, spidx)
spidx=spidx+1;
b = barh([combined_report(1:3:end,1)...
     combined_report(2:3:end,1)...
     combined_report(3:3:end,1)],...
     'FaceColor','flat');
for k = 1:length(marker_color)
    b(k).CData = marker_color_num(k,:);
end
xlabel(report_labels{plotcrit(3,1)})
set(gca, 'YTick', 1:length(solution_types))
set(gca, 'YTickLabel', solution_types)
axis square
legend(C(1:3), 'Location', 'best')

subplot(spx, spy, spidx)
spidx=spidx+1;
plotdata = mean([combined_report(1:3:end,1)...
     combined_report(2:3:end,1)...
     combined_report(3:3:end,1)],2);
[plotdata, plotsort] = sort(plotdata, 'ascend'); 
b = barh(plotdata,...
     'FaceColor','flat');
xlabel(['Mean ' report_labels{plotcrit(3,1)}])
set(gca, 'YTick', 1:length(solution_types))
set(gca, 'YTickLabel', solution_types(plotsort))
axis square

subplot(spx, spy, spidx)
spidx=spidx+1;
b = barh([combined_report(1:3:end,5)...
     combined_report(2:3:end,5)...
     combined_report(3:3:end,5)],...
     'FaceColor','flat');
for k = 1:length(marker_color)
    b(k).CData = marker_color_num(k,:);
end
xlabel(report_labels{plotcrit(4,1)})
set(gca, 'YTick', 1:length(solution_types))
set(gca, 'YTickLabel', solution_types)
axis square
legend(C(1:3), 'Location', 'best')

subplot(spx, spy, spidx)
spidx=spidx+1;
plotdata = mean([combined_report(1:3:end,5)...
     combined_report(2:3:end,5)...
     combined_report(3:3:end,5)],2);
[plotdata, plotsort] = sort(plotdata, 'ascend'); 
b = barh(plotdata,...
     'FaceColor','flat');
xlabel(['Mean ' report_labels{plotcrit(4,1)}])
set(gca, 'YTick', 1:length(solution_types))
set(gca, 'YTickLabel', solution_types(plotsort))
xlim([0 1])
axis square

orient landscape

print(fig, '-painters', '-dpdf', '-r600', '-bestfit',...
    [figureFolder, 'figSX_Accuracy_F1_classification_drugs',...
    '_class', strrep(num2str(classthreshold),'.','_'),...
    '_corr', strrep(num2str(corrthreshold),'.','_'),...
    '_strictclass', num2str(flag_strictclass)])
