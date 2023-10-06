function [confmat, classlabels] = compare_two_model_results(met_names1, met_bestsols1,...
                                   met_names2, met_bestsols2,...
                                   figureFolder)
% assess two modelling results (e.g. to two datasets,
% or two models to the same dataset)
% figureFolder is the folder to store plots

% use classification criteria to 
% substrate (-1)
% product (1)
% neither (0)
% of the model coefficients

% corr threshold for reliable solutions
corrthreshold = 0.7;
classthreshold = 0.5;%1;
% make a flag for positive or negative class additional requirement of 
% LI PCC passing the threshold as well
flag_strictclass = 1;%0;%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% intersect modelled metabolite names
[met_names, idx1, idx2] = intersect(lower(met_names1.CompoundName),...
                                    lower(met_names2.CompoundName));
% set classes for each model
% get classification by each model
model1_classes = zeros(length(met_names),1);
model2_classes = zeros(length(met_names),1);
% record also the "sources" - diet, host or microbe
model1_sources = zeros(length(met_names),1);
model2_sources = zeros(length(met_names),1);

if length(met_bestsols1.coefvalues)>5
    model1_classes = zeros(length(met_names),2);
    model1_sources = zeros(length(met_names),2);
end
if length(met_bestsols2.coefvalues)>5
    model2_classes = zeros(length(met_names),2);
    model2_sources = zeros(length(met_names),2);
end

model1_corr = met_bestsols1.x_sel_CorrRev(idx1);
model1_corrLI = met_bestsols1.x_sel_CorrRevLI(idx1);
for i=1:length(met_names)
    bestsol = met_bestsols1.x(idx1(i),:);
    %normalize without the first coefficient (f)
    bestsol = bestsol(2:end)/max(abs(bestsol(2:end)));
    if length(met_bestsols1.coefvalues)>5
        model1_classes(i,:) = bestsol(end-1:end);
    else
        model1_classes(i) = bestsol(end);
    end
    model1_sources(i) = determine_met_source(bestsol);
end

model2_corr = met_bestsols2.x_sel_CorrRev(idx2);
model2_corrLI = met_bestsols2.x_sel_CorrRevLI(idx2);
for i=1:length(met_names)
    bestsol = met_bestsols2.x(idx2(i),:);
    %normalize without the first coefficient (f)
    bestsol = bestsol(2:end)/max(abs(bestsol(2:end)));
    if length(met_bestsols2.coefvalues)>5
        model2_classes(i,:) = bestsol(end-1:end);
    else
        model2_classes(i) = bestsol(end);
    end
    model2_sources(i) = determine_met_source(bestsol);
end

% keep only metabolites that pass the corr threshold
select_mets = (model1_corr >= corrthreshold) & ...
              (model2_corr >= corrthreshold); 
model1_classes = model1_classes(select_mets,:);
model2_classes = model2_classes(select_mets,:);
model1_corr = model1_corr(select_mets,:);
model1_corrLI = model1_corrLI(select_mets,:);
model2_corr = model2_corr(select_mets,:);
model2_corrLI = model2_corrLI(select_mets,:);
met_names = met_names(select_mets);

% set model classifiers based on threshold
model1_classes(model1_classes<=-classthreshold)=-1;
model1_classes(model1_classes>=classthreshold)=1;
model1_classes((model1_classes>-classthreshold) &...
              (model1_classes<classthreshold))=0;

model2_classes(model2_classes<=-classthreshold)=-1;
model2_classes(model2_classes>=classthreshold)=1;
model2_classes((model2_classes>-classthreshold) &...
              (model2_classes<classthreshold))=0;

          
model1_classes = model1_classes(:,1);

if flag_strictclass
    % keep only positions that pass the LI corr threshold
    % for nozero classes
    unreliable_classes = ((model1_classes~=0) & (model1_corrLI<corrthreshold)) |...
                         ((model2_classes~=0) & (model2_corrLI<corrthreshold));
    model1_classes(unreliable_classes)=[];
    model1_corrLI(unreliable_classes)=[];
    model1_corr(unreliable_classes)=[];
    model2_classes(unreliable_classes)=[];
    model2_corrLI(unreliable_classes)=[];
    model2_corr(unreliable_classes)=[];
end
    
          
% get the confusion matrix
[confmat, classlabels] = confusionmat(model1_classes(:,1),...
                                      model2_classes(:,1));
                                  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot confusion matrix
fig = figure;%('units','normalized','outerposition',[0 0 1 1]);

heatmap(confmat, ...
        'YDisplayLabels', classlabels,...
        'XDisplayLabels', classlabels)
xlabel([met_bestsols2.modelname ' class'])
ylabel([met_bestsols1.modelname ' class'])
suptitle('Confusion matrix')
orient landscape

print(fig, '-painters', '-dpdf', '-r600', '-bestfit',...
    [figureFolder, 'figSX_confusion_matrix_'...
    met_bestsols2.modelname, '_vs_',...
    met_bestsols1.modelname,...
    '_class', strrep(num2str(classthreshold),'.','_'),...
    '_corr', strrep(num2str(corrthreshold),'.','_'),...
    '_strictclass', num2str(flag_strictclass)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the confusion matrix for metabolite sources
% select_mets = (model1_corr >= corrthreshold) & ...
%               (model2_corr >= corrthreshold); 
% model1_sources = model1_sources(select_mets,:);
% model2_sources = model2_sources(select_mets,:);
% 
% [confmat, classlabels] = confusionmat(model1_sources(:,1),...
%                                       model2_sources(:,1));
% 
% heatmap(confmat, ...
%         'YDisplayLabels', classlabels,...
%         'XDisplayLabels', classlabels)
% xlabel([met_bestsols2.modelname ' source'])
% ylabel([met_bestsols1.modelname ' source'])
% suptitle('Confusion matrix')
% orient landscape
% 
% % create a matrix with 
% % accuracy, precision, recall, specificity, F1, support 
% report_labels = {'accuracy', 'precision', 'recall',...
%                  'specificity', 'F1', 'support'}; 
% classification_report = zeros(length(classlabels),length(report_labels));
% for i = 1:length(classlabels)
%     TP = confmat(i,i);
%     FP = sum(confmat(:, i), 1) - TP;
%     FN = sum(confmat(i, :), 2) - TP;
%     TN = sum(confmat(:)) - TP - FP - FN;
% 
%     Accuracy = (TP+TN)./(TP+FP+TN+FN);
%     classification_report(i,1) = Accuracy;
% 
%     PPV = TP./ (TP + FP); % tp / predicted positive PRECISION
%     if isnan(PPV)
%         PPV = 0;
%     end
%     classification_report(i,2) = PPV;
% 
%     TPR = TP./(TP + FN);%tp/actual positive  RECALL SENSITIVITY
%     if isnan(TPR)
%         TPR = 0;
%     end
%     classification_report(i,3) = TPR;
% 
%     TNR = TN./ (TN+FP); %tn/ actual negative  SPECIFICITY
%     if isnan(TNR)
%         TNR = 0;
%     end
%     classification_report(i,4) = TNR;
% 
% %     FPR = FP./ (TN+FP);
% %     if isnan(FPR)
% %         FPR = 0;
% %     end
%     FScore = (2*(PPV * TPR)) / (PPV+TPR);
%     if isnan(FScore)
%         FScore = 0;
%     end
%     classification_report(i,5) = FScore;
% 
%     % number of class instances
%     classification_report(i,6) = sum(confmat(i,:));
% 
% end
% classification_reports{sol_type} = classification_report;
% confusion_matrices{sol_type} = confmat;
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % plot classification report for different solutions
% ST = reshape(repmat(solution_types, length(classlabels),1),[],1); 
% C = repmat(arrayfun(@(x) num2str(x), classlabels, 'unif', 0),length(solution_types),1);
% combined_report = cell2mat(vertcat(classification_reports));
% figure
% heatmap(combined_report, ...
%     'YDisplayLabels', strcat(ST, ': class ', C),...
%     'XDisplayLabels', report_labels)
% caxis([0 1]) 
% title('Classification report')
% 
% print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
%     [figureFolder, 'figSX_classification_report_drugs',...
%     '_class', strrep(num2str(classthreshold),'.','_'),...
%     '_corr', strrep(num2str(corrthreshold),'.','_')])
%     
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % plot confusion matrices for different solutions
% fig = figure('units','normalized','outerposition',[0 0 1 1]);
% 
% spx = 3;
% spy = 3;
% 
% for sol_type=1:length(solution_types)
%     subplot(spx, spy, sol_type)
%     confmat = confusion_matrices{sol_type};
%     heatmap(confmat, ...
%         'YDisplayLabels', classlabels,...
%         'XDisplayLabels', classlabels)
%     xlabel('Predicted class')
%     ylabel('True class')
%     title(solution_types{sol_type})
% end
% suptitle('Confusion matrices')
% orient landscape
% 
% print(fig, '-painters', '-dpdf', '-r600', '-bestfit',...
%     [figureFolder, 'figSX_confusion_matrices_classification_drugs',...
%     '_class', strrep(num2str(classthreshold),'.','_'),...
%     '_corr', strrep(num2str(corrthreshold),'.','_')])
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % plot all solutions in precision recall and specificity sensitivity
% % metrics
% marker_shape = {'o', '+', '*', 'x', 's', 'd', '^'};
% marker_color = {'b','k','r'};
% marker_color_num = [0 0 1;...
%                     0 0 0;...
%                     1 0 0];
% plotcrit = [3 2;
%             4 3;
%             1 1];
% 
% fig = figure('units','normalized','outerposition',[0 0 1 1]);
% 
% subplot(1,3,1)
% hold on
% for class_i = 1:length(classlabels)
%     gscatter(combined_report(class_i:3:end,plotcrit(1,1)),...
%              combined_report(class_i:3:end,plotcrit(1,2)),...
%         solution_types',...
%         repmat(marker_color{class_i},1,length(marker_shape)),...
%         strjoin(marker_shape,''))
% end
% xlim([0 1])
% ylim([0 1])
% axis square
% xlabel(report_labels{plotcrit(1,1)})
% ylabel(report_labels{plotcrit(1,2)})
% 
% subplot(1,3,2)
% hold on
% for class_i = 1:length(classlabels)
%     gscatter(1-combined_report(class_i:3:end,plotcrit(2,1)),...
%              combined_report(class_i:3:end,plotcrit(2,2)),...
%         solution_types',...
%         repmat(marker_color{class_i},1,length(marker_shape)),...
%         strjoin(marker_shape,''))
% end
% xlim([0 1])
% ylim([0 1])
% axis square
% xlabel(report_labels{plotcrit(2,1)})
% ylabel(report_labels{plotcrit(2,2)})
% 
% subplot(1,3,3)
% b = barh([combined_report(1:3:end,1)...
%      combined_report(2:3:end,1)...
%      combined_report(3:3:end,1)],...
%      'FaceColor','flat');
% for k = 1:length(marker_color)
%     b(k).CData = marker_color_num(k,:);
% end
% xlabel(report_labels{plotcrit(3,1)})
% set(gca, 'YTick', 1:length(solution_types))
% set(gca, 'YTickLabel', solution_types)
% axis square
% legend(C(1:3), 'Location', 'best')
% orient landscape
% 
% print(fig, '-painters', '-dpdf', '-r600', '-bestfit',...
%     [figureFolder, 'figSX_selection_criteria_classification_drugs',...
%     '_class', strrep(num2str(classthreshold),'.','_'),...
%     '_corr', strrep(num2str(corrthreshold),'.','_')])
% 
