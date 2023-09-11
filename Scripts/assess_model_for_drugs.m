function assess_model_for_drugs(met_names, met_bestsols)
% assess the model fitting to drugs data 
% based on information of bacterial drug metabolism

% use classification criteria to 
% substrate (-1)
% product (1)
% neither (0)
% of the model coefficients

% corr threshold for reliable solutions
corrthreshold = 0.7;
classthreshold = 0.5;

solution_types = met_bestsols{1}.selection_criterion;
classification_reports = cell(size(solution_types))';

for sol_type = 1:length(solution_types)
    picksol = find(ismember(met_bestsols{1}.selection_criterion,...
                      solution_types(sol_type)));% {'total PCC'}));%
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
    model_classes = zeros(size(met_classes));
    model_corr = zeros(size(met_classes));
    for i=1:length(met_bestsols)
        if size(met_bestsols{i}.x,2)>=picksol
            bestsol = met_bestsols{i}.x(:, picksol);
            model_corr(i) = met_bestsols{i}.selection_value(picksol);
        else
            bestsol = zeros(size(met_bestsols{i}.x,1));
            model_corr(i) = 0;
        end
        %normalize without the first coefficient (f)
        bestsol = bestsol(2:end)/max(abs(bestsol(2:end)));
        model_classes(i) = bestsol(end);

    end

    % keep only posutions that pass the corr threshold
    model_classes(model_corr<corrthreshold)=[];
    met_classes(model_corr<corrthreshold)=[];

    % set model classifiers based on threshold
    model_classes(model_classes<=-classthreshold)=-1;
    model_classes(model_classes>=classthreshold)=1;
    model_classes((model_classes>-classthreshold) &...
                  (model_classes<classthreshold))=0;



    % get the confusion matrix
    [confmat, classlabels] = confusionmat(met_classes,model_classes);
    % create a matrix with 
    % accuracy, precision, recall, specificity, F1, support 
    report_labels = {'accuracy', 'precision', 'recall',...
                     'specificity', 'F1', 'support'}; 
    classification_report = zeros(length(classlabels),length(report_labels));
    for i = 1:length(classlabels)
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
    classification_reports{sol_type} = classification_report;
end

ST = reshape(repmat(solution_types, length(classlabels),1),[],1); 
C = repmat(arrayfun(@(x) num2str(x), classlabels, 'unif', 0),length(solution_types),1);
combined_report = cell2mat(vertcat(classification_reports));
figure
heatmap(combined_report, ...
    'YDisplayLabels', strcat(ST, ': class ', C),...
    'XDisplayLabels', report_labels)
caxis([0 1]) 

marker_shape = {'o', '+', '*', 'x', 's', 'd', '^'};
marker_color = {'b','k','r'};
marker_color_num = [0 0 1;...
                    0 0 0;...
                    1 0 0];
plotcrit = [3 2;
            4 3;
            1 1];

fig = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(1,3,1)
hold on
for class_i = 1:length(classlabels)
    gscatter(combined_report(class_i:3:end,plotcrit(1,1)),...
             combined_report(class_i:3:end,plotcrit(1,2)),...
        solution_types',...
        repmat(marker_color{class_i},1,length(marker_shape)),...
        strjoin(marker_shape,''))
end
xlim([0 1])
ylim([0 1])
axis square
xlabel(report_labels{plotcrit(1,1)})
ylabel(report_labels{plotcrit(1,2)})

subplot(1,3,2)
hold on
for class_i = 1:length(classlabels)
    gscatter(1-combined_report(class_i:3:end,plotcrit(2,1)),...
             combined_report(class_i:3:end,plotcrit(2,2)),...
        solution_types',...
        repmat(marker_color{class_i},1,length(marker_shape)),...
        strjoin(marker_shape,''))
end
xlim([0 1])
ylim([0 1])
axis square
xlabel(report_labels{plotcrit(2,1)})
ylabel(report_labels{plotcrit(2,2)})

subplot(1,3,3)
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