function classification_report = create_classification_report(confmat, classlabels)

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
% convert class labels to string if necessary
if isnumeric(classlabels)
    classlabels = arrayfun(@(x) num2str(x), classlabels, 'unif', 0);
end
% convert classification report to a table
classification_report = array2table(classification_report,...
    'VariableNames', report_labels,...
    'RowNames', classlabels);