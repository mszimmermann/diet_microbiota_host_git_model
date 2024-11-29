function [model_x_normalized] = normalize_model_coefs(met_bestsols)
% normalize modelling results and assess the
% confidence of bacterial products/substrates (confidentclassflag)

% use classification criteria to 
% substrate (-1)
% product (1)
% neither (0)
% of the model coefficients

% corr threshold for reliable solutions
corrthreshold = 0.7;
corrthresholdLI = 0.7;
classthreshold = 0.5;%1;
% make a flag for positive or negative class additional requirement of 
% LI PCC passing the threshold as well
%flag_strictclass = 1;%0;%


% normalize mode coefficients
model_corr = met_bestsols.x_sel_CorrRev;
model_corrLI = met_bestsols.x_sel_CorrRevLI;
model_x_normalized = zeros(size(met_bestsols.x,1), size(met_bestsols.x,2)-1);
confidentclassflag = ones(size(met_bestsols.x,1),1);

for i=1:size(met_bestsols.x,1)
    bestsol = met_bestsols.x(i,:);
    %normalize without the first coefficient (f)
    bestsol = bestsol(2:end)/max(abs(bestsol(2:end)));
    model_x_normalized(i,:) = bestsol;

    % check if metabolite is a bacterial substrate or product
    if length(met_bestsols.coefvalues)>5
        model_classes = bestsol(end-1:end);
    else
        model_classes = bestsol(end);
    end
    
    % set model classifiers based on threshold
    model_classes(model_classes<=-classthreshold)=-1;
    model_classes(model_classes>=classthreshold)=1;
    model_classes((model_classes>-classthreshold) &...
              (model_classes<classthreshold))=0;

    confidentclassflag(i) = (sum(abs(model_classes))==0) | ...
                            ((sum(abs(model_classes))>0) &...
                             (model_corrLI(i)>corrthresholdLI));
end

% convert matrix to table
model_x_normalized = array2table(model_x_normalized,...
            'VariableNames',met_bestsols.coefvalues(2:end));
model_x_normalized.confidentclassflag = confidentclassflag;
model_x_normalized.x_sel_CorrRev = model_corr;
model_x_normalized.x_sel_CorrRevLI = model_corrLI;
