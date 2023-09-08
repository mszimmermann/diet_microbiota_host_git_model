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
picksol = ismember(met_bestsols{1}.selection_criterion,...
                                {'total PCC'});
model_classes = zeros(size(met_classes));
model_corr = zeros(size(met_classes));
for i=1:length(met_bestsols)
    bestsol = met_bestsols{i}.x(:, picksol);
    %normalize without the fisr coefficient (f)
    bestsol = bestsol(2:end)/max(abs(bestsol(2:end)));
    model_classes(i) = bestsol(end);
    model_corr(i) = met_bestsols{i}.selection_value(picksol);
end

% keep only posutions that pass the corr threshold
model_classes(model_corr<corrthreshold)=[];
met_classes(model_corr<corrthreshold)=[];

% set model classifiers based on threshold
model_classes(model_classes<=-classthreshold)=-1;
model_classes(model_classes>=classthreshold)=1;
model_classes((model_classes>-classthreshold) &...
              (model_classes<classthreshold))=0;



confmat = confusionmat(met_classes,model_classes);


    
