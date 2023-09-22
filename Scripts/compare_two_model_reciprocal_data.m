function compare_two_model_reciprocal_data(met_names1, met_bestsols1,...
                                   met_names2, met_bestsols2,...
                                   figureFolder)
% assess two modelling results (e.g. to two datasets,
% or two models to the same dataset)
% figureFolder is the folder to store plots
% plot original and restored data for both models

% use classification criteria to 
% substrate (-1)
% product (1)
% neither (0)
% of the model coefficients

% corr threshold for reliable solutions
corrthreshold = 0.7;
classthreshold = 0.5;%1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% intersect modelled metabolite names
[met_names, idx1, idx2] = intersect(lower(met_names1.CompoundName),...
                                    lower(met_names2.CompoundName));
% set classes for each model
% get classification by each model
model1_classes = zeros(length(met_names),1);
model2_classes = zeros(length(met_names),1);

if length(met_bestsols1.coefvalues)>5
    model1_classes = zeros(length(met_names),2);
end
if length(met_bestsols2.coefvalues)>5
    model2_classes = zeros(length(met_names),2);
end

model1_corr = met_bestsols1.x_sel_CorrRev(idx1);
for i=1:length(met_names)
    bestsol = met_bestsols1.x(idx1(i),:);
    %normalize without the first coefficient (f)
    bestsol = bestsol(2:end)/max(abs(bestsol(2:end)));
    if length(met_bestsols1.coefvalues)>5
        model1_classes(i,:) = bestsol(end-1:end);
    else
        model1_classes(i) = bestsol(end);
    end
end

model2_corr = met_bestsols2.x_sel_CorrRevLI(idx2);
%model2_corr = met_bestsols2.x_sel_CorrRev(idx2);
for i=1:length(met_names)
    bestsol = met_bestsols2.x(idx2(i),:);
    %normalize without the first coefficient (f)
    bestsol = bestsol(2:end)/max(abs(bestsol(2:end)));
    if length(met_bestsols2.coefvalues)>5
        model2_classes(i,:) = bestsol(end-1:end);
    else
        model2_classes(i) = bestsol(end);
    end
end

% keep only metabolites that pass the corr threshold
select_mets = find((model1_corr >= corrthreshold) & ...
                   (model2_corr >= corrthreshold)); 
% set model classifiers based on threshold
model1_classes(model1_classes<=-classthreshold)=-1;
model1_classes(model1_classes>=classthreshold)=1;
model1_classes((model1_classes>-classthreshold) &...
              (model1_classes<classthreshold))=0;

model2_classes(model2_classes<=-classthreshold)=-1;
model2_classes(model2_classes>=classthreshold)=1;
model2_classes((model2_classes>-classthreshold) &...
              (model2_classes<classthreshold))=0;

% add second colun in class variable if needed
if size(model1_classes,2)==1
    model1_classes(:,2) = model1_classes(:,1);
end
if size(model2_classes,2)==1
    model2_classes(:,2) = model2_classes(:,1);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot data and reciprocal data for all metabolites
fig = figure('units','normalized','outerposition',[0 0 1 1]);

spi = 2;
spj = 3;

idx1 = idx1(select_mets);
idx2 = idx2(select_mets);

for met_i = 1:length(select_mets)

    clf(fig);
    spidx = 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot original data for model1
    subplot(spi,spj,spidx)
    x = met_bestsols1.kmeanMatrix_joint_orig(idx1(met_i),:);
    x = reshape(x,4,6);
    plot(x')
    xlim([1,6])
    title(sprintf('Original data %s', met_bestsols1.modelname))
    lh=legend(cellfun(@(x) strrep(x, '_','-'),...
                met_bestsols1.kmeanMatrix_joint_names(1:6:end), 'unif', 0),...
                'Location', 'West');
    set(lh,'position',[.0 .75 .1 .1])
    spidx = spidx+1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot selected solution
    subplot(spi,spj,spidx)
    % take only metabolic coefs
    x = met_bestsols1.x(idx1(met_i),:);
    % remove first coefficient (f)
    x = x(2:end);
    % normalize by max in column
    x = x/max(abs(x));
    bar(x')
    title(sprintf('PCC=%.2f PCCLI=%.2f Class=%d/%d', ...
                    met_bestsols1.x_sel_CorrRev(idx1(met_i)),...
                    met_bestsols1.x_sel_CorrRevLI(idx1(met_i)),...
                    model1_classes(select_mets(met_i),1),...
                    model1_classes(select_mets(met_i),2)))
    ylim([-1.2 1.2])
    spidx = spidx + 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot restored data for selected solutions
    dataR = met_bestsols1.x_sel_dataR(idx1(met_i),:);
    dataR = reshape(dataR,6,4)';
    subplot(spi, spj, spidx)
    plot(dataR')
    xlim([1,6])
    title(sprintf('%s Rdata', met_bestsols1.selection_criterion))
    spidx = spidx+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % second model
    % plot original data for model2
    subplot(spi,spj,spidx)
    x = met_bestsols2.kmeanMatrix_joint_orig(idx2(met_i),:);
    x = reshape(x,4,6);
    plot(x')
    xlim([1,6])
    title(sprintf('Original data %s', met_bestsols2.modelname))
    lh=legend(cellfun(@(x) strrep(x, '_','-'),...
                met_bestsols2.kmeanMatrix_joint_names(1:4), 'unif', 0),...
                'Location', 'West');
    set(lh,'position',[.0 0.25 .1 .1])
    spidx = spidx+1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot selected solution
    subplot(spi,spj,spidx)
    % take only metabolic coefs
    x = met_bestsols2.x(idx2(met_i),:);
    % remove first coefficient (f)
    x = x(2:end);
    % normalize by max in column
    x = x/max(abs(x));
    bar(x')
    title(sprintf('PCC=%.2f PCCLI=%.2f, Class=%d/%d', ...
                    met_bestsols2.x_sel_CorrRev(idx2(met_i)),...
                    met_bestsols2.x_sel_CorrRevLI(idx2(met_i)),...
                    model2_classes(select_mets(met_i),1),...
                    model2_classes(select_mets(met_i),2)))
    ylim([-1.2 1.2])
    spidx = spidx + 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot restored data for selected solutions
    dataR = met_bestsols2.x_sel_dataR(idx2(met_i),:);
    dataR = reshape(dataR,4,6);
    subplot(spi, spj, spidx)
    plot(dataR')
    xlim([1,6])
    title(sprintf('%s Rdata', met_bestsols2.selection_criterion))
    
    
    suptitle(strrep(met_names(select_mets(met_i)),'_','-'));

    orient landscape

    print(fig, '-painters', '-dpsc', '-r600', '-append', '-bestfit',...
        [figureFolder, 'figSX_dataR_comparison_'...
        met_bestsols2.modelname, '_vs_',...
        met_bestsols1.modelname,...
        '_class', strrep(num2str(classthreshold),'.','_'),...
        '_corr', strrep(num2str(corrthreshold),'.','_')])

end
