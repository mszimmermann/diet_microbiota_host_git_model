
addpath(genpath('.\'));
% load data from Meier et al and reformat for the model
figureFolder = '.\Figures\';
outputFolder = '.\ProcessedData\';
inputFolder = '.\InputData\';

fileName = '2023_Meier_NatMet_42255_2023_802_MOESM3_ESM.xlsx';
 
 meier_data = readtable([inputFolder, fileName],...
                        'Sheet', 'intensities');
 %meier_data = readtable(fileName, 'Sheet', 'concentration');
 % select only cecum measurements
 select_content = ismember(meier_data.habitat, 'c'); 
 meier_matrix = meier_data{select_content, 5:end}';
 meier_mets = meier_data.Properties.VariableNames(5:end)';
 meier_condition = meier_data.Var1(select_content);
 meier_location = meier_data.site(select_content);
 meier_replicate = meier_data.replicate(select_content);
 
 meier_condition_unique = unique(meier_condition);
 meier_location_unique = unique(meier_location);
 select_locations = [3 6 9 12 14 15];
 
 % make a boxplot of all data
 fig = figure('units','normalized','outerposition',[0 0 1 1]);
 boxplot(log10(meier_matrix))
 set(gca, 'XTick', 1:length(meier_location))
 set(gca, 'XTickLabel', meier_location)
 
 % perform quantile normalization
 meier_matrix_norm =  meier_matrix;
 % replace 0 with nan
 meier_isnan = (meier_matrix_norm==0);
 meier_matrix_norm = quantilenorm(meier_matrix_norm);
 meier_matrix_norm(meier_isnan) = nan;
  
 fig = figure('units','normalized','outerposition',[0 0 1 1]);
 boxplot(log10(meier_matrix_norm))
 set(gca, 'XTick', 1:length(meier_location))
 set(gca, 'XTickLabel', meier_location)

% define colors for different nouse groups and tissue names
mycolors = [204 227 240;...%light blue
            0 115 178]/256;%dark blue
git_labels = {'Du', 'Je', 'Il', 'Cec', 'Col', 'Fec'};

meier_figure_file_name = 'Figures_final\Meier_data_plots_smooth_by_sum.ps';

% smooth first metabolite to create placeholder matrices
[smoothloc, smoothcond, smoothrep, smoothdata] = ...
        smooth_logitudinal_data(meier_location,...
        meier_condition,...
        meier_replicate,...
        meier_matrix_norm(1, :));
meier_data_smooth = zeros(length(meier_mets), length(smoothloc));   
 % plot all datapoints separately
printflag=0;
if printflag
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
end
for met_i=1:length(meier_mets)
     % smooth data for the current metabolite
     [smoothloc, smoothcond, smoothrep, smoothdata] = ...
        smooth_logitudinal_data(meier_location,...
        meier_condition,...
        meier_replicate,...
        meier_matrix_norm(met_i, :));
     smoothdata_normbymax = smoothdata./max(max(smoothdata));
     
     if printflag
         clf
         hold on
         for cond_i = 1:length(meier_condition_unique)
             subplot(2,2,1)
             % original data
             hold on
             curdata = meier_matrix(met_i, ...
                 ismember(meier_condition, meier_condition_unique{cond_i}));
             curloc = meier_location(ismember(meier_condition, meier_condition_unique{cond_i}));
             scatter(curloc,...
                 curdata, 'filled');
             title('Raw intensities')

             subplot(2,2,2)
             % quantile normalized
             hold on
             curdata = meier_matrix_norm(met_i, ...
                 ismember(meier_condition, meier_condition_unique{cond_i}));
             curloc = meier_location(ismember(meier_condition, meier_condition_unique{cond_i}));
             scatter(curloc,...
                 curdata, 'filled');
             title('Quantile-normalized')

             subplot(2,2,3)
             % selected sections, normalized
             hold on
             curdata = smoothdata(ismember(smoothcond, meier_condition_unique{cond_i}));
             curloc = smoothloc(ismember(smoothcond, meier_condition_unique{cond_i}));
             plot(curloc, curdata,...
                 'o', 'MarkerSize', 5, 'MarkerEdgeColor', mycolors(cond_i,:),...
                                     'MarkerFaceColor', mycolors(cond_i,:));

             kmeanMatrix = zeros(length(unique(curloc)),1);
             kstdMatrix = zeros(length(unique(curloc)),1);
             for tissue_i=1:length(unique(curloc))
                kmeanMatrix(tissue_i) = nanmean(curdata(curloc==tissue_i));
                kstdMatrix(tissue_i) = nanstd(curdata(curloc==tissue_i));
             end
             plot(1:length(unique(curloc)), kmeanMatrix,...
                  'LineWidth', 2, 'Color', mycolors(cond_i,:));
              errorbar(1:length(unique(curloc)), kmeanMatrix, kstdMatrix,...
                  '.','LineWidth', 2, 'Color', mycolors(cond_i,:)) 
              axis square
              set(gca, 'XTick', 1:length(git_labels))
              set(gca, 'XTickLabel', git_labels)
              xlim([1,length(git_labels)])  
              title('Smoothed')
             % normalize by max
             subplot(2,2,4)
             hold on
             % selected sections, normalized
             curdata = smoothdata_normbymax(ismember(smoothcond,...
                 meier_condition_unique{cond_i}));
             plot(curloc,curdata,...
                 'o', 'MarkerSize', 5, 'MarkerEdgeColor', mycolors(cond_i,:),...
                                     'MarkerFaceColor', mycolors(cond_i,:));

             kmeanMatrix = zeros(length(unique(curloc)),1);
             kstdMatrix = zeros(length(unique(curloc)),1);
             for tissue_i=1:length(unique(curloc))
                kmeanMatrix(tissue_i) = nanmean(curdata(curloc==tissue_i));
                kstdMatrix(tissue_i) = nanstd(curdata(curloc==tissue_i));
             end
             lg1(cond_i) = plot(1:length(unique(smoothloc)), kmeanMatrix,...
                  'LineWidth', 2, 'Color', mycolors(cond_i,:));
             errorbar(1:length(unique(smoothloc)), kmeanMatrix, kstdMatrix,...
                  '.','LineWidth', 2, 'Color', mycolors(cond_i,:)) 
            axis square
            set(gca, 'XTick', 1:length(git_labels))
            set(gca, 'XTickLabel', git_labels)
            xlim([1,length(git_labels)])
            title('Normalized to max')
         end
         suptitle(meier_mets{met_i});
         legend(lg1,meier_condition_unique, 'location', 'eastoutside');
         print(gcf, '-painters', '-dpsc2', '-r600', '-append', '-bestfit',...
             meier_figure_file_name);
     end
     % save smoothed data
     meier_data_smooth(met_i,:) = smoothdata;
 end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% apply modelling framework to Meier data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sampleType_unique = {'SPF', 'GF'};
sampleDiet_unique = {'Chow1', 'Chow2'};

% calculate mean across replicates
smoothloc_unique = unique(smoothloc);
meier_smooth_mean = zeros(size(meier_data_smooth,1),...
    length(smoothloc_unique)*length(meier_condition_unique));
meier_smooth_mean_condition = cell(length(smoothloc_unique)*length(meier_condition_unique),1);
idx=1;
for cond_i = 1:length(meier_condition_unique)
    for loc_i = 1:length(smoothloc_unique)
        meier_smooth_mean(:,idx) = nanmean(meier_data_smooth(:,...
            ismember(smoothcond, meier_condition_unique{cond_i}) &...
            (smoothloc == smoothloc_unique(loc_i))),2);
        meier_smooth_mean_condition{idx} = [meier_condition_unique{cond_i},...
            '_', git_labels{loc_i}];
    idx = idx+1;
    end
end

meanConditions = [strcat(meier_smooth_mean_condition, '_Chow1');
                  strcat(meier_smooth_mean_condition, '_Chow2')];              
% duplicate since we have only one diet
meanMatrix = [meier_smooth_mean meier_smooth_mean];
meanMatrix_mets = meier_mets;

mycolors = [0 115 178;... %dark blue
            204 227 240;...%light blue
            211 96 39;... %dark orange
            246 223 212]/256;%light orange
[bb,aa] = ndgrid(sampleType_unique, sampleDiet_unique); 
condLabels = strcat(aa(:),'-', bb(:));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the mean normalized value across replicates

fileNameFigure = [figureFolder...
    'fig_meierInt_diag_modelSMOOTH_2LIcoefHost_1LIbact.ps'];

diag_plot_flag = 0; % diagnostic plotting flag
if diag_plot_flag
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
end


% % % selected_mets are metabolites detected along the GI tract
selected_mets = 1:length(meier_mets)-1; % last metabolite is not measured

met_gitfits = cell(length(selected_mets),1);
met_bestsols = cell(length(selected_mets),1);

for met_i = 33:33%1:length(selected_mets)
   
        cmpd_interest_idx = selected_mets(met_i);
%         % set volume to CV or GF/WT
%         if contains(meanMatrix_mets{cmpd_interest_idx}, '_CV')
%             volumeMatrix = volumeMatrix_CVR;
%         else
%             volumeMatrix = volumeMatrix_GF;
%         end
        % calculate mean profiles            
        idx=1;
        kmeanMatrix_joint = zeros(4,6);
        kmeanMatrix_joint_names = cell(4,6);
        for diet_i = 1:length(sampleDiet_unique)
            for type_i = 1:length(sampleType_unique)
                selectDiet = sampleDiet_unique{diet_i};
                selectMouse = sampleType_unique{type_i};
                kmeanMatrix = meanMatrix(:,cellfun(@(x) contains(x,selectMouse) & ...
                                        contains(x,selectDiet),meanConditions));
                kmeanMatrix = kmeanMatrix(cmpd_interest_idx,1:6);
                kmeanMatrix_joint(idx,:) = kmeanMatrix;
                kmeanMatrix_joint_names(idx,:) = strcat(selectDiet,...
                                                        '_',...
                                                        selectMouse,...
                                                        '_',...
                                                        git_labels);
                idx = idx+1;
            end
        end
        
%         if use_volume_flag
%             % multiply by volume to get amounts
%             kmeanMatrix_joint = kmeanMatrix_joint.*volumeMatrix; 
%         end
        % normalize by max intensity
        kmeanMatrix_joint = kmeanMatrix_joint./max(max(kmeanMatrix_joint));
        kmeanMatrix_joint(isnan(kmeanMatrix_joint))=0;

        if nnz(kmeanMatrix_joint)
            %[gitfit] = fitGITmodel(kmeanMatrix_joint, ncond, shuffle_flag)
            [gitfit] = fitGITmodel(kmeanMatrix_joint, kmeanMatrix_joint_names, 2, 1);
            % get best solution
            [bestsol] = select_gitfit_sol(gitfit);

            % save gitfit and best solution for the current metabolite
            met_gitfits{met_i} = gitfit;
            met_bestsols{met_i} = bestsol;

            if diag_plot_flag
                gitfit_diag_plot(gitfit,meanMatrix_mets{met_i},fig)

                orient landscape
                %print to figure
                print(fig, '-painters', '-dpsc2', '-r600', '-append', '-bestfit',...
                        fileNameFigure);
                clf('reset')
            end    
        end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot correlation of restored and original data
% calculate differentce in corr distrbutions
filename = [figureFolder,...
    'Fig4a_histogram_MAXcorr_model_2LIcoefHost1LIcoefbact_meier'];
plot_gitfit_model_corr(met_gitfits, filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select LI within best total solution
x_selected = zeros(size(met_bestsols{1}.x,1), length(met_bestsols));
x_data_corr = zeros(1, length(met_bestsols));
x_data_corr_SI = zeros(1, length(met_bestsols));
x_data_corr_LI = zeros(1, length(met_bestsols));
x_data_corr_mean = zeros(1, length(met_bestsols));

for i=1:length(met_bestsols)
    picksol = ismember(met_bestsols{i}.selection_criterion,...
                                {'LI PCC within high total'});
    if size(met_bestsols{i}.x,2)==length(picksol)
        x_selected(:,i) = met_bestsols{i}.x(:,picksol);
        x_data_corr(i) = met_bestsols{i}.x_sel_CorrRev(picksol);
        x_data_corr_SI(i) = met_bestsols{i}.x_sel_CorrRevSI(picksol);
        x_data_corr_LI(i) = met_bestsols{i}.x_sel_CorrRevLI(picksol);
        x_data_corr_mean(i) = mean([x_data_corr(i) x_data_corr_SI(i) x_data_corr_LI(i)]);
    else
        % there was no total PCC above thresold to select LI PCC,
        % get the best possible total PCC
        picksol = ismember(met_bestsols{i}.selection_criterion,...
                                {'total PCC'});
        picksol = picksol(1:size(met_bestsols{i}.x,2));
        x_selected(:,i) = met_bestsols{i}.x(:,picksol);
        x_data_corr(i) = met_bestsols{i}.x_sel_CorrRev(picksol);
        x_data_corr_SI(i) = met_bestsols{i}.x_sel_CorrRevSI(picksol);
        x_data_corr_LI(i) = met_bestsols{i}.x_sel_CorrRevLI(picksol);
        x_data_corr_mean(i) = mean([x_data_corr(i) x_data_corr_SI(i) x_data_corr_LI(i)]);
    end        
end


% nrand=100;
% % call calculateAmatrix_final with zero matrices of correct size
% % to get coefvalues output, which is used to allocate variables further
% [~, coefvalues] = calculateAmatrix_final(zeros(length(condLabels),...
%                                 length(git_labels)));
% 
% % allocate variables for model prediction results
% x_met_mean = zeros(length(coefvalues),length(selected_mets));
% x_met_std = zeros(length(coefvalues),length(selected_mets));
% x_met_smooth = zeros(length(coefvalues),length(selected_mets));
% 
% % allocate variables for correlation calculated from values restored with
% % reverse problem with original and randomly shuffled parameters
% x_data_corr = zeros(length(selected_mets),1);
% x_resid = zeros(length(selected_mets),1);
% x_data_corr_shuffled = zeros(length(selected_mets),1);
% 
% % calculate reverse A matrix with zero matrix of the right size to allocate
% % variables according to Ra size
% [Ra] = calculateRAmatrix_final(zeros(size(coefvalues)));
% x_rdata = zeros(size(Ra,2),length(selected_mets));
% x_rdata_shuffled = zeros(size(Ra,2),length(selected_mets));
% 
% minval = 0.01;
% 
% fprintf('Starting parameter estimation for %d metabolites\n',length(selected_mets));
% 
% for met_i = 1:length(selected_mets)
%    
%         cmpd_interest_idx = selected_mets(met_i);
%         % calculate mean profiles            
%         idx=1;
%         kmeanMatrix_joint = zeros(4,6);
%         for diet_i = 1:length(sampleDiet_unique)
%             for type_i = 1:length(sampleType_unique)
%                 selectDiet = sampleDiet_unique{diet_i};
%                 selectMouse = sampleType_unique{type_i};
%                 kmeanMatrix = meanMatrix(:,cellfun(@(x) contains(x,selectMouse) & ...
%                                         contains(x,selectDiet),meanConditions));
%                 kmeanMatrix = kmeanMatrix(cmpd_interest_idx,1:6);
%                 kmeanMatrix_joint(idx,:) = kmeanMatrix;
%                 idx = idx+1;
%             end
%         end
%         % normalize by max intensity
%         kmeanMatrix_joint = kmeanMatrix_joint./max(max(kmeanMatrix_joint));
%         kmeanMatrix_joint(isnan(kmeanMatrix_joint))=0;
%         % replace small values with noise
%         kmeanMatrix_joint_orig = kmeanMatrix_joint;
%         kmeanMatrix_joint(kmeanMatrix_joint<minval)=...
%         kmeanMatrix_joint(kmeanMatrix_joint<minval)+rand(nnz(kmeanMatrix_joint<minval),1)*minval;
% 
%         % calculate fluxes for average profile
%         Aorig = calculateAmatrix_final(kmeanMatrix_joint_orig);
%         [A,coefvalues] = calculateAmatrix_final(kmeanMatrix_joint);
%         A(Aorig(:,1)==0,:)=[];
%         b = zeros(size(A,1),1);
%         options = optimoptions(@lsqlin,'Display', 'off');
% 
%         xlowerlim = zeros(size(A,2),1);
%         xlowerlim(2:end)=-Inf;
% 
%         %%%%%%%%%%%%%%%%%%%%%%%%%
%         %set bact to zero
%         xupperlim = Inf*ones(size(A,2),1);
%         %%%%%%%%%%%%%%%%%%%%%%%%%
%         % check if A is empty an if yes move to the next metabolite
%         if isempty(A)
%             continue;
%         end
% 
%         [x, xres] = lsqlin(A,b,[],[],[],[],xlowerlim, xupperlim, [], options);
%         x_met_smooth(:, met_i) = x;
%         x_resid(met_i) = xres;
% 
% 
%     [Ra,rb] = calculateRAmatrix_final(x*1000);
%     dataR = lsqlin(Ra,rb,[],[],[],[],zeros(1,size(Ra,2)), [], [], options);
%     % save reconstructed data
%     x_rdata(:,met_i) = dataR;
% 
%     %rsshape to matrix form
%     dataR = reshape(dataR,[],4)';
% 
%     x_data_corr(met_i) = corr(kmeanMatrix_joint_orig(:), dataR(:));
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % calculate reverse problem for shuffled x
%     x_shuffled = x(randperm(length(x)));
%     [Ra_shuffled,rb_shuffled] = calculateRAmatrix_final(x_shuffled*1000);
%     dataR_shuffled = lsqlin(Ra_shuffled,rb_shuffled,[],[],[],[],...
%                             zeros(1,size(Ra_shuffled,2)), [], [], options);
%     % save reconstructed data
%     x_rdata_shuffled(:,met_i) = dataR_shuffled;
% 
%     %rsshape to matrix form
%     dataR_shuffled = reshape(dataR_shuffled,[],4)';
% 
%     x_data_corr_shuffled(met_i) = corr(kmeanMatrix_joint_orig(:), dataR_shuffled(:));
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end
% 
% % get the original metabolomics data as vectors to calculate correlations
% % between original and restored measurements
% kmean_vector_joint_orig = zeros(size(Ra,2),length(selected_mets));
% kmean_vector_joint = zeros(size(Ra,2),length(selected_mets));
% for met_i = 1:length(selected_mets)
%    
%     cmpd_interest_idx = selected_mets(met_i);
%     % calculate mean profiles            
%     idx=1;
%     kmeanMatrix_joint = zeros(4,6);
%     for diet_i = 1:length(sampleDiet_unique)
%         for type_i = 1:length(sampleType_unique)
%             selectDiet = sampleDiet_unique{diet_i};
%             selectMouse = sampleType_unique{type_i};
%             kmeanMatrix = meanMatrix(:,cellfun(@(x) contains(x,selectMouse) & contains(x,selectDiet),meanConditions));
%             % get serum and liver data
%             kmeanMatrix = kmeanMatrix(cmpd_interest_idx,1:6);
%             kmeanMatrix_joint(idx,:) = kmeanMatrix;
%             idx = idx+1;
%         end
%     end
%     % normalize by max intensity
%     kmeanMatrix_joint = kmeanMatrix_joint./max(max(kmeanMatrix_joint));
%     % replace small values with noise
%     kmeanMatrix_joint_orig = kmeanMatrix_joint;
%     kmeanMatrix_joint(kmeanMatrix_joint<minval)=...
%     kmeanMatrix_joint(kmeanMatrix_joint<minval)+rand(nnz(kmeanMatrix_joint<minval),1)*minval;
%     
%     kmean_vector_joint_orig(:, met_i) = kmeanMatrix_joint_orig(:);
%     kmean_vector_joint(:, met_i) = kmeanMatrix_joint(:);
% end
% 
% % calculate Spearman correlations
% x_corr_Spearman = zeros(size(x_data_corr));
% corr_Mouse = zeros(length(x_data_corr),4);
% corr_Mouse_Spearman = zeros(length(x_data_corr),4);
% 
% for i = 1:length(x_data_corr)
%     dataR = reshape(x_rdata(:,i),[],4)';
%     dataOrig = reshape(kmean_vector_joint_orig(:,i),4,[]);
% 
%     x_corr_Spearman(i) = corr(dataOrig(:), dataR(:), 'type', 'Spearman');
%     corr_Mouse(i,:) = diag(corr(dataOrig', dataR'));
%     corr_Mouse_Spearman(i,:) = diag(corr(dataOrig', dataR', 'type', 'Spearman'));
% 
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot data and restored intensities and model coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define colora and GIT section names for plotting
mycolors = [0 115 178;... %dark blue
            204 227 240;...%light blue
            211 96 39;... %dark orange
            246 223 212]/256;%light orange
git_labels = {'Du', 'Je', 'Il', 'Cec', 'Col', 'Fec'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting file name
fileNameprofiles = 'fig_profiles_figureselected_modelSMOOTH_2LIcoefHost_1LIbact_MeierDataNORMsmoothbysum.ps';
          
curData_cols = reshape(meanConditions, 6, 4);
compoundsInterest = 1:length(meanMatrix_mets);

fig = figure('units','normalized','outerposition',[0 0 1 1]);

fprintf('Plotting profiles for %d metabolites\n',length(compoundsInterest));
for cpdix=1:length(compoundsInterest)
    
    testidx = compoundsInterest(cpdix);

    spx=1;
    spy=3;
    spidx = 1;
    coloridx = 1;
    idx=1;
    curmat = zeros(4,6);
    cur_data = reshape(kmean_vector_joint_orig(:,testidx)', 4, 6)';
    %cur_data = reshape(meanMatrix(testidx,:), 6, 4);
    cur_rdata = reshape(x_rdata(:,testidx)', 6, 4);
    
    %normalize rdata to max
    cur_rdata = (cur_rdata-0.8);
    cur_rdata = cur_rdata/max(max(cur_rdata));
    
    legend_entries = cell(4,1);
    for i = 1:size(cur_data,2)

        subplot(spx,spy,idx)
        hold on
        h(i) = plot(cur_data(:,i),...
             'LineWidth', 2,...
             'Color', mycolors(i,:));
        ylabel('Original normbymax')

        subplot(spx,spy,idx+1);
        hold on
        plot((cur_rdata(:,i)),...
                 'LineWidth', 2,...
                 'Color', mycolors(i,:));
        ylabel('Restored normbymax')

    end
    title(sprintf('%s PCC=%.2f',...
                  meanMatrix_mets{testidx},...
                        x_data_corr(testidx)),...
                        'Interpreter', 'none');

    legend(h, curData_cols(1,:),...
             'Interpreter', 'none');%, 'Location', 'bestoutside');
    
    subplot(spx,spy,idx+2);
    
    curcoefs = x_met_smooth(2:end, testidx);
    barh(curcoefs./max(abs(curcoefs)))
    set(gca, 'YTick', 1:length(curcoefs));
    set(gca, 'YTickLabel', coefvalues(2:end));
    set(gca, 'YDir','reverse')
    ylim([0.5 length(curcoefs)+0.5])
    xlim([-1 1]);
    axis square
    
    for spi = 1:(spx*spy)-1
        subplot(spx,spy,spi)
        set(gca, 'XTick', 1:6)
        xlim([1 6])
        %ylim([0 1])
        set(gca, 'XTick', 1:length(git_labels))
        set(gca, 'XTickLabel', git_labels)
        
        axis square
    end
%     spt = suptitle({sprintf('MZ=%.3f',testmz(idx,1)),...
%                                         testannID{1},...
%                                         testann{1}});
%     set(spt,'FontSize',8,'FontWeight','normal')
    orient landscape
    %print to figure
    print(gcf, '-painters', '-dpsc2', '-r600', '-append', '-bestfit',...
            fileNameprofiles);%[figureFolder,...
            % fileNameprofiles])
    clf('reset')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cluster only metabolites with corr>0.7
% uncomment if using whole metabolite table
clustidx = 1:size(x_selected,2); % for testing purposed plot all 
clustdata = x_selected(2:end,clustidx);
clustrows = met_bestsols{1}.coefvalues(2:end);
% filter by resiprocal corr
clustcorr = x_data_corr(clustidx);
clustdata = clustdata(:,clustcorr>=0.7);
clustidx = clustidx(clustcorr>=0.7);

for i=1:size(clustdata,2)
    clustdata(:,i) = clustdata(:,i)/max(abs(clustdata(:,i)));
end
clustnan = isnan(sum(clustdata,1));
clustdata(:,clustnan) = [];
clustidx(clustnan)=[];

clustdist ='cityblock';%'correlation';%'euclidean';%
cgo = clustergram(clustdata,...
            'RowLabels', clustrows,...
            'ColumnLabels', clustidx,...
            'ColumnPdist',clustdist,...
            'RowPdist', clustdist,...
            'DisplayRange', 1,...
            'colormap', redbluecmap);
% at this point colorbar and print to figure have to be done manually in
% the clustergram window
% print to figure
fig = cgo.plot;

orient landscape
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
    [figureFolder,...
    'fig_clustergram_2LIhos1LIbact_model_coefs_MeierData_Rmodel_LIwithinTotcorr_0_7_smoothbysum'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print best solutions to files
filename = [outputFolder ...
            'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_Meier'];
% create met_info object needed for the printing function
met_info.CompoundID = meier_mets;
met_info.CompoundName = meier_mets;          
% print solutions to files
print_bestsol_to_files(met_info, met_bestsols, filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % save model results to file - raw coefficients
% 
% fid = fopen([outputFolder ...
%             'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_MeierData_smoothbysum.csv'], 'w');
%              %'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_allions.csv'], 'w');
% fprintf(fid, 'CompoundName\tReciprocalCorr');
% for i=1:length(coefvalues)
%     fprintf(fid, '\t%s', coefvalues{i});
% end
% fprintf(fid, '\n');
% for i=1:size(x_met_smooth,2)
%     fprintf(fid, '%s', meier_mets{i});
%     fprintf(fid, '\t%.3f', x_data_corr(i));
%     for j=1:size(x_met_smooth,1)
%         fprintf(fid, '\t%e', x_met_smooth(j,i));
%     end
%     fprintf(fid, '\n');
% end
% fclose(fid);
% 
% % save model results to file - normalized by max coefficient
% fid = fopen([outputFolder...
%     'model_results_SMOOTH_normbyabsmax_2LIcoefHost1LIcoefbact_MeierData_smoothbysum.csv'], 'w');
%     %'model_results_SMOOTH_normbyabsmax_2LIcoefHost1LIcoefbact_allions.csv'], 'w');
% fprintf(fid, 'CompoundName\tReciprocalCorr');
% for i=1:length(coefvalues)
%     fprintf(fid, '\t%s', coefvalues{i});
% end
% fprintf(fid, '\n');
% for i=1:size(x_met_smooth,2)
%     fprintf(fid, '%s', meier_mets{i});
%     fprintf(fid, '\t%.3f', x_data_corr(i));
%     for j=1:size(x_met_smooth,1)
%         fprintf(fid, '\t%.3f', x_met_smooth(j,i)./max(abs(x_met_smooth(:,i))));
%     end
%     fprintf(fid, '\n');
% end
% fclose(fid);
% 
% % save model results to file - metabolism_coefficients_only
% fid = fopen([outputFolder...
%     'model_results_SMOOTH_normbyabsmax_ONLYMETCOEF_2LIcoefHost1LIcoefbact_MeierData_smoothbysum.csv'], 'w');
%     %'model_results_SMOOTH_normbyabsmax_ONLYMETCOEF_2LIcoefHost1LIcoefbact_allions.csv'], 'w');
% fprintf(fid, 'CompoundName\tReciprocalCorr');
% for i=2:length(coefvalues)
%     fprintf(fid, '\t%s', coefvalues{i});
% end
% fprintf(fid, '\n');
% for i=1:size(x_met_smooth,2)
%     fprintf(fid, '%s', meier_mets{i});
%     fprintf(fid, '\t%.3f', x_data_corr(i));
%     for j=2:size(x_met_smooth,1)
%         fprintf(fid, '\t%.3f', x_met_smooth(j,i)./max(abs(x_met_smooth(2:end,i))));
%     end
%     fprintf(fid, '\n');
% end
% fclose(fid);
%     
% % save model results to file - reciprocal data restoration
% fid = fopen([outputFolder...
%     'model_results_SMOOTH_normbyabsmax_reciprocal_problem_MeierData_smoothbysum.csv'], 'w');
%     %'model_results_SMOOTH_normbyabsmax_reciprocal_problem_allions.csv'], 'w');
% columnNames = cell(size(x_rdata,1),1);
% idx=1;
% for diet_i = 1:length(sampleDiet_unique)
%     for type_i = 1:length(sampleType_unique)
%         for tiss_i = 1:length(git_labels)
%                 selectDiet = sampleDiet_unique{diet_i};
%                 selectMouse = sampleType_unique{type_i};
%                 selectTissue = git_labels{tiss_i};
%                 columnNames{idx} = [selectDiet '_' selectMouse '_' selectTissue];
%                 idx = idx+1;
%         end
%     end
% end
%          
% fprintf(fid, 'CompoundName\tReciprocalCorr\tRandomCorr');
% for i=1:length(columnNames)
%     fprintf(fid, '\t%s', columnNames{i});
% end
% for i=1:length(columnNames)
%     fprintf(fid, '\tRecip_%s', columnNames{i});
% end
% for i=1:length(columnNames)
%     fprintf(fid, '\tRandom_%s', columnNames{i});
% end
% fprintf(fid, '\n');
% for i=1:size(x_rdata,2)
%     fprintf(fid, '%s', meier_mets{i});
%     fprintf(fid, '\t%.3f', x_data_corr(i));
%     fprintf(fid, '\t%.3f', x_data_corr_shuffled(i));
%     for j=1:size(kmean_vector_joint_orig(:,i),1)
%         fprintf(fid, '\t%.3f', kmean_vector_joint_orig(j,i));
%     end
%     for j=1:size(x_rdata(:,i),1)
%         fprintf(fid, '\t%.3f', x_rdata(j,i));
%     end
%     for j=1:size(x_rdata_shuffled(:,i),1)
%         fprintf(fid, '\t%.3f', x_rdata_shuffled(j,i));
%     end
%     fprintf(fid, '\n');
% end
% fclose(fid);

