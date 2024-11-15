%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare modelling results to differential analysis to aid interpretation
% of the model coefficients

% call script defining file dependenciesand global variables
addpath(genpath('.\'));
add_global_and_file_dependencies

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Required files:
% 'metabolites_allions_combined_norm_intensity.csv'
% 'metabolites_allions_combined_formulas_with_metabolite_filters.csv'
% 'table_diff_abundance_metabolite_ions_removed2outliers.csv'
% 'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_bestDC.csv'
%
% Output: 
% Files:
% Figures:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
metaboliteFilters = readtable([inputFolder,...
    'metabolites_allions_combined_formulas_with_metabolite_filters.csv']);

metaboliteDiff = readtable([resultsFolder,...
    'table_diff_abundance_metabolite_ions_removed2outliers.csv']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = [outputFolder ...
            'model_results_SMOOTH_raw_2LIcoefHost1LIcoefbact_bestDC'];
corrthreshold = 0.7;
% try combining solutions
sel_crit1 = 'IP';
sel_crit2 = 'LI_PCC_within_high_total';

[met_info_combined, met_bestsol_combined] = ...
                combine_bestsols_from_file(filename, ...
                sel_crit1, sel_crit2,...
                corrthreshold);

% analyze only annotated metabolites for which solution was possible
selected_mets = find( (metaboliteFilters.MetaboliteFilter==1) &...
                      ~isnan(met_bestsol_combined.x_sel_CorrRev) &...
                      (sum(abs(met_bestsol_combined.x),2)>0));

met_bestsol_combined = filter_bestsols_by_index(met_bestsol_combined, selected_mets);

% filter differential analysis table with the same indeces
metaboliteDiff = metaboliteDiff(selected_mets,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get information from differential analysis on whether substrates and
% products are changing in at least one low intestine tissue in at least
% one diet between GF and colonized mice
liTissues = {'Cecum', 'Colon', 'Feces'};
dietChangeNames = {'HFD_DCGF', 'CTR_DCGF'};
fcThreshold = log2(1.5);
pfdrThreshold = 0.05;
pvalColumn = '_pFDR_'; % choose between _p_ and _pFDR_

% normalize coefficient matrix
model_x_norm = met_bestsol_combined.x(:,2:end);
model_x_norm = model_x_norm./max(abs(model_x_norm), [], 2);

res_norm_linHCD = zeros(length(dietChangeNames), length(liTissues));
res_norm_tanhHCD = zeros(length(dietChangeNames), length(liTissues));
res_norm_linHFD = zeros(length(dietChangeNames), length(liTissues));
res_norm_tanhHFD = zeros(length(dietChangeNames), length(liTissues));
corrS_HCD = zeros(length(dietChangeNames), length(liTissues));
corrS_HFD = zeros(length(dietChangeNames), length(liTissues));

mycolors = [211 96 39;...   %dark orange 
            0 115 178]/256; %dark blue

% make one figure for HCD coef and one for HFD
figHCD = figure('units','normalized','outerposition',[0 0 1 1]);
figHFD = figure('units','normalized','outerposition',[0 0 1 1]);

set(figHCD, 'PaperUnits', 'centimeters');
set(figHCD, 'PaperPosition', [0 0 10 25]); %x_width=10cm y_width=15cm

set(figHFD, 'PaperUnits', 'centimeters');
set(figHFD, 'PaperPosition', [0 0 10 25]); %x_width=10cm y_width=15cm
 
spidx = 1;
for diet_i=1:length(dietChangeNames)
        for tissue_i=1:length(liTissues)
  
            curfc = log2(metaboliteDiff{:,...
                 ismember(metaboliteDiff.Properties.VariableNames,...
                 {[liTissues{tissue_i} '_fc_', dietChangeNames{diet_i}]})});

            % calculate correlation between fold changes and coefficients
            [corrcoef, corrp] = corr(model_x_norm(:,end-1), curfc);
            [corrScoef, corrSp] = corr(model_x_norm(:,end-1), curfc, 'type', 'spearman');
            
            axHCD = subplot(2,3,spidx, 'Parent', figHCD);
            scatter(axHCD, curfc, model_x_norm(:,end-1), 'MarkerEdgeColor',mycolors(diet_i,:))
                
            xlabel(axHCD, met_bestsol_combined.coefvalues(end-1))
            yl = ylabel(axHCD, [liTissues{tissue_i} '_fc_', dietChangeNames{diet_i}]);
            set(yl, 'Interpreter', 'none')
            title(axHCD, {sprintf("PCC=%.2f, p=%.2e", corrcoef, corrp),...
                   sprintf("SCC=%.2f, p=%.2e", corrScoef, corrSp)})

            % add linear fit
            [P, polyS] = polyfit(curfc,model_x_norm(:,end-1),1);
            yfit = polyval(P,curfc);
            hold(axHCD, 'on');
            plot(axHCD, curfc,yfit,'k--', 'LineWidth', 3);
            eqn = string("y=" + P(1)) + "x+" + string(P(2));
            text(axHCD, min(curfc),max(model_x_norm(:,end-1)),eqn,...
                "HorizontalAlignment","left","VerticalAlignment","top")

            % fit sigmoid tannh function
            [a,res] = lsqcurvefit(@myfunsigmoid,[1 1],curfc,model_x_norm(:,end-1));
            plot(axHCD, sort(curfc),myfunsigmoid(a,sort(curfc)), 'k:', 'LineWidth', 3)
            text(axHCD, min(curfc),min(model_x_norm(:,end)+0.3),...
                sprintf('y=tanh(%.3f x)(x>=0) + tanh(%.3f x)(x<0)', a(1), a(2)),...
                "HorizontalAlignment","left","VerticalAlignment","top")

            xlim(axHCD, [-20 20])
            ylim(axHCD, [-1 1])
            
            res_norm_linHCD(diet_i, tissue_i) = polyS.normr;
            res_norm_tanhHCD(diet_i, tissue_i) = sqrt(res);
            corrS_HCD(diet_i, tissue_i) = corrScoef;
            
            % axis square

            axHFD = subplot(2,3,spidx, 'Parent', figHFD);
            % calculate correlations between fold changes and model
            % coefficients
            [corrcoef, corrp] = corr(model_x_norm(:,end), curfc);
            [corrScoef, corrSp] = corr(model_x_norm(:,end), curfc, 'type', 'spearman');
           
          
            scatter(axHFD, curfc, model_x_norm(:,end), 'MarkerEdgeColor',mycolors(diet_i,:))
            xlabel(axHFD,met_bestsol_combined.coefvalues(end))
            yl = ylabel(axHFD,[liTissues{tissue_i} '_fc_', dietChangeNames{diet_i}]);
            set(yl, 'Interpreter', 'none')
            xlim(axHFD, [-20 20])
            ylim(axHFD, [-1 1])
            %axis square

            % add linear fit
            [P, polyS] = polyfit(curfc,model_x_norm(:,end),1);
            yfit = polyval(P,curfc);
            hold(axHFD, 'on');
            plot(axHFD, curfc,yfit,'k--', 'LineWidth', 3);
            eqn = string("y=" + P(1)) + "x+" + string(P(2));
            text(axHFD, min(curfc),max(model_x_norm(:,end)),eqn,...
                "HorizontalAlignment","left","VerticalAlignment","top")

            % fit tannh function
            [a,res] = lsqcurvefit(@myfunsigmoid,[1 1],curfc,model_x_norm(:,end));
            plot(axHFD, sort(curfc),myfunsigmoid(a,sort(curfc)),'k:', 'LineWidth', 3)
            text(axHFD, min(curfc),min(model_x_norm(:,end))+0.3,...
                sprintf('y=tanh(%.3f x)(x>=0) + tanh(%.3f x)(x<0)', a(1), a(2)),...
                "HorizontalAlignment","left","VerticalAlignment","top")

            title(axHFD,{sprintf("PCC=%.2f, p=%.2e", corrcoef, corrp),...
                   sprintf("SCC=%.2f, p=%.2e", corrScoef, corrSp)});
                   %sprintf('res=%.3f (l) %.3f (s)', polyS.normr, sqrt(res))})
           
            res_norm_linHFD(diet_i, tissue_i) = polyS.normr;
            res_norm_tanhHFD(diet_i, tissue_i) = sqrt(res);
            corrS_HFD(diet_i, tissue_i) = corrScoef;
            
            spidx=spidx+1;
            res_idx = res_idx+1;

        end
end
sgtitle(sprintf('Correlation between HCD model coefficients and metabolite fold changes, n=%d',...
            size(model_x_norm,1)), 'Parent', figHCD);
sgtitle(sprintf('Correlation between HFD model coefficients and metabolite fold changes, n=%d',...
            size(model_x_norm,1)), 'Parent', figHFD);

orient(figHCD, 'landscape')
print(figHCD, '-vector', '-dpdf', '-r600', ...
    [figureFolder,...
    'figED_scatter_modelHCD_vs_met_foldchanges'])

orient(figHFD, 'landscape')
print(figHFD, '-vector', '-dpdf', '-r600', ...
    [figureFolder,...
    'figED_scatter_modelHFD_vs_met_foldchanges'])


% plot spearman corr
fig = figure('units','normalized','outerposition',[0 0 1 1]);
set(fig, 'PaperUnits', 'centimeters');
set(fig, 'PaperPosition', [0 0 10 15]); %x_width=10cm y_width=15cm

subplot(2,3,1)
    h = heatmap(corrS_HCD);
    ax = gca;
    ax.XData = liTissues;
    ax.YData = dietChangeNames;
    clim([0 0.7]);
    title('HCD model Spearman coefficients')
    h.CellLabelFormat = '%.2f';

subplot(2,3,2)
    h = heatmap(res_norm_linHCD); 
    ax = gca;
    ax.XData = liTissues;
    ax.YData = dietChangeNames;
    clim([25 29]);
    title('HCD residual norm linear')
    h.CellLabelFormat = '%.2f';

subplot(2,3,3)
    h = heatmap(res_norm_tanhHCD); 
    ax = gca;
    ax.XData = liTissues;
    ax.YData = dietChangeNames;
    clim([25 29]);
    title('HCD residual norm tanh')
    h.CellLabelFormat = '%.2f';

subplot(2,3,4)
    h = heatmap(corrS_HFD);
    ax = gca;
    ax.XData = liTissues;
    ax.YData = dietChangeNames;
    clim([0 0.7]);
    title('HFD model Spearman coefficients')
    h.CellLabelFormat = '%.2f';
    
subplot(2,3,5)
    h = heatmap(res_norm_linHFD);
    ax = gca;
    ax.XData = liTissues;
    ax.YData = dietChangeNames;
    clim([22 27]);
    title('HFD residual norm linear')
    h.CellLabelFormat = '%.2f';

subplot(2,3,6)
    h = heatmap(res_norm_tanhHFD);
    ax = gca;
    ax.XData = liTissues;
    ax.YData = dietChangeNames;
    clim([22 27]);
    title('HFD residual norm tanh')
    h.CellLabelFormat = '%.2f';    

orient(fig, 'landscape')
print(fig, '-vector', '-dpdf', '-r600', ...
    [figureFolder,...
    'figED_heatmap_model_vs_met_corr_and_residuals'])




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kegg_sub_li_changes = zeros(size(shortestPathTable,1), length(liTissues)*length(dietChangeNames));
kegg_prod_li_changes = zeros(size(shortestPathTable,1), length(liTissues)*length(dietChangeNames));
% record also liver and serum changes
kegg_sub_liver_changes = zeros(size(shortestPathTable,1), length(dietChangeNames));
kegg_prod_liver_changes = zeros(size(shortestPathTable,1), length(dietChangeNames));
kegg_sub_serum_changes = zeros(size(shortestPathTable,1), length(dietChangeNames));
kegg_prod_serum_changes = zeros(size(shortestPathTable,1), length(dietChangeNames));
for i=1:size(shortestPathTable,1)
    idx=1;
    % sub
    cur_sub_idx = kegg_ids_unique_index(ismember(kegg_ids_unique,...
                                        shortestPathTable.Substrate{i}));
    % prod
    cur_prod_idx = kegg_ids_unique_index(ismember(kegg_ids_unique,...
                                        shortestPathTable.Product{i}));
    for diet_i=1:length(dietChangeNames)
        for tissue_i=1:length(liTissues)
            % sub
            kegg_sub_li_changes(i,idx) = (abs(log2(metaboliteDiff{cur_sub_idx,...
                 ismember(metaboliteDiff.Properties.VariableNames,...
                 {[liTissues{tissue_i} '_fc_', dietChangeNames{diet_i}]})})) >= fcThreshold).*...
                 (metaboliteDiff{cur_sub_idx,...
                 ismember(metaboliteDiff.Properties.VariableNames,...
                 {[liTissues{tissue_i} pvalColumn, dietChangeNames{diet_i}]})} <= pfdrThreshold);
      
            % prod
            kegg_prod_li_changes(i,idx) = (abs(log2(metaboliteDiff{cur_prod_idx,...
                 ismember(metaboliteDiff.Properties.VariableNames,...
                 {[liTissues{tissue_i} '_fc_', dietChangeNames{diet_i}]})})) >= fcThreshold).*...
                 (metaboliteDiff{cur_prod_idx,...
                 ismember(metaboliteDiff.Properties.VariableNames,...
                 {[liTissues{tissue_i} pvalColumn, dietChangeNames{diet_i}]})} <= pfdrThreshold);
            idx = idx+1;
        end
        % liver
        % sub (record sign of change)
        kegg_sub_liver_changes(i,diet_i) = (abs(log2(metaboliteDiff{cur_sub_idx,...
             ismember(metaboliteDiff.Properties.VariableNames,...
             {['Liver_fc_', dietChangeNames{diet_i}]})})) >= fcThreshold).*...
             (metaboliteDiff{cur_sub_idx,...
             ismember(metaboliteDiff.Properties.VariableNames,...
             {['Liver' pvalColumn, dietChangeNames{diet_i}]})} <= pfdrThreshold).*...
             sign(log2(metaboliteDiff{cur_sub_idx,...
             ismember(metaboliteDiff.Properties.VariableNames,...
             {['Liver_fc_', dietChangeNames{diet_i}]})}));

        % prod (record sign of change)
        kegg_prod_li_changes(i,diet_i) = (abs(log2(metaboliteDiff{cur_prod_idx,...
             ismember(metaboliteDiff.Properties.VariableNames,...
             {['Liver_fc_', dietChangeNames{diet_i}]})})) >= fcThreshold).*...
             (metaboliteDiff{cur_prod_idx,...
             ismember(metaboliteDiff.Properties.VariableNames,...
             {['Liver' pvalColumn, dietChangeNames{diet_i}]})} <= pfdrThreshold).*...
             sign(log2(metaboliteDiff{cur_prod_idx,...
             ismember(metaboliteDiff.Properties.VariableNames,...
             {['Liver_fc_', dietChangeNames{diet_i}]})}));
        % serum
        % sub (record sign of change)
        kegg_sub_serum_changes(i,diet_i) = (abs(log2(metaboliteDiff{cur_sub_idx,...
             ismember(metaboliteDiff.Properties.VariableNames,...
             {['Serum_fc_', dietChangeNames{diet_i}]})})) >= fcThreshold).*...
             (metaboliteDiff{cur_sub_idx,...
             ismember(metaboliteDiff.Properties.VariableNames,...
             {['Serum' pvalColumn, dietChangeNames{diet_i}]})} <= pfdrThreshold).*...
             sign(log2(metaboliteDiff{cur_sub_idx,...
             ismember(metaboliteDiff.Properties.VariableNames,...
             {['Serum_fc_', dietChangeNames{diet_i}]})}));

        % prod (record sign of change)
        kegg_prod_serum_changes(i,diet_i) = (abs(log2(metaboliteDiff{cur_prod_idx,...
             ismember(metaboliteDiff.Properties.VariableNames,...
             {['Serum_fc_', dietChangeNames{diet_i}]})})) >= fcThreshold).*...
             (metaboliteDiff{cur_prod_idx,...
             ismember(metaboliteDiff.Properties.VariableNames,...
             {['Serum' pvalColumn, dietChangeNames{diet_i}]})} <= pfdrThreshold).*...
             sign(log2(metaboliteDiff{cur_prod_idx,...
             ismember(metaboliteDiff.Properties.VariableNames,...
             {['Serum_fc_', dietChangeNames{diet_i}]})}));
             
    end
    
    if mod(i,1000)==0
        fprintf('Extracted fold changes for %d metabolites\n',i);
    end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate correlations between RNA and metabolites
kegg_sub_prod_EC_sub_corr = zeros(size(kegg_substrate_product_1enzyme,1), length(species_list));
kegg_sub_prod_EC_sub_corrP = ones(size(kegg_substrate_product_1enzyme,1), length(species_list));
kegg_sub_prod_EC_sub_sum_corr = zeros(size(kegg_substrate_product_1enzyme,1), length(species_list));
kegg_sub_prod_EC_sub_sum_corrP = ones(size(kegg_substrate_product_1enzyme,1), length(species_list));
kegg_sub_prod_EC_sub_corr_id = cell(size(kegg_substrate_product_1enzyme,1), length(species_list));
kegg_sub_prod_EC_sub_corr_geneidx = cell(size(kegg_substrate_product_1enzyme,1), length(species_list));
% prepare matrices for linear model fits and r adjusted
kegg_sub_prod_EC_sub_mdlX = zeros(size(kegg_substrate_product_1enzyme,1), length(species_list));
kegg_sub_prod_EC_sub_mdlP = ones(size(kegg_substrate_product_1enzyme,1), length(species_list));
kegg_sub_prod_EC_sub_mdlRsqadj = zeros(size(kegg_substrate_product_1enzyme,1), length(species_list));

% substrates
for i=1:size(kegg_substrate_product_1enzyme,1)
    curenzymes = strsplit(kegg_substrate_product_1enzyme.EC_path{i},'|');
    if ~isempty(curenzymes)
        % intersect enzymes with the ones from the gene table
        [~, idxEC, idxGene] = intersect(curenzymes, geneTable_EC_unique, 'stable');
        if ~isempty(idxEC)
            curenzymes_geneidx = geneTable_EC_species(idxGene,:);
            % prepare tables for correlation values etc
            curenzymes_mat = zeros(length(idxEC), length(species_list));
            curenzymes_p = ones(length(idxEC), length(species_list));
            curenzymes_ec = cell(length(idxEC), length(species_list));
            % for linear model
            curenzymes_mdlX = zeros(length(idxEC), length(species_list));
            curenzymes_mdlXp = ones(length(idxEC), length(species_list));
            curenzymes_mdlRsqadj = zeros(length(idxEC), length(species_list));
            
            for k=1:size(curenzymes_geneidx,2)
                for j=1:size(curenzymes_geneidx,1)
                    if ~isempty(curenzymes_geneidx{j,k})
                        cur_expression = table2array(countsMatrixGetMM_RNA_table(curenzymes_geneidx{j,k},gene_idx))';
                        cur_expression = sum(cur_expression, 2);

                        if sum(cur_expression)>0
                            % corr with substrate
                            y = kegg_sub_sum_intensities(select_path_idx(i),metabolomics_idx)';                          

                            [curcorr, curcorrP] = corr(cur_expression,y);
                            % add to matrix
                            curenzymes_mat(j,k) = curcorr;
                            curenzymes_p(j,k) = curcorrP;
                            
                             % calculate mdl coefficients and rsquared
                            mdl = fitlm(cur_expression,y);
                            curenzymes_mdlX(j,k) = mdl.Coefficients.Estimate(2);
                            curenzymes_mdlXp(j,k) = mdl.Coefficients.pValue(2);
                            curenzymes_mdlRsqadj(j,k) = mdl.Rsquared.Adjusted;
                        end
                        curenzymes_ec{j, k} = curec;
                        
                    end
                end
                % calculate corr with sum of all enzymes
                if ~isempty(horzcat(curenzymes_geneidx{:,k}))
                    cur_expression = table2array(countsMatrixGetMM_RNA_table(horzcat(curenzymes_geneidx{:,k}),gene_idx))';
                    cur_expression = sum(cur_expression, 2);
                    
                    if sum(cur_expression)>0
                        % corr with substrate
                        y = kegg_sub_sum_intensities(select_path_idx(i),metabolomics_idx)';                          

                        [curcorr, curcorrP] = corr(cur_expression, y);
                        % add to matrix
                        kegg_sub_prod_EC_sub_sum_corr(i,k) = curcorr;
                        kegg_sub_prod_EC_sub_sum_corrP(i,k) = curcorrP;
                        
                    end
                end
            end
            % keep max absolute value per species across all relevant enzymes
            [~, idx] = max(abs(curenzymes_mat), [], 1, 'linear');
            kegg_sub_prod_EC_sub_corr(i,:) = curenzymes_mat(idx);
            kegg_sub_prod_EC_sub_corrP(i,:) = curenzymes_p(idx);
            kegg_sub_prod_EC_sub_corr_id(i,:) = curenzymes_ec(idx);
            kegg_sub_prod_EC_sub_corr_geneidx(i,:) = curenzymes_geneidx(idx);
            kegg_sub_prod_EC_sub_mdlX(i,:)  = curenzymes_mdlX(idx);
            kegg_sub_prod_EC_sub_mdlP(i,:)  = curenzymes_mdlXp(idx);
            kegg_sub_prod_EC_sub_mdlRsqadj(i,:) = curenzymes_mdlRsqadj(idx);
        end
    end
end