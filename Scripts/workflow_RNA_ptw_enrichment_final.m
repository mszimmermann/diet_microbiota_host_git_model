%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot differential analysis results and perform pathway enrichment

% call script defining file dependenciesand global variables
addpath(genpath('.\'));
add_global_and_file_dependencies

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Required files:
% 'edgeR_gene_fold_changes_and_ann.csv'
% 'gene_annotation_kegg.tsv'
% 'gene_annotation_COG.tsv'
% 'gene_annotation_EC.tsv'
% Output: 
% Files:
% 'edgeR_gene_fold_changes_and_ann_filtered_FDRrecalculated.csv'
% 'ptwenr_recalc_KEGG_DOWN_updfiltered_eggnog_per_species.csv'
% 'ptwenr_recalc_KEGG_DOWN_updfiltered_eggnog_total.csv'
% 'ptwenr_recalc_KEGG_UP_updfiltered_eggnog_per_species.csv'
% 'ptwenr_recalc_KEGG_UP_updfiltered_eggnog_total.csv'
% 'ptwNchangingGenes_KEGG_DOWN_updfiltered_eggnog.csv'
% 'ptwNchangingGenes_KEGG_UP_updfiltered_eggnog.csv'
% Figures:
% 'fig_1e_bar_fraction_changing_genes_per_species_phylum_order.pdf'
% 'fig_1f_clustergram_ptwenr_KEGGDOWN_eggnog_updfilter_ANY_FDR0_1_ngenes3_perspeciesNgenes_noSum.pdf'
% 'fig_1f_clustergram_ptwenr_KEGGUP_eggnog_updfilter_ANY_FDR0_1_ngenes3_perspeciesNgenes_noSum.pdf'
% 'fig_sup_volcanos_edgeR_LRT_fc_HFD_CTR_RNA_per_species_detected_filtered_n5_rnaRESEQ.pdf'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read updated table with edgeR and deseq2 results instead of edgeR_gene_fold_changes_and_ann.csv
edgeRTable = readtable([outputFolder, ...
                    'edgeR_deseq2_gene_fold_changes_and_ann.csv']);

% filter out genes that are not detected in >5 samples
edgeRTable_filtered = edgeRTable(edgeRTable.geneFilter==1,:);

% recalculate FDR
edgeRTable_filtered.fdrHFDCTR_DNA = mafdr(edgeRTable_filtered.pHFDCTR_DNA, 'bhfdr', 1);
edgeRTable_filtered.fdrHFDCTR_RNA = mafdr(edgeRTable_filtered.pHFDCTR_RNA, 'bhfdr', 1);

% recalculate FDR for deseq2
edgeRTable_filtered.fdrDeseq2HFDCTR_DNA = mafdr(edgeRTable_filtered.pDeseq2HFDCTR_DNA, 'bhfdr', 1);
edgeRTable_filtered.fdrDeseq2HFDCTR_RNA = mafdr(edgeRTable_filtered.pDeseq2HFDCTR_RNA, 'bhfdr', 1);

% write filtered edgeR table with updated FDRs to file
writetable(edgeRTable_filtered, [outputFolder, ...
    'edgeR_deseq2_gene_fold_changes_and_ann_filtered_FDRrecalculated.csv']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read pathway-gene information
keggPathways_geneidx = readtable([inputFolderSeq,...
                            'gene_annotation_kegg.tsv'],...
                        'FileType','text');
keggPathways = unique(keggPathways_geneidx.PathwayID);

cogPathways_geneidx = readtable([inputFolderSeq,...
                            'gene_annotation_COG.tsv'],...
                        'FileType','text');
cogPathways = unique(cogPathways_geneidx.PathwayID);

ecPathways_geneidx = readtable([inputFolderSeq,...
                            'gene_annotation_EC.tsv'],...
                        'FileType','text');
ecPathways = unique(ecPathways_geneidx.PathwayID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get species names
edgeRFileNames_abbr = unique(edgeRTable_filtered.abbrSpecies);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot volcano for each species
fcThreshold = log2(1.5);
fdrThreshold = 0.05;

fig = figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:length(edgeRFileNames_abbr)
    curidx = cellfun(@(x) isequal(x,edgeRFileNames_abbr{i}),...
                            edgeRTable_filtered.abbrSpecies);
    sp = subplot(3,5,i);
  
    plotVolcano(edgeRTable_filtered.fcHFDCTR_RNA(curidx),...
                edgeRTable_filtered.fdrHFDCTR_RNA(curidx), [],...
                fcThreshold, fdrThreshold, [],sp);

    title(edgeRFileNames_abbr{i})
end
orient landscape
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
    [figureFolder, ...
    'fig_sup_volcanos_edgeR_LRT_fc_HFD_CTR_RNA_per_species_detected_filtered_n5_rnaRESEQ.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% cap FDR at -10
fig = figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:length(edgeRFileNames_abbr)
    curidx = cellfun(@(x) isequal(x,edgeRFileNames_abbr{i}),...
                            edgeRTable_filtered.abbrSpecies);
    sp = subplot(3,5,i);
    plotFDR = edgeRTable_filtered.fdrHFDCTR_RNA(curidx);
    % cap at -log10(x)=10;
    plotFDR(plotFDR<10^(-10)) = 10^(-10);
    plotVolcano(edgeRTable_filtered.fcHFDCTR_RNA(curidx),...
                plotFDR, [],...
                fcThreshold, fdrThreshold, [],sp);
    %ylim([0 80])
    title(edgeRFileNames_abbr{i})
end
orient landscape
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
    [figureFolder, ...
    'fig_sup_volcanos_edgeR_LRT_fc_HFD_CTR_RNA_per_species_detected_filtered_n5_rnaRESEQ_capped.pdf'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot one volcano for all species
% figure
% plotVolcano(edgeRTable_filtered.fcHFDCTR_RNA,...
%                 edgeRTable_filtered.fdrHFDCTR_RNA,...
%                 [], fcThreshold, fdrThreshold, [],[]);
% orient landscape
% print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
%     [figureFolder, 'volcanos_edgeR_LRT_fc_HFD_CTR_RNA_combined_species_detected_filtered_n5_rnaRESEQ.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot number of changing genes per species
changingGenesUP = zeros(length(edgeRFileNames_abbr),1);
changingGenesDOWN = zeros(length(edgeRFileNames_abbr),1);
% calculate number of changing genes per deseq2
changingGenesDeseq2UP = zeros(length(edgeRFileNames_abbr),1);
changingGenesDeseq2DOWN = zeros(length(edgeRFileNames_abbr),1);
% calculate number of changing genes per deseq2 and edgeR
changingGenesEdgeRDeseq2UP = zeros(length(edgeRFileNames_abbr),1);
changingGenesEdgeRDeseq2DOWN = zeros(length(edgeRFileNames_abbr),1);
% save total number of genes
changingGenesTotalNum = zeros(length(edgeRFileNames_abbr),1);
for i=1:length(edgeRFileNames_abbr)
    curidx = cellfun(@(x) isequal(x,edgeRFileNames_abbr{i}),...
                            edgeRTable_filtered.abbrSpecies);
    changingGenesUP(i) = nnz( (edgeRTable_filtered.fcHFDCTR_RNA(curidx) >= fcThreshold) &...
                           (edgeRTable_filtered.fdrHFDCTR_RNA(curidx) <= fdrThreshold));

    changingGenesDOWN(i) = nnz( (edgeRTable_filtered.fcHFDCTR_RNA(curidx) <= -fcThreshold) &...
                           (edgeRTable_filtered.fdrHFDCTR_RNA(curidx) <= fdrThreshold));
    % deseq2
    changingGenesDeseq2UP(i) = nnz( (edgeRTable_filtered.fcDeseq2HFDCTR_RNA(curidx) >= fcThreshold) &...
                           (edgeRTable_filtered.fdrDeseq2HFDCTR_RNA(curidx) <= fdrThreshold));

    changingGenesDeseq2DOWN(i) = nnz( (edgeRTable_filtered.fcDeseq2HFDCTR_RNA(curidx) <= -fcThreshold) &...
                           (edgeRTable_filtered.fdrDeseq2HFDCTR_RNA(curidx) <= fdrThreshold));

    % deseq2 and edgeR
    changingGenesEdgeRDeseq2UP(i) = nnz( (edgeRTable_filtered.fcHFDCTR_RNA(curidx) >= fcThreshold) &...
                           (edgeRTable_filtered.fdrHFDCTR_RNA(curidx) <= fdrThreshold) & ...
                           (edgeRTable_filtered.fcDeseq2HFDCTR_RNA(curidx) >= fcThreshold) &...
                           (edgeRTable_filtered.fdrDeseq2HFDCTR_RNA(curidx) <= fdrThreshold));

    changingGenesEdgeRDeseq2DOWN(i) = nnz( (edgeRTable_filtered.fcHFDCTR_RNA(curidx) <= -fcThreshold) &...
                           (edgeRTable_filtered.fdrHFDCTR_RNA(curidx) <= fdrThreshold) & ...
                            (edgeRTable_filtered.fcDeseq2HFDCTR_RNA(curidx) <= -fcThreshold) &...
                           (edgeRTable_filtered.fdrDeseq2HFDCTR_RNA(curidx) <= fdrThreshold));

    % total unfiltered number of genes
    curidx = cellfun(@(x) isequal(x,edgeRFileNames_abbr{i}),...
                            edgeRTable_filtered.abbrSpecies);
    changingGenesTotalNum(i) = nnz(curidx);
end
changingGenesTotalNum(changingGenesTotalNum==0)=1;
% plot number of changing genes
% [~, sortidx] = sort(changingGenesUP+changingGenesDOWN);
% figure
% barh(-changingGenesDOWN(sortidx))
% hold on
% barh(changingGenesUP(sortidx))
% set(gca, 'YTick', 1:length(sortidx))
% set(gca, 'YTickLabel', edgeRFileNames_abbr(sortidx))
% xlabel('Number of changing genes (DOWN/UP)')
% orient landscape
% print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
%     [figureFolder 'bar_number_changing_genes_per_species.pdf'])
% % plot fraction of changing genes
% [~, sortidx] = sort((changingGenesUP+changingGenesDOWN)./changingGenesTotalNum);
% figure
% barh((-changingGenesDOWN(sortidx))./(changingGenesTotalNum(sortidx)))
% hold on
% barh(changingGenesUP(sortidx)./changingGenesTotalNum(sortidx))
% set(gca, 'YTick', 1:length(sortidx))
% set(gca, 'YTickLabel', edgeRFileNames_abbr(sortidx))
% xlabel('Fraction of changing genes (DOWN/UP)')
% xlim([-0.1 0.1])
% orient landscape
% print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
%     [figureFolder 'bar_fraction_changing_genes_per_species.pdf'])

% plot fraction of changing genes in the order of abundance clustergram
sortidx = fliplr([3 13 14 9 10 1 8 11 12 5 4 7 6 2]);
figure
barh((-changingGenesDOWN(sortidx))./(changingGenesTotalNum(sortidx)))
hold on
barh(changingGenesUP(sortidx)./changingGenesTotalNum(sortidx))
set(gca, 'YTick', 1:length(sortidx))
set(gca, 'YTickLabel', edgeRFileNames_abbr(sortidx))
xlabel('Fraction of changing genes (DOWN/UP)')
xlim([-0.1 0.1])
ylim([0.5 length(sortidx)+0.5])
orient landscape
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
    [figureFolder ...
    'fig_1e_bar_fraction_changing_genes_per_species_clustergram_order.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%
% plot deseq2 and edgeR comparison
sortidx = fliplr([3 13 14 9 10 1 8 11 12 5 4 7 6 2]);
figure
barh((-[changingGenesDOWN(sortidx) changingGenesDeseq2DOWN(sortidx) changingGenesEdgeRDeseq2DOWN(sortidx)])./...
      ([changingGenesTotalNum(sortidx) changingGenesTotalNum(sortidx) changingGenesTotalNum(sortidx)]))
hold on
barh([changingGenesUP(sortidx) changingGenesDeseq2UP(sortidx) changingGenesEdgeRDeseq2UP(sortidx)]./...
    [changingGenesTotalNum(sortidx) changingGenesTotalNum(sortidx) changingGenesTotalNum(sortidx)])
set(gca, 'YTick', 1:length(sortidx))
set(gca, 'YTickLabel', edgeRFileNames_abbr(sortidx))
xlabel('Fraction of changing genes (DOWN/UP)')
xlim([-0.1 0.1])
ylim([0.5 length(sortidx)+0.5])
legend({"EdgeR Down", "Deseq2 Down", "Both Down", "EdgeR UP", "Deseq2 UP", "Both Up"})
orient landscape
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
    [figureFolder ...
    'fig_1e_bar_fraction_changing_genes_per_species_clustergram_order_edgeRdeseq.pdf'])
%%%%%%%%%%%%%%%%%%%%%%%%

% plot fraction of changing genes in the order of phylum
%sortidx = fliplr([12 10 3 9 14 1 8 4 2 6 7 5 13 11]);
figure
barh((-changingGenesDOWN(sortidx))./(changingGenesTotalNum(sortidx)))
hold on
barh(changingGenesUP(sortidx)./changingGenesTotalNum(sortidx))
set(gca, 'YTick', 1:length(sortidx))
set(gca, 'YTickLabel', edgeRFileNames_abbr(sortidx))
xlabel('Fraction of changing genes (DOWN/UP)')
xlim([-0.1 0.1])
orient landscape
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
    [figureFolder ...
    'fig_1e_bar_fraction_changing_genes_per_species_phylum_order.pdf'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pathway enrichment analysis for all genes simultaneously
% kegg pathway enrichment
for ptw_i = 1:3
    switch ptw_i
        case 1
            pathwayNames = keggPathways; % goPathwaysNumbers;% ecPathways;% 
            pathway_genes = keggPathways_geneidx; % goPathwaysNumbers_geneidx;% ecPathways_geneidx; % 
        case 2
            pathwayNames = cogPathways;% goPathwaysNumbers;% ecPathways;% 
            pathway_genes = cogPathways_geneidx;%  goPathwaysNumbers_geneidx;% ecPathways_geneidx; % 
        case 3
            pathwayNames = ecPathways;% 
            pathway_genes = ecPathways_geneidx; % 
        case 4
            pathwayNames = figfamPathways; 
            pathway_genes = figfamPathways_geneidx;
        case 5
            pathwayNames = plfPathways;
            pathway_genes = plfPathways_geneidx;
        case 6
            pathwayNames = pgfPathways;
            pathway_genes = pgfPathways_geneidx;
    end
    % define change direction
    for change_i = 1:2
        %for functional group enrichment
        enrichmentScores = zeros(length(pathwayNames), 1);
        nenrichCol = 9;
        enrichmentTable = cell(length(pathwayNames), nenrichCol);
 
        switch change_i
            case 1
                changedir = 'DOWN';%'UP';
            case 2
                changedir = 'UP';
        end
        % calculate the enrichment score for each pathway for each group
        if isequal(changedir, 'UP')
            % changing_group = (edgeRTable_filtered.fcHFDCTR_RNA>=fcThreshold) &...
            %                  (edgeRTable_filtered.fdrHFDCTR_RNA<=fdrThreshold);
            % select genes changing both in edgeR and deseq2
            changing_group = (edgeRTable_filtered.fcHFDCTR_RNA>=fcThreshold) &...
                             (edgeRTable_filtered.fdrHFDCTR_RNA<=fdrThreshold) &...
                             (edgeRTable_filtered.fcDeseq2HFDCTR_RNA>=fcThreshold) &...
                             (edgeRTable_filtered.fdrDeseq2HFDCTR_RNA<=fdrThreshold);
            changing_fc = sort(unique(edgeRTable_filtered.fcHFDCTR_RNA(changing_group)),...
                            'descend');
        else
            % changing_group = (edgeRTable_filtered.fcHFDCTR_RNA<=-fcThreshold) &...
            %                  (edgeRTable_filtered.fdrHFDCTR_RNA<=fdrThreshold);
            % select genes changing both in edgeR and deseq2
            changing_group = (edgeRTable_filtered.fcHFDCTR_RNA<=-fcThreshold) &...
                              (edgeRTable_filtered.fdrHFDCTR_RNA<=fdrThreshold) &...
                              (edgeRTable_filtered.fcDeseq2HFDCTR_RNA<=-fcThreshold) &...
                              (edgeRTable_filtered.fdrDeseq2HFDCTR_RNA<=fdrThreshold);
            changing_fc = sort(unique(edgeRTable_filtered.fcHFDCTR_RNA(changing_group)),...
                            'ascend');
        end    
        changing_fdr = sort(unique(edgeRTable_filtered.fdrHFDCTR_RNA(changing_group)),...
                            'ascend');
        % this is N - total number of genes
        nDetectedGenes = length(changing_group);
        % comparison of very small numbers
        scoreround = 1000000000;
        for iPw = 1:length(pathwayNames)
            pw_genes = zeros(size(changing_group,1),1);
            cur_pathway_genes = pathway_genes.GeneIDXfiltered{iPw};
            cur_pathway_genes = cellfun(@(x) str2num(x), strsplit(cur_pathway_genes,';'), 'unif', 0);
            cur_pathway_genes = cell2mat(cur_pathway_genes(cellfun(@(x) ~isempty(x), cur_pathway_genes)));
            pw_genes(cur_pathway_genes)=1;
                
            bestScore = 2;
            %tic
            for fdrthres = length(changing_fdr):-1:1
                if isequal(changedir, 'UP')
                    %changing_group = (edgeRTable_filtered.fcHFDCTR_RNA>=fcThreshold) &...
                    %                 (edgeRTable_filtered.fdrHFDCTR_RNA<=changing_fdr(fdrthres));
                    % both edgeR and deseq2
                    changing_group = (edgeRTable_filtered.fcHFDCTR_RNA>=fcThreshold) &...
                                     (edgeRTable_filtered.fdrHFDCTR_RNA<=changing_fdr(fdrthres)) &...
                                     (edgeRTable_filtered.fcDeseq2HFDCTR_RNA>=fcThreshold) &...
                                     (edgeRTable_filtered.fdrDeseq2HFDCTR_RNA<=changing_fdr(fdrthres));
                else
                    % changing_group = (edgeRTable_filtered.fcHFDCTR_RNA<=-fcThreshold) &...
                    %                  (edgeRTable_filtered.fdrHFDCTR_RNA<=changing_fdr(fdrthres));
                    % both edgeR and deseq2
                    changing_group = (edgeRTable_filtered.fcHFDCTR_RNA<=-fcThreshold) &...
                                     (edgeRTable_filtered.fdrHFDCTR_RNA<=changing_fdr(fdrthres)) &...
                                     (edgeRTable_filtered.fcDeseq2HFDCTR_RNA<=-fcThreshold) &...
                                     (edgeRTable_filtered.fdrDeseq2HFDCTR_RNA<=changing_fdr(fdrthres));
                end
                % this is n - size of the group
                group_size = sum(changing_group);
                % this is K - number of pathway genes
                nPwDetectedGenes = length(cur_pathway_genes);
                % this is k - number of pathway genes which are in the group
                genesFound = sum(changing_group & pw_genes);

                % calculate P(x<=k) - hypergeometric distribution f(k,N,K,n)
                % where N - total number of genes, n - size of the group
                score = sum( hygepdf(genesFound:nPwDetectedGenes, nDetectedGenes,...
                             nPwDetectedGenes, group_size) );
                if round(score*scoreround)<=round(bestScore*scoreround)
                    enrichmentScores(iPw) = score;
                    enrichmentTable{iPw, 1} = score;
                    enrichmentTable{iPw, 3} = pathwayNames{iPw};
                    enrichmentTable{iPw, 4} = genesFound;
                    enrichmentTable{iPw, 5} = group_size;
                    enrichmentTable{iPw, 6} = nPwDetectedGenes;
                    enrichmentTable{iPw, 7} = nDetectedGenes;
                    enrichmentTable{iPw, 8} = changing_fdr(fdrthres);
                    enrichmentTable{iPw, 9} = fcThreshold*(isequal(changedir, 'UP')-(1-isequal(changedir, 'UP')));
                    bestScore=score;
                end
                % if there are no genes in the changing group, break
                if genesFound==0
                    break
                end
            end
            for fcthres = length(changing_fc):-1:1
                if isequal(changedir, 'UP')
                    % changing_group = (edgeRTable_filtered.fcHFDCTR_RNA>=changing_fc(fcthres)) &...
                    %                  (edgeRTable_filtered.fdrHFDCTR_RNA<=fdrThreshold);
                    % both edgeR and deseq2
                    changing_group = (edgeRTable_filtered.fcHFDCTR_RNA>=changing_fc(fcthres)) &...
                                     (edgeRTable_filtered.fdrHFDCTR_RNA<=fdrThreshold) &...
                                     (edgeRTable_filtered.fcDeseq2HFDCTR_RNA>=changing_fc(fcthres)) &...
                                     (edgeRTable_filtered.fdrDeseq2HFDCTR_RNA<=fdrThreshold);
                else
                    % changing_group = (edgeRTable_filtered.fcHFDCTR_RNA<=changing_fc(fcthres)) &...
                    %                  (edgeRTable_filtered.fdrHFDCTR_RNA<=fdrThreshold);
                    % both edgeR and deseq2
                    changing_group = (edgeRTable_filtered.fcHFDCTR_RNA<=changing_fc(fcthres)) &...
                                     (edgeRTable_filtered.fdrHFDCTR_RNA<=fdrThreshold) &...
                                     (edgeRTable_filtered.fcDeseq2HFDCTR_RNA<=changing_fc(fcthres)) &...
                                     (edgeRTable_filtered.fdrDeseq2HFDCTR_RNA<=fdrThreshold);
                end
                % this is n - size of the group
                group_size = sum(changing_group);
                % this is K - number of pathway genes
                nPwDetectedGenes = length(cur_pathway_genes);
                % this is k - number of pathway genes which are in the group
                genesFound = sum(changing_group & pw_genes);

                % calculate P(x<=k) - hypergeometric distribution f(k,N,K,n)
                % where N - total number of genes, n - size of the group
                score = sum( hygepdf(genesFound:nPwDetectedGenes, nDetectedGenes,...
                             nPwDetectedGenes, group_size) );
                if round(score*scoreround)<=round(bestScore*scoreround)
                    enrichmentScores(iPw) = score;
                    enrichmentTable{iPw, 1} = score;
                    enrichmentTable{iPw, 3} = pathwayNames{iPw};
                    enrichmentTable{iPw, 4} = genesFound;
                    enrichmentTable{iPw, 5} = group_size;
                    enrichmentTable{iPw, 6} = nPwDetectedGenes;
                    enrichmentTable{iPw, 7} = nDetectedGenes;
                    enrichmentTable{iPw, 8} = fdrThreshold;
                    enrichmentTable{iPw, 9} = changing_fc(fcthres);
                    bestScore=score;
                end
                % if there are no genes in the changing group, break
                if genesFound==0
                    break
                end
            end
            %fprintf('Pathway %d of %d\n', iPw, length(pathwayNames));
            %toc
        end

        % multiple hypothesis adjustment
        padj = my_bhfdr(cell2mat(enrichmentTable(:,1)));
        for i=1:length(padj)
            enrichmentTable{i,2} = padj(i);
        end

        switch ptw_i
            case 1
                switch changedir
                    case 'DOWN'
                        enrichmentTable_KEGG_DOWN = enrichmentTable;
                    case 'UP'
                        enrichmentTable_KEGG_UP = enrichmentTable;
                end
            case 2
                switch changedir
                    case 'DOWN'
                        %enrichmentTable_GO_DOWN = enrichmentTable;
                        enrichmentTable_COG_DOWN = enrichmentTable;
                    case 'UP'
                        %enrichmentTable_GO_UP = enrichmentTable;
                        enrichmentTable_COG_UP = enrichmentTable;
                end
            case 3
                switch changedir
                    case 'DOWN'
                        enrichmentTable_EC_DOWN = enrichmentTable;
                    case 'UP'
                        enrichmentTable_EC_UP = enrichmentTable;
                end
            case 4
                switch changedir
                    case 'DOWN'
                        enrichmentTable_FIGFAM_DOWN = enrichmentTable;
                    case 'UP'
                        enrichmentTable_FIGFAM_UP = enrichmentTable;
                end
            case 5
                switch changedir
                    case 'DOWN'
                        enrichmentTable_PLF_DOWN = enrichmentTable;
                    case 'UP'
                        enrichmentTable_PLF_UP = enrichmentTable;
                end
            case 6
                switch changedir
                    case 'DOWN'
                        enrichmentTable_PGF_DOWN = enrichmentTable;
                    case 'UP'
                        enrichmentTable_PGF_UP = enrichmentTable;
                end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % pathway enrichment analysis per species
        ncol = 8;
        enrichmentTable = cell(length(pathwayNames), ncol*length(edgeRFileNames_abbr));
     
        for abbr_i = 1:length(edgeRFileNames_abbr)
            % calculate the enrichment score for each pathway for each group
            if isequal(changedir, 'UP')
                % changing_group = (edgeRTable_filtered.fcHFDCTR_RNA>=fcThreshold) &...
                %                  (edgeRTable_filtered.fdrHFDCTR_RNA<=fdrThreshold);
                % both edgeR and deseq2
                changing_group = (edgeRTable_filtered.fcHFDCTR_RNA>=fcThreshold) &...
                                 (edgeRTable_filtered.fdrHFDCTR_RNA<=fdrThreshold) &...
                                 (edgeRTable_filtered.fcDeseq2HFDCTR_RNA>=fcThreshold) &...
                                 (edgeRTable_filtered.fdrDeseq2HFDCTR_RNA<=fdrThreshold);
                changing_fc = sort(unique(edgeRTable_filtered.fcHFDCTR_RNA(changing_group)),...
                            'descend');
            else
                % changing_group = (edgeRTable_filtered.fcHFDCTR_RNA<=-fcThreshold) &...
                %                  (edgeRTable_filtered.fdrHFDCTR_RNA<=fdrThreshold);
                % both edgeR and deseq2
                changing_group = (edgeRTable_filtered.fcHFDCTR_RNA<=-fcThreshold) &...
                                 (edgeRTable_filtered.fdrHFDCTR_RNA<=fdrThreshold);
                                 (edgeRTable_filtered.fcDeseq2HFDCTR_RNA<=-fcThreshold) &...
                                 (edgeRTable_filtered.fdrDeseq2HFDCTR_RNA<=fdrThreshold);
                changing_fc = sort(unique(edgeRTable_filtered.fcHFDCTR_RNA(changing_group)),...
                            'ascend');
            end    
            % select genes of one species
            species_genes = ismember(edgeRTable_filtered.abbrSpecies, edgeRFileNames_abbr{abbr_i});
            changing_group = changing_group & species_genes;

            changing_fdr = sort(unique(edgeRTable_filtered.fdrHFDCTR_RNA(changing_group)),...
                                'ascend');
            % this is N - total number of genes
            nDetectedGenes = nnz(species_genes);
            % comparison of very small numbers
            scoreround = 1000000000;
            for iPw = 1:length(pathwayNames)
                pw_genes = zeros(size(changing_group,1),1);
                cur_pathway_genes = pathway_genes.GeneIDXfiltered{iPw};
                cur_pathway_genes = cellfun(@(x) str2num(x), strsplit(cur_pathway_genes,';'), 'unif', 0);
                cur_pathway_genes = cell2mat(cur_pathway_genes(cellfun(@(x) ~isempty(x), cur_pathway_genes)));
                pw_genes(cur_pathway_genes)=1;
                
                bestScore = 2;
                %tic
                for fdrthres = length(changing_fdr):-1:1
                    if isequal(changedir, 'UP')
                        % changing_group = (edgeRTable_filtered.fcHFDCTR_RNA>=fcThreshold) &...
                        %                  (edgeRTable_filtered.fdrHFDCTR_RNA<=fdrThreshold);
                        % both edgeR and deseq2
                        changing_group = (edgeRTable_filtered.fcHFDCTR_RNA>=fcThreshold) &...
                                         (edgeRTable_filtered.fdrHFDCTR_RNA<=fdrThreshold) &...
                                         (edgeRTable_filtered.fcDeseq2HFDCTR_RNA>=fcThreshold) &...
                                         (edgeRTable_filtered.fdrDeseq2HFDCTR_RNA<=fdrThreshold);
                    else
                        % changing_group = (edgeRTable_filtered.fcHFDCTR_RNA<=-fcThreshold) &...
                        %                  (edgeRTable_filtered.fdrHFDCTR_RNA<=changing_fdr(fdrthres));
                        % both edgeR and deseq2
                        changing_group = (edgeRTable_filtered.fcHFDCTR_RNA<=-fcThreshold) &...
                                         (edgeRTable_filtered.fdrHFDCTR_RNA<=changing_fdr(fdrthres)) & ...
                                         (edgeRTable_filtered.fcDeseq2HFDCTR_RNA<=-fcThreshold) &...
                                         (edgeRTable_filtered.fdrDeseq2HFDCTR_RNA<=changing_fdr(fdrthres));
                    end
                    changing_group = changing_group & species_genes;
                    % this is n - size of the group
                    group_size = sum(changing_group);
                    % this is K - number of pathway genes
                    nPwDetectedGenes = length(intersect(cur_pathway_genes, find(species_genes)));
                    % this is k - number of pathway genes which are in the group
                    genesFound = sum(changing_group & pw_genes);

                    % calculate P(x<=k) - hypergeometric distribution f(k,N,K,n)
                    % where N - total number of genes, n - size of the group
                    score = sum( hygepdf(genesFound:nPwDetectedGenes, nDetectedGenes,...
                                 nPwDetectedGenes, group_size) );
                    if round(score*scoreround)<round(bestScore*scoreround)
                        enrichmentTable{iPw, (abbr_i-1)*ncol+1} = score;
                        enrichmentTable{iPw, (abbr_i-1)*ncol+3} = genesFound;
                        enrichmentTable{iPw, (abbr_i-1)*ncol+4} = group_size;
                        enrichmentTable{iPw, (abbr_i-1)*ncol+5} = nPwDetectedGenes;
                        enrichmentTable{iPw, (abbr_i-1)*ncol+6} = nDetectedGenes;
                        enrichmentTable{iPw, (abbr_i-1)*ncol+7} = changing_fdr(fdrthres);
                        enrichmentTable{iPw, (abbr_i-1)*ncol+8} = fcThreshold*(isequal(changedir, 'UP')-(1-isequal(changedir, 'UP')));
                        bestScore=score;
                    end
                    % if there are no genes in the changing group, break
                    if genesFound==0
                        break;
                    end
                end
                for fcthres = length(changing_fdr):-1:1
                    if isequal(changedir, 'UP')
                        % changing_group = (edgeRTable_filtered.fcHFDCTR_RNA>=fcThreshold) &...
                        %                  (edgeRTable_filtered.fdrHFDCTR_RNA<=fdrThreshold);
                        % both edgeR and deseq2
                        changing_group = (edgeRTable_filtered.fcHFDCTR_RNA>=fcThreshold) &...
                                         (edgeRTable_filtered.fdrHFDCTR_RNA<=fdrThreshold) &...
                                         (edgeRTable_filtered.fcDeseq2HFDCTR_RNA>=fcThreshold) &...
                                         (edgeRTable_filtered.fdrDeseq2HFDCTR_RNA<=fdrThreshold);
                    else
                        % changing_group = (edgeRTable_filtered.fcHFDCTR_RNA<=-fcThreshold) &...
                        %                  (edgeRTable_filtered.fdrHFDCTR_RNA<=changing_fdr(fdrthres));
                        % both edgeR and deseq2
                        changing_group = (edgeRTable_filtered.fcHFDCTR_RNA<=-fcThreshold) &...
                                         (edgeRTable_filtered.fdrHFDCTR_RNA<=changing_fdr(fdrthres)) &...
                                         (edgeRTable_filtered.fcDeseq2HFDCTR_RNA<=-fcThreshold) &...
                                         (edgeRTable_filtered.fdrDeseq2HFDCTR_RNA<=changing_fdr(fdrthres));
                    end
                    changing_group = changing_group & species_genes;
                    % this is n - size of the group
                    group_size = sum(changing_group);
                    % this is K - number of pathway genes
                    nPwDetectedGenes = length(intersect(cur_pathway_genes, find(species_genes)));
                    % this is k - number of pathway genes which are in the group
                    genesFound = sum(changing_group & pw_genes);

                    % calculate P(x<=k) - hypergeometric distribution f(k,N,K,n)
                    % where N - total number of genes, n - size of the group
                    score = sum( hygepdf(genesFound:nPwDetectedGenes, nDetectedGenes,...
                                 nPwDetectedGenes, group_size) );
                    if round(score*scoreround)<round(bestScore*scoreround)
                        enrichmentTable{iPw, (abbr_i-1)*ncol+1} = score;
                        enrichmentTable{iPw, (abbr_i-1)*ncol+3} = genesFound;
                        enrichmentTable{iPw, (abbr_i-1)*ncol+4} = group_size;
                        enrichmentTable{iPw, (abbr_i-1)*ncol+5} = nPwDetectedGenes;
                        enrichmentTable{iPw, (abbr_i-1)*ncol+6} = nDetectedGenes;
                        enrichmentTable{iPw, (abbr_i-1)*ncol+7} = fdrThreshold;
                        enrichmentTable{iPw, (abbr_i-1)*ncol+8} = changing_fc(fcthres);
                        bestScore=score;
                    end
                    % if there are no genes in the changing group, break
                    if genesFound==0
                        break;
                    end
                end
                %fprintf('Pathway %d of %d\n', iPw, length(pathwayNames));
                %toc
            end

            % multiple hypothesis adjustment
            padj = my_bhfdr(cell2mat(enrichmentTable(:,(abbr_i-1)*ncol+1)));
            for i=1:length(padj)
                enrichmentTable{i,(abbr_i-1)*ncol+2} = padj(i);
            end
            fprintf('Finished enrichment for %s\n', edgeRFileNames_abbr{abbr_i});
        end
        % add pathway names
        enrichmentTable = [pathwayNames enrichmentTable];
         switch ptw_i
            case 1
                switch changedir
                    case 'DOWN'
                        enrichmentTable_KEGG_DOWN_perspecies = enrichmentTable;
                    case 'UP'
                        enrichmentTable_KEGG_UP_perspecies = enrichmentTable;
                end
            case 2
                switch changedir
                    case 'DOWN'
                        enrichmentTable_COG_DOWN_perspecies = enrichmentTable;
                    case 'UP'
                        enrichmentTable_COG_UP_perspecies = enrichmentTable;
                end
            case 3
                switch changedir
                    case 'DOWN'
                        enrichmentTable_EC_DOWN_perspecies = enrichmentTable;
                    case 'UP'
                        enrichmentTable_EC_UP_perspecies = enrichmentTable;
                end
            case 4
                switch changedir
                    case 'DOWN'
                        enrichmentTable_FIGFAM_DOWN_perspecies = enrichmentTable;
                    case 'UP'
                        enrichmentTable_FIGFAM_UP_perspecies = enrichmentTable;
                end
            case 5
                switch changedir
                    case 'DOWN'
                        enrichmentTable_PLF_DOWN_perspecies = enrichmentTable;
                    case 'UP'
                        enrichmentTable_PLF_UP_perspecies = enrichmentTable;
                end
            case 6
                switch changedir
                    case 'DOWN'
                        enrichmentTable_PGF_DOWN_perspecies = enrichmentTable;
                    case 'UP'
                        enrichmentTable_PGF_UP_perspecies = enrichmentTable;
                end
         end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save enrichment tables to file
enrichmentColumns = {'pvalue', 'FDR',...
                     'n_pw_changin', 'n_changing', 'n_pw', 'n_total',...
                     'p_threshold', 'fc_threshold'};
enrichmentColumns_species = strcat(edgeRFileNames_abbr(kron(1:length(edgeRFileNames_abbr),...
                                                            ones(1,length(enrichmentColumns))))',...
                                   '_',...
                                   enrichmentColumns(kron(ones(1,length(edgeRFileNames_abbr)),...
                                                          1:length(enrichmentColumns))));

enrichmentColumns = [{'pathwayID'}, enrichmentColumns];
enrichmentColumns_species = [{'pathwayID'}, enrichmentColumns_species];

for ptw_i = 1:3
    switch ptw_i
            case 1
                enrichmentTable_DOWN_perspecies = enrichmentTable_KEGG_DOWN_perspecies;
                enrichmentTable_UP_perspecies = enrichmentTable_KEGG_UP_perspecies;
                enrichmentTable_DOWN = enrichmentTable_KEGG_DOWN(:,[3 1 2 4:end]);
                enrichmentTable_UP = enrichmentTable_KEGG_UP(:,[3 1 2 4:end]);
                enrichmentFileName = 'ptwenr_recalc_KEGG_';
            case 2
                enrichmentTable_DOWN_perspecies = enrichmentTable_COG_DOWN_perspecies;% GO_DOWN_perspecies;
                enrichmentTable_UP_perspecies = enrichmentTable_COG_UP_perspecies;% GO_UP_perspecies;
                enrichmentTable_DOWN = enrichmentTable_COG_DOWN(:,[3 1 2 4:end]);
                enrichmentTable_UP = enrichmentTable_COG_UP(:,[3 1 2 4:end]);
                enrichmentFileName = 'ptwenr_recalc_COG_';%GO_';
            case 3
                enrichmentTable_DOWN_perspecies = enrichmentTable_EC_DOWN_perspecies;
                enrichmentTable_UP_perspecies = enrichmentTable_EC_UP_perspecies;
                enrichmentTable_DOWN = enrichmentTable_EC_DOWN(:,[3 1 2 4:end]);
                enrichmentTable_UP = enrichmentTable_EC_UP(:,[3 1 2 4:end]);
                enrichmentFileName = 'ptwenr_recalc_EC_';
    end
    enrichmentTable_DOWN_perspecies = cell2table(enrichmentTable_DOWN_perspecies,...
                                                'VariableNames', enrichmentColumns_species);
    enrichmentTable_UP_perspecies = cell2table(enrichmentTable_UP_perspecies,...
                                                'VariableNames', enrichmentColumns_species);
    enrichmentTable_DOWN = cell2table(enrichmentTable_DOWN,...
                                      'VariableNames', enrichmentColumns);
    enrichmentTable_UP = cell2table(enrichmentTable_UP,...
                                    'VariableNames', enrichmentColumns);

    % save to file
    writetable(enrichmentTable_DOWN_perspecies,...
               [outputFolder, enrichmentFileName, 'DOWN_updfiltered_eggnog_edgerdeseq_per_species.csv']);
    writetable(enrichmentTable_UP_perspecies,...
               [outputFolder, enrichmentFileName, 'UP_updfiltered_eggnog_edgerdeseq_per_species.csv']);
    writetable(enrichmentTable_DOWN,...
               [outputFolder, enrichmentFileName, 'DOWN_updfiltered_eggnog_edgerdeseq_total.csv']);
    writetable(enrichmentTable_UP,...
               [outputFolder, enrichmentFileName, 'UP_updfiltered_eggnog_edgerdeseq_total.csv']);
           
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate total number of changing genes per pathway (with the same
% threshold)
for ptw_i = 1:3
    switch ptw_i
        case 1
            pathwayNames = keggPathways; 
            pathway_genes = keggPathways_geneidx; 
        case 2
            pathwayNames = cogPathways;
            pathway_genes = cogPathways_geneidx;
        case 3
            pathwayNames = ecPathways;% 
            pathway_genes = ecPathways_geneidx; % 
        case 4
            pathwayNames = figfamPathways; 
            pathway_genes = figfamPathways_geneidx;
        case 5
            pathwayNames = plfPathways;
            pathway_genes = plfPathways_geneidx;
        case 6
            pathwayNames = pgfPathways;
            pathway_genes = pgfPathways_geneidx;
    end
    % define change direction
    for change_i = 1:2
        %for functional group enrichment
        pathwayNchagingGenes = zeros(length(pathwayNames), 1);
        
        switch change_i
            case 1
                changedir = 'DOWN';%'UP';
            case 2
                changedir = 'UP';
        end
        % calculate the enrichment score for each pathway for each group
        if isequal(changedir, 'UP')
            changing_group = (edgeRTable_filtered.fcHFDCTR_RNA>=fcThreshold) &...
                             (edgeRTable_filtered.fdrHFDCTR_RNA<=fdrThreshold);
        else
            changing_group = (edgeRTable_filtered.fcHFDCTR_RNA<=-fcThreshold) &...
                             (edgeRTable_filtered.fdrHFDCTR_RNA<=fdrThreshold);
        end    
        
        for iPw = 1:length(pathwayNames)
            pw_genes = zeros(size(changing_group,1),1);
            cur_pathway_genes = pathway_genes.GeneIDXfiltered{iPw};
            cur_pathway_genes = cellfun(@(x) str2num(x), strsplit(cur_pathway_genes,';'), 'unif', 0);
            cur_pathway_genes = cell2mat(cur_pathway_genes(cellfun(@(x) ~isempty(x), cur_pathway_genes)));
            pw_genes(cur_pathway_genes)=1;
            
            pathwayNchagingGenes(iPw) = sum(changing_group & pw_genes);
        end
        switch ptw_i
            case 1
                switch changedir
                    case 'DOWN'
                        pathwayNchagingGenes_KEGG_DOWN = pathwayNchagingGenes;
                    case 'UP'
                        pathwayNchagingGenes_KEGG_UP = pathwayNchagingGenes;
                end
            case 2
                switch changedir
                    case 'DOWN'
                        %enrichmentTable_GO_DOWN = enrichmentTable;
                        pathwayNchagingGenes_COG_DOWN = pathwayNchagingGenes;
                    case 'UP'
                        %enrichmentTable_GO_UP = enrichmentTable;
                        pathwayNchagingGenes_COG_UP = pathwayNchagingGenes;
                end
            case 3
                switch changedir
                    case 'DOWN'
                        pathwayNchagingGenes_EC_DOWN = pathwayNchagingGenes;
                    case 'UP'
                        pathwayNchagingGenes_EC_UP = pathwayNchagingGenes;
                end
            case 4
                switch changedir
                    case 'DOWN'
                        pathwayNchagingGenes_FIGFAM_DOWN = pathwayNchagingGenes;
                    case 'UP'
                        pathwayNchagingGenes_FIGFAM_UP = pathwayNchagingGenes;
                end
            case 5
                switch changedir
                    case 'DOWN'
                        pathwayNchagingGenes_PLF_DOWN = pathwayNchagingGenes;
                    case 'UP'
                        pathwayNchagingGenes_PLF_UP = pathwayNchagingGenes;
                end
            case 6
                switch changedir
                    case 'DOWN'
                        pathwayNchagingGenes_PGF_DOWN = pathwayNchagingGenes;
                    case 'UP'
                        pathwayNchagingGenes_PGF_UP = pathwayNchagingGenes;
                end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % number of changing genes per species
        pathwayNchagingGenes = zeros(length(pathwayNames), length(edgeRFileNames_abbr));
     
        for abbr_i = 1:length(edgeRFileNames_abbr)
            % calculate the enrichment score for each pathway for each group
            if isequal(changedir, 'UP')
                changing_group = (edgeRTable_filtered.fcHFDCTR_RNA>=fcThreshold) &...
                                 (edgeRTable_filtered.fdrHFDCTR_RNA<=fdrThreshold);
            else
                changing_group = (edgeRTable_filtered.fcHFDCTR_RNA<=-fcThreshold) &...
                                 (edgeRTable_filtered.fdrHFDCTR_RNA<=fdrThreshold);
            end    
            % select genes of one species
            species_genes = ismember(edgeRTable_filtered.abbrSpecies, edgeRFileNames_abbr{abbr_i});
            changing_group = changing_group & species_genes;

            for iPw = 1:length(pathwayNames)
                pw_genes = zeros(size(changing_group,1),1);
                cur_pathway_genes = pathway_genes.GeneIDXfiltered{iPw};
                cur_pathway_genes = cellfun(@(x) str2num(x), strsplit(cur_pathway_genes,';'), 'unif', 0);
                cur_pathway_genes = cell2mat(cur_pathway_genes(cellfun(@(x) ~isempty(x), cur_pathway_genes)));
                pw_genes(cur_pathway_genes)=1;
                
                pathwayNchagingGenes(iPw, abbr_i) = sum(changing_group & pw_genes);
            end

        end
        % add pathway names
        switch ptw_i
            case 1
                switch changedir
                    case 'DOWN'
                        pathwayNchagingGenes_KEGG_DOWN_perspecies = pathwayNchagingGenes;
                    case 'UP'
                        pathwayNchagingGenes_KEGG_UP_perspecies = pathwayNchagingGenes;
                end
            case 2
                switch changedir
                    case 'DOWN'
                        %enrichmentTable_GO_DOWN_perspecies = enrichmentTable;
                        pathwayNchagingGenes_COG_DOWN_perspecies = pathwayNchagingGenes;
                    case 'UP'
                        %enrichmentTable_GO_UP_perspecies = enrichmentTable;
                        pathwayNchagingGenes_COG_UP_perspecies = pathwayNchagingGenes;
                end
            case 3
                switch changedir
                    case 'DOWN'
                        pathwayNchagingGenes_EC_DOWN_perspecies = pathwayNchagingGenes;
                    case 'UP'
                        pathwayNchagingGenes_EC_UP_perspecies = pathwayNchagingGenes;
                end
            case 4
                switch changedir
                    case 'DOWN'
                        pathwayNchagingGenes_FIGFAM_DOWN_perspecies = pathwayNchagingGenes;
                    case 'UP'
                        pathwayNchagingGenes_FIGFAM_UP_perspecies = pathwayNchagingGenes;
                end
            case 5
                switch changedir
                    case 'DOWN'
                        pathwayNchagingGenes_PLF_DOWN_perspecies = pathwayNchagingGenes;
                    case 'UP'
                        pathwayNchagingGenes_PLF_UP_perspecies = pathwayNchagingGenes;
                end
            case 6
                switch changedir
                    case 'DOWN'
                        pathwayNchagingGenes_PGF_DOWN_perspecies = pathwayNchagingGenes;
                    case 'UP'
                        pathwayNchagingGenes_PGF_UP_perspecies = pathwayNchagingGenes;
                end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save number of genes per pathway to file
pathwayNgenesColumns_species = [{'pathwayID'}, {'Total'}, edgeRFileNames_abbr'];

for ptw_i = 1:3
    switch ptw_i
            case 1
                pathwayNchangingGenes_DOWN_perspecies = [pathwayNchagingGenes_KEGG_DOWN,...
                                                         pathwayNchagingGenes_KEGG_DOWN_perspecies];
                pathwayNchangingGenes_UP_perspecies = [ pathwayNchagingGenes_KEGG_UP,...
                                                        pathwayNchagingGenes_KEGG_UP_perspecies];
                pathwayNames = keggPathways;
                enrichmentFileName = 'ptwNchangingGenes_KEGG_';
            case 2
                pathwayNchangingGenes_DOWN_perspecies = [pathwayNchagingGenes_COG_DOWN,...
                                                         pathwayNchagingGenes_COG_DOWN_perspecies];
                pathwayNchangingGenes_UP_perspecies = [pathwayNchagingGenes_COG_UP,...
                                                       pathwayNchagingGenes_COG_UP_perspecies];
                pathwayNames = cogPathways;
                enrichmentFileName = 'ptwNchangingGenes_COG_';%GO_';
            case 3
                pathwayNchangingGenes_DOWN_perspecies = [pathwayNchagingGenes_EC_DOWN,...
                                                         pathwayNchagingGenes_EC_DOWN_perspecies];
                pathwayNchangingGenes_UP_perspecies = [pathwayNchagingGenes_EC_UP,...
                                                       pathwayNchagingGenes_EC_UP_perspecies];
                pathwayNames = ecPathways;
                enrichmentFileName = 'ptwNchangingGenes_EC_';
    end
    pathwayNchangingGenes_DOWN_perspecies = [cell2table(pathwayNames, ...
                                                'VariableNames', pathwayNgenesColumns_species(1)),...
                                                 array2table(pathwayNchangingGenes_DOWN_perspecies,...
                                                'VariableNames', pathwayNgenesColumns_species(2:end))];
    pathwayNchangingGenes_UP_perspecies = [cell2table(pathwayNames, ...
                                                'VariableNames', pathwayNgenesColumns_species(1)),...
                                                 array2table(pathwayNchangingGenes_UP_perspecies,...
                                                'VariableNames', pathwayNgenesColumns_species(2:end))];

    % save to file
    writetable(pathwayNchangingGenes_DOWN_perspecies,...
               [outputFolder, enrichmentFileName, 'DOWN_updfiltered_eggnog_edgerdeseq.csv']);
    writetable(pathwayNchangingGenes_UP_perspecies,...
               [outputFolder, enrichmentFileName, 'UP_updfiltered_eggnog_edgerdeseq.csv']);
           
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot enrichment results
% get enrichment results
enrichmentFileName = 'ptwenr_recalc_KEGG_';
enrichmentTable_KEGG_DOWN_perspecies = readtable([outputFolder, enrichmentFileName, 'DOWN_updfiltered_eggnog_edgerdeseq_per_species.csv']);
enrichmentTable_KEGG_UP_perspecies = readtable([outputFolder, enrichmentFileName, 'UP_updfiltered_eggnog_edgerdeseq_per_species.csv']);
enrichmentTable_KEGG_DOWN = readtable([outputFolder, enrichmentFileName, 'DOWN_updfiltered_eggnog_edgerdeseq_total.csv']);
enrichmentTable_KEGG_UP = readtable([outputFolder, enrichmentFileName, 'UP_updfiltered_eggnog_edgerdeseq_total.csv']);
enrichmentFileName = 'ptwNchangingGenes_KEGG_';
pathwayNchangingGenes_KEGG_DOWN_perspecies = readtable([outputFolder, enrichmentFileName, 'DOWN_updfiltered_eggnog_edgerdeseq.csv']);
pathwayNchangingGenes_KEGG_UP_perspecies = readtable([outputFolder, enrichmentFileName, 'UP_updfiltered_eggnog_edgerdeseq.csv']);
edgeRFileNames_abbr = pathwayNchangingGenes_KEGG_DOWN_perspecies.Properties.VariableNames';
edgeRFileNames_abbr(cellfun(@(x) contains(x, 'pathway') | contains(x, 'Total'), edgeRFileNames_abbr))=[];

% get pathway names
kegg_ptw_names = table2cell(readtable([rawdataFolder...
                                    'kegg_ptw_names_and_bact_flag.csv'],...
                                    'ReadVariableNames', 1));
kegg_ptw_names(:,1) = cellfun(@(x) strrep(x, 'path:',''),kegg_ptw_names(:,1), 'unif',0);

pathwayNames = cellfun(@(x) strcat(x,': ', kegg_ptw_names(ismember(kegg_ptw_names(:,1),x),2)), keggPathways, 'unif',0);
pathwayNames(cellfun(@(x) ~isempty(x),pathwayNames)) = cellfun(@(x) x{1},pathwayNames(cellfun(@(x) ~isempty(x),pathwayNames)), 'unif', 0);
pathwayNames(cellfun(@(x) isempty(x),pathwayNames)) = keggPathways(cellfun(@(x) isempty(x),pathwayNames));

pathwayFlag = zeros(size(pathwayNames));
for i = 1:length(pathwayFlag)
    if nnz(ismember(kegg_ptw_names(:,1),keggPathways{i}))
        pathwayFlag(i) = kegg_ptw_names{ismember(kegg_ptw_names(:,1),keggPathways{i}),3};
    end
end

% define threshold for visualisation
fdrThresholdplot = 0.1;% 0.01;%
fdrThresholdplot_all = 0.1;% 0.01;%
ngenesThreshold = 3;
clustercoef = 1;
ncol = 8;
cgo = cell(2,1);
for dir_i = 1:2
    switch dir_i
        case 1
            plotenr = [enrichmentTable_KEGG_DOWN(:,2) ...
                enrichmentTable_KEGG_DOWN_perspecies(:, 3:ncol:end)];
            plotenr_ngenes = [enrichmentTable_KEGG_DOWN(:,4)...
                enrichmentTable_KEGG_DOWN_perspecies(:, 4:ncol:end)];
            plotenr_ngenes_total = pathwayNchangingGenes_KEGG_DOWN_perspecies(:,2:end);
            clustercoef = -1;
        case 2
            plotenr = [enrichmentTable_KEGG_UP(:,2) ...
                enrichmentTable_KEGG_UP_perspecies(:, 3:ncol:end)];
            plotenr_ngenes = [enrichmentTable_KEGG_UP(:,4)...
                enrichmentTable_KEGG_UP_perspecies(:, 4:ncol:end)];
            plotenr_ngenes_total = pathwayNchangingGenes_KEGG_UP_perspecies(:,2:end);
            clustercoef = 1;
    end

    
    % replace empty with 1
    %plotenr(cellfun(@(x) isempty(x), plotenr))={1};
    % convert to matrix
    plotenr = table2array(plotenr);
    
    % replace empty genes with 0
    %plotenr_ngenes(cellfun(@(x) isempty(x), plotenr_ngenes))={0};
    plotenr_ngenes = table2array(plotenr_ngenes);
    
    % add rows and column names
    plotenrrows = pathwayNames;
    plotenrcols = ["All"; edgeRFileNames_abbr(1:end)];
    if size(plotenr,2) > length(plotenrcols)
        plotenr(:,end) = [];
    end

    %keep_rows = sum(plotenr <=fdrThresholdplot,2)>0;
    % keep only pathways enriched in all
%     keep_rows = (plotenr(:,1) <=fdrThresholdplot) &...
%                 (plotenr_ngenes(:,1) >= ngenesThreshold);
    % species-specific enrichment
%     keep_rows = (plotenr(:,1) > fdrThresholdplot_all) &...
%                 sum((plotenr <=fdrThresholdplot).*...
%                     (plotenr_ngenes >= ngenesThreshold),2)>0;
% 
    % keep both all and species-specific
%     keep_rows = (sum((plotenr <=fdrThresholdplot).*...
%                     (plotenr_ngenes >= ngenesThreshold),2)>0) &...
%                     (pathwayFlag==0);
    % no total
    keep_rows = (sum((plotenr(:,2:end) <=fdrThresholdplot).*...
                    (plotenr_ngenes(:,2:end) >= ngenesThreshold),2)>0) &...
                    (pathwayFlag==0);


    % (un)comment to plot FDR
    %clustermat = -log10(plotenr);
    %clustermat(clustermat<-log10(fdrThresholdplot))=0;
    %clustermat = clustermat*clustercoef;
    
    % (un)comment to plot number of genes
    clustermat = plotenr_ngenes*clustercoef;
    clusterrows = plotenrrows;
    %clustersort = [1 1+[12 11 8 1 10 9 14 13 3 4 5 7 6 2]];
    %clustersort = [1 1+[12 10 3 9 14 1 8 4 2 6 7 5 13 11]];
    % with total
    %clustersort = [1 1+[3 13 14 9 10 1 8 11 12 5 4 7 6 2]];
    % without total
    clustersort = [1+[3 13 14 9 10 1 8 11 12 5 4 7 6 2]];
    clustercols = plotenrcols;
    clusterkeep = keep_rows;
    

    cgo{dir_i} = clustergram(clustermat(clusterkeep,clustersort),...
                'ColumnLabels', clustercols(clustersort),...
                'RowLabels', clusterrows(clusterkeep), ...
                'Symmetric', 1,...
                'Colormap', redbluecmap,...parula,...gray,...
                'Cluster', 'column',...'all',...
                'DisplayRange', 10,...
                'Annotate', 'on');%bone)
end
fig = cgo{1}.plot;
orient landscape
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
    [ figureFolder, ...
    'fig_1f_clustergram_ptwenr_KEGGDOWN_eggnog_edgerdeseq_updfilter_ANY_FDR0_1_ngenes3_perspeciesNgenes_noSum.pdf']);
fig = cgo{2}.plot;
orient landscape
print(gcf, '-painters', '-dpdf', '-r600', '-bestfit',...
    [ figureFolder, ...
    'fig_1f_clustergram_ptwenr_KEGGUP_eggnog_edgerdeseq_updfilter_ANY_FDR0_1_ngenes3_perspeciesNgenes_noSum.pdf']);
