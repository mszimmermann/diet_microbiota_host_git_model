org_dir <- '.\\InputData\\sequencing_data\\ballgown_RNA\\eggNOGann\\'
output_dir <- '.\\ProcessedData\\DEGoutput\\'

org_files = list.files(org_dir)
#setwd(org_dir)
# get raw read files

if(grepl('DNA',org_dir)){
  
  org_files = org_files[grepl("stringtie_DNA", org_files)]
  org_files = org_files[startsWith(org_files, "stringtie_DNA")]
}else{
  org_files = org_files[grepl("stringtie_RNA100", org_files)]
  org_files = org_files[startsWith(org_files, "stringtie_RNA100")]
}
# remove all zero files
org_files = org_files[!grepl("stringtie_DNA_Bhyd_t_e_data_merged_", org_files)]
org_files = org_files[!grepl("stringtie_DNA_Blon_t_e_data_merged_", org_files)]
org_files = org_files[!grepl("stringtie_DNA_Erec_t_e_data_merged_", org_files)]

org_files = org_files[!grepl("stringtie_RNA100_Bhyd_t_e_data_merged_", org_files)]
org_files = org_files[!grepl("stringtie_RNA100_Blon_t_e_data_merged_", org_files)]
org_files = org_files[!grepl("stringtie_RNA100_Erec_t_e_data_merged_", org_files)]


########################################################
# perform deseq2 analysis and normalization
library("DESeq2")

mincount_threshold <- 20000

# for each organism, perform the same analysis
for (infilename in org_files){
  outfilename <- paste("compare_methods_deseq_results_hisatEB_hfd_vs_ctr_", infilename, sep="")
  normalizedfilename <- paste("compare_methods_deseq_normalized_", infilename, sep="")
  # read gene counts
  countData <- as.matrix(read.csv(paste0(org_dir, infilename)))#, row.names="gene_id"))
  rownames <- paste(countData[,colnames(countData)=="chr"], countData[,colnames(countData)=="gene_id"], sep="_")
  rownames(countData) <- rownames
  # replace NA with 0
  countData[is.na(countData)] <- 0
  # leave only data columns
  countData_counts <- countData[,grepl("MSZ", colnames(countData))]
  # convert to numeric
  storage.mode(countData_counts) <- "numeric"
  print(colnames(countData_counts))
  diet <- factor(c('HFD', 'HFD', 'HFD', 'HFD', 'HFD', 'CTR', 'CTR', 'CTR', 'CTR', 'CTR'))
  #################################################
  # remove samples where bacteria were not detected
  removeCols <- colSums(countData_counts) < mincount_threshold
  countData_counts <- countData_counts[,!removeCols]
  diet <- diet[!removeCols]
  #################################################
  if(length(diet) > 3){
    dds <- DESeqDataSetFromMatrix(countData = countData_counts, colData = DataFrame(diet), ~ diet)
    dds = DESeq(dds)
    res = results(dds)
    table_counts_normalized <- counts(dds, normalized=TRUE)
    # add gene annotations from the original file to results
    countData_ann <- countData[,!grepl("MSZ", colnames(countData))]
    res = cbind(res,countData_ann)
    table_counts_normalized = cbind(table_counts_normalized, countData_ann)
    # write to files
    write.table(res,paste0(output_dir,outfilename), sep="\t")
    write.table(table_counts_normalized,paste0(output_dir, normalizedfilename), sep="\t")
  }else{
    print(paste('Warning: File does not contain enough count data: ', infilename))
  }
}

########################################################
# perform edgeR analysis and normalization

library(edgeR)

for (infilename in org_files){
  outfilename1 <- paste("compare_methods_edgeR_QLF_results_hisatEB_hfd_vs_ctr_", infilename, sep="")
  outfilename2 <- paste("compare_methods_edgeR_LRT_results_hisatEB_hfd_vs_ctr_", infilename, sep="")
  outfilename3 <- paste("compare_methods_edgeR_Exact_results_hisatEB_hfd_vs_ctr_", infilename, sep="")
  normalizedfilename <- paste("compare_methods_edgeR_cpm_normalized_", infilename, sep="")
  # read gene counts
  countData <- as.matrix(read.csv(paste0(org_dir,infilename)))#, row.names="gene_id"))
  rownames <- paste(countData[,colnames(countData)=="chr"], countData[,colnames(countData)=="gene_id"], sep="_")
  rownames(countData) <- rownames
  # replace NA with 0
  countData[is.na(countData)] <- 0
  # leave only data columns
  countData_counts <- countData[,grepl("MSZ", colnames(countData))]
  # convert to numeric
  storage.mode(countData_counts) <- "numeric"
  print(colnames(countData_counts))
  diet <- factor(c('HFD', 'HFD', 'HFD', 'HFD', 'HFD', 'CTR', 'CTR', 'CTR', 'CTR', 'CTR'))
  #################################################
  # remove samples where bacteria were not detected
  removeCols <- colSums(countData_counts) < mincount_threshold
  countData_counts <- countData_counts[,!removeCols]
  diet <- diet[!removeCols]
  #################################################
  if(length(diet) > 3){
    # Creating a DGEList object
    dgList <- DGEList(counts=countData_counts, genes=rownames(countData_counts), group=diet)
    #normalise RNA-seq with the rimmed mean of M-values (TMM) method.
    dgList <- calcNormFactors(dgList, method="TMM")
    design <- model.matrix(~diet)
    dgList <- estimateDisp(dgList,design)
    # To perform quasi-likelihood F-tests
    fit <- glmQLFit(dgList,design)
    qlf <- glmQLFTest(fit,coef=2)
    edgeR_result_QLF <- topTags(qlf,n=1000000)
    # To perform likelihood ratio tests:
    fit <- glmFit(dgList,design)
    lrt <- glmLRT(fit,coef=2)
    edgeR_result_LRT <- topTags(lrt,n=1000000)
    # To perform exact test
    de <- exactTest(dgList)
    edgeR_result_Exact <- topTags(de,n=1000000)
    table_counts_normalized <- cpm(dgList)
    # add gene annotations from the original file to results
    countData_ann <- countData[,!grepl("MSZ", colnames(countData))]
    edgeR_result_QLF$table = cbind(edgeR_result_QLF$table,
                                   countData_ann[match(rownames(edgeR_result_QLF$table),
                                                       rownames(countData_ann)),])
    edgeR_result_LRT$table = cbind(edgeR_result_LRT$table,
                                   countData_ann[match(rownames(edgeR_result_LRT$table),
                                                       rownames(countData_ann)),])
    edgeR_result_Exact$table = cbind(edgeR_result_Exact$table,
                                     countData_ann[match(rownames(edgeR_result_Exact$table),
                                                         rownames(countData_ann)),])
    table_counts_normalized = cbind(table_counts_normalized,
                                    countData_ann[match(rownames(table_counts_normalized),
                                                        rownames(countData_ann)),])
    write.table(edgeR_result_QLF$table,paste0(output_dir,outfilename1), sep="\t")
    write.table(edgeR_result_LRT$table,paste0(output_dir,outfilename2), sep="\t")
    write.table(edgeR_result_Exact$table,paste0(output_dir,outfilename3), sep="\t")
    write.table(table_counts_normalized,paste0(output_dir,normalizedfilename), sep="\t")
  }else{
    print(paste('Warning: File does not contain enough count data: ', infilename))
    }
}

########################################################
# perform getMM normalization

# for each organism, perform the same analysis
for (infilename in org_files){
  normalizedfilenameEdgeR <- paste("compare_methods_edgeR_nogroup_cpm_normalized_", infilename, sep="")
  normalizedfilenameGeTMM <- paste("compare_methods_getmm_normalized_", infilename, sep="")
  # read gene counts
  countData <- as.matrix(read.csv(paste0(org_dir,infilename)))#, row.names="gene_id"))
  rownames <- paste(countData[,colnames(countData)=="chr"], countData[,colnames(countData)=="gene_id"], sep="_")
  rownames(countData) <- rownames
  # replace NA with 0
  countData[is.na(countData)] <- 0
  # leave only data columns
  countData_counts <- countData[,grepl("MSZ", colnames(countData))]
  # convert to numeric
  storage.mode(countData_counts) <- "numeric"
  #################################################
  # remove samples where bactreia were not detected
  removeCols <- colSums(countData_counts) < mincount_threshold
  countData_counts <- countData_counts[,!removeCols]
  #################################################
  
  # get gene length
  #countData_genelength <- countData[,grepl("length", colnames(countData))]
  countData_genelength <- countData[,which(colnames(countData) == "length")]
  # convert to numeric
  countData_genelengthN <- as.numeric(gsub("\\[", "", gsub("\\]", "", countData_genelength)))
  # replace NA with 1
  countData_genelengthN[is.na(countData_genelengthN)] <- 1
  # calculate RPK
  rpk <- (countData_counts/countData_genelengthN)
  # for normalization purposes, no grouping of samples
  diet <- c(rep("A",ncol(countData_counts)))
  #################################################
  # remove samples where bactreia were not detected
  diet <- diet[!removeCols]
  #################################################
  if(length(diet) > 3){
    #EdgeR
    x.norm.edger <- DGEList(counts=countData_counts,group=diet)
    x.norm.edger <- calcNormFactors(x.norm.edger)
    norm.counts.edger <- cpm(x.norm.edger)
    #GeTMM
    rpk.norm <- DGEList(counts=rpk,group=diet)
    rpk.norm <- calcNormFactors(rpk.norm)
    norm.counts.rpk_edger <- cpm(rpk.norm)
    # add gene annotations from the original file to results
    countData_ann <- countData[,!grepl("MSZ", colnames(countData))]
    norm.counts.edger = cbind(norm.counts.edger,
                              countData_ann[match(rownames(norm.counts.edger),
                                                  rownames(countData_ann)),])
    norm.counts.rpk_edger = cbind(norm.counts.rpk_edger,
                                  countData_ann[match(rownames(norm.counts.rpk_edger),
                                                      rownames(countData_ann)),])
    write.table(norm.counts.edger,paste0(output_dir,normalizedfilenameEdgeR), sep="\t")
    write.table(norm.counts.rpk_edger,paste0(output_dir,normalizedfilenameGeTMM), sep="\t")
  }else{
    print(paste('Warning: File does not contain enough count data: ', infilename))
  }
}

########################################################
# perform getMM normalization on merged read matrices


org_files_subset = org_files[grepl("ucount", org_files)]
# remove merged count file from the file list (it does not include species name between RNA100 and _t_e)
org_files_subset = org_files_subset[!grepl("RNA100_t_e_", org_files_subset)]


rpk_all = matrix(, nrow = 0, ncol = 10)
rpk_ann = matrix(, nrow = 0, ncol = ncol(countData_ann))
colnames(rpk_ann) = colnames(countData_ann)
rpk_genes = c("")
for (infilename in org_files_subset){
  normalizedfilenameEdgeR <- paste("compare_methods_edgeR_nogroup_cpm_normalized_", infilename, sep="")
  normalizedfilenameGeTMM <- paste("compare_methods_getmm_normalized_", infilename, sep="")
  # read gene counts
  countData <- as.matrix(read.csv(paste0(org_dir,infilename), row.names="gene_id"))
  # replace NA with 0
  countData[is.na(countData)] <- 0
  # leave only data columns
  countData_counts <- countData[,grepl("MSZ", colnames(countData))]
  # convert to numeric
  storage.mode(countData_counts) <- "numeric"
  # get gene length
  #countData_genelength <- countData[,grepl("length", colnames(countData))]
  countData_genelength <- countData[,which(colnames(countData) == "length")]
  # convert to numeric
  countData_genelengthN <- as.numeric(gsub("\\[", "", gsub("\\]", "", countData_genelength)))
  # replace NA with 1
  countData_genelengthN[is.na(countData_genelengthN)] <- 1
  # calculate RPK
  rpk <- (countData_counts/countData_genelengthN)
  # add org to gene names
  curgenes <- paste(infilename, rownames(rpk), sep="_")
  rpk_all <- rbind(rpk_all,rpk)
  rpk_genes <- c(rpk_genes, curgenes)
  # add annotation matrix
  countData_ann <- countData[,!grepl("MSZ", colnames(countData))]
  rpk_columns_intersect <- intersect(colnames(rpk_ann), colnames(countData_ann))
  rpk_ann <- rbind(rpk_ann[,rpk_columns_intersect], countData_ann[, rpk_columns_intersect])
}
nrow(rpk_all)
nrow(rpk_ann)
if (rpk_genes[1]==""){
  rpk_genes <- rpk_genes[2:length(rpk_genes)]
}
# for normalization purposes, no grouping of samples
diet <- c(rep("A",ncol(rpk_all)))
x.norm.edger <- DGEList(counts=rpk_all,group=diet)
x.norm.edger <- calcNormFactors(x.norm.edger)
norm.counts.edger <- cpm(x.norm.edger)
# add gene names with organism info
rownames(norm.counts.edger) <- rpk_genes
norm.counts.edger = cbind(norm.counts.edger, rpk_ann)
#GeTMM
rpk.norm <- DGEList(counts=rpk_all,group=diet)
rpk.norm <- calcNormFactors(rpk.norm)
norm.counts.rpk_edger <- cpm(rpk.norm)
# add gene names with organism info
rownames(norm.counts.rpk_edger) <- rpk_genes
norm.counts.rpk_edger = cbind(norm.counts.rpk_edger, rpk_ann)
if(grepl('DNA',org_dir)){
  normalizedfilenameEdgeR <- "compare_methods_edgeR_nogroup_cpm_normalized_stringtie_DNA_t_e_extendedIDs_EGGNOGann.txt"
  normalizedfilenameGeTMM <- "compare_methods_getmm_normalized_stringtie_DNA_t_e_extendedIDs_EGGNOGann.txt"
}else{
  normalizedfilenameEdgeR <- "compare_methods_edgeR_nogroup_cpm_normalized_stringtie_RNA_t_e_extendedIDs_EGGNOGann.txt"
  normalizedfilenameGeTMM <- "compare_methods_getmm_normalized_stringtie_RNA_t_e_extendedIDs_EGGNOGann.txt"
}
write.table(norm.counts.edger,paste0(output_dir,normalizedfilenameEdgeR), sep="\t")
write.table(norm.counts.rpk_edger,paste0(output_dir,normalizedfilenameGeTMM), sep="\t")
# write r session info to file
writeLines(capture.output(sessionInfo()), "Scripts\\r_session_info.txt")