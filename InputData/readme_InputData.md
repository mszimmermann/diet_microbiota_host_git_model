#######################################################
# Folder InputData
- mouse_info.csv - information on Mouse_number,Gender,Colonization state,Diet,and Weight at day 2 (end of experiment)
- HMDBclass072021_items.csv - list of HMDB classes and metabolite IDs belonging to the classess (accessed 07 2021)                  
- HMDBsubclass072021_items.csv - list of HMDB subclasses and metabolite IDs belonging to the classess (accessed 07 2021)               
- HMDBsuperclass072021_items.csv - list of HMDB superclasses and metabolite IDs belonging to the classess (accessed 07 2021)             
- kegg_agora_09_2020_hmdb_06_2021_table.csv - list of metabolites used for ion annotation 
- kegg_ptw_names_and_bact_flag.csv - names of KEGG pathways and flag of being present in bacteria           

## Folder KEGGreaction_path
Contains matlab files with metabolite-metabolite paths calculated from KEGG reaction-pair information (Each matrix contains a subset of paths). These files are used by the script workflow_extract_keggECpathes_for_SPpairs_final.m.

## Folder metabolomics_data
Contains raw metabolomics data from six methods (three LC columns: C08, C18 and HILIC, and positive and negative acquisition modes) and file tissue_weights.txt with tissue weight information used for normalization. 

## Folder sequencing_data
Contains folders ballgown_DNA and ballgown_RNA with results of metagenomic and metatranscriptomic data analysis (raw counts, GetMM normalized counts, EdgeR and DeSeq2 analysis).   

## Folder genomes14_eggnog
Contains folders eggNOG annotations of the 14 bacterial genomes.   

## Folder human_microbiome_data
Should contain files hub.cellcount.motu.Phylum.v2.data.frame.R, hub.cellcount.motu.Species.v2.data.frame.R and hub.adjusted_KEGG_module.down.10000000.v3.data.frame.R from the depository https://zenodo.org/records/6242715 file input_features.tar.gz

