\#######################################################

# Folder InputData

Is available at Zenodo: https://zenodo.org/records/7040318

* mouse\_info.csv - information on Mouse\_number,Gender,Colonization state,Diet,and Weight at day 2 (end of experiment)
* mouse\_fecal\_weights.csv - information on the weights of fecal pellets collected at the end of the experiment
* mouse\_CFU\_DC.csv - colony forming units in fecal samples from mice colonized with the defined community
* mouse\_CFU\_CVR.csv - colony forming units in fecal samples from conventional mice  
* HMDBclass072021\_items.csv - list of HMDB classes and metabolite IDs belonging to the classess (accessed 07 2021)
* HMDBsubclass072021\_items.csv - list of HMDB subclasses and metabolite IDs belonging to the classess (accessed 07 2021)
* HMDBsuperclass072021\_items.csv - list of HMDB superclasses and metabolite IDs belonging to the classess (accessed 07 2021)
* kegg\_agora\_09\_2020\_hmdb\_06\_2021\_table.csv - list of metabolites used for ion annotation
* kegg\_ptw\_names\_and\_bact\_flag.csv - names of KEGG pathways and flag of being present in bacteria
* kegg\_hmdb\_06\_2021\_mimedb\_03\_2024\_table.csv - updated list of metabolites used for ion annotation which includes MiMeDB (microbial metabolite database) and does not include metabolite IDs from AGORA models
* hmdb\_07\_2021\_mimedb\_03\_2024\_classification\_table.csv - joint table of HMDB classes and metabolite IDs from HMDB and MiMeDB daabases (accessed 07 2021 and 08 2025 correspondingly)



## Folder KEGGreaction\_path

Contains matlab files with metabolite-metabolite paths calculated from KEGG reaction-pair information (Each matrix contains a subset of paths). These files are used by the script workflow\_extract\_keggECpathes\_for\_SPpairs\_final.m.

## Folder metabolomics\_data

Contains raw metabolomics data from six methods (three LC columns: C08, C18 and HILIC, and positive and negative acquisition modes) and file tissue\_weights.txt with tissue weight information used for normalization.

## Folder sequencing\_data

Contains folders ballgown\_DNA and ballgown\_RNA with results of metagenomic and metatranscriptomic data analysis (raw counts, GetMM normalized counts, EdgeR and DeSeq2 analysis).

## Folder genomes14\_eggnog

Contains folders eggNOG annotations of the 14 bacterial genomes.

## Folder human\_microbiome\_data

Should contain files hub.cellcount.motu.Phylum.v2.data.frame.R, hub.cellcount.motu.Species.v2.data.frame.R and hub.adjusted\_KEGG\_module.down.10000000.v3.data.frame.R from the depository https://zenodo.org/records/6242715 file input\_features.tar.gz

## Folder public\_data

Should contain two files:

* 2023\_Meier\_NatMet\_42255\_2023\_802\_MOESM3\_ESM.xlsx which is Supplementary Data file (Processed and quantified metabolite data per sample) from Meier et al 2023 (https://doi.org/10.1038/s42255-023-00802-1)
* aat9931\_tables\_s1-26.xlsx which is Supplementary Tables file from Zimmermann et al. 2019 (https://doi.org/10.1126/science.aat9931)



