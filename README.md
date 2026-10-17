# diet\_microbiota\_host\_git\_model

Pipelines of the Project
*Multiomics and quantitative modelling disentangle diet, host, and microbiota contributions to the host metabolome*

By Zimmermann-Kogadeeva, Maria; Henseler, Zachariah M.; Bencivenga-Barry, Natasha A; Gross, Matthias W.; Zimmermann, Michael; Bartmanski, Bartosz J.; Bork, Peer; Goodman, Andrew L.

*Contents:*

1. *Coarse-grained modelling workflow: a coarse-grained intestinal flux model that describes metabolite propagartion along the gastrointestinal tract (GIT). In each GIT compartment, the changes in metabolite levels are described with an equation including compartment-specific consuming and producing metabolic fluxes. The model connects GIT compartments via intestinal flux which represents metabolite transit. In each compartment, the change in metabolite abundance due to intestinal transit follows the law of mass action with a condition- and compartment-independent parameter f. In the small intestine, parameter f\_SI\_host describes host metabolism flux, and in the large intestine host metabolism fluxes are described by parameters f\_Cecum\_host and f\_Colon\_host for cecum and colon, respectively. In the cecum of colonized mice, metabolite levels are also influenced by microbiota metabolism flux parameter f\_LI\_microbe. The model assumption is that host metabolism parameters are independent of mouse colonization state, while all host and microbiota metabolism parameters are diet-specific. The resulting model contains 20 equations (for jejunum, ileum, cecum, colon and feces for each of the four conditions) with nine parameters (the nonnegative intestinal flux parameter f and eight diet-specific host and microbiota metabolism parameters which can have either negative or positive values corresponding to metabolite consumption or production).* 
2. *Scripts to reproduce the analysis and figures of the project.* 



*Software requirements: The code is developed with MatLab2023a Optimization toolbox and requires MatLab version 2023a, although most scripts should be compatible with MatLab version 2019.*



*License: The code is distributed under the terms of the GNU General Public License (please read copyright\_and\_license and LICENSE files for details).*



*Installation: There is no specific installation required to run the code.* 



Folder *Scripts* contains pipelines of the project:

* add\_global\_and\_file\_dependencies.m

Script defining global variables and file dependencies for other scripts.



Transcriptomics analysis



1\) workflow\_preprocess\_RNADNA\_stringtie\_data\_eggnog\_final.m

preprocess tables from RNA/DNA and merge

get counts matrices

get KEGG/COG/GO terms from the annotations

get full and filtered gene annotation tables



2\) workflow\_RNA\_ptw\_enrichment\_final.m

Perform pathway enrichment on the merged tables



3\) workflow\_merge\_rnaseq\_tables.m

Merge read files and pathway enrichment files created in 1 and 2 and make tables S4 and S5



4\) workflow\_PCA\_analysis\_metatranscriptomics\_final.m

Makes PCA plots for supplementary figures and count matrix files with species names that are needed for subsequent analysis



5\) workflow\_analysize\_metaGmetaTcomposition\_final.m

Gets counts tables as input (with species, generates in the script in point 4), calculates metaG and metaT species abundances, compares directly species compositions across metaG, metaT and metaphlan OTUs, makes plots for RNAseq species composition and RNA/DNA ratios. 



Metabolomics analysis



1\) workflow\_merge\_untargeted\_diet\_data\_final.m

Takes raw metabolomics files and tissue weights and produces quantile-normalizes intensity table and full (unfiltered) annotation table

% 'metabolites\_allions\_combined\_norm\_intensity.csv'

% 'metabolites\_allions\_combined\_norm\_intensity\_with\_CVR.csv'

% 'metabolites\_allions\_combined\_formulas.csv'



2\) workflow\_filter\_metabolite\_ions.m

Takes annotation table and normalized ions and filtered based on a set of criteria. Filtered ion list as output: 'metabolites\_allions\_combined\_formulas\_with\_metabolite\_filters\_updated\_filtering\_0925.csv'

Also combined intensities with annotation filtering table to make table S6: 'table\_s6\_metabolites\_allions\_norm\_intensity\_with\_filters\_updated\_0925.csv'



3\) workflow\_spatial\_clustering\_final.m

Performs spatial clustering and HMDB group enrichment. Independent of others, requires annotation table and normalized intensities. 



4\) workflow\_metabolomics\_diff\_analysis.m

Performs differential analysis per tissue between groups



5\) workflow\_PCA\_analysis\_metabolomics\_final.m

Plots PCA analysis results and statistics on mouse group ion overlap and tissue overlap for supplementary figures. Requires annotated and filtered metabolite tables. 



6\) workflow\_plot\_metabolite\_profiles.m

Plots profiles of selected metabolites with standard deviations and data points (independent of modelling)



Modelling scripts



7\) workflow\_fit\_steadystate\_fluxes\_to\_drugs\_final.m

prepare drug data for modelling and perform modelling and model type selection

Ideally should be run before all metabolite model fitting as this is the model selection procedure, but in principle independent

* calculateAmatrix\_final.m

Script that creates A matrix for the forward problem given metabolite measurements.

* calculateRAmatrix\_final.m,

Script that creates RA matrix and vector b for the reverse problem given intestinal flux and metabolic flux parameter estimates.



8\) workflow\_fit\_steadystate\_fluxes\_final.m

Fit the model to the metabolites (DC and CVR separately) and save all solutions and best solutions to file. Note for all ions (\~57000) model run takes \~2.5 hours. 

Pipeline that fits parameters of the coarse-grained intestinal flux model to metabolomics measurements along the GIT. Both forward problem (intestinal flux and metabolic flux parameters are estimated given normalized metabolite abundances) and reverse problem (metabolite abundances are estimated given intestinal flux and metabolic flux parameters estimated with the forward problem as coefficients) are solved.



9\) workflow\_plot\_and\_analyse\_modeling\_results\_final.m

Reads combined solution from file and makes a clustergram. Clustergram is manually split into subclusters and pathway enrichment per cluster is performed. Output tables for model solutions are produced and correlation plots for reciprocal solution. 



10\) plot\_modeling\_profiles\_and\_coefficients.m

Plots model coefficients and original and reciprocal restored data per metabolite



11\) workflow\_interpretation\_of\_model\_coefficients.m

Read model coefficients and metabolite fold changes and compare them to identify how to best interpret model coefficients. 



12\) workflow\_analyze\_meier\_data.m

Prepares Meier et al 2023 data for modelling and calculates coefficients



13\) workflow\_model\_pairwise\_comparison.m

Compares modelling results between CVR and DC groups and between DC and CVR group and Meier et al modelling results

Note: file workflow\_compare\_with\_meier\_study.m is outdated, all comparison is done within the workflow\_model\_pairwise\_comparison.m script. 



14\) workflow\_extract\_keggECpathes\_for\_SPpairs\_final.m

Extracts enzyme pairs from KEGG pathways.



15\) workflow\_get\_ec\_pathes\_short.m

Filters enzyme paths for paths of short length. 



16\) workflow\_combine\_paths\_with\_rnaseq\_final.m

Combined paths analysis between substrates and products with information of gene expression from metatranscriptomics measurements. 



Folders *InputData* and *ProcessedData* should contain files required to run the pipelines. These files can be downloaded from Zenodo: https://doi.org/10.5281/zenodo.6993160



Expected output: Expected output is provided in the ProcessedData/output folders (model parameter estimates, normalized and analyzed metabolomics and metatranscriptomics datasets).





