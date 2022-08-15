# diet_microbiota_host_git_model
Pipelines of the Project 
*Multiomics and quantitative modelling disentangle diet, host, and microbiota contributions to the host metabolome*

By Zimmermann-Kogadeeva, Maria; Bencivenga-Barry, Natasha A; Zimmermann, Michael; Bork, Peer; Goodman, Andrew L.

*Contents:*

Folder _Scripts_ contains pipelines of the project: 

- add_global_and_file_dependencies.m 

Script defining global variables and file dependencies for other scripts. 
- workflow_fit_steadystate_fluxes_final.m

Pipeline that fits parameters of the coarse-grained intestinal flux model to metabolomics measurements along the GIT. Both forward problem (intestinal flux and metabolic flux parameters are estimated given normalized metabolite abundances) and reverse problem (metabolite abundances are estimated given intestinal flux and metabolic flux parameters estimated with the forward problem as coefficients) are solved. 
- calculateAmatrix_final.m 

Script that creates A matrix for the forward problem given metabolite measurements. 
- calculateRAmatrix_final.m 

Script that creates RA matrix and vector b for the reverse problem given intestinal flux and metabolic flux parameter estimates. 
- workflow_plot_and_analyse_modeling_results_final.m

Script that plots modeling results and performs hierarchical clustering of normalized mdel coefficients. 

Folder _ProcessedData_ should contain files required to run the pipelines. These files can be downloaded from Zenodo: https://doi.org/10.5281/zenodo.6993160 
