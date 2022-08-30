%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create matlab object with HMDB pathways from tables prepared in Python

% call script defining file dependenciesand global variables
addpath(genpath('.\'));
add_global_and_file_dependencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data requirements: 
% 'hmdbV4_072021_ptw.mat'
% 'metabolites_allions_combined_formulas_with_metabolite_filters_spatial100clusters.csv'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output:
% Files:
% 'hmdbPTWtables.mat'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

annotationTable = readtable([inputFolder...
'metabolites_allions_combined_formulas_with_metabolite_filters_spatial100clusters.csv']);

hmdbPTW = load([processeddataFolder filesep 'util' filesep ...
        'hmdbV4_072021_ptw.mat']);

%prepare hmdb ptw-met tables
keggCompoundsID = annotationTable.CompoundID;
for hmdbtype = 1:3
    switch hmdbtype
        case 1
            hmdbFG = hmdbPTW.hmdbClass_ptw;
            [hmdbPTWclassTable, hmdbPTWclassNames] = ...
                    create_pathway_ion_table(hmdbFG, keggCompoundsID);
        case 2
            hmdbFG = hmdbPTW.hmdbSubClass_ptw;
            [hmdbPTWSubClassTable, hmdbPTWSubClassNames] = ...
                    create_pathway_ion_table(hmdbFG, keggCompoundsID);

        case 3
            hmdbFG = hmdbPTW.hmdbSuperClass_ptw;
            [hmdbPTWSuperClassTable, hmdbPTWSuperClassNames] = ...
                    create_pathway_ion_table(hmdbFG, keggCompoundsID);

    end
end

save([outputFolder 'hmdbPTWtables.mat'], 'hmdbPTWclassNames', 'hmdbPTWclassTable',...
    'hmdbPTWSubClassNames', 'hmdbPTWSubClassTable',...
    'hmdbPTWSuperClassNames', 'hmdbPTWSuperClassTable')
