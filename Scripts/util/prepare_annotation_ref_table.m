% utility workflow to prepare metabolite annotation reference file
% get old reference table and remove AGORA annotation since they are not
% relevant anymore
annKEGGAGORA = readtable([rawdataFolder 'kegg_agora_09_2020_hmdb_06_2021_table.csv']);
% remove AGORA and leave only KEGG and HMDB
remove_agora = cellfun(@(x) (contains(x, 'cpd:C') | contains(x, 'HMDB')),...
    annKEGGAGORA.KEGGID);
annKEGGAGORA(remove_agora==0,:) = [];
% round mass
annKEGGAGORA.EXACT_MASS = round(annKEGGAGORA.EXACT_MASS*10000)/10000;
writetable(annKEGGAGORA, [rawdataFolder 'kegg_hmdb_06_2021_table.csv']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read MiMeDB database file downloaded 31/08/2025, release 2024-03-19
annMiMeDB = readtable([rawdataFolder 'mimedb_metabolites_v1_2024_03_19.csv']);
% leave only compound ID, name, formula, mol weight and mol mass columns to
% merge with KEGG and HMDB annotations
select_cols = {'mime_id', 'name', 'moldb_formula', 'moldb_mono_mass', 'moldb_average_mass'};
annMiMeDB_selected = annMiMeDB(:,select_cols);
% rename columns to match the ones from annKEGGAGORA
annMiMeDB_selected.Properties.VariableNames = annKEGGAGORA.Properties.VariableNames;
% merge tables
annKEGGAGORA_merged = [annKEGGAGORA; annMiMeDB_selected];
% remove compounds without exact mass
annKEGGAGORA_merged(isnan(annKEGGAGORA_merged.EXACT_MASS),:)=[];

% get unique formulas
[FORMULA_unique, idx_formula, idx_formula_all] = unique(annKEGGAGORA_merged.FORMULA);

% recalculate exact mass of formulas containing CHNOPS and no special
% characters

FORMULA_unique_recalculated_mass = calculateExactMass(FORMULA_unique, 0);

% add recalculates mass to the annotation table
idx_all_recalculated = cellfun(@(x) find(ismember(FORMULA_unique,x)),...
                                annKEGGAGORA_merged.FORMULA);
FORMULA_recalculated_mass_all = FORMULA_unique_recalculated_mass(idx_all_recalculated);
annKEGGAGORA_merged.EXACT_MASS_RECALCULATED = FORMULA_recalculated_mass_all;

% check non-matching masses
test = annKEGGAGORA_merged((abs(annKEGGAGORA_merged.EXACT_MASS_RECALCULATED-...
    annKEGGAGORA_merged.EXACT_MASS) > 0.0001) & (annKEGGAGORA_merged.EXACT_MASS_RECALCULATED>0),:);

% update mass to contain recalculated exact mass if it is nonzero and the
% original one otherwise
annKEGGAGORA_merged.EXACT_MASS_ORIGINAL = annKEGGAGORA_merged.EXACT_MASS;
annKEGGAGORA_merged.EXACT_MASS(annKEGGAGORA_merged.EXACT_MASS_RECALCULATED>0) = ...
    annKEGGAGORA_merged.EXACT_MASS_RECALCULATED(annKEGGAGORA_merged.EXACT_MASS_RECALCULATED>0);
% round mass
annKEGGAGORA_merged.EXACT_MASS = round(annKEGGAGORA_merged.EXACT_MASS*10000)/10000;

writetable(annKEGGAGORA_merged, [rawdataFolder 'kegg_hmdb_06_2021_mimedb_03_2024_table.csv']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make a table of Classyfire class and MimeDB compounds
mimedbclassification = annMiMeDB(:, {'classification', 'mime_id'});
writetable(mimedbclassification, [rawdataFolder 'mimedb_classification_03_2024_table.csv']);
% merge with HMDB classification
hmdbclassification = readtable([rawdataFolder 'HMDBclass072021_items.csv']);
% rename mimedbclassification columns and merge
mimedbclassification.Properties.VariableNames = hmdbclassification.Properties.VariableNames;
hmdb_mimedb_classification = [hmdbclassification; mimedbclassification];
writetable(hmdb_mimedb_classification, [rawdataFolder 'hmdb_07_2021_mimedb_03_2024_classification_table.csv']);

