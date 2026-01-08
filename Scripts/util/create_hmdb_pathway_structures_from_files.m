%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create matlab object with HMDB pathways from tables prepared in Python

% call script defining file dependenciesand global variables
addpath(genpath('.\'));
add_global_and_file_dependencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data requirements: 
% 'HMDBclass072021_items.csv'
% 'HMDBsubclass072021_items.csv'
% 'HMDBsuperclass072021_items.csv'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output:
% Files:
% 'hmdbV4_072021_ptw.mat'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%hmdbClassFile = [rawdataFolder 'HMDBclass072021_items.csv'];
hmdbClassFile = [rawdataFolder 'hmdb_07_2021_mimedb_03_2024_classification_table.csv'];%'HMDBclass072021_items.csv'];
hmdbSubClassFile = [rawdataFolder 'HMDBsubclass072021_items.csv'];
hmdbSuperclassFile = [rawdataFolder 'HMDBsuperclass072021_items.csv'];

hmdbClass = readtable(hmdbClassFile, 'delim', ',');
hmdbSubClass = readtable(hmdbSubClassFile, 'delim', ',');
hmdbSuperclass = readtable(hmdbSuperclassFile, 'delim', ',');

hmdbClass_unique = unique(hmdbClass.HMDBclass);
hmdbSubClass_unique = unique(hmdbSubClass.HMDBsubclass);
hmdbSuperclass_unique = unique(hmdbSuperclass.HMDBsuperclass);

for i=1:length(hmdbClass_unique)
    hmdbClass_ptw(i).ptw.name = hmdbClass_unique{i};
    hmdbClass_ptw(i).ptw.cmpdID = hmdbClass.HMDBID(ismember(hmdbClass.HMDBclass,...
                                                            hmdbClass_unique{i}));
end
                                                        
for i=1:length(hmdbSubClass_unique)
    hmdbSubClass_ptw(i).ptw.name = hmdbSubClass_unique{i};
    hmdbSubClass_ptw(i).ptw.cmpdID = hmdbSubClass.HMDBID(ismember(hmdbSubClass.HMDBsubclass,...
                                                            hmdbSubClass_unique{i}));
end
         
for i=1:length(hmdbSuperclass_unique)
    hmdbSuperClass_ptw(i).ptw.name = hmdbSuperclass_unique{i};
    hmdbSuperClass_ptw(i).ptw.cmpdID = hmdbSuperclass.HMDBID(ismember(hmdbSuperclass.HMDBsuperclass,...
                                                            hmdbSuperclass_unique{i}));
end

%save([processeddataFolder filesep 'util' filesep 'hmdbV4_072021_ptw.mat'],...
%        'hmdbClass_ptw', 'hmdbSuperClass_ptw', 'hmdbSubClass_ptw')
save([processeddataFolder filesep 'util' filesep 'hmdbV4_072021_mimedb_0324_ptw.mat'],...
        'hmdbClass_ptw', 'hmdbSuperClass_ptw', 'hmdbSubClass_ptw')
