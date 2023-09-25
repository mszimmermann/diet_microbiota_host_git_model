function [met_info, gitfits] = read_allsols_from_files(filename)
% met_info contains information about metabolites
% gitfits contains all solutions
% filename is the input file name without extension
% WARNING! column from which coefs starts is hardcoded
datastartidx = 7;
% WARNING! number of solutions is hardcoded
nsols = 100;

curfilename = [filename...
               '.csv'];
soldata = readtable(curfilename,...
    'ReadVariableNames',1, 'HeaderLines',0);
% get metabolite info
met_info.MZ = soldata.MZ;
met_info.RT = soldata.RT;
met_info.CompoundID = soldata.CompoundID;
met_info.CompoundName = soldata.CompoundName;
met_info.MetaboliteFilter = soldata.MetaboliteFilter;
met_info.SumGITclusters = soldata.SumGITclusters;

% get solution correlation info and solution as matrix
columnNames = soldata.Properties.VariableNames;
columnNames_coefs = columnNames(datastartidx:end);
% save the filed names in the same way as in bestsol
testx = soldata{:,columnNames_coefs};

coefvalues = reshape(columnNames_coefs,[],nsols);
coefvalues = coefvalues(:,1)';
coefvalues = cellfun(@(x) x(1:strfind(x, '_')-1), coefvalues, 'unif',0);

%create array of gitfits
gitfits = cell(size(testx,1),1);
for i=1:size(testx,1)
    gitfit = struct; 
    gitfit.coefvalues = coefvalues;
    gitfit.testx = reshape(testx(i,:), [], 100);
    gitfits{i} = gitfit;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% read restored data from reciprocal file name
curfilenameR = [filename...
               '_reciprocal.csv'];
soldataR = readtable(curfilenameR,...
    'ReadVariableNames',1, 'HeaderLines',0);

% check that the metabolites are the same in soldata and soldataR
if ~isequal(soldata.CompoundID, soldataR.CompoundID)
    fprintf(['Warning: metabolite IDs do not match between:\n%s\nand\n%s\n'...
             'omit loading reciprocal data\n'], curfilename, curfilenameR)
else
    % load original data
    % get columns with original data
    columnNamesR = soldataR.Properties.VariableNames;
    columnNamesR_reciprocal_data_idx = find(cellfun(@(x) contains(lower(x), ...
                    'recip_'), columnNamesR));
    columnNamesR_data = columnNamesR(datastartidx:columnNamesR_reciprocal_data_idx(1)-1);
    kmeanMatrix_joint_orig = soldataR{:,datastartidx:columnNamesR_reciprocal_data_idx(1)-1};
    kmeanMatrix_joint_names = reshape(columnNamesR_data, [],6);
    dataR = soldataR{:,columnNamesR_reciprocal_data_idx};
    for i=1:size(dataR,1)
        cur_kmeanMatrix_joint_orig = reshape(kmeanMatrix_joint_orig(i,:),[],6);
        cur_dataR = reshape(dataR(i,:), nsols,[]);
        gitfits{i}.kmeanMatrix_joint_orig = cur_kmeanMatrix_joint_orig;
        gitfits{i}.kmeanMatrix_joint_names = kmeanMatrix_joint_names;
        gitfits{i}.testdataR = cur_dataR;
        %calculate correlations between original and restored data
        testCorrRev = zeros(1,nsols);
        testCorrRevSI = zeros(1,nsols);
        testCorrRevLI = zeros(1,nsols);
        for k=1:nsols
            cur_dataRi = reshape(cur_dataR(k,:),[],6);
            testCorrRev(k) = corr(kmeanMatrix_joint_orig(k,:)',...
                                      cur_dataR(k,:)');
            testCorrRevSI(k) = corr(reshape(cur_kmeanMatrix_joint_orig(:,1:3),[],1),...
                                    reshape(cur_dataRi(:,1:3),[],1));
            testCorrRevLI(k) = corr(reshape(cur_kmeanMatrix_joint_orig(:,4:6),[],1),...
                                    reshape(cur_dataRi(:,4:6),[],1)); 
        end
        gitfits{i}.testCorrRev = testCorrRev;
        gitfits{i}.testCorrRevSI = testCorrRevSI;
        gitfits{i}.testCorrRevLI = testCorrRevLI;
    end
end
        
        
        
        
        