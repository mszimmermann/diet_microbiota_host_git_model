function print_allsols_to_files(met_info, gitfits, filename)
% met_info contains information about metabolites
% gitfits contains all solutions
% filename is the output file name

% check whether met_info contains MZ and RT fiels
% if not (e.g. drug info) set them to 0
if isstruct(met_info)
    if isfield(met_info, 'MZ') == 0
        met_info.MZ = zeros(size(gitfits,1),1);
    end
    if isfield(met_info, 'RT') == 0
        met_info.RT = zeros(size(gitfits,1),1);
    end
    % check whether met_info contains MetaboliteFilter - flag indicating
    % annotation filtering
    if isfield(met_info, 'MetaboliteFilter') == 0
        met_info.MetaboliteFilter = ones(size(gitfits,1),1);
    end
    
    % check whether met_info contains gut_filter - flag indicating
    % whether metabolites were detected in the GIT
    if isfield(met_info, 'gut_filter') == 0
        met_info.gut_filter = ones(size(gitfits,1),1);
    end
elseif istable(met_info)
    if ~ismember('MZ', met_info.Properties.VariableNames)
        met_info.MZ = zeros(size(gitfits,1),1);
    end
    if ~ismember('RT', met_info.Properties.VariableNames)
        met_info.RT = zeros(size(gitfits,1),1);
    end
     % check whether met_info contains MetaboliteFilter - flag indicating
    % annotation filtering
    if ~ismember('MetaboliteFilter', met_info.Properties.VariableNames)
        met_info.MetaboliteFilter = ones(size(gitfits,1),1);
    end
    
    % check whether met_info contains gut_filter - flag indicating
    % whether metabolites were detected in the GIT
    if ~ismember('gut_filter', met_info.Properties.VariableNames)
        met_info.gut_filter = ones(size(gitfits,1),1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save model results to file - raw coefficients

% save each type of selected solution to a separate file
curfilename = [filename...
               '.csv'];
fid = fopen(curfilename, 'w');
fprintf(fid,'MZ\tRT\tCompoundID\tCompoundName\tMetaboliteFilter\tSumGITclusters');
% create column names
columnnames = strcat(repmat(gitfits{1}.coefvalues, 1, size(gitfits{1}.testx,2)),...
    '_',...
    arrayfun(@(x) sprintf('%03d',x), reshape(repmat((1:size(gitfits{1}.testx,2))',...
                                             1,length(gitfits{1}.coefvalues))',...
                                      1,[]),...
                              'unif',0));
                                 
for i=1:length(columnnames)
    fprintf(fid, '\t%s', columnnames{i});
end
fprintf(fid, '\n');
for i=1:length(gitfits)
    fprintf(fid, '%.3f\t%3f',   met_info.MZ(i),...
                                met_info.RT(i));
    fprintf(fid, '\t%s\t%s\t%d',met_info.CompoundID{i},...
                                met_info.CompoundName{i},...
                                met_info.MetaboliteFilter(i));
    fprintf(fid, '\t%d', met_info.gut_filter(i));


    if ~isempty(gitfits{i})
        testx = gitfits{i}.testx;
        testx = reshape(testx, [],1);
        for j=1:length(testx)
            fprintf(fid, '\t%e', testx(j));
        end
    else
        for j=1:numel(gitfits{1}.testx)
            fprintf(fid, '\t0');
        end
    end
        
    fprintf(fid, '\n');
end
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% save model results to file - reciprocal data restoration
% save each type of selected solution to a separate file
curfilename = [filename...
               '_reciprocal.csv'];
fid = fopen(curfilename, 'w');

% column names for the data
columnNames = gitfits{1}.kmeanMatrix_joint_names(:);

fprintf(fid, 'MZ\tRT\tCompoundID\tCompoundName\tMetaboliteFilter\tSumGITclusters');
for i=1:length(columnNames)
    fprintf(fid, '\t%s', columnNames{i});
end
for k=1:size(gitfits{1}.testdataR,2)
    for i=1:length(columnNames)
        fprintf(fid, '\tRecip_%s_%03d', columnNames{i},k);
    end
end
fprintf(fid, '\n');
for i=1:length(gitfits)
    fprintf(fid, '%.3f\t%3f',   met_info.MZ(i),...
                                met_info.RT(i));
    fprintf(fid, '\t%s\t%s\t%d',met_info.CompoundID{i},...
                                met_info.CompoundName{i},...
                                met_info.MetaboliteFilter(i));
    fprintf(fid, '\t%d', met_info.gut_filter(i));

    if ~isempty(gitfits{i})
        kmean_vector_joint_orig = gitfits{i}.kmeanMatrix_joint_orig(:);
        for j=1:length(kmean_vector_joint_orig)
            fprintf(fid, '\t%.3f', kmean_vector_joint_orig(j));
        end
        dataR = reshape(gitfits{i}.testdataR,[],1); 
        for j=1:length(dataR)
            fprintf(fid, '\t%.3f', dataR(j));
        end
    else
        %"original data"
        for j=1:numel(gitfits{1}.kmeanMatrix_joint_orig)
            fprintf(fid, '\t0');
        end
        % "reciprocal data"
        for j=1:numel(gitfits{1}.testdataR)
            fprintf(fid, '\t0');
        end
    end
    fprintf(fid, '\n');
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save model results to file - random coefficients (for comparison)

% save each type of selected solution to a separate file
curfilename = [filename...
               '_random_coef_correv.csv'];
fid = fopen(curfilename, 'w');
fprintf(fid,'MZ\tRT\tCompoundID\tCompoundName\tMetaboliteFilter\tSumGITclusters');
% create column names
columnnames = arrayfun(@(x) sprintf('randomCorr_%03d',x), (1:length(gitfits{1}.testCorrRev_shuffled)), 'unif', 0);
                                 
for i=1:length(columnnames)
    fprintf(fid, '\t%s', columnnames{i});
end
fprintf(fid, '\n');
for i=1:length(gitfits)
    fprintf(fid, '%.3f\t%3f',   met_info.MZ(i),...
                                met_info.RT(i));
    fprintf(fid, '\t%s\t%s\t%d',met_info.CompoundID{i},...
                                met_info.CompoundName{i},...
                                met_info.MetaboliteFilter(i));
    fprintf(fid, '\t%d', met_info.gut_filter(i));


    if ~isempty(gitfits{i})
        testx = gitfits{i}.testCorrRev_shuffled;
        for j=1:length(testx)
            fprintf(fid, '\t%e', testx(j));
        end
    else
        for j=1:numel(gitfits{1}.testCorrRev_shuffled)
            fprintf(fid, '\t0');
        end
    end
        
    fprintf(fid, '\n');
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% save random coefficients as well
curfilename = [filename...
               '_random_coefs.csv'];
fid = fopen(curfilename, 'w');
fprintf(fid,'MZ\tRT\tCompoundID\tCompoundName\tMetaboliteFilter\tSumGITclusters');
% create column names
if length(gitfits{1}.testx_shuffled) > length(gitfits{1}.coefvalues)
    % coef names should be doubled because there was only one diet
    % condition
    coefvalues_expanded = {gitfits{1}.coefvalues{1},...
                           gitfits{1}.coefvalues{2}, [gitfits{1}.coefvalues{2} '2'],...
                           gitfits{1}.coefvalues{3}, [gitfits{1}.coefvalues{3} '2'],...
                           gitfits{1}.coefvalues{4}, [gitfits{1}.coefvalues{4} '2'],...
                           gitfits{1}.coefvalues{5}, [gitfits{1}.coefvalues{5} '2']};
    columnnames = strcat(repmat(coefvalues_expanded, 1, size(gitfits{1}.testx_shuffled,2)),...
        '_',...
        arrayfun(@(x) sprintf('%03d',x), reshape(repmat((1:size(gitfits{1}.testx_shuffled,2))',...
                                             1,length(coefvalues_expanded))',...
                                      1,[]),...
                              'unif',0));                       
else
    columnnames = strcat(repmat(gitfits{1}.coefvalues, 1, size(gitfits{1}.testx_shuffled,2)),...
        '_',...
        arrayfun(@(x) sprintf('%03d',x), reshape(repmat((1:size(gitfits{1}.testx_shuffled,2))',...
                                             1,length(gitfits{1}.coefvalues))',...
                                      1,[]),...
                              'unif',0));
end                             
for i=1:length(columnnames)
    fprintf(fid, '\t%s', columnnames{i});
end
fprintf(fid, '\n');
for i=1:length(gitfits)
    fprintf(fid, '%.3f\t%3f',   met_info.MZ(i),...
                                met_info.RT(i));
    fprintf(fid, '\t%s\t%s\t%d',met_info.CompoundID{i},...
                                met_info.CompoundName{i},...
                                met_info.MetaboliteFilter(i));
    fprintf(fid, '\t%d', met_info.gut_filter(i));


    if ~isempty(gitfits{i})
        testx = gitfits{i}.testx_shuffled;
        testx = reshape(testx, [],1);
        for j=1:length(testx)
            fprintf(fid, '\t%e', testx(j));
        end
    else
        for j=1:numel(gitfits{1}.testx)
            fprintf(fid, '\t0');
        end
    end
        
    fprintf(fid, '\n');
end
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
