function [ann, annFlag] = annotateIonsByMZ(MzRT, cpdMASS, cpdID, cpdNAME, annThres, varargin)
% annotate ions in MzRT based on theit m/z ratio given a list of 
% cpdMASS - compound masses
% cpdID - compound ID to save in annotation
% cpdNAME - compound name to save in annotation
% annThres - annotation threshold
% varargin - additional argument containing a list of reference Cpd RT to 
% use for annotation
% RETURN annotation structure ann and annotation flag annFlag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ann = struct;
annFlag = zeros(length(MzRT),1);
refRT = {'',0,0};
refRTthreshold = 1;
cpdFORMULA = [];
if nargin >= 6
    if ~isempty(varargin{1})
        refRT = varargin{1};
    end
    if nargin >=7
        if ~isempty(varargin{2})
            refRTthreshold = varargin{2};
        end
    end
    if nargin >=8
        cpdFORMULA = varargin{3};
    end
end
for i=1:length(MzRT)
    ann(i).annID = [];
    ann(i).annNames = [];
    ann(i).annMzdelta = [];
    ann(i).annIDX = [];
    if ~isempty(cpdFORMULA)
        ann(i).annFormula = [];
    end 
    ann(i).byReference = 0;
    annFlag(i) = 0;
                    
    closestMatch = find(abs(cpdMASS-MzRT(i,1)) <= annThres);
    if ~isempty(closestMatch)
        for j=1:length(closestMatch)
            cur_refRT = refRT(ismember(refRT(:,1),cpdID(closestMatch(j))),3);
            if ~isempty(cur_refRT)
                if abs(MzRT(i,2)-cell2mat(cur_refRT)) <= refRTthreshold
                    ann(i).annID = [ann(i).annID;cpdID(closestMatch(j))];
                    ann(i).annNames = [ann(i).annNames; cpdNAME(closestMatch(j))];
                    ann(i).annMzdelta = [ann(i).annMzdelta cpdMASS(closestMatch(j))-MzRT(i,1)];
                    ann(i).annIDX = [ann(i).annIDX; closestMatch(j)];
                    ann(i).byReference = 1;
                    annFlag(i) = 1;
                    if ~isempty(cpdFORMULA)
                         ann(i).annFormula = [ann(i).annFormula cpdFORMULA(closestMatch(j))];
                    end 
                end
            else
                ann(i).annID = [ann(i).annID;cpdID(closestMatch(j))];
                ann(i).annNames = [ann(i).annNames; cpdNAME(closestMatch(j))];
                ann(i).annMzdelta = [ann(i).annMzdelta cpdMASS(closestMatch(j))-MzRT(i,1)];
                ann(i).annIDX = [ann(i).annIDX; closestMatch(j)];
                annFlag(i) = 1;
                if ~isempty(cpdFORMULA)
                    ann(i).annFormula = [ann(i).annFormula cpdFORMULA(closestMatch(j))];
                end 
            end
        end
    end
end

% delete alternative annotations from ions annotated by reference
references = find([ann.byReference] == 1);
for i=1:length(references)
    [~, refidx] = intersect(ann(references(i)).annID, refRT(:,1));
    ann(references(i)).annID = ann(references(i)).annID(refidx);
    ann(references(i)).annNames = ann(references(i)).annNames(refidx);
    ann(references(i)).annMzdelta = ann(references(i)).annMzdelta(refidx);
    ann(references(i)).annIDX = ann(references(i)).annIDX(refidx);
    if ~isempty(cpdFORMULA)
         ann(references(i)).annFormula = ann(references(i)).annFormula(refidx);
    end 
end

    
    