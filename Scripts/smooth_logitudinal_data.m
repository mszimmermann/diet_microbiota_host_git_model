function [smoothloc, smoothcond, smoothrep, smoothdata] = ...
    smooth_logitudinal_data(curloc,curcond, currep, curdata)
% take average between locations to smooth metabolite profiles
smooth_sections = { [1 2 3 4 5],...
                    [6 7 8],...
                    [9 10 11],...
                    [12],...
                    [13 14],...
                    [15]};
smoothloc_unique = 1:length(smooth_sections);
curcond_unique = unique(curcond);
currep_unique = unique(currep);

numsamples = length(smoothloc_unique)*...
             length(curcond_unique)*...
             length(currep_unique);

smoothloc = zeros(numsamples,1);
smoothcond = cell(numsamples,1);
smoothrep = cell(numsamples,1);
smoothdata = zeros(numsamples,1);

idx = 1;
for sec_i=1:length(smooth_sections)
    cursec = smooth_sections{sec_i};
    for cond_i = 1:length(curcond_unique)
        for rep_i = 1:length(currep_unique)
            cursmoothdata = curdata(ismember(curloc,cursec)&...
                ismember(curcond, curcond_unique(cond_i))&...
                ismember(currep, currep_unique(rep_i)));
            % take sume of intensities in several sections of the same SI
            smoothdata(idx) = nansum(cursmoothdata);%nanmean(cursmoothdata);
            smoothloc(idx) = sec_i;
            smoothcond{idx} = curcond_unique{cond_i};
            smoothrep{idx} = currep_unique{rep_i};
            idx = idx+1;
        end
    end
end
            
            
        


