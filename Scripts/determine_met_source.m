function model_source = determine_met_source(bestsol)

% normalize just in case it was not 
bestsol = bestsol/max(abs(bestsol));
[~, maxidx] = max(abs(bestsol));
maxval = bestsol(maxidx);
model_source = 0;

% determine the source
if length(bestsol)==4
    switch maxidx
        case 1
            model_source = 0; % diet
        case {2,3}
            model_source = 1; % host production
        case 4
            model_source = 2; % bacterial production
    end
else
    if length(bestsol)==8
        switch maxidx
        case {1,4}
            model_source = 0; % diet
        case {2,3,5,6}
            model_source = 1; % host production
        case {7,8}
            model_source = 2; % bacterial production
        end
    else
        fprintf(['Warning: solution in the wrong format, expect length either 4 or 8\n'...
                 'Provided solution has length %d\n'], length(bestsol));
    end
end
model_source = model_source*maxval; % determine positive or negative
