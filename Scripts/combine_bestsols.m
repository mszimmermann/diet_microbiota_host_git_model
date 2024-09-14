function [bestsol, picksol] = combine_bestsols(cursol, picksol1, picksol2, corrthreshold)

if size(cursol.x,2)>=max(picksol1,picksol2)
    if (cursol.x_sel_CorrRev(picksol1)>=corrthreshold) &&...
       (cursol.x_sel_CorrRevLI(picksol1)>=corrthreshold)
       picksol = picksol1;
    else
        % second solution type passess both threshlds - take
        % this one
        if (cursol.x_sel_CorrRev(picksol2)>=corrthreshold) &&...
           (cursol.x_sel_CorrRevLI(picksol2)>=corrthreshold)
            picksol = picksol2;
        else
            % first solution passes total corr threshold but
            % not LI corr threshold
            if (cursol.x_sel_CorrRev(picksol1)>=corrthreshold) 
                picksol = picksol1;
            else
                % second solution passes the threshold
                if (cursol.x_sel_CorrRev(picksol2)>=corrthreshold) 
                    picksol = picksol2;
                else
                    % second solution is not better - keep the
                    % first one
                    picksol = picksol1;
                end
            end
        end
    end
     
    bestsol = cursol.x(:, picksol);
else
    % only one type of solution exist for this metabolite,
    % pick this one
    if size(cursol.x,2)>=min(picksol1,picksol2)
        bestsol = cursol.x(:, min(picksol1,picksol2));
        % pick total correlation and not selection value
        %model_corr(i) = met_bestsols{i}.selection_value(min(picksol1,picksol2));
        picksol = min(picksol1,picksol2);
    else
        % both types do not exist - set to 0
        bestsol = zeros(size(cursol.x,1));
        picksol = 0;
    end
end