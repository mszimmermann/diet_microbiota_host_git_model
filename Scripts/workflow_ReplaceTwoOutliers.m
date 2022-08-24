function [totDataMat_withnan] = workflow_ReplaceTwoOutliers(totDataMat,...
                                                               totDataDiet,...
                                                               totDataCondition)

diets = {'HFD', 'CTR'};
conditions = {'GF', 'DC'};

totDataMat_withnan = totDataMat;

for diet_i = 1:length(diets)
    for cond_i = 1:length(conditions)
        curcols = ismember(totDataDiet, diets(diet_i)) &...
                  ismember(totDataCondition, conditions(cond_i));
        tempmat = totDataMat(:, curcols);
        for i=1:size(tempmat,1)
             % remove two 5000
             indx5000 = find(tempmat(i,:)==5000);
             if length(indx5000)<=2
                tempmat(i,indx5000) = nan;
             end
        end
        totDataMat_withnan(:,curcols) = tempmat;
    end
end
        
