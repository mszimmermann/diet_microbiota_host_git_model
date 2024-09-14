function gitfit_diag_plot(gitfit, metname, fig)

clf(fig)

spi = 2;
spj = 7;
spidx = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot original data
subplot(spi,spj,spidx)
plot(gitfit.kmeanMatrix_joint_orig')
xlim([1,size(gitfit.kmeanMatrix_joint_orig,2)])
title('Original data')
lh=legend(cellfun(@(x) strrep(x, '_','-'),...
            gitfit.kmeanMatrix_joint_names(:,1), 'unif', 0),...
            'Location', 'West');
set(lh,'position',[.0 .75 .1 .1])
spidx = spidx+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot distribution of solutions
subplot(spi,spj,spidx)
testx_norm = gitfit.testx./max(abs(gitfit.testx));
boxplot(testx_norm')
hold on
% plot interior point solution
testx_norm = gitfit.x_ip./max(abs(gitfit.x_ip));
plot((1:length(testx_norm)), testx_norm', 'g*')
title('All solutions')
ylim([-1.2 1.2])
spidx = spidx+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot distribution of correlations of trust-region and shuffled
subplot(spi,spj,spidx)
histogram(gitfit.testCorrRev)
hold on
histogram(gitfit.testCorrRev_shuffled)
plot([gitfit.x_ip_CorrRev,gitfit.x_ip_CorrRev], [0, 50], 'g')
legend({'Total TR', 'Shuffled', 'IP'})
title({'Corr rec vs orig',...
        sprintf('mn=%.2f md=%.2f n7=%d',...
                 mean(gitfit.testCorrRev),...
                 median(gitfit.testCorrRev),...
                 nnz(gitfit.testCorrRev>0.7)),...
        sprintf('SH mn=%.2f md=%.2f n7=%d',...
                 mean(gitfit.testCorrRev_shuffled),...
                 median(gitfit.testCorrRev_shuffled),...
                 nnz(gitfit.testCorrRev_shuffled>0.7))})
xlim([-1 1])
spidx = spidx+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot distribution of correlations of of total, SI and LI
subplot(spi,spj,spidx)
histogram(gitfit.testCorrRev)
hold on
histogram(gitfit.testCorrRevSI)
histogram(gitfit.testCorrRevLI)
plot([gitfit.x_ip_CorrRevSI,gitfit.x_ip_CorrRevSI], [0, 50], 'g--') 
plot([gitfit.x_ip_CorrRevLI,gitfit.x_ip_CorrRevLI], [0, 50], 'g.-') 
legend({'Total', 'SI', 'LI'})
title('Corr SI/LI')
xlim([-1 1])
spidx = spidx+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot distribution of selected solutions
% get the best x
x_selected = select_gitfit_sol(gitfit);
% get the best reciprocals
x_selected_dataR = zeros(numel(gitfit.kmeanMatrix_joint_orig),...
                         size(x_selected.x,2));
options = optimoptions(@lsqlin,'Display', 'off',...
    'Algorithm','interior-point','MaxIterations',1500);
for i=1:size(x_selected.x,2)
    x = x_selected.x(:,i);
    if length(x) < 9 % duplicate to restore values 
        x = [x(1); x(2:4); x(2:4); x(5); x(5)];
    end
    % calculate reverse problem
    [Ra,rb] = calculateRAmatrix_final(x*1000);
    dataR = lsqlin(Ra,rb,[],[],[],[],zeros(1,size(Ra,2)), [], [], options);
    x_selected_dataR(:,i) = dataR;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot distribution of selected solutions
subplot(spi,spj,spidx)
% take only metabolic coefs
x = x_selected.x(2:end,:);
% normalize by max in column
x = bsxfun(@rdivide,x,max(abs(x),[],1));
boxplot(x')
title('Best solutions')
ylim([-1.2 1.2])
spidx = spidx + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pick sum of corr solution
subplot(spi,spj,spidx)
% take only metabolic coefs
x = x_selected.x(2:end,ismember(x_selected.selection_criterion,...
                                {'sum PCC'}));
% normalize by max in column
x = x/max(abs(x));
bar(x')
title({'Sel sum PCC',...
    sprintf('Mean PCC=%.2f', ...
    x_selected.selection_value(ismember(x_selected.selection_criterion,...
                                {'sum PCC'}))/3)})
ylim([-1.2 1.2])
spidx = spidx + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pick LI within total corr solution
subplot(spi,spj,spidx)
% take only metabolic coefs
picksol = ismember(x_selected.selection_criterion,...
                                {'LI PCC within high total'});
if size(x_selected.x,2)==length(picksol)
    x = x_selected.x(2:end,picksol);
    % normalize by max in column
    x = x/max(abs(x));
    bar(x')
    title({'LI within total PCC',...
        sprintf('Mean PCC=%.2f', ...
        mean([x_selected.x_sel_CorrRev(picksol)...
              x_selected.x_sel_CorrRevSI(picksol)...
              x_selected.x_sel_CorrRevLI(picksol)]))})
else
    x = zeros(size(x_selected.x,1));
    bar(x')
    title('No sol selected')
end

ylim([-1.2 1.2])
spidx = spidx + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot restored data for selected solutions
spidx = spj+1;
for i=1:size(x_selected_dataR,2)
    dataR = x_selected_dataR(:,i);
    dataR = reshape(dataR,[],4)';
    subplot(spi, spj, spidx)
    plot(dataR')
    xlim([1,size(gitfit.kmeanMatrix_joint_orig,2)])
    title({x_selected.selection_criterion{i}
        sprintf('%.2f %.2f %.2f', ...
            x_selected.x_sel_CorrRev(i),...
            x_selected.x_sel_CorrRevSI(i),...
            x_selected.x_sel_CorrRevLI(i))})
    spidx = spidx+1;
end

% subplot(spi,spj,3)
% %testx_norm = testx(2:end,:)./max(abs(testx(2:end,:)));
% testx_norm = gitfit.testx(:,gitfit.testCorrRev>0.7)./...
%     max(abs(gitfit.testx(:,gitfit.testCorrRev>0.7)));
% boxplot(testx_norm', 'PlotStyle','compact')

sgtitle(strrep(metname,'_','-'));
