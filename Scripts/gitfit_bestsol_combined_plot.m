function gitfit_bestsol_combined_plot(met_info_combined, met_bestsol_combined,...
                                      figFileName)

fig = figure('units','normalized','outerposition',[0 0 1 1]);
spi = 1;
spj = 3;

% define colora and GIT section names for plotting
mycolors = [0 115 178;... %dark blue
            211 96 39]/256;%dark orange
%mycolors = [0 115 178;... %dark blue
%            ]/256;%dark orange
git_labels = {'Du', 'Je', 'Il', 'Cec', 'Col', 'Fec'};
mylinestyles = ["-"; "--"];

for i=1:length(met_info_combined.CompoundName) 
    clf(fig)

    spidx = 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot original data
    plot_data = reshape(met_bestsol_combined.kmeanMatrix_joint_orig(i,:)', [],6)';
    plot_data_names = reshape(met_bestsol_combined.kmeanMatrix_joint_names, [], 6)';
    subplot(spi,spj,spidx)
 
    plot(plot_data)

    ax = gca; 
    ax.ColorOrder = mycolors;
    ax.LineStyleOrder = mylinestyles;
    ax.LineStyleCyclingMethod = "beforecolor";

    xlim([1,size(plot_data,1)])
    xticks(gca, 1:length(git_labels))
    xticklabels(gca, git_labels)
    title('Original data')
    ylabel('Intensity, normalized')
    ylim([0 1])
    lh=legend(cellfun(@(x) strrep(x, '_','-'),...
                plot_data_names(1,:), 'unif', 0),...
                'Location', 'West');
    set(lh,'position',[.0 .75 .1 .1])
    axis square
    spidx = spidx+1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % pick sum of corr solution
    subplot(spi,spj,spidx)
    % take only metabolic coefs
    x = met_bestsol_combined.x(i, 2:end);
    % normalize by max in column
    x = x/max(abs(x));
    x_coeflabels = met_bestsol_combined.coefvalues(2:end); 
    barh(x')
    yticks(gca, 1:length(x_coeflabels))
    yticklabels(gca, x_coeflabels)
    title({met_bestsol_combined.selection_criteria{i},...
        sprintf('Total PCC=%.2f', ...
        met_bestsol_combined.x_sel_CorrRev(i))})
    xlim([-1 1])
    set(gca, 'YDir', 'reverse')
    axis square
    spidx = spidx + 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot restored data for selected solutions
    dataR = met_bestsol_combined.x_sel_dataR(i,:);
    dataR = reshape(dataR,[],6)';

    %normalize rdata to max
    dataR = (dataR-0.8);
    dataR = dataR/max(max(dataR));

    subplot(spi, spj, spidx)

    plot(dataR)

    ax = gca; 
    ax.ColorOrder = mycolors;
    ax.LineStyleOrder = mylinestyles;
    ax.LineStyleCyclingMethod = "beforecolor";

    xlim([1,size(plot_data_names,1)])
    xticks(gca, 1:length(git_labels))
    xticklabels(gca, git_labels)
    title(sprintf('CorRev=%.2f SI=%.2f LI=%.2f', ...
            met_bestsol_combined.x_sel_CorrRev(i),...
            met_bestsol_combined.x_sel_CorrRevSI(i),...
            met_bestsol_combined.x_sel_CorrRevLI(i)))
    ylabel('Restored intensity, normalized')
    ylim([0 1])
    axis square

    sgtitle(strrep(met_info_combined.CompoundName{i},'_','-'));

    % print to file
     orient landscape
    %print to figure
    print(gcf, '-vector', '-dpsc2', '-r600', '-append', '-bestfit',...
            figFileName);

end
