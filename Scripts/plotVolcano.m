% volcano plot with colorcode
% input: logFC, pvalues, data labels, logFC treshold, pvalue treshold
% speclabel is the ids of points to label specially
function [upreg, downreg, fig] = plotVolcano(logFC, pt, labels, logFCtres, ptres, speclabel, fig)
% some of the plot parameters
fsize = 10; %axis font size
textfontsize = 10; %text font size
markersize = 3; %marker size
ylimValue = 5;
% take the logarithm 10 of the p-values
% set min value instead of 0
pt(pt==0) = min(pt(pt>0));
logpval = -log10(pt);

showlabels = 0;
if isempty(fig)
    fig = figure;
    showlabels = 1;
end
% plot everything in grey
plot(logFC, logpval, 'o', 'MarkerEdgeColor',[.5 .5 .5],...
                            'MarkerFaceColor', [.5 .5 .5],...
                            'MarkerSize', markersize)
hold on

%up-regulated in red
upreg = (logpval > -log10(ptres)) & (logFC > logFCtres);
plot(logFC(upreg), logpval(upreg), 'o', 'MarkerEdgeColor','r',...
                                        'MarkerFaceColor', 'r',...
                                        'MarkerSize', markersize)

% plot boundaries in dashed grey
%plot([min(logFC), max(logFC)], [-log10(ptres), -log10(ptres)], '--', 'Color', [.5 .5 .5])
%text(max(logFC), -log10(ptres), sprintf('P-value<%.2f', ptres), 'FontSize', textfontsize, 'horizontalAlignment', 'right', 'verticalAlignment', 'top')
plot([-10 10], [-log10(ptres), -log10(ptres)], '--', 'Color', [.5 .5 .5])
%text(3, -log10(ptres), sprintf('P-value<%.2f', ptres), 'FontSize', textfontsize, 'horizontalAlignment', 'right', 'verticalAlignment', 'top')
text(1, -log10(ptres), sprintf('Adj.P-value<%.2f', ptres), 'FontSize', textfontsize, 'horizontalAlignment', 'right', 'verticalAlignment', 'top')

plot([-logFCtres, -logFCtres], [0, ylimValue], '--', 'Color', [.5 .5 .5])
plot([logFCtres, logFCtres], [0, ylimValue], '--', 'Color', [.5 .5 .5])

plot([-logFCtres, -logFCtres], [0, max(logpval)], '--', 'Color', [.5 .5 .5])
plot([logFCtres, logFCtres], [0, max(logpval)], '--', 'Color', [.5 .5 .5])
% text(-logFCtres, max(logpval), sprintf('Log2(FC)<%.2f', -logFCtres), 'FontSize', textfontsize, 'horizontalAlignment', 'right')
% text(logFCtres, max(logpval), sprintf('Log2(FC)>%.2f', logFCtres), 'FontSize', textfontsize, 'horizontalAlignment', 'left')
% text(-logFCtres, max(logpval), sprintf('FC<%.2f', -2^logFCtres), 'FontSize', textfontsize, 'horizontalAlignment', 'right')
% text(logFCtres, max(logpval), sprintf('FC>%.2f', 2^logFCtres), 'FontSize', textfontsize, 'horizontalAlignment', 'left')
text(-logFCtres, ylimValue-0.3, sprintf('FC<%.2f', -2^logFCtres), 'FontSize', textfontsize, 'horizontalAlignment', 'right')
text(logFCtres, ylimValue-0.3, sprintf('FC>%.2f', 2^logFCtres), 'FontSize', textfontsize, 'horizontalAlignment', 'left')


% down-regulated in blue
downreg = (logpval > -log10(ptres)) & (logFC < -logFCtres);
plot(logFC(downreg), logpval(downreg), 'o', 'MarkerEdgeColor','b',...
                                            'MarkerFaceColor', 'b',...
                                            'MarkerSize', markersize)

% if there are special labels
if ~isempty(speclabel)
    plot(logFC(speclabel' & (downreg | upreg)), logpval(speclabel' & (downreg | upreg)), 'o', 'MarkerEdgeColor',[0 0.4 0.2],...
                                                'MarkerFaceColor', [0 0.4 0.2],...
                                                'MarkerSize', markersize)
end
% if ~isempty(speclabel2)
%     plot(logFC(speclabel2' & (downreg | upreg)), logpval(speclabel2' & (downreg | upreg)), 'o', 'MarkerEdgeColor','r',...
%                                                 'MarkerFaceColor', 'r',...
%                                                 'MarkerSize', markersize)
% end

if showlabels
    % redefine datacursormode show function
    global POINTLABELS POINTXY
    POINTXY = [logFC, logpval];
    POINTLABELS = labels;

    dcm_obj = datacursormode(fig);
    set(dcm_obj,'UpdateFcn',@myupdatefcn)
    datacursormode on
end

% add x and y labels
xlabel('Log2(Fold Change)', 'fontSize', fsize)
ylabel('-Log10(P-value)', 'fontSize', fsize)
set(gca, 'fontSize', fsize)

% set plot limits
if nnz(logFC)
    xlimabs = max(abs(min(logFC)), max(logFC));
    xlim([-xlimabs, xlimabs])
end