close all; 

filename = 'Bleed+Cancellation+listening+test_June+24,+2021_00.05.csv';
opts = detectImportOptions(filename);
% preview(filename,opts)

%%
progressID = 4;
opts.SelectedVariableNames = progressID;
progress = readmatrix(filename, opts);
notCompleted = find(progress < 100);
% notCompleted = [2; notCompleted];   %eliminating subject because of wrong rating in training question
nSubjects = length(progress) - length(notCompleted);


snames = ["Kick", "Snare", "Rack tom", "Floor tom"];
startId = 17;        %columns for kick, snare, floor, tom
randAPS = [true, true, true, true];       % have choices in the question been randomized?
randIPS = [false, true, true, true];
nChoices = 6;
plotLabels = {'Recorded', 'MCWF', 'MLE $\sigma$=1', 'MLE $\sigma$=100', ...
    'MAP $\sigma$=1','MAP $\sigma$=100'};
totalAPS = zeros(length(snames)*nSubjects, nChoices);
totalIPS = zeros(length(snames)*nSubjects, nChoices);

for k = 1:length(snames)
    
    % read data for APS
    apsVars = startId + (0:nChoices-1);     %columns corresponding to drum part APS
    opts.SelectedVariableNames = apsVars;
    APS = readmatrix(filename, opts);
    % remove participants who did not complete the test
    APS(notCompleted,:) = [];
    
    if randAPS(k)
        choiceOrderAPSVars = apsVars(end)+(1:nChoices);
        opts.SelectedVariableNames = choiceOrderAPSVars;
        choiceOrderAPS = readmatrix(filename, opts);
        choiceOrderAPS(notCompleted,:) = [];
        sortedAPS = sort_by_choice_order(APS, choiceOrderAPS);
    else
        choiceOrderAPSVars = apsVars(end);
        sortedAPS = APS; 
    end
    sortedAPS(isnan(sortedAPS)) = 0;
    totalAPS((k-1)*nSubjects+1:k*nSubjects,:) = sortedAPS;
    
    % plot
    figure(1);
    boxplot(sortedAPS, 'Labels', plotLabels);
    bp = gca;
    bp.XAxis.TickLabelInterpreter = 'latex';
    ylim([0, 110]);
    title((snames(k) + ' APS'));
    
    
    % read data for IPS
    ipsVars = choiceOrderAPSVars(end)+(1:nChoices);     %columns corresponding to drum part IPS
    opts.SelectedVariableNames = ipsVars;
    IPS = readmatrix(filename, opts);
    % remove participants who did not complete the test
    IPS(notCompleted,:) = [];
    
   
    if randIPS(k)
        choiceOrderIPSVars = ipsVars(end) + (1:nChoices);
        opts.SelectedVariableNames = choiceOrderIPSVars;
        choiceOrderIPS = readmatrix(filename, opts);
        choiceOrderIPS(notCompleted,:) = [];
        sortedIPS = sort_by_choice_order(IPS, choiceOrderIPS);
    else
        choiceOrderIPSVars = ipsVars(end);
        sortedIPS = IPS;  
    end
    sortedIPS(isnan(sortedIPS)) = 0;
    totalIPS((k-1)*nSubjects+1:k*nSubjects,:) = sortedIPS;

    
    % plot
    figure(2);
    boxplot(sortedIPS, 'Labels', plotLabels);
    bp = gca;
    bp.XAxis.TickLabelInterpreter = 'latex';
    ylim([0, 110]);
    title((snames(k) + ' IPS'));
    
    startId = choiceOrderIPSVars(end)+1;

end

%% TO-DO - ANOVA test to ensure all methods give significantly different results
% conditions for anova - normal distributions, same sample sizes,
% independent samples, equal population variances (?)

close all;
% figure('Units','inches', 'Position',[0 0 2.9 3.3],'PaperPositionMode','auto');
% set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',7, 'FontName','Times');
[p1,tbl1,stats1] = anova1(totalAPS);
bp = gca;
bp.XAxis.TickLabels = plotLabels;
bp.XAxis.TickLabelInterpreter = 'latex';
ylim([0, 100]);ylabel('Rating');
title('Perceptual quality of target signal');
set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',10, 'FontName','Times');
print('listening_test_APS.eps', '-depsc');

    
%%
close all;
% figure('Units','inches', 'Position',[0 0 2.9 3.3],'PaperPositionMode','auto');
% set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',7, 'FontName','Times');
[p2,tbl2,stats2] = anova1(totalIPS);
bp = gca;
bp.XAxis.TickLabels = plotLabels;
bp.XAxis.TickLabelInterpreter = 'latex';
ylim([0, 100]);ylabel('Rating');
title('Amount of interference cancellation');
set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',10, 'FontName','Times');
print('listening_test_IPS.eps', '-depsc');

%% Results - p value signficantly small, all estimators have different means, 
% i.e, they give significantly different results. The notches do not overlap,
% i.e, none of the group medians overlap.