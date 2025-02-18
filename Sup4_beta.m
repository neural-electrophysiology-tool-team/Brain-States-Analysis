%% Supplamentary figure 4: 

% plot the over night beta power for regular and stim nights - red and LED
SA.setCurrentRecording('Animal=PV161,recNames=Night18');
SA.plotDelta2BetaRatio('stim',1,'stimCh',15)
SA.setCurrentRecording('Animal=PV161,recNames=Night17');
SA.plotDelta2BetaRatio

SA.setCurrentRecording('Animal=PV157,recNames=Night18');
SA.plotDelta2BetaRatio('stim',1,'stimCh',15)

%% beta during full cycle + beta during SWS.

%% plot one night: Beta during full cycle
 
% for i = 1:height(stimTable)
i =22 ;
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
SA.setCurrentRecording(recName);

dbStimData = {1./stimTable.betaFull(i).Pre 1./stimTable.betaFull(i).Stim 1./stimTable.betaFull(i).Post};

fdbDec = figure;
groupNames = ["Pre", "Stim","Post"];
colors = [0.2, 0.6, 0.8; 0.9, 0.4, 0.3; 0.5, 0.8, 0.5];

yswarm = [dbStimData{1} dbStimData{2} dbStimData{3}];
xswarm=[1*ones(1,length(dbStimData{1})), 2*ones(1,length(dbStimData{2})), 3*ones(1,length(dbStimData{3}))];
colorswarm = [repmat([0.2, 0.6, 0.8;], length(dbStimData{1}), 1);  % Red for Array 1
    repmat([0.9, 0.4, 0.3], length(dbStimData{2}), 1);  % Green for Array 2
    repmat([0.5, 0.8, 0.5], length(dbStimData{3}), 1)]; % Blue for Array 3
swarmchart(xswarm,yswarm,10,colorswarm,"filled",'o','XJitterWidth',0.5)
ylabel('1/beta power during full cycle')
% ylim([0 450])
set(gca, 'XTick', 1:numel(dbStimData), 'XTickLabel', groupNames);
title('1/beta power during full cycles - one night')

% statistics:

% Group names
groupNames = {'Pre', 'During', 'After'};

% Combine data into a single vector and create a grouping variable
allData = [dbStimData{1},dbStimData{2},dbStimData{3}]'; % Concatenate all data
groupLabels = [...
    repmat({'Pre'}, numel(dbStimData{1}), 1); ...
    repmat({'During'}, numel(dbStimData{2}), 1); ...
    repmat({'After'}, numel(dbStimData{3}), 1)]';

% Kruskal-Wallis test
[pKruskal, tblKruskal, statsKruskal] = kruskalwallis(allData, groupLabels, 'off');
fprintf('Kruskal-Wallis test p-value: %.4f\n', pKruskal);
annotation('textbox', [0.1, 0.8, 0.2, 0.1], 'String', ...
    sprintf('Kruskal-Wallis test p-value: %.4f\n', pKruskal), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');


if pKruskal < 0.05
    disp('Significant differences found in Kruskal-Wallis test. Proceeding with pairwise comparisons...');
    pMannWhitneyAll = [];
    % Pairwise comparisons with Mann-Whitney U test
    pairwiseComparisons = combnk(1:3, 2); % All pairwise combinations
    numComparisons = height(pairwiseComparisons);
    alpha = 0.05;
    correctedAlpha = alpha / numComparisons;

    for i = 1:size(pairwiseComparisons, 1)
        group1 = dbStimData{pairwiseComparisons(i, 1)};
        group2 = dbStimData{pairwiseComparisons(i, 2)};
        
        % Mann-Whitney U test (rank-sum test)
        pMannWhitney = ranksum(group1, group2);
        fprintf('Mann-Whitney U test (%s vs %s): p-value = %.4f\n', ...
                groupNames{pairwiseComparisons(i, 1)}, ...
                groupNames{pairwiseComparisons(i, 2)}, ...
                pMannWhitney);
        pMannWhitneyAll = [pMannWhitneyAll, pMannWhitney]; %for later in the caption

        annotation('textbox', [0.1+0.2*i, 0.8, 0.2, 0.1], 'String', ...
            sprintf('MW (%s vs %s): p = %.4f\n', ...
            groupNames{pairwiseComparisons(i, 1)}, ...
            groupNames{pairwiseComparisons(i, 2)}, ...
            pMannWhitney), 'EdgeColor', 'none', 'HorizontalAlignment', ...
            'right', 'VerticalAlignment', 'middle');

    end
else
    disp('No significant differences found in Kruskal-Wallis test.');
end


%savefigure
set(fdbDec,'PaperPositionMode','auto');
fileName=[SA.currentPlotFolder filesep 'betafullcycpreStimPost'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);
set(fdbDec,'PaperPosition',[1 1 2.7 2.3]);
fileName=[analysisFolder filesep 'BetafullcycOneNight'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);


%% plot one night: Beta during full cycle
 
% for i = 1:height(stimTable)
i =22 ;
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
SA.setCurrentRecording(recName);

dbStimData = {1./stimTable.betaSW(i).Pre 1./stimTable.betaSW(i).Stim 1./stimTable.betaSW(i).Post};

fdbDec = figure;
groupNames = ["Pre", "Stim","Post"];
colors = [0.2, 0.6, 0.8; 0.9, 0.4, 0.3; 0.5, 0.8, 0.5];

yswarm = [dbStimData{1} dbStimData{2} dbStimData{3}];
xswarm=[1*ones(1,length(dbStimData{1})), 2*ones(1,length(dbStimData{2})), 3*ones(1,length(dbStimData{3}))];
colorswarm = [repmat([0.2, 0.6, 0.8;], length(dbStimData{1}), 1);  % Red for Array 1
    repmat([0.9, 0.4, 0.3], length(dbStimData{2}), 1);  % Green for Array 2
    repmat([0.5, 0.8, 0.5], length(dbStimData{3}), 1)]; % Blue for Array 3
swarmchart(xswarm,yswarm,10,colorswarm,"filled",'o','XJitterWidth',0.5)
ylabel('beta power during SWS')
% ylim([0 450])
set(gca, 'XTick', 1:numel(dbStimData), 'XTickLabel', groupNames);
title('Beta power during SWS - one night')

% statistics:

% Group names
groupNames = {'Pre', 'During', 'After'};

% Combine data into a single vector and create a grouping variable
allData = [dbStimData{1},dbStimData{2},dbStimData{3}]'; % Concatenate all data
groupLabels = [...
    repmat({'Pre'}, numel(dbStimData{1}), 1); ...
    repmat({'During'}, numel(dbStimData{2}), 1); ...
    repmat({'After'}, numel(dbStimData{3}), 1)]';

% Kruskal-Wallis test
[pKruskal, tblKruskal, statsKruskal] = kruskalwallis(allData, groupLabels, 'off');
fprintf('Kruskal-Wallis test p-value: %.4f\n', pKruskal);
annotation('textbox', [0.1, 0.8, 0.2, 0.1], 'String', ...
    sprintf('Kruskal-Wallis test p-value: %.4f\n', pKruskal), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');


if pKruskal < 0.05
    disp('Significant differences found in Kruskal-Wallis test. Proceeding with pairwise comparisons...');
    pMannWhitneyAll = [];
    % Pairwise comparisons with Mann-Whitney U test
    pairwiseComparisons = combnk(1:3, 2); % All pairwise combinations
    numComparisons = height(pairwiseComparisons);
    alpha = 0.05;
    correctedAlpha = alpha / numComparisons;

    for i = 1:size(pairwiseComparisons, 1)
        group1 = dbStimData{pairwiseComparisons(i, 1)};
        group2 = dbStimData{pairwiseComparisons(i, 2)};
        
        % Mann-Whitney U test (rank-sum test)
        pMannWhitney = ranksum(group1, group2);
        fprintf('Mann-Whitney U test (%s vs %s): p-value = %.4f\n', ...
                groupNames{pairwiseComparisons(i, 1)}, ...
                groupNames{pairwiseComparisons(i, 2)}, ...
                pMannWhitney);
        pMannWhitneyAll = [pMannWhitneyAll, pMannWhitney]; %for later in the caption

        annotation('textbox', [0.1+0.2*i, 0.8, 0.2, 0.1], 'String', ...
            sprintf('MW (%s vs %s): p = %.4f\n', ...
            groupNames{pairwiseComparisons(i, 1)}, ...
            groupNames{pairwiseComparisons(i, 2)}, ...
            pMannWhitney), 'EdgeColor', 'none', 'HorizontalAlignment', ...
            'right', 'VerticalAlignment', 'middle');

    end
else
    disp('No significant differences found in Kruskal-Wallis test.');
end


%savefigure
set(fdbDec,'PaperPositionMode','auto');
fileName=[SA.currentPlotFolder filesep 'betaSWpreStimPost'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);
set(fdbDec,'PaperPosition',[1 1 2.7 2.3]);
fileName=[analysisFolder filesep 'BetaSWOneNight'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);



%% plot all nights: Beta full cycle - only red

type = 'Red';
wavelength = '635';
curTrials = contains(stimTable.Remarks,wavelength)...
        & ~contains(stimTable.Remarks,'Ex') ...
        & ~any(isnan(stimTable.betaMeans),2); %& contains(stimTable.Animal,curAni);
n = sum(curTrials);
N = length(unique(stimTable.Animal(curTrials)));
groupNames = {'Pre', 'During', 'After'};
curData = 1./stimTable.betaMeans(curTrials,:);


% statistics:

% Assuming data in columns where each row is a subject and each column is a timepoint
% bdSWDataStat = stimTable.dbSWMeans(curTrials,:);
% cleanDBSW = bdSWDataStat(~any(isnan(bdSWDataStat), 2), :);
% n=height(cleanDBSW);

[p, tbl, stats] = friedman(curData, 1,'off'); % Here, 1 indicates within-subjects design
fprintf('p-value for freidman ANOVA test: %.5f\n',p)
% p-valure is very low, post hoc:
% Example data for three groups
before = curData(:,1);
during = curData(:,2);
after = curData(:,3);

% Bonferroni-corrected alpha level
alpha = 0.05 / 3;

% Pairwise Wilcoxon signed-rank tests
[p_before_during, ~, stats_before_during] = signrank(before, during);
[p_during_after, ~, stats_during_after] = signrank(during, after);
[p_after_before, ~, stats_after_before] = signrank(after, before);

% Display results with Bonferroni correction
fprintf('Wilcoxon signed-rank test results with Bonferroni correction:\n');
fprintf('Before vs During: p-value = %.4f (Significant if < %.4f)\n', p_before_during, alpha);
fprintf('During vs After: p-value = %.4f (Significant if < %.4f)\n', p_during_after, alpha);
fprintf('After vs Before: p-value = %.4f (Significant if < %.4f)\n', p_after_before, alpha);

%plot
fdb = figure;
% plot(stimTable.dbSWMeans(curTrials,:)','Color',[0.5 0.5 0.5],'Marker','.','MarkerSize',10)
%create color code:

[~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
curColorMat = animalsColors(animalIndices, :); 
hold on; 
for i = 1:height(curData)
    plot(curData(i,:),'Color',curColorMat(i,:),'Marker','.','MarkerSize',10)
end
plot(mean(curData,1,'omitnan'),'Color','k','LineWidth',2,'Marker','.','MarkerSize',10)
xlim([0.5, 3.5])
% grid on;
% ylim([0 450])
xticks(1:3)
xticklabels(groupNames)
ylabel('1/Beta means during full cycle')

annotation('textbox', [0.8, 0.85, 0.03, 0.1], 'String', ...
    sprintf('n=%i,N=%i',n,N), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

annotation('textbox', [0.1, 0.8, 0.4, 0.1], 'String', ...
    sprintf('p-value for Friedman ANOVA test: %.5f',p), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

annotation('textbox', [0.15, 0.65, 0.25, 0.1], 'String', ...
    sprintf('p-value = %.4f (Significant if < %.4f)', p_before_during, alpha), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

annotation('textbox', [0.55, 0.65, 0.25, 0.1], 'String', ...
    sprintf('p-value = %.4f', p_during_after), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

annotation('textbox', [0.3, 0.1, 0.25, 0.1], 'String', ...
    sprintf('p-value = %.4f', p_after_before), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

% savefigure
set(fdb,'PaperPositionMode','auto');
fileName=[analysisFolder filesep 'BetachangeFullCycle'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);
%% plot all nights: Beta full cycle - only red

type = 'Red';
wavelength = '635';
curTrials = contains(stimTable.Remarks,wavelength)...
        & ~contains(stimTable.Remarks,'Ex') ...
        & ~any(isnan(stimTable.betaMeans),2); %& contains(stimTable.Animal,curAni);
n = sum(curTrials);
N = length(unique(stimTable.Animal(curTrials)));
groupNames = {'Pre', 'During', 'After'};
curData = 1./stimTable.betaSWMeans(curTrials,:);


% statistics:

% Assuming data in columns where each row is a subject and each column is a timepoint
% bdSWDataStat = stimTable.dbSWMeans(curTrials,:);
% cleanDBSW = bdSWDataStat(~any(isnan(bdSWDataStat), 2), :);
% n=height(cleanDBSW);

[p, tbl, stats] = friedman(curData, 1,'off'); % Here, 1 indicates within-subjects design
fprintf('p-value for freidman ANOVA test: %.5f\n',p)
% p-valure is very low, post hoc:
% Example data for three groups
before = curData(:,1);
during = curData(:,2);
after = curData(:,3);

% Bonferroni-corrected alpha level
alpha = 0.05 / 3;

% Pairwise Wilcoxon signed-rank tests
[p_before_during, ~, stats_before_during] = signrank(before, during);
[p_during_after, ~, stats_during_after] = signrank(during, after);
[p_after_before, ~, stats_after_before] = signrank(after, before);

% Display results with Bonferroni correction
fprintf('Wilcoxon signed-rank test results with Bonferroni correction:\n');
fprintf('Before vs During: p-value = %.4f (Significant if < %.4f)\n', p_before_during, alpha);
fprintf('During vs After: p-value = %.4f (Significant if < %.4f)\n', p_during_after, alpha);
fprintf('After vs Before: p-value = %.4f (Significant if < %.4f)\n', p_after_before, alpha);

%plot
fdb = figure;
% plot(stimTable.dbSWMeans(curTrials,:)','Color',[0.5 0.5 0.5],'Marker','.','MarkerSize',10)
%create color code:

[~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
curColorMat = animalsColors(animalIndices, :); 
hold on; 
for i = 1:height(curData)
    plot(curData(i,:),'Color',curColorMat(i,:),'Marker','.','MarkerSize',10)
end
plot(mean(curData,1,'omitnan'),'Color','k','LineWidth',2,'Marker','.','MarkerSize',10)
xlim([0.5, 3.5])
% grid on;
% ylim([0 450])
xticks(1:3)
xticklabels(groupNames)
ylabel('1/Beta means during full cycle')

annotation('textbox', [0.8, 0.85, 0.03, 0.1], 'String', ...
    sprintf('n=%i,N=%i',n,N), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

annotation('textbox', [0.1, 0.8, 0.4, 0.1], 'String', ...
    sprintf('p-value for Friedman ANOVA test: %.5f',p), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

annotation('textbox', [0.15, 0.65, 0.25, 0.1], 'String', ...
    sprintf('p-value = %.4f (Significant if < %.4f)', p_before_during, alpha), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

annotation('textbox', [0.55, 0.65, 0.25, 0.1], 'String', ...
    sprintf('p-value = %.4f', p_during_after), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

annotation('textbox', [0.3, 0.1, 0.25, 0.1], 'String', ...
    sprintf('p-value = %.4f', p_after_before), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

% savefigure
set(fdb,'PaperPositionMode','auto');
fileName=[analysisFolder filesep 'BetachangeSWS'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);




%% plot D/B decrease full cycles - all nights

stimType = ["Blue","Green","Red","WhiteEx"];
stimWaveL = ["47","532","635","LED"];
% plotColors = {[0 0.586 0.9766],[0.05 0.81 0.379],[1 0.27 0.27], [0.5 0.5 0.5]};
numType = length(stimType);
x=1:3;
%plot Period Times:
f=figure;
set(f, 'Position', [100, 100, 800, 400]);
hold on
for type = 1:numType
    %plot the data
    subplot(1,4,type)
    curType = stimWaveL(type);
    curName = stimType(type);
    curTrials = contains(stimTable.Remarks,curType) &...
                ~contains(stimTable.Remarks,"Ex") &...
                all(~isnan(stimTable.dbMeans), 2);
    n = sum(curTrials);
    N = numel(unique(stimTable.Animal(curTrials)));
    % curCol = plotColors{type};
    curData = stimTable.dbMeans(curTrials,:);
    curMean = mean(stimTable.dbMeans(curTrials,:),1,'omitnan');
    statsdb = struct();
    %statistics:

    [p, tbl, stats] = friedman(curData, 1,'off'); % Here, 1 indicates within-subjects design
    fprintf('p-value for freidman ANOVA test: %.5f\n',p)
    % % p-valure is very low, post hoc:
    % data for the four groups
    beforedb = curData(:,1);
    duringdb = curData(:,2);
    afterdb = curData(:,3);

    % Bonferroni-corrected alpha level
    alpha = 0.05 / 3;

    % Pairwise Wilcoxon signed-rank tests
    [p_pre_during, ~, stats_before_during] = signrank(beforedb, duringdb);
    [p_during_post, ~, stats_during_after] = signrank(duringdb, afterdb);
    [p_pre_post, ~, stats_wake_during] = signrank(beforedb, afterdb);

    % Display results with Bonferroni correction
    fprintf('Wilcoxon signed-rank test results with Bonferroni correction:\n');
    fprintf('Before vs During: p-value = %.4f (Significant if < %.4f)\n', p_pre_during, alpha);
    fprintf('During vs After: p-value = %.4f (Significant if < %.4f)\n', p_during_post, alpha);
    fprintf('Pre vs. Post: p-value = %.4f (Significant if < %.4f)\n', p_pre_post, alpha);

    %save statistics:
    statsdb.(curName).alpha = alpha;
    statsdb.(curName).pAnova = p;
    statsdb.(curName).p_pre_post = p_pre_post;
    statsdb.(curName).p_pre_during = p_pre_during;
    statsdb.(curName).p_during_post = p_during_post;


    if n>0
        
        [~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
        curColorMat = animalsColors(animalIndices, :); 
        hold on
        for j = 1:height(curData)
            plot(x,curData(j,:),'Color',curColorMat(j,:),'Marker', '.')
        end
            plot(x,curMean,'color','k','LineWidth',3,'Marker', '.')
        
        %     % plot ExLight:
        %     ExTrials = contains(stimTable.Remarks,curType) & ...
        %         contains(stimTable.Remarks,'Ex') &...
        %         all(~isnan(stimTable.dbSWMeans),2);
        %     if sum(ExTrials)>0
        %         plot(x,stimTable.dbSWMeans(ExTrials,:),'Color',[0.2 0.2 0.2],'LineStyle','--','Marker','.')
        %     end
        % hold off

        annotation('textbox', [.05 + 0.202*type 0.85, 0.03, 0.1], 'String', ...
            sprintf('n=%i,N=%i',n,N), 'EdgeColor', 'none', 'HorizontalAlignment', ...
            'right', 'VerticalAlignment', 'middle');
    end
    ylabel('Time[s]')
    ylim([0 600])
    xticklabels({'Pre','During','Post'})
    xticks(1:3); xlim([0.5 3.5])
end

sgtitle ('D/B decrease full cycles according to wavelangth ')


% savefigure
set(gcf,'PaperPositionMode','auto');
fileName=[analysisFolder filesep 'DBfullcycleAllnigthscolors'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);

%% plot D/B decrease full cycles - all nights

stimType = ["Blue","Green","Red","WhiteEx"];
stimWaveL = ["47","532","635","LED"];
% plotColors = {[0 0.586 0.9766],[0.05 0.81 0.379],[1 0.27 0.27], [0.5 0.5 0.5]};
numType = length(stimType);
x=1:3;
%plot Period Times:
f=figure;
set(f, 'Position', [100, 100, 800, 400]);
hold on
for type = 1:numType
    %plot the data
    subplot(1,4,type)
    curType = stimWaveL(type);
    curName = stimType(type);
    curTrials = contains(stimTable.Remarks,curType) &...
                ~contains(stimTable.Remarks,"Ex") &...
                all(~isnan(stimTable.dbMeans), 2);
    n = sum(curTrials);
    N = numel(unique(stimTable.Animal(curTrials)));
    % curCol = plotColors{type};
    curData = stimTable.betaMeans(curTrials,:);
    curMean = mean(stimTable.betaMeans(curTrials,:),1,'omitnan');
    statsbeta = struct();
    %statistics:

    [p, tbl, stats] = friedman(curData, 1,'off'); % Here, 1 indicates within-subjects design
    fprintf('p-value for freidman ANOVA test: %.5f\n',p)
    % % p-valure is very low, post hoc:
    % data for the four groups
    beforedb = curData(:,1);
    duringdb = curData(:,2);
    afterdb = curData(:,3);

    % Bonferroni-corrected alpha level
    alpha = 0.05 / 3;

    % Pairwise Wilcoxon signed-rank tests
    [p_pre_during, ~, stats_before_during] = signrank(beforedb, duringdb);
    [p_during_post, ~, stats_during_after] = signrank(duringdb, afterdb);
    [p_pre_post, ~, stats_wake_during] = signrank(beforedb, afterdb);

    % Display results with Bonferroni correction
    fprintf('Wilcoxon signed-rank test results with Bonferroni correction:\n');
    fprintf('Before vs During: p-value = %.4f (Significant if < %.4f)\n', p_pre_during, alpha);
    fprintf('During vs After: p-value = %.4f (Significant if < %.4f)\n', p_during_post, alpha);
    fprintf('Pre vs. Post: p-value = %.4f (Significant if < %.4f)\n', p_pre_post, alpha);

    %save statistics:
    statsbeta.(curName).alpha = alpha;
    statsbeta.(curName).pAnova = p;
    statsbeta.(curName).p_pre_post = p_pre_post;
    statsbeta.(curName).p_pre_during = p_pre_during;
    statsbeta.(curName).p_during_post = p_during_post;


    if n>0
        
        [~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
        curColorMat = animalsColors(animalIndices, :); 
        hold on
        for j = 1:height(curData)
            plot(x,curData(j,:),'Color',curColorMat(j,:),'Marker', '.')
        end
            plot(x,curMean,'color','k','LineWidth',3,'Marker', '.')
        
        %     % plot ExLight:
        %     ExTrials = contains(stimTable.Remarks,curType) & ...
        %         contains(stimTable.Remarks,'Ex') &...
        %         all(~isnan(stimTable.dbSWMeans),2);
        %     if sum(ExTrials)>0
        %         plot(x,stimTable.dbSWMeans(ExTrials,:),'Color',[0.2 0.2 0.2],'LineStyle','--','Marker','.')
        %     end
        % hold off

        annotation('textbox', [.05 + 0.202*type 0.85, 0.03, 0.1], 'String', ...
            sprintf('n=%i,N=%i',n,N), 'EdgeColor', 'none', 'HorizontalAlignment', ...
            'right', 'VerticalAlignment', 'middle');
    end
    ylabel('Time[s]')
    ylim([0 600])
    xticklabels({'Pre','During','Post'})
    xticks(1:3); xlim([0.5 3.5])
end

sgtitle ('Beta decrease full cycles according to wavelangth ')


% savefigure
set(gcf,'PaperPositionMode','auto');
fileName=[analysisFolder filesep 'BetafullcycleAllnigthscolors'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);
