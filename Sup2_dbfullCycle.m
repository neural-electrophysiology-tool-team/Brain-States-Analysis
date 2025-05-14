%% Sup2 - d/b for full cycles
%% plot one night: D/B decrease full cycle
 
% for i = 1:height(stimTable)
i =22 ;
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
SA.setCurrentRecording(recName);

dbStimData = {stimTable.dbFull(i).Pre stimTable.dbFull(i).Stim stimTable.dbFull(i).Post};

fdbDec = figure;
groupNames = ["Pre", "Stim","Post"];
colors = [0.2, 0.6, 0.8; 0.9, 0.4, 0.3; 0.5, 0.8, 0.5];

yswarm = [stimTable.dbFull(i).Pre stimTable.dbFull(i).Stim stimTable.dbFull(i).Post];
xswarm=[1*ones(1,length(stimTable.dbFull(i).Pre)), 2*ones(1,length(stimTable.dbFull(i).Stim)), 3*ones(1,length(stimTable.dbFull(i).Post))];
colorswarm = [repmat([0.2, 0.6, 0.8;], length(stimTable.dbSW(i).Pre), 1);  % Red for Array 1
    repmat([0.9, 0.4, 0.3], length(stimTable.dbSW(i).Stim), 1);  % Green for Array 2
    repmat([0.5, 0.8, 0.5], length(stimTable.dbSW(i).Post), 1)]; % Blue for Array 3
swarmchart(xswarm,yswarm,10,colorswarm,"filled",'o','XJitterWidth',0.5)
ylabel('D/B power during full cycle')
ylim([0 300])
set(gca, 'XTick', 1:numel(dbStimData), 'XTickLabel', groupNames);
title('D/B power during full cycles - one night')

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
% fprintf('Kruskal-Wallis test p-value: %.4f\n', pKruskal);
% annotation('textbox', [0.1, 0.8, 0.2, 0.1], 'String', ...
%     sprintf('Kruskal-Wallis test p-value: %.4f\n', pKruskal), 'EdgeColor', 'none', 'HorizontalAlignment', ...
%     'right', 'VerticalAlignment', 'middle');


if pKruskal < 0.05
    disp('Significant differences found in Kruskal-Wallis test. Proceeding with pairwise comparisons...');

    if isstring(groupLabels)
        groupLabels = cellstr(groupLabels);
    end

    % Get unique group names
    [groupNames, ~, groupIdx] = unique(groupLabels);
    numGroups = numel(groupNames);

    % Initialize
    comparisons = {};
    raw_pvals = [];
    idx = 1;

    % Loop through all group pairs
    for i = 1:numGroups-1
        for j = i+1:numGroups
            % Extract data for group i and j
            data_i = allData(groupIdx == i);
            data_j = allData(groupIdx == j);

            % Wilcoxon rank-sum (Mann-Whitney U)
            [p, ~] = ranksum(data_i, data_j);

            % Store results
            comparisons{idx,1} = [groupNames{i} ' vs ' groupNames{j}];
            raw_pvals(idx,1) = p;
            idx = idx + 1;
        end
    end

    % Bonferroni correction
    corrected_pvals = min(raw_pvals * length(raw_pvals), 1);

    % Display results
    fprintf('\nPairwise Wilcoxon Rank-Sum Test (Mann-Whitney U):\n');
    for i = 1:length(raw_pvals)
        fprintf('%s:\t raw p = %.4f,\t Bonferroni-corrected p = %.4f\n', ...
            comparisons{i}, raw_pvals(i), corrected_pvals(i));
    end
    
else
    disp('No significant differences found in Kruskal-Wallis test.');
end


%savefigure
set(fdbDec,'PaperPositionMode','auto');
fileName=[SA.currentPlotFolder filesep 'DBfullcycpreStimPost'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);
set(fdbDec,'PaperPosition',[1,1,2.7,2.3]);
fileName=[analysisFolder filesep 'DBfullcycOneNight'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);


%% plot all nights: D/B full cycle - only red

type = 'Red';
wavelength = '635';
curTrials = contains(stimTable.Remarks,wavelength)...
        & ~contains(stimTable.Remarks,'Ex') ...
        & ~any(isnan(stimTable.dbMeans),2); %& contains(stimTable.Animal,curAni);
n = sum(curTrials);
N = length(unique(stimTable.Animal(curTrials)));
groupNames = {'Pre', 'During', 'After'};
curData = stimTable.dbMeans(curTrials,:);


% statistics:

% Assuming data in columns where each row is a subject and each column is a timepoint
% bdSWDataStat = stimTable.dbSWMeans(curTrials,:);
% cleanDBSW = bdSWDataStat(~any(isnan(bdSWDataStat), 2), :);
% n=height(cleanDBSW);

[p, tbl, stats] = friedman(curData, 1); % Here, 1 indicates within-subjects design
fprintf('p-value for freidman ANOVA test: %.5f\n',p)
% p-valure is very low, post hoc:
% Example data for three groups
before = stimTable.dbMeans(curTrials,1);
during = stimTable.dbMeans(curTrials,2);
after = stimTable.dbMeans(curTrials,3);

% Bonferroni-corrected alpha level
alpha = 0.05 / 3;

% Pairwise Wilcoxon signed-rank tests
[p_before_during, ~, stats_before_during] = signrank(before, during);
[p_during_after, ~, stats_during_after] = signrank(during, after);
[p_after_before, ~, stats_after_before] = signrank(after, before);
raw_pvals = [p_before_during,p_during_after,p_after_before];
num_comparisons = 3;

corrected_pvals_bonferroni = min(raw_pvals * num_comparisons, 1);
fprintf('Wilcoxon signed-rank test results with Bonferroni correction:\n');
fprintf('Before vs During: p-value = %.4f \n', corrected_pvals_bonferroni(1));
fprintf('During vs After: p-value = %.4f\n', corrected_pvals_bonferroni(2));
fprintf('After vs Before: p-value = %.4f\n ', corrected_pvals_bonferroni(3));
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
plot(mean(stimTable.dbMeans(curTrials,:),1,'omitnan'),'Color','k','LineWidth',2,'Marker','.','MarkerSize',10)
xlim([0.5, 3.5])
% grid on;
ylim([0 350])
xticks(1:3)
xticklabels(groupNames)
ylabel('D/B means during full cycle')

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
set(fdb,'PaperPosition',[1 1 2.7 2.3]);
fileName=[analysisFolder filesep 'DBchangeFullCycle'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);
