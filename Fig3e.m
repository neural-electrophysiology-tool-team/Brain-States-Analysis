%% fig3F
load([analysisFolder filesep 'LMData.mat'])
% plot head movements - All nights!

stimType = ["Blue","Green","Red","LED"];
stimWaveL = ["47","532","635","LED"];
numType = length(stimType);
Groups = ["Wake","Pre","Stim"];
% statsHeadMovAll = struct();

diffStimPre = [];
groupNum = [];
colorMat = [];
diffStimPreM =[];
ns = [];
Ns = [];
% Head movements:
fLMdiff=figure;
for type = 1:numType
    subplot(2,1,1)
    %plot the data
    %curAni = animals{animal};
    curName = stimType(type);
    curType = stimWaveL(type);
    curTrials = contains(stimTable.Remarks,curType) & ...
             ~contains(stimTable.Remarks,'Ex') &...
             all(~isnan(LMData.LMpre),2); 
    n = sum(curTrials);
    N = length(unique(stimTable.Animal(curTrials)));
    ns = [ns; n];
    Ns = [Ns; N];
    LMpre = LMData.LMpre(curTrials, :);
    LMwake = LMData.LMwake(curTrials, :);
    LMstimbinM = cell2mat(LMData.LMstimbin); % takes out the nan val
    LMstimbintrialM = mean(LMstimbinM(curTrials),2); % mean for each night
    
    curDiffStimPre = LMstimbintrialM-LMpre;
    diffStimPre = [diffStimPre;curDiffStimPre];
    diffStimPreM = [diffStimPreM; mean(curDiffStimPre)];
    groupNum = [groupNum; repmat(type,length(curDiffStimPre),1)];
    curDiffPreWake = LMpre;
    curDdiffStimWake = LMstimbintrialM;
    curData = [curDiffPreWake,curDdiffStimWake];
    % plot normelized to wake:

    [~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
    curColorMat = animalsColors(animalIndices, :); 
    colorMat = [colorMat; curColorMat];    

    hold on
    for i = 1:length(curData)
        plot([type*2-1, type*2],curData(i,:),'Color',curColorMat(i,:), 'Marker','.','LineWidth',0.5)
    end
    plot([type*2-1, type*2],mean(curData),'Color','k', 'Marker','.','LineWidth',1)

   
    % check the statistics:
    % Assuming data in columns where each row is a subject and each column is a timepoint
  
    % Pairwise Wilcoxon signed-rank tests

    p = signrank(curData(:,1), curData(:,2));
    fprintf('%s: p-value = %.5f\n',stimType(type), p);
   % 
   %  annotation('textbox', [.095 + 0.195*type, 0.85, 0.03, 0.1], 'String', ...
   %      sprintf('n=%i,N=%i',n,N), 'EdgeColor', 'none', 'HorizontalAlignment', ...
   %      'right', 'VerticalAlignment', 'middle');
   % % ylim([-40 160])
  
    % add titles. labels...
    if type==1
        ylabel('Head Movements')
    end
end
xlim([0.7 8.2])
xticks(1:8);xticklabels(repmat(["Pre","Stim"],1,numType))

% plot only stim-pre diff:
subplot(2,1,2)
swarmchart(groupNum,diffStimPre,15,colorMat,'filled','XJitterWidth',0.5);
hold on;
scatter(1:4,diffStimPreM,15,'k','Marker','+')
xticks(1:4), xticklabels(stimType); xlim([0.7 4.3]); %ylim([0 10]);
yline(0,'--','Color',[0.4 0.4 0.4])
[pKW, tbl, stats] = kruskalwallis(diffStimPre,groupNum,'off');

% if pKW < 0.05
%     c = multcompare(stats, 'CType', 'dunn-sidak','off');
% end

psFromZero = zeros(1,numType);
for i= 1:numType
    curData = diffStimPre(groupNum==i);
    [pfromZero, h] = signrank(curData, 0);
    psFromZero(1,i) = pfromZero; 

end

% savefigure
set(fLMdiff,'PaperPosition',[1 5 3.5 4]);
fileName=[analysisFolder filesep 'LMdiffs'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);
