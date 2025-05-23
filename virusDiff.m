%% virus Figure
% plot a figure that check the differences in reactivnes to blue light.
% animals with virus: PV161, 149
% animals with no virus: 157 159 162 126
%% Sup 1F
viAnimals = ["161","149"];
curType = "47";

f=figure;

%plot the data
% plot nights without :
nonViTrials = contains(stimTable.Remarks,curType) &...
    ~contains(stimTable.Remarks,'Ex') & ...
    any(~contains(stimTable.Animal,viAnimals),2);
nNon = sum(nonViTrials);
Nnon= length(unique(stimTable.Animal(nonViTrials)));
ViTrials = contains(stimTable.Remarks,curType) &...
    ~contains(stimTable.Remarks,'Ex') & ...
    any(contains(stimTable.Animal,viAnimals),2);
nVi = sum(ViTrials);
Nvi = length(unique(stimTable.Animal(ViTrials)));
 
% Statistics:
[pMW_sham,h_sham] = ranksum(stimTable.dbDiffShamM(ViTrials), stimTable.dbDiffShamM(nonViTrials));
    fprintf('Mann-Whitney Signed-Rank Test p-value - sham: %.4f\n', pMW_sham);
    fprintf('h = %i\n',h_sham)

[pMW_stim, h_stim] = ranksum(stimTable.dbDiffStimM(ViTrials), stimTable.dbDiffStimM(nonViTrials));
    fprintf('Mann-Whitney Signed-Rank Test p-value - stim: %.4f\n', pMW_stim);
    fprintf('h = %i\n',h_stim)

%plot
x = 1:2;
h1 = plot(x,[stimTable.dbDiffShamM(nonViTrials) stimTable.dbDiffStimM(nonViTrials)], ...
    'Marker','.','LineStyle','-','Color',[1 0.65 0]);
hold on
h2 = plot(x,[stimTable.dbDiffShamM(ViTrials) stimTable.dbDiffStimM(ViTrials)], ...
    'Marker','.','LineStyle','-','Color','g');
xlim([0.9 2.1])
xticks([x]);xticklabels(["Sham","Stim"])

legend([h1(1), h2(1)],"Animals not injected","Animals injected with ChR2")


% save figue

set(f,'PaperPosition',[1 1 2 3]);
fileName=[analysisFolder filesep 'virusChangeStimShamBlue'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);



%% Traces of stimulations with virus:

i = 6;
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
SA.setCurrentRecording(recName);
DB = SA.getDelta2BetaRatio;

% stimulations timings:
t_ch = stimTable.StimTrighCh(i);
T=SA.getDigitalTriggers;
stims = T.tTrig{t_ch};
stimStartT = stimTable.stimStartT(i);
stimEndT = stimTable.stimEndT(i);
firstTrig=stims(1:8:end-2);
endStim=stims(8:8:end)+200;
trial = reshape(stims,[8,length(stims)/8])';

pre=20000;
post=100000;
ch = 17;
% get the raw data and LFP:
% for j = 10:40
 j = 39;
 VipTmp=find(DB.t_ms>(firstTrig(j)-pre) & DB.t_ms<(firstTrig(j)+post));
 VidbTmp=DB.bufferedDelta2BetaRatio(VipTmp);
 [vlfp,vlfp_t] = SA.currentDataObj.getData(ch,firstTrig(j)-pre,pre+post);
    

  f = figure;
    h1 = subplot(2,1,1);
    % set(f, 'Position', [100, 100, 1200, 400]);
    % yyaxis left
    plot(vlfp_t/1000,squeeze(vlfp),'k'); hold on;
    curstims = trial(j,:) -trial(j,1) +pre ;
    xline(curstims/1000,'r','LineWidth',1.5)
    % yyaxis right
    plot(VidbTmp,'Color','b','LineWidth',2)
    title('With virus');
    % sgtitle(sprintf('Trial num: %i',j))
    hold off;
% end

% Traces of stimulations without virus:

i = 9;
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
SA.setCurrentRecording(recName);
DB = SA.getDelta2BetaRatio;

% stimulations timings:
t_ch = stimTable.StimTrighCh(i);
T=SA.getDigitalTriggers;
stims = T.tTrig{t_ch};
stimStartT = stimTable.stimStartT(i);
stimEndT = stimTable.stimEndT(i);
firstTrig=stims(1:8:end-2);
endStim=stims(8:8:end)+400;
trial = reshape(stims,[8,length(stims)/8])';

pre=20000;
post=100000;
ch = 17;

% get the raw data and LFP:
% for j = 10:40
 j = 28;
 pTmp=find(DB.t_ms>(firstTrig(j)-pre) & DB.t_ms<(firstTrig(j)+post));
 dbTmp=DB.bufferedDelta2BetaRatio(pTmp);
 [lfp,lfp_t] = SA.currentDataObj.getData(ch,firstTrig(j)-pre,pre+post);
    
  % f = figure;
    h2 = subplot(2,1,2);
    % set(f, 'Position', [100, 100, 1200, 400]);
    % yyaxis left
    plot(lfp_t/1000,squeeze(lfp),'k'); hold on;
    curstims = trial(j,:) -trial(j,1) +pre ;
    xline(curstims/1000,'r','LineWidth',1.5)
    % yyaxis right
    plot(dbTmp,'Color','b','LineWidth',2)
    title('Without Virus');
    % sgtitle(sprintf('W/O Trial num:%i',j))
    hold off;
%end


%save figue
set(f,'PaperPosition', [1 1 4 2.5]);
fileName=[analysisFolder filesep 'tracesVirusNovirus'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);

%% Third eye:
% plot a figure that shows the stim sham of the 
curTrialsThird = stimTable.Mani ==2;
animalNum = numel(unique(stimTable.Animal(curTrialsThird)));
stimTypes = unique(stimTable.Remarks(curTrialsThird));
recNum = sum(curTrialsThird);


%% plot the bar plot - D/B change - stim sham -  all colors
animals = unique(stimTable.Animal);
stimType = ["Green","Red","LED"];
stimWaveL = ["532","635","LED"];
% plotColors = {[0 0.586 0.9766],[0.05 0.81 0.379],[1 0.27 0.27], [0.5 0.5 0.5]};
numAnimal = length(animals);
numType = length(stimType);
statsStimSham = struct();

% mean normelized change in D/B.
f=figure;
sgtitle('mean normelized D/B')

for type = 1:numType
    h = subplot(1,numType,type);
    %plot the data
    %curAni = animals{animal};
    curType = stimWaveL(type);
    curName = stimType(type);
    curTrialsThird = contains(stimTable.Remarks,curType) &...
                     ~contains(stimTable.Remarks,'Ex') & ...
                     all(~isnan(stimTable.dbDiffStimM),2) & ...
                     all(~isnan(stimTable.dbDiffShamM),2) & ...
                     stimTable.Mani==1;
    nThird = sum(curTrialsThird);
    NThird = length(unique(stimTable.Animal(curTrialsThird)));
    curDataT = [stimTable.dbDiffShamM(curTrialsThird) stimTable.dbDiffStimM(curTrialsThird)];

    curTrialsNoThird = contains(stimTable.Remarks,curType) &...
                     ~contains(stimTable.Remarks,'Ex') & ...
                     all(~isnan(stimTable.dbDiffStimM),2) & ...
                     all(~isnan(stimTable.dbDiffShamM),2) & ...
                     stimTable.Mani==2;
    nNoThird = sum(curTrialsNoThird);
    NNoThird = numel(unique(stimTable.Animal(curTrialsNoThird)));
    curDataNT = [stimTable.dbDiffShamM(curTrialsNoThird) stimTable.dbDiffStimM(curTrialsNoThird)];

    % curCol = plotColors{type};

    % curData = [stimTable.dbDiffShamM(curTrialsThird), stimTable.dbDiffStimM(curTrialsThird)];
    % %statistics:
    % % 1. Wilcoxon Signed-Rank Test
    % [pWilcoxon, ~, statsWilcoxon] = signrank(curData(:,1), curData(:,2));
    % fprintf('Wilcoxon Signed-Rank Test p-value: %.4f\n', pWilcoxon);
    % disp(statsWilcoxon)

    if n>0
        x=[1,2];
        [~, animalIndices] = ismember(stimTable.Animal(curTrialsThird), uniqueAnimals);
        curColorMatT = animalsColors(animalIndices, :);
        [~, animalIndices] = ismember(stimTable.Animal(curTrialsNoThird), uniqueAnimals);
        curColorMatNT = animalsColors(animalIndices, :);
        hold on
        for j = 1:height(curDataT)
            plot(x,curDataT(j,:),'Color',curColorMatT(j,:),'Marker','.')
        end
        hold on
        for j = 1:height(curDataNT)
            plot(x,curDataNT(j,:),'Color',curColorMatNT(j,:),'LineStyle','--','Marker','.')

        end
        % plot(x,[curMeanNdbSham,curMeanNdbStim],'Color','k','LineWidth',3,'Marker','.')
       
        % 
        % %plot Ex
        % curExI = contains(stimTable.Remarks,curType) &...
        %         contains(stimTable.Remarks,'Ex') & ...
        %         all(~isnan(stimTable.dbDiffStimM),2) & ...
        %         all(~isnan(stimTable.dbDiffShamM),2);
        % curEx = [stimTable.dbDiffShamM(curExI) stimTable.dbDiffStimM(curExI)];
        % if ~isempty(curEx )
        %     plot(x,curEx,'Color',[0.3 0.3 0.3],'LineStyle','--','Marker','.')
        % end
        hold off
    end
    xticks([1, 2]); % Position of the x-ticks
    xticklabels({'Sham', 'Stim'}); % Labels for the x-ticks
    xlim([0.5, 2.5]);
    annotation('textbox', [.095 + 0.195*type, 0.85, 0.03, 0.1], 'String', ...
        sprintf('n=%i,N=%i',n,N), 'EdgeColor', 'none', 'HorizontalAlignment', ...
        'right', 'VerticalAlignment', 'middle');
   annotation('textbox', [.15 + 0.195*type, 0.85, 0.03, 0.1], 'String', ...
    sprintf('Wilcoxon test p =%.3f',pWilcoxon), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');
    
    ylim([-55 160])

    if n==0
        plot(0,0)
    end

    % add titles. labels...
    ylabel('D2B power')
    title(stimType(type))
 
end

% savefigure
set(f,'PaperPositionMode','auto');
fileName=[analysisFolder filesep 'meanNormBDStimShamThirdEye'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);


% clearvars -except stimTable SA analysisFolder
