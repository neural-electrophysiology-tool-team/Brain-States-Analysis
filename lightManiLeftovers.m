 % light manipulation left overs
 %stuff i dont need or use for the paper:



%% get the stim sham for single rec
SA.setCurrentRecording('Animal=PV162,recNames=Night23');
getStimSham(SA,11);
SA.getDelta2BetaRatio('tStart', 0 ,'win',1000*60*60*11);
SA.getDelta2BetaAC;
SA.plotDelta2BetaRatio
SA.plotDelta2BetaAC;
SA.plotDelta2BetaSlidingAC;
plotStimSham(SA)
 %% Check for spikes in recording:
stimTable.spikes = zeros(height(stimTable),1);

for i = 1:height(stimTable)
    recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
    SA.setCurrentRecording(recName);
    timeSeriesViewer(SA.currentDataObj)
    spikes = 0;
    % Some computations or operations
    disp(['Rec: ', recName])

    % Pause for manual input
    new_value = input('Does this rec have spikes? (or press Enter for no): ', 's');
    
    % If the user provides input, update the parameter
    if ~isempty(new_value)
        spikes = str2double(new_value);  % Convert to number
    end
    stimTable.spikes(i) = spikes;
    
end

%% plot imagesc spike rate
% figure;
% imagesc(spikeRateM);
% colormap(flipud(gray))  % Set colormap to grayscale
% colorbar;        % Optional: show colorbar
% 
% xlabel('Time[100 ms]')
% ylabel('Trial #')
% % Hold the plot to add lines
% hold on;
% 
% % Define X positions for the red lines (customize as needed)
% stimDiff = mean(mean(diff(trial,[],2)));
% xPositions = (pre-spikeRateT(1) + (0:7)*(stimDiff))/100;  % Example positions for lines
% 
% % Add red vertical lines
% for i = 1:length(xPositions)
%     line([xPositions(i), xPositions(i)], [0.5, size(spikeRateM, 1)+1], 'Color', 'red', 'LineWidth', 0.5);
% end
% 
% hold off;
% title('spikeRate (spike per sec) average per stim trial, PV161N18')
% % 
% saveas (gcf, [analysisFolder filesep 'spikeRateAllTrialsPV161N18.pdf']);


%% create rec list for batch analysis:
SA=sleepAnalysis('/media/sil1/Data/Pogona Vitticeps/brainStatesWake.xlsx'); 
% animals = unique(stimTable.Animal);                        
recList = {};
for i = 1:height(stimTable)
    animal = stimTable.Animal{i};
    recName = stimTable.recNames{i};
    curRec=['Animal=' animal ',recNames=' recName];
    recList{end+1} = curRec;
    
end
% disp(recList);

% run batch analysis on all nights
SA.batchProcessData('getLizardMovements',recList)
%% plot all stim sham, SlidingAc and D2B for all records (indevidually)

load([analysisFolder filesep 'stimTable.mat'])
for i = 1:height(stimTable)
    % set te current rec:
    recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
    SA.setCurrentRecording(recName);
    SA.plotDelta2BetaRatio;
    SA.plotDelta2BetaSlidingAC;
    getStimSham(SA,stimTable.StimTrighCh(i),1);
    plotStimSham(SA);
%     keyboard;
end
%% plot D2B general decrease for each part of the stimulation - peaks analysis:
% get and orgenize the data: 

% calculate the D2B average for each 30 sec for each part of the stimulation
p = 30*60*1000; %some time diff for the cycle to change.
DBwin = 2*60*60*1000; % 2 hrs in ms

DBpeaksPre = zeros(height(stimTable),1);
DBpeaksStim = zeros(height(stimTable),1);
DBpeaksPost = zeros(height(stimTable),1);

% for i = 1:height(stimTable)
i =22 ;
    recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
    SA.setCurrentRecording(recName);
    DB = SA.getDelta2BetaRatio;

    % stimulations timings:
    t_ch = stimTable.StimTrighCh(i);
    T=SA.getDigitalTriggers;
    stimStartT = T.tTrig{t_ch}(1);
    stimEndT = T.tTrig{t_ch}(end);
    firstTrig=T.tTrig{t_ch}(1:8:end-2);
    endStim=T.tTrig{t_ch}(8:8:end)+400;

   % get the DB 
    preind = find(DB.t_ms < stimStartT & DB.t_ms > (stimStartT-DBwin)) ;
    DBpre = DB.bufferedDelta2BetaRatio(preind);
    stimind = find(DB.t_ms > stimStartT+p & DB.t_ms < (stimStartT+ p+DBwin)) ;
    DBstim = DB.bufferedDelta2BetaRatio(stimind);
    postind = find(DB.t_ms > stimEndT+p & DB.t_ms < (stimEndT +p + DBwin)) ;
    DBpost = DB.bufferedDelta2BetaRatio(postind);

    % peak analysis:
    [peakPre, peakLocPre] = findpeaks(DBpre, 'MinPeakHeight', 100);
    [peakstim, peakLocStim] = findpeaks(DBstim, 'MinPeakHeight', 100);
    [peakPost, peakLocPost] = findpeaks(DBpost, 'MinPeakHeight', 100);

 

%% Create a figure for the plot - violin
figure; hold on;
vioData = {peakPre, peakstim, peakPost };
% Define colors for each group
colors = [0.2, 0.6, 0.8; 0.9, 0.4, 0.3; 0.5, 0.8, 0.5];

violin(vioData,'xlable',groupNames, 'facecolor',colors,'edgecolor',[0.6 0.6 0.6],'medc',[])

% Overlay individual data points
jitterAmount = 0.5;  % Amount to jitter the points along the x-axis
for i = 1:numel(vioData)
    % Generate x-coordinates with jitter for each group
    x = i + (rand(size(vioData{i})) - 0.5) * jitterAmount;
    
    % Plot data points for this group
    scatter(x, vioData{i}, 20, colors(i,:), 'filled', 'MarkerFaceAlpha', 0.5); % Adjust size and transparency
end

% Customize the plot
xlim([0.5, numel(vioData) + 0.5]);  % Adjust x-axis limits
set(gca, 'XTick', 1:numel(vioData), 'XTickLabel', groupNames);
ylabel('D/B peak');
title('D/B peaks changes during stimulations');

hold off;
% savefigure
set(gcf,'PaperPosition',[.25 3 8 6])
saveas (gcf, [analysisFolder filesep 'DBpeaksrec22.pdf']);


%% plot all nights:
figure;

%plot using plot - to show which of the dots are the same night. 
peaksAllrecs = [DBpeaksPre,DBpeaksStim, DBpeaksPost];
means = [mean(DBpeaksPre,'omitnan'),mean(DBpeaksStim,'omitnan'),mean(DBpeaksPost,'omitnan')];
plot(peaksAllrecs','-o','Color',[0.5 0.5 0.5]);
hold on
plot(means,'-o', 'Color','black', 'LineWidth',3)

hold off;
% Add labels and title

xlim([0.5, 3.5])
xticks([1 2 3]);  % Set x-axis ticks to 1, 2, 3
xticklabels({'Pre', 'stimulations', 'Post'});  % Group labels
ylabel('Peak D/B');
title('Scatter Plot of Peak D/B Values - All nights');
grid on;  % Optional: Add grid lines
% savefigure
set(gcf,'PaperPosition',[.25 3 8 6])
saveas (gcf, [analysisFolder filesep 'DBpeaksAllRecs.pdf']);
clearvars -except stimTable SA analysisFolder
%% plot AC - Basic
ACcomPer = [ACprePer ACstimPer ACpostPer];
ACcomP2V = [ACpreP2V ACstimP2V ACpostP2V];
stimTable.ACcomPer = ACcomPer;
stimTable.ACcomP2V = ACcomP2V;

figure;
plot(ACcomPer'/1000, '-o')
xlim([0.5 3.5])
figure; plot(ACcomP2V','-o')


%% PLOT stim activation according to stim type
% assuming the table is in the workspace
% load([analysisFolder filesep 'stimTable.mat'])

% loop on table
f = figure;
set(f, 'Position', [100, 100, 1200, 800]);

times = stimTable.times{1};
pre=50000;
% STIM DUTRATION NEEDS A THINK!!!!!!!
mStimDur = mean(stimTable.stimDuration(3:end));
post=100000;

%plot the blues:
ax1 = subplot(2,2,1);
hold on
blues = contains(stimTable.Remarks,'47');
meanBlue = mean(cell2mat(stimTable.StimAvg(blues)),1,'omitnan');
meanSham = mean(cell2mat(stimTable.StimAvgSham(blues)),1,'omitnan');
for i = 1:height(stimTable)
    if blues(i) && stimTable.Mani(i) == 1
        plot(times,stimTable.StimAvg{i},'blue')
    end 
    if blues(i) && stimTable.Mani(i) == 2
        plot(times,stimTable.StimAvg{i},'blue','LineStyle','--')
    end

end
plot(times,meanBlue,'color','blue','LineWidth',4)
plot(times,meanSham,'color','black','LineWidth',2)
xline(pre/1000, 'LineWidth',1,'Color','black')
xline((pre+mStimDur)/1000,'LineWidth',1,'Color','black')
title('blue laser stimulation');
xlabel('Time[S]'), ylabel('D2B power')
hold off

%plot the greens:
ax2 = subplot(2,2,2);
hold on
greens = contains(stimTable.Remarks,'532');
meangreen = mean(cell2mat(stimTable.StimAvg(greens)),1,'omitnan');


for i = 1:height(stimTable)
    if greens(i) && stimTable.Mani(i) == 1
        plot(times,stimTable.StimAvg{i},'green')
    end 
    if greens(i) && stimTable.Mani(i) == 2
        plot(times,stimTable.StimAvg{i},'green','LineStyle','--')
    end

end
meanSham = mean(cell2mat(stimTable.StimAvgSham(greens)),1,'omitnan');
plot(times,meanSham,'color','black','LineWidth',2)
plot(times,meangreen,'color','green','LineWidth',4)
xline(pre/1000, 'LineWidth',1,'Color','black')
xline((pre+mStimDur)/1000,'LineWidth',1,'Color','black')
title('Green laser stimulation');
xlabel('Time[S]'), ylabel('D2B power')
hold off


% plot the reds
ax3 = subplot(2,2,3);
hold on
reds = contains(stimTable.Remarks,'635');
meanreds = mean(cell2mat(stimTable.StimAvg(reds)),1,'omitnan');
for i = 1:height(stimTable)
    if reds(i) && stimTable.Mani(i) == 1
        plot(times,stimTable.StimAvg{i},'r')
    end 
    if reds(i) && stimTable.Mani(i) == 2
        plot(times,stimTable.StimAvg{i},'r','LineStyle','--')
    end

end
meanSham = mean(cell2mat(stimTable.StimAvgSham(reds)),1,'omitnan');
plot(times,meanSham,'color','black','LineWidth',2)
plot(times,meanreds,'color','r','LineWidth',4)
xline(pre/1000, 'LineWidth',1,'Color','black')
xline((pre+mStimDur)/1000,'LineWidth',1,'Color','black')
title('Red laser stimulation');

hold off
xlabel('Time[S]'), ylabel('D2B power')


% plot the LEDs
ax4 = subplot(2,2,4);
hold on
leds = contains(stimTable.Remarks,'LED');
meanleds = mean(cell2mat(stimTable.StimAvg(leds)),1,'omitnan');
gray = [0.5 0.5 0.5];
for i = 1:height(stimTable)
    if leds(i) && stimTable.Mani(i) == 1
        plot(times,stimTable.StimAvg{i},'color',gray)
    end 
    if leds(i) && stimTable.Mani(i) == 2
        plot(times,stimTable.StimAvg{i},'color',gray,'LineStyle','--')
    end

end
meanSham = mean(cell2mat(stimTable.StimAvgSham(leds)),1,'omitnan');
plot(times,meanSham,'color','black','LineWidth',2)
plot(times,meanleds,'color',gray,'LineWidth',4)
xline(pre/1000, 'LineWidth',1,'Color','black')
xline((pre+mStimDur)/1000,'LineWidth',1,'Color','black')
title('LED stimulation');

hold off


% general proprties
linkaxes([ax1,ax2,ax3,ax4],'xy')
ylims([0 200])
xlabel('Time[S]'), ylabel('D2B power')

% savefigure
saveas (gcf, [analysisFolder filesep 'all_colors_sham.jpg']);
% clearvars -except stimTable SA analysisFolder
%% figure according to time:

f=figure;
set(f, 'Position', [100, 100, 1200, 400]);
sgtitle('mean normelized D/B')
% stimTable.nightn    
colors=jet(max(stimTable.nightnum));

for type = 1:numType
    h = subplot(1,numType,type);
    %plot the data
    %curAni = animals{animal};
    curType = stimWaveL(type);
    curTrials = contains(stimTable.Remarks,curType); %& contains(stimTable.Animal,curAni);
    n = sum(curTrials);
    N = length(unique(stimTable.Animal(curTrials)));
    curCol = plotColors{type};
    curMeanNdbStim = mean(stimTable.dbDiffStimM(curTrials),1,'omitnan');
    curMeanNdbSham = mean(stimTable.dbDiffShamM(curTrials),1,'omitnan');
    if n>0
        x=[1,2];
        plot(x,[curMeanNdbSham,curMeanNdbStim],'-o','color','black','LineWidth',3)
        hold on
        for j=1:height(stimTable)
            if curTrials(j)==1
                plot(x,[stimTable.dbDiffShamM(j), stimTable.dbDiffStimM(j)] ...
                ,'-o','Color',colors(stimTable.nightnum(j),:))
            end
        end
        hold off
    end
    

    xticks([1, 2]); % Position of the x-ticks
    xticklabels({'Sham', 'Stim'}); % Labels for the x-ticks
    xlim([0.5, 2.5]);
    annotation('textbox', [.095 + 0.195*type, 0.85, 0.03, 0.1], 'String', ...
        sprintf('n=%i,N=%i',n,N), 'EdgeColor', 'none', 'HorizontalAlignment', ...
        'right', 'VerticalAlignment', 'middle');
   ylim([-40 160])
    if n==0
        plot(0,0)
    end

    % add titles. labels...
    ylabel('D2B power')
    title(stimType(type))
 
end
c=colorbar;    
clim([min(stimTable.nightnum), max(stimTable.nightnum)]);

% Define custom tick positions for the colorbar
c.Ticks = linspace(min(stimTable.nightnum), max(stimTable.nightnum), 6);

c.Label.String = 'Night Number';  % Add a descriptive label
c.TickLabels = arrayfun(@num2str, round(linspace(5, 40, 6)), ...
    'UniformOutput', false);  % Ticks between 5 and 40

% savefigure
set(gcf,'PaperPosition',[.25 3 8 6])
saveas (gcf, [analysisFolder filesep 'meanNormBDStimShamTimeColor.pdf']);

clearvars -except stimTable SA analysisFolder LM



% 
% % stimulations timings:
% t_ch = stimTable.StimTrighCh(i);
% T=SA.getDigitalTriggers;
% stimStartT = T.tTrig{t_ch}(1);
% stimEndT = T.tTrig{t_ch}(end);
% firstTrig=T.tTrig{t_ch}(1:8:end-2);
% endStim=T.tTrig{t_ch}(8:8:end)+400;
% stimDuration=(endStim(1)-firstTrig(1));
% %% plot the movement binned to the whole night. 
% LM_DBt = stimTable.LM_DBt{i};
% figure; 
% plot(DB.t_ms/1000, LM_DBt); xlabel('Time[s]'), ylabel('MovCount');
% % save fig
% set(gcf,'PaperPosition',[.25 3 8 6])
% saveas (gcf, [analysisFolder filesep 'wholeNightmove.pdf']);
% 
% %% calculate the mean of events in the hour before stim, during the stims,
% % and an hour after. 
% befind = find(DB.t_ms>stimStartT-1000-(60*60*1000)&DB.t_ms<stimStartT-1000);
% befMov = mean(LM_DBt(befind));
% wakeMov = mean(LM_DBt(1:3600));
% postind = find(DB.t_ms>stimEndT+(60*60*1000)&DB.t_ms<stimEndT+(2*60*60*1000));
% aftMov = mean(LM_DBt(postind));
% stimlength = 150;
% stimLM = zeros(numel(firstTrig),stimlength);
% post =150*1000;
% binSize = 10;
% stimLMbin = zeros(numel(firstTrig),stimlength/binSize);
% % Add the last bin edge 
% for i=1:numel(firstTrig)
%     pTmp=find(DB.t_ms>(firstTrig(i)) & DB.t_ms<=(firstTrig(i)+post));
%         %StimDB(i,:)=1./DB.bufferedDelta2BetaRatio(pTmp);
%         if length(pTmp) ~=stimlength
%             pTmp = pTmp(1:stimlength);
%         end
%         curstimLM = LM_DBt(pTmp);
%         stimLM(i,:) = curstimLM;
% 
%         % Reshape the array to 15x10
%         curReshLM = reshape(curstimLM, binSize,stimlength/binSize);
% 
%         % Sum along the rows (dim=1) to get a 1x15 array
%         curStimLMBin = sum(curReshLM, 1);
%         stimLMbin(i,:) = curStimLMBin/binSize;
% 
% end
% stimLMmean = mean(stimLM,1);
% stimLMbinmean = mean(stimLMbin,1);
% 
% %% plot the data from the stimulation time:
% figure; subplot (2,1,1)
% plot(stimLM')
% hold on; plot(stimLMmean,'black','LineWidth',3); 
% ylabel('Event count'), xlabel('Time [s]'); xline(0,'color','r');xline(38,'color','r');
% hold off
% 
% subplot (2,1,2);plot(stimLMbin')
% hold on;
% plot(stimLMbinmean,'black','LineWidth',3); 
% ylabel('Event count'), xlabel('Time [10s]'); xline(0,'color','r');xline(3.8,'color','r');
% hold off
% % save fig
% set(gcf,'PaperPosition',[.25 3 8 6])
% saveas (gcf, [analysisFolder filesep 'movDurStimsinglenight.pdf']);

%% D/B stim sham old calculation 
f=figure;
% set(f, 'Position', [100, 100, 1200, 400]);
% sgtitle('Max average D/B')
% for type = 1:numType
%     h = subplot(1,numType,type);
%     %plot the data
%     %curAni = animals{animal};
%     curType = stimWaveL(type);
%     curTrials = contains(stimTable.Remarks,curType); %& contains(stimTable.Animal,curAni);
%     n = sum(curTrials);
%     N = length(unique(stimTable.Animal(curTrials)));
%     curCol = plotColors{type};
%     curMeanStim = mean(stimTable.maxStim(curTrials),1,'omitnan');
%     curMeanSham = mean(stimTable.maxSham(curTrials),1,'omitnan');
%     if n>0
%         x=[1,2];
%         plot(x,[curMeanSham,curMeanStim],'-o','color',curCol,'LineWidth',3)
%         hold on
%         plot(x,[stimTable.maxSham(curTrials), stimTable.maxStim(curTrials)] ...
%             ,'-o','Color',curCol)
%           hold off
%     end
%     xticks([1, 2]); % Position of the x-ticks
%     xticklabels({'Sham', 'Stim'}); % Labels for the x-ticks
%     xlim([0.5, 2.5]);
%     annotation('textbox', [.095 + 0.195*type, 0.85, 0.03, 0.1], 'String', ...
%         sprintf('n=%i,N=%i',n,N), 'EdgeColor', 'none', 'HorizontalAlignment', ...
%         'right', 'VerticalAlignment', 'middle');
% 
%     if n==0
%         plot(0,0)
%     end
% 
%     % add titles. labels...
%     ylabel('D2B power')
%     title(stimType(type))
%     ylim([0 450])
% end
% 
% % savefigure
% set(gcf,'PaperPosition',[.25 3 8 6])
% saveas (gcf, [analysisFolder filesep 'maxStimShamAll1.pdf']);

%%
	
		



% look at the angeles changesin the videos:
% check 5 recs:
% 1. PV161, night18 (probably not the right properties of the headstage)
recName = 'Animal=PV161,recNames=Night18';
SA.setCurrentRecording(recName);
LM = SA.getLizardMovements;
figure; plot(LM.t_static_ms, LM.angles(1,:)); 
figure; plot(LM.t_static_ms, LM.angles(2,:)); 
figure; plot(LM.t_static_ms, LM.angles(3,:)); 
diff1 = 0.9765-0.976;
diff2 = 2.1942-2.195;
diff3 = 2.1995-2.198; %biggest diff. tis is the critical channel.
% small bump in the middle of the sleep - invisible in the recording. (back
% camera). the head lift in the end is seen in the graph. at the same
% timings.
 %% 2. PV149, night12
recName = 'Animal=PV149,recNames=Night12';
SA.setCurrentRecording(recName);
LM = SA.getLizardMovements;
figure; plot(LM.t_static_ms, LM.angles(1,:)); 
figure; plot(LM.t_static_ms, LM.angles(2,:)); 
figure; plot(LM.t_static_ms, LM.angles(3,:)); 
% diff1 = ;
% diff2 = 2.1942-2.195;
% diff3 = 2.1995-2.198; %biggest diff between wake and sleep times. 
% .
% small bump in the middle of the sleep - invisible in the recording. (back
% camera). the head lift in the end is seen in the graph. at the same
% timings.
%% 3. PV159, night10 (probably not the right properties of the headstage)
recName = 'Animal=PV159,recNames=Night10';
SA.setCurrentRecording(recName);
LM = SA.getLizardMovements;
figure; plot(LM.t_static_ms, LM.angles(1,:)); 
figure; plot(LM.t_static_ms, LM.angles(2,:)); 
figure; plot(LM.t_static_ms, LM.angles(3,:)); 
% diff1 =  extreamlt small
diff2 = 2.1945-2.1955;% same and small
diff3 = 2.1995-2.1985;
% i don't understand which of the axis is the "right" one. 

%% 4. PV161, night13 - run the LM anallysis again. right parameters for headstage 
recName = 'Animal=PV161,recNames=Night13';
SA.setCurrentRecording(recName);
SA.getLizardMovements('overwrite',1)
LM = SA.getLizardMovements;
figure; plot(LM.t_static_ms, LM.angles(1,:)); 
figure; plot(LM.t_static_ms, LM.angles(2,:)); 
figure; plot(LM.t_static_ms, LM.angles(3,:)); % this is the right angle!
%  =  extreamlt small

%% 5. PV157, night16 -seems like the right headstage. 
recName = 'Animal=PV157,recNames=Night16';
SA.setCurrentRecording(recName);
% SA.getLizardMovements('overwrite',1)
LM = SA.getLizardMovements;
figure; plot(LM.t_static_ms, LM.angles(1,:)); 
figure; plot(LM.t_static_ms, LM.angles(2,:)); 
figure; plot(LM.t_static_ms, LM.angles(3,:)); % this is the right angle

%% 6. PV159, night34 -seems like the right headstage. 
recName = 'Animal=PV159,recNames=Night34';
SA.setCurrentRecording(recName);
SA.getLizardMovements('overwrite',1)
LM = SA.getLizardMovements;
figure; plot(LM.t_static_ms, LM.angles(1,:)); 
figure; plot(LM.t_static_ms, LM.angles(2,:)); 
figure; plot(LM.t_static_ms, LM.angles(3,:)); % this is the right angle

%% might be useful to find the timings between stimulations:
ISIt = extract_spike_intervals(spikeRateT, xPositions);



%% PLOT ISI vs ITI: violin plot with swarm plot 
N = 4;
groupNames = ["Pre"; "Post"];
colors = [0.6 0.8 0.6; 1 0.7 0.4];

fall = figure; % all units:
subplot(1,2,1)
IIdata = [ITIsamples, ISIsamples];
x = [ones(length(ITIsamples),1) 2*ones(length(ITIsamples),1)];
swarmchart(x,IIdata,15,colors,'filled','XJitterWidth',0.2);
ylabel('Spikes/S')

subplot(1,2,2)
plot([1,1.2] ,IIdata,'Color',[0.5 0.5 0.5], 'Marker','.','MarkerSize',10);
hold on;
plot([1,1.2] ,mean(IIdata), 'Color','k','Marker','.', 'LineWidth',2);
xlim([0.7 1.5])
% xticks([1,1.2]);xticklabels(["Inter-Trial Interval", "Inter-Stim Interval"]);
ylabel('Spikes/S')
sgtitle('all units')
n = length(ITIsamples);
[pWilcoxon, ~, statsWilcoxon] = signrank(ITIsamples, ISIsamples);

annotation('textbox', [0.8, 0.85, 0.03, 0.1], 'String', ...
    sprintf('n=%i,N=%i',n,N), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');
annotation('textbox', [0.8, 0.7, 0.03, 0.1], 'String', ...
    sprintf('Wilcoxon  p-va: %.4f\n', pWilcoxon), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');
% savefigures
    set(fall,'PaperPositionMode','auto');
    fileName=[analysisFolder filesep 'ITIISIallunits'];
    print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);


    %% PLOT - Movement during stimulation - All Nights 
% change LMbin to mat:
LMstimbinM = cell2mat(LMdata.LMstimbin);
LMstimbinMean = mean(LMstimbinM,1);
n = height(LMstimbinM);

figure;
hold on;
xWak = 1;
xBef = 2;
xdur = linspace(3,5,length(LMstimbinMean));
xaft = 6;

plot(xWak,LMdata.LMwake,'.','Color','black', 'MarkerSize',20);
plot(xBef,LMdata.LMpre,'.','Color','black','MarkerSize',20);
plot(xdur,LMstimbinM','-','Color',[0.75, 0.75, 0.75]);
plot(xdur,LMstimbinMean,'-','Color','black','LineWidth',2);
plot(xaft,LMdata.LMpost,'.','Color','black','MarkerSize',20);

xlim([0.5,xaft+0.5])
% ylim([0 100])
Groups = {'Wake',' Sleep Before','during stimulation','sleep After'};
xticklabels(Groups)
xticks([xWak,xBef,xdur(round(length(LMstimbinMean)/2)),xaft])
ylabel('Mov/s')

title('Mean movement during stimulation, All Nights, norm')
annotation('textbox', [0.5, 0.85, 0.03, 0.1], 'String', ...
    sprintf('n=%i',n), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

xline(xdur(1),'Color','r','LineWidth',2);
xline(xdur(4),'Color','r','LineWidth',2)
    
% save fig
set(gcf,'PaperPosition',[.25 3 8 6])
saveas (gcf, [analysisFolder filesep 'lizMovAllNights_norm.pdf']);

% clearvars -except stimTable SA analysisFolder


%% plots all nights combined:

% STATISTICAL TEST:

[p, tbl, stats] = friedman(headAngAvgTh, 1); % Here, 1 indicates within-subjects design
fprintf('p-value for freidman ANOVA test: %.5f\n',p)
% p-valure is very low, post hoc:
% Example data for three groups
wakeAng = headAngAvgTh(:,1);
beforeAng = headAngAvgTh(:,2);
duringAng = headAngAvgTh(:,3);
afterAng = headAngAvgTh(:,4);

% Bonferroni-corrected alpha level
alpha = 0.05 / 4;

% Pairwise Wilcoxon signed-rank tests
[p_wake_before, ~, stats_wake_before] = signrank(wakeAng, beforeAng);
[p_before_during, ~, stats_before_during] = signrank(beforeAng, duringAng);
[p_during_after, ~, stats_during_after] = signrank(duringAng, afterAng);
[p_wake_during, ~, stats_wake_during] = signrank(wakeAng, duringAng);

% Display results with Bonferroni correction
fprintf('Wilcoxon signed-rank test results with Bonferroni correction:\n');
fprintf('Wake vs Before: p-value = %.4f (Significant if < %.4f)\n', p_wake_before, alpha);
fprintf('Before vs During: p-value = %.4f (Significant if < %.4f)\n', p_before_during, alpha);
fprintf('During vs After: p-value = %.4f (Significant if < %.4f)\n', p_during_after, alpha);
fprintf('Wake vs During: p-value = %.4f (Significant if < %.4f)\n', p_wake_during, alpha);

% plot head angles:
figure;
x1 = 1:4;
plot(x1,headAngAvgTh,'Color',[0.5 0.5 0.5],'Marker','.'); hold on;
plot(x1, mean(headAngAvgTh), 'Color','k','Marker','.','LineWidth',1.5)
xlim([0.7,4.2]); xticks(x1); xticklabels(["Wake","Sleep","Stim","Sleep after"])
ylabel('Avg Head Angles (Deg)')
yline(0,'--','Headstage penpendicular to floor')
title('Head Angle avg - all')
n= height(headAngAvgTh);
N = numel(unique(stimTable.Animal(inds)));


annotation('textbox', [0.8, 0.85, 0.03, 0.1], 'String', ...
    sprintf('n=%i,N=%i',n,N), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

annotation('textbox', [0.1, 0.8, 0.4, 0.1], 'String', ...
    sprintf('p-value for Friedman ANOVA test: %.5f',p), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

annotation('textbox', [0.1, 0.55, 0.25, 0.1], 'String', ...
    sprintf('wake-before p = %.4f (Significant if < %.4f)', p_wake_before, alpha), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

annotation('textbox', [0.15, 0.2, 0.25, 0.1], 'String', ...
    sprintf('before-during p-value = %.4f ', p_before_during, alpha), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

annotation('textbox', [0.55, 0.65, 0.25, 0.1], 'String', ...
    sprintf('during after p-value = %.4f', p_during_after), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

annotation('textbox', [0.3, 0.1, 0.25, 0.1], 'String', ...
    sprintf('wake during p-value = %.4f', p_wake_during), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');


% savefigure
set(gcf,'PaperPositionMode','auto')
saveas (gcf, [analysisFolder filesep 'HeadliftsAllNights.pdf']);


%% plot stim according to animal and color
% assuming the table is in the workspace
% load([analysisFolder filesep 'stimTable.mat'])

animals = unique(stimTable.Animal);
stimType = ["Blue","Green","Red","LED"];
stimWaveL = ["47","532","635","LED"];
plotColors = {"blue","green","red", [0.5 0.5 0.5]};
numAnimal = length(animals);
numType = length(stimType);

times = stimTable.times{1};
pre=50000;
% STIM DUTRATION NEEDS A THINK!!!!!!!
mStimDur = mean(stimTable.stimDuration(3:end));
post=100000;

f=figure;
set(f, 'Position', [100, 100, 1200, 800]);

ylims = [0 300];
for animal = 1:numAnimal
    for type = 1:numType
        h= subplot(numAnimal,numType, (animal-1)*numType+type);
        
        %plot the data
        curAni = animals{animal};
        curType = stimWaveL(type);
        curTrials = contains(stimTable.Remarks,curType) & contains(stimTable.Animal,curAni);
        n = sum(curTrials);
        curCol = plotColors{type};
        curMean = mean(cell2mat(stimTable.StimAvg(curTrials)),1,'omitnan');
        
        
        
        
        if n>0
            plot(times,curMean,'color',curCol,'LineWidth',4)
            hold on
            for i = 1:height(stimTable)
                if curTrials(i) && stimTable.Mani(i) == 1
                    plot(times,stimTable.StimAvg{i},'Color',curCol)
                end
                if curTrials(i) && stimTable.Mani(i) == 2
                    plot(times,stimTable.StimAvg{i},'Color',curCol,'LineStyle','--')
                end

            end

            xline(pre/1000, 'LineWidth',1,'Color','black')
            xline((pre+mStimDur)/1000,'LineWidth',1,'Color','black')
            hold off
            
            annotation('textbox', [.095 + 0.195*type, 0.98 - 0.142*animal, 0.03, 0.1], 'String', ...
                sprintf('n=%i',n), 'EdgeColor', 'none', 'HorizontalAlignment', ...
                'right', 'VerticalAlignment', 'middle');

        end
        if n==0
            plot(0,0)
        end

        % add titles. labels...
        ylim(ylims);xlabel('Time[S]'), ylabel('D2B power')
        if animal == 1
            title(stimType(type))
        end
%         legend()
      
        pos = get(h, 'Position');
        pos(1) = pos(1) + 0.03; % Shift subplot right to make space for row titles
        set(h, 'Position', pos);
        
        % Add row titles using text function
        if type == 1
            ypos = (ylims(1) + ylims(2)) / 2; % Center vertically

            annotation('textbox', [0.1, 0.98 - 0.142*animal, 0.03, 0.1], 'String', animals{animal}, ...
                'EdgeColor', 'none', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
                'FontWeight', 'bold');
        
        end
    end 
end


% savefigure
set(gcf,'PaperPosition',[.25 3 8 6])
saveas (gcf, [analysisFolder filesep 'animal_type.pdf']);
% clearvars -except stimTable SA analysisFolder
