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
%% plot all raw traces, 1 nights:
% recNums = [5 17 19 26 37]; %one for each animal. 149 159 161 157 162 correpondingly.
% animals = ["PV149","PV159","PV161","PV157","PV162"];
recNums = [27 38 51];
animals = ["157","162","126"];

figure;
for j=1:3
    i = recNums(j);
    recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
    SA.setCurrentRecording(recName);
    
    % get the stimualtions time stamps:
    t_ch = stimTable.StimTrighCh(i);
    T=SA.getDigitalTriggers;
    stims = T.tTrig{t_ch}; %stimulations timimng in ms
    % first_stims = stims(1:8:end-2);
    
    % get the raw data for the second around each stim:
    pre = 200;
    post = 1000;
    [mV, mV_t] = SA.currentDataObj.getData(SA.recTable.defaulLFPCh(SA.currentPRec),stims-pre,(pre+post));
    % try abs for all data:
    mVa = abs(mV);

    % plot
    subplot(3,1,j)
   
    x = linspace(-pre,post,length(mV_t));
    plot(x,squeeze(mVa),'Color',[0.7 0.7 0.7])
    hold on
    plot(x,mean(squeeze(mVa),1),'Color','k','LineWidth',2)
    ylims = [0 300];
    ylim(ylims)
    title(animals(j))
    x_shade = [0 200 200 0];  % X-coordinates of the shaded region
    y_shade = [ylims(1) ylims(1) ylims(2) ylims(2)]; % Y-coordinates covering the full y-range
    patch(x_shade, y_shade, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    xlabel('Time(ms)');ylabel('mV')

end

%% nnight with artifact but no synchronization of sleep cycle:
i = 37;
    recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
    SA.setCurrentRecording(recName);
    plotStimSham(SA);
    SA.plotDelta2BetaRatio



    %% plot Head Angle - RED NIGHTS ONLY"
headAngDiff = diff(HeadAngleAvg,[],2);
% set the zero to 90 Deg, according to accelerometer data ( this is the z
% axis, when it is 90 the accelerometer is penpendicular to the ground)
HeadAngleAvgP = HeadAngleAvg -90;


wavelength = '635';
curTrials = contains(stimTable.Remarks,wavelength) & ~contains(stimTable.Remarks,'Ex') ...
    & (headAngDiff(:,1)>3 | headAngDiff(:,1)<-3); %& contains(stimTable.Animal,curAni);
n = sum(curTrials);
N = length(unique(stimTable.Animal(curTrials)));
groupNames = {'Pre', 'During', 'After'};
curHeadAvg = HeadAngleAvgP(curTrials,1:3);

% STATISTICAL TEST:

[p, tbl, stats] = friedman(curHeadAvg, 1); % Here, 1 indicates within-subjects design
fprintf('p-value for freidman ANOVA test: %.5f\n',p)
% p-valure is very low, post hoc:
% Example data for three groups
wakeAng = curHeadAvg(:,1);
beforeAng = curHeadAvg(:,2);
duringAng = curHeadAvg(:,3);
% afterAng = curHeadAvg(:,4);

% Bonferroni-corrected alpha level
alpha = 0.05 / 3;

% Pairwise Wilcoxon signed-rank tests
[p_wake_before, ~, stats_wake_before] = signrank(wakeAng, beforeAng);
[p_before_during, ~, stats_before_during] = signrank(beforeAng, duringAng);
% [p_during_after, ~, stats_during_after] = signrank(duringAng, afterAng);
[p_wake_during, ~, stats_wake_during] = signrank(wakeAng, duringAng);

% Display results with Bonferroni correction
fprintf('Wilcoxon signed-rank test results with Bonferroni correction:\n');
fprintf('Wake vs Before: p-value = %.4f (Significant if < %.4f)\n', p_wake_before, alpha);
fprintf('Before vs During: p-value = %.4f (Significant if < %.4f)\n', p_before_during, alpha);
% fprintf('During vs After: p-value = %.4f (Significant if < %.4f)\n', p_during_after, alpha);
fprintf('Wake vs During: p-value = %.4f (Significant if < %.4f)\n', p_wake_during, alpha);

% plot head angles:
figure;
x1 = 1:width(curHeadAvg);
[~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
curColorMat = animalsColors(animalIndices, :); 
for i= 1:height(curHeadAvg)
    plot(x1,curHeadAvg(i,:),'Color',curColorMat(i,:),'Marker','.'); hold on;
end
plot(x1, mean(curHeadAvg), 'Color','k','Marker','.','LineWidth',1.5)
xlim([0.7,3.2]); xticks(x1); xticklabels(["Wake","Sleep","Stim"])
ylabel('Avg Head Angles (Deg)')
yline(0,'--','Headstage penpendicular to floor')
title('Head Angle avg - red nights')


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

% annotation('textbox', [0.55, 0.65, 0.25, 0.1], 'String', ...
%     sprintf('during after p-value = %.4f', p_during_after), 'EdgeColor', 'none', 'HorizontalAlignment', ...
%     'right', 'VerticalAlignment', 'middle');

annotation('textbox', [0.3, 0.1, 0.25, 0.1], 'String', ...
    sprintf('wake during p-value = %.4f', p_wake_during), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');


% savefigure
set(gcf,'PaperPositionMode','auto')
saveas (gcf, [analysisFolder filesep 'HeadliftsREDNights.pdf']);
%% plot the bar plot - D/B change - stim sham -  all colors
animals = unique(stimTable.Animal);
stimType = ["Blue","Green","Red","LED"];
stimWaveL = ["47","532","635","LED"];
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
    curTrials = contains(stimTable.Remarks,curType) &...
                ~contains(stimTable.Remarks,'Ex') & ...
                all(~isnan(stimTable.dbDiffStimM),2) & ...
                all(~isnan(stimTable.dbDiffShamM),2); %& contains(stimTable.Animal,curAni);
    n = sum(curTrials);
    N = length(unique(stimTable.Animal(curTrials)));
    % curCol = plotColors{type};
    curMeanNdbStim = mean(stimTable.dbDiffStimM(curTrials),1);
    curMeanNdbSham = mean(stimTable.dbDiffShamM(curTrials),1);
    curData = [stimTable.dbDiffShamM(curTrials), stimTable.dbDiffStimM(curTrials)];
    %statistics:
    % 1. Wilcoxon Signed-Rank Test
    [pWilcoxon, ~, statsWilcoxon] = signrank(curData(:,1), curData(:,2));
    fprintf('Wilcoxon Signed-Rank Test p-value: %.4f\n', pWilcoxon);
    disp(statsWilcoxon)

    if n>0
        x=[1,2];
        [~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
        curColorMat = animalsColors(animalIndices, :);
        hold on
        for j = 1:height(curData)
            plot(x,curData(j,:),'Color',curColorMat(j,:),'Marker','.')
        end
        hold on
        plot(x,[curMeanNdbSham,curMeanNdbStim],'Color','k','LineWidth',3,'Marker','.')
        

        %plot Ex
        curExI = contains(stimTable.Remarks,curType) &...
                contains(stimTable.Remarks,'Ex') & ...
                all(~isnan(stimTable.dbDiffStimM),2) & ...
                all(~isnan(stimTable.dbDiffShamM),2);
        curEx = [stimTable.dbDiffShamM(curExI) stimTable.dbDiffStimM(curExI)];
        if ~isempty(curEx )
            plot(x,curEx,'Color',[0.3 0.3 0.3],'LineStyle','--','Marker','.')
        end
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
    
    ylim([-40 160])

    if n==0
        plot(0,0)
    end

    % add titles. labels...
    ylabel('D2B power')
    title(stimType(type))
 
end

% savefigure
set(f,'PaperPositionMode','auto');
fileName=[analysisFolder filesep 'meanNormBDStimSham'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);
%% PLOT P2V
    % subplot(2,4,type+4)
    % curP2V = stimTable.ACcomP2V(curTrials,:);
    % curP2Vmean = mean(curP2V,1);
    % hold on
    % ylim([0 1.3])
    % if type ==1
    %     ylabel('P2V')
    % end
    % xticks(1:3); xlim([0.9 3.1])
    % xticklabels({'Pre','During','Post'})
    % for j = 1:height(curP2V)
    %     plot(x,curP2V(j,:),'Color',curColorMat(j,:),'Marker','.')
    % end
    % plot(x,curP2Vmean,'Color','k','LineWidth',2,'Marker','.')
    % notTrials = contains(stimTable.Remarks,curType) & ...
    %             ~contains(stimTable.Remarks,'Ex') &...
    %             all(~isnan(stimTable.ACcomPer),2) &...
    %             any(stimTable.ACcomP2V < 0.15,2);
    % yline(0.15,'r')
    % if any(notTrials~=0)
    %     plot(x,stimTable.ACcomP2V(notTrials,:),'Color',[0.3 0.3 0.3])
    % end
    % hold off

    %% plot P2V Times:
f=figure;
plotColors = {[0 0.586 0.9766],[0.05 0.81 0.379],[1 0.27 0.27], [0.5 0.5 0.5]};

set(f, 'Position', [100, 200, 800, 300]);
hold on
for type = 1:numType
    %plot the data
    subplot(1,4,type)
    curType = stimWaveL(type);
    curTrials = contains(stimTable.Remarks,curType);
    n = sum(curTrials);
    curCol = plotColors{type};
    curMean = mean(stimTable.ACcomP2V(curTrials,:)/1000,1,'omitnan');
    
    if n>0

        plot(x,curMean,'color',curCol,'LineWidth',4)
        hold on
        for animal = 1:numAnimal
            curAni = animals{animal};
            curcurT = contains(stimTable.Remarks,curType) & contains(stimTable.Animal,curAni);
            if any(curcurT)
            plot(x, stimTable.ACcomP2V(curcurT,:)/1000,'Marker',markers{animal},'Color',curCol ...
                ,'MarkerFaceColor',curCol)
            end
        end
        hold off

        annotation('textbox', [.07 + 0.202*type 0.85, 0.03, 0.1], 'String', ...
            sprintf('n=%i',n), 'EdgeColor', 'none', 'HorizontalAlignment', ...
            'right', 'VerticalAlignment', 'middle');
    end
    ylabel('P2V')
    ylim([0 0.0015])
    xticklabels({'Pre','During','Post'})
    xticks(1:3); xlim([0.5 3.5])
end

sgtitle ('P2V According to stim wavelangth and animal')

% savefigure
set(gcf,'PaperPosition',[.25 3 6 2])
saveas (gcf, [analysisFolder filesep 'ACP2VstimColor.pdf']);
% clearvars -except stimTable SA analysisFolder

 %% plot D/B decrease - all nights

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
                all(~isnan(stimTable.dbSWMeans), 2);
    n = sum(curTrials);
    N = numel(unique(stimTable.Animal(curTrials)));
    % curCol = plotColors{type};
    curData = stimTable.dbSWMeans(curTrials,:);
    curMean = mean(stimTable.dbSWMeans(curTrials,:),1,'omitnan');
    statsdbSW = struct();
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
    statsdbSW.(curName).alpha = alpha;
    statsdbSW.(curName).pAnova = p;
    statsdbSW.(curName).p_pre_post = p_pre_post;
    statsdbSW.(curName).p_pre_during = p_pre_during;
    statsdbSW.(curName).p_during_post = p_during_post;


    if n>0
        
        [~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
        curColorMat = animalsColors(animalIndices, :); 
        hold on
        for j = 1:height(curData)
            plot(x,curData(j,:),'Color',curColorMat(j,:),'Marker', '.')
        end
            plot(x,curMean,'color','k','LineWidth',3,'Marker', '.')
        
            % plot ExLight:
            ExTrials = contains(stimTable.Remarks,curType) & ...
                contains(stimTable.Remarks,'Ex') &...
                all(~isnan(stimTable.dbSWMeans),2);
            if sum(ExTrials)>0
                plot(x,stimTable.dbSWMeans(ExTrials,:),'Color',[0.2 0.2 0.2],'LineStyle','--','Marker','.')
            end
        hold off

        annotation('textbox', [.05 + 0.202*type 0.85, 0.03, 0.1], 'String', ...
            sprintf('n=%i,N=%i',n,N), 'EdgeColor', 'none', 'HorizontalAlignment', ...
            'right', 'VerticalAlignment', 'middle');
    end
    ylabel('Time[s]')
    ylim([0 600])
    xticklabels({'Pre','During','Post'})
    xticks(1:3); xlim([0.5 3.5])
end

sgtitle ('D/B decrease during SWS bouts according to wavelangth ')


% savefigure
set(gcf,'PaperPositionMode','auto');
fileName=[analysisFolder filesep 'DBdecreaseAllnigthscolors'];
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
xticks(1:3)
xticklabels(groupNames)
ylabel('1/Beta means during full cycle')
ylims = [0 0.14];
ylim(ylims)
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
fileName=[analysisFolder filesep 'BetachangeFullCycle'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);

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
ylims = [0 0.25];
ylim(ylims)
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


%% plot head movements - Red nights!
% swarm plot - smae color. 
wavelength = '635';
curTrials = contains(stimTable.Remarks,wavelength) & ...
    ~contains(stimTable.Remarks,'Ex') & ...
     ~any(isempty(stimTable.LM_DBt), 2) &...
     ~cellfun(@isempty, stimTable.LM_DBt);

n = sum(curTrials);
N = length(unique(stimTable.Animal(curTrials)));
colors = [0.5 0.5 0.5;0.2, 0.6, 0.8; 0.9, 0.4, 0.3; 0.5, 0.8, 0.5];

LMpre = LMData.LMpre(curTrials, :);
LMwake = LMData.LMwake(curTrials, :);
LMpost = LMData.LMpost(curTrials, :);
LMstimbinM = cell2mat(LMData.LMstimbin); % takes out the nan val
LMstimbintrialM = mean(LMstimbinM(curTrials),2); % mean for each night


LMplotData = [LMwake, LMpre, LMstimbintrialM, LMpost];

% check the statistics:
% Assuming data in columns where each row is a subject and each column is a timepoint

[p, tbl, stats] = friedman(LMplotData, 1); % Here, 1 indicates within-subjects design
fprintf('p-value for freidman ANOVA test: %.5f\n',p)
% p-valure is very low, post hoc:
% Bonferroni-corrected alpha level
alpha = 0.05 / 6;
% Pairwise Wilcoxon signed-rank tests
[p_wake_pre, ~, stats_wake_pre] = signrank(LMwake, LMpre);
[p_wake_during, ~, stats_wake_during] = signrank(LMwake, LMstimbintrialM);
[p_wake_after, ~, stats_wake_after] = signrank(LMwake, LMpost);
[p_pre_during, ~, stats_pre_during] = signrank(LMstimbintrialM,LMpre);
[p_during_after, ~, stats_during_after] = signrank(LMstimbintrialM, LMpost);
[p_pre_after, ~, stats_pre_after] = signrank(LMpre, LMpost);

% Display results with Bonferroni correction
fprintf('Wilcoxon signed-rank test results with Bonferroni correction:\n');
fprintf('Wake vs Pre: p-value = %.5f (Significant if < %.4f)\n', p_wake_pre, alpha);
fprintf('Wake vs During: p-value = %.5f (Significant if < %.4f)\n', p_wake_during, alpha);
fprintf('Wake vs Post: p-value = %.5f (Significant if < %.4f)\n', p_wake_after, alpha);
fprintf('During vs Pre: p-value = %.5f (Significant if < %.4f)\n', p_pre_during, alpha);
fprintf('During vs Post: p-value = %.5f (Significant if < %.4f)\n', p_during_after, alpha);
fprintf('Pre vs Post: p-value = %.5f (Significant if < %.4f)\n', p_pre_after, alpha);

fileName = [analysisFolder filesep 'postHocPvalLMRedNights.mat'];
save(fileName, 'alpha','p_wake_pre','p_wake_during','p_wake_after','p_pre_during', ...
    'p_during_after','p_pre_after')


fLMr = figure;
set(fLMr, 'PaperPositionMode','auto')
Groups = ["Wake","Pre","During Stim","Post"];
% colors = [0.5 0.5 0.5;0.2, 0.6, 0.8; 0.9, 0.4, 0.3; 0.5, 0.8, 0.5];
% violin(LMplotData, 'xlabel',Groups, 'facecolor',colors,'edgecolor',[0.6 0.6 0.6],'medc',[]);

x = [ones(length(LMplotData),1) 2*ones(length(LMplotData),1) 3*ones(length(LMplotData),1) 4*ones(length(LMplotData),1)];
swarmchart(x,LMplotData,15,colors,'filled','XJitterWidth',0.1);
xticks(1:4), xticklabels(Groups); xlim([0.7 4.3]);
ylim([0 5.5]);

annotation('textbox', [0.4, 0.85, 0.03, 0.1], 'String', ...
    sprintf('n=%i,N=%i',n,N), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

annotation('textbox', [0.5, 0.65, 0.3, 0.1], 'String', ...
    sprintf(['Wake vs Pre: p-value = %.5f (Significant if < %.4f)\nWake vs During: ' ...
    'p-value = %.5f \nWake vs Post: p-value = %.5f \nDuring vs Pre: p-value = ' ...
    '%.5f \nDuring vs Post: p-value = %.5f \nPre vs Post: p-value = %.5f \n'],...
    p_wake_pre, alpha,p_wake_during,p_wake_after,p_pre_during, p_during_after,p_pre_after), ...
    'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

% savefigure
set(fLMr,'PaperPositionMode','auto');
fileName=[analysisFolder filesep 'LMRednights2'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);

% clearvars -except stimTable SA analysisFolder LMdata

%% plot P2V diff between stim-pre
numType =4;
P2Vdiffs = [];
groupNums = [];
colorMat = [];
stimWaveL = ["47","532","635","LED"];
stimType = ["Blue","Green","Red","LED"];

for type = 1:numType
    curType = stimWaveL(type);
    curTrials = contains(stimTable.Remarks,curType) & ...
                ~contains(stimTable.Remarks,'Ex') &...
                all(~isnan(stimTable.ACcomPer),2);
    curP2Vdiff = stimTable.ACcomP2V(curTrials,2)-stimTable.ACcomP2V(curTrials,3);
    P2Vdiffs = [P2Vdiffs; curP2Vdiff];
    groupNums = [groupNums; repmat(type, length(curP2Vdiff), 1)];
    [~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
    curColorMat = animalsColors(animalIndices, :);
    colorMat = [colorMat; curColorMat];
end


[p, tbl, statsP2V] = kruskalwallis(P2Vdiffs,groupNums);

if p < 0.05
    cP2V = multcompare(statsP2V, 'CType', 'dunn-sidak');
end

f = figure;
swarmchart(groupNums,P2Vdiffs,10,colorMat,'filled','XJitterWidth',0.5);
ylabel('P2V diff (Stim-Pre)')
xticks(1:4);xticklabels(stimType)
yline(0);
% savefigure
set(f,'PaperPosition',[1 1 1.8 1.2]);
fileName=[analysisFolder filesep 'P2vDiffallcolors'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);