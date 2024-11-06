%% light manipulation analysis 08/05/2024
% This is the new version. 

SA=sleepAnalysis('/media/sil1/Data/Pogona Vitticeps/brainStatesWake.xlsx');
analysisFolder = '/media/sil3/Data/Pogona_Vitticeps/NitzanAnalysisFiles';
load([analysisFolder filesep 'stimTable.mat'])
%% analysis folder
% analysisFolder = '/media/sil3/Data/Pogona_Vitticeps/NitzanAnalysisFiles';
% SA.batchProcessData('getDelta2BetaRatio',{})

%% get the stim sham for single rec
SA.setCurrentRecording('Animal=PV162,recNames=Night23');
getStimSham(SA,11);
SA.getDelta2BetaRatio('tStart', 0 ,'win',1000*60*60*11);
SA.getDelta2BetaAC;
SA.plotDelta2BetaRatio
SA.plotDelta2BetaAC;
SA.plotDelta2BetaSlidingAC;
plotStimSham(SA)


%% get all the stim sham avg from all recs and put in a new table
% this part goes over all the records in SA.
%  for every recoerd that is tagged (1/2/3..) 
SA=sleepAnalysis('/media/sil1/Data/Pogona Vitticeps/brainStatesWake.xlsx');
% SA.setCurrentRecording('Animal=PV162,recNames=Night27');
maniRecs = SA.recTable.Mani>0; % taking all the rows with manipulation
stimTable = SA.recTable(maniRecs,{'Animal','recNames','Remarks','Mani','LizMov','StimTrighCh'});  % creating new table
stimTable.StimAvg = cell(height(stimTable),1);
stimTable.StimAvgSham = cell(height(stimTable),1);
stimTable.times = cell(height(stimTable),1);
stimTable.stimDuration = zeros(height(stimTable),1);
stimTable.ACpre = cell(height(stimTable),1);
stimTable.ACstim = cell(height(stimTable),1);
stimTable.ACpost = cell(height(stimTable),1);

%% geting the AND AC stim sham data from all records 

for i = 1:height(stimTable)

    % set te current rec:
    recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
    SA.setCurrentRecording(recName);
    % run all required analysis:
    SA.getDelta2BetaRatio;
    SA.getDelta2BetaAC;
    SA.getDigitalTriggers
    % get and save the stim avg
    s = getStimSham(SA,stimTable.StimTrighCh(i),1);
    disp('got Stim for this rec')
    stimTable.StimAvg(i) = {mean(s.StimDB,1)};
    stimTable.StimAvgSham(i) = {mean(s.StimDBSham,1)};
    stimTable.times(i) = {s.ts};
    stimTable.stimDuration(i) = s.stimDur;
    disp('stimsham in table')

    % stimulations timings:
    t_ch = stimTable.StimTrighCh(i);
    T=SA.getDigitalTriggers;
    stimStartT = T.tTrig{t_ch}(1);
    stimEndT = T.tTrig{t_ch}(end);
    firstTrig=T.tTrig{t_ch}(1:8:end-2);
    endStim=T.tTrig{t_ch}(8:8:end)+400;
%     stimDuration=(endStim(1)-firstTrig(1));

     % calculate the AC and the P2V for each part of the stimulation
    p = 30*60*1000; %some time diff for the cycle to change.
    ACwin = 2*60*60*1000; % 2 hrs in ms

    SA.getDelta2BetaAC('tStart',stimStartT-ACwin, 'win',ACwin,'overwrite',1);
    stimTable.ACpre(i) = {SA.getDelta2BetaAC('tStart',stimStartT-ACwin, 'win',ACwin)};
    SA.getDelta2BetaAC('tStart',stimStartT+p, 'win',ACwin,'overwrite',1);
    stimTable.ACstim(i) = {SA.getDelta2BetaAC('tStart',stimStartT+p, 'win',ACwin)};
    SA.getDelta2BetaAC('tStart',stimEndT+p,'win',ACwin,'overwrite',1);
    stimTable.ACpost(i) = {SA.getDelta2BetaAC('tStart',stimEndT+p,'win',ACwin)};
    disp('AC in table')
    
      
end
clear recName
clear s
%% GETTING THE LIZ MOV FOR ALL RECS + add to table
SA=sleepAnalysis('/media/sil1/Data/Pogona Vitticeps/brainStatesWake.xlsx');
stimTable.LM_DBt = cell(height(stimTable),1);
for i = 1:height(stimTable)
    % set te current rec:
    recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
    SA.setCurrentRecording(recName);
    % run all required analysis:
    if stimTable.LizMov(i) ==1
        SA.getLizardMovements
        LM = SA.getLizardMovements;

        DB = SA.getDelta2BetaRatio;
        

        % calculate the number of movements to each DB bin
        LM_DBt = zeros(size(DB.t_ms));
        % Loop through each bin in DB and count the events in LM that fall within each bin
        for j = 1:length(DB.t_ms)-1
            % Count events from LM that fall within the current bin (DB(i) to DB(i+1))
            LM_DBt(j) = sum(LM.t_mov_ms >= DB.t_ms(j)& LM.t_mov_ms < DB.t_ms(j+1));
        end
        % Count any events at the last bin edge
        LM_DBt(end) = sum(LM.t_mov_ms >= DB.t_ms(end));

        %put in stimTable:
        stimTable.LM_DBt(i) = {LM_DBt};
        disp('LM in stimTabl')
    end
end


%% save stimTable
save([analysisFolder filesep 'stimTable.mat'], "stimTable",'-mat');

%% calculate the AC analysis again on all night for further analysis. 

for i = 1:height(stimTable)
    % set te current rec:
    recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
    SA.setCurrentRecording(recName);
    SA.getDelta2BetaAC('tStart', 0, 'overwrite',1)
end
%% look at single traces from one night:

i = 22;
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
% recName2 = 'Animal=PV157,recNames=Night37';
SA.setCurrentRecording(recName);
DB = SA.getDelta2BetaRatio;
% SA.getDelta2BetaRatio
% SA.plotDelta2BetaRatio

% stimulations timings:
t_ch = stimTable.StimTrighCh(i);
T=SA.getDigitalTriggers;
stims = T.tTrig{t_ch};
stimStartT = stims(1);
stimEndT = stims(end);
firstTrig=stims(1:8:end-2);
endStim=stims(8:8:end)+400;
trial = reshape(stims,[8,length(stims)/8])';

pre=20000;
post=100000;
ch = 17;

%% plot:

 % for j=22:numel(firstTrig)
    j = 36;
    pTmp=find(DB.t_ms>(firstTrig(j)-pre) & DB.t_ms<(firstTrig(j)+post));
    %StimDB(i,:)=1./DB.bufferedDelta2BetaRatio(pTmp);
    dbTmp=DB.bufferedDelta2BetaRatio(pTmp);
    [lfp,lfp_t] = SA.currentDataObj.getData(ch,firstTrig(j)-pre,pre+post);
    
    % add a raster plot according to stimulation
    % Load spike times and cluster IDs (adjust file names as needed)
    curPhyFoler = [SA.currentDataObj.recordingDir filesep 'spikeSorting' filesep 'Kilosort4'];
    spikes=load([curPhyFoler filesep 'spike_data.mat']);

    % spikeTimes = readNPY([curPhyFoler filesep 'spike_times.npy'] ); % spike times in samples.
    spikeTimes_ms = spikes.spike_times/ (SA.currentDataObj.samplingFrequency(ch)/1000); %spike times in ms.
    spikeClusters = spikes.spike_clusters; % cluster ID per spike
    cluster_info = readtable([curPhyFoler filesep 'cluster_info.tsv'], 'FileType', 'text', 'Delimiter', '\t');

    % Parameters for the plot
    goodClusters  = cluster_info.cluster_id(find(contains(cluster_info.group,'good')|contains(cluster_info.group,'mua')));
    % uniqueClusters = unique(spikeClusters); % all neuron/unit ID
    nClusters = length(goodClusters);
   
%% plot
    figure;
    h1 = subplot(2,1,1);
    % set(f, 'Position', [100, 100, 1200, 400]);
    % yyaxis left
    plot(lfp_t/1000,squeeze(lfp),'k'); hold on;
    curstims = trial(j,:) -trial(j,1) +pre ;
    xline(curstims/1000,'r')
    % yyaxis right
    plot(dbTmp,'Color','b','LineWidth',2)
    sgtitle(sprintf('Trial num: %i',j));
    hold off;

    % Plot RASTER!! spikes for each cluster, for each trig time!
    h2 = subplot(2,1,2);
    hold on;
    for l = 1:nClusters
        clusterID = goodClusters(l);
        
        clusterSpikeTimes = spikeTimes_ms(spikeClusters == clusterID); % times for this unit
        curClusterSpikeT_ms=clusterSpikeTimes(clusterSpikeTimes>(firstTrig(j)-pre) & clusterSpikeTimes<(firstTrig(j)+post));
        TClusterSpikeT_ms = (curClusterSpikeT_ms-(firstTrig(j)-pre));
        % Plot each spike as a tick at y = cluster number
        for k = 1:length(TClusterSpikeT_ms)
            line([TClusterSpikeT_ms(k), TClusterSpikeT_ms(k)], [l - 0.4, l + 0.4], 'Color', 'k'); % tick mark
        end

    end
    curstims = trial(j,:) -trial(j,1) +pre ;
    xline(curstims,'r')
    box on
    % Convert x-axis labels from ms to s by setting the x-axis ticks and labels
    xticks = get(gca, 'XTick');          % Get current x-axis tick values in ms
    set(gca, 'XTick', xticks);           % Set the same ticks
    set(gca, 'XTickLabel', xticks / 1000); % Display tick labels in seconds
    xlabel('Time (s)');
    % xlabel('Time (ms)');
    ylabel('Neuron/Unit');
    ylim([1,nClusters+1])
    title('Raster Plot');
    hold off;

    % saveas(gcf, [SA.currentPlotFolder filesep sprintf('DBRaterLFPT%i.pdf',j)])
    % waitforbuttonpress;

% end


% saveas (gcf, [analysisFolder filesep 'singleTrialRasterPV161N13t36.pdf']);


%% spike rates for each trigger: 

%get the mean spike rate:
spikeRateT = 0:OL:win-meanWin; %time in ms
spikeRateM = zeros(numel(firstTrig),length(spikeRateT));
win = pre+post;
meanWin = 1000;
OL =100;
for j = 1:numel(firstTrig)
    [curSpikeRate] = getSpikeRate(spikes,goodClusters,firstTrig(j)-pre,win,meanWin,OL);
    spikeRateM(j,:)= mean(curSpikeRate,1);
end
% save([curPhyFoler filesep 'spikeRateTrigers.mat'],"spikeRateM","spikeRateT",'-mat')
load([curPhyFoler filesep 'spikeRateTrigers.mat'],"spikeRateM")

%% plot
figure;
imagesc(spikeRateM);
colormap(flipud(gray))  % Set colormap to grayscale
colorbar;        % Optional: show colorbar
% Set x-tick values and labels for time
% timeLabels = (0:10:(length(spikeRateT)-1)/10);

%Set x-tick values and labels for every 100 columns (every 10 seconds)
% xTickPositions = 1:100:length(spikeRateT);  % Positions at every 10 seconds
% xticks(xTickPositions);  % Set x-ticks
% xticklabels(arrayfun(@(t) sprintf('%.1f', t), timeLabels, 'UniformOutput', false));  % Set labels

xlabel('Time[100 ms]')
ylabel('Trial #')
% Hold the plot to add lines
hold on;

% Define X positions for the red lines (customize as needed)
stimDiff = mean(mean(diff(trial,[],2)));
xPositions = (pre + (0:7)*(stimDiff))/100;  % Example positions for lines

% Add red vertical lines
for i = 1:length(xPositions)
    line([xPositions(i), xPositions(i)], [1, size(spikeRateM, 1)], 'Color', 'red', 'LineWidth', 0.5);
end

hold off;
title('spikeRate (spike per sec) average per stim trial, PV161N18')
% 
saveas (gcf, [analysisFolder filesep 'spikeRateAllTrialsPV161N18.pdf']);

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


%% PLOTINGS

%% plot D2B general decrease for each part of the stimulation
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
    [peakPre, peakLocPre] = findpeaks(DBpre, 'MinPeakHeight', 50);
    [peakstim, peakLocStim] = findpeaks(DBstim, 'MinPeakHeight', 50);
    [peakPost, peakLocPost] = findpeaks(DBpost, 'MinPeakHeight', 50);

    % % peak analysis, no min peaks: 
    % [peakPre, peakLocPre] = findpeaks(DBpre);
    % [peakstim, peakLocStim] = findpeaks(DBstim);
    % [peakPost, peakLocPost] = findpeaks(DBpost);



%     DBpeaksPre(i) = mean(peakPre);
%     DBpeaksStim(i) = mean(peakstim);
%     DBpeaksPost(i) = mean(peakPost);
% end


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

%% plot sliding AC sith stimulations:
%set the recording:
i = 22;
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
SA.setCurrentRecording(recName);
AC = SA.getDelta2BetaAC;
SA.plotDelta2BetaSlidingAC

%% AC - plot specific before during and after AC - for a specific rec/
i = 22;
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
SA.setCurrentRecording(recName);
ACpre = stimTable.ACpre{i};
ACstim = stimTable.ACstim{i};
ACpost = stimTable.ACpost{i};

ACstructs = {ACpre,ACstim,ACpost};
labels = {'preStim', 'Stim', 'postStim'};

for j = 1:3
    struct2vars(ACstructs{j});
    fAC = figure;
    h = axes;
    %plot:
    lineHandles = stem(autocorrTimes/1000,real(xcf),'filled','r-o');
    ylim([-0.4 1])
    set(lineHandles(1),'MarkerSize',4);
    grid('on');
    xlabel('Period [s]');
    ylabel('Auto corr.');
    hold('on');

    plot(period/1000,real(xcf(pPeriod)),'o','MarkerSize',5,'color','k');
    text(period/1000,0.05+real(xcf(pPeriod)),num2str(period/1000));

    a = axis;
    plot([a(1) a(1); a(2) a(2)],[xcf_bounds([1 1]) xcf_bounds([2 2])],'-b');
    plot([a(1) a(2)],[0 0],'-k');
    hold('off');
    
    % save fig:

    set(fAC,'PaperPositionMode','auto');
    fileName=[SA.currentPlotFolder filesep 'dbAC_ch' num2str(parDbAutocorr.ch) '_t' num2str(parDbAutocorr.tStart) '_w' num2str(parDbAutocorr.win) labels{j} '.pdf'];
    saveas(fAC,fileName)
    % print(fileName,'-pdf',['-r' num2str(SA.figResJPG)]);
    % if printLocalCopy
    %     fileName=[cd filesep obj.recTable.Animal{obj.currentPRec} '_Rec' num2str(obj.currentPRec) '_dbAC_ch' num2str(parDbAutocorr.ch) '_t' num2str(parDbAutocorr.tStart) '_w' num2str(parDbAutocorr.win)];
    %     print(fileName,'-pdf',['-r' num2str(obj.figResJPG)]);
    % end

end





%% AC - get the Data 

ACprePer = zeros(height(stimTable),1);
ACpreP2V = zeros(height(stimTable),1);
ACstimPer = zeros(height(stimTable),1);
ACstimP2V = zeros(height(stimTable),1);
ACpostPer = zeros(height(stimTable),1);
ACpostP2V = zeros(height(stimTable),1);

for i = 1:height(stimTable)
   curACpre = stimTable.ACpre{i};
   curACstim = stimTable.ACstim{i};
   curACpost = stimTable.ACpost{i};
   ACprePer(i) = curACpre.period;
   ACpreP2V(i) = curACpre.peak2VallyDiff;
   ACstimPer(i) = curACstim.period;
   ACstimP2V(i) = curACstim.peak2VallyDiff;
   ACpostPer(i) = curACpost.period;
   ACpostP2V(i) = curACpost.peak2VallyDiff;
   
end

%% plot AC - Basic
ACcomPer = [ACprePer ACstimPer ACpostPer];
ACcomP2V = [ACpreP2V ACstimP2V ACpostP2V];
stimTable.ACcomPer = ACcomPer;
stimTable.ACcomP2V = ACcomP2V;

figure;
plot(ACcomPer'/1000, '-o')
xlim([0.5 3.5])
figure; plot(ACcomP2V','-o')
%% plot AC - only Red nights:
type = 'Red';
wavelength = '635';
curTrials = contains(stimTable.Remarks,wavelength); %& contains(stimTable.Animal,curAni);
n = sum(curTrials);
N = length(unique(stimTable.Animal(curTrials)));

% statistical tests:
% first, we need to check the differences in genral, in Freidman test,
% which is a a-parametrical ANOVA test. then we can use wilcoxon post-hoc
% to check where is the differenc (with benforoni corection)

% Assuming data in columns where each row is a subject and each column is a timepoint
[p, tbl, stats] = friedman(stimTable.ACcomPer(curTrials,:), 1); % Here, 1 indicates within-subjects design
fprintf('p-value for freidman ANOVA test: %.5f\n',p)
% p-valure is very low, post hoc:
% Example data for three groups
before = stimTable.ACcomPer(curTrials,1);
during = stimTable.ACcomPer(curTrials,2);
after = stimTable.ACcomPer(curTrials,3);

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
figure;
x = 1:3;
curMean = mean(stimTable.ACcomPer(curTrials,:)/1000,1,'omitnan');
plot(x, stimTable.ACcomPer(curTrials,:)/1000,'Color',[0.5 0.5 0.5],'Marker','.')
hold on
plot(x, curMean,'Color','k','LineWidth',2)
xticklabels({'Pre','During','Post'})
xticks(1:3); xlim([0.5 3.5])
grid on
ylabel('Period Time[s]')
title ('Perios Times changes - all red nights')
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
set(gcf,'PaperPositionMode','auto')
saveas (gcf, [analysisFolder filesep 'ACperiodReds.pdf']);
%% plot AC - According to color and animal

animals = unique(stimTable.Animal);
stimType = ["Blue","Green","Red","LED"];
stimWaveL = ["47","532","635","LED"];
plotColors = {"blue","green","red", [0.5 0.5 0.5]};
numAnimal = length(animals);
numType = length(stimType);
markers = {'o', 's', 'd', '^', 'v','x'};
x=1:3;
%plot Period Times:
f=figure;
set(f, 'Position', [100, 100, 800, 400]);
hold on
for type = 1:numType
    %plot the data
    subplot(1,4,type)
    curType = stimWaveL(type);
    curTrials = contains(stimTable.Remarks,curType);
    n = sum(curTrials);
    curCol = plotColors{type};
    curMean = mean(stimTable.ACcomPer(curTrials,:)/1000,1,'omitnan');
    
    if n>0

        plot(x,curMean,'color',curCol,'LineWidth',4)
        hold on
        for animal = 1:numAnimal
            curAni = animals{animal};
            curcurT = contains(stimTable.Remarks,curType) & contains(stimTable.Animal,curAni);
            if any(curcurT)
            plot(x, stimTable.ACcomPer(curcurT,:)/1000,'Marker',markers{animal},'Color',curCol ...
                ,'MarkerFaceColor',curCol)
            end
        end
        hold off

        annotation('textbox', [.05 + 0.202*type 0.85, 0.03, 0.1], 'String', ...
            sprintf('n=%i',n), 'EdgeColor', 'none', 'HorizontalAlignment', ...
            'right', 'VerticalAlignment', 'middle');
    end
    ylabel('Time[s]')
    ylim([40 250])
    xticklabels({'Pre','During','Post'})
    xticks(1:3); xlim([0.5 3.5])
end

sgtitle ('Perios Times According to stim wavelangth and animal')

% savefigure
set(gcf,'PaperPosition',[.25 3 8 6])
saveas (gcf, [analysisFolder filesep 'ACperiodstimColor.pdf']);

%% plot P2V Times:
f=figure;
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
clearvars -except stimTable SA analysisFolder

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
clearvars -except stimTable SA analysisFolder

%% plot the statistics:
% i want a bar plot (with points) of the max values of the D/B during
% stim/Sham
% also, plot the max slop, i think that will show better results.

% assuming the table is in the workspace
% load([analysisFolder filesep 'stimTable.mat'])

% calculate max value and max slop
% max:
stimTable.maxStim = cellfun(@max,stimTable.StimAvg);
stimTable.maxSham = cellfun(@max,stimTable.StimAvgSham);

% stimtime minus pre stim - change in D/B:
% take the last 30 sec of stim and substract the 30 sec before the start of
% stim, ans avrage that. 
stimTable.dbDiffStim = cell(height(stimTable),1);
stimTable.dbDiffSham = cell(height(stimTable),1);
firstStimInd = 50000/1000;
win = 30;
for i=1:height(stimTable)
    CurStimDur = round(stimTable.stimDuration(i)/1000);
    % calc substracted stimulation from baseline
    befStim = stimTable.StimAvg{i}(firstStimInd-win:firstStimInd-1);
    durStim = stimTable.StimAvg{i}(firstStimInd+CurStimDur-win:firstStimInd+CurStimDur-1);
    stimTable.dbDiffStim(i) = {durStim-befStim};
    % calc substracted stimulation from baseline
    befSham = stimTable.StimAvgSham{i}(firstStimInd-win:firstStimInd-1);
    durSham = stimTable.StimAvgSham{i}(firstStimInd+CurStimDur-win:firstStimInd+CurStimDur-1);
    stimTable.dbDiffSham(i) = {durSham-befSham};
end
stimTable.dbDiffStimM = cellfun(@(x) mean(x,'omitnan'),stimTable.dbDiffStim);
stimTable.dbDiffShamM = cellfun(@(x) mean(x,'omitnan'),stimTable.dbDiffSham);

save([analysisFolder filesep 'stimTable.mat'], "stimTable",'-mat');
 
%% plot D/B diff stimSham bar plot - only Red nights
% animals = unique(stimTable.Animal);
% numAnimal = length(animals);
type = 'Red';
wavelength = '635';
curTrials = contains(stimTable.Remarks,wavelength); %& contains(stimTable.Animal,curAni);
n = sum(curTrials);
N = length(unique(stimTable.Animal(curTrials)));

% statistical tests for that figure - diff StimSham only red nights:
% trying different methods: 
groupSham = stimTable.dbDiffShamM(curTrials);
groupStim = stimTable.dbDiffStimM(curTrials);

% 1. Wilcoxon Signed-Rank Test
[pWilcoxon, ~, statsWilcoxon] = signrank(groupSham, groupStim);
fprintf('Wilcoxon Signed-Rank Test p-value: %.4f\n', pWilcoxon);
disp(statsWilcoxon)

% 2. Paired t-test (check normality first)
% Normality test (Shapiro-Wilk) - used
[hNorm1, pNorm1] = adtest(groupSham);
[hNorm2, pNorm2] = adtest(groupStim);

if pNorm1 > 0.05 && pNorm2 > 0.05
    [hTtest, pTtest] = ttest(groupSham, groupStim);
    fprintf('Paired t-test p-value: %.4f\n', pTtest);
else
    fprintf('Data is not normally distributed; Paired t-test may not be appropriate.\n');
end


%plot:
figure;
title('Change in mean D/B norm across trials - All nights')
x=[1,2];
plot(x,[stimTable.dbDiffShamM(curTrials), stimTable.dbDiffStimM(curTrials)] ...
    ,'-o','Color',[0.5, 0.5, 0.5])
hold on
curMeanStim = mean(stimTable.dbDiffStimM(curTrials),1,'omitnan');
curMeanSham = mean(stimTable.dbDiffShamM(curTrials),1,'omitnan');
plot(x,[curMeanSham,curMeanStim],'-o','color','k','LineWidth',2)
hold off

grid on
xlim([0.5 2.5]);
xticks([1, 2]); % Position of the x-ticks
xticklabels({'Sham', 'Stim'}); % Labels for the x-ticks
annotation('textbox', [0.85, 0.85, 0.03, 0.1], 'String', ...
    sprintf('n=%i,N=%i',n,N), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');
annotation('textbox', [0.15, 0.8, 0.3, 0.2], 'String', ...
    sprintf('Wilcoxon test p =%.3f',pWilcoxon), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');
annotation('textbox', [0.15, 0.75, 0.3, 0.2], 'String', ...
    sprintf('T-test p =%.3f',pTtest), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle')

set(gcf,'PaperPositionMode','auto')
saveas (gcf, [analysisFolder filesep 'DBdiffStimShamRedNights.pdf']);


clearvars -except stimTable SA analysisFolder

%% plot the bar plot - D/B change - stim sham -  all colors
animals = unique(stimTable.Animal);
stimType = ["Blue","Green","Red","LED"];
stimWaveL = ["47","532","635","LED"];
plotColors = {"blue","green","red", [0.5 0.5 0.5]};
numAnimal = length(animals);
numType = length(stimType);


f=figure;
set(f, 'Position', [100, 100, 1200, 400]);
sgtitle('Max average D/B')
for type = 1:numType
    h = subplot(1,numType,type);
    %plot the data
    %curAni = animals{animal};
    curType = stimWaveL(type);
    curTrials = contains(stimTable.Remarks,curType); %& contains(stimTable.Animal,curAni);
    n = sum(curTrials);
    N = length(unique(stimTable.Animal(curTrials)));
    curCol = plotColors{type};
    curMeanStim = mean(stimTable.maxStim(curTrials),1,'omitnan');
    curMeanSham = mean(stimTable.maxSham(curTrials),1,'omitnan');
    if n>0
        x=[1,2];
        plot(x,[curMeanSham,curMeanStim],'-o','color',curCol,'LineWidth',3)
        hold on
        plot(x,[stimTable.maxSham(curTrials), stimTable.maxStim(curTrials)] ...
            ,'-o','Color',curCol)
          hold off
    end
    xticks([1, 2]); % Position of the x-ticks
    xticklabels({'Sham', 'Stim'}); % Labels for the x-ticks
    xlim([0.5, 2.5]);
    annotation('textbox', [.095 + 0.195*type, 0.85, 0.03, 0.1], 'String', ...
        sprintf('n=%i,N=%i',n,N), 'EdgeColor', 'none', 'HorizontalAlignment', ...
        'right', 'VerticalAlignment', 'middle');
   
    if n==0
        plot(0,0)
    end

    % add titles. labels...
    ylabel('D2B power')
    title(stimType(type))
    ylim([0 450])
end

% savefigure
set(gcf,'PaperPosition',[.25 3 8 6])
saveas (gcf, [analysisFolder filesep 'maxStimShamAll1.pdf']);
%% mean normelized change in D/B.
f=figure;
set(f, 'Position', [100, 100, 1200, 400]);
sgtitle('mean normelized D/B')
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
        plot(x,[curMeanNdbSham,curMeanNdbStim],'-o','color',curCol,'LineWidth',3)
        hold on
        plot(x,[stimTable.dbDiffShamM(curTrials), stimTable.dbDiffStimM(curTrials)] ...
            ,'-o','Color',curCol)
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

% savefigure
set(gcf,'PaperPosition',[.25 3 8 6])
saveas (gcf, [analysisFolder filesep 'meanNormBDStimSham.pdf']);

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


%% Accelerometer Data analysis:

% for one night:
i = 29;
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
stimDuration=(endStim(1)-firstTrig(1));
%% plot the movement binned to the whole night. 
LM_DBt = stimTable.LM_DBt{i};
figure; 
plot(DB.t_ms/1000, LM_DBt); xlabel('Time[s]'), ylabel('MovCount');
% save fig
set(gcf,'PaperPosition',[.25 3 8 6])
saveas (gcf, [analysisFolder filesep 'wholeNightmove.pdf']);

%% calculate the mean of events in the hour before stim, during the stims,
% and an hour after. 
befind = find(DB.t_ms>stimStartT-1000-(60*60*1000)&DB.t_ms<stimStartT-1000);
befMov = mean(LM_DBt(befind));
wakeMov = mean(LM_DBt(1:3600));
aftind = find(DB.t_ms>stimEndT+(60*60*1000)&DB.t_ms<stimEndT+(2*60*60*1000));
aftMov = mean(LM_DBt(aftind));
stimlength = 150;
stimLM = zeros(numel(firstTrig),stimlength);
post =150*1000;
binSize = 10;
stimLMbin = zeros(numel(firstTrig),stimlength/binSize);
% Add the last bin edge 
for i=1:numel(firstTrig)
    pTmp=find(DB.t_ms>(firstTrig(i)) & DB.t_ms<=(firstTrig(i)+post));
        %StimDB(i,:)=1./DB.bufferedDelta2BetaRatio(pTmp);
        if length(pTmp) ~=stimlength
            pTmp = pTmp(1:stimlength);
        end
        curstimLM = LM_DBt(pTmp);
        stimLM(i,:) = curstimLM;
        
        % Reshape the array to 15x10
        curReshLM = reshape(curstimLM, binSize,stimlength/binSize);

        % Sum along the rows (dim=1) to get a 1x15 array
        curStimLMBin = sum(curReshLM, 1);
        stimLMbin(i,:) = curStimLMBin/binSize;
              
end
stimLMmean = mean(stimLM,1);
stimLMbinmean = mean(stimLMbin,1);

%% plot the data from the stimulation time:
figure; subplot (2,1,1)
plot(stimLM')
hold on; plot(stimLMmean,'black','LineWidth',3); 
ylabel('Event count'), xlabel('Time [s]'); xline(0,'color','r');xline(38,'color','r');
hold off

subplot (2,1,2);plot(stimLMbin')
hold on;
plot(stimLMbinmean,'black','LineWidth',3); 
ylabel('Event count'), xlabel('Time [10s]'); xline(0,'color','r');xline(3.8,'color','r');
hold off
% save fig
set(gcf,'PaperPosition',[.25 3 8 6])
saveas (gcf, [analysisFolder filesep 'movDurStimsinglenight.pdf']);

%% plot the full movement data for a night
figure;
hold on;
xWak = 1;
xBef = 2;
xdur = 3:3+14;
xaft = 18;

plot(xWak,wakeMov,'.','Color','black', 'MarkerSize',20);
plot(xBef,befMov,'.','Color','black','MarkerSize',20);
plot(xdur,stimLMbinmean,'-o','Color','black','MarkerFaceColor','black');
plot(xaft,aftMov,'.','Color','black','MarkerSize',20);
xlim([0,19])
xticklabels({'Wake','sleep Before','during stimulation','sleep After'})
xticks([1,2,10,18])
title('Mean movement during stimulation, PV161,Night18')
% save fig
set(gcf,'PaperPosition',[.25 3 8 6])
saveas (gcf, [analysisFolder filesep 'lizMovWholeNightPV161N18.pdf']);

%% Movement during stimulation - All Nights:
% get data:
LMwake = zeros(height(stimTable),1);
LMpre = zeros(height(stimTable),1);
LMstim = cell(height(stimTable),1);
LMstimbin = cell(height(stimTable),1);
LMpost = zeros(height(stimTable),1);
LMallMean = zeros(height(stimTable),1);

stimlength = 150;
post =150*1000;
binSize = 10;

for i = 1:height(stimTable)
    %set the recording:
    recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
    SA.setCurrentRecording(recName);
    DB = SA.getDelta2BetaRatio;
    LM_DBt = stimTable.LM_DBt{i};
    % AC = SA.getDelta2BetaAC;
    
    %check if LM analysis was already done for this Rec
    if isempty(LM_DBt)
        disp('run Lizard movement on this recording. Moving to next rec.');
        continue; % Skip to the next iteration;
    end

    % stimulations timings:
    t_ch = stimTable.StimTrighCh(i);
    T=SA.getDigitalTriggers;
    stimStartT = T.tTrig{t_ch}(1);
    stimEndT = T.tTrig{t_ch}(end);
    firstTrig=T.tTrig{t_ch}(1:8:end-2);
    endStim=T.tTrig{t_ch}(8:8:end)+400;
    stimDuration=(endStim(1)-firstTrig(1));

    % get and save mov for each part
    LMallMean(i) = mean(LM_DBt(LM_DBt<450));
    
    % deivided byt the general mean of the recording:
    befind = find(DB.t_ms>stimStartT-1000-(60*60*1000)&DB.t_ms<stimStartT-1000);
    LMpre(i) = mean(LM_DBt(befind))/LMallMean(i);
    LMwake(i) = mean(LM_DBt(1:3600))/LMallMean(i);
    aftind = find(DB.t_ms>stimEndT+(60*60*1000)&DB.t_ms<stimEndT+(2*60*60*1000));
    LMpost(i) = mean(LM_DBt(aftind))/LMallMean(i);
    
    stimLM = zeros(numel(firstTrig),stimlength);
    stimLMbin = zeros(numel(firstTrig),stimlength/binSize);
    % Add the last bin edge
    for j=1:numel(firstTrig)
        pTmp=find(DB.t_ms>(firstTrig(j)) & DB.t_ms<=(firstTrig(j)+post));
        %StimDB(i,:)=1./DB.bufferedDelta2BetaRatio(pTmp);
        if length(pTmp) ~=stimlength
            pTmp = pTmp(1:stimlength);
        end
        curstimLM = LM_DBt(pTmp);
        stimLM(i,:) = curstimLM;

        % Reshape the array to 15x10
        curReshLM = reshape(curstimLM, binSize,stimlength/binSize);

        % Sum along the rows (dim=1) to get a 1x15 array
        curStimLMBin = sum(curReshLM, 1);
        stimLMbin(j,:) = curStimLMBin/(binSize*LMallMean(i)); %normelized to sec - to bin size.

    end
    LMstim(i) = {mean(stimLM,1)};
    LMstimbin(i) = {mean(stimLMbin,1)};
    
end



%% PLOT - Movement during stimulation - All Nights 
figure;
hold on;
xWak = 1;
xBef = 2;
xdur = 3:3+14;
xaft = 18;

% change LMbin to mat:
LMstimbinM = cell2mat(LMstimbin);
LMstimbinMean = mean(LMstimbinM,1);
n = height(LMstimbinM);

plot(xWak,LMwake,'.','Color','black', 'MarkerSize',20);
plot(xBef,LMpre,'.','Color','black','MarkerSize',20);
plot(xdur,LMstimbinM','-','Color',[0.75, 0.75, 0.75]);
plot(xdur,LMstimbinMean,'-','Color','black','LineWidth',2);
plot(xaft,LMpost,'.','Color','black','MarkerSize',20);

xlim([0,19])
% ylim([0 100])

xticklabels({'Wake','Sleep Before','during stimulation','sleep After'})
xticks([1,2,10,18])
title('Mean movement during stimulation, All Nights, norm')
annotation('textbox', [0.5, 0.85, 0.03, 0.1], 'String', ...
    sprintf('n=%i',n), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

xline(3,'Color','r','LineWidth',2);
xline(3+4,'Color','r','LineWidth',2)
    

% save fig
set(gcf,'PaperPosition',[.25 3 8 6])
saveas (gcf, [analysisFolder filesep 'lizMovAllNights_norm.pdf']);

clearvars -except stimTable SA analysisFolder



