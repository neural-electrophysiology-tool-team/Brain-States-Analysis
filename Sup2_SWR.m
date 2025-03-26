%% sup 3 - SWR:

%% plot all raw traces, 1 nights:

i = 17;
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
SA.setCurrentRecording(recName);
    
% get the stimualtions time stamps:
t_ch = stimTable.StimTrighCh(i);
T=SA.getDigitalTriggers;
stims = T.tTrig{t_ch}; %stimulations timimng in ms
    
% get the raw data for the second around each stim:
pre = 200;
post = 1000;
[mV, mV_t] = SA.currentDataObj.getData(SA.recTable.defaulLFPCh(SA.currentPRec),stims-pre,(pre+post));
% try abs for all data:
% mVa = abs(mV);

% plot
fraw = figure;
x = linspace(-pre,post,length(mV_t));
% plot(x,squeeze(mV),'Color',[0.7 0.7 0.7])
mV1 = squeeze(mV);
[d,h] = hist2(repmat(x,height(mV1),1),mV1,'dX1',3,'dX2',5);
hold on
plot(downsample(x,100),downsample(mean(mV1,1),100),'Color','r','LineWidth',2)
ylims =[-750 500];
ylim(ylims)
% title('Raw traces for 1 night')
x_shade = [0 200 200 0];  % X-coordinates of the shaded region
y_shade = [ylims(1) ylims(1) ylims(2) ylims(2)]; % Y-coordinates covering the full y-range
patch(x_shade, y_shade, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
xlabel('Time(ms)');ylabel('Voltage[uV]')

% save figure
set(fraw ,'PaperPosition',[1 2 4.3 1.7]);
fileName=[analysisFolder filesep 'SWRrawPV159N34'];
% print(fileName,'-depsc','-vector');
print(fileName, '-dpdf', '-r300');

%% plot cros corr 1 night:
% detect SWR using SWR detection method.

% for 1 night first:
i = 22;
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
SA.setCurrentRecording(recName);
SA.getSharpWaves('detectOnlyDuringSWS',false)
% gettig the SWR timestamps
SW = SA.getSharpWaves('detectOnlyDuringSWS',false);
swrT = SW.tSW'; %SWR timings in ms
swrN = length(swrT); 
% get the stimualtions time stamps:
t_ch = stimTable.StimTrighCh(i);
T=SA.getDigitalTriggers;
stims = T.tTrig{t_ch}; %stimulations timimng in ms
first_stims = stims(1:8:end-2);
% cross cor: Short, one night 
% Convert Timestamps to Binary Time Series:
bin_size = 200; % ms 
winShort = 10*1000; % ten seconds

ylims = [-0.04 0.1];
t_min = min([stims swrT]);
t_max = max([stims swrT]);
time_bins = t_min:bin_size:t_max;

% Convert timestamps to binary spike trains
swrBin = histcounts(swrT, [time_bins, time_bins(end)+bin_size]);
stimsBin = histcounts(stims, [time_bins, time_bins(end)+bin_size]);
max_lag = winShort/bin_size;
[c, lags,bounds] = crosscorr(stimsBin,swrBin,"NumLags",max_lag); % Normalized cross-correlation. remove 'coeff' for not-normalized

% plt:
lag_times = (lags * bin_size)/1000; % Convert lag indices to time (and moving from ms to s
f = figure;
subplot (1,2,1)
plot(lag_times, c,'Color','k','Marker','.');
yline(bounds,'Color',[0.7 0.7 0.7],'LineStyle','--')
yline(0,'Color',[0.7 0.7 0.7])
xlabel(['Lag (s)']);
ylim(ylims)
ylabel('Cross-correlation');
title('Cross-Correlation SWR and stimulations');

% cross cor: Long, one night 
% Convert Timestamps to Binary Time Series:
bin_size = 5000; % ms 
winLong = 200*1000; % ten seconds

t_min = min([stims swrT]);
t_max = max([stims swrT]);
time_bins = t_min:bin_size:t_max;

% Convert timestamps to binary spike trains
swrBin = histcounts(swrT, [time_bins, time_bins(end)+bin_size]);
% stimsBin = histcounts(stims, [time_bins, time_bins(end)+bin_size]);
firstStimsBin = histcounts(first_stims, [time_bins, time_bins(end)+bin_size]);
max_lag = winLong/bin_size;
[c, lags,bounds] = crosscorr(firstStimsBin,swrBin,"NumLags",max_lag); % Normalized cross-correlation. remove 'coeff' for not-normalized

% plt:
lag_times = (lags * bin_size)/1000; % Convert lag indices to time (and moving from ms to s
% figure;
subplot(1,2,2)
plot(lag_times, c,'Color','k','Marker','.');
yline(bounds,'Color',[0.7 0.7 0.7],'LineStyle','--')
yline(0,'Color',[0.7 0.7 0.7])
xlabel(['Lag (s)']);
ylim(ylims)
% ylabel('Cross-correlation');
title('Cross-Correlation SWR and stimulations');

% save figure:
set(f,'PaperPosition',[1 2 7 2]);
fileName=[analysisFolder filesep 'SWRcrosPV161N18'];
print(fileName,'-dpdf','-r300');

%% run get sharp waves on all nights with stim.
% or i = 25:height(stimTable)
% i = 25;  
% recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
    % SA.setCurrentRecording(recName);
    % SA.getSharpWaves('detectOnlyDuringSWS',false,'startEnds',[5*60*1000 SA.currentDataObj.recordingDuration_ms]) % if you change to true, make sure you run SlowWave detction first.

    
    %% calc and plot cross all nights
% calc the long and short, keep bounds for each 
% setting general variables, true for all recs and for plot!!!
binSizeS = 200; % ms 
binSizeL = 5000; % 5 sec in ms
winShort = 10*1000; % ten seconds
winLong = 200*1000; %200 se in ms
maxLagShort = winShort/binSizeS;
maxLagLong = winLong/binSizeL;
lagTimeS = ([-maxLagShort:1:maxLagShort] * binSizeS)/1000; % Convert lag indices to time (and moving from ms to s
lagTimeL = ([-maxLagLong:1:maxLagLong] * binSizeL)/1000; 
%% creating variables for saving the data of each experiment:
cShort = zeros(height(stimTable),(2*maxLagShort +1));
boundsShort = zeros(height(stimTable),2);
cLong = zeros(height(stimTable),(2*maxLagLong +1));
boundsLong =zeros(height(stimTable),2);


for i = 1:height(stimTable)
    recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
    SA.setCurrentRecording(recName);

    % gettig the SWR timestamps
    SW = SA.getSharpWaves('detectOnlyDuringSWS',false);
    swrT = SW.tSW'; %SWR timings in ms
    % get the stimualtions time stamps:
    t_ch = stimTable.StimTrighCh(i);
    T=SA.getDigitalTriggers;
    stims = T.tTrig{t_ch}; %stimulations timimng in ms
    firstStims = stims(1:8:end-2);

    % Calculate bin times:
    t_min = min([stims swrT]);
    t_max = max([stims swrT]);
    timeBinsS = t_min:binSizeS:t_max;
    timeBinsL = t_min:binSizeL:t_max;

    % Convert timestamps to binary spike trains - Short
    swrBinS = histcounts(swrT, [timeBinsS, timeBinsS(end)+timeBinsS]);
    stimsBinS = histcounts(stims, [timeBinsS, timeBinsS(end)+timeBinsS]);
    % Calculate short cros cor.
    [curCS,lagsS,curBoundsS] = crosscorr(stimsBinS,swrBinS,"NumLags",maxLagShort); % Normalized cross-correlation. remove 'coeff' for not-normalized
    % save to vector:
    cShort(i,:) = curCS;
    boundsShort(i,:) =  curBoundsS;

    % Convert timestamps to binary spike trains - Long
    swrBinL = histcounts(swrT, [timeBinsL, timeBinsL(end)+timeBinsL]);
    fStimsBinL = histcounts(firstStims, [timeBinsL, timeBinsL(end)+timeBinsL]);
    % Calculate short cros cor.
    [curCL,lagsL,curBoundsL] = crosscorr(fStimsBinL,swrBinL,"NumLags",maxLagLong); % Normalized cross-correlation. remove 'coeff' for not-normalized
    % save to vector:
    cLong(i,:) = curCL;
    boundsLong(i,:) =  curBoundsL;

end


%% SHW save and load
fName = [analysisFolder filesep 'ShWcrossCorrData.mat'];
% save(fName,"cShort","cLong","boundsLong","boundsShort",'-mat')
load(fName)
%% plot: only red nights:

wavelength = '635';
curTrials = contains(stimTable.Remarks,wavelength) ...
            & ~contains(stimTable.Remarks,'Ex');
n = sum(curTrials);
N = length(unique(stimTable.Animal(curTrials)));
curDataShort = cShort(curTrials,:);
curDataLong = cLong(curTrials,:);

fr = figure;
% plot Short according to animaL
subplot(1,2,1)
[~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
curColorMat = animalsColors(animalIndices, :); 
hold on; 
for i = 1:height(curDataShort)
    plot(lagTimeS, curDataShort(i,:),'Color',curColorMat(i,:))%,'Marker','.','MarkerSize',4)
end
plot(lagTimeS, mean(curDataShort,1)','Color','k','LineWidth',1.5)%,'Marker','.','MarkerSize',5)
yline(mean(boundsShort(curTrials,:)),'Color',[0.4 0.4 0.4],'LineStyle','--')
yline(0,'Color',[0.4 0.4 0.4])
ylims = [-0.02 .12];
ylim(ylims)

% plot Long according to animaL
subplot(1,2,2)
[~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
curColorMat = animalsColors(animalIndices, :); 
hold on; 
for i = 1:height(curDataLong)
    plot(lagTimeL, curDataLong(i,:),'Color',curColorMat(i,:))%,'Marker','.','MarkerSize',4)
end
plot(lagTimeL, mean(curDataLong,1)','Color','k','LineWidth',1.5)%,'Marker','.','MarkerSize',5)
yline(mean(boundsLong(curTrials,:)),'Color',[0.4 0.4 0.4],'LineStyle','--')
yline(0,'Color',[0.4 0.4 0.4])
ylim(ylims)

% save figure:
set(fr,'PaperPosition',[1 2 7 2]);
fileName=[analysisFolder filesep 'SWRcrosRedNights'];
print(fileName,'-dpdf','-r300');