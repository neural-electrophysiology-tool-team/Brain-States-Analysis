%% sup 3 - SWR:

%% plot all raw traces, 1 nights:
recNums = [5 17 19 26 37]; %one for each animal. 149 159 161 157 162 correpondingly. 
figure;
for j=1:5
    i = recNums(j);
    recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
    SA.setCurrentRecording(recName);
    
    % get the stimualtions time stamps:
    t_ch = stimTable.StimTrighCh(i);
    T=SA.getDigitalTriggers;
    stims = T.tTrig{t_ch}; %stimulations timimng in ms
    % first_stims = stims(1:8:end-2);
    
    % get the raw data for the second around each stim:
    pre = 100;
    post = 300;
    [mV, mV_t] = SA.currentDataObj.getData(SA.recTable.defaulLFPCh(SA.currentPRec),stims-pre,(pre+post));
        
    % plot
    subplot(5,1,j)
    x = linspace(-pre,post,length(mV_t));
    plot(x,squeeze(mV),'Color',[0.7 0.7 0.7])
    hold on
    plot(x,mean(squeeze(mV),1),'Color','k','LineWidth',2)
    ylims = [-1100 450];
    ylim(ylims)
    x_shade = [0 200 200 0];  % X-coordinates of the shaded region
    y_shade = [ylims(1) ylims(1) ylims(2) ylims(2)]; % Y-coordinates covering the full y-range
    patch(x_shade, y_shade, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    xlabel('Time(ms)');ylabel('mV')

end
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
% cross cor: 
% Convert Timestamps to Binary Time Series:
bin_size = 200; % ms 
t_min = min([stims swrT]);
t_max = max([stims swrT]);
time_bins = t_min:bin_size:t_max;

% Convert timestamps to binary spike trains
x = histcounts(swrT, [time_bins, time_bins(end)+bin_size]);
y = histcounts(stims, [time_bins, time_bins(end)+bin_size]);
[c, lags] = xcorr(x, y, 'coeff'); % Normalized cross-correlation. remove 'coeff' for not-normalized

%% plt results:
lag_times = (lags * bin_size)/1000; % Convert lag indices to time (and moving from ms to s
pTmp = lag_times<10 & lag_times>-10;

figure;
stem(lag_times(pTmp), c(pTmp),'filled','MarkerSize',3);
xlabel(['Lag (s)']);
ylabel('Cross-correlation');
title('Cross-Correlation SWR and stimulations');


%% long

% cross cor: 
% Convert Timestamps to Binary Time Series:
bin_size = 5000; % ms 
t_min = min([stims swrT]);
t_max = max([stims swrT]);
time_bins = t_min:bin_size:t_max;

% Convert timestamps to binary spike trains
x = histcounts(swrT, [time_bins, time_bins(end)+bin_size]);
y = histcounts(first_stims, [time_bins, time_bins(end)+bin_size]);
[c, lags] = xcorr(x, y, 'coeff'); % Normalized cross-correlation. remove 'coeff' for not-normalized

% plt results:
lag_times = (lags * bin_size)/1000; % Convert lag indices to time (and moving from ms to s
pTmp = lag_times<200 & lag_times>-200;

figure;
stem(lag_times(pTmp), c(pTmp),'filled','MarkerSize',3);
xlabel(['Lag (s)']);
ylabel('Cross-correlation');
title('Cross-Correlation SWR and stimulations');

%% run get sharp waves on all nights with stim.
% or i = 25:height(stimTable)
i = 25;  
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
    SA.setCurrentRecording(recName);
    SA.getSharpWaves('detectOnlyDuringSWS',false,'startEnds',[5*60*1000 SA.currentDataObj.recordingDuration_ms]) % if you change to true, make sure you run SlowWave detction first.

