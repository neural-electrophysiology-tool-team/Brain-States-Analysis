%% sup 3 - SWR:
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

% cross cor: 
% Convert Timestamps to Binary Time Series:
bin_size = 200; % ms 
t_min = min([stims swrT]);
t_max = max([stims swrT]);
time_bins = t_min:bin_size:t_max;

% Convert timestamps to binary spike trains
x = histcounts(stims, [time_bins, time_bins(end)+bin_size]);
y = histcounts(swrT, [time_bins, time_bins(end)+bin_size]);
[c, lags] = xcorr(x, y, 'coeff'); % Normalized cross-correlation. remove 'coeff' for not-normalized

%% plt results:
lag_times = (lags * bin_size)/1000; % Convert lag indices to time (and moving from ms to s
figure;
stem(lag_times, c,'filled','MarkerSize',3);
xlabel(['Lag (s)']);
ylabel('Cross-correlation');
title('Cross-Correlation SWR and stimulations');

%% run get sharp waves on all nights with stim.
for i = 1:height(stimTable)
    recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
    SA.setCurrentRecording(recName);
    SA.getSharpWaves('detectOnlyDuringSWS',false) % if you change to true, make sure you run SlowWave detction first.

end

