% check trials
WA = wakeAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWake.xlsx');
analysisFolderWake = '/media/sil1/Data/Nitzan/WakeStatesPaper/plots';

%%
WA.setCurrentRecording(['Animal=PV157,recNames=Wake13']);
arena=WA.getArenaCSVs('trigStop',2000,'overwrite',1);
arena.videoFPS
WA.getDelta2BetaRatio( ...
    'movWin',5000, ...
    'movOLWin',4500);
WA.getFreqBandDetection;
% WA.getLizardMovements;

%% 
WA.setCurrentRecording('Animal=PV157,recNames=Night11');
arena=WA.getArenaCSVs('trigStop',2000);
arena.videoFPS
WA.getDelta2BetaRatio( ...
    'movWin',5000, ...
    'movOLWin',4500);
WA.getFreqBandDetection;
WA.getDelta2BetaAC;
WA.plotDelta2BetaSlidingAC;
% WA.getStatesAC;
% WA.plotStatesAC;