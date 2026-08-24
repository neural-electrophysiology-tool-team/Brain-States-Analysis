% check trials
WA = wakeAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWakeN.xlsx');
analysisFolderWake = '/media/sil1/Data/Nitzan/WakeStatesPaper/plots';

%%
WA.setCurrentRecording(['Animal=PV153,recNames=Wake35']);
arena=WA.getArenaCSVs('trigStop',2000,'overwrite',1);
% arena.videoFPS
WA.getDelta2BetaRatio( 'overwrite',1,...
    'movWin',5000, ...
    'movOLWin',4500);
WA.getFreqBandDetection;
% WA.getLizardMovements;

%% 
WA.setCurrentRecording('Animal=PV153,recNames=Night4');
arena=WA.getArenaCSVs('trigStop',2000,'overwrite',1);
% arena.videoFPS
WA.getDelta2BetaRatio( ...
    'movWin',5000, ...
    'movOLWin',4500);
WA.getFreqBandDetection;
WA.getDelta2BetaAC;
WA.plotDelta2BetaSlidingAC;
% WA.getStatesAC;
% WA.plotStatesAC;