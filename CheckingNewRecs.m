%% check new nights recording:
SA = wakeAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWake.xlsx');
SA.setCurrentRecording('Animal=PV153,recNames=Night11');
SA.getDelta2BetaRatio('overwrite',1);
SA.getDelta2BetaAC;
SA.plotDelta2BetaRatio;
SA.plotDelta2BetaSlidingAC;

%% check new hunter rec:
WA = wakeAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWake.xlsx');
WA.setCurrentRecording('Animal=PV106,recNames=Hunter1');
WA.getArenaCSVs(1);
WA.getBeta2GammaRatio;
WA.plotTrailsBG;
A = WA.getArenaCSVs;
