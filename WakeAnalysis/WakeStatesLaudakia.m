%% Wake Analysis in Laudakia

WAs = wakeAnalysis('/media/sil2/Data/Lizard/Stellagama/brainStatesSS.xlsx');
analysisFolderWake = '/media/sil1/Data/Nitzan/WakeStatesPaper/plots';
%%
WAs.setCurrentRecording('Animal=SA07,recNames=sleepNight14')
