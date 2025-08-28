
% plot the nights for specific animal
SA=sleepAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWake.xlsx');
animal = 'PV106';
T = SA.recTable;
curTrials = T(strcmp(T.Animal,animal) & T.Mani == 1,:);

for i = 1:height(curTrials)
    recName = ['Animal=' curTrials.Animal{i} ',recNames=' curTrials.recNames{i}];
    SA.setCurrentRecording(recName);
    SA.plotDelta2BetaSlidingAC;

end



