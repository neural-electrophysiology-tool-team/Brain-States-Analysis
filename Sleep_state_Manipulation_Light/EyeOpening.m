%% eye Opening
eyeInd = contains(stimTable.Remarks,'white')&stimTable.eyeMov==1;
animals =stimTable.Animal(eyeInd);recnames = stimTable.recNames(eyeInd);
recList = cellfun(@(x,y) ['Animal=' x ',recNames=' y], animals, recnames, 'UniformOutput', false);
eyeOpening = cell(length(recList),1);
for i = 2:length(recList)
    SA.setCurrentRecording(recList{i})
    stims = SA.getStimTriggers;
    firstTrig = stims(1:8:end);
    curEyeOpen = [];
    eyeOpening{i} = {curEyeOpen};
end
%% eye opening:
% insert the eyeOpenFrames manually.
% translate to OE times:
SA.setCurrentRecording(recList{5});
camTrigs = SA.getCameraTriggers;
eyeOpenStim = camTrigs(eyeOpenFrames);
filename = [SA.currentAnalysisFolder filesep 'eyeOpenTimes'];
save(filename,"eyeOpenStim","eyeOpenFrames")

% plot stim db cycle, and the eye open events/
SC = SA.getSlowCycles;
stims = SA.getStimTriggers;
firstTrig = stims(1:8:end);
%% plt
p=120*1000;
pTmp = SC.t_ms>stims(1) & SC.t_ms<(stims(end)+p);
ms2h = 1000*60*60;
% plot(SC.t_ms(pTmp), SC.HAng(pTmp));hold on;

figure;
h1 = plot(SC.t_ms(pTmp)/ms2h, SC.DBRatioMedFilt(pTmp));
hold on;

h2 = xline(firstTrig/ms2h, 'Color',[0.5 0.5 0.5]);
h3 = xline(eyeOpenStim/ms2h, 'r');

legend([h1(1), h2(1), h3(1)], "\delta/\beta", "firstStim", "EyeOpen");