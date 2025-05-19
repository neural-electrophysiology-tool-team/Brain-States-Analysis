%% Supplemental Materials
% initiating stimTable and other essentials. 

SA=sleepAnalysis('/media/sil1/Data/Pogona Vitticeps/brainStatesWake.xlsx');
analysisFolder = '/media/sil1/Data/Nitzan/Light Manipulation paper/NitzanAnalysisFiles';
load([analysisFolder filesep 'stimTable.mat'])
% load([analysisFolder filesep 'LMdata.mat'])
animalsColors = [
    255/255, 142/255, 71/255;% HEX:  FF8E47 - orange  - PV126
    28/255, 144/255, 217/255;  % HEX: 1C90D9 - blue - PV149
    148/255, 120/255, 186/255; % HEX: 9478BA - perpule - PV157
    217/255, 62/255, 67/255; % HEX: D93E43 - red - PV159
    255/255, 202/255, 58/255; % HEX: FFCA3A - yellow -  PV161
    97/255, 184/255, 68/255;  % HEX:61B844 - Green -PV162
];
uniqueAnimals = unique(stimTable.Animal);


%% Supplementary Figure 1


%% Supplementary Figure 2

%% Supplementary Figure 3


%% Supplementary Figure 4

%% Supplementary Figure 5

i = 17;
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
SA.setCurrentRecording(recName);
LM = SA.getLizardMovements;%('overwrite',1);
pitchAngles = -LM.angles(2,:);
[angleF, angleF_t] = getHeadLifts(pitchAngles,LM.t_static_ms,100,5);

wakeEnd = 1000*60*60*1;
if wakeEnd>stimTable.sleepStartT(i)
    wakeEnd = stimTable.sleepStartT(i);
end

parts ={[0, wakeEnd],[stimTable.sleepStartT(i),stimTable.stimStartT(i)],...
    [stimTable.stimStartT(i),stimTable.stimEndT(i)], ...
    [stimTable.stimEndT(i),stimTable.sleepEndT(i)]};
numParts = numel(parts);

% plot
f = figure;
plot(angleF_t/(1000*60*60), angleF,'k');
xlabel('Time (hours)'); ylabel('Head Angle')
xline(stimTable.sleepStartT(i)/(1000*60*60),'b')
xline(stimTable.stimStartT(i)/(1000*60*60),'r')
xline(stimTable.stimEndT(i)/(1000*60*60),'Color','r');
xline(stimTable.sleepEndT(i)/(1000*60*60),'b')
yline(0,'--','color',[0.5 0.5 0.5]); ylim([-30 70])
legend({'';'Start Sleep';'Start Stimulations';'End Stimulations';'End Sleep'})

%save figure:
set(f,'PaperPositionMode','auto');
fileName=[analysisFolder filesep 'headAnglesPV159N34'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);

%% Supplementary Figure 6


