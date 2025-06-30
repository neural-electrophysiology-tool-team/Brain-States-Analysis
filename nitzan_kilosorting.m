%% setting up the spike sorting:

%% Light Manipulation - spike sorting
   % this section is for the same thing but using the stimTable as ref
WA=wakeAnalysis('/media/sil1/Data/Pogona Vitticeps/brainStatesWake.xlsx');
WA.setCurrentRecording('Animal=PV126,recNames=Hunter54');


WA.currentDataObj.generateChannelMapFile('40_16x2_FlexLin'); % 120_32x1_H4_CamNeuro for PV24. all layouts in TSV/electrode layouts 
binaryileName = [WA.currentDataObj.recordingDir filesep 'spikeSorting' filesep 'ch1_32.bin'];
WA.currentDataObj.convert2Binary(binaryileName);
%% loop of all the recoring to create binary files for them
% look for only the ones in stim table + have above 0 in spikes. 

SA=sleepAnalysis('/media/sil1/Data/Pogona Vitticeps/brainStatesWake.xlsx');

for i = 25:height(stimTable)
    if stimTable.spikes(i) >0 
        recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
        SA.setCurrentRecording(recName);
        SA.currentDataObj.generateChannelMapFile('40_16x2_FlexLin');
        binaryileName = [SA.currentDataObj.recordingDir filesep 'spikeSorting' filesep 'ch1_32.bin'];
        SA.currentDataObj.convert2Binary(binaryileName);
    end
end

%% upload the data to the kilosort
addpath(genpath('/home/nitzan/Documents/MATLAB/Kilosort-main')) % path to kilosort folder
addpath(genpath('/home/nitzan/Documents/MATLAB/npy-matlab-master')) % for converting to Phy


WA=sleepAnalysis('/media/sil1/Data/Pogona Vitticeps/brainStatesWake.xlsx');
% bRecName = [recName 'b'];
bRecName = 'Animal=PV126,recNames=Hunter23b';
WA.setCurrentRecording(bRecName);

WA.currentDataObj.getKiloSort('/home/nitzan/tempKilosort','tStart',30*60*1000); % this is a new format where
WA.currentDataObj.convertPhySorting2tIc([],2*60*60*1000); % is there manual annotation





%% create the binary file for kilosorting: (only if not yet excisting)
% add another line in the excel, to the binary: sleepNight12b for example,
% and add it in a folder namde spikeSorting in the original recording
% folder. 
Sb=sleepAnalysis('/media/sil2/Data/Lizard/Stellagama/brainStatesSS.xlsx'); 
Sb.setCurrentRecording('Animal=SA07,recNames=sleepNight11');
SA.currentDataObj.generateChannelMapFile('40_16x2_FlexLin'); % 120_32x1_H4_CamNeuro for PV24. all layouts in TSV/electrode layouts 
SA.currentDataObj.convert2Binary([SA.currentDataObj.recordingDir filesep 'spikeSorting' filesep 'ch1_32.bin'],1:32,1:32);

%% upload the data to the kilosort
addpath(genpath('/home/nitzan/Documents/MATLAB/Kilosort-main')) % path to kilosort folder
addpath(genpath('/home/nitzan/Documents/MATLAB/npy-matlab-master')) % for converting to Phy


SA=sleepAnalysis('/media/sil2/Data/Lizard/Stellagama/brainStatesSS.xlsx'); 
Sb.setCurrentRecording('Animal=SA15,recNames=sleepNight5b');
SA.currentDataObj.getKiloSort('/home/nitzan/tempKilosort','tStart',30*60*1000); % this is a new format where
SA.currentDataObj.convertPhySorting2tIc([],2*60*60*1000); % is there manual annotation




%[out,envs] = system(['conda env list'])
%conda('activate','phy2')
%set phy2 environment
%run phy on current analysis object
system(['phy template-gui "' fullfile(SA.currentDataObj.recordingDir,['kiloSortResults_',SA.currentDataObj.recordingName]) filesep 'params.py"'])
%Convert to t-ic - this should run after performeing the phy manual curation analysis
SA.currentDataObj.convertPhySorting2tIc; % 
                      



