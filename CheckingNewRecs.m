%% check new nights recording:
SA = sleepAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWake.xlsx');
SA.setCurrentRecording('Animal=PV143,recNames=Night41');
SA.getDelta2BetaRatio;
SA.getDelta2BetaAC;
SA.plotDelta2BetaRatio;
SA.plotDelta2BetaSlidingAC(stim=1);
SA.getStimDiodeTrig;
% trig = SA.getStimDiodeTrig
% tTrig = SA.getDigitalTriggers
%% check new hunter rec:
WA = wakeAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWake.xlsx');
WA.setCurrentRecording('Animal=PV106,recNames=Hunter17');
WA.getArenaCSVs(1);
WA.getBeta2GammaRatio;
WA.plotTrailsBG;
A = WA.getArenaCSVs;


%% run on all and do...?
SA = wakeAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWake.xlsx');

Animals = ["PV106"; "PV143"; "PV153"];
% Animals="PV106";
% plot the D/B together

% figure;
% tiledlayout('flow'); 
% ax1 = nexttile;
for i = 1:length(Animals)
    for j = 1:42
        try  
            txt = evalc("SA.setCurrentRecording(sprintf('Animal=%s,recNames=Night%d',Animals{i},j));");
            if contains(txt,'Selected recording/s were not found')
                continue
            end
                       
            % SA.setCurrentRecording(sprintf('Animal=%s,recNames=Night%d',Animals{i},j));
            % if SA.currentDataObj==[]
            %     continue
            % end

            % SA.plotDelta2BetaRatio;
            % SA.plotDelta2BetaSlidingAC;
            if SA.recTable.Mani(SA.currentPRec) ==4 & contains(SA.recTable.Remarks(SA.currentPRec),'blue')
                plotStimSham(SA)
                % SA.plotDelta2BetaSlidingAC;

                % ax1 = nexttile; % left tile
                % getStimSham(SA,[],1)
                % ax1 = plotStimSham(SA);
                % SA.getDelta2BetaRatio;
                % SA.getDelta2BetaAC;
               
                % SA.plotDelta2BetaRatio('h',ax1,saveFigure = 0);
                % title(ax1,sprintf('%s, stimType = %s',Animals{i},SA.recTable.Remarks{SA.currentPRec}))

                % ax2 = nexttile;
                % SA.plotDelta2BetaSlidingAC('h',ax2,'saveFigures',0);
            %     [lfp,lfp_ms] = SA.currentDataObj.getData(SA.recTable.defaulLFPCh(SA.currentPRec),1000*60*10,1000);
            %     plot(lfp_ms,squeeze(lfp))

            end
        catch ME
            fprintf('Error with %s Night%d: %s\n', Animals(i), j, ME.message);

        end
    end
end


%% run on all and do...?
SA = wakeAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWake.xlsx');

Animals = ["PV106"; "PV143"; "PV153"];
% Animals="PV106";
% plot the D/B together
figure;
tiledlayout('flow'); 
% ax1 = nexttile;
for i = 1:length(Animals)
    for j = 1:39
        try  
            txt = evalc("SA.setCurrentRecording(sprintf('Animal=%s,recNames=Night%d',Animals{i},j));");
            if contains(txt,'Selected recording/s were not found')
                continue
            end
                       
            % SA.setCurrentRecording(sprintf('Animal=%s,recNames=Night%d',Animals{i},j));
            % if SA.currentDataObj==[]
            %     continue
            % end

            % SA.plotDelta2BetaRatio;
            % SA.plotDelta2BetaSlidingAC;
            if SA.recTable.Mani(SA.currentPRec) ==4 && contains(SA.recTable.Remarks(SA.currentPRec),'red') && SA.recTable.spikes(SA.currentPRec)==1
                % SA.plotDelta2BetaSlidingAC;
                % SA.getEyeMovements('startTime',2*60*60,'endTime',12*60*60)
                % SA.getRespirationMovements('startTime',2*60*60,'endTime',12*60*60)
                % ax1 = nexttile; % left tile
                % getStimSham(SA,15)
                % ax1 = plotStimSham(SA);
                % SA.getDelta2BetaRatio;
                % SA.getDelta2BetaAC;
                % SA.plotDelta2BetaRatio('h',ax1,saveFigure = 0,stim=1);clim(ax1,[0 450])
                % title(ax1,sprintf('%s, stimType = %s',Animals{i},SA.recTable.Remarks{SA.currentPRec}))

                % ax2 = nexttile;
                % SA.plotDelta2BetaSlidingAC('h',ax2,'saveFigures',0,stim=1);clim(ax2,[0 450])
                %  [lfp,lfp_ms] = SA.currentDataObj.getData(SA.recTable.defaulLFPCh(SA.currentPRec),1000*60*10,1000);
                %  plot(lfp_ms,squeeze(lfp))

                % rec to binary: kilosort
                if ~contain(SA.recTable.Animal(SA.currentPRec),'PV153') 
                    SA.currentDataObj.generateChannelMapFile('40_16x2_FlexLin');
                end
                binaryileName = [SA.currentDataObj.recordingDir filesep 'spikeSorting' filesep 'ch1_32.bin'];
                SA.currentDataObj.convert2Binary(binaryileName);
                
                
            end
        catch ME
            fprintf('Error with %s Night%d: %s\n', Animals(i), j, ME.message);

        end
    end
end




%% print
analysisFolder = '/media/sil1/Data/Nitzan/Light Manipulation paper/NitzanAnalysisFiles';
fileName=[analysisFolder filesep 'ExtRedNights'];
print(fileName,'-djpeg');


%% updating recs: 
%% 
% first, I want to see all the plots according to the stimulations.
analysisFolder = '/media/sil1/Data/Nitzan/Light Manipulation paper/NitzanAnalysisFiles';
load([analysisFolder filesep 'stimTableAll.mat'])

%%
recList=[];
for i = 1:height(stimTable)

    recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
    SA.setCurrentRecording(recName);
    stimData = SA.getStimTriggers;
    
    if isempty(stimData)
        recList = [recList; recName];
    end
end
%% get the trigger fixed

% SA = wakeAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWake.xlsx');
% SA.setCurrentRecording('Animal=PV106,recNames=Night18');
vidTrigs = SA.currentDataObj.getCamerasTrigger(7);
%

figure; plot(diff(vidTrigs),'.');title('Original triggers')

% Step 1: find small diffs
dt = diff(vidTrigs);
thresh = 15; 
badIdx = find(dt < thresh);   % indices of bad diffs

% Step 2: detect cluster starts
isStart = [false, diff(badIdx) == 1]; 

% Step 3: keep the indx if it have 1 before and after
keeps = find(diff(isStart)==0);
keeps2 = reshape(keeps,2,length(keeps)/2);
real_keeps = keeps2(2,:);
isStart(real_keeps)=0;

removeIdx = badIdx(isStart);  

% Step 4: clean timings
t_clean = vidTrigs;
t_clean(removeIdx) = [];

figure; plot(diff(t_clean),'.');title('after cleanup')



%% run Eye movement analysis #1
SA.setCurrentRecording('Animal=PV106,recNames=Night10');
SA.plotDelta2BetaSlidingAC(stim=1);
SA.getEyeMovements('startTime',3*60*60,'endTime',10*60*60,'overwrite',1,'opticFlowNoiseThreshold',0.002,'loadInitialConditions',0);
trigs = SA.getStimTriggers;
% SA.getEyeMovements('startTime',trigs(1)/1000,'endTime',(trigs(end)/1000)+(30*60),'minTrackingPoints',25,'overwrite',1, 'opticFlowNoiseThreshold',0.002);
SA.getCameraTriggers;
SA.getSyncedDBEyeMovements('digitalVideoSyncCh',7,'useRobustFloatingAvg',0,'overwrite',1);
SA.plotSyncedDBEyeMovementsRaster
%% run Eye movement analysis #2
SA.setCurrentRecording('Animal=PV106,recNames=Night14');
% SA.plotDelta2BetaSlidingAC(stim=1);
SA.getEyeMovements('startTime',3*60*60,'endTime',10*60*60,'overwrite',1,'opticFlowNoiseThreshold',0.002,'loadInitialConditions',0);
trigs = SA.getStimTriggers;
% SA.getEyeMovements('startTime',trigs(1)/1000,'endTime',(trigs(end)/1000)+(30*60),'minTrackingPoints',25,'overwrite',1, 'opticFlowNoiseThreshold',0.002);
SA.getCameraTriggers;
SA.getSyncedDBEyeMovements('digitalVideoSyncCh',7,'useRobustFloatingAvg',0,'overwrite',1);
SA.plotSyncedDBEyeMovementsRaster



%% run behavioral:
SA.setCurrentRecording('Animal=PV106,recNames=Night16');
SA.plotDelta2BetaSlidingAC;
SA.getCameraTriggers;
SA.getStimTriggers;
SA.getSyncedDBEyeMovements;
SA.plotSyncedDBEyeMovementsRaster;


%% check for spiking activity in recs:

%% spike sorting - revision
SA=sleepAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWake.xlsx');
analysisFolder = '/media/sil1/Data/Nitzan/Light Manipulation paper/NitzanAnalysisFiles';
% load([analysisFolder filesep 'stimTableAll.mat'])

recInd = (contains(SA.recTable.Remarks,'white')| contains(SA.recTable.Remarks,'DayTime')) & SA.recTable.Mani >0;% SA.recTable.spikes ==1
animals =SA.recTable.Animal(recInd);recnames = SA.recTable.recNames(recInd);
subset = SA.recTable(recInd,["Animal","recNames","Remarks"]);
recList = cellfun(@(x,y) ['Animal=' x ',recNames=' y], animals, recnames, 'UniformOutput', false);
% tStart=1000*60*60*5;
figure; tiledlayout('flow'); 
for i = 1:length(recList)
    SA.setCurrentRecording(recList{i});
    % [lfp,lfp_t]=SA.currentDataObj.getData(SA.recTable.defaulLFPCh(SA.currentPRec),tStart,40000);
    % nexttile;
    % plot(lfp_t,squeeze(lfp)); title(recList{i})
    if ~contains(recList{i},'PV153')
        SA.currentDataObj.generateChannelMapFile('40_16x2_FlexLin'); % 120_32x1_H4_CamNeuro for PV24. all layouts in TSV/electrode layouts
    end
    binaryileName = [SA.currentDataObj.recordingDir filesep 'spikeSorting' filesep 'ch1_32.bin'];
    SA.currentDataObj.convert2Binary(binaryileName);
end


%% find best whits - revision
SA=sleepAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWake.xlsx');
analysisFolder = '/media/sil1/Data/Nitzan/Light Manipulation paper/NitzanAnalysisFiles';
% load([analysisFolder filesep 'stimTableAll.mat'])

recInd = (contains(SA.recTable.Remarks,'white')| contains(SA.recTable.Remarks,'DayTime')) & SA.recTable.Mani >0;% SA.recTable.spikes ==1
animals =SA.recTable.Animal(recInd);recnames = SA.recTable.recNames(recInd);
subset2 = SA.recTable(recInd,["Animal","recNames","Remarks","spikes","folder"]);
recList = cellfun(@(x,y) ['Animal=' x ',recNames=' y], animals, recnames, 'UniformOutput', false);
%% tStart=1000*60*60*5;
figure; tiledlayout('flow'); 
for i = 1:length(recList)
    SA.setCurrentRecording(recList{i});
    ax1 = nexttile; % left tile
    SA.plotDelta2BetaRatio('h',ax1,saveFigure = 0,stim=1);clim(ax1,[0 450])
    title(ax1,sprintf('%s, stimType = %s',recList{i},SA.recTable.Remarks{SA.currentPRec}))

    ax2 = nexttile;
    SA.plotDelta2BetaSlidingAC('h',ax2,'saveFigures',0,stim=1);
    % clim(ax2,[0 450]);
end

%%
indskilo = subset2.kilosort==0&subset2.spikes==1&strcmp(subset2.Animal,'PV153');
recsFolders = subset2.folder(indskilo);
folders = string(recsFolders);          % convert to string array
X = 17;
folders = "W:/" + extractAfter(folders, X) + "/spikeSorting/ch1_32.bin";
folders = replace(folders,'/','\');