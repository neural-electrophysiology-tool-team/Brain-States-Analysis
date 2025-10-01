%% check new nights recording:
% SA = wakeAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWake.xlsx');
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

% Animals = ["PV106"; "PV143"; "PV153"];
Animals="PV106";
% plot the D/B together
figure;
tiledlayout('flow'); 
% ax1 = nexttile;
for i = 1:length(Animals)
    for j = 7:35
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
            if SA.recTable.Mani(SA.currentPRec) ==4
                SA.plotDelta2BetaSlidingAC;
                SA.getEyeMovements('startTime',2*60*60,'endTime',12*60*60)
                SA.getRespirationMovements('startTime',2*60*60,'endTime',12*60*60)
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
            if SA.recTable.Mani(SA.currentPRec) ==4 & contains(SA.recTable.Remarks(SA.currentPRec),'daylight')
                % SA.plotDelta2BetaSlidingAC;
                % SA.getEyeMovements('startTime',2*60*60,'endTime',12*60*60)
                % SA.getRespirationMovements('startTime',2*60*60,'endTime',12*60*60)
                ax1 = nexttile; % left tile
                getStimSham(SA,15)
                ax1 = plotStimSham(SA);
                % SA.getDelta2BetaRatio;
                % SA.getDelta2BetaAC;
               
                % SA.plotDelta2BetaRatio('h',ax1,saveFigure = 0);
                title(ax1,sprintf('%s, stimType = %s',Animals{i},SA.recTable.Remarks{SA.currentPRec}))

                ax2 = nexttile;
                SA.plotDelta2BetaSlidingAC('h',ax2,'saveFigures',0);
            %     [lfp,lfp_ms] = SA.currentDataObj.getData(SA.recTable.defaulLFPCh(SA.currentPRec),1000*60*10,1000);
            %     plot(lfp_ms,squeeze(lfp))

            end
        catch ME
            fprintf('Error with %s Night%d: %s\n', Animals(i), j, ME.message);

        end
    end
end




%% print
analysisFolder = '/media/sil1/Data/Nitzan/Light Manipulation paper/NitzanAnalysisFiles';
fileName=[analysisFolder filesep 'PV106Nights'];
print(fileName,'-djpeg');


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



%% run Eye movement analysis
SA.setCurrentRecording('Animal=PV106,recNames=Night18');
SA.plotDelta2BetaSlidingAC;
SA.getEyeMovements('startTime',2*60*60,'endTime',12*60*60);