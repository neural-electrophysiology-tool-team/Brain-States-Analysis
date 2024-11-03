
%% getStimSham
% this functio is getting the data rellevant for the stimulation and the
% "sham" stimulation for the nights with stimulation. It should save the
% data in the analysis folder. 

function data = getStimSham(SA, t_ch, overwrite)
% SA is an instance of sleep analysis class,with a record currently
% selected
if nargin ==2
    overwrite = 0;
end
  %check if analysis was already done done
    SA.files.stimSham=[SA.currentAnalysisFolder filesep 'stimSham.mat'];
    if exist(SA.files.stimSham,'file') & ~overwrite
        if nargout==1
            data=load(SA.files.stimSham);
        else
            disp('stim sham analysis already exists for this recording');
        end
        return;
    end

    DB=SA.getDelta2BetaRatio;
    AC=SA.getDelta2BetaAC;
    T=SA.getDigitalTriggers;
    firstTrig=T.tTrig{t_ch}(1:8:end-2);
    endStim=T.tTrig{t_ch}(8:8:end)+400;
    stimDuration=(endStim(1)-firstTrig(1));
    pre=50000;
    post=100000;
%     clear StimDB; %change to zeros
%     StimDB = zeros(1,numel(firstTrig));
    for i=1:numel(firstTrig)
        pTmp=find(DB.t_ms>(firstTrig(i)-pre) & DB.t_ms<=(firstTrig(i)+post));
        %StimDB(i,:)=1./DB.bufferedDelta2BetaRatio(pTmp);
        if length(pTmp) ~=150
            pTmp = ones([1,150]);
        end
        StimDB(i,:)=DB.bufferedDelta2BetaRatio(pTmp);
    end

    ts=(DB.t_ms(pTmp)-DB.t_ms(pTmp(1)))/1000;
    meadStimInterval=mean(diff(firstTrig));
    firstTrigSham=(AC.tStartSleep:meadStimInterval:(firstTrig(1)-post))+10000;
    endStimSham=firstTrigSham+max(endStim-firstTrig);
    
%     clear StimDBSham;
    for i=1:numel(firstTrigSham)
        pTmp=find(DB.t_ms>(firstTrigSham(i)-pre) & DB.t_ms<=(firstTrigSham(i)+post));
        %StimDBSham(i,:)=1./DB.bufferedDelta2BetaRatio(pTmp);
        StimDBSham(i,:)=DB.bufferedDelta2BetaRatio(pTmp);
    end

   % save the data
save(SA.files.stimSham,'StimDBSham','ts','StimDB','stimDuration','pre','post')
data.StimDBSham = StimDBSham;
data.ts = ts;
data.StimDB = StimDB;
data.stimDur = stimDuration;
data.pre = pre;
data.post = post;

end

function plotStimSham(SA)
    % SA is an instance of sleep analysis class,with a record currently
    % selected
    stimShamFile=[SA.currentAnalysisFolder filesep  'stimSham.mat'];
    SA.checkFileRecording(stimShamFile,'stim Sham file missing, please first run getStimSham');
    load(stimShamFile); %load data
    
    
    colorLim=[0 600];
    f=figure;
    subplot(4,2,[1:2:6]);imagesc(StimDBSham,colorLim);ylabel('Trial #');title('Sham');hold on;set(gca,'XTick',[]);
    cb=colorbar('Position',[0.47 0.76 0.013 0.17]);ylabel(cb,'\delta/\beta');
    line([pre/1000 pre/1000],ylim,'color','r');
    subplot(4,2,7);plot(ts-pre/1000,nanmean(StimDBSham));xlabel(['Time [s]']);ylabel('Avg.');ylim(colorLim/3);
    line([0 0],ylim,'color','r');
    line([stimDuration/1000 stimDuration/1000],ylim,'color','r');
    subplot(4,2,[2:2:6]);imagesc(StimDB,colorLim);ylabel('Trial #');title('Stim');set(gca,'XTick',[]);
    cb=colorbar('Position',[ 0.91 0.76 0.013 0.17]);ylabel(cb,'\delta/\beta');
    line([pre/1000 pre/1000],ylim,'color','r');
    subplot(4,2,8);plot(ts-pre/1000,nanmean(StimDB));xlabel(['Time [s]']);ylabel('Avg.');ylim(colorLim/3);
    line([0 0],ylim,'color','r');
    line([stimDuration/1000 stimDuration/1000],ylim,'color','r');

%   save?
    fileName=[SA.currentPlotFolder filesep 'stim_sham_activation.jpg'];
    saveas (f, fileName);

    
end

function [spikeRate, spikeRate_t] = getSpikeRate(data, units, tStart, win, meanWin, OL)
    % This function computes the spike rate for multiple specified units (neurons) in a dataset
    % across time windows with specified duration, averaging within a window size (meanWin)
    % with a specified overlap (OL).
    %
    % Inputs:
    %   - data: A structure containing spikes and timings, output from Phy after conversion to .mat
    %   - units: Array of unit IDs (e.g., neuron IDs) for which spike rates are calculated
    %   - tStart: The start time in ms (default = 0)
    %   - win: Duration of time window in ms (default = 2 hours)
    %   - meanWin: Duration of the averaging window in ms (default = 1000 ms)
    %   - OL: Overlap duration between consecutive windows in ms (default = 500 ms)
    %
    % Outputs:
    %   - spikeRate: Matrix of spike rates (spikes per second) with each row corresponding to a unit
    %   - spikeRate_t: Time array corresponding to each spike rate value

    % Assign parameters if needed
    if nargin < 6
        OL = 100; % ms
    end
    if nargin < 5
        meanWin = 1000; % ms
    end
    if nargin < 4
        win = 2 * 60 * 60 * 1000; % 2 hours in ms
    end
    if nargin < 3
        tStart = 0;
    end
    if nargin < 2 
        units = 1; % Default to unit 1 if no units are specified
    end

    % Extract spike data
    spikeClusters = data.spikes_clusters;
    spikeTimes_ms = data.spikes_times / 1000; % Convert to ms if needed

    % Define the end time of the analysis window
    tEnd = tStart + win;

    % Generate time windows
    spikeRate_t = tStart:OL:(tEnd - meanWin);
    numWindows = length(spikeRate_t);
    numUnits = length(units);

    % Pre-allocate spike rate matrix
    spikeRate = zeros(numUnits, numWindows);

    % Loop over each unit
    for u = 1:numUnits
        % Filter spikes for the current unit
        unitSpikes = spikeTimes_ms(spikeClusters == units(u));

        % Calculate spike rate for each window for the current unit
        for i = 1:numWindows
            % Define the start and end of the current window
            winStart = spikeRate_t(i);
            winEnd = winStart + meanWin;

            % Count spikes within the current window
            spikesInWindow = unitSpikes(unitSpikes >= winStart & unitSpikes < winEnd);
            spikeCount = length(spikesInWindow);

            % Calculate the spike rate (spikes per second)
            spikeRate(u, i) = spikeCount / (meanWin / 1000); % Convert ms to s
        end
    end
end

