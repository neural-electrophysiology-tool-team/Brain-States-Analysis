%% getStimMicrophonerig
function [microphoneTriggers]=getStimMicrophoneTrig(obj,templateS, templateF,varargin)

microphoneTriggers = table; 

obj.checkFileRecording;

channels=obj.currentDataObj.analogChannelNames; %get the order of the analog channels
microphoneOrdinaryCh=find(cellfun(@(x) strcmp(x,'ADC2'),channels)); %get the ordinary number of the microphone channel
fs = obj.currentDataObj.samplingFrequencyAnalog(microphoneOrdinaryCh);

%get the timeStamp of the sound given before the reward, using the blockLog
arena = obj.getArenaCSVs(1);
rewardIdx=find(contains(arena.blockLog.Var5,"Reward given") & contains(arena.blockLog.Var5,"playing audio"));
if isempty(rewardIdx)
    return
end
event_frame = obj.getVideoFrames(arena.videoFrames, arena.blockLog.Timestamps(rewardIdx));
event_frame = event_frame - arena.IRFrames(1);
event_oe_time = arena.oeCamTrigs(event_frame); %get the sound event in time (ms)

timeStartS=[]; timeStartF=[]; %find the start of the sound's and fedder's voices
for event=1:length(event_oe_time)
    parseObj = inputParser;
    addParameter(parseObj,'microphoneCh',microphoneOrdinaryCh,@isnumeric); 
    addParameter(parseObj,'tStart',event_oe_time(event)-6000,@isnumeric); %window begins 1 sec prior to the sound event in the blockLog. If 0 begins from the start of the recording
    addParameter(parseObj,'win',12000,@isnumeric); %if 0 uses the whole recording duration
    addParameter(parseObj,'overwrite',0,@isnumeric);
     %event_oe_time(event)-2000,8000       0,obj.currentDataObj.recordingDuration_ms
    addParameter(parseObj,'inputParams',false,@isnumeric);
    parseObj.parse(varargin{:});
    if parseObj.Results.inputParams
        disp(parseObj.Results);
        return;
    end
    
    %evaluate all input parameters in workspace
    for i=1:numel(parseObj.Parameters)
        eval([parseObj.Parameters{i} '=' 'parseObj.Results.(parseObj.Parameters{i});']);
    end
    
    %make parameter structure
    parTrigger=parseObj.Results;
    
    %check if analysis was already done
    obj.files.microphoneTriggers=[obj.currentAnalysisFolder filesep 'microphoneTriggers_ch' num2str(microphoneCh) '.mat'];
    if exist(obj.files.microphoneTriggers,'file') & ~overwrite
        disp('MicrophoneTriggers analysis already exists for this recording');
        microphoneTriggers=load(obj.files.microphoneTriggers);
        return;
    end
    
    [m,m_ms] = obj.currentDataObj.getAnalogData(microphoneCh,tStart,win);
    %[m,m_ms] = obj.currentDataObj.getAnalogData(microphoneCh,0,obj.currentDataObj.recordingDuration_ms);
    m = squeeze(m);
    % m = m - mean(m);
    templateS=templateS(1:20000);
    matchS = convn(m,flip(templateS.'),'same'); %convolution with the sound template
    matchF = convn(m,flip(templateF.'),'same'); %convolution with the feeder template
    %figure; plot(m); hold on; xline(floor(event_oe_time),'k');
    % figure; subplot(3,1,1); plot(m); subplot(3,1,2); plot(template); subplot(3,1,3); plot(match)
    [~, maxIdxS] = max(abs(matchS)); [~, maxIdxF] = max(abs(matchF));
    idxS = round(maxIdxS - length(templateS)/2);
    idxF = round(maxIdxF - length(templateF)/2);
    timeStartS = [timeStartS; m_ms(idxS) + event_oe_time(event) - 6000]; %get the start time
    timeStartF = [timeStartF; m_ms(idxF) + event_oe_time(event) - 6000] ;

    % idxStartS = [idxStartS; idxS + (event_oe_time(event)*fs/1e3) - 120001]; %get the start sample from the recording beginning
    % idxStartF = [idxStartF; idxF + (event_oe_time(event)*fs/1e3) - 120001];

    figure; plot(m_ms',m); hold on; xline(m_ms(120001),'k'); xline(m_ms(idxS),'r'); xline(m_ms(idxF),'g');

    % figure; plot(m_ms',m); hold on; xline(m_ms(120001),'k'); xline(m_ms(idxS),'r'); xline(timeStartF,'g');
    
    % timeStartF =[]; timeStartF = [8156+ event_oe_time(event) - 6000];
     % timeStartF = [8250+ event_oe_time(event) - 6000; 8226+ event_oe_time(event) - 6000]; 
    % timeStartF = [NaN; NaN; NaN; NaN; NaN; NaN; NaN; NaN; NaN; NaN; NaN; NaN; 8147+ event_oe_time(event) - 6000; NaN; NaN; 8358+ event_oe_time(event) - 6000; 8339+ event_oe_time(event) - 6000; 9044+ event_oe_time(event) - 6000];
    % timeStartF = [9170+ event_oe_time(event) - 6000; 8399 + event_oe_time(event) - 6000; 8682 + event_oe_time(event) - 6000; NaN; 7607 + event_oe_time(event) - 6000; 8408 + event_oe_time(event) - 6000; 8425 + event_oe_time(event) - 6000; 8414 + event_oe_time(event) - 6000; 8403 + event_oe_time(event) - 6000; 8420 + event_oe_time(event) - 6000; 8396 + event_oe_time(event) - 6000; 7836 + event_oe_time(event) - 6000; 8375 + event_oe_time(event) - 60000];
    % timeStartS=[]; timeStartF=[];
    % if ~any(m(idxS:idxS+20000) > 7e4) % remove sound when sound wasn't above threshold
    %     timeStartS(timeStartS==(m_ms(idxS) + event_oe_time(event) - 6000))=[];
    % end
    % delete=0;
    % if idxF <= idxS % remove feeder when feeder was detected before sound
    %     % add: find if it doesn't belong to the right event, but the previous. maybe draw it before removing
    % 
    %     % timeStartF(timeStartF==(m_ms(idxF) + event_oe_time(event) - 6000))=[];
    %     delete=1;
    % elseif idxF-idxS < 15000 % remove feeder when feeder was detected as a part from the sound
    %     timeStartF(timeStartF==(m_ms(idxF) + event_oe_time(event) - 6000))=[];
    %     delete=1;
    % end
    % if delete==0
    %     if sum(m(idxF:idxF+10000) >2e4) > 50 % remove feeder when the signal detected is too long than a typical size of a feeder
    %         timeStartF(timeStartF==(m_ms(idxF) + event_oe_time(event) - 6000))=[];
    %     end
    % end
    a=1;
end

a=1;
% microphoneTriggers=table(idxStartS,idxStartF,'VariableNames',{'Sound','Feeder'});
microphoneTriggers=struct();
microphoneTriggers.Sound = timeStartS;
microphoneTriggers.Feeder = timeStartF;
disp('yes')
save(obj.files.microphoneTriggers,'microphoneTriggers','parTrigger');
% 

end