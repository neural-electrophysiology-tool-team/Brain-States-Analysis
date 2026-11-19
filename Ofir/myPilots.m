%% A new table - 22.12.25 - all rewards

% reward = table([], [], [], [], [], [], [], [], [], [], [], [], 'VariableNames', {'Animal','Record','Idx','timeSound','timeFeeder','mpSound','mpFeeder','signalSound','signalFeeder', 'typeSound', 'typeTouch','Notes'});
reward = table([], [], [], [], 'VariableNames', {'Animal','Record','Idx','typeSound'});
records=SA.recTable;
animals=["PV143","PV106","PV153"];
idxTrial=ismember(records.Animal,animals);
records_filtered=records(idxTrial, :);
% winPrev=6000; win=20000;

for i=1:size(records_filtered,1) 
    i
    curRecName=string(records_filtered.recNames(i));
    if contains (curRecName,"Hunter")
        curAnimal=string(records_filtered.Animal(i));
        row=sprintf('Animal=%s,recNames=%s', curAnimal, curRecName);
        current=SA.setCurrentRecording(row);
        arena = SA.getArenaCSVs(1);
        rewardIdx = find(contains(arena.blockLog.Var5,"Reward given") & contains(arena.blockLog.Var5,"playing audio"));

        if ~isempty(rewardIdx)
            % analog_channels=SA.currentDataObj.analogChannelNames; %get the order of the analog channels
            % microphoneOrdinaryCh=find(cellfun(@(x) strcmp(x,'ADC2'),analog_channels)); %get the ordinary number of the microphone channel
            % % fs_analog = obj.currentDataObj.samplingFrequencyAnalog(microphoneOrdinaryCh);
            % 
            % LFPchannel=records_filtered.defaulLFPCh(i); %get the order of the analog channels
            % fs_LFP = SA.currentDataObj.samplingFrequency(LFPchannel);
    
            microphoneTriggers=getStimMicrophoneTrig(SA,mean_template_sound_aligned,mean_template_feeder_aligned);
            % s=microphoneTriggers.microphoneTriggers.Sound;
            % f=microphoneTriggers.microphoneTriggers.Feeder;
            % mp_sound=[]; mp_feeder=[]; signal_sound=[]; signal_feeder=[];
            % for s_num=1:length(s)
            %     [m_s,~] = SA.currentDataObj.getAnalogData(microphoneOrdinaryCh,s(s_num)-winPrev,win);
            %     m_s = squeeze(m_s);
            %     mp_sound=[mp_sound,m_s];
            %     sig_s = SA.currentDataObj.getData(LFPchannel, s(s_num)-winPrev, win);
            %     sig_s = squeeze(sig_s);
            %     signal_sound=[signal_sound,sig_s];
            % end
            % 
            % for f_num=1:length(f)
            %     [m_f,~] = SA.currentDataObj.getAnalogData(microphoneOrdinaryCh,f(f_num)-winPrev,win);
            %     m_f = squeeze(m_f);
            %     mp_feeder=[mp_feeder,m_f];
            %     sig_f = SA.currentDataObj.getData(LFPchannel, f(f_num)-winPrev, win);
            %     sig_f = squeeze(sig_f);
            %     signal_feeder=[signal_feeder,sig_f];
            % end

            touchIdx = find(contains(arena.blockLog.Var5,"touch event"));
            % touchMat = [];
            nSound = numel(rewardIdx);
            soundType = strings(nSound,1);   % "hit" or "sound_only"
            soundUsed = false(nSound,1);     % prevent double matching

            if ~isempty(touchIdx)
                % touchMat = arena.screenTouch.is_hit;

                maxDelay = 3;   % seconds
                nTouch = numel(touchIdx);
                

                for t = 1:nTouch
                    % find sounds AFTER this touch
                    idx = find(rewardIdx > touchIdx(t) & ~soundUsed, 1, 'first');

                    t1 = arena.blockLog.DateTime(touchIdx(t));
                    t2 = arena.blockLog.DateTime(idx);
                    dt1 = datetime(t1,'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
                    dt2 = datetime(t2,'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
                    
                    delta = seconds(dt2 - dt1);
                
                    if ~isempty(idx) && delta <= maxDelay
                        soundType(idx) = "touch";
                        soundUsed(idx) = true;
                    end
                end
            end

            soundType(soundType == "") = "sound_only";
            a=1;
            
            if rewardIdx > 1
                for id = 1:length(rewardIdx)
                    newRow = table(curAnimal, curRecName, rewardIdx(id), soundType(id), ...
                       'VariableNames', reward.Properties.VariableNames);
                    reward = [reward; newRow];
                    % reward = [reward; {curAnimal, curRecName, rewardIdx, s, f, mp_sound, mp_feeder, signal_sound, signal_feeder, soundType, touchMat}];
                end
            else
                newRow = table(curAnimal, curRecName, rewardIdx, soundType, ...
                       'VariableNames', reward.Properties.VariableNames);
                reward = [reward; newRow];
            end
        end
    end
end

%% Find touch screen
% reward = table([], [], [], 'VariableNames', {'animal','record','idx'});
records=SA.recTable;
animals=["PV153","PV106","PV143"]; %"PV143","PV106",
idx=ismember(records.Animal,animals);
records_filtered=records(idx, :);
for i=1:size(records_filtered,1)
    curRecName=string(records_filtered.recNames(i));
    if contains (curRecName,"Hunter")
        curAnimal=string(records_filtered.Animal(i));
        row=sprintf('Animal=%s,recNames=%s', curAnimal, curRecName);
        current=SA.setCurrentRecording(row);
        arena = SA.getArenaCSVs(1);
        rewardIdx = find(contains(arena.blockLog.Var5,"Reward given") & contains(arena.blockLog.Var5,"playing audio"));
        if ~isempty(rewardIdx)
            reward = [reward; {curAnimal, curRecName, rewardIdx}];
        end
        % touchIdx = find(contains(arena.blockLog.Var5,"touch event"));
        % if ~isempty(touchIdx)
        %     touch = [touch; {curAnimal, curRecName, touchIdx}];
        % end

    end
end


%% Find touch with sound
touch_sound = table([], [], [], [], 'VariableNames', {'animal','record','timeS','timeF'});
for i=1:size(touch,1)
    i
    curAnimal=string(touch.animal(i));
    curRecName=string(touch.record(i));
    row=sprintf('Animal=%s,recNames=%s', curAnimal, curRecName);
    current=SA.setCurrentRecording(row);
    % microphoneTriggers=getStimMicrophoneTrig(SA,mean_template_sound_aligned,mean_template_feeder_aligned);
    % if ~isempty(microphoneTriggers)
    %     touch_sound = [touch_sound; {curAnimal, curRecName, microphoneTriggers.Sound, microphoneTriggers.Feeder}];
    % end
    arena = SA.getArenaCSVs(1);
    hitIdx = find(contains(arena.blockLog.Var5,"Reward given") & contains(arena.blockLog.Var5,"playing audio"));
    if ~isempty(hitIdx)
        touch_sound = [touch_sound; {curAnimal, curRecName, microphoneTriggers.microphoneTriggers.Sound, microphoneTriggers.microphoneTriggers.Feeder}];
    end
    a=1;
end

%% Get mp vectors for touch with sound
row=sprintf('Animal=%s,recNames=%s', "PV106", "Hunter38");
current=SA.setCurrentRecording(row);
s=cell2mat(touch_sound.timeS(2));
f=cell2mat(touch_sound.timeF(2));
touchIdx=[24;36;67];
% event_frame = SA.getVideoFrames(arena.videoFrames, arena.blockLog.Timestamps(touchIdx));
% event_frame = event_frame - arena.IRFrames(1);
% event_oe_timeTouch= arena.oeCamTrigs(event_frame); %get the sound event in time (ms)
% m_s_all=[]; m_f_all=[];
for sID=1:length(s)
[m_s,m_ms_s] = SA.currentDataObj.getAnalogData(5,s(sID)-6000,12000);
m_s=squeeze(m_s);
m_s_all=[m_s_all,m_s];
end
for fID=1:length(f)
[m_f,m_ms_f] = SA.currentDataObj.getAnalogData(5,f(fID)-6000,12000);
m_f=squeeze(m_f);
m_f_all=[m_f_all,m_f];
end

%% Collect and inspect rapidly the data
% records=SA.recTable;
% animals=["PV153","PV106","PV143"]; %"PV143","PV106",
% idx=ismember(records.Animal,animals);
% records_filtered=records(idx, :);
% mpRec=[];
% v=53:129; %PV153:38:108;
% v(ismember(v,[57,106,110,117:120,122,124]))=[]; %strange records in PV153 [52,68,69,80,89,99,108], strange records in PV106 [57,106,110,117:120,122,124]
for i=24:25%:size(rewardNewWOtouch,1) %38:108, 1:size(records_filtered,1), 
    % PV153 skip: hunter15 (52),hunter31 (68),hunter32(69), hunter45(80), hunter55(89), hunter65(99), hunter74(108)
    % PV106 skip: hunter5 (57), hunter54 (106), hunter58 (110), hunter65-70
    % (117-120), hunter70 (122), hunter72 (124)
    curAnimal=rewardNewWOtouch.animal(i);
    curRecName=rewardNewWOtouch.record(i);
     % curRecName=string(records_filtered.recNames(i));
    % curRecName
    % if contains (curRecName,"Hunter")
        % curAnimal=string(records_filtered.Animal(i));
        row=sprintf('Animal=%s,recNames=%s', curAnimal, curRecName);
        row
        % row=sprintf('Animal=%s,recNames=%s', "PV106", "Hunter68");
        current=SA.setCurrentRecording(row);
        microphoneTriggers=getStimMicrophoneTrig(SA,mean_template_sound_aligned,mean_template_feeder_aligned);
        timeSeriesViewer(SA.currentDataObj)
        s=cell2mat(rewardNewWOtouch.Sound(i))
        f=cell2mat(rewardNewWOtouch.Feeder(i))
        a=1;
end
        % arena = SA.getArenaCSVs(1);
        % log=arena.blockLog;
        % find(contains(arena.blockLog.Var5,"touch event"))
        % 
        % channels=SA.currentDataObj.analogChannelNames; %get the order of the analog channels
        % microphoneOrdinaryCh=find(cellfun(@(x) strcmp(x,"ADC2"),channels)); %get the ordinary number of the microphone channel
        % SA.files.microphoneTriggers=[SA.currentAnalysisFolder filesep 'microphoneTriggers_ch' num2str(microphoneOrdinaryCh) '.mat'];
        % if exist(SA.files.microphoneTriggers,'file') 
        %     microphoneTriggers=load(SA.files.microphoneTriggers);
        %     if (~isempty(microphoneTriggers.microphoneTriggers.Sound) && ~isempty(microphoneTriggers.microphoneTriggers.Feeder))
        %         mpRec=[mpRec,curRecName];
        %     end
        % end
      
    % end
% end

%% DONE: Create the sound template for getStimMicrophoneTrig
t1=template_PV153_Hunter42; t2=template_PV153_Hunter44; t3=template_PV153_Hunter51_feeder1;
t4=template_PV153_Hunter52; t5=template_PV153_Hunter54; t6=template_PV153_Hunter90;
maxlen = max([numel(t1) numel(t2) numel(t3) numel(t4) numel(t5) numel(t6)]);
M = nan(6, maxlen);

M(1,1:numel(t1)) = t1;
M(2,1:numel(t2)) = t2;
M(3,1:numel(t3)) = t3;
M(4,1:numel(t4)) = t4;
M(5,1:numel(t5)) = t5;
M(6,1:numel(t6)) = t6;
mean_template_sound = mean(M,1,'omitnan');
figure; plot(M.'); hold on; plot(mean_template_sound,'k','LineWidth',2)

ref=M(1,:);
[acor2, lag2] = xcorr(t2, ref);
[~, I2] = max(acor2);
shift2 = lag2(I2);

[acor3, lag3] = xcorr(t3, ref);
[~, I3] = max(acor3);
shift3 = lag3(I3);

[acor4, lag4] = xcorr(t4, ref);
[~, I4] = max(acor4);
shift4 = lag4(I4);

[acor5, lag5] = xcorr(t5, ref);
[~, I5] = max(acor5);
shift5 = lag5(I5);

[acor6, lag6] = xcorr(t6, ref);
[~, I6] = max(acor6);
shift6 = lag6(I6);

t2_aligned = circshift(t2, -shift2);
t3_aligned = circshift(t3, -shift3);
t4_aligned = circshift(t4, -shift4);
t5_aligned = circshift(t5, -shift5);
t6_aligned = circshift(t6, -shift6);

maxlen = max([numel(ref) numel(t2_aligned) numel(t3_aligned) numel(t4_aligned) numel(t5_aligned) numel(t6_aligned)]);
M = nan(6, maxlen);

M(1,1:numel(t1)) = t1;
M(2,1:numel(t2_aligned)) = t2_aligned;
M(3,1:numel(t3_aligned)) = t3_aligned;
M(4,1:numel(t4_aligned)) = t4_aligned;
M(5,1:numel(t5_aligned)) = t5_aligned;
M(6,1:numel(t6_aligned)) = t6_aligned;
mean_template_sound_aligned = mean(M,1,'omitnan');
figure; plot(ref); hold on; plot(t2_aligned); plot(t3_aligned); plot(t4_aligned); plot(t5_aligned); plot(t6_aligned); plot(mean_template_sound_aligned,'k','LineWidth',2)

%%
t1=template_PV228_Wake23_1; t2=template_PV228_Wake23_2; t3=template_PV228_Wake23_3; t4=template_PV228_Wake23_4; t5=template_PV228_Wake23_5;
t6=template_PV167_Wake5_1; t7=template_PV167_Wake5_2; t8=template_PV167_Wake5_3; t9=template_PV167_Wake5_4; 
t10=template_PV228_Wake36_1; t11=template_PV228_Wake36_2; t12=template_PV228_Wake36_3; t13=template_PV228_Wake36_4; t14=template_PV228_Wake36_5;

maxlen = max([numel(t1) numel(t2) numel(t3) numel(t4) numel(t5) numel(t6) numel(t7) numel(t8) numel(t9) numel(t10) numel(t11) numel(t12) numel(t13) numel(t14)]);
M = nan(14, maxlen);

M(1,1:numel(t1)) = t1;
M(2,1:numel(t2)) = t2;
M(3,1:numel(t3)) = t3;
M(4,1:numel(t4)) = t4;
M(5,1:numel(t5)) = t5;
M(6,1:numel(t6)) = t6;
M(7,1:numel(t7)) = t7;
M(8,1:numel(t8)) = t8;
M(9,1:numel(t9)) = t9;
M(10,1:numel(t10)) = t10;
M(11,1:numel(t11)) = t11;
M(12,1:numel(t12)) = t12;
M(13,1:numel(t13)) = t13;
M(14,1:numel(t14)) = t14;


% mean_template_sound_4 = mean(M,1,'omitnan');
% figure; plot(M.'); hold on; plot(mean_template_sound_4,'k','LineWidth',2)
%%
t1=template_PV228_Wake36_1; t2=template_PV228_Wake36_2; t3=template_PV228_Wake36_3; t4=template_PV228_Wake36_4; t5=template_PV228_Wake36_5;
t6=again_template_PV167_Wake5_1; t7=again_template_PV167_Wake5_2; t8=again_template_PV167_Wake5_3; t9=again_template_PV167_Wake5_4; 

ref=t1;
[acor2, lag2] = xcorr(t2, ref);
[~, I2] = max(acor2);
shift2 = lag2(I2);

[acor3, lag3] = xcorr(t3, ref);
[~, I3] = max(acor3);
shift3 = lag3(I3);

[acor4, lag4] = xcorr(t4, ref);
[~, I4] = max(acor4);
shift4 = lag4(I4);

[acor5, lag5] = xcorr(t5, ref);
[~, I5] = max(acor5);
shift5 = lag4(I5);

[acor6, lag6] = xcorr(t6, ref);
[~, I6] = max(acor6);
shift6 = lag4(I6);

[acor7, lag7] = xcorr(t7, ref);
[~, I7] = max(acor7);
shift7 = lag4(I7);

[acor8, lag8] = xcorr(t8, ref);
[~, I8] = max(acor8);
shift8 = lag4(I8);

[acor9, lag9] = xcorr(t9, ref);
[~, I9] = max(acor9);
shift9 = lag4(I9);
% 
% [acor10, lag10] = xcorr(t10, ref);
% [~, I10] = max(acor10);
% shift10 = lag4(I10);
% 
% [acor11, lag11] = xcorr(t11, ref);
% [~, I11] = max(acor11);
% shift11 = lag4(I11);
% 
% [acor12, lag12] = xcorr(t12, ref);
% [~, I12] = max(acor12);
% shift12 = lag4(I12);
% 
% [acor13, lag13] = xcorr(t13, ref);
% [~, I13] = max(acor13);
% shift13 = lag4(I13);
% 
% [acor14, lag14] = xcorr(t14, ref);
% [~, I14] = max(acor14);
% shift14 = lag4(I14);

t2_aligned = circshift(t2, -shift2);
t3_aligned = circshift(t3, -shift3);
t4_aligned = circshift(t4, -shift4);
t5_aligned = circshift(t5, -shift5);
t6_aligned = circshift(t6, -shift6);
t7_aligned = circshift(t7, -shift7);
t8_aligned = circshift(t8, -shift8);
t9_aligned = circshift(t9, -shift9);
% t10_aligned = circshift(t10, -shift10);
% t11_aligned = circshift(t11, -shift11);
% t12_aligned = circshift(t12, -shift12);
% t13_aligned = circshift(t13, -shift13);
% t14_aligned = circshift(t14, -shift14);


maxlen = max([numel(ref) numel(t2_aligned) numel(t3_aligned) numel(t4_aligned) numel(t5_aligned) numel(t6_aligned) numel(t7_aligned) numel(t8_aligned) numel(t9_aligned)]);% numel(t10_aligned) numel(t11_aligned) numel(t12_aligned) numel(t13_aligned) numel(t14_aligned)]);
M = nan(9, maxlen);

M(1,1:numel(t1)) = t1;
M(2,1:numel(t2_aligned)) = t2_aligned;
M(3,1:numel(t3_aligned)) = t3_aligned;
M(4,1:numel(t4_aligned)) = t4_aligned;
M(5,1:numel(t5_aligned)) = t5_aligned;
M(6,1:numel(t6_aligned)) = t6_aligned;
M(7,1:numel(t7_aligned)) = t7_aligned;
M(8,1:numel(t8_aligned)) = t8_aligned;
M(9,1:numel(t9_aligned)) = t9_aligned;
% M(10,1:numel(t10_aligned)) = t10_aligned;
% M(11,1:numel(t11_aligned)) = t11_aligned;
% M(12,1:numel(t12_aligned)) = t12_aligned;
% M(13,1:numel(t13_aligned)) = t13_aligned;
% M(14,1:numel(t14_aligned)) = t14_aligned;

mean_template_sound_aligned = mean(M,1,'omitnan');
% figure; plot(ref); hold on; plot(t2_aligned); plot(t3_aligned); plot(t4_aligned); plot(t5_aligned); plot(t6_aligned); plot(t7_aligned); plot(t8_aligned);  plot(t9_aligned); plot(t10_aligned); plot(t11_aligned); plot(t12_aligned); plot(t13_aligned); plot(t14_aligned);
% plot(mean_template_sound_aligned,'k','LineWidth',2) 
mean_template_sound_aligned_7= mean_template_sound_aligned-mean(mean_template_sound_aligned);
% mean_template_sound_aligned_3 = mean_template_sound_aligned_2(1:21730);
%% DONE: Create an audio from the mp signal
% Example audio (you can replace this with your own)
m2audio= m / max(abs(m));   % scale to [-1, 1]
audiowrite('PV143Wake10.WAV', m2audio, fs);
[y, Fs] = audioread('PV143Wake10.wav');  


%% Get the mp vectors for mpRec matrix
records=SA.recTable;
% animals=["PV153"]; %"PV143","PV106",
% curAnimal=animals;
% m_s_all=[]; m_f_all=[];
for r=1:size(rewardNew,1)
    r
    curAnimal=rewardNew.animal(r);
    curRecName=rewardNew.record(r);
    row=sprintf('Animal=%s,recNames=%s', curAnimal, curRecName);
    current=SA.setCurrentRecording(row);
    channels=SA.currentDataObj.analogChannelNames; %get the order of the analog channels
    microphoneOrdinaryCh=find(cellfun(@(x) strcmp(x,"ADC2"),channels)); %get the ordinary number of the microphone channel
    fs = SA.currentDataObj.samplingFrequencyAnalog(microphoneOrdinaryCh);
    microphoneTriggers=getStimMicrophoneTrig(SA,mean_template_sound_aligned,mean_template_feeder_aligned);
    s=microphoneTriggers.microphoneTriggers.Sound;
    m_s_Record=[];
    for s_num=1:length(s)
        [m_s,m_ms_s] = SA.currentDataObj.getAnalogData(5,s(s_num)-6000,12000);
        m_s = squeeze(m_s);
        % m_s = m_s - mean(m_s);
        % m_s_all=[m_s_all,m_s];
        m_s_Record=[m_s_Record,m_s];
        rewardNew.mpSound{r} = m_s_Record;
    end
    f=microphoneTriggers.microphoneTriggers.Feeder;
    m_f_Record=[];
    for f_num=1:length(f)
        [m_f,m_ms_f] = SA.currentDataObj.getAnalogData(5,f(f_num)-6000,12000);
        m_f = squeeze(m_f);
        % m_f = m_f - mean(m_f);
        % m_f_all=[m_f_all,m_f];
         m_f_Record=[m_f_Record,m_f];
        rewardNew.mpFeeder{r} = m_f_Record;
    end
    a=1;
end

%% Get the mp vectors sorted by feeder-sound intehttps://www.nature.com/articles/nature08010rval for mpRec matrix
% rows = [15 4 23];
% others = setdiff(1:height(rewardNewWOtouch), rows, 'stable');
% rewardNewWOtouch = [rewardNewWOtouch(others, :) ; rewardNewWOtouch(rows, :)];

diffs = []; diffsIDX=[];
iterate=0; soundOnly=[];
for r=1:size(rewardNewWOtouch)
    r
    curAnimal=rewardNewWOtouch.animal(r);
    curRecName=rewardNewWOtouch.record(r);
    row=sprintf('Animal=%s,recNames=%s', curAnimal, curRecName);
    current=SA.setCurrentRecording(row);

    microphoneTriggers=getStimMicrophoneTrig(SA,mean_template_sound_aligned,mean_template_feeder_aligned);
    s=microphoneTriggers.microphoneTriggers.Sound;
    f=microphoneTriggers.microphoneTriggers.Feeder;

    % if isempty(f)
    %     iterate=iterate+1;
    %     soundOnly=[soundOnly,iterate];
    %     continue;
    % end

    if length(s)>length(f)
        iterate=iterate+2;
        if f >= s(1) && f <= s(2)
            diffs = [diffs, f-s(1)];
            diffsIDX = [diffsIDX, iterate-1];
            % soundOnly=[soundOnly,iterate];
        else
            diffs = [diffs, f-s(2)];
            diffsIDX = [diffsIDX, iterate];
            % soundOnly=[soundOnly,iterate-1];
        end

    elseif length(s) > 1 && length(s) == length(f)
        for l=1:length(s)
           iterate=iterate+l;
           diffs = [diffs, f(l)-s(l)];
           diffsIDX = [diffsIDX, iterate];
        end
    else
        iterate=iterate+1;
        diffs = [diffs, f-s];
        diffsIDX = [diffsIDX, iterate];
    end
end
% diffs4sort=[diffs,1:length(diffs)];
[~,sortIDX] = sort(diffs);

m_s_all=rewardNewWOtouch.mpSound;
m_f_all=rewardNewWOtouch.mpFeeder;

m_s_all=cell2mat(m_s_all);
m_s_allSort=m_s_all((1:end-1),sortIDX);
m_f_allSort=m_f_all((1:end-1),sortIDX);

% diffs = [diffs;diffsIDX];
% [~,sortIDX] = sort(diffs(1,:));
% diffsSorted = diffs(:,sortIDX);
% id2sort = [diffsSorted(2,:),soundOnly];https://www.nature.com/articles/nature08010
% 
% m_s_allSort=m_s_all(:,id2sort);
% m_f_allSort=m_f_all(:,diffsSorted(2,:));

%% Plot mp vectors

time=linspace(-6000,5999,20*12000)';
time=repmat(time,1,size(m_s_allSort,2));
figure; plotShifted(time,m_s_allSort); xlabel('Time from Sound (ms)'); ylabel('Voltage (uV)');
time=linspace(-6000,5999,20*12000)';
time=repmat(time,1,size(m_f_allSort,2));
figure; plotShifted(time,m_f_allSort); xlabel('Time from Feeder (ms)'); ylabel('Voltage (uV)');
% plot(t_ms,signal); hold on; xline(t_ms(idxS),'r'); xline(t_ms(idxF),'g'); xlabel('Time (ms)'); ylabel('Voltage (uV)'); %Samples (20 kHz)

%% Get the signal vectors for mpRec matrix
ch=17; fs=20000; winPrev=6000; win=12000; %winPrev=2000; win=4000;
records=SA.recTable;
% animals=["PV153"]; %"PV143","PV106",
% curAnimal=animals;
% signal_s_allRes=[]; signal_f_allRes=[];
for r=1:size(rewardNew)
    r
    signal_s_allRes=[]; signal_f_allRes=[];
    curAnimal=rewardNew.animal(r);
    curRecName=rewardNew.record(r);
    row=sprintf('Animal=%s,recNames=%s', curAnimal, curRecName);
    current=SA.setCurrentRecording(row);
    microphoneTriggers=getStimMicrophoneTrig(SA,mean_template_sound_aligned,mean_template_feeder_aligned);
    s=microphoneTriggers.microphoneTriggers.Sound;
    for s_num=1:length(s)
        tStart=s(s_num)-winPrev;
        signalS = SA.currentDataObj.getData(ch, tStart, win);
        signalS = squeeze(signalS);
        % signalS = signalS - mean(signalS);
        signal_s_allRes=[signal_s_allRes,signalS];
        rewardNew.sigSound{r} = signal_s_allRes;
    end
    f=microphoneTriggers.microphoneTriggers.Feeder;
    for f_num=1:length(f)
        tStart=f(f_num)-winPrev;
        signalF = SA.currentDataObj.getData(ch, tStart, win);
        signalF = squeeze(signalF);
        % signalF = signalF - mean(signalF);
        signal_f_allRes=[signal_f_allRes,signalF];
        rewardNew.sigFeeder{r} = signal_f_allRes;
    end

end


%% Get the signal vectors for touch with sound
ch=17; fs=20000; winPrev=10000; win=20000; %winPrev=2000; win=4000;
records=SA.recTable;
row=sprintf('Animal=%s,recNames=%s', "PV106", "Hunter38");
current=SA.setCurrentRecording(row);
s=cell2mat(touch_sound.timeS(2));
f=cell2mat(touch_sound.timeF(2));
touchIdx=24;%[24;36;67];
% signal_s_allRes=[]; signal_f_allRes=[];
for s_num=1:length(s)
    tStart=s(s_num)-winPrev;
    signalS = SA.currentDataObj.getData(ch, tStart, win);
    signalS = squeeze(signalS);
    % signalS = signalS - mean(signalS);
    signal_s_allRes=[signal_s_allRes,signalS];
end

for f_num=1:length(f)
    tStart=f(f_num)-winPrev;
    signalF = SA.currentDataObj.getData(ch, tStart, win);
    signalF = squeeze(signalF);
    % signalF = signalF - mean(signalF);
    signal_f_allRes=[signal_f_allRes,signalF];
end


%% Plot signal vectors
winPrev=6000; win=20000;
% signal_s_allSortRes=allsigSSort;
% signal_s_allSortRes=allsigSSort;
fs_ms = 20;                 % samples per ms
Tms   = 20000;              % total duration in ms

time = (-winPrev : 1/fs_ms : Tms-winPrev - 1/fs_ms)';
% time=linspace(-winPrev,win-winPrev-1,2*win)';
% time=repmat(time,1,size(signal_s_allSortRes,2));
figure; plotShifted(time,signal_s_all_sound_Sort); xlabel('Time from Sound (ms)'); ylabel('Voltage (uV)');
hold on;
fTimes = round(sortedDiffsfinal, -2)';
% sortFeeder=
idx_timeFall=[];
responseFeeder = nan(size(signal_s_allSortRes));
for r=nNaNs+1:length(fTimes)
    idx_timeF = find(time==sortedDiffsfinal(r));
    idx_timeF
    idx_timeFall = [idx_timeFall, idx_timeF];
    responseFeeder(idx_timeF,r) = signal_s_allSortRes(idx_timeF,r);
end
% responseFeeder = signal_s_allSortRes(fTimes(4:end),4:size(signal_s_allSortRes,2));
% plotShifted(time, responseFeeder, 'Color','k*');%,'LineWidth',4);
plotShifted(time, responseFeeder, 'k*');%,'LineWidth',4);

% time=linspace(-winPrev,win-winPrev-1,20*win)';
% time=repmat(time,1,size(signal_f_allSortRes,2));
% figure; plotShifted(time,signal_f_allSortRes); xlabel('Time from Feeder (ms)'); ylabel('Voltage (uV)');

%% Spectrogram
% mean_signal_s_allRes=mean(signal_s_allRes,2);
% mean_signal_f_allRes=mean(signal_f_allRes,2);
% mean_signal_s_allRes_reshape = reshape(mean_signal_s_allRes, 1, 1, 14000);
% mean_signal_f_allRes_reshape = reshape(mean_signal_f_allRes, 1, 1, 14000);
nTrials_s=size(signalAll,2);
data_reshape = reshape(signalAll', 1, nTrials_s, 300000);

F=filterData(obj.currentDataObj.samplingFrequency(1));
F.downSamplingFactor=100;
F=F.designDownSample;
[filtSig,tFilt]=F.getFilteredData(data_reshape);

% figure; spectrogram(squeeze(filtSig),[],[],[],100,'yaxis');
% [S,F,T] = spectrogram(squeeze(signalS),[],[],[],20000,'yaxis');
% figure; spectrogram(squeeze(signalS),[],[],[],20000,'yaxis');
% xticks(T); xticklabels(time(:,1));
% 
% figure; pspectrum(signalS,20000,"spectrogram",TimeResolution=0.1,OverlapPercent=99,Leakage=0.85)

sig4spec=squeeze(filtSig);
% [sspec,fspec,tspec] = pspectrum(sig4spec(1,:),100,"spectrogram",TimeResolution=0.1,OverlapPercent=99,Leakage=0.85);
% figure; imagesc(tspec, fspec, 10*log10(sspec)); axis xy; colorbar;
S=[];
for r=1:size(sig4spec,1)
   % [Pspec,fspec,tspec] = pspectrum(sig4spec(r,:),200,"spectrogram",TimeResolution=0.5,OverlapPercent=90,Leakage=0.85);
   [Pspec,fspec,tspec] = pspectrum(sig4spec(r,:),200,"spectrogram",TimeResolution=0.1,OverlapPercent=90,Leakage=0.85);
    S(:,:,r) = Pspec;     % freq * time * trial
end


SdB = 10*log10(S);
% clim = [ ...
%     min(SdB(:),[],'omitnan'), ...
%     max(SdB(:),[],'omitnan') ];
clim = [ ...
    0, ...
    max(SdB(:),[],'omitnan') ];
figure;
tiledlayout(nTrials_s,1,'TileSpacing','compact','Padding','compact') %nTrials_s

for k = 1:nTrials_s
    nexttile
    imagesc(tspec,fspec, SdB(:,:,nTrials_s-k+1))
    axis xy
    caxis(clim)        %  lock color range
    % ylabel(sprintf('Trial %d',k))

    if k <nTrials_s
        set(gca,'XTickLabel',[])
    else
        eventTime = 6;   % seconds
        ax = gca;
        
        ax.XTickLabel = ax.XTick - eventTime;
        xlabel('Time from Sound (s)')

    end
end

colormap turbo
colorbar


% figure; imagesc(tspec, fspec, 10*log10(sspec)); axis xy; colorbar;
% ax = gca;
% numTicks=numel(ax.XTick);
% % xticks(t); 
% idxt = round( linspace(1, length(time(:,1)), numTicks) );
% t_sampled = time(idxt);% original time vector
% t_sampled = t_sampled*1e-3;
% xticklabels(t_sampled);  % your custom labels
% xlabel('Time from feeder (sec)')

%% Spectrogram try
[S,F,T] = spectrogram(squeeze(signalS));

S = abs(S);
figure; imagesc(time(:,1), F, S);
axis xy
xlabel('Custom Time')
ylabel('Frequency (Hz)')

%get the spectro
% [s] = spectrogram (signal);

[S,F,T] = spectrogram(squeeze(filtSig));
figure;
imagesc(abs(S)'); axis xy; colorbar;
imagesc(T, F, abs(S)'); axis xy; colorbar;
%% Plot a single signal
signal = squeeze(signal);
figure; 
subplot(2,1,1); plot(t_ms,signal); hold on; xline(t_ms(idxS),'r'); xline(t_ms(idxF),'g'); xlabel('Time (ms)'); ylabel('Voltage (uV)'); %Samples (20 kHz)
subplot(2,1,2); plot (t_ms,m); hold on; xline(t_ms(idxS),'r'); xline(t_ms(idxF),'g'); xlabel('Time (ms)'); ylabel('Voltage (uV)');

%% Beta2GammaRatio
% getDelta2BetaRatio(obj)
% [h]=plotDelta2BetaRatio(obj);
% event_oe_time=4.154870000000000e+04;



% for i=137:size(reward240526,1)
%     i
    for i=1:size(records,1) 
    i
    curRecName=string(records.recNames(i));
    if (contains(curRecName, "Hunter"))

    curAnimal=string(records.Animal(i));
    % curRecName=reward240526.Record(i);
    row=sprintf('Animal=%s,recNames=%s', curAnimal, curRecName); %curAnimal, curRecName
        % animals=["PV143","PV106","PV153"];
        animals=["PV126"];
   if ismember(curAnimal,animals)
        % SA=sleepAnalysis('/media/sil3/Data/Ofir/brainStatesWake.xlsx');
    % else
    %     SA=sleepAnalysis('/media/sil3/Data/Ofir/BrainStatesExperimentsOfir.xlsx');
    % end
    % records=SA.recTable;
     % if reward.TypeSound(i) == "sound_only"
     % if i == 1 || ...
     %   (curAnimal ~= records.Animal(i-1) || curRecName ~= records.recNames(i-1))
       
        current=SA.setCurrentRecording(row);
        getBeta2GammaRatio(SA);
        % getDelta2BetaRatioMed(SA);
        SA.getDelta2BetaRatio;
     end
       end
  end
    % end
% ch=SA.recTable.defaulLFPCh(SA.currentPRec);
% bgRatioFile=[SA.currentAnalysisFolder filesep 'bgRatio_ch' num2str(ch) '.mat'];
% SA.checkFileRecording(bgRatioFile,'beta to gamma file missing, please first run getBGRatio');
% load(bgRatioFile); %load data

%% Load b2g and create a mat for imagesc
% load('/home/ofir/Documents/MATLAB/sound.mat', 'mean_template_sound_aligned')
% load('/home/ofir/Documents/MATLAB/feeder.mat', 'mean_template_feeder_aligned')
% time_mat_all1 = {}; BG_mat_all1 = {};
% BG_mat_all_up80_106= [];
% BG_mat_all_delta= []; BG_mat_all_beta= []; BG_mat_all_betaHigh= []; BG_mat_all_gamma= [];
% BG_mat_all_208 = [];

% 228 rs: 92:96,98:99,103,105:108    ns:    97,100:102,104    167 rs:   75:77,80:82,86,88,90:91 - "ns"     ns:  78:79,83:85,87 - "rs" 

 BG_mat_Rsound_1673 = []; DB_mat_Rsound_1673 = []; F_mat_Rsound_1673 = []; %[110,113,115:122] 126:136
for r=[75:77,80:82,86,88,90:91]     %30:48 %8:16%17:29%30:48%17:29%1:7%8:16%17:29%30:size(screenDiode,1)%17:30 %108:159 %60:107 %31:59%1:16 1:size(reward,1)
     r
%     % if reward.TypeSound(r) == "sound_only"
    curAnimal=screenDiode.Animal(r);
    curRecName=screenDiode.Record(r);
    row=sprintf('Animal=%s,recNames=%s', curAnimal, curRecName);
    animals=["PV143","PV106","PV153"];
    if ismember(curAnimal,animals)
        SA=sleepAnalysis('/media/sil3/Data/Ofir/brainStatesWake.xlsx');
    else
        SA=sleepAnalysis('/media/sil3/Data/Ofir/BrainStatesExperimentsOfir.xlsx');
    end
    records=SA.recTable;
    current=SA.setCurrentRecording(row); 

ch=SA.recTable.defaulLFPCh(SA.currentPRec);

dbRatioFile=[SA.currentAnalysisFolder filesep 'dbRatio_ch' num2str(ch) '.mat'];
SA.checkFileRecording(dbRatioFile,'Delta to beta analysis missing, please first run getDelta2BetaRatio');
load(dbRatioFile); %load data
bgRatioFile=[SA.currentAnalysisFolder filesep 'bgRatio_ch' num2str(ch) '.mat'];
SA.checkFileRecording(bgRatioFile,'beta to gamma file missing, please first run getBGRatio');
load(bgRatioFile); %load data

% zeroTimes=1;

%%% FOR SCREEN DIODE
%%%%%[screenDiodeTriggers]=getScreenDiodeTrig(SA);
% stat = screenDiodeTriggers.screenDiodeTriggers.Stat;
%%%%%move = screenDiodeTriggers.screenDiodeTriggers.Move;
%%%%%zeroTimes = move -26000; %6 seconds before trial started
% zeroTimes = stat -26000; %6 seconds before trial started
%%%%%fullwin = 46000;
%matLen = length(find((t_ms>zeroTimes(1))&(t_ms<(zeroTimes(1)+fullwin))));
% % BG_mat = nan([length(zeroTimes),460]); %matLen

% microphoneTriggers=[SA.currentAnalysisFolder filesep 'microphoneTriggers_ch' num2str(3) '.mat'];
% load(microphoneTriggers); %load data
% % microphoneTriggers=load(SA.files.microphoneTriggers);
% s=microphoneTriggers.Sound;

template = template_NS_PV167_RL3;
reward = 1;
channels=SA.currentDataObj.analogChannelNames; %get the order of the analog channels
microphoneOrdinaryCh=find(cellfun(@(x) strcmp(x,'ADC2'),channels));
microphoneTriggers=[SA.currentAnalysisFolder filesep 'microphoneTriggers_ch' num2str(microphoneOrdinaryCh) '.mat'];
load(microphoneTriggers); %load data
s=microphoneTriggers.Sound;
fullwin = 24000; win = 12000;
zeroTimes = s -win;
BG_mat = []; DB_mat=[]; F_mc_up280=[];
% % BG_mc_up280 = nan(length(zeroTimes),328,460); %matLen
for i=1:length(zeroTimes)
    BG_ind = find((t_ms>zeroTimes(i))&(t_ms<(zeroTimes(i)+fullwin)));
    BG_mc = beta2gammaRatio(BG_ind);
        % BG_mat= [BG_mat;BG_mc'];
    BG_mat(i,1:length(BG_mc)) = BG_mc';


    DB_mc = bufferedDelta2BetaRatio(BG_ind);
    DB_mat(i,1:length(DB_mc)) = DB_mc';


    up280 = find(f>=0 & f<80);
    % F_mc =  (nan([length(zeroTimes),460])); %matLen
    F_mc =  (allFreqProfiles(up280,BG_ind)); 
    if ndims(F_mc) == 2
        F_mc = reshape(F_mc,1,size(F_mc,1),size(F_mc,2));
    end
    F_mc_up280(i,:,1:length(BG_ind)) = F_mc;


end
BG_mat_Rsound_1673  = [BG_mat_Rsound_1673 ;BG_mat];
DB_mat_Rsound_1673  = [DB_mat_Rsound_1673 ;DB_mat];
F_mat_Rsound_1673=cat(1,F_mat_Rsound_1673,F_mc_up280);
%%%

% for l=1:length(s)
    % sRnd=round(s, -2); % fRnd=round(f(l)/1000); %round to hundreds ms
    % startPointIdx = find(t_ms==sRnd)-6*5; %6 seconds before trial started
    % endPoint = fRnd*1000+6000;
    
    % endPoint = startPoint+20000-6000;
    % matLen = length(find((t_ms>=startPoint)&(t_ms<=endPoint)));
    % BG_ind_start = find(t_ms>=sRnd, 1)-6*10;
    % BG_mc = beta2gammaRatio(BG_ind_start:BG_ind_start+(6*10-1)); %window of 12 s
        % % BG_mc = beta2gammaRatio(BG_ind_start:BG_ind_start+(20*5-1)); %window of 20 s
        % timeVec = t_ms(BG_ind_start:BG_ind_start+199) - sRnd*1000;
        % BG_mat = zeros([length(zeroTimes),matLen]);
    % BG_mat = [];
    % for i=1:length(zeroTimes)
    % BG_ind_start = find((t_ms>=startPoint)); %&(t_ms<=endPoint));

        %     BG_mat= [BG_mat;BG_mc'];
    % BG_mat(i,1:length(BG_mc)) = BG_mc';
    
% end
    % time_mat_all1{end+1,1} = timeVec;
    % BG_mat_all1{end+1,1}= BG_mc;


%%% FOR SOUND
% s=reward.Sound(r);
% microphoneTriggers=load(SA.files.microphoneTriggers);
% s=microphoneTriggers.microphoneTriggers.Sound;
% fullwin = 24000;
% for i=1:length(s)
% BG_start = s-12*1000;
% matLen = length(find((t_ms>BG_start)&(t_ms<(BG_start+fullwin))));
% BG_ind = find((t_ms>BG_start)&(t_ms<(BG_start+fullwin)));
% BG_mc =  beta2gammaRatio(BG_ind);  %  bufferedDelta2BetaRatio
% BG_mat= [BG_mat;BG_mc'];
% BG_mat(i,1:length(BG_mc)) = BG_mc';
% BG_mat_all_sound_213=[BG_mat_all_sound_213;BG_mc'];
% end

% up280 = find(f>=0 & f<80);
% BG_mc =  (nan([length(zeroTimes),460])); %matLen
% BG_mc =  (allFreqProfiles(up280,BG_ind)); 
% if ndims(BG_mc) == 2
%     BG_mc = reshape(BG_mc,1,size(BG_mc,1),size(BG_mc,2));
% end
% BG_mc_up280(i,:,1:length(BG_ind)) = BG_mc;
% % cat(3,BG_mc_up280,BG_mc');
% end
% BG_mat_diode_228_up80=cat(1,BG_mat_diode_228_up80,BG_mc_up280);

% delta = find(f>=0 & f<5);
% beta = find(f>=10 & f<30);
% betaHigh = find(f>=30 & f<50);
% gamma = find(f>=60 & f<80);
% 
% BG_mc_delta =  allFreqProfiles(delta,BG_ind); 
% BG_mc_beta =  allFreqProfiles(beta,BG_ind); 
% BG_mc_betaHigh =  allFreqProfiles(betaHigh,BG_ind); 
% BG_mc_gamma =  allFreqProfiles(gamma,BG_ind); 
% BG_mat_all_delta=cat(3,BG_mat_all_delta,BG_mc_delta');
% BG_mat_all_beta=cat(3,BG_mat_all_beta,BG_mc_beta');
% BG_mat_all_betaHigh=cat(3,BG_mat_all_betaHigh,BG_mc_betaHigh');
% BG_mat_all_gamma=cat(3,BG_mat_all_gamma,BG_mc_gamma');

%%%

    
    a=1;
 end
% end
% BG_mat_diode_228_up80=permute(BG_mat_diode_228_up80,[3 2 1]);
% BG_mat_all(1:3,:)= [];
%% Plot b2g ratio
% % % BG_ind_downSampled = BG_ind(1:10:end); matLen_downSampled = matLen/10;
% % % time = (t_ms(BG_ind_downSampled)-s)/1000;
figure; 
% t = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
% --- Top: imagesc ---

% ax1 =subplot(2,3,1);
ax1 = nexttile;
% ax1 = subplot('Position', [0.1, 0.3, 0.7, 0.6]);
% tdiff = (par.movWin-par.movOLWin);
tdiff = 1000; %100
% timeAxis = (((0:matLen_downSampled-1) * tdiff-6000))/1000; %in s %(((0:matLen-1)-60) * tdiff)
% 
% % % fs1 = 10;                     % samples per second
% % % n = 240;  %120 sound old %460 screen diode      240 new sound             % number of samples
% % % time1 = (0:n-1)/fs1;              % time in seconds ? 0 to 11.9
% % % time1_centered = time1 - 12;     % now goes from -6 to +5.9    -12 new sound

% BG_mat_all2 = BG_mat_all;
% nan_mask = isnan(BG_mat_all2);
% BG_mat_all2(nan_mask) = min(BG_mat_all(:), [], 'omitnan') - 1;  % put below data range
% 
BG_mat_all_figutr = BG_mat_Rsound_1673;
% BG_mat_all_figutr([209],:) =[];
% ax1 = subplot('Position', [0.1, 0.3, 0.7, 0.6]); 
tdiff = (par.movWin-par.movOLWin);
timeAxis = (((0:size(BG_mat_all_figutr,2)-1) * tdiff)-win)/1000; %in s
imagesc(timeAxis ,1:size(BG_mat_all_figutr), BG_mat_all_figutr) %1:175 for PV213

% set(h, 'AlphaData', ~isnan(BG_mat_all))   % NaNs become transparent
% set(gca, 'Color', [1 1 1])       % background color = NaNs
% imagesc(time1_centered,1:length(zeroTimes), BG_mat)
% % imagesc(timeAxis ,[1:11 13:size(BG_mat_allSort)], BG_mat_allSort([1:11 13:size(BG_mat_allSort)],:))
axis xy % important: row 1 at top
title('PV167 Neutral Sound','FontSize',18)
% xlabel("Time from Sound Onset (s)");
ylabel('Trial #')
set(gca,'FontSize',18)
hold on
% xticks(-6:5)
% xticklabels(xticks)
% fTimes=(round(sortedDiffsfinal, -2)/1000);
% plot(fTimes(4:end),1:size(BG_mat_allSort),  'Color','g','LineWidth',2);
% % plot(fTimes([4:11 13:end]),[4:11 13:size(BG_mat_allSort)],  'Color','g','LineWidth',2);
% offset = time(7) - floor(time(7));   % ~0.59265
% % xticks_vals = (-6:5);   % or whatever range you want
% 
% xticks(timeAxis - offset)
% xticklabels(timeAxis)
% xline(0 - offset,'r--','LineWidth',2);

% xline(0,'r','LineWidth',2)
% xline(0,'g','LineWidth',2 , 'Label', 'Movement', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
% xline(-20,'r--','LineWidth',2 , 'Label', 'Static', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left')
xline(0,'g','LineWidth',2 , 'Label', 'Sound', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left')

  %add colorbar
c = colorbar(ax1);

% 
% caxis([min(BG_mat_all(:),[],'omitnan')-1, max(BG_mat_all(:))])
% c.Position = [0.905 0.585 0.02 0.35];
% % Manually adjust its position so it only covers the top tile
% ax1_pos = ax1.Position;    % [left bottom width height]
% c.Position = [ax1_pos(1)+ax1_pos(3)+0.01, ...  % right of ax1
%               ax1_pos(2), ...                  % same bottom as ax1
%               0.02, ...                        % width of colorbar
%               ax1_pos(4)];                      % same height as ax1

% % c.Position = [0.82 0.3 0.02 0.6];
% % Set vertical label for colorbar

% thr = prctile(BG_mat_all_figutr(:),90);
% clim([min(BG_mat_all_figutr(:)),thr]);

c = colorbar;
 c.Label.String = '\beta/\gamma'; % Set label text   \delta/\beta   \beta/\gamma
% % c.Label.Position = [4, 80, 0]; % Adjust label position (relative to colorbar)
% % c.Label.Rotation = 90; % Set label rotation to 0 degrees for vertical alignment


% ax2 = subplot(2,3,4);
ax2 = nexttile;
% ax2 = subplot('Position', [0.1, 0.1, 0.7, 0.15]);
 % yMean = nanmean(BG_mat);
yMean = nanmean(BG_mat_all_figutr);  %1:175 for PV213
% BG_mat_all0 = BG_mat_all; 
% BG_mat_all0(isnan(BG_mat_all)) = 0;
% yMean0 = mean(BG_mat_all0);
% % yMean = nanmean(BG_mat_allSort([1:11 13:size(BG_mat_allSort)],:));%, 1, 'omitnan');
plot(timeAxis,yMean,'k','Color',[0.5 0.5 0.5],'LineWidth',1.5); hold on; %1:175 for PV213
% xticks(timeAxis - offset)
% xticklabels(timeAxis)
% xticks(-6:5)
% xticklabels(xticks)
% xline(0,'g','LineWidth',2, 'Label', 'Movement', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left')
% xline(-20,'r--','LineWidth',2 , 'Label', 'Static', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
xline(0,'g','LineWidth',2 , 'Label', 'Sound', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left')


xlabel("Time from Sound (s)");  ylabel('\beta/\gamma '); %B2G  \delta/\beta     \beta/\gamma
% xlim([timeAxis(1) timeAxis(end)]);
set(gca,'FontSize',18)
linkaxes([ax1,ax2],'x');
% ax = findall(gcf,'type','axes');

%% Single trial analysis

trial=30; feeder=0;
figure; 
tiledlayout(8,1,'TileSpacing','compact','Padding','compact')

nexttile
timeSec=time/1000;
plot(timeSec,signal_s_all_sound_Sort(:,trial)); xlabel('Time from Sound (ms)'); ylabel('Voltage (uV)');
sigLim = [min(signal_s_all_sound_Sort(:),[],'omitnan'), max(signal_s_all_sound_Sort(:),[],'omitnan') ];
ylim(sigLim);
hold on;
xline(0,'r--','LineWidth',2)
% if feeder
%     xline(fTimes(trial),'g--','LineWidth',2)
% end

SdB=SdB_sound; tspec=tspec_s; fspec=fspec_s;

nexttile
eventTime = 6;   % seconds
tspecNew = tspec - eventTime;
imagesc(tspecNew,fspec, SdB(:,:,trial))
axis xy
caxis(clim)        %  lock color range
% ylabel(sprintf('Trial %d',k))

colormap turbo
c1=colorbar;
c1.Label.String = 'Power (dB)';

xlabel('Time from Sound (s)')
ylabel('Frequency (Hz)')
hold on;
xline(0,'w--','LineWidth',2)
% if feeder
%     xline(fTimes(trial),'w--','LineWidth',2)
% end

nexttile
fspecIdx= find(fspec>=60 & fspec<80);
% fspecCurs= fspec(fspecIdx(2:end));
plot(tspecNew, nanmean(SdB(fspecIdx,:,trial),1),'b')
ax = gca;
xlabel('Time from Sound (s)')
ylabel('Power (dB)')
text(ax, 1.02, 0.5, '60-80 Hz','Units','normalized','Rotation', -90,'HorizontalAlignment','center','VerticalAlignment','middle')
ylim(clim) 
hold on;
xline(0,'r--','LineWidth',2)
% if feeder
%     xline(fTimes(trial),'g--','LineWidth',2)
% end

nexttile
fspecIdx= find(fspec>=50 & fspec<60);
% fspecCurs= fspec(fspecIdx(2:end));
plot(tspecNew, nanmean(SdB(fspecIdx,:,trial),1),'c')
ax = gca;
xlabel('Time from Sound (s)')
ylabel('Power (dB)')
text(ax, 1.02, 0.5, '50-60 Hz','Units','normalized','Rotation', -90,'HorizontalAlignment','center','VerticalAlignment','middle')
ylim(clim) 
hold on;
xline(0,'r--','LineWidth',2)
% if feeder
%     xline(fTimes(trial),'g--','LineWidth',2)
% end

nexttile
fspecIdx= find(fspec>=30 & fspec<50);
% fspecCurs= fspec(fspecIdx(2:end));
plot(tspecNew, nanmean(SdB(fspecIdx,:,trial),1),'m')
ax = gca;
xlabel('Time from Sound (s)')
ylabel('Power (dB)')
text(ax, 1.02, 0.5, '30-50 Hz','Units','normalized','Rotation', -90,'HorizontalAlignment','center','VerticalAlignment','middle')
ylim(clim) 
hold on;
xline(0,'r--','LineWidth',2)
% if feeder
%     xline(fTimes(trial),'g--','LineWidth',2)
% end

nexttile
fspecIdx= find(fspec>=10 & fspec<=30);
% fspecCurs= fspec(fspecIdx(2:end));
plot(tspecNew, nanmean(SdB(fspecIdx,:,trial),1),'g')
ax = gca;
xlabel('Time from Sound (s)')
ylabel('Power (dB)')
text(ax, 1.02, 0.5, '10-30 Hz','Units','normalized','Rotation', -90,'HorizontalAlignment','center','VerticalAlignment','middle')
ylim(clim) 
hold on;
xline(0,'r--','LineWidth',2)
% if feeder
%     xline(fTimes(trial),'g--','LineWidth',2)
% end

nexttile
fspecIdx= find(fspec>=2 & fspec<=10);
% fspecCurs= fspec(fspecIdx(2:end));
plot(tspecNew, nanmean(SdB(fspecIdx,:,trial),1),'r')
ax = gca;
xlabel('Time from Sound (s)')
ylabel('Power (dB)')
text(ax, 1.02, 0.5, '2-10 Hz','Units','normalized','Rotation', -90,'HorizontalAlignment','center','VerticalAlignment','middle')
ylim(clim) 
hold on;
xline(0,'r--','LineWidth',2)
% if feeder
%     xline(fTimes(trial),'g--','LineWidth',2)
% end

nexttile
fspecIdx= find(fspec>=0 & fspec<=2);
% fspecCurs= fspec(fspecIdx(2:end));
plot(tspecNew, nanmean(SdB(fspecIdx,:,trial),1),'k')
ax = gca;
xlabel('Time from Sound (s)')
ylabel('Power (dB)')
text(ax, 1.02, 0.5, '0-2 Hz','Units','normalized','Rotation', -90,'HorizontalAlignment','center','VerticalAlignment','middle')
ylim(clim) 
hold on;
xline(0,'r--','LineWidth',2)
% if feeder
%     xline(fTimes(trial),'g--','LineWidth',2)
% end

%% Specific frequencies analysis

%spectrogram
signalNow=signal_s_all_sound_Sort;
nTrials_s=size(signalNow,2);
data_reshape = reshape(signalNow', 1, nTrials_s, 400000);
F=filterData(SA.currentDataObj.samplingFrequency(1));
F.downSamplingFactor=100;
F=F.designDownSample;
[filtSig,tFilt]=F.getFilteredData(data_reshape);
sig4spec=squeeze(filtSig);
S_sound=[];
for r=1:size(sig4spec,1)
   [Pspec_s,fspec_s,tspec_s] = pspectrum(sig4spec(r,:),200,"spectrogram",TimeResolution=0.1,OverlapPercent=90,Leakage=0.85);
    S_sound(:,:,r) = Pspec_s;     % freq * time * trial
end
SdB_touch = 10*log10(S_sound);
% clim = [ ...
%     0, ...
%     max(SdB(:),[],'omitnan') ]; 

%%
eventTime = 6;   % seconds
tspecNew = tspec_s - eventTime;
fspecIdx= find(fspec_s>=10 & fspec_s<=30);
betaTrace = nanmean(SdB_touch(fspecIdx,:,:),1);
betaTrace = permute(betaTrace, [2 3 1]);
tspecNew1 = tspecNew(496:696);
betaTrace1=betaTrace(496:696,:,:);
figure;
tiledlayout(2,1,'TileSpacing','compact','Padding','compact')

nexttile
imagesc(tspecNew1, 1:nTrials_s, betaTrace1')
axis xy
clim=([-20 60]);
 caxis(clim)        %  lock color range


colormap turbo
c1=colorbar;
c1.Label.String = 'Power (dB)';
title('10-30 Hz');
xlabel('Time from Sound (s)')
ylabel('Trial')
hold on;
xline(0,'k--','LineWidth',2)

nexttile
plot(tspecNew, betaTrace)
hold on; plot(tspecNew1, mean(betaTrace1,2),'k','LineWidth',2)
xlabel('Time from Sound (s)')
ylabel('Power (dB)')
ylim([-20 60]);
xline(0,'K--','LineWidth',2)



%%
figure; subplot(2,1,1);
% ax1 = subplot('Position', [0.1, 0.3, 0.7, 0.6]); 
imagesc(timeCommon, 1:nIter, Y)
axis xy                     % important: row 1 at top
%add colorbar
c = colorbar;
% c.Position = [0.82 0.3 0.02 0.6];
% % Set vertical label for colorbar
% c.Label.String = 'B/G ratio'; % Set label text

xline(0,'k--','LineWidth',1.5)

subplot(2,1,2);
yMean = mean(Y, 1, 'omitnan');
plot(timeCommon,yMean,'k','LineWidth',2); hold on;
xline(0,'k--','LineWidth',1.5)




%% add cols to the table

for r=1:size(rewardNewWOtouch,1)
    r
    curAnimal=rewardNewWOtouch.animal(r);
    curRecName=rewardNewWOtouch.record(r);
    row=sprintf('Animal=%s,recNames=%s', curAnimal, curRecName);
    current=SA.setCurrentRecording(row);
    microphoneTriggers=getStimMicrophoneTrig(SA,mean_template_sound_aligned,mean_template_feeder_aligned);
    rewardNewWOtouch.Sound{r}=[microphoneTriggers.microphoneTriggers.Sound];
    rewardNewWOtouch.Feeder{r}=[microphoneTriggers.microphoneTriggers.Feeder];
end

%% Sort the rows by feeder delay from sound

s=cell2mat(rewardNewWOtouch.Sound);
f=cell2mat(rewardNewWOtouch.Feeder);
diffs=f-s;
nNaNs = length(find(isnan(diffs)));
diffsDivided = diffs([find(~isnan(diffs)); find(isnan(diffs))]);
[sortedDiffs,sortIDX] = sort(diffsDivided);
sortedDiffsfinal=[sortedDiffs(end-nNaNs+1:end);sortedDiffs(1:end-nNaNs)];
sortIDXfinal=[sortIDX(end-nNaNs+1:end);sortIDX(1:end-nNaNs)];
% sortS=s(sortIDXfinal); sortF=f(sortIDXfinal);
allsigSSort=allsigS(:,sortIDXfinal);
allsigFSort=allsigF(:,sortIDXfinal(nNaNs+1:end));
% BG_mat_allSort=BG_mat_all(sortIDXfinal,:);
BG_mat_all_B2GlowSort=BG_mat_all_B2Glow(sortIDXfinal,:);

%% Get signals from reward table
% allsigS=cell2mat(rewardNewWOtouch.sigSound);
allsigS = [];
for i = 1:numel(rewardNewWOtouch.sigSound)
    signal=cell2mat(rewardNewWOtouch.sigSound(i));
    allsigS = [allsigS, signal];
end
allsigF = [];
for i = 1:numel(rewardNewWOtouch.sigFeeder)
    signal=cell2mat(rewardNewWOtouch.sigFeeder(i));
    allsigF = [allsigF, signal];
end

%% All Hunter sounds
[~, idx] = unique(reward(:, {'Animal','Record'}), 'rows', 'stable');
reward_unique = reward(idx,:);
for i = 80:size(reward_unique,1)
    curAnimal=reward_unique.Animal(i);
    curRecName=reward_unique.Record(i);
    row=sprintf('Animal=%s,recNames=%s', curAnimal, curRecName); %curAnimal, curRecName
    current=SA.setCurrentRecording(row);
    microphoneTriggers=getStimMicrophoneTrig(SA,mean_template_sound_aligned,mean_template_feeder_aligned);
    s=microphoneTriggers.microphoneTriggers.Sound;
    % f={microphoneTriggers.microphoneTriggers.Feeder};
    % if isempty(f)
    %     f= {NaN};
    % else
    %     for k = 1:numel(f)
    %         if isempty(f{k})
    %             f{k} = NaN;
    %         end
    %     end
    % end

    % find rows in table with this ID
    idxRows = find((reward.Animal == curAnimal) & (reward.Record == curRecName));

    % assign each element of s to the table rows
    reward.Sound(idxRows) = s(:);  
    % reward.Feeder(idxRows) = cell2mat(f(:));  
end

%%
winPrev=6000; win=20000;
signal_s_all_sound=[]; signal_f_all_sound=[]; 
signal_s_all_touch=[]; signal_f_all_touch=[]; 
diffs_sound=[]; diffs_touch=[];

for i=1:size(reward,1)
    i
    curAnimal=reward.Animal(i);
    curRecName=reward.Record(i);
    row=sprintf('Animal=%s,recNames=%s', curAnimal, curRecName);
    current=SA.setCurrentRecording(row);
    s=reward.Sound(i);
    f=reward.Feeder(i);

    % for s_num=1:length(s)
        tStart=s-winPrev;
        signalS = SA.currentDataObj.getData(17, tStart, win);
        signalS = squeeze(signalS);
        if reward.typeSound(i) == "sound_only"
            signal_s_all_sound=[signal_s_all_sound,signalS];
        else
            signal_s_all_touch=[signal_s_all_touch,signalS];
        end
    % end
    % for f_num=1:length(f)
    if ~isnan(f)
        tStart=f-winPrev;
        signalF = SA.currentDataObj.getData(17, tStart, win);
        signalF = squeeze(signalF);
        if reward.typeSound(i) == "sound_only"
             signal_f_all_sound=[signal_f_all_sound,signalF];
        else
             signal_f_all_touch=[signal_f_all_touch,signalF];
        end
    else
        if reward.typeSound(i) == "sound_only"
            signal_f_all_sound=[signal_f_all_sound,NaN(length(signalF),1)];
        else
            signal_f_all_touch=[signal_f_all_touch,NaN(length(signalF),1)];
        end
    end
    % end 

     if reward.typeSound(i) == "sound_only"
         diffs_sound=[diffs_sound,reward.Diffs(i)];
     else
         diffs_touch=[diffs_touch,reward.Diffs(i)];
     end
end



%%

% dibetaTrace1=betaTrace(496:696,:,:);ffs=reward.Feeder-reward.Sound;
% id_sound=[]; id_touch=[];
% for a=1:size(reward,1)
%     if reward.typeSound(a) == "sound_only"
%          id_sound=[id_sound,reward.Animal(a)];
%      else
%          id_touch=[id_touch,reward.Animal(a)];
%     end
% end


nNaNs_sound = length(find(isnan(diffs_sound)));
nNaNs_touch = length(find(isnan(diffs_touch)));
diffsDivided_sound = diffs_sound([find(~isnan(diffs_sound)), find(isnan(diffs_sound))]);
diffsDivided_touch = diffs_touch([find(~isnan(diffs_touch)), find(isnan(diffs_touch))]);
diffsDivided_sound = [diffsDivided_sound; [find(~isnan(diffs_sound)), find(isnan(diffs_sound))]];
diffsDivided_touch = [diffsDivided_touch; [find(~isnan(diffs_touch)), find(isnan(diffs_touch))]];
[sortedDiffs_sound,sortIDX_sound] = sort(diffsDivided_sound(1,:));
[sortedDiffs_touch,sortIDX_touch] = sort(diffsDivided_touch(1,:));
diffsDivided_sound_Sort = diffsDivided_sound(:,sortIDX_sound);
diffsDivided_touch_Sort = diffsDivided_touch(:,sortIDX_touch);
sortedDiffsfinal_sound=[sortedDiffs_sound(end-nNaNs_sound+1:end),sortedDiffs_sound(1:end-nNaNs_sound)];
sortedDiffsfinal_touch=[sortedDiffs_touch(end-nNaNs_touch+1:end),sortedDiffs_touch(1:end-nNaNs_touch)];
sortIDXfinal_sound=[diffsDivided_sound_Sort(2,end-nNaNs_sound+1:end),diffsDivided_sound_Sort(2,1:end-nNaNs_sound)];
sortIDXfinal_touch=[diffsDivided_touch_Sort(2,end-nNaNs_touch+1:end),diffsDivided_touch_Sort(2,1:end-nNaNs_touch)];
signal_s_all_sound_Sort=signal_s_all_sound(:,sortIDXfinal_sound);
signal_f_all_sound_Sort=signal_f_all_sound(:,sortIDXfinal_sound(nNaNs_sound+1:end));
signal_s_all_touch_Sort=signal_s_all_touch(:,sortIDXfinal_touch);
signal_f_all_touch_Sort=signal_f_all_touch(:,sortIDXfinal_touch(nNaNs_touch+1:end));

id_sound=id_sound'; id_touch=id_touch';
id_sound_Sort=id_sound(sortIDXfinal_sound);
id_touch_Sort=id_touch(sortIDXfinal_touch);

signal_s_all_sound_Sort1 = signal_s_all_sound_Sort((10000:14000),:);
figure; plotShifted(signal_s_all_sound_Sort1);
% figure; 
% subplot(2,2,1); plot (signal_s_all_sound_Sort); title('sound_only'); xlabel('Time from Sound (ms)'); ylabel('Voltage (uV)');
% subplot(2,2,2); plot (signal_f_all_sound_Sort); title('touch'); xlabel('Time from Feeder (ms)'); ylabel('Voltage (uV)');
% subplot(2,2,3); plot (signal_s_all_touch_Sort); title('sound_only'); xlabel('Time from Sound (ms)'); ylabel('Voltage (uV)');
% subplot(2,2,4); plot (signal_f_all_touch_Sort); title('touch'); xlabel('Time from Feeder (ms)'); ylabel('Voltage (uV)');

%% Wake trials

% reward = table([], [], [], [], 'VariableNames', {'Animal','Record','Idx','typeSound'});
records=SA.recTable;
animals=["PV143","PV106","PV153"];
idxTrial=ismember(records.Animal,animals);
records_filtered=records(idxTrial, :);
% winPrev=6000; win=20000;

for i=426:size(records_filtered,1) 
    i
    curRecName=string(records_filtered.recNames(i));
    if contains (curRecName,"Wake")
        curAnimal=string(records_filtered.Animal(i));
        row=sprintf('Animal=%s,recNames=%s', curAnimal, curRecName);  
        current=SA.setCurrentRecording(row);
        arena = SA.getArenaCSVs(1);
        rewardIdx = find(contains(arena.blockLog.Var5,"Reward given") & contains(arena.blockLog.Var5,"playing audio"));

        if ~isempty(rewardIdx)
            % analog_channels=SA.currentDataObj.analogChannelNames; %get the order of the analog channels
            % microphoneOrdinaryCh=find(cellfun(@(x) strcmp(x,'ADC2'),analog_channels)); %get the ordinary number of the microphone channel
            % % fs_analog = obj.currentDataObj.samplingFrequencyAnalog(microphoneOrdinaryCh);
            % 
            % LFPchannel=records_filtered.defaulLFPCh(i); %get the order of the analog channels
            % fs_LFP = SA.currentDataObj.samplingFrequency(LFPchannel);
    
            microphoneTriggers=getStimMicrophoneTrig(SA,mean_template_sound_aligned,mean_template_feeder_aligned,overwrite=0);
            s=microphoneTriggers.microphoneTriggers.Sound;
            f={microphoneTriggers.microphoneTriggers.Feeder};
            if isempty(f{:})
                f=nan(size(s));
            else
                for k = 1:numel(f)
                    if isempty(f{k})
                        f{k} = NaN;
                    end
                end
                f=cell2mat(f);
            end
        

            diff=f-s;
            % % find rows in table with this ID
            % idxRows = find((reward.Animal == curAnimal) & (reward.Record == curRecName));
            % 
            % % assign each element of s to the table rows
            % reward.Sound(idxRows) = s(:);  
            % reward.Feeder(idxRows) = cell2mat(f(:));  

            soundType = "sound_only";

            if rewardIdx > 1
                for id = 1:length(rewardIdx)
                    newRow = table(curAnimal, curRecName, rewardIdx(id), soundType, s(id), f(id) , diff(id),...
                       'VariableNames', reward.Properties.VariableNames);
                    reward = [reward; newRow];
                    % reward = [reward; {curAnimal, curRecName, rewardIdx, s, f, mp_sound, mp_feeder, signal_sound, signal_feeder, soundType, touchMat}];
                end
            else
                newRow = table(curAnimal, curRecName, rewardIdx, soundType, s, f , diff, ...
                       'VariableNames', reward.Properties.VariableNames);
                reward = [reward; newRow];
            end

        end
    end
end


%% Feeder check
m_f_all=[];
for i=1:size(reward,1) 
    i
    curAnimal=string(reward.Animal(i));
    curRecName=string(reward.Record(i));
    row=sprintf('Animal=%s,recNames=%s', curAnimal, curRecName);  
    current=SA.setCurrentRecording(row);
    if ~isnan(reward.Feeder(i))
        [m,~] = SA.currentDataObj.getAnalogData(5,reward.Feeder(i)-6000,12000);
        m = squeeze(m);
        m_f_all = [m_f_all,m];
    end
end


%%
x=d;

win = 4000;   % choose based on time scale (samples)
rmsSig = movstd(x, win);   % or movvar / movrms

thr = median(rmsSig);   % robust threshold
isOsc = rmsSig > thr;

diffs = diff(isOsc);
dropIdx = find(diffs == -1, 1, 'first');  % oscillations to quiet
riseIdx = find(diffs ==  -1, 1, 'last');   % quiet to oscillations

figure
plot(x), hold on
xline(dropIdx, 'r', 'LineWidth', 2)
xline(riseIdx, 'g', 'LineWidth', 2)
legend('signal','drop','rise')

d_MT=d(dropIdx:riseIdx);
d_MT_ms = d_ms(dropIdx:riseIdx);
%%
thresh = (mean(d_MT)-3*std(d_MT));
thresh_signal = d_MT < thresh;
% thresh_signal_changes = find(thresh_signal);
diffs = [0; diff(thresh_signal)];
dropIdx = find(diffs == -1);
% figure; plot(d_MT); hold on; xline(dropIdx,'r');

transitions = [0; diff(dropIdx)];
stratTrialsIdx = find(transitions > 18000 & transitions < 22000);
stratTrialsSamples = dropIdx(stratTrialsIdx-1);
screenDiodeTriggers = d_MT_ms(stratTrialsSamples);
figure; plot(d_MT_ms,d_MT); hold on; xline(screenDiodeTriggers);

%% TRY: screen diode

[screenDiodeTriggers]=getScreenDiodeTrig(SA);

%% Obtaining the neural data - screenDiode

winPrev=2000; win=44000;
% signal_143=[]; signal_106=[]; signal_153=[];
signal_208=[];
% signal_143_M=[]; signal_143_S=[];
% signal_106_M=[]; signal_106_S=[];
% signal_153_M=[]; signal_153_S=[];
% signal_208_M=[]; signal_208_S=[];


for i = 30:31 %size(screenDiode,1) %1:184

        curAnimal=screenDiode.Animal(i);
        curRecName=screenDiode.Record(i);
        row=sprintf('Animal=%s,recNames=%s', curAnimal, curRecName); %curAnimal, curRecName
        current=SA.setCurrentRecording(row);
        idxRow = find((string(records.Animal) == curAnimal) & (string(records.recNames) == curRecName));
        LFPchannel = records.defaulLFPCh(idxRow); %get the order of the analog channels

        [screenDiodeTriggers]=getScreenDiodeTrig(SA);
        % move = screenDiodeTriggers.screenDiodeTriggers.Move;
        stat = screenDiodeTriggers.screenDiodeTriggers.Stat;

        % for m=1:length(move)
        %     tStart=move(m)-winPrev;
        %     signalS = SA.currentDataObj.getData(LFPchannel, tStart, win);
        %     signalS = squeeze(signalS);
        %     switch screenDiode.Animal(i)
        %         case "PV143"
        %                 signal_143_M=[signal_143_M,signalS];
        %         case "PV106"
        %                 signal_106_M=[signal_106_M,signalS];
        %         case "PV153"
        %                 signal_153_M=[signal_153_M,signalS];
        %         case "PV208"
        %                 signal_208_M=[signal_208_M,signalS];
        %     end
        % end

        for s=1:length(stat)
            tStart=stat(s)-winPrev;
            signalS = SA.currentDataObj.getData(LFPchannel, tStart, win);
            signalS = squeeze(signalS);
            switch screenDiode.Animal(i)
                case "PV143"
                        signal_143=[signal_143,signalS];
                case "PV106"
                        signal_106=[signal_106,signalS];
                case "PV153"
                        signal_153=[signal_153,signalS];
                case "PV208"
                        signal_208=[signal_208,signalS];
            end
        end

end


%%
fs_ms = 20; Tms=42000;
time = (-winPrev : 1/fs_ms : Tms-winPrev - 1/fs_ms)';
trials = [1:50];
figure; plotShifted(time, signal_106_M(1:20000*5,:),'verticalShift',200); xlabel('Time from Meadow Stop (ms)'); ylabel('Voltage (uV)'); title('PV208, 3rd group (50 trials)');%(sprintf('PV153, n= %s trials', size(signal_153_M,2)));

% ylim([-0.2e4 5.7e04]);

%% Spectral analysis screenDiode

signal = [signal_s_153_H,signal_s_153_W];
data_reshape = reshape(signal', 1, size(signal,2), size(signal,1));

F=filterData(20000);
F.downSamplingFactor=100;
F=F.designDownSample;
[filtSig,tFilt]=F.getFilteredData(data_reshape);
sig4spec=squeeze(filtSig);

S153=[];%S213nall=[];
for r=1:size(signal,2)
%1.2e3+1-400:1.2e3+1+800-1
   [Pspec,fspec,tspec] = pspectrum(sig4spec(r,:),200,"spectrogram",TimeResolution=0.1,OverlapPercent=90,Leakage=0.5);
  % Pnorm = Pspec ./ sum(Pspec,1);
   S153(:,:,r) = Pspec;     % freq * time * trial
  % S213nall(:,:,r) = Pnorm;
end


% SdB213all = 10*log10(S213all); 

beta = find(fspec>=10 & fspec<30); gamma = find(fspec>=60 & fspec<80); 
b2g153 = [];
for i=1:size(S153,3)
    b2g153_trial = nanmean(S153(beta,:,i),1)./nanmean(S153(gamma,:,i),1);
    b2g153 = [b2g153;b2g153_trial];
end
% tdiff = 100;
% timeAxis = ([-60:200-60-1] * tdiff)/1000; %in s
% ax1 = subplot('Position', [0.1, 0.3, 0.7, 0.6]); 
% imagesc(timeAxis ,1:size(BG_mat_allSort), BG_mat_allSort)
 figure; imagesc(b2g153);
 clim([0 4e4]);
  h(2)=colorbar;
            set(h(2),'position',[0.9115    0.7040    0.0129    0.2220]);
            ylabel(h(2),'\beta/\gamma');

%%
% reward_sound_only = reward080326(reward080326.TypeSound == "sound_only", :);
% reward = reward_sound_only;

for i=62%145:152%72:size(reward,1) %1:184

    curAnimal=reward.Animal(i);
    curRecName=reward.Record(i);
    row=sprintf('Animal=%s,recNames=%s', curAnimal, curRecName); %curAnimal, curRecName
    current=SA.setCurrentRecording(row);
    s = reward.Sound(i);

    % run for first row automatically
    % OR if current row matches previous row
    if i == 62 || ...
       (curAnimal ~= reward.Animal(i-1) || curRecName ~= reward.Record(i-1))
        % close all;
        timeSeriesViewer(SA.currentDataObj) %for not SA object:rec
    end
    a=1;
end

%%
%% plot for all trial from 1 animal - get data:

SA = wakeAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWake.xlsx');
analysisFolder = '/media/sil1/Data/Nitzan/WakeStatesPaper/hunterAnalysis';
Animal = "PV153";
movWin = 1600;
movOLWin = 1200;
segWel = 1000;
fullwin = 50000;
bug_t = 6000;

curHunTrials = contains(SA.recTable.recNames,'Hunter') & contains(SA.recTable.Animal,Animal) & isnan(SA.recTable.Exclude);
hunterTable = SA.recTable(curHunTrials,:);
%%
% Pre-allocate struct array to collect data across recordings
% allData = struct('recName', {}, 'BG_mat', {}, 'endTrialT', {}, ...
                 % 'strikeTrialNum', {}, 'strikeT', {});
%
for i = 61:height(hunterTable)
    SA.setCurrentRecording(['Animal=' hunterTable.Animal{i} ',recNames=' hunterTable.recNames{i}]);
    ArenaCSVs = SA.getArenaCSVs('overwrite', 1);

    % Check for not-sync Hunters:
    frameDiff = length(ArenaCSVs.videoFrames) - length(ArenaCSVs.oeCamTrig);
    if frameDiff > 10 || frameDiff < -10
        continue
    end

    % Get BG:
    % BG = SA.getBeta2GammaRatio('maxVoltage', 1000, 'tStart', 0, 'win', 0, ...
    %     'movLongWin', 1000*60*30, 'movWin', movWin, 'movOLWin', movOLWin, ...
    %     'segmentWelch', segWel); %, 'overwrite', 1

    ch=SA.recTable.defaulLFPCh(SA.currentPRec);
    bgRatioFile=[SA.currentAnalysisFolder filesep 'bgRatio_ch' num2str(ch) '.mat'];
    BG = load(bgRatioFile); %load data


    zeroTimes = ArenaCSVs.startTrigSh - bug_t;
    endTrialT161 = ArenaCSVs.endTrigSh - zeroTimes;

    % Fix order: get strikeTrialNum before using it
    if isfield(ArenaCSVs, 'strikeTrialNum') && ~isempty(ArenaCSVs.strikeTrialNum)
        strikeTrialNum = ArenaCSVs.strikeTrialNum;
        strikeT = ArenaCSVs.strikeTrigSh - zeroTimes(strikeTrialNum);
    else
        strikeTrialNum = [];
        strikeT = [];
    end

    % Build BG_mat for this recording
    matLen = length(find((BG.t_ms > zeroTimes(1)) & (BG.t_ms < (zeroTimes(1) + fullwin))));
    BG_mat = zeros(length(zeroTimes), matLen);
    for j = 1:length(zeroTimes)
        BG_ind = find((BG.t_ms > zeroTimes(j)) & (BG.t_ms < (zeroTimes(j) + fullwin)));
        BG_mc = BG.beta2gammaRatio(BG_ind);
        BG_mat(j, 1:length(BG_mc)) = BG_mc';
    end

    % Append to allData struct array
    newEntry.recName       = hunterTable.recNames{i};
    newEntry.BG_mat        = BG_mat;
    newEntry.endTrialT     = endTrialT161;       % vector, one per trial
    newEntry.strikeTrialNum = strikeTrialNum; % indices of strike trials (local to this rec)
    newEntry.strikeT       = strikeT;         % times, same length as strikeTrialNum
    allData(end+1) = newEntry;
end

% Save for later plotting
% save(fullfile(analysisFolder, ['allData_' char(Animal) '.mat']), 'allData', 'movWin', 'movOLWin', 'segWel', 'bug_t', 'fullwin');
%% plot all trials
% Load and concatenate all recordings into one big matrix
% load(fullfile(analysisFolder, ['allData_' char(Animal) '.mat']));

tdiff = (movWin - movOLWin);
timeAxis = ((0:size(allData(1).BG_mat, 2)-1) * tdiff) / 1000;

% Concatenate BG_mat and adjust trial indices per recording
big_BG_mat    = [];
big_endTrialT = [];
big_strikeT   = [];
big_strikeTrial = [];
trialOffset   = 0;

for i = 1:length(allData)
    big_BG_mat    = [big_BG_mat;    allData(i).BG_mat];
    big_endTrialT = [big_endTrialT; allData(i).endTrialT(:)];

    % Offset strike trial indices to match position in big matrix
    if ~isempty(allData(i).strikeTrialNum)
        big_strikeTrial = [big_strikeTrial; allData(i).strikeTrialNum(:) + trialOffset];
        big_strikeT     = [big_strikeT;     allData(i).strikeT(:)];
    end

    trialOffset = trialOffset + size(allData(i).BG_mat, 1);
end

% Plot (same structure as before)

figure;


ax1 = subplot('Position', [0.1, 0.3, 0.7, 0.6]);
imagesc(timeAxis, 1:size(big_BG_mat,1), big_BG_mat);
hold on
xline(bug_t/1000, 'r', 'LineWidth', 1.5, 'Label', 'Bug Appearance', ...
      'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');

h1 = []; h2 = [];
for i = 1:size(big_BG_mat, 1)
    h1 = plot(big_endTrialT(i)/1000, i, '.', 'Color', 'k', 'MarkerSize', 5);
    if ismember(i, big_strikeTrial)
        strikeIdx = find(big_strikeTrial == i);
        h2 = plot(big_strikeT(strikeIdx)/1000, i, 'pentagram', 'Color', 'g', ...
                  'MarkerFaceColor', 'g', 'MarkerSize', 5);
    end
end
hold off

% Legend
legendHandles = h1(1);
legendLabels  = {'End Trial'};
if ~isempty(h2) && isgraphics(h2(1))
    legendHandles = [h2(1), h1(1)];
    legendLabels  = {'Strike', 'End Trial'};
end
legend(legendHandles, legendLabels, 'Position', [0.85 0.15 0.04 0.08], 'Box', 'off');
titlestr ='Mean BG for all trials - PV153'; 
title(titlestr)

%add colorbar
c = colorbar;
c.Position = [0.82 0.3 0.02 0.6];
% Set vertical label for colorbar
c.Label.String = 'G/B ratio'; % Set label text
c.Label.Position = [4, 80, 0]; % Adjust label position (relative to colorbar)
c.Label.Rotation = 90; % Set label rotation to 0 degrees for vertical alignment


ax2 = subplot('Position', [0.1, 0.1, 0.7, 0.15]);
plot(timeAxis, mean(big_BG_mat, 1,'omitnan'), 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('avg. B2G');
xline(bug_t/1000, 'r', 'LineWidth', 1.5);
linkaxes([ax1, ax2], 'x');

set(gcf,'PaperPosition',[.25 3 8 6])
% saveas(gcf,strcat([analysisFolder, filesep 'MeanBGallTrialsPV153.pdf']));

%%
  %% plotTrialsGB
 SA = wakeAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWake.xlsx');
% SA = wakeAnalysis('/media/sil3/Data/Ofir/BrainStatesExperimentsOfir.xlsx');
% analysisFolder = '/media/sil1/Data/Nitzan/WakeStatesPaper/hunterAnalysis';
Animal = "PV161";
% movWin = 1600;
% movOLWin = 1200;
% segWel = 1000;
% fullwin = 50000;
% bug_t = 6000;

curHunTrials = contains(SA.recTable.recNames,'Hunter') & contains(SA.recTable.Animal,Animal) & isnan(SA.recTable.Exclude);
hunterTable = SA.recTable(curHunTrials,:);
%%
% endTrialT161 = []; strikeTimeVec161=[];
% BG_mat_hunter_161 = []; DB_mat_hunter_161 = [];
% F_mat_hunter_161 = [];

for n = 65:height(hunterTable)

    SA.setCurrentRecording(['Animal=' hunterTable.Animal{n} ',recNames=' hunterTable.recNames{n}]);
    arenaCSVs = SA.getArenaCSVs('overwrite', 1);
    ch = SA.recTable.defaulLFPCh(SA.currentPRec);  
    win = 5000;
    fullwin = 60000;
    % bug_t = 6000;
    
    dbRatioFile=[SA.currentAnalysisFolder filesep 'dbRatio_ch' num2str(ch) '.mat'];
    SA.checkFileRecording(dbRatioFile,'delta to beta file missing, please first run getDBRatio');
    load(dbRatioFile); %load data

    bgRatioFile=[SA.currentAnalysisFolder filesep 'bgRatio_ch' num2str(ch) '.mat'];
    SA.checkFileRecording(bgRatioFile,'beta to gamma file missing, please first run getBGRatio');
    load(bgRatioFile); %load data
        
    arenaCSVFile=[SA.currentAnalysisFolder filesep 'ArenaCSV.mat'];
    SA.checkFileRecording(arenaCSVFile,'arenaCSV file missing, please first run getArenaCSV');
    load(arenaCSVFile); %load data
    
    
    zeroTimes = arenaCSVs.startTrigSh - win; %5 seconds before trial started

    strikeTimeVec = nan(length(zeroTimes),1);
    try
        strikeTrialNum =arenaCSVs.strikeTrialNum;
        strikeT = arenaCSVs.strikeTrigSh -arenaCSVs.startTrigSh(strikeTrialNum);
        strikeTimeVec(strikeTrialNum) = strikeT';
    catch
        strikeTimeVec = nan(length(zeroTimes),1);
    end
    strikeTimeVec161 = [strikeTimeVec161;strikeTimeVec];

    endTrialT = arenaCSVs.endTrigSh - arenaCSVs.startTrigSh;
    endTrialT161 = [endTrialT161;endTrialT];

    % matLen = length(find((t_ms>zeroTimes(1))&(t_ms<(zeroTimes(1)+fullwin))));
    % BG_mat = zeros([length(zeroTimes),matLen]);
    BG_mat=[]; DB_mat=[]; 
    F_mc_up280=[];
    for i=1:length(zeroTimes)
        BG_ind = find((t_ms>zeroTimes(i))&(t_ms<(zeroTimes(i)+fullwin)));
        BG_mc = beta2gammaRatio(BG_ind);
            % BG_mat= [BG_mat;BG_mc'];
        BG_mat(i,1:length(BG_mc)) = BG_mc';


        DB_mc = bufferedDelta2BetaRatio(BG_ind);
        DB_mat(i,1:length(DB_mc)) = DB_mc';
    
    
        up280 = find(f>=0 & f<80);
        % F_mc =  (nan([length(zeroTimes),460])); %matLen
        F_mc =  (allFreqProfiles(up280,BG_ind)); 
        if ndims(F_mc) == 2
            F_mc = reshape(F_mc,1,size(F_mc,1),size(F_mc,2));
        end
        F_mc_up280(i,:,1:length(BG_ind)) = F_mc;


    end
    BG_mat_hunter_161  = [BG_mat_hunter_161 ;BG_mat];
    DB_mat_hunter_161  = [DB_mat_hunter_161 ;DB_mat];
    F_mat_hunter_161=cat(1,F_mat_hunter_161,F_mc_up280);

a=1;
end
%%

idx_strikes = find(~isnan(strikeTimeVec161));
idx_rests = find(isnan(strikeTimeVec161));
% idx_2plot = idx_rests(endTrialT153(idx_rests)<10000); %([1:402,540:683],:)
%106:[919:end] LHN 143:[403:539,684:end] LHN [1:402,540:683] CF
idx_2plot = idx_strikes;
% idx_2plot = idx_strikes;
BG_mat_2plot=BG_mat_hunter_161(idx_2plot,:); endTrialT_2plot=endTrialT161(idx_2plot); strikeTimeVec_2plot=strikeTimeVec161(idx_2plot);
% BG_mat_2plot=BG_mat_hunter_161; endTrialT_2plot=endTrialT161; strikeTimeVec_2plot=strikeTimeVec161;
% rows = all(~isnan(BG_mat_2plot(:,51:70)), 2);
% BG_mat_2plot=BG_mat_2plot(rows,:);
% BG_mat_2plot=BG_mat_all106; endTrialT_2plot=endTrialT106; strikeTimeVec_2plot=strikeTimeVec106;
% BG_mat_2plot(82,:)=[]; endTrialT_2plot(82,:)=[]; strikeTimeVec_2plot(82,:)=[];
% BG_mat_2plot(403:404,:)=[];
% plot
figure;
%add subplot for the B/G ratio
ax1 = subplot('Position', [0.1, 0.3, 0.7, 0.6]);
tdiff = (par.movWin-par.movOLWin);
% tdiff = (parDBRatio.movWin-parDBRatio.movOLWin);
timeAxis = (((0:size(BG_mat_2plot,2)-1) * tdiff)-win)/1000; %in s
% timeAxis = (0:500:(length(BG_mat)-1)*500)/1000; % 0 ms to 19500 ms
imagesc(timeAxis,1:size(BG_mat_2plot),BG_mat_2plot)
axis xy % important: row 1 at top
xlabel("Time (s)")
hold on
title("PV228 No strikes, circles");
% Plot start time horizontal line:
xline(0,'r' ,'LineWidth', 1.5, 'Label', 'Bug Apperance', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
ylabel('Trial #')
% Plot end times and strike times
for i = 1:size(BG_mat_2plot, 1)
    h1 = plot(endTrialT_2plot(i)/1000, i, '.', 'Color', 'k','MarkerSize',20); % 'ko' for black circles
    % if ismember (i,strikeTrialNum)
        h2 = plot(strikeTimeVec_2plot(i)/1000, i, 'pentagram', 'Color', 'g','MarkerFaceColor','g','MarkerSize',8); % 'ko' for black circles
        % str_ind = str_ind+1;
    % end
end
%add legend and title:
% if ~isempty(strikeTrialNum)
    legend([h2(1), h1(1)], {'Strike','End Trial'},'Position',[0.85 0.15 0.04 0.08],'Box','off');%[x, y, width, height]
% else
    % legend(h1(1), 'End Trial','Position',[0.85 0.15 0.04 0.08],'Box','off');%[x, y, width, height]

% end

thr = prctile(BG_mat_2plot(:),95);
clim([min(BG_mat_2plot(:)),thr]);

%add colorbar
c = colorbar;
c.Position = [0.82 0.3 0.02 0.6];
% Set vertical label for colorbar
c.Label.String = 'B/G ratio'; % Set label text
c.Label.Position = [4, 80, 0]; % Adjust label position (relative to colorbar)
c.Label.Rotation = 90; % Set label rotation to 0 degrees for vertical alignment

% Plot Average
ax2 = subplot('Position', [0.1, 0.1, 0.7, 0.15]);
plot(timeAxis, nanmean(BG_mat_2plot,1),'LineWidth', 1.5)
xlabel('Time [s]'); ylabel('avg. B2G')
% ylim([10 85])
xline(0,'r' ,'LineWidth', 1.5)
linkaxes([ax1,ax2],'x');

% % Save plot:
% set(gcf,'PaperPosition',[.25 3 8 6])
% saveas(gcf,strcat(obj.currentPlotFolder, '/TrialsBG.pdf'));


%%
video = VideoReader('/media/sil1/Data/Ofir/regular_meadow_silent.mp4');

opticFlow = opticalFlowFarneback;

startTime = 20;   % seconds
endTime   = 40;   % seconds

video.CurrentTime = startTime;
% figure; 
k=1;
while hasFrame(video) && video.CurrentTime <= endTime
    frame = readFrame(video);

        % Convert to grayscale
    grayFrame = rgb2gray(frame);

    % Estimate flow
    flow = estimateFlow(opticFlow, grayFrame);
     % flow = estimateFlow(opticFlow,frame);
     motion(k) = sum(flow.Magnitude(:));

    k = k+1;
    % % Display image
    % imshow(frame)
    % hold on
    % 
    % % Display flow vectors
    % plot(flow,'DecimationFactor',[8 8],'ScaleFactor',10)
    % 
    % hold off
    % drawnow
end

%%
x = motion;   % your per-frame motion signal
fs = 30;             % frames per second
x = detrend(x);
z = hilbert(x);
env = abs(z);
env_smooth = movmean(env, round(fs*1)); % ~1 sec smoothing
Y = abs(fft(env_smooth));
f = (0:length(Y)-1)*(fs/length(Y));

plot(f, Y);
xlim([0 1.5]); % low frequency range


%%
t1 = linspace(0,20,600); t2 = linspace(0,20,200); t3 = linspace(0,40,2600); t4 = linspace(0,40,400);
figure;

% Left y-axis
% yyaxis left
% stimulus = arenaCSVs.bugs.x(1:2600); %env_smooth(1:end-1)

plot(t1, env_smooth_meadow/max(abs(env_smooth_meadow)), 'Color',[0.5 0.5 0.5], 'LineWidth', 1.5);
% plot(t3, stimulus/max(abs(stimulus)), 'Color',[0.5 0.5 0.5], 'LineWidth', 1.5)
hold on;
ylabel('Optic flow meadow')
% ylabel('Bug movement')

% b2g2plot153s=nanmean(BG_mat_all153(idx_strikes,201:600),1);
% b2g2plot153r=nanmean(BG_mat_all153(idx_rests,201:600),1);
f_band = gammaHigh;
f_g2plot106before = squeeze(nanmean(nanmean(F_mat_diode_106_up80(:,f_band,61:260),2),1)); f_g2plot106after = squeeze(nanmean(nanmean(F_mat_diode_106_up80(:,f_band,261:460),2),1)); 
f_g2plot143before = squeeze(nanmean(nanmean(F_mat_diode_143_up80(:,f_band,61:260),2),1)); f_g2plot143after = squeeze(nanmean(nanmean(F_mat_diode_143_up80(:,f_band,261:460),2),1)); 
f_g2plot153before = squeeze(nanmean(nanmean(F_mat_diode_153_up80(:,f_band,61:260),2),1)); f_g2plot153after = squeeze(nanmean(nanmean(F_mat_diode_153_up80(:,f_band,261:460),2),1)); 
f_g2plot167before = squeeze(nanmean(nanmean(F_mat_diode_167_up80(:,f_band,61:260),2),1)); f_g2plot167after = squeeze(nanmean(nanmean(F_mat_diode_167_up80(:,f_band,261:460),2),1)); 
f_g2plot167before_w = squeeze(nanmean(nanmean(F_mat_diode_167_up80_w(:,f_band,61:260),2),1)); f_g2plot167after_w = squeeze(nanmean(nanmean(F_mat_diode_167_up80_w(:,f_band,261:460),2),1)); 
f_g2plot208before = squeeze(nanmean(nanmean(F_mat_diode_208_up80(:,f_band,61:260),2),1)); f_g2plot208after = squeeze(nanmean(nanmean(F_mat_diode_208_up80(:,f_band,261:460),2),1)); 
f_g2plot228before = squeeze(nanmean(nanmean(F_mat_diode_228_up80(:,f_band,61:260),2),1)); f_g2plot228after = squeeze(nanmean(nanmean(F_mat_diode_228_up80(:,f_band,261:460),2),1)); 
f_g2plot228before_w = squeeze(nanmean(nanmean(F_mat_Rsound_1673(:,f_band,61:260),2),1)); f_g2plot228after_w = squeeze(nanmean(nanmean(F_mat_Rsound_1673(:,f_band,261:460),2),1)); 

f_g2plot = [f_g2plot106before,f_g2plot143before,f_g2plot153before,f_g2plot167before,f_g2plot167before_w,f_g2plot208before,f_g2plot228before,f_g2plot228before_w,f_g2plot106after,f_g2plot143after,f_g2plot153after,f_g2plot167after,f_g2plot167after_w,f_g2plot208after,f_g2plot228after,f_g2plot228after_w];

for c = 1:size(f_g2plot,2)
    f_g2plot(:,c) = f_g2plot(:,c)/max(abs(f_g2plot(:,c)));
end


d2b2plot106before = nanmean(DB_mat_diode_106(:,61:260),1); d2b2plot106after = nanmean(DB_mat_diode_106(:,261:460),1); 
d2b2plot143before = nanmean(DB_mat_diode_143(:,61:260),1); d2b2plot143after = nanmean(DB_mat_diode_143(:,261:460),1); 
d2b2plot153before = nanmean(DB_mat_diode_153(:,61:260),1); d2b2plot153after = nanmean(DB_mat_diode_153(:,261:460),1); 
d2b2plot167before = nanmean(DB_mat_diode_167(:,61:260),1); d2b2plot167after = nanmean(DB_mat_diode_167(:,261:460),1); 
d2b2plot167before_w = nanmean(DB_mat_diode_167_w(:,61:260),1); d2b2plot167after_w = nanmean(DB_mat_diode_167_w(:,261:460),1); 
d2b2plot208before = nanmean(DB_mat_diode_208(:,61:260),1); d2b2plot208after = nanmean(DB_mat_diode_208(:,261:460),1); 
d2b2plot228before = nanmean(DB_mat_diode_228(:,61:260),1); d2b2plot228after = nanmean(DB_mat_diode_228(:,261:460),1); 
d2b2plot228before_w = nanmean(DB_mat_Rsound_1673(:,61:260),1); d2b2plot228after_w = nanmean(DB_mat_Rsound_1673(:,261:460),1); 

b2g2plot106 = nanmean(BG_mat_hunter_161(idx_2plot,51:450),1); d2b2plot106 = nanmean(DB_mat_hunter_161(idx_2plot,51:450),1);
b2g2plot143 = nanmean(BG_mat_hunter_143(idx_2plot,51:450),1); d2b2plot143 = nanmean(DB_mat_hunter_143(idx_2plot,51:450),1);
b2g2plot153 = nanmean(BG_mat_hunter_153(idx_2plot,51:450),1); d2b2plot153 = nanmean(DB_mat_hunter_153(idx_2plot,51:450),1);
b2g2plot167 = nanmean(BG_mat_hunter_167(idx_2plot,51:450),1); d2b2plot167 = nanmean(DB_mat_hunter_167(idx_2plot,51:450),1);
b2g2plot208 = nanmean(BG_mat_hunter_208(idx_2plot,51:450),1); d2b2plot208= nanmean(DB_mat_hunter_208(idx_2plot,51:450),1);
b2g2plot213 = nanmean(BG_mat_Rsound_1673(:,51:450),1); d2b2plot213= nanmean(DB_mat_Rsound_1673(idx_2plot,51:450),1);
b2g2plot228 = nanmean(BG_mat_hunter_161(:,51:450),1); d2b2plot228= nanmean(DB_mat_hunter_161(idx_2plot,51:450),1);

f_d2plot213 = squeeze(nanmean(nanmean(F_mat_Rsound_1673(idx_2plot,delta,51:450),2),1));
f_b2plot213 = squeeze(nanmean(nanmean(F_mat_Rsound_1673(idx_2plot,beta,51:450),2),1));
f_g2plot213 = squeeze(nanmean(nanmean(F_mat_Rsound_1673(idx_2plot,gammaHigh,51:450),2),1));

f_g2plot = [f_g2plot106,f_g2plot143,f_g2plot153,f_g2plot167,f_g2plot208,f_g2plot213,f_g2plot228];
for c = 1:size(f_g2plot,2)
    f_g2plot(:,c) = f_g2plot(:,c)/max(abs(f_g2plot(:,c)));
end


% % Right y-axis
% % yyaxis right
% plot(t4,  b2g2plot153r/max(abs(b2g2plot153r)), 'LineWidth', 1.5); hold on;
% plot(t4,  b2g2plot153s/max(abs(b2g2plot153s)), 'LineWidth', 1.5); hold on;
d2b2plot = [d2b2plot106before;d2b2plot143before;d2b2plot153before;d2b2plot167before;d2b2plot167before_w;d2b2plot208before;d2b2plot228before;d2b2plot228before_w;d2b2plot106after;d2b2plot143after;d2b2plot153after;d2b2plot167after;d2b2plot167after_w;d2b2plot208after;d2b2plot228after;d2b2plot228after_w]';
for c = 1:size(d2b2plot,2)
    d2b2plot(:,c) = d2b2plot(:,c)/max(abs(d2b2plot(:,c)));
end
plot(t2,  b2g2plot106/max(abs(b2g2plot106)), 'LineWidth', 1.5); hold on;
% plot(t2, b2g2plot208/max(abs(b2g2plot208)), 'LineWidth', 1.5); 
% plot(t2, b2g2plot106/max(abs(b2g2plot106)), 'LineWidth', 1.5); 
ylabel('Norm')

xlabel('Time (s)')

[pks, locs] = findpeaks(env_smooth_meadow,'MinPeakDistance',90); %stimulus
t = locs * 20 / 600;           % if using a direct proportional conversion   % * 40 / 2600
xlim([0 20])     % [0 40]
xticks(t)
grid on
% legend('Bug movement','B/G 228 rest','B/G 228 strike')
% legend('Optic flow meadow','B/G 208')
% legend('B/G 153','B/G 208','B/G 106')
% title('Comparison of two signals')

%%
 figure;
% plotShifted(t2,f_g2plot(:,1:8), 'Color',[0.6 0 0], 'LineWidth', 1.5);
% plotShifted(t2,f_g2plot(:,9:16), 'Color',[0 0.5 0], 'LineWidth', 1.5);
plotShifted(t4, f_b2plot, 'Color',[0 0.5 0], 'LineWidth', 1.5); % 'verticalShift', 1,
xlabel('Time from Bug Appearance (s)'); ylabel('Norm \beta');  %  \delta/\beta    \beta/\gamma
% xlabel('Time from Meadow Movement (s)'); ylabel('Norm \gamma');  %\beta/\gamma
hold on; 
plot(t3, stimulus/max(abs(stimulus)), 'Color',[0.5 0.5 0.5], 'LineWidth', 1.5); %env_smooth_meadow/max(abs(env_smooth_meadow))
% plot(t1, env_smooth_meadow/max(abs(env_smooth_meadow)), 'Color',[0.5 0.5 0.5], 'LineWidth', 1.5);
% legend('106 Stat','143 Stat','153 Stat','167 Stat','167 White Stat','208 Stat','228 Stat','228 White Stat','106 Move','143 Move','153 Move','167 Move','167 White Move','208 Move','228 Move','228 White Move','Optic flow meadow');
[pks, locs] = findpeaks(stimulus,'MinPeakDistance',90); %env_smooth_meadow
t = locs * 40 / 2600;           % if using a direct proportional conversion   %  stimulus:* 40 / 2600 meadow:   * 20 / 600
xlim([0 40])     % stimulus:[0 40] meadow: [0 20]
xticks(t)
grid on
title('Beta of Bug Appearance');
% ylim([0.7 2.7])

%% D2B sleep analysis
vec = bufferedDelta2BetaRatio;
threshold = 600;
% idx = find(vec(1:end-1) <= threshold & vec(2:end) > threshold) + 1;
nSamples = 50;

% [pks, locs] = findpeaks(vec,'MinPeakDistance',50,'MinPeakProminence',50); %,'MinPeakDistance',90
[valls, locs] = findpeaks(-vec,'MinPeakDistance',50,'MinPeakProminence',50); %,'MinPeakDistance',90
valls = -valls;
% Extract segments (each row is one event)
% segments = zeros(numel(idx), nSamples);
segments=[];
for k = 10:numel(locs)
    segments(k,:) = vec(locs(k)-nSamples:locs(k)+nSamples-1);
end


tdiff = (parDBRatio.movWin-parDBRatio.movOLWin);
timeAxis = (((0:size(segments,2)-1) * tdiff)-50000)/1000; %in s

% Plot
figure;
imagesc(timeAxis,1:length([50:150]),segments(50:150,:)); %40:470
xlabel('Time from D2B peak (s)');
ylabel('Cycle');
clim([0 800]);
colorbar;
