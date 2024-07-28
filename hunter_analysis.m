%% wake state - bug appearnce
% here I want to shwo traces of the bug apperance in PV149 during circle
% trials/ 
%
% TO DO:
%   -- change the variables of the paths so i can change it easily / it
%       will read from the excel.
%   -- Take the trigger channels from the excel - atimulation channel and
%       camera triggers chanel.
%   -- triggers are by defenition will be different for this animal. adjust
%   the code accordingly. problem expained near "num_missing_frames".
%   
% WHY ARE THE TIMESTAMPS IN S and not MS?!
%
% trigger chanels for PV149: cams - 7-8, stimulations - 15-16/

%% getting the data for this sepecific experimeent:
% notice - now everything should be changes manualy.
SA = sleepAnalysis('/media/sil1/Data/Pogona Vitticeps/brainStatesWake.xlsx');
SA.setCurrentRecording(['Animal=PV157,recNames=Hunter5']);
% get the arena timings, and translte to oe triger times:
ArenaCSVs = SA.getArenaCSVs(1);

analysisFolder = '/media/sil3/Data/Pogona_Vitticeps/NitzanAnalysisFiles/hunterAnalysis';
%% changing the triggers according to the new time shift:
% for PV159, Hunter5:
% startFrameSh = ArenaCSVs.startFrame - ArenaCSVs.IRFrames(1);
% EndFramSh = ArenaCSVs.endFrame - ArenaCSVs.IRFrames(1)+2;
% strikeFrameSh = ArenaCSVs.strikesFrame - ArenaCSVs.IRFrames(1)+1;
% 
% % for PV161, Hunter9: and PV162,Hunter27:
% startFrameSh = ArenaCSVs.startFrame - ArenaCSVs.IRFrames(1);
% EndFramSh = ArenaCSVs.endFrame - ArenaCSVs.IRFrames(1);
% strikeFrameSh = ArenaCSVs.strikesFrame - ArenaCSVs.IRFrames(1);
% 
% startTrigSh = ArenaCSVs.oeCamTrigs(startFrameSh) ;
% endTrigSh = ArenaCSVs.oeCamTrigs(EndFramSh) ;
% strikeTrigSh = ArenaCSVs.oeCamTrigs(strikeFrameSh ) ;



%% getting the triggers from the OERecording
% correction of frames for the behavioral data:
% for this recording - PV149, Hunter61, the frames are shifted 4 triggers
% (the videos is from trigger number 4 until the end of num of frames).
% the CSVs are 15 frames from the times of the triggers (11 after
% correction)
% shifts = SA.recTable.blockShift(SA.currentPRec);
% shiftsnum = str2num(shifts{1});
% 
% % timestamps for the events:
% if isempty(shiftsnum) ==1
%     oeStartTrial = OEcamTrig(sTrialFrame);
%     oeEndTrial = OEcamTrig(eTrialFrame);
%     oeStrikes = OEcamTrig(strikesFrames);
% elseif numel(shiftsnum) ==1
%     oeStartTrial = OEcamTrig(sTrialFrame-shiftsnum);
%     oeEndTrial = OEcamTrig(eTrialFrame-shiftsnum);
%     oeStrikes = OEcamTrig(strikesFrames-shiftsnum);
% else 
%     oeStartTrial = OEcamTrig(sTrialFrame-shiftsnum(1));
%     oeEndTrial = OEcamTrig(eTrialFrame-shiftsnum(1));
%     oeStrikes = OEcamTrig(strikesFrames-shiftsnum(2));
% end

%% check synchronization:
videoFile=SA.recTable.VideoFiles{SA.currentPRec};
%check if its a top video and change to back video:
videoFileparts = strsplit(videoFile,'/');
vidName = videoFileparts{end};
if contains(vidName,'top')
    %replace top with back:
    vidName = strrep(vidName,'top','left');
    videoFileparts{end} = vidName; %replace them
    videoFile = strjoin(videoFileparts,'/');
end


VR=VideoReader(videoFile); %creates class for the video file
b=1;
    for u=1:length(startFrameSh) %trigger_frame = frame of e.g. bug appearance - a few frames to check video before bug should appear
        for i=0:10% or
            frame1=read(VR, startFrameSh(u)+i+5);
            frame1 = insertText(frame1,[1010,0],string(i), 'FontSize',70, 'TextColor','white', 'BoxColor', 'black', 'BoxOpacity', 1); %add text to image with framenumber
            frame1 = insertText(frame1,[1210,0],string(u), 'FontSize',70, 'TextColor','white', 'BoxColor', 'black', 'BoxOpacity', 1); %add text to image with framenumber
            imshow(frame1)
            pause() %any button will continue to next frame
            if get(gcf,'CurrentCharacter')=='q' %if you want to endends the loop with q
                break
            end 
            if get(gcf,'CurrentCharacter')=='b' %if you want to endends the whole thing
                b=5;
                break
            end 
        end
        if b==5
            break
        end
    end


%% plot stikes with imagesc:
%getting the data from 4 seconds around the event:
%samples per sec
defCh = SA.recTable.defaulLFPCh(SA.currentPRec);

strikeData = squeeze(SA.currentDataObj.getData(defCh,strikeTrigSh-2000,4000));
% plot
figure; plot(strikeData');
% imagesc:
figure; imagesc(strikeData, [-400, 500]);colorbar;xline(2*20000, 'LineWidth',1.5,'Color', 'r')
figure;imagesc(strikeData(:,[1:(1.75*20000), 2.1*20000:end]),[-400,500]); colorbar;
xline(1.75*20000, 'LineWidth',1.5,'Color', 'r')
%% plot start trial with the imagesc + Plot Shift:
%getting the data from 4 seconds around the event:
dateWin = 10000; % in ms
% defCh = 3;
% fs = SA.currentDataObj.samplingFrequency(1); %samples per sec
startTrigSh = ArenaCSVs.startTrigSh;
[startData, startDataTimes]= SA.currentDataObj.getData(defCh,startTrigSh-(dateWin/2),dateWin);
startData = squeeze(startData);
% imagesc:
% figure; imagesc(startData, [-400, 500]);
% colorbar;xline(dateWin/2*20, 'LineWidth',1.5,'Color', 'r')
% figure;imagesc(startData(:,[1:(1.75*20000), 2.1*20000:end]),[-400,500]); colorbar;
% xline(1.75*20000, 'LineWidth',1.5,'Color', 'r')


% plot shifted - start trial :
% good trials for PV157,Hunter5:
goodTrials = [1,3:4, 6:13, 15:20];
% good trials for PV161,Hunter9:
% goodTrials = [1,3:15,18:20];
% good trials for PV162,Hunter27:
% goodTrials = [2:5,8:10,12,14:20];
% 12,18 - didn't want to go and i think the signal is different. 
figure;
set(gcf,'PaperPosition',[.25 3 8 6])
plotShifted(startDataTimes, startData(goodTrials,:).','verticalShift',-800); % exuluding trials: 2 (didn't look) , 5 (didn't look at first second), 14 (going to get the real bugs)
xline(dateWin/2,'r'); ylabel('uV');xlabel ('Time[ms]');
title('Bug apperance, for one block')
saveas(gcf,strcat(SA.currentPlotFolder, '/bugApperanceFull.pdf'));


%% get behavioral data
% % read the parquet 
videoPath = SA.recTable.VideoFiles(SA.currentPRec);
[videosFolderPath,vidName,~] = fileparts(videoPath{1});
locFilename = strcat(videosFolderPath,'/predictions/front_head_ephys_resnet_101__',vidName,'.parquet');

% check file exists:
if ~isfile(locFilename)
    error('File does not exist: %s', locFilename);
end
locTable = parquetread(locFilename);
% Display the first few rows and summary of the data
disp(locTable(1:5, :));
% summary(locTable);
%% plot example trial - nose y-time, LFPs
% set the start of example timefrom rec + zeroFrame
zeroTime = 82812.4; 
inds = find(ArenaCSVs.oeCamTrigs<zeroTime); 
zeroFrame = inds(end)+1;

% get the recdata 
fullwin = 14000;
[recData_M,recData_T] = SA.currentDataObj.getData(defCh,zeroTime,fullwin);
recData_M = squeeze(recData_M);

% get the behavioral data:
frameWin = round(fullwin/(1000/ArenaCSVs.videoFPS));
nose_y = locTable.x__nose___y__(zeroFrame:zeroFrame+frameWin-1); 
nose_x = locTable.x__nose___y__(zeroFrame:zeroFrame+frameWin-1);
nose_t = linspace(0,fullwin/1000,frameWin);

%% plot
figure;
set(gcf,'Position',[100 20 800 600])
ax1 = subplot(3,3,[1:3]);
%plot the nose location
plot(nose_t,nose_y);
ylabel('Nose from screen [cm]');
ylim([-20 90])
% plot the LFPs
ax2 = subplot(3,3,[4:6]);
plot(recData_T/1000,recData_M, 'Color','black');
ylabel('uV');
xlabel('Time[s]');

linkaxes([ax1,ax2],'x');

% add the timings of the strik and end and start trial. 
startTime = (ArenaCSVs.startFrameSh(1)-zeroFrame)/(1000/30);
strikeFrame = (ArenaCSVs.strikeFrameSh(1)+2);
strikeTime = (ArenaCSVs.oeCamTrigs(strikeFrame)-zeroTime)/1000;
endtime = (ArenaCSVs.endTrigSh(1)-zeroTime)/1000;
xline(ax1,startTime,'b--' ,'LineWidth', 2, 'Label', 'Bug Apperance', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'right');
xline(ax2,startTime,'b--','LineWidth', 2);
xline(ax1,strikeTime,'r--' ,'LineWidth', 2, 'Label', 'Strike', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'right');
xline(ax2,strikeTime,'r--','LineWidth', 2);
xline(ax1,endtime,'g--' ,'LineWidth', 2, 'Label', 'Bug Disapearance', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'right');
xline(ax2,endtime,'g--','LineWidth', 2);
hold on
fs = SA.currentDataObj.samplingFrequency(defCh);
cutoff1 = 2.5*fs:3.5*fs;
yLimits = ylim;
patch([cutoff1(1)/fs cutoff1(end)/fs  cutoff1(end)/fs cutoff1(1)/fs], ...
      [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
      'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Highlight with a semi-transparent red
cutoff2 = 5*fs:6*fs;
patch([cutoff2(1)/fs cutoff2(end)/fs  cutoff2(end)/fs cutoff2(1)/fs], ...
      [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
      'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); 

cutoff3 = 12*fs:13*fs;
patch([cutoff3(1)/fs cutoff3(end)/fs  cutoff3(end)/fs cutoff3(1)/fs], ...
      [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
      'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Highlight with a semi-transparent red
hold off

% add close up graphs:
ax3 = subplot(3,3,7);
plot(recData_T(cutoff1)/1000,recData_M(cutoff1), 'Color','black');
xlabel('Time[s]');
ylabel('uV');

ax4 = subplot(3,3,8);
plot(recData_T(cutoff2)/1000,recData_M(cutoff2), 'Color','black');
xlabel('Time[s]');

ax5 = subplot(3,3,9);
plot(recData_T(cutoff3)/1000,recData_M(cutoff3), 'Color','black');

xlabel('Time[s]');
linkaxes([ax3,ax4,ax5],'y');


set(gcf,'PaperPosition',[.25 3 8 6])
saveas(gcf,strcat(SA.currentPlotFolder, '/TrialNoiseHL.pdf'));
% saveas(gcf,strcat(SA.currentPlotFolder, '/TrialNoiseHL.jpg'));
%% get the dendogram of the signal - only the best parts
% using Mark,s function 
tStarts = startTrigSh(goodTrials)-5000;
% tStarts = startTrigSh-5000;
% Freq = getFreqBandDetectionN(SA,'multiStart',1,'tStart',tStarts,'fMax',80);
getFreqBandDetectionN(SA,'overwrite',1,'multiStart',1,'tStart',tStarts,'fMax',100,'win',10000, 'WelchOL',0.5,...
   'binDuration',2500,'segmentLength',2000);
SA.plotFreqBandDetection('savePlots',false);

%             addParameter(parseObj,'segmentLength',1000);
%             addParameter(parseObj,'WelchOL',0.5);
%             addParameter(parseObj,'binDuration',10000)

%% re-write the pwelch and dendrogram:

% get the data:
win = 5000; % in ms
fs = SA.currentDataObj.samplingFrequency(1); %samples per sec
% Bug apperance: before and after
% defCh = 3;
postStart = SA.currentDataObj.getData(defCh,startTrigSh(goodTrials)+400,win);
preStart = SA.currentDataObj.getData(defCh,startTrigSh(goodTrials)-win, win);
prePost = [squeeze(preStart);squeeze(postStart)];

% prePost= transpose(reshape(startData(goodTrials,:),[win*fs/2000,length(goodTrials)*2]));


welchWin = 1*fs; 
welchOL = 0.5; % the percent of overlap
samplesOL = welchOL*welchWin;
%calculate the ppx:
[pxx,f] = pwelch(prePost',welchWin,samplesOL,[],fs);

fMax = 110;
% times=(tStart+binDuration/2):binDuration:(tStart+win);
% dendrogram:          
maxDendroClusters = 2;

p=find(f<fMax);
pp=find(sum(pxx(p,:))<0.4e6); %reject signals with very high amplitudes (probably noise)
            
sPxx=pxx(p,pp);
freqHz=f(p);
normsPxx=bsxfun(@rdivide,sPxx,mean(sPxx,2));
corrMat=corrcoef(normsPxx);
% times=times(pp);


if maxDendroClusters==2

    [DC,order,clusters]=DendrogramMatrix(corrMat,'linkMetric','euclidean','linkMethod','ward','maxClusters',maxDendroClusters);

    S1=mean(normsPxx(:,clusters==1),2);
    S2=mean(normsPxx(:,clusters==2),2);
    if mean(S1(1:3))>mean(S2(1:3))
        crossFreq=freqHz(find(S2-S1>=0,1,'first'));
    else
        crossFreq=freqHz(find(S1-S2>=0,1,'first'));
    end
else
    [DC,order,clusters]=DendrogramMatrix(corrMat,'linkMetric','euclidean','linkMethod','ward','maxClusters',maxDendroClusters);

    for i=1:maxDendroClusters
        S(:,i)=mean(normsPxx(:,clusters==i),2);
    end
    crossFreq=[];
end

%not normelized PSD:
sPxx1 = mean(sPxx(:,clusters==1),2);
sPxx2 = mean(sPxx(:,clusters==2),2);



% plot

k = length(goodTrials); % Threshold

% Display the correlation matrix
fig1 = figure;
set(fig1, 'Position', [100, 100, 900, 800]); % [x, y, width, height]

imagesc(DC);
colorbar;
% axis equal;
hold on;

% Set axis labels
% xticks(1:length(order));
% yticks(1:length(order));

% Loop through each segment number and display it
for i = 1:length(order)
    if order(i) <= k
        text_color = 'black';
        background_color = [0.7, 0.7, 1]; % Light purple background
    else
        text_color = 'black';
        background_color = [1, 0.85, 0.7]; % Light orange background
    end
    
    % Display the segment number with background
    text(-0.5, i, num2str(order(i)), 'HorizontalAlignment', 'center', ...
        'BackgroundColor', background_color, 'Color',text_color , ...
        'FontWeight', 'bold');
    text(i, -0.5, num2str(order(i)), 'HorizontalAlignment', 'center', ...
        'BackgroundColor', background_color, 'Color',text_color , ...
        'FontWeight', 'bold');

end
% Adjust limits
% xlim([-1, length(order)]);
% ylim([0.5, length(order) + 0.5]);
set(gca, 'XTick', [], 'YTick', []);
hold off;
set(fig1,'PaperPosition',[.25 3 8 6])
saveas(gcf,strcat(SA.currentPlotFolder, '/corelogram.pdf'));

%plot the freqs:
figure;
set(gcf, 'Position', [300, 300, 700, 500]); % [x, y, width, height]
plot(freqHz,S1,'Color',[0.7, 0.7, 1], 'LineWidth',2);
hold on
plot(freqHz,S2,'Color',[1, 0.85, 0.7],'LineWidth',2);
legend('Quiet','Attention')
%plot(freqHz,normsPxx(:,clusters==1), 'Color', [0.7, 0.7, 1, 0.5])
%plot(freqHz,normsPxx(:,clusters==2), 'Color', [1, 0.85, 0.7, 0.5])

set(gcf,'PaperPosition',[.25 3 8 6])
saveas(gcf,strcat(SA.currentPlotFolder, '/freqBands.pdf'));

%plot the freqs - not norelized::
figure;
set(gcf, 'Position', [300, 400, 700, 500]); % [x, y, width, height]
plot(freqHz,sPxx1,'Color',[0.7, 0.7, 1], 'LineWidth',2);
hold on
plot(freqHz,sPxx2,'Color',[1, 0.85, 0.7],'LineWidth',2);
set(gca, 'YScale', 'log');
title(sprintf('Ch%i',defCh))
legend('Quiet','Attention')

%% WAKE SESSIONS
SA = sleepAnalysis('/media/sil1/Data/Pogona Vitticeps/brainStatesWake.xlsx');
SA.setCurrentRecording(['Animal=PV157,recNames=Wake2']);

% SA.getFreqBandDetection('fMax', 110, )



%% FFT
welchWin = 1*fs; 
% get the data:
win = 5000; % in ms
fs = SA.currentDataObj.samplingFrequency(1); %samples per sec
% Bug apperance: before and after

postStart = SA.currentDataObj.getData(defCh,startTrigSh,win);
preStart = SA.currentDataObj.getData(defCh,startTrigSh-win, win);
%calculate the ppx:
[ppxPre,f1] = pwelch(squeeze(preStart(1,goodTrials,:)).',welchWin,welchWin/2,[],fs);
[ppxPost,f2] = pwelch(squeeze(postStart(1,goodTrials,:)).',welchWin,welchWin/2,[],fs);
meanPPXpre = mean(ppxPre,2);
meanPPXpost = mean(ppxPost,2);

%% normalizing each trial then average:
ntPPXpre = ppxPre./sum(ppxPre,1);
ntPPXpost = ppxPost./sum(ppxPost,1);
meanNTpre = mean(ntPPXpre,2);
meanNTpost = mean(ntPPXpost,2);
freqCutoff = 150;
%plot it:
figure;
plot(f1(1:freqCutoff),meanNTpre(1:freqCutoff,:),'Color','black',LineWidth=2);
hold on
plot(f2(1:freqCutoff),meanNTpost(1:freqCutoff,:),'Color','r',LineWidth=2)
xlabel ('Freq[Hz]'); ylabel('nPSD');
legend('Before bug apperance','After bug apperance')

titlestr = sprintf(['Normalized trials,data length: %d sec,' ...
    ' welch window: %d sec'],win/1000,welchWin/fs);
title(titlestr);
saveas(gcf,strcat(SA.currentPlotFolder, '/ntFreqBandsBugAppeare.jpg'));

%
ntmaenPPX = (meanNTpre + meanNTpost)/2;
normPPXPre = meanNTpre - ntmaenPPX;
normPPXPost = meanNTpost - ntmaenPPX;
figure;
plot(f1(1:freqCutoff),normPPXPre(1:freqCutoff,:),'Color','black',LineWidth=2);
hold on
plot(f2(1:freqCutoff),normPPXPost(1:freqCutoff,:),'Color','r',LineWidth=2)
xlabel ('Freq[Hz]'); ylabel('nPSD');
legend('Before bug apperance','After bug apperance')

titlestr = sprintf(['Normalized trials to mean,data length: %d sec,' ...
    ' welch window: %d sec'],win/1000,welchWin/fs);
title(titlestr);
saveas(gcf,strcat(SA.currentPlotFolder, '/tnFreqBandsBugAppeare.jpg'));
%% mark's analysis - not sure about it/
fMax=40;
p=find(f1<fMax);
ppPre=find(sum(ppxPre(p,:))<0.4e6); %reject signals with very high amplitudes (probably noise)
ppPost=find(sum(ppxPost(p,:))<0.4e6); %reject signals with very high amplitudes (probably noise)

sPxxPre=ppxPre(p,ppPre);
sPxxPost=ppxPost(p,ppPost);
sPxx=[sPxxPre sPxxPost];

freqHz=f1(p);
normsPxx=bsxfun(@rdivide,sPxx,mean(sPxx,2));
corrMat=corrcoef(normsPxx);

%maxDendroClusters=2;
%[DC,order,clusters]=DendrogramMamediantrix(corrMat,'linkMetric','euclidean','linkMethod','ward','maxClusters',maxDendroClusters);

figure;
plot(freqHz,median(normsPxx(:,1:25),2),'r');hold on;
plot(freqHz,median(normsPxx(:,26:end),2),'b')

%%
%plot:
fig = figure;

subplot (2,1,1)
plot(f1(1:60),meanPPXpre(1:60,:),'Color','r',LineWidth=2);
hold on
plot(f1(1:40),ppxPre(1:40,:))
title ('Before Bug Apperance')
legend('Average')
hold off

subplot (2,1,2)
plot(f2(1:60),meanPPXpost(1:60,:),'Color','black',LineWidth=2)
hold on
plot(f2(1:20),ppxPost(1:20,:))
title ('After Bug Apperance')
legend('Average')
hold off


han = axes(fig,'Visible','off');
han.XLabel.Visible = 'on';
xlabel(han,'Frequency[Hz]')
titlestr = sprintf(['Welch transform, time win of data = %d sec,' ...
    ' welch window = %d sec'],win/1000,welchWin/fs);
%han.Title.Visible = 'on';
%title(han, titlestr,);
sgtitle(titlestr);

% plot only averegaes on on figure:
figure;
plot(f1(1:60),meanPPXpre(1:60,:),'Color','black',LineWidth=2);
hold on
plot(f2(1:60),meanPPXpost(1:60,:),'Color','r',LineWidth=2)
legend('Before bug apperance','After bug apperance')
xlabel ('Freq[Hz]'); ylabel('PSD');
title('un-normelized');
saveas(gcf,strcat(SA.currentPlotFolder, '/FreqBandsBugAppeare.jpg'));
%% normelized: 
maenPPX = (meanPPXpre + meanPPXpost)/2;
normPPXPre = meanPPXpre - maenPPX;
normPPXPost = meanPPXpost - maenPPX;
figure;
plot(f1(1:100),normPPXPre(1:100,:),'Color','black',LineWidth=2);
hold on
plot(f2(1:100),normPPXPost(1:100,:),'Color','r',LineWidth=2)
xlabel ('Freq[Hz]'); ylabel('nPSD');
legend('Before bug apperance','After bug apperance')

titlestr = sprintf(['Normalized,data length: %d sec,' ...
    ' welch window: %d sec'],win/1000,welchWin/fs);
title(titlestr);
saveas(gcf,strcat(SA.currentPlotFolder, '/nFreqBandsBugAppeare.jpg'));
%% D2B ratio - swarm
%beta = [13, 30]; %hz
%delta = [0.1,3.5]; %hz
deltaBandCutoff = 5; 
betaBandLowcutoff = 10;
betaBandHighcutoff = 30;
betaFs = find(f1>=betaBandLowcutoff & f1<=betaBandHighcutoff);
deltaFs = find(f1<=deltaBandCutoff);

dbPre = [sum(ntPPXpre(deltaFs,:),1)] ./ [sum(ntPPXpre(betaFs,:),1)];
dbPost = [sum(ntPPXpost(deltaFs,:),1)] ./ [sum(ntPPXpost(betaFs,:),1)];

figure;
swarmchart(ones(1,length(dbPre)),dbPre, "filled");
hold on
swarmchart(2*ones(1,length(dbPost)),dbPost, "filled");
legend({'before','After'})
title('db ratio for each trial')
ylabel('ratio')
hold off
saveas(gcf,strcat(SA.currentPlotFolder, '/dbswarm.jpg'))
%% find the effect: FFT for only the 11 trials. 

welchWin = 3*fs; 

%calculate the ppx:
[ppxPre,f1] = pwelch(squeeze(preStartSlim).',welchWin,welchWin/2,[],fs);
[ppxPost,f2] = pwelch(squeeze(postStartSlim).',welchWin,welchWin/2,[],fs);
meanPPXpre = median(ppxPre,2);
meanPPXpost = median(ppxPost,2);

fig = figure;

subplot (2,1,1)
plot(f1(1:60),meanPPXpre(1:60,:),'Color','r',LineWidth=2);
hold on
plot(f1(1:40),ppxPre(1:40,:))
title ('Before Bug Apperance')
legend('Average')
hold off

subplot (2,1,2)
plot(f2(1:60),meanPPXpost(1:60,:),'Color','black',LineWidth=2)
hold on
plot(f2(1:20),ppxPost(1:20,:))
title ('After Bug Apperance')
legend('Average')
hold off


han = axes(fig,'Visible','off');
han.XLabel.Visible = 'on';
xlabel(han,'Frequency[Hz]')
titlestr = sprintf(['Welch transform, time win of data = %d sec,' ...
    ' welch window = %d sec'],win/1000,welchWin/fs);
%han.Title.Visible = 'on';
%title(han, titlestr,);
sgtitle(titlestr);


%% plot bug disappearing:

% get the data:
% instacnese:
% window- 3 seconds (for each part)
% channels - for now, I'm taking one that DO have spikes.
%   to have channels with spikes: we can chose 24,25,8,11,9. but many have. 

win = 3000; % in ms
fs = rec.samplingFrequency(1); %samples per sec
times = linspace(1,win,(win/1000)*fs); % time vector
defCh = 9;
% Bug apperance: before and after

postEnd = rec.getData(defCh,oeEndTrig,win);
preEnd = rec.getData(defCh,oeEndTrig-win, win);


%% plot traces:
% plot the averave trace of each trail
% avarge: IF ONLY ONE CH - NOTHING WILL HAPPEN
meanPostEnd = mean(postEnd,1);
meanPreEnd = mean(preEnd,1);

%ploting
figure
subplot(2,1,1);
plotShifted(times, squeeze(meanPreEnd).');
title('Before Bug Disapperance')
subplot(2,1,2);
plotShifted(times, squeeze(meanPostEnd).');
title('After Bug Disapperance')
xlabel ('Time(ms)')
sgtitle (sprintf("Traces of 25 trial, before/after bug Disapperance. Ch %d",defCh))

%% Get the Max amplitude
maxPreEnd = max(meanPreEnd,[],3);
maxPostEnd = max(meanPostEnd,[],3); 
meanMaxPre = mean(maxPreEnd);
meanMaxPost = mean(maxPostEnd);
stdMaxPre = std(maxPreEnd);
stdMaxPost = std(maxPostEnd);

% swarmchart:
figure

swarmchart(ones(1,length(maxPreEnd)),maxPreEnd, "filled");
hold on
swarmchart(2*ones(1,length(maxPostEnd)),maxPostEnd, "filled");
legend({'before','After'})
title('Max Values for each trial')
ylabel('voltage')
hold off


%% functions de
%% getFreqBandDetection
        function [data]=getFreqBandDetectionN(SA,varargin)
            SA.checkFileRecording;

            parseObj = inputParser;
            addParameter(parseObj,'ch',SA.recTable.defaulLFPCh(SA.currentPRec),@isnumeric);
            addParameter(parseObj,'fMax',45,@isnumeric); %max freq. to examine
            addParameter(parseObj,'dftPoints',2^10,@isnumeric);
            addParameter(parseObj,'tStart',0,@isnumeric);
            addParameter(parseObj,'win',1000*6,@isnumeric);
            addParameter(parseObj,'maxDendroClusters',2,@isnumeric);
            addParameter(parseObj,'saveFile',[]);
            addParameter(parseObj,'remove50HzArtifcats',false);
            addParameter(parseObj,'overwrite',0,@isnumeric);
            addParameter(parseObj,'segmentLength',1000);
            addParameter(parseObj,'WelchOL',0.5);
            addParameter(parseObj,'binDuration',10000);
            addParameter(parseObj,'inputParams',false,@isnumeric);
            addParameter(parseObj,'multiStart',false,@isnumeric);
            
            parseObj.parse(varargin{:});
            if parseObj.Results.inputParams
                disp(parseObj.Results);
                return;
            end
            
            %evaluate all input parameters in workspace
            for i=1:numel(parseObj.Parameters)
                eval([parseObj.Parameters{i} '=' 'parseObj.Results.(parseObj.Parameters{i});']);
            end
            
            if isnan(ch)
                error('LFP channel not define, either define in database as ''defaulLFPCh'' or as input to method ,eg ''ch'',''1''');
            end
            %make parameter structure
            parFreqBandDetection=parseObj.Results;
            
            %check if analysis was already done done
            if isempty(saveFile)
                SA.files.spectralClustering=[SA.currentAnalysisFolder filesep 'spectalClustering_ch' num2str(ch) '.mat'];
            else
                SA.files.spectralClustering=[saveFile '.mat'];
            end
            
            if exist(SA.files.spectralClustering,'file') & ~overwrite
                if nargout==1
                    data=load(SA.files.spectralClustering);
                else
                    disp('Spectral clustering analysis already exists for this recording');
                end
                return;
            end
            SA.getFilters;
            
            if ~multiStart
                if win>SA.currentDataObj.recordingDuration_ms-tStart
                    win=SA.currentDataObj.recordingDuration_ms-tStart;
                    fprintf('Window larger than recordings length, cutting window to %f [ms]\n',win);
                end
                win=floor(win/binDuration)*binDuration; %making win an integer number of segment length
            
                MLong=SA.currentDataObj.getData(ch,tStart,win);
                
            else
                if win>SA.currentDataObj.recordingDuration_ms-tStart(end)
                    win=SA.currentDataObj.recordingDuration_ms-tStart(end);
                    fprintf('Window larger than recordings length, cutting window to %f [ms]\n',win);
                end
                MLongM = SA.currentDataObj.getData(ch,tStart,win);
                MLong = reshape(MLongM,1,1,size(MLongM,2)*size(MLongM,3));
                win = length(MLong);
            
            end

            %filter data
            FMLong=SA.filt.F.getFilteredData(MLong);
    

            if remove50HzArtifcats
                SA.filt.notch=filterData(SA.filt.F.filteredSamplingFrequency);
                SA.filt.notch.filterDesign='cheby1';
                SA.filt.notch=SA.filt.notch.designNotch;
                SA.filt.notch.padding=true;
                FMLong=SA.filt.notch.getFilteredData(FMLong);
            end
            
            times=(tStart+binDuration/2):binDuration:(tStart+win);
            
            %calculate initial parameters
            segmentSamples = round(segmentLength/1000*SA.filt.FFs);
            samplesOL = round(segmentSamples*WelchOL);
            samplesBin = binDuration/1000*SA.filt.FFs;
            
            nBins=numel(FMLong)/samplesBin;

            FMLongB=reshape(FMLong,[samplesBin,nBins]);
            
            if (numel(FMLong)/samplesBin)~=round(numel(FMLong)/samplesBin)
                nBins=nBins-1;
                FMLong=FMLong(1:(samplesBin*nBins));
                disp('Last bin in recording not included due to a missmatch between recording duration and binDuration');
            end
                
            [pxx,f] = pwelch(FMLongB,segmentSamples,samplesOL,dftPoints,SA.filt.FFs);
            %plot(10*log10(pxx))
            p=find(f<fMax);
            pp=find(sum(pxx(p,:))<0.4e6); %reject signals with very high amplitudes (probably noise)
            
            sPxx=pxx(p,pp);
            freqHz=f(p);
            normsPxx=bsxfun(@rdivide,sPxx,mean(sPxx,2));
            corrMat=corrcoef(normsPxx);
            times=times(pp);
            if maxDendroClusters==2
                
                [DC,order,clusters]=DendrogramMatrix(corrMat,'linkMetric','euclidean','linkMethod','ward','maxClusters',maxDendroClusters);
                
                S1=mean(normsPxx(:,clusters==1),2);
                S2=mean(normsPxx(:,clusters==2),2);
                if mean(S1(1:3))>mean(S2(1:3))
                    crossFreq=freqHz(find(S2-S1>=0,1,'first'));
                else
                    crossFreq=freqHz(find(S1-S2>=0,1,'first'));
                end
            else
                [DC,order,clusters]=DendrogramMatrix(corrMat,'linkMetric','euclidean','linkMethod','ward','maxClusters',maxDendroClusters);
                
                for i=1:maxDendroClusters
                    S(:,i)=mean(normsPxx(:,clusters==i),2);
                end
                crossFreq=[];
            end
            
            save(SA.files.spectralClustering,'times','corrMat','sPxx','normsPxx','freqHz','parFreqBandDetection','order','clusters','crossFreq');
        end

