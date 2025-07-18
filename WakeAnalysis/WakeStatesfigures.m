%% Wake States figures. 
% 8.6.2025
%% Inition parameters:
% SA=sleepAnalysis('/media/sil1/Data/Pogona Vitticeps/brainStatesWake.xlsx');
WA = wakeAnalysis('/media/sil1/Data/Pogona Vitticeps/brainStatesWake.xlsx');
analysisFolderWake = '/media/sil1/Data/Nitzan/WakeStatesPaper/plots';

%% Hunter example - all parameters:
% 			i. Body movement
% 			ii. Head movement (head angle)
% 			iii. Accelerometer (over threshold movement bins per second).
% 			iv. Visual stimulation
% 			v. LFP trace
% 			vi. High pass data
% 			vii. Raster
% 			viii. Pupil diameter
% *** Maybe show long epoch with zoom in on 3 segments.

% experiment : PV126, Hunter54
WA.setCurrentRecording('Animal=PV126,recNames=Hunter54');
ArenaCSV = WA.getArenaCSVs(WA.recTable.camTriggerCh(WA.currentPRec));


% Trial timings:
pre = 2000;
post = 8000;
trialNum = 15;
bugApp = ArenaCSV.startTrigSh(trialNum);
bugDis = ArenaCSV.endTrigSh(trialNum);
start_t = bugApp-pre;

% Get all data
% LFP:
LFPch = WA.recTable.defaulLFPCh(WA.currentPRec);
[LFP_uV, LFP_t] = WA.currentDataObj.getData(LFPch,bugApp-pre, pre+post);

% High pass - using the SA filters. 
% took the one with 80-2000 Hz. (FH)
WA.getFilters;
hplfp = WA.filt.FH.getFilteredData(LFP_uV); 

% getMovementdata
% locTable = WA.getHeadPosition();
% pTmp = find(locTable.t_ms >= bugApp-pre & locTable.t_ms <= bugApp+post);
% headLocX = locTable.x__nose___x__(pTmp);
% headLocY = locTable.x__nose___y__(pTmp);
% headLoc_t = locTable.t_ms(pTmp)-(start_t);

% pupil diameter. 
[rightData, leftData] = WA.getPupilData();
pTmp = find(rightData.ms_axis >= bugApp-pre & rightData.ms_axis <= bugApp+post);
rightPupil = rightData.major_ax(pTmp);
leftPupil = leftData.major_ax(pTmp);
pupil_t = leftData.ms_axis(pTmp)-start_t;

% accelerometer:
WA.getLizardMovements;
LM = WA.getLizardMovements;
%over threshold, bin to seconds. 
pTmp = find(LM.t_mov_ms >= bugApp-pre & LM.t_mov_ms <= bugApp+post);
mov = LM.movAll(pTmp);
mov_t = LM.t_mov_ms(pTmp)-start_t;

% raster data: 
%% plot Traces - 1 example
plotNum = 4;
f = figure;

bugAppT = pre; 
bugDisT = (bugDis-bugApp)+pre;

% LFP
ax1 = subplot(plotNum,1,1);
plot(LFP_t/1000,squeeze(LFP_uV),'k');
xlabel('Time [s]'); ylabel('Voltage [uV]')
xline([bugAppT/1000,bugDisT/1000],'r');

% High pass
ax2 = subplot(plotNum,1,2);
plot(LFP_t/1000,squeeze(hplfp),'k');
xlabel('Time [s]'); ylabel('Voltage [uV]')
xline([bugAppT/1000,bugDisT/1000],'r');

% pupil diamters
ax3 = subplot(plotNum,1,3);
pupilColors = [0 0.52 0.32;0.09 0.50 0.98];
plot(pupil_t/1000,rightPupil, Color=pupilColors(1,:));
hold on
plot(pupil_t/1000,leftPupil,Color=pupilColors(2,:));
xlabel('Time [s]'); ylabel('Pupil Diameter (px)')
xline([bugAppT/1000,bugDisT/1000],'r');
ylim([20 60])
legend({'Right Eye','Left Eye'},"Location",'southwest')

% acceleromerter
ax4 = subplot(plotNum,1,4);
plot(mov_t/1000,mov,Color='k',LineStyle='none',Marker='.',MarkerSize=4)
xlabel('Time [s]'); ylabel('Movement Amp.')
xline([bugAppT/1000,bugDisT/1000],'r');

% head location - nose_y
% subplot(plotNum,1,3)
% plot(headLoc_t/1000,headLocY,'b');
% xlabel('Time [s]'); ylabel('Head from screen (cm)')
% xline([bugAppT/1000,bugDisT/1000],'r');

% raster plot

linkaxes([ax1,ax2,ax3],'x')

% save plot
set(f,'PaperPosition',[1,1,5,4]);
fileName=[analysisFolderWake filesep 'ExampleTraces-PV126H54'];
print(fileName,'-dpdf','-r300')


%% example plot - old - only LFP and behavior - PV157,H5
% plot example trial - nose y-time, LFPs
% set the start of example timefrom rec + zeroFrame
WA.setCurrentRecording('Animal=PV157,recNames=Hunter5')
ArenaCSVs = WA.getArenaCSVs(1);
locTable = WA.getHeadPosition;


zeroTime = 82812.4; 
inds = find(ArenaCSVs.oeCamTrigs<zeroTime); 
zeroFrame = inds(end)+1;

% get the recdata 
fullwin = 14000;
[recData_M,recData_T] = WA.currentDataObj.getData(WA.recTable.defaulLFPCh(WA.currentPRec),zeroTime,fullwin);
recData_M = squeeze(recData_M);

% get the behavioral data:
frameWin = round(fullwin/(1000/ArenaCSVs.videoFPS));
nose_y = locTable.x__nose___y__(zeroFrame:zeroFrame+frameWin-1); 
nose_x = locTable.x__nose___y__(zeroFrame:zeroFrame+frameWin-1);
nose_t = linspace(0,fullwin/1000,frameWin);

% plot
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
fs = WA.currentDataObj.samplingFrequency(WA.recTable.defaulLFPCh(WA.currentPRec));
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
title('PV157,H5 - Example trial')

%saveFigure:
set(gcf,'PaperPosition',[.25 3 8 6])

saveas(gcf,strcat(WA.currentPlotFolder, '/TrialNoiseHL.pdf'));
fileName=[analysisFolderWake filesep 'ExampleTrailPV157,H5'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);





%% get beta to gamma ratio:
WA.setCurrentRecording('Animal=PV157,recNames=Hunter5');

movWin = 2000;
movOLWin = 1000;
segWel = 1000;
% the ratio is band2 to band1
WA.getMultiBandSpectralAnalysis('band1Low', 60,'band1High',80,'band2Low', ...
    10,'band2High',30, ...
    'maxVoltage',1000, 'tStart',0,'win',0, ...
    'movLongWin',1000*60*30,'movWin',movWin,'movOLWin',movOLWin,'segmentWelch' ...
    ,segWel,'overwrite',1);

BG = WA.getMultiBandSpectralAnalysis('band1Low', 60,'band1High',80,'band2Low', ...
    10,'band2High',30, ...
    'maxVoltage',1000, 'tStart',0,'win',0, ...
    'movLongWin',1000*60*30,'movWin',movWin,'movOLWin',movOLWin,'segmentWelch' ...
    ,segWel,'overwrite',0);

%get only non-Nan values.
p = find(~isnan(BG.band2to1Ratio));
BG.band2to1Ratio = BG.band2to1Ratio(p);
BG.t_ms = BG.t_ms(p);

%% Create the plot
figure;

% Plot the time series as an image
% imagesc(BG.t_ms/1000, 1, BG.band2to1Ratio');
% colorbar;
plot(BG.t_ms/1000,BG.band2to1Ratio, color='k',LineStyle='-')
hold on; % Important: allows overlaying additional plots


n = length(ArenaCSV.endTrigSh);
y_height = ylim();
x_all = [ArenaCSV.startTrigSh'; ArenaCSV.endTrigSh'; ArenaCSV.endTrigSh'; ArenaCSV.startTrigSh'; ArenaCSV.startTrigSh']/1000;  % 5�n matrix
y_all = [repmat(y_height(1),1,n);repmat(y_height(1),1,n); repmat(y_height(2),1,n); repmat(y_height(2),1,n); repmat(y_height(1),1,n)];

fill(x_all, y_all, 'red', 'FaceAlpha', 0.5, 'EdgeColor', 'none');


% Labels and formatting
xlabel('Time [s]');

% Optional: Add legend
legend('B/G', 'Trial Time', 'Location', 'best');

hold off;





%% Check the video annotations:
videoPath = WA.recTable.VideoFiles{WA.currentPRec};
matPath = WA.files.ArenaLocation;
xColumnName = 'x__nose___cam_x__';
yColumnName = 'x__nose___cam_y__';

overlayPointsOnVideo(videoPath, matPath, xColumnName, yColumnName)
