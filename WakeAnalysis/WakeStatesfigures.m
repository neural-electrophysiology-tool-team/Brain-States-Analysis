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


%% Check the video annotations:
videoPath = WA.recTable.VideoFiles{WA.currentPRec};
matPath = WA.files.ArenaLocation;
xColumnName = 'x__nose___cam_x__';
yColumnName = 'x__nose___cam_y__';

overlayPointsOnVideo(videoPath, matPath, xColumnName, yColumnName)

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

ArenaCSV = WA.getArenaCSVs;


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