%% Wake States figures. 
% 8.6.2025
%% Inition parameters:
% SA=sleepAnalysis('/media/sil1/Data/Pogona Vitticeps/brainStatesWake.xlsx');
WA = wakeAnalysis('/media/sil1/Data/Pogona Vitticeps/brainStatesWake.xlsx');
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
LFPch = WA.recTable.defaulLFPCh(SA.currentPRec);
[LFP_uV, LFP_t] = WA.currentDataObj.getData(LFPch,bugApp-pre, pre+post);

% High pass - using the SA filters. 
% took the one with 80-2000 Hz. (FH)
WA.getFilters;
hplfp = WA.filt.FH.getFilteredData(LFP_uV); 

% getMovementdata
locTable = WA.getHeadPosition();
pTmp = find(locTable.t_ms >= bugApp-pre & locTable.t_ms <= bugApp+post);
headLocX = locTable.x__nose___x__(pTmp);
headLocY = locTable.x__nose___y__(pTmp);
headLoc_t = locTable.t_ms(pTmp)-(start_t);

%% Check the video annotations:
videoPath = WA.recTable.VideoFiles{WA.currentPRec};
matPath = WA.files.ArenaLocation;
xColumnName = 'x__nose___cam_x__';
yColumnName = 'x__nose___cam_y__';

overlayPointsOnVideo(videoPath, matPath, xColumnName, yColumnName)


%% plot
plotNum = 3;
figure;

bugAppT = pre; 
bugDisT = (bugDis-bugApp)+pre;

% LFP
subplot(plotNum,1,1)
plot(LFP_t/1000,squeeze(LFP_uV),'k');
xlabel('Time [s]'); ylabel('Voltage [uV]')
xline([bugAppT/1000,bugDisT/1000],'r');

% High pass
subplot(plotNum,1,2)
plot(LFP_t/1000,squeeze(hplfp),'k');
xlabel('Time [s]'); ylabel('Voltage [uV]')
xline([bugAppT/1000,bugDisT/1000],'r');

% head location - nose_y
subplot(plotNum,1,3)
plot(headLoc_t/1000,headLocY,'b');
xlabel('Time [s]'); ylabel('Head from screen (cm)')
xline([bugAppT/1000,bugDisT/1000],'r');
