% Awake States Figures - figure 1, states

WA = wakeAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWake.xlsx');
analysisFolderWake = '/media/sil1/Data/Nitzan/WakeStatesPaper/plots';
clearvars -except WA analysisFolderWake wakeData wakeSubset clustersAll nPxxDataAll

%% Fig.1 A  - example raw data:
WA.setCurrentRecording('Animal=PV126,recNames=Wake5');
% options: 
% Animal=PV126,recNames=Wake5: tStart4=94208, win =8000 / tStart = 1075153,
% win = 16000
% tStart = 1016000, win = 16000, no mov/ 
%%
tStart = 1074000; win = 20000;
ds = 100;
analog_ch = WA.recTable.accelerometerCh{WA.currentPRec};


% raw LFP:
[V_uV,t_ms]=WA.currentDataObj.getData(WA.recTable.defaulLFPCh(WA.currentPRec),tStart,win);
V_uVd = downsample(squeeze(V_uV),ds);
t_msd = downsample(t_ms,ds);

% ACC: 
% [lm, l_t_ms] = WA.currentDataObj.getAnalogData(str2num(analog_ch),tStart,win);
% lmd = downsample(squeeze(lm),ds);
% l_t_msd = downsample(l_t_ms,ds);
LM = WA.getLizardMovements;
moves = zeros(win,1);
pTmpMov = LM.t_mov_ms>tStart & LM.t_mov_ms<(tStart+win);
mov_t = LM.t_mov_ms(pTmpMov)-tStart;
moves(mov_t)=1;
moves_t = linspace(0,win/1000,win);
% pTmpSt = LM.t_static_ms>tStart & LM.t_static_ms<(tStart+win);
% static = LM.staticAll(pTmpSt);
% static_t = LM.t_static_ms(pTmpSt)-tStart;


% d/b: 
db=WA.getDelta2BetaRatio;%('overwrite',1,'movWin',2000,'movOLWin',1500);
pTmp = db.t_ms>tStart & db.t_ms<(tStart+win);
dbRatio = db.bufferedDelta2BetaRatio(pTmp);
dbRatio_t = db.t_ms(pTmp)-tStart;

% pupil diameter
r_eye_path ='/media/sil2/Data/Nimrod/experiments/PV_126/2024_07_18/block_009/analysis/right_eye_data_degrees_raw_verified.csv';
l_eye_path ='/media/sil2/Data/Nimrod/experiments/PV_126/2024_07_18/block_009/analysis/left_eye_data_degrees_raw_verified.csv';
r_eye = readtable(r_eye_path);
l_eye = readtable(l_eye_path);
pTmpREye = r_eye.ms_axis>tStart & r_eye.ms_axis<(tStart+win);
pTmpLEye = l_eye.ms_axis>tStart & l_eye.ms_axis<(tStart+win);
r_pup = r_eye.major_ax(pTmpREye);
l_pup = l_eye.major_ax(pTmpLEye);
r_pup_t = r_eye.ms_axis(pTmpREye)-tStart;
l_pup_t = l_eye.ms_axis(pTmpLEye)-tStart;

% head loc and angle:

%% Plot
f = figure;
set(f, 'Position', [500, 200, 600, 800]);

left  = 0.13;
width = 0.78;

% --- manual vertical layout (bottom, height) for each subplot ---
ax4_bottom = 0.08; ax4_h = 0.18;
gap1       = 0.05;
ax3_bottom = ax4_bottom + ax4_h + gap1;   ax3_h = 0.06;   % narrow: half height
gap2       = 0.05;
ax2_bottom = ax3_bottom + ax3_h + gap2;   ax2_h = 0.20;
gap3       = 0.09;                        % extra room reserved for colorbar
ax1_bottom = ax2_bottom + ax2_h + gap3;   ax1_h = 0.20;

% raw
ax1 = subplot('Position',[left ax1_bottom width ax1_h]);
plot(t_msd/1000,V_uVd,'k'); hold on;
ylabel('\muV','Color','k')
box off
ax1.YColor = 'k';
ax1.XColor = 'k';

% d/b
ax2 = subplot('Position',[left ax2_bottom width ax2_h]);
imagesc(dbRatio',"XData",dbRatio_t/1000)
ylabel('\delta/\beta','Color','k')
box off
ax2.YColor = 'k';
ax2.XColor = 'k';

cb = colorbar(ax2, 'northoutside');
drawnow;
cbWidth = width * 0.25;
cb.Position = [left + width - cbWidth, ...      % right-aligned
               ax2_bottom + ax2_h + 0.015, ...    % sits in the gap above ax2
               cbWidth, ...
               0.015];
cb.Label.String = '\delta/\beta';

% restore ax2 width in case colorbar() call shrank it
ax2.Position = [left ax2_bottom width ax2_h];
ax2.YColor = 'k';
ax2.XColor = 'k';
ax2.YTick = [];   % remove y-axis ticks

% Acc (narrow)
ax3 = subplot('Position',[left ax3_bottom width ax3_h]);
se = ones(1,15);
moves_disp = conv(double(moves(:)'), se, 'same') > 0;
imagesc(moves_t,[0 1],moves_disp);
colormap(ax3,[0.6 0.6 0.6; 0 0 0])   % grey = static, black = movement
ylim(ax3,[0 1])
ylabel('Acc','Color','k')
box off
ax3.YColor = 'k';
ax3.XColor = 'k';
ax3.YTick = [];   % remove y-axis ticks

% legend via dummy patches (imagesc has no native legend support)
h_static = patch(ax3, nan, nan, [0.6 0.6 0.6]);
h_moving = patch(ax3, nan, nan, [0 0 0]);
lg3 = legend(ax3, [h_static, h_moving], {'Static','Movement'}, ...
    'Orientation','horizontal', 'Box','off');
drawnow;
lg3.Position(1) = left + (width - lg3.Position(3))/2;   % center horizontally
lg3.Position(2) = ax3_bottom - lg3.Position(4) - 0.015;  % just below ax3
ax3.Position = [left ax3_bottom width ax3_h];  % restore width after legend


% pupils
ax4 = subplot('Position',[left ax4_bottom width ax4_h]);
plot(r_pup_t/1000,r_pup,'r'); hold on;
plot(l_pup_t/1000,l_pup,'b')
legend(["Right","Left"],"Location","southeast")
ylabel('Pupil Diameter','Color','k')
xlabel('Time [s]','Color','k')
box off
ax4.YColor = 'k';
ax4.XColor = 'k';

linkaxes([ax1,ax2,ax3,ax4],'x')
%save figure:
set(f,'PaperPositionMode','auto');
fileName=[analysisFolderWake filesep 'ExampleTraces_PV126W5'];
print(fileName,'-dpdf','-r300');
%% raw examples - short 
% WA.setCurrentRecording('Animal=PV106,recNames=Wake1');
% fs = WA.currentDataObj.samplingFrequency(1);
% s1raw = squeeze(V_uV(fs*4:fs*6-1));t1 = t_ms(fs*4:fs*6-1);
% s2raw = squeeze(V_uV(fs*12:fs*14-1));t2 = t_ms(fs*12:fs*14-1);
% 
% f2 = figure;
% set(f2, 'Position', [100, 100, 400, 100]);
% subplot(1,2,1)
% plot(t1/1000,s1raw,'k'); hold on;
% ylabel('\muV','Color','k')
% xlabel('Time[s]','Color','k')
% box off
% ax = gca; ax.YColor = 'k';ax.XColor = 'k';
% 
% subplot(1,2,2)
% plot(t2/1000,s2raw,'k'); hold on;
% ylabel('\muV','Color','k')
% xlabel('Time[s]','Color','k')
% box off
% ax = gca;ax.YColor = 'k';ax.XColor = 'k';
% %save figure:
% set(f2,'PaperPositionMode','auto');
% fileName=[analysisFolderWake filesep 'ExampleTraceszoomin_PV106W1'];
% print(fileName,'-dpdf','-r300');

%% correlogram + frq bands - 1 example
WA.setCurrentRecording('Animal=PV126,recNames=Wake5');
WA.getFreqBandDetection(...
'fMax',          80, ...
'segmentLength', 1000, ...
'binDuration',   5000, ...
'WelchOL',       0.75, ...
'overwrite',     1,...
'dftPoints',     2^10 ...
);
WA.plotFreqBandDetection

%% transtion pint and avg frq:
%% Plot transient point for the states.2 states.  (high/low)
crossFreqsAll=zeros(height(wakeData),1);
S1avg = zeros(height(wakeData),1);
S2avg = zeros(height(wakeData),1);
% freqdete = WA.getFreqBandDetection;
% freqHz = freqdete.freqHz;
for i= 1:height(wakeData)
    recName = ['Animal=' wakeData.Animal{i} ',recNames=' wakeData.recNames{i}];
    WA.setCurrentRecording(recName);
    % re-calc the 2 states:
    curFreq = WA.getFreqBandDetection(fMax=80,overwrite=1);
    crossFreqsAll(i) = curFreq.crossFreq;
    S1= mean(curFreq.normsPxx(:,curFreq.clusters==1),2);
    S2= mean(curFreq.normsPxx(:,curFreq.clusters==2),2);
    % get the avg of values above 1
    S1avg(i) = mean(curFreq.freqHz(S1>1));
    S2avg(i) = mean(curFreq.freqHz(S2>1));
end

%% plot - transition and avg for both states:
f = figure('Color','w');
colors = lines(4);
ax1=subplot(1,3,1:2);
x = ones(length(S1avg),1);
swarmchart(x,S1avg,'filled',color= colors(1,:));
hold on;
swarmchart(x*2,S2avg,'filled',color=colors(2,:));
ylabel("Frequency (Hz)")
xlim([0.5 2.5])
xticks([1,2]); 
xticklabels({'Cluster 1','Cluster 2'})
ax1.XTickLabelRotation=0;
ylims=[0,60];
box(ax1,'off');
ylim(ylims);box(ax1,'off');
% title('Clusters Avgs above 1')

ax2=subplot(1,3,3);
swarmchart(ones(size(crossFreqsAll)),crossFreqsAll,'filled','k');
ylim(ylims)
xticks(1);
xticklabels({"Cross Freq"})
xlim([0.5 1.5]);
% title('Cross Freq')
box(ax2,'off');
ax2.XTickLabelRotation=0;

% save:
set(f,'PaperPosition',[2 2 3 2]);
fileName=[analysisFolderWake filesep 'AvgCrossFreqAll'];
print(fileName,'-dpdf','-r300');

%% Silhouette histogram for 2 states:

%% evalute the clustering
k=2;
expNum = height(wakeData);
evaSilAll=zeros(expNum,1);

for i = 1:expNum
    curNormPxx=wakeData.pxxDataAll{i,2}; % all segments normPxx (nSegXnFreq)
    curCorrMat=corrcoef(curNormPxx);
    [DC,order,clusters]=DendrogramMatrix(curCorrMat,'linkMetric','euclidean','linkMethod','ward','maxClusters',k);
    s = silhouette(curCorrMat, clusters);  % or just silhouette(Xeval,idx)
    evaSilAll(i) = mean(s);
end

%% plot silhouette
f=figure('color','W');
histogram(evaSilAll,FaceColor=[0.5 0.5 0.5])
xlabel('Mean Silhouette');
ylabel('# of sessions')

% save:
set(f,'PaperPosition',[2 2 1.5 1.5]);
fileName=[analysisFolderWake filesep 'meanSilhistK2'];
print(fileName,'-dpdf','-r300');