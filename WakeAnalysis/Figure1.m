% Awake States Figures - figure 1, states

WA = wakeAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWake.xlsx');
analysisFolderWake = '/media/sil1/Data/Nitzan/WakeStatesPaper/plots';
clearvars -except WA analysisFolderWake wakeData wakeSubset clustersAll nPxxDataAll



%% Fig.1 A  - example raw data:
WA.setCurrentRecording('Animal=PV106,recNames=Wake1');
[V_uV,t_ms]=WA.currentDataObj.getData(WA.recTable.defaulLFPCh(WA.currentPRec),871400,16000);
V_uVd = downsample(squeeze(V_uV),100);
t_msd = downsample(t_ms,100);

f = figure;
set(f, 'Position', [100, 100, 400, 100]);

plot(t_msd/1000,V_uVd,'k'); hold on;
ylabel('\muV','Color','k')
xlabel('Time[s]','Color','k')
box off
ax = gca;
ax.YColor = 'k';
ax.XColor = 'k';

% Limits = [-1000 500]; % Get left y-axis limits

%save figure:
set(f,'PaperPositionMode','auto');
fileName=[analysisFolderWake filesep 'ExampleTraces_PV106W1'];
print(fileName,'-dpdf','-r300');
%%
WA.setCurrentRecording('Animal=PV106,recNames=Wake1');
fs = WA.currentDataObj.samplingFrequency(1);
s1raw = squeeze(V_uV(fs*4:fs*6-1));t1 = t_ms(fs*4:fs*6-1);
s2raw = squeeze(V_uV(fs*12:fs*14-1));t2 = t_ms(fs*12:fs*14-1);

f2 = figure;
set(f2, 'Position', [100, 100, 400, 100]);
subplot(1,2,1)
plot(t1/1000,s1raw,'k'); hold on;
ylabel('\muV','Color','k')
xlabel('Time[s]','Color','k')
box off
ax = gca; ax.YColor = 'k';ax.XColor = 'k';

subplot(1,2,2)
plot(t2/1000,s2raw,'k'); hold on;
ylabel('\muV','Color','k')
xlabel('Time[s]','Color','k')
box off
ax = gca;ax.YColor = 'k';ax.XColor = 'k';
%save figure:
set(f2,'PaperPositionMode','auto');
fileName=[analysisFolderWake filesep 'ExampleTraceszoomin_PV106W1'];
print(fileName,'-dpdf','-r300');

%% correlogram + frq bands - 1 example
WA.setCurrentRecording('Animal=PV106,recNames=Wake1');
WA.getFreqBandDetection(fMax=80, overwrite=1);
WA.plotFreqBandDetection;


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