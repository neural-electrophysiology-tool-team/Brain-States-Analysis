%% Spontaneous Wake analysis, wakefulness, 06.01.2026
WA = wakeAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWake.xlsx');
analysisFolderWake = '/media/sil1/Data/Nitzan/WakeStatesPaper/plots';

%% all wake trials 
wakeTrials = contains(WA.recTable.recNames,'Wake') & ...
             cellfun(@isempty,WA.recTable.Remarks);
wakeTable = WA.recTable(wakeTrials,:);

%% subset to analyse for the effect:
rng('shuffle')   % optional: for randomness

animals = unique(wakeTable.Animal);
rowsToKeep = [];

for i = 1:numel(animals)
    idx = find(strcmp(wakeTable.Animal,animals{i}));   % rows of this animal
    nPick = min(3, numel(idx));   % in case some animals have < 3 rows
    rowsToKeep = [rowsToKeep; idx(randperm(numel(idx), nPick))];
end

wakeSubset = wakeTable(rowsToKeep, :);
%%
wakeSubset = readtable([analysisFolderWake filesep 'wakeSubset.xlsx']);

%% checking correlogram and freq bands:

cali_result = load(['/media/sil3/Data/accelerometer_calibrations/' ...
    'headtagse_cali_recs/calibration_results.mat']).cali_result;


for i = 2:height(wakeSubset)
    recName = ['Animal=' wakeSubset.Animal{i} ',recNames=' wakeSubset.recNames{i}];
    WA.setCurrentRecording(recName);
    % WA.getFreqBandDetection(fMax=80,maxDendroClusters=3,overwrite=1)
    % WA.plotFreqBandDetection
    hs = WA.recTable.Headstage{WA.currentPRec};
    curSens = cali_result.(hs).sensetivity';
    curZeroG = cali_result.(hs).zeroGbais';
    curLM = WA.getLizardMovements(sensitivity=curSens,zeroGBias=curZeroG);
    DB = WA.getDelta2BetaRatio;


    %%plot ACC to DB
    convWin = 60;
    th= 100;
    dbInd = ~isnan(DB.bufferedDelta2BetaRatio);
    dbtimes= DB.t_ms(dbInd);
    dbratio = DB.bufferedDelta2BetaRatio(dbInd);
    binSize = DB.parDBRatio.movWin -DB.parDBRatio.movOLWin;
    edges = [dbtimes-(binSize/2), dbtimes(end)+(binSize/2)];
    [counts,e] =  histcounts(curLM.t_mov_ms,edges);
    
    h_thresh= counts>th;
    y = ones(1,convWin);
    conv = convn(h_thresh,y,'same');
    f=figure; 
    subplot(2,1,1);plot(conv);hold on;plot(dbratio); legend("mov","DB")
    subplot(2,1,2); scatter(conv,dbratio,'.');
    xlabel('Movement');ylabel("D/B");
    set(f,'PaperPositionMode','auto');
    fileName=[WA.currentPlotFolder filesep 'dbtoACC'];
    print(fileName,'-djpeg',['-r' num2str(WA.figResJPG)]);
    


    % k1 = 2*i-1; k2=2*1;
    % ax1 = nexttile(t,k1);   % dendrogram tile
    % % ax2 = nexttile(t,k2);   % spectra tile
    % WA.plotFreqBandDetection('plotDendrogram',true,'plotSpectralBands',true, ...
    %                      'hDendro',ax1,'hSpectra',ax2, ...
    %                      'savePlots',false);

end
%% figure out the states parameters:
WA.setCurrentRecording('Animal=PV126,recNames=Wake18');
data = WA.getFreqBandDetection(fMax=80,maxDendroClusters=3);
WA.plotFreqBandDetection;
% take the 50 most type 1 or type 2 from each and plot their pxx
type1i = data.order(1:50);
type2i = data.order(end-49:end);
%for timeseriesviewer:
type1t = data.times(type1i);
type2t = data.times(type2i);
%%
type1pxx = mean(data.sPxx(:,type1i),2);
type2pxx = mean(data.sPxx(:,type2i),2);
figure; plot(data.freqHz,type1pxx);hold on; plot(data.freqHz,type2pxx);legend(["quiet","active"])
title('pxx');xlabel('Freq[Hz]'),ylabel('Power')

type1npxx = mean(data.normsPxx(:,type1i),2);
type2npxx = mean(data.normsPxx(:,type2i),2);
figure; plot(data.freqHz,type1npxx);hold on; plot(data.freqHz,type2npxx);legend(["quiet","active"])
title('normelized pxx');xlabel('Freq[Hz]'),ylabel('nPower')
%% correlation with mov. - acc
% channels = 1:3; startT = 0;win = 1000*60*60*4;
% [V_uV,t_ms] = WA.currentDataObj.getAnalogData(channels,startT,win);
% V_uV = squeeze(V_uV)';


%
cali_result = load(['/media/sil3/Data/accelerometer_calibrations/' ...
    'headtagse_cali_recs/calibration_results.mat']).cali_result;
hs = WA.recTable.Headstage{WA.currentPRec};
curSens = cali_result.(hs).sensetivity';
curZeroG = cali_result.(hs).zeroGbais';
LM = WA.getLizardMovements(overwrite=1,sensitivity=curSens,zeroGBias=curZeroG);

%%
bin = 1000;
convWin = 60;
th= 100;
[counts, edges] =  histcounts(LM.t_mov_ms,BinWidth=bin);
h_thresh= counts>th;
% h_times = h.BinEdges;
y = ones(1,convWin);
conv = convn(h_thresh,y,'same');
figure; plot(conv);hold on;plot(DB.bufferedDelta2BetaRatio)
figure; scatter(conv(1:length(DB.bufferedDelta2BetaRatio)),DB.bufferedDelta2BetaRatio,'.');
xlabel('Movement');ylabel("D/B")
% move_conv = conv>0;
% db_stativ = DB.bufferedDelta2BetaRatio;

%% behaviors:
arena = WA.getArenaCSVs;
