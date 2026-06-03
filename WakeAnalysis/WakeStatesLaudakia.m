%% Wake Analysis in Laudakia

WAs = wakeAnalysis('/media/sil2/Data/Lizard/Stellagama/brainStatesSS.xlsx');
analysisFolderWake = '/media/sil1/Data/Nitzan/WakeStatesPaper/plots';
%%
WAs.setCurrentRecording('Animal=SA07,recNames=sleepNight5');
% WAs.plotDelta2BetaRatio;

% claculate LM for all rec
LM= WAs.getLizardMovements;
% get the sleep timings
DB = WAs.getDelta2BetaRatio;
AC = WAs.getDelta2BetaAC;
p=30*60*1000; % 30 mins after wake

wakeTstart = AC.tEndSleep+p;
%% get the DB and ACC for the wake time: 
    %%plot ACC to DB
    % convWin = DB.parDBRatio.movWin/1000;
    convWin=60;
    th= 0;
    tStart = wakeTstart;
    wakeOnly=1;


    bin = DB.parDBRatio.movWin - DB.parDBRatio.movOLWin;
    DBt = DB.t_ms(~isnan(DB.bufferedBetaRatio)); % remove nan
    buffDB=DB.bufferedDelta2BetaRatio(~isnan(DB.bufferedBetaRatio));
    %take only wake times:
    if wakeOnly
        pTmp= DBt>tStart;
        DBt = DBt(pTmp);
        buffDB=buffDB(pTmp);

    end
    edges = [DBt-bin/2, DBt(end)+bin/2];
    [counts, edges] =  histcounts(LM.t_mov_ms,edges);
    % h_thresh= counts>th;

    h_thresh=zeros(size(counts));
    h_thresh(counts>th)= counts(counts>th);

    % h_times = h.BinEdges;
    y = ones(1,convWin);
    conv = convn(h_thresh,y,'same');
    % figure; plot(conv);hold on;plot(DB.bufferedDelta2BetaRatio)
    
    xdata=conv;
    ydata=buffDB;
% hist 2 plot:
figure;
[d,h] = hist2(xdata,ydata,'dX1',50,'dX2',10);
% colorbar
colormap(parula)
xlabel('Mov Count')
ylabel('D/B')



%% for all animals anfd trials:
