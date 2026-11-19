%% ============================================================
%  MEAN BAND POWER IN 2-s WINDOWS (5 freq. bands), per animal + mean
%
%  Windows (all winLen = 2 s, end point NOT included, relative to the reference event):
%     pre    : [preStart,    preStart    + winLen)  -> [-3,-1) s  ("baseline")
%     onset  : [onsetStart,  onsetStart  + winLen)  -> [-1, 1) s
%     offset : [offsetStart, offsetStart + winLen)  -> [ 1, 3) s
%     pre-offset (hunters) = the "pre" window taken relative to the OFFSET -> [-3,-1) s from offset
%
%  Z-scoring ("window", PAIRED): all windows of the SAME trial are concatenated in time, one mean/std
%  per band over all of them, then every window is averaged -> one value per band per window.
%  (A single window z-scored by itself would always average to exactly 0, so every figure uses the
%   joint z-score of its own trial; figures 4/5 simply re-use the onset/offset values of those trials.)
%
%  Per trial windows:
%     HUNTER circle           : [pre, onset, pre-offset, offset]   (pre-offset only for circle)
%     HUNTER other            : [pre, onset, offset]
%     Hunters are run for NO-STRIKE trials (blue) and STRIKE trials (orange) separately.
%     MEADOW (same table row) : [pre-move, onset-move, pre-stat, onset-stat]
%           move = Move_ms of the row;  stat = Triger_ms of the SAME row
%           pre-move = "pre" window before the move of this row  (inside this row's STAT period)
%           pre-stat = "pre" window before the stat of this row  (inside the previous row's MOVE period)
%     SOUND (wake only)       : [pre, onset, offset]; trial dropped if another sound (its 1-s interval)
%                               overlaps ANY of the 3 windows (same record)
%
%  Figures:
%   1. hunters, circle (row 1) / other (row 2): baseline | onset | [pre-offset (circle only)] | offset ; baseline subtracted
%   2. meadow: pre-move | onset-move | pre-stat | onset-stat ; pre-move subtracted
%   3. sound: baseline | onset | offset ; baseline subtracted
%   4. ONSET: circle bug | other bug | meadow move | sound ; no subtraction
%   5. OFFSET: circle bug | other bug | sound ; no subtraction (no meadow)
%  Columns = the 5 bands. No-strike: animals blue (dots joined by a line), mean black.
%  Strike: animals orange, mean gray. (Meadow / sound have no strike variant -> single set.)
% =============================================================
clear; clc;

%% ---------------- parameters ----------------
hunterTables = {'hunter106','hunter126','hunter143','hunter153','hunter157','hunter159','hunter161', ...
                'hunter162','hunter167','hunter208','hunter213','hunter228'};
strikeMinTime_s = 1;          % strike trials: Strike_ms > this (s)
minTrials      = 15;          % min valid trials per animal (no-strike hunters / meadow)
minTrialsStrike = 15;         % min valid strike trials per animal
minTrialsSound = 3;           % min valid sound events per animal
winLen = 2; preStart = -3; onsetStart = -1; offsetStart = 1;     % s
requirePreInPeriod = true;    % meadow: pre windows must lie inside the stat / previous-move period
normToSum = false;            % divide each band by the sum of the 5 bands before z-scoring
colNo = [0.20 0.45 0.80];  colStrike = [0.90 0.45 0.10];     % animals: no-strike blue, strike orange
colMeanNo = [0 0 0];       colMeanStrike = [0.5 0.5 0.5];    % means: no-strike black, strike gray

meadowVarFmt = 'reg_meadow%d';  soundVarName = 'reward_sound';  soundDur_ms = 1000;
soundRecordFilter = @(rec) contains(lower(string(rec)),'wake') & ~contains(lower(string(rec)),'hunter');
extraMatFiles = {};

bandRanges = [0 4; 5 12; 13 30; 31 40; 60 80];
bandNames  = ["Delta (0-4 Hz)","Theta+Alpha (5-12 Hz)","Beta (13-30 Hz)","GammaLow (31-40 Hz)","GammaHigh (60-80 Hz)"];
nBands = size(bandRanges,1);

loadData = load('/home/ofir/Documents/MATLAB/data4paper.mat');
for k = 1:numel(extraMatFiles)
    ex = load(extraMatFiles{k});
    for fn = fieldnames(ex)', loadData.(fn{1}) = ex.(fn{1}); end
end

recCache = containers.Map('KeyType','char','ValueType','any');
saCache  = containers.Map('KeyType','char','ValueType','any');

nA = numel(hunterTables);
VC = cell(nA,1); VO = cell(nA,1); VCs = cell(nA,1); VOs = cell(nA,1); VM = cell(nA,1); VS = cell(nA,1);   % [nBands x nWin x nTrials] per animal
animalLabels = cell(1,nA);

%% ---------------- extraction ----------------
for a = 1:nA
    hName = hunterTables{a};
    animalID = str2double(regexp(hName,'\d+','match','once'));
    animalLabels{a} = hName;

    if ismember(animalID,[208 213 228 167])
        excelPath = '/media/sil3/Data/Ofir/BrainStatesExperimentsOfir.xlsx';
    else
        excelPath = '/media/sil1/Data/Nitzan/Experiments/brainStatesWakeN.xlsx';
    end
    if ~isKey(saCache,excelPath), saCache(excelPath) = wakeAnalysis(excelPath); end
    X = struct('SA',saCache(excelPath),'bandRanges',bandRanges, ...
               'recCache',recCache,'normToSum',normToSum);

    % ---- hunters: circle / other ----
    if isfield(loadData,hName)
        H = loadData.(hName);
        for isS = [false true]
            for cType = {'circle','other'}
                tbl = bugTrials(H,cType{1},isS,strikeMinTime_s);
                if strcmp(cType{1},'circle')
                    % keep only trials whose pre-offset window does NOT overlap the onset window:
                    % pre-offset start (offset + preStart) >= onset-window end (onset + onsetStart + winLen)
                    okPO = (tbl.End_ms + preStart*1000) >= (tbl.Triger_ms  + (onsetStart + winLen)*1000);
                    fprintf('%s | circle | strike=%d : %d of %d trials removed (pre-offset overlaps onset window)\n', ...
                            hName,isS,sum(~okPO),height(tbl));
                    tbl = tbl(okPO,:);
                end
                E   = table(string(tbl.Animal),string(tbl.Record),'VariableNames',{'Animal','Record'});
                on  = tbl.Triger_ms ;  off = tbl.End_ms;
                if strcmp(cType{1},'circle')      % pre-offset window only for circle
                    refs   = [on on off off];
                    starts = [preStart onsetStart preStart offsetStart];
                else
                    refs   = [on on off];
                    starts = [preStart onsetStart offsetStart];
                end
                V = keepValid(computeWindows(E,refs,starts,winLen,X));
                fprintf('%s | %s | strike=%d : %d valid trials\n',hName,cType{1},isS,size(V,3));
                if strcmp(cType{1},'circle') && ~isS, VC{a} = V;
                elseif strcmp(cType{1},'circle'),     VCs{a} = V;
                elseif ~isS,                          VO{a} = V;
                else,                                 VOs{a} = V; end
            end
        end
    end

    % ---- meadow ----
    mName = sprintf(meadowVarFmt,animalID);
    if isfield(loadData,mName)
        [E,refs,starts] = meadowTrials(loadData.(mName),animalID,preStart,onsetStart,requirePreInPeriod);
        VM{a} = keepValid(computeWindows(E,refs,starts,winLen,X));
        fprintf('%s | meadow : %d valid trials\n',hName,size(VM{a},3));
    end

    % ---- sound (wake only) ----
    if isfield(loadData,soundVarName)
        [E,refs,starts] = soundTrials(loadData.(soundVarName),animalID,soundRecordFilter,soundDur_ms, ...
                                      [preStart onsetStart offsetStart],winLen);
        VS{a} = keepValid(computeWindows(E,refs,starts,winLen,X));
        fprintf('%s | sound : %d valid trials\n',hName,size(VS{a},3));
    end
end

%% ---------------- figures ----------------
% window columns: circle [1 pre,2 onset,3 pre-offset,4 offset] | other [1 pre,2 onset,3 offset]
% meadow [1 pre-move,2 on-move,3 pre-stat,4 on-stat] | sound [1 pre,2 onset,3 offset]
% each figure input: A{row} = {no-strike array, strike array ([] if none)}, arrays [nAnimals x nBands x nCond]
yl = 'Mean z-scored power (2-s window)';
col = struct('no',colNo,'strike',colStrike,'meanNo',colMeanNo,'meanStrike',colMeanStrike);
sk = @(V,cols,sub) stackAnimals(V,cols,sub,minTrialsStrike,nBands);
nk = @(V,cols,sub) stackAnimals(V,cols,sub,minTrials,nBands);
noD = nan(nA,nBands,1);

% Fig 1 : hunters, baseline subtracted
A1 = {{nk(VC,1:4,1), sk(VCs,1:4,1)}, {nk(VO,1:3,1), sk(VOs,1:3,1)}};
figureBands('Fig1 hunters (baseline subtracted)',{'CIRCLE','OTHER'},A1, ...
    {{'Baseline','Onset','Pre-offset','Offset'},{'Baseline','Onset','Offset'}}, ...
    animalLabels,bandNames,[yl ' - baseline'],col);

% Fig 2 : meadow (no strike distinction)
A2 = {{stackAnimals(VM,1:4,1,minTrials,nBands), []}};
figureBands('Fig2 meadow (pre-move subtracted)',{'MEADOW'},A2, ...
    {{'Pre move','Onset move','Pre stat','Onset stat'}},animalLabels,bandNames,[yl ' - pre move'],col);

% Fig 3 : sound (no strike distinction)
A3 = {{stackAnimals(VS,1:3,1,minTrialsSound,nBands), []}};
figureBands('Fig3 sound (baseline subtracted)',{'SOUND'},A3, ...
    {{'Baseline','Onset','Offset'}},animalLabels,bandNames,[yl ' - baseline'],col);

% Fig 4 : onsets, no subtraction (strike only exists for the two bug conditions)
A4no = cat(3,nk(VC,2,[]),nk(VO,2,[]),stackAnimals(VM,2,[],minTrials,nBands),stackAnimals(VS,2,[],minTrialsSound,nBands));
A4s  = cat(3,sk(VCs,2,[]),sk(VOs,2,[]),noD,noD);
figureBands('Fig4 onsets',{'ONSET'},{{A4no,A4s}},{{'Circle bug','Other bug','Meadow move','Sound'}}, ...
    animalLabels,bandNames,yl,col);

% Fig 5 : offsets, no subtraction, no meadow
A5no = cat(3,nk(VC,4,[]),nk(VO,3,[]),stackAnimals(VS,3,[],minTrialsSound,nBands));
A5s  = cat(3,sk(VCs,4,[]),sk(VOs,3,[]),noD);
figureBands('Fig5 offsets',{'OFFSET'},{{A5no,A5s}},{{'Circle bug','Other bug','Sound'}}, ...
    animalLabels,bandNames,yl,col);

save('lfp_window_power_results.mat','A1','A2','A3','A4','A5','animalLabels','bandNames', ...
     'winLen','preStart','onsetStart','offsetStart');


%% ============================================================
%  HELPERS
% =============================================================
function V = keepValid(V)
    if isempty(V), return; end
    ok = squeeze(~any(any(isnan(V),1),2));
    V = V(:,:,ok(:));
end

function A = stackAnimals(cellV,cols,subCol,minN,nBands)
% -> [nAnimals x nBands x numel(cols)] ; per-trial subtraction of window subCol (if given), then mean over trials
    nA = numel(cellV);  A = nan(nA,nBands,numel(cols));
    for a = 1:nA
        V = cellV{a};
        if isempty(V) || size(V,3) < minN, continue; end
        if ~isempty(subCol), V = V - V(:,subCol,:); end
        A(a,:,:) = mean(V(:,cols,:),3,'omitnan');
    end
end

function V = computeWindows(E,refs,starts,winLen,X)
% PAIRED window z-score: for every trial all windows (refs(i,w)+starts(w) ... +winLen s) are
% concatenated, one mean/std per band over them, then each window is averaged.  V: [nBands x nWin x nEv]
    nBands = size(X.bandRanges,1);  nW = numel(starts);  nEv = height(E);
    V = nan(nBands,nW,nEv);
    for i = 1:nEv
        animal = char(E.Animal(i));  record = char(E.Record(i));
        if any(isnan(refs(i,:))), continue; end
        R = getRec(X.recCache,X.SA,animal,record,X.bandRanges,X.normToSum);
        if isempty(R), continue; end
        idx = cell(1,nW);  ok = true;
        for w = 1:nW
            rel = (R.t_ms - refs(i,w))/1000;
            idx{w} = find(rel >= starts(w) & rel < starts(w)+winLen);
            if isempty(idx{w}), ok = false; break; end
        end
        if ~ok, continue; end
        both = R.bandRaw(:,[idx{:}]);
        mu = mean(both,2,'omitnan');  sd = std(both,0,2,'omitnan');  sd(sd==0) = NaN;
        for w = 1:nW
            V(:,w,i) = mean((R.bandRaw(:,idx{w}) - mu)./sd,2,'omitnan');
        end
    end
end

function [E,refs,starts] = meadowTrials(M,animalID,preStart,onsetStart,requireIn)
% windows per meadow row i (moveOn = Move_ms of row i, statOn = Triger_ms of the SAME row i):
%   [pre-move, onset-move, pre-stat, onset-stat]
    n = height(M);
    if ismember('Animal',M.Properties.VariableNames), an = string(M.Animal); else, an = repmat("PV"+animalID,n,1); end
    rec = string(M.Record);
    statOn = M.Triger_ms ;  moveOn = M.Move_ms;
    prevMove = [NaN; moveOn(1:end-1)];
    samePrev = [false; rec(2:end)==rec(1:end-1) & an(2:end)==an(1:end-1)];
    valid = ~isnan(statOn) & ~isnan(moveOn);
    if requireIn
        valid = valid & samePrev & (prevMove <= statOn + preStart*1000) ...   % pre-stat window inside previous MOVE period
                      & (statOn <= moveOn + preStart*1000);                    % pre-move window inside this STAT period
    end
    E = table(an(valid),rec(valid),'VariableNames',{'Animal','Record'});
    refs = [moveOn(valid) moveOn(valid) statOn(valid) statOn(valid)];
    starts = [preStart onsetStart preStart onsetStart];
end

function [E,refs,starts] = soundTrials(S,animalID,recFilter,dur_ms,starts,winLen)
% wake sounds of one animal; windows [pre, onset, offset]; trial dropped if another sound interval
% [onset, onset+dur] of the same record overlaps any of the windows
    num = cellfun(@(x) str2double(regexp(char(x),'\d+','match','once')),cellstr(string(S.Animal)));
    S = S(num(:)==animalID & logical(recFilter(S.Record)),:);
    rec = string(S.Record);  tOn = S.Triger_ms ;
    refsAll = [tOn tOn tOn+dur_ms];
    keep = true(height(S),1);
    for i = 1:height(S)
        oth = find(rec==rec(i));  oth(oth==i) = [];
        for w = 1:numel(starts)
            wS = refsAll(i,w) + starts(w)*1000;  wE = wS + winLen*1000;
            if any(tOn(oth) < wE & tOn(oth)+dur_ms > wS), keep(i) = false; break; end
        end
    end
    E = table(string(S.Animal(keep)),rec(keep),'VariableNames',{'Animal','Record'});
    refs = refsAll(keep,:);
end

function tbl = bugTrials(H,cType,isStrike,strikeMin_s)
    isCircle = (H.Type=="circle" | H.Type=="circle_flip");
    if strcmp(cType,'circle'), selType = isCircle; else, selType = ~isCircle; end
    if isStrike, tbl = H(selType & H.Strike_ms > strikeMin_s*1000,:);
    else,        tbl = H(selType & isnan(H.Strike_ms),:); end
end

function R = getRec(recCache,SA,animal,record,bandRanges,normToSum)
    key = [animal '|' record];
    if ~isKey(recCache,key), recCache(key) = loadBand(SA,animal,record,bandRanges,normToSum); end
    R = recCache(key);
end

function R = loadBand(SA,animal,record,bandRanges,normToSum)
    n = find(strcmp(string(SA.recTable.Animal),string(animal)) & strcmp(string(SA.recTable.recNames),string(record)));
    if isempty(n), warning('Recording %s / %s not found.',animal,record); R = []; return; end
    folder = fullfile(SA.recTable.folder{n},'analysis',['Animal=' animal ',recNames=' record]);
    ch = SA.recTable.defaulLFPCh(n);
    if strcmp(animal,'PV106') && ismember(record,{'Hunter19','Hunter20','Hunter21','Hunter34'}), ch = 5;
    elseif strcmp(animal,'PV157') && ismember(record,{'Hunter40','Hunter41'}), ch = 17; end
    fileName = [folder filesep 'bgRatio_ch' num2str(ch) '.mat'];
    if ~isfile(fileName), warning('File not found: %s',fileName); R = []; return; end
    S = load(fileName,'f','t_ms','allFreqProfiles');
    nB = size(bandRanges,1);  B = nan(nB,numel(S.t_ms));
    for fb = 1:nB
        fi = S.f >= bandRanges(fb,1) & S.f <= bandRanges(fb,2);
        if any(fi), B(fb,:) = mean(S.allFreqProfiles(fi,:),1,'omitnan'); end
    end
    if normToSum, B = B ./ sum(B,1,'omitnan'); B(~isfinite(B)) = NaN; end
    R.bandRaw = B;  R.t_ms = S.t_ms(:)';
end

function figureBands(figName,rowNames,A,condLabels,animalNames,bandNames,yLab,col)
% A{r} = {noStrike, strike} ; each [nAnimals x nBands x nCond] (strike may be [])
% no-strike: animals blue, mean black | strike: animals orange, mean gray (small x offset so they do not overlap)
    nR = numel(A);  nB = numel(bandNames);
    figure('Name',figName,'Color','w','Position',[50 50 1750 330*nR+120]);
    tl = tiledlayout(nR,nB,'TileSpacing','compact');
    hasStrike = false;
    for r = 1:nR
        nC = size(A{r}{1},3);
        for b = 1:nB
            ax = nexttile; hold(ax,'on');
            sets = {A{r}{1}, A{r}{2}};
            cA = {col.no, col.strike};  cM = {col.meanNo, col.meanStrike};  dx = [-0.06 0.06];
            if isempty(sets{2}), dx = [0 0]; end
            for k = 1:2
                if isempty(sets{k}), continue; end
                if k == 2, hasStrike = true; end
                nAn = size(sets{k},1);
                M = reshape(sets{k}(:,b,:),nAn,nC);
                for a = 1:nAn
                    plot(ax,(1:nC)+dx(k),M(a,:),'-o','Color',cA{k},'MarkerFaceColor',cA{k},'LineWidth',1.2,'MarkerSize',6);
                end
                if sum(any(~isnan(M),2)) > 1
                    plot(ax,(1:nC)+dx(k),mean(M,1,'omitnan'),'-s','Color',cM{k},'MarkerFaceColor',cM{k}, ...
                         'LineWidth',3,'MarkerSize',8);
                end
            end
            yline(ax,0,':','Color',[.5 .5 .5]);
            xlim(ax,[0.5 nC+0.5]); box(ax,'on'); grid(ax,'on');
            xticks(ax,1:nC); xticklabels(ax,condLabels{r}); xtickangle(ax,30); ax.TickLabelInterpreter = 'none';
            if r == 1, title(ax,bandNames(b),'FontSize',9); end
            if b == 1, ylabel(ax,rowNames{r},'FontWeight','bold'); end
        end
    end
    ylabel(tl,yLab); title(tl,figName,'Interpreter','none');
    h = plot(nan,nan,'-o','Color',col.no,'MarkerFaceColor',col.no);
    h(2) = plot(nan,nan,'-s','Color',col.meanNo,'MarkerFaceColor',col.meanNo,'LineWidth',3);
    lab = {'Animals, no strike','Mean, no strike'};
    if hasStrike
        h(3) = plot(nan,nan,'-o','Color',col.strike,'MarkerFaceColor',col.strike);
        h(4) = plot(nan,nan,'-s','Color',col.meanStrike,'MarkerFaceColor',col.meanStrike,'LineWidth',3);
        lab = [lab {'Animals, strike','Mean, strike'}];
    end
    lgd = legend(h,lab);  lgd.Layout.Tile = 'east';
end