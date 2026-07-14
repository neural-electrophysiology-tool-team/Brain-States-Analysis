%% Sleep-Wake Analysis
clearvars
WA = wakeAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWake.xlsx');
analysisFolderWake = '/media/sil1/Data/Nitzan/WakeStatesPaper/plots';

%% transition points of states - WAKE
movWin = 5000;
movOLWin = 4500;

WA.setCurrentRecording('Animal=PV126,recNames=Wake5');

% DB=WA.getDelta2BetaRatio('overwrite',1, ...
    % 'movWin',5000, ...
    % 'movOLWin',4500, ...
    % 'win', 0, ...
    % 'tStart',0); %for shorther proceccing, until we decide on parameters. 

WA.getDelta2BetaAC
wSC = WA.getSlowCycles;
WA.plotDelta2BetaRatio;
WA.plotSlowCycles;

%% Transition point of states - SLEEP
WA.setCurrentRecording('Animal=PV126,recNames=Night3');
% tStartSleep = 4.5*60*60*1000;
% win= 1*60*60*1000;

DB=WA.getDelta2BetaRatio('overwrite',1, ...
    'movWin',5000, ...
    'movOLWin',4500);

WA.plotDelta2BetaRatio;
WA.getDelta2BetaAC('overwrite',1);
SC = WA.getSlowCycles;
WA.plotSlowCycles

%% get slow cycles - old
% parameters: 
medianFiltWin= 1000*5; %default = 1000*20
longOrdFiltWin = 1000*100; %defualt = 1000*1000
longOrdFiltOrd = 0.6; %default = 0.6
maxCycle= 1000*250; % defualt = 140*1000


WA.setCurrentRecording('Animal=PV126,recNames=Wake5');
wDB = WA.getDelta2BetaRatio;
wSCo=WA.getSlowCyclesOld( ...
    'overwrite',    1,...
    'medianFiltWin',medianFiltWin, ...
    'longOrdFiltWin', longOrdFiltWin, ...
    'longOrdFiltOrd', longOrdFiltOrd,...
    'maxCycle',     maxCycle ...
    );

% WA.plotSlowCycles;

% plot 

% Build state vector from cycle timestamps
t_ms = wDB.t_ms(:)';           % ensure row vector
dbRatio = wDB.bufferedDelta2BetaRatio(:)';

% Trim trailing NaNs (from both dbRatio and the corresponding time/state vectors)
lastValid = find(~isnan(dbRatio), 1, 'last');
t_ms      = t_ms(1:lastValid);
dbRatio   = dbRatio(1:lastValid);

nSamples = numel(t_ms);

% State codes: 1 = Unknown, 2 = SW, 3 = RE
stateVec = ones(1, nSamples);  % default = Unknown

nCycles = numel(wSCo.TcycleOnset);
for iCyc = 1:nCycles
    onset  = wSCo.TcycleOnset(iCyc);
    mid    = wSCo.TcycleMid(iCyc);
    offset = wSCo.TcycleOffset(iCyc);

    swIdx = t_ms >= onset & t_ms < mid;
    reIdx = t_ms >= mid   & t_ms < offset;

    stateVec(swIdx) = 2;
    stateVec(reIdx) = 3;
end

% Plot
figure('Color', 'w');

left  = 0.13;
width = 0.75;

% --- manual vertical layout (bottom, height) ---
ax2_bottom = 0.18; ax2_h = 0.18;   % smaller height to leave room for colorbar below
gap        = 0.12;                 % reserved space for colorbar between plots
ax1_bottom = ax2_bottom + ax2_h + gap; ax1_h = 0.35;

% --- Subplot 1: delta/beta ratio ---
ax1 = subplot('Position',[left ax1_bottom width ax1_h]);
plot(t_ms/1000, dbRatio, 'k', 'LineWidth', 1);
ylabel('\delta/\beta ratio');
title('Delta/Beta Ratio');
box off;

% --- Subplot 2: state imagesc ---
ax2 = subplot('Position',[left ax2_bottom width ax2_h]);
imagesc(t_ms/1000, 1, stateVec);
axis tight;
yticks([]);
xlabel('Time (s)');
title('Sleep State');

% Custom colormap: Unknown = light blue, SW = one color, RE = another
stateColors = [ ...
    0.68 0.85 0.90;   % 1 = Unknown - light blue
    0.85 0.33 0.10;   % 2 = SW - orange/red
    0.20 0.55 0.85];  % 3 = RE - blue
colormap(ax2, stateColors);
clim(ax2, [1 3]);

% Colorbar with labels, placed below ax2 (horizontal, doesn't shrink ax2 sideways)
cb = colorbar(ax2, 'southoutside', 'Ticks', [1.33 2 2.67], 'TickLabels', {'Unknown','SW','RE'});
drawnow;
cbWidth = width * 0.6;
cb.Position = [left + (width - cbWidth)/2, ...   % centered horizontally
               ax2_bottom - 0.09, ...              % sits below ax2
               cbWidth, ...
               0.02];

% restore ax2 position in case colorbar() shrank/moved it
ax2.Position = [left ax2_bottom width ax2_h];

% Link x-axes so zoom/pan stay synced
linkaxes([ax1, ax2], 'x');
xlim(ax1, [t_ms(1)/1000 t_ms(end)/1000]);

%% sleep - get slow cycle, old, same parameters as the wake. 

WA.setCurrentRecording('Animal=PV126,recNames=Night3');
SCo=WA.getSlowCyclesOld( ...
    'overwrite',    1,...
    'medianFiltWin',medianFiltWin, ...
    'longOrdFiltWin', longOrdFiltWin, ...
    'longOrdFiltOrd', longOrdFiltOrd,...
    'maxCycle',     maxCycle ...
    );

% WA.plotSlowCycles;
% Build state vector from cycle timestamps
t_ms = wDB.t_ms(:)';           % ensure row vector
dbRatio = wDB.bufferedDelta2BetaRatio(:)';

% Trim trailing NaNs (from both dbRatio and the corresponding time/state vectors)
lastValid = find(~isnan(dbRatio), 1, 'last');
t_ms      = t_ms(1:lastValid);
dbRatio   = dbRatio(1:lastValid);

nSamples = numel(t_ms);

% State codes: 1 = Unknown, 2 = SW, 3 = RE
stateVec = ones(1, nSamples);  % default = Unknown

nCycles = numel(wSCo.TcycleOnset);
for iCyc = 1:nCycles
    onset  = wSCo.TcycleOnset(iCyc);
    mid    = wSCo.TcycleMid(iCyc);
    offset = wSCo.TcycleOffset(iCyc);

    swIdx = t_ms >= onset & t_ms < mid;
    reIdx = t_ms >= mid   & t_ms < offset;

    stateVec(swIdx) = 2;
    stateVec(reIdx) = 3;
end

% Plot
figure('Color', 'w');

left  = 0.13;
width = 0.75;

% --- manual vertical layout (bottom, height) ---
ax2_bottom = 0.18; ax2_h = 0.18;   % smaller height to leave room for colorbar below
gap        = 0.12;                 % reserved space for colorbar between plots
ax1_bottom = ax2_bottom + ax2_h + gap; ax1_h = 0.35;

% --- Subplot 1: delta/beta ratio ---
ax1 = subplot('Position',[left ax1_bottom width ax1_h]);
plot(t_ms/1000, dbRatio, 'k', 'LineWidth', 1);
ylabel('\delta/\beta ratio');
title('Delta/Beta Ratio');
box off;

% --- Subplot 2: state imagesc ---
ax2 = subplot('Position',[left ax2_bottom width ax2_h]);
imagesc(t_ms/1000, 1, stateVec);
axis tight;
yticks([]);
xlabel('Time (s)');
title('Sleep State');

% Custom colormap: Unknown = light blue, SW = one color, RE = another
stateColors = [ ...
    0.68 0.85 0.90;   % 1 = Unknown - light blue
    0.85 0.33 0.10;   % 2 = SW - orange/red
    0.20 0.55 0.85];  % 3 = RE - blue
colormap(ax2, stateColors);
clim(ax2, [1 3]);

% Colorbar with labels, placed below ax2 (horizontal, doesn't shrink ax2 sideways)
cb = colorbar(ax2, 'southoutside', 'Ticks', [1.33 2 2.67], 'TickLabels', {'Unknown','SW','RE'});
drawnow;
cbWidth = width * 0.6;
cb.Position = [left + (width - cbWidth)/2, ...   % centered horizontally
               ax2_bottom - 0.09, ...              % sits below ax2
               cbWidth, ...
               0.02];

% restore ax2 position in case colorbar() shrank/moved it
ax2.Position = [left ax2_bottom width ax2_h];

% Link x-axes so zoom/pan stay synced
linkaxes([ax1, ax2], 'x');
xlim(ax1, [t_ms(1)/1000 t_ms(end)/1000])

%% state transition according to th (based on getslocycles old)
% parameters states:: 
medianFiltWin= 1000*5; %default = 1000*5
longOrdFiltWin = 1000*200; %defualt = 1000*200
longOrdFiltOrd = 0.4; %default = 0.4
minStateDur= 1000*7; % defualt = 7*1000

%parameters D/b
movWin = 5000;
movOLWin = 4500;

WA.setCurrentRecording('Animal=PV126,recNames=Wake5');
wDB=WA.getDelta2BetaRatio('overwrite',1, ...
    'movWin',movWin, ...
    'movOLWin',movOLWin);
WA.getDelta2BetaAC('overwrite',1)
wSS=WA.getStatesAC( ...
    'overwrite',    1,...
    'medianFiltWin',medianFiltWin, ...
    'longOrdFiltWin', longOrdFiltWin, ...
    'longOrdFiltOrd', longOrdFiltOrd,...
    'minStateDur',minStateDur...
    );

WA.plotStatesAC

%% sleep
WA.setCurrentRecording('Animal=PV126,recNames=Night3');
% DB=WA.getDelta2BetaRatio('overwrite',1, ...
%     'movWin',movWin, ...
%     'movOLWin',movOLWin);

% WA.getDelta2BetaAC('overwrite',1)
SS=WA.getStatesAC( ...
    'overwrite',    1,...
    'medianFiltWin',medianFiltWin, ...
    'longOrdFiltWin', longOrdFiltWin, ...
    'longOrdFiltOrd', longOrdFiltOrd,...
    'minStateDur',minStateDur...
    );
WA.plotStatesAC

