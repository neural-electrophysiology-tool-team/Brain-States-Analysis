function Rout = plotShiftedFigure(figName,rows,align,tGrid,minAtT,bandNames,yLab,col,shiftFactor,globalShift)
% PLOTSHIFTEDFIGURE  band-power time courses, all lines of a subplot shifted one above the other.
%
%   Rout = plotShiftedFigure(figName,rows,align,tGrid,minAtT,bandNames,yLab,col,shiftFactor,globalShift)
%
%   rows         : cell, one struct per row (fields name,no,st,minNo,minSt,ref,labels) - as built in the main script
%   align        : 'onset' or 'offset'
%   shiftFactor  : vertical distance between adjacent lines = shiftFactor * (largest peak-to-peak range of the column)
%                  (lower = the lines overlap more but look bigger).  default 0.5
%   globalShift  : false (default) -> one shift per band (column), identical in all rows
%                  true            -> one shift for all bands (then the ylim of all subplots of a row is identical)
%
%   Layout: the HEIGHT of each row is proportional to its number of lines, so that the distance between
%   adjacent lines (in z units AND on the screen) is the same throughout a whole column (and the same on the whole figure).
%   Each band has its own z scale (shift written in the column title).
%   y ticks (first column) = zero level of each line, labelled with the animal.
%   Dot: onset-aligned figure -> offset time ; offset-aligned figure -> onset time (only if inside the window).

    if nargin < 9 || isempty(shiftFactor), shiftFactor = 0.5; end
    if nargin < 10 || isempty(globalShift), globalShift = false; end
    figPos = [10 10 1850 1300];                         % figure size (pixels) - enlarge to see more detail

    nB = numel(bandNames);  nR = numel(rows);  nT = numel(tGrid);
    isOff = strcmp(align,'offset');
    sgn = 1 - 2*isOff;                                  % dot x = sgn * duration
    cA = {col.no, col.strike};  cM = {col.meanNo, col.meanStrike};  sfx = {'',' (strike)'};

    %% ---- 1. lines of every row (bottom -> top: animals, mean; then the same for strike) ----
    Ls = cell(1,nR);  anyStrike = false;
    for r = 1:nR
        row = rows{r};
        sets = {row.no, row.st};  mins = [row.minNo row.minSt];
        L = struct('data',{},'col',{},'lw',{},'ev',{},'label',{});
        for k = 1:2
            if isempty(sets{k}), continue; end
            Ak  = stackAligned(sets{k},mins(k),minAtT,nB,nT,align);
            evk = animalEv(sets{k},mins(k));
            haveA = find(reshape(any(any(~isnan(Ak),2),3),1,[]));
            for a = haveA
                L(end+1) = struct('data',reshape(Ak(a,:,:),nB,nT),'col',cA{k},'lw',1.2,'ev',evk(a), ...
                                  'label',[strrep(row.labels{a},'hunter','') sfx{k}]); %#ok<AGROW>
            end
            if ~isempty(haveA)
                if k == 2, anyStrike = true; end
                Mk = reshape(mean(Ak,1,'omitnan'),nB,nT);
                Mk(reshape(sum(~isnan(Ak),1),nB,nT) < 2) = NaN;        % mean over animals needs >= 2 animals
                if any(~isnan(Mk(:)))
                    L(end+1) = struct('data',Mk,'col',cM{k},'lw',3,'ev',median(evk,'omitnan'), ...
                                      'label',['mean' sfx{k}]); %#ok<AGROW>
                end
            end
        end
        Ls{r} = L;
    end

    %% ---- 2. shift of every column (same for all rows) ----
    rngB = zeros(1,nB);
    for r = 1:nR
        for i = 1:numel(Ls{r})
            d = Ls{r}(i).data;
            rngB = max(rngB,(max(d,[],2) - min(d,[],2))');
        end
    end
    shiftB = shiftFactor*rngB;
    shiftB(~isfinite(shiftB) | shiftB == 0) = 1;
    if globalShift, shiftB(:) = max(shiftB); end

    %% ---- 3. y-range of every row in units of the shift, row heights proportional to it ----
    lowM = zeros(1,nR);  hiM = zeros(1,nR);  span = 4*ones(1,nR);
    for r = 1:nR
        L = Ls{r};  nL = numel(L);
        if nL == 0, continue; end
        lm = 0;  hm = 0;
        for b = 1:nB
            lm = max(lm,-min(L(1).data(b,:),[],'omitnan')/shiftB(b));
            hm = max(hm, max(L(nL).data(b,:),[],'omitnan')/shiftB(b));
        end
        lowM(r) = max(0.3,lm+0.2);  hiM(r) = max(0.3,hm+0.2);
        span(r) = (nL-1) + lowM(r) + hiM(r);
    end

    %% ---- 4. figure with manual axes positions ----
    figure('Name',figName,'Color','w','Units','pixels','Position',figPos);
    mL = 0.09;  mR = 0.10;  mT = 0.08;  mB = 0.06;  gX = 0.012;  gY = 0.04;
    wAx  = (1 - mL - mR - (nB-1)*gX)/nB;
    avail = 1 - mT - mB - (nR-1)*gY;
    hAx  = avail*span/sum(span);
    annotation('textbox',[0 0.965 1 0.03],'String',[figName '     |     y: ' yLab], ...
               'EdgeColor','none','HorizontalAlignment','center','FontWeight','bold','Interpreter','none');

    Rout = cell(1,nR);  lastAx = [];  yTop = 1 - mT;
    for r = 1:nR
        row = rows{r};  L = Ls{r};  nL = numel(L);
        if isOff, refLab = row.ref{2}; else, refLab = row.ref{1}; end
        yBot = yTop - hAx(r);
        rowAx = gobjects(nB,1);
        for b = 1:nB
            ax = axes('Position',[mL+(b-1)*(wAx+gX), yBot, wAx, hAx(r)]);  hold(ax,'on');
            rowAx(b) = ax;  lastAx = ax;
            sh = shiftB(b);
            if nL > 0
                for i = 1:nL                                          % zero level of each line
                    plot(ax,[tGrid(1) tGrid(end)],(i-1)*sh*[1 1],':','Color',[.85 .85 .85]);
                end
                % ---- shifted lines (this loop can be replaced by your plotShifted call) ----
                for i = 1:nL
                    y = L(i).data(b,:) + (i-1)*sh;
                    plot(ax,tGrid,y,'-','Color',L(i).col,'LineWidth',L(i).lw);
                    plotDot(ax,tGrid,y,sgn*L(i).ev,L(i).col);
                end
                ylim(ax,[-lowM(r), (nL-1)+hiM(r)]*sh);
                yticks(ax,(0:nL-1)*sh);
                if b == 1, yticklabels(ax,{L.label}); ax.YAxis.FontSize = 7; else, ax.YTickLabel = {}; end
            else
                text(ax,0.5,0.5,'no data','Units','normalized','HorizontalAlignment','center');
                ax.YTick = [];
            end
            xline(ax,0,'-','Color',[.3 .3 .3]);
            xlim(ax,[tGrid(1) tGrid(end)]);  box(ax,'on');
            ax.XGrid = 'on';  ax.YGrid = 'off';  ax.FontSize = 8;
            if r == 1, title(ax,{char(bandNames(b)),sprintf('shift = %.2g z',sh)},'FontSize',9,'Interpreter','none'); end
            if r == nR, xlabel(ax,['Time from ' align ' (s)']); end
            if b == 1
                ylabel(ax,{row.name,['t = 0: ' refLab]},'FontWeight','bold','FontSize',9,'Interpreter','none');
            end
        end
        linkaxes(rowAx,'x');
        Rout{r} = struct('name',row.name,'t',tGrid,'shift',shiftB,'labels',{{L.label}},'data',cat(3,L.data));
        yTop = yBot - gY;
    end

    %% ---- legend ----
    h = [plot(lastAx,nan,nan,'-','Color',col.no), plot(lastAx,nan,nan,'-','Color',col.meanNo,'LineWidth',3)];
    lab = {'Animals, no strike','Mean, no strike'};
    if anyStrike
        h(3) = plot(lastAx,nan,nan,'-','Color',col.strike);
        h(4) = plot(lastAx,nan,nan,'-','Color',col.meanStrike,'LineWidth',3);
        lab = [lab {'Animals, strike','Mean, strike'}];
    end
    h(end+1) = plot(lastAx,nan,nan,'o','MarkerFaceColor','w','MarkerEdgeColor','k');
    if isOff, lab{end+1} = 'Onset time (median)'; else, lab{end+1} = 'Offset time (median)'; end
    lgd = legend(lastAx,h,lab);  lgd.Units = 'normalized';  lgd.FontSize = 8;
    lgd.Position = [1-mR+0.004, 0.45, mR-0.008, 0.12];
end


%% ---------------- local helpers ----------------
function A = stackAligned(cellS,minN,minAtT,nBands,nT,align)
% -> [nAnimals x nBands x nT] : per-animal mean over trials
    A = nan(numel(cellS),nBands,nT);
    for a = 1:numel(cellS)
        S = cellS{a};
        if numel(S.dur) < minN, continue; end
        if strcmp(align,'offset'), Zc = S.off; else, Zc = S.on; end
        P = cat(3,Zc{:});
        cnt = sum(any(~isnan(P),1),3);                            % 1 x nT trials contributing
        m = mean(P,3,'omitnan');  m(:,cnt < minAtT) = NaN;
        A(a,:,:) = m;
    end
end

function ev = animalEv(cellS,minN)
% median onset->offset duration per animal (NaN for animals below the minimum number of trials)
    ev = nan(numel(cellS),1);
    for a = 1:numel(cellS)
        if numel(cellS{a}.dur) >= minN, ev(a) = median(cellS{a}.dur); end
    end
end

function plotDot(ax,t,y,x,c)
% dot on the line at time x, only if x lies inside the plotted window
    if ~isfinite(x) || x < t(1) || x > t(end), return; end
    yv = interp1(t,y,x,'nearest');
    if isnan(yv), return; end
    plot(ax,x,yv,'o','MarkerFaceColor',c,'MarkerEdgeColor','k','MarkerSize',6,'LineWidth',0.8);
end