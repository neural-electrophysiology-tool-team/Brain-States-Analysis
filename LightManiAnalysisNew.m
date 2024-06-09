%% light manipulation analysis 08/05/2024
% This is the new version. 

% SA=sleepAnalysis('/media/sil1/Data/Pogona Vitticeps/brainStatesWake.xlsx');
% SA.batchProcessData('getDelta2BetaRatio',{})
%% analysis folder
analysisFolder = '/media/sil3/Data/Pogona_Vitticeps/NitzanAnalysisFiles';
%% geting the stim sham data from all records 

% new column in recTable - mani/ if 1 - manipulation night. if 2,
% manipulation night after third eye removal.
% 1. think of how to insert the different colors. 
% the ides is to: 1. go over all records with 1 or 2, get the stim sham.
% save them in a new variable (cell array?) .
% get also the anima number. once I have this structure I can start
% painting a pic


%% get all the stim sham avg from all recs and put in a new table
% this part goes over all the records in SA.
%  for every recoerd that is tagged (1/2/3..) 
SA=sleepAnalysis('/media/sil1/Data/Pogona Vitticeps/brainStatesWake.xlsx');
% SA.setCurrentRecording('Animal=PV162,recNames=Night27');
maniRecs = SA.recTable.Mani>0; % taking all the rows with manipulation
stimTable = SA.recTable(maniRecs,{'Animal','recNames','Remarks','Mani','StimTrighCh'});  % creating new table
stimTable.StimAvg = cell(height(stimTable),1);
stimTable.StimAvgSham = cell(height(stimTable),1);
stimTable.times = cell(height(stimTable),1);
stimTable.stimDuration = zeros(height(stimTable),1);

stimTable.ACbefore = cell(height(stimTable),1);
stimTable.ACduring = cell(height(stimTable),1);
stimTable.ACafter = cell(height(stimTable),1);

%% AC - under construction
i = 1;
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
SA.setCurrentRecording(recName);

% get the stimulation times:
T = SA.getDigitalTriggers();
trigs = T.tTrig{stimTable.StimTrighCh(i)};
firstStim = trigs(1);
lastStim = trigs(end);
ACbefore = SA.getDelta2BetaAC('tStart',10*1000*60,'win',3*60*60*1000,'overwrite',1);


%%

for i = 1:height(stimTable)
    % set te current rec:
    recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
    SA.setCurrentRecording(recName);
    % run all required analysis:
    SA.getDelta2BetaRatio;
    SA.getDelta2BetaAC;
    SA.getDigitalTriggers
    % get and save the stim avg
    s = getStimSham(SA,stimTable.StimTrighCh(i),1);
    disp('got Stim for this rec')
    stimTable.StimAvg(i) = {mean(s.StimDB,1)};
    stimTable.StimAvgSham(i) = {mean(s.StimDBSham,1)};
    stimTable.times(i) = {s.times};
    stimTable.stimDuration(i) = s.stimDur;
    disp('stimsham in table')
    

    
end
clear recName
clear s

%% save stimTable
save([analysisFolder filesep 'stimTable.mat'], "stimTable",'-mat');

%% PLOTINGS
%% plot the AC for nights manipulation. 
SA.setCurrentRecording('Animal=PV157,recNames=Night20')
SA.plotDelta2BetaSlidingAC


%% plot all stim sham, SlidingAc and D2B for stim recs:

load([analysisFolder filesep 'stimTable.mat'])
for i = 1:height(stimTable)
    % set te current rec:
    recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
    SA.setCurrentRecording(recName);
    SA.plotDelta2BetaRatio;
    SA.plotDelta2BetaSlidingAC;
    getStimSham(SA,stimTable.StimTrighCh(i),1);
    plotStimSham(SA);
%     keyboard;
end

%% get and plot D2B+SlidingAC for all nights for the five animals
SA=sleepAnalysis('/media/sil1/Data/Pogona Vitticeps/brainStatesWake.xlsx'); 
animals = unique(stimTable.Animal);                        
recList = {};
for i = 1:height(SA.recTable)
    animal = SA.recTable.Animal{i};
    recName = SA.recTable.recNames{i};
    if any(strcmp(animal,animals)) && contains(recName,'Night') && isnan(SA.recTable.Exclude(i))
        curRec=['Animal=' animal ',recNames=' recName];
        recList{end+1} = curRec;
    end
end
% disp(recList);
%% run batch analysis on all nights
SA.batchProcessData('getDelta2BetaRatio',recList)



%% PLOT stim activation according to stim type
% assuming the table is in the workspace
% load([analysisFolder filesep 'stimTable.mat'])

% loop on table
f = figure;
set(f, 'Position', [100, 100, 1200, 800]);

times = stimTable.times{1};
pre=50000;
% STIM DUTRATION NEEDS A THINK!!!!!!!
mStimDur = mean(stimTable.stimDuration(3:end));
post=100000;

%plot the blues:
ax1 = subplot(2,2,1);
hold on
blues = contains(stimTable.Remarks,'47');
meanBlue = mean(cell2mat(stimTable.StimAvg(blues)),1,'omitnan');
meanSham = mean(cell2mat(stimTable.StimAvgSham(blues)),1,'omitnan');
for i = 1:height(stimTable)
    if blues(i) && stimTable.Mani(i) == 1
        plot(times,stimTable.StimAvg{i},'blue')
    end 
    if blues(i) && stimTable.Mani(i) == 2
        plot(times,stimTable.StimAvg{i},'blue','LineStyle','--')
    end

end
plot(times,meanBlue,'color','blue','LineWidth',4)
plot(times,meanSham,'color','black','LineWidth',2)
xline(pre/1000, 'LineWidth',1,'Color','black')
xline((pre+mStimDur)/1000,'LineWidth',1,'Color','black')
title('blue laser stimulation');
xlabel('Time[S]'), ylabel('D2B power')
hold off

%plot the greens:
ax2 = subplot(2,2,2);
hold on
greens = contains(stimTable.Remarks,'532');
meangreen = mean(cell2mat(stimTable.StimAvg(greens)),1,'omitnan');


for i = 1:height(stimTable)
    if greens(i) && stimTable.Mani(i) == 1
        plot(times,stimTable.StimAvg{i},'green')
    end 
    if greens(i) && stimTable.Mani(i) == 2
        plot(times,stimTable.StimAvg{i},'green','LineStyle','--')
    end

end
meanSham = mean(cell2mat(stimTable.StimAvgSham(greens)),1,'omitnan');
plot(times,meanSham,'color','black','LineWidth',2)
plot(times,meangreen,'color','green','LineWidth',4)
xline(pre/1000, 'LineWidth',1,'Color','black')
xline((pre+mStimDur)/1000,'LineWidth',1,'Color','black')
title('Green laser stimulation');
xlabel('Time[S]'), ylabel('D2B power')
hold off


% plot the reds
ax3 = subplot(2,2,3);
hold on
reds = contains(stimTable.Remarks,'635');
meanreds = mean(cell2mat(stimTable.StimAvg(reds)),1,'omitnan');
for i = 1:height(stimTable)
    if reds(i) && stimTable.Mani(i) == 1
        plot(times,stimTable.StimAvg{i},'r')
    end 
    if reds(i) && stimTable.Mani(i) == 2
        plot(times,stimTable.StimAvg{i},'r','LineStyle','--')
    end

end
meanSham = mean(cell2mat(stimTable.StimAvgSham(reds)),1,'omitnan');
plot(times,meanSham,'color','black','LineWidth',2)
plot(times,meanreds,'color','r','LineWidth',4)
xline(pre/1000, 'LineWidth',1,'Color','black')
xline((pre+mStimDur)/1000,'LineWidth',1,'Color','black')
title('Red laser stimulation');

hold off
xlabel('Time[S]'), ylabel('D2B power')


% plot the LEDs
ax4 = subplot(2,2,4);
hold on
leds = contains(stimTable.Remarks,'LED');
meanleds = mean(cell2mat(stimTable.StimAvg(leds)),1,'omitnan');
gray = [0.5 0.5 0.5];
for i = 1:height(stimTable)
    if leds(i) && stimTable.Mani(i) == 1
        plot(times,stimTable.StimAvg{i},'color',gray)
    end 
    if leds(i) && stimTable.Mani(i) == 2
        plot(times,stimTable.StimAvg{i},'color',gray,'LineStyle','--')
    end

end
meanSham = mean(cell2mat(stimTable.StimAvgSham(leds)),1,'omitnan');
plot(times,meanSham,'color','black','LineWidth',2)
plot(times,meanleds,'color',gray,'LineWidth',4)
xline(pre/1000, 'LineWidth',1,'Color','black')
xline((pre+mStimDur)/1000,'LineWidth',1,'Color','black')
title('LED stimulation');

hold off


% general proprties
linkaxes([ax1,ax2,ax3,ax4],'xy')
ylims([0 200])
xlabel('Time[S]'), ylabel('D2B power')

% savefigure
saveas (gcf, [analysisFolder filesep 'all_colors_sham.jpg']);
clearvars -except stimTable SA analysisFolder

%% plot stim according to animal and color
% assuming the table is in the workspace
load([analysisFolder filesep 'stimTable.mat'])



animals = unique(stimTable.Animal);
stimType = ["Blue","Green","Red","LED"];
stimWaveL = ["47","532","635","LED"];
plotColors = {"blue","green","red", [0.5 0.5 0.5]};
numAnimal = length(animals);
numType = length(stimType);

times = stimTable.times{1};
pre=50000;
% STIM DUTRATION NEEDS A THINK!!!!!!!
mStimDur = mean(stimTable.stimDuration(3:end));
post=100000;

f=figure;
set(f, 'Position', [100, 100, 1200, 800]);

ylims = [0 300];
for animal = 1:numAnimal
    for type = 1:numType
        h= subplot(numAnimal,numType, (animal-1)*numType+type);
        
        %plot the data
        curAni = animals{animal};
        curType = stimWaveL(type);
        curTrials = contains(stimTable.Remarks,curType) & contains(stimTable.Animal,curAni);
        n = sum(curTrials);
        curCol = plotColors{type};
        curMean = mean(cell2mat(stimTable.StimAvg(curTrials)),1,'omitnan');
        if n>0
            plot(times,curMean,'color',curCol,'LineWidth',4)
            hold on
            for i = 1:height(stimTable)
                if curTrials(i) && stimTable.Mani(i) == 1
                    plot(times,stimTable.StimAvg{i},'Color',curCol)
                end
                if curTrials(i) && stimTable.Mani(i) == 2
                    plot(times,stimTable.StimAvg{i},'Color',curCol,'LineStyle','--')
                end

            end

            xline(pre/1000, 'LineWidth',1,'Color','black')
            xline((pre+mStimDur)/1000,'LineWidth',1,'Color','black')
            hold off
            
            annotation('textbox', [.07 + 0.21*type, 1.02 - 0.17*animal, 0.03, 0.1], 'String', ...
                sprintf('n=%i',n), 'EdgeColor', 'none', 'HorizontalAlignment', ...
                'right', 'VerticalAlignment', 'middle');

        end
        if n==0
            plot(0,0)
        end

        % add titles. labels...
        ylim(ylims);xlabel('Time[S]'), ylabel('D2B power')
        if animal == 1
            title(stimType(type))
        end
%         legend()
      
        pos = get(h, 'Position');
        pos(1) = pos(1) + 0.03; % Shift subplot right to make space for row titles
        set(h, 'Position', pos);
        
        % Add row titles using text function
        if type == 1
            ypos = (ylims(1) + ylims(2)) / 2; % Center vertically

            annotation('textbox', [0.1, 0.98 - 0.17*animal, 0.03, 0.1], 'String', animals{animal}, ...
                'EdgeColor', 'none', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
                'FontWeight', 'bold');
        
        end
    end 
end


% savefigure
saveas (gcf, [analysisFolder filesep 'animal_type.jpg']);
clearvars -except stimTable SA analysisFolder





%% getStimSham
% this functio is getting the data rellevant for the stimulation and the
% "sham" stimulation for the nights with stimulation. It should save the
% data in the analysis folder. 

function data = getStimSham(SA, t_ch, overwrite)
% SA is an instance of sleep analysis class,with a record currently
% selected
if nargin ==2
    overwrite = 0;
end
  %check if analysis was already done done
    SA.files.stimSham=[SA.currentAnalysisFolder filesep 'stimSham.mat'];
    if exist(SA.files.stimSham,'file') & ~overwrite
        if nargout==1
            data=load(SA.files.stimSham);
        else
            disp('stim sham analysis already exists for this recording');
        end
        return;
    end

    DB=SA.getDelta2BetaRatio;
    AC=SA.getDelta2BetaAC;
    T=SA.getDigitalTriggers;
    firstTrig=T.tTrig{t_ch}(1:8:end);
    endStim=T.tTrig{t_ch}(8:8:end)+400;
    stimDuration=(endStim(1)-firstTrig(1));
    pre=50000;
    post=100000;
%     clear StimDB; %change to zeros
%     StimDB = zeros(1,numel(firstTrig));
    for i=1:numel(firstTrig)
        pTmp=find(DB.t_ms>(firstTrig(i)-pre) & DB.t_ms<=(firstTrig(i)+post));
        %StimDB(i,:)=1./DB.bufferedDelta2BetaRatio(pTmp);
        StimDB(i,:)=DB.bufferedDelta2BetaRatio(pTmp);
    end

    ts=(DB.t_ms(pTmp)-DB.t_ms(pTmp(1)))/1000;
    meadStimInterval=mean(diff(firstTrig));
    firstTrigSham=(AC.tStartSleep:meadStimInterval:(firstTrig(1)-post))+10000;
    endStimSham=firstTrigSham+max(endStim-firstTrig);
    
%     clear StimDBSham;
    for i=1:numel(firstTrigSham)
        pTmp=find(DB.t_ms>(firstTrigSham(i)-pre) & DB.t_ms<=(firstTrigSham(i)+post));
        %StimDBSham(i,:)=1./DB.bufferedDelta2BetaRatio(pTmp);
        StimDBSham(i,:)=DB.bufferedDelta2BetaRatio(pTmp);
    end

   % save the data
save(SA.files.stimSham,'StimDBSham','ts','StimDB','stimDuration','pre','post')
data.StimDBSham = StimDBSham;
data.ts = ts;
data.StimDB = StimDB;
data.stimDur = stimDuration;
data.pre = pre;
data.post = post;

end

function plotStimSham(SA)
    % SA is an instance of sleep analysis class,with a record currently
    % selected
    stimShamFile=[SA.currentAnalysisFolder filesep  'stimSham.mat'];
    SA.checkFileRecording(stimShamFile,'stim Sham file missing, please first run getStimSham');
    load(stimShamFile); %load data
    
    
    colorLim=[0 600];
    f=figure;
    subplot(4,2,[1:2:6]);imagesc(StimDBSham,colorLim);ylabel('Trial #');title('Sham');hold on;set(gca,'XTick',[]);
    cb=colorbar('Position',[0.47 0.76 0.013 0.17]);ylabel(cb,'\delta/\beta');
    line([pre/1000 pre/1000],ylim,'color','r');
    subplot(4,2,7);plot(ts-pre/1000,nanmean(StimDBSham));xlabel(['Time [s]']);ylabel('Avg.');ylim(colorLim/3);
    line([0 0],ylim,'color','r');
    line([stimDuration/1000 stimDuration/1000],ylim,'color','r');
    subplot(4,2,[2:2:6]);imagesc(StimDB,colorLim);ylabel('Trial #');title('Stim');set(gca,'XTick',[]);
    cb=colorbar('Position',[ 0.91 0.76 0.013 0.17]);ylabel(cb,'\delta/\beta');
    line([pre/1000 pre/1000],ylim,'color','r');
    subplot(4,2,8);plot(ts-pre/1000,nanmean(StimDB));xlabel(['Time [s]']);ylabel('Avg.');ylim(colorLim/3);
    line([0 0],ylim,'color','r');
    line([stimDuration/1000 stimDuration/1000],ylim,'color','r');

%   save?
    fileName=[SA.currentPlotFolder filesep 'stim_sham_activation.jpg'];
    saveas (f, fileName);

    
end
