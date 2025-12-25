%% is d.b increase is related to d/b before?
% figure;
% tiledlayout('flow')
HighLowTrialsAll=cell(height(stimTable),4);

for i = 1:height(stimTable)
    wavelength = stimTable.Remarks{i};
    if contains(wavelength, "Ex")
        continue
    end
    % % nexttile;
    % plotDBLowHighStim(SA,stimTable,i)
% end

%
% figure
% plotDBLowHighStim(SA,stimTable,i)
% function plotDBLowHighStim(SA,stimTable,i)
    % i = 82; % set the recording to PV153,N11
    recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
    HighLowTrialsAll(i,1)={recName};
    HighLowTrialsAll(i,2)=stimTable.Remarks(i);
    SA.setCurrentRecording(recName);
    stims = SA.getStimTriggers;
    firstTrigs = stims(1:8:end);
    pre = 30*1000;
    post= 100*1000;
    DB =SA.getDelta2BetaRatio;
    
    DBTrials =zeros([length(firstTrigs),(pre+post)/1000]);
    for j=1:length(firstTrigs)
        tStart = firstTrigs(j)-pre;
        tEnd = firstTrigs(j)+post;
        pTmp = DB.t_ms>tStart&DB.t_ms<tEnd;
        DBTrials(j,:) = DB.bufferedDelta2BetaRatio(pTmp);
    end
    
    % figure;
    % imagesc(DBTrials)
    % xline(30,'r')
    
    meanprewin = 10;%sec
    dbPre = mean(DBTrials(:,pre/1000-meanprewin:pre/1000),2);
    % figure; plot(dbPre)
    thresh = mean(dbPre);
    LowTrials=dbPre<thresh;
    HighTrials = dbPre>thresh;
    stimDur = round(stimTable.stimDuration(i)/1000);
    % win =30;
    dbLow = mean(DBTrials(LowTrials,(pre/1000):stimDur+pre/1000),1,'omitmissing');
    dbHigh = mean(DBTrials(HighTrials,(pre/1000):stimDur+pre/1000),1,'omitmissing');
    HighLowTrialsAll(i,3)={dbLow};
    HighLowTrialsAll(i,4)={dbHigh};
    
    % figure;
    % x = [ones(length(dbLow),1) ;2*ones(length(dbHigh),1)];
    % colors = [repmat([1 0.5 0],[length(dbLow),1]);repmat([0.5 0 0.5],[length(dbHigh),1])];
    % swarmchart(x,[dbLow;dbHigh],20,colors,'filled','XJitterWidth',0.5); %hold on;
    % scatter(2,dbHigh,"magenta",'filled','XJitterWidth',0.8);
    colors = [
    % Light RGB        Dark RGB
    0.80 0.80 0.80     0.40 0.40 0.40 ;   % Gray
    0.70 0.85 1.00     0.10 0.35 0.70 ;   % Blue
    0.70 0.90 0.70     0.20 0.55 0.20 ;   % Green
    0.95 0.70 0.70     0.60 0.10 0.10     % Red
    ];
    wavelength = stimTable.Remarks{i};

    if contains(wavelength,'47')| contains(wavelength,'lue')
        LowColor = colors(2,1:3);
        highColor = colors(2,4:6);
    elseif contains(wavelength,'53')| strcmp(wavelength,'green')
        LowColor = colors(3,1:3);
        highColor = colors(3,4:6);
    elseif contains(wavelength,'63')| strcmp(wavelength,'red')
        LowColor = colors(4,1:3);
        highColor = colors(4,4:6);
    elseif contains(wavelength,'DayTime') || contains(wavelength,'white')
        LowColor = colors(1,1:3);
        highColor = colors(1,4:6);
    end
    
    plot(dbLow,color=LowColor);hold on; plot(dbHigh,color=highColor);
    % xlabel('Time[s]');ylabel('d/b')
    % legend(["Low","High"])

    % xticks(1:2);xticklabels(["Low D2B","High D2B"])
    % xlim([0.5 2.5]); hold on;

    %stats and print:
    % if ~isempty(dbLow)& ~isempty(dbHigh)
    %     [p, ~] = ranksum(dbLow, dbHigh);
    %     if p<0.05
    %         scatter(1.5,200,50,"black","Marker","*")
    %     end
    %     % dim = [.2 .5 .3 .3];
    % 
    %     % annotation('textbox',[],'String',sprintf('%s\np value:%.4f\n',recName,p),FitBoxToText='on')
    %     % fprintf('p value diff between high and low trials:%.4f\n',p)
    %     fprintf('%s\np value:%.4f\n',recName,p)
    % 
    % end
    % set(gcf,'PaperPosition',[1 1 3.5 3]);
    % fileName=[SA.currentPlotFolder filesep 'LowHighDBbeforeStim'];
    % print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);
end
%%
wavelength = stimTable.Remarks;
redTrials= contains(wavelength,'63')| strcmp(wavelength,'red');
whiteTrials= contains(wavelength,'DayTime')| strcmp(wavelength,'white');

meanRedLow = mean(cell2mat(HighLowTrialsAll(redTrials,3)));