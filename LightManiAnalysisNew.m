%% light manipulation analysis 08/05/2024
% This is the new version. 

% SA=sleepAnalysis('/media/sil1/Data/Pogona Vitticeps/brainStatesWake.xlsx');

%%

%% geting the stim sham data from all records 

% new column in recTable - mani/ if 1 - manipulation night. if 2,
% manipulation night after third eye removal.
% 1. think of how to insert the different colors. 
% the ides is to: 1. go over all records with 1 or 2, get the stim sham.
% save them in a new variable (cell array?) .
% get also the anima number. once I have this structure I can start
% painting a pic

%% get all the stim sham avg from al recs anf put in a new var
SA=sleepAnalysis('/media/sil1/Data/Pogona Vitticeps/brainStatesWake.xlsx');

maniRecs = SA.recTable.Mani>0; % taking all the rows with manipulation
stimTable = SA.recTable(maniRecs,{'Animal','recNames','Remarks','Mani','StimTrighCh'});  % creating new table
stimTable.StimAvg = cell(height(stimTable),1);
stimTable.StimAvgSham = cell(height(stimTable),1);

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
    s = getStimSham(SA,stimTable.StimTrighCh(i));
    disp('got Stim for this rec')
    stimTable.StimAvg(i) = {mean(s.StimDB,1)};
    stimTable.StimAvgSham(i) = {mean(s.StimDBSham,1)};
    disp('stimsham in table')
    
end
clear recName
clear s

%% save stimTable
fileName = '/media/sil1/Data/Pogona Vitticeps/stimTable.m';
save(fileName, "stimTable",'-mat');




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
%     stimDuation=(endStim(1)-firstTrig(1));
    pre=50000;
    post=100000;
%     clear StimDB; %change to zeros
%     StimDB = zeros(1,numel(firstTrig));
    for i=1:numel(firstTrig)
        pTmp=find(DB.t_ms>(firstTrig(i)-pre) & DB.t_ms<=(firstTrig(i)+post));
        %StimDB(i,:)=1./DB.bufferedDelta2BetaRatio(pTmp);
        StimDB(i,:)=DB.bufferedDelta2BetaRatio(pTmp);
    end

    times=(DB.t_ms(pTmp)-DB.t_ms(pTmp(1)))/1000;
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
save(SA.files.stimSham,'StimDBSham','times','StimDB')
data.StimDBSham = StimDBSham;
data.times = times;
data.StimDB = StimDB;

end

function plotStimSham (rec)
    

    
end
