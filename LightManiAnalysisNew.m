%% light manipulation analysis 08/05/2024
% This is the new version. 

SA=sleepAnalysis('/media/sil1/Data/Pogona Vitticeps/brainStatesWake.xlsx');
SA.setCurrentRecording('Animal=PV157,recNames=Night29');
getStimSham(SA);

%% getStimSham
% this functio is getting the data rellevant for the stimulation and the
% "sham" stimulation for the nights with stimulation. It should save the
% data in the analysis folder. 

function data = getStimSham(SA, overwrite)
% SA is an instance of sleep analysis class,with a record currently
% selected
if nargin ==1
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
    firstTrig=T.tTrig{15}(1:8:end);
    endStim=T.tTrig{15}(8:8:end)+400;
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
