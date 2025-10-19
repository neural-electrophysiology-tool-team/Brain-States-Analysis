%% getStimSham
% this functio is getting the data rellevant for the stimulation and the
% "sham" stimulation. It should save the data in the table. 

function data = getStimSham(SA, overwrite)
    % SA is an instance of sleep analysis class,with a record currently
    % selected
    if nargin ==2 
        overwrite = 0;
    end
    %check if analysis was already done
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
    
    stimData = SA.getStimTriggers;
    firstTrig = stimData(1:8:end-2);
    endStim = stimData(8:8:end)+200;
    stimDuration=(endStim(1)-firstTrig(1));

     
    pre=50000; %ms
    post=100000; %ms
%     clear StimDB; %change to zeros
    
    StimDB = zeros(numel(firstTrig),150);
    for i=1:numel(firstTrig)
        pTmp=find(DB.t_ms>(firstTrig(i)-pre) & DB.t_ms<=(firstTrig(i)+post));
        %StimDB(i,:)=1./DB.bufferedDelta2BetaRatio(pTmp);
        if length(pTmp) ~=150
            pTmp = ones([1,150]);
        end
        StimDB(i,:)=DB.bufferedDelta2BetaRatio(pTmp);
    end

    ts=(DB.t_ms(pTmp)-DB.t_ms(pTmp(1)))/1000;% creating Time array
    meadStimInterval=mean(diff(firstTrig)); %array of times between stims
    firstTrigSham=(AC.tStartSleep:meadStimInterval:(firstTrig(1)-post))+10000;
    endStimSham=firstTrigSham+max(endStim-firstTrig);
    
%     clear StimDBSham;
    for i=1:numel(firstTrigSham)
        pTmp=find(DB.t_ms>(firstTrigSham(i)-pre) & DB.t_ms<=(firstTrigSham(i)+post));
        %StimDBSham(i,:)=1./DB.bufferedDelta2BetaRatio(pTmp);
        if length(pTmp) ~=150
            continue
        end
        StimDBSham(i,:)=DB.bufferedDelta2BetaRatio(pTmp);
    end

   % save the data
save(SA.files.stimSham,'StimDBSham','ts','StimDB','stimDuration','pre','post')
data.StimDBSham = StimDBSham;
data.ts = ts;
data.StimDB = StimDB;
data.stimDuration = stimDuration;
data.pre = pre;
data.post = post;

end