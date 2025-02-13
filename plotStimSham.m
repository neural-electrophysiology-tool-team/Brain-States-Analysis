

function plotStimSham(SA)
    % SA is an instance of sleep analysis class,with a record currently
    % selected
    stimShamFile=[SA.currentAnalysisFolder filesep  'stimSham.mat'];
    SA.checkFileRecording(stimShamFile,'stim Sham file missing, please first run getStimSham');
    load(stimShamFile); %load data
    
    trialsham = height(StimDBSham)-70+1;
    
    colorLim=[0 600];
    f=figure;
    subplot(4,2,[1:2:6]);imagesc(StimDBSham(trialsham:end,:),colorLim);ylabel('Trial #');title('Sham');hold on;set(gca,'XTick',[]);
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
    fileName=[SA.currentPlotFolder filesep 'stim_sham_activation.pdf'];
    saveas (f, fileName);
    % saveas (gcf, [analysisFolder filesep 'lizMovAllNights_norm.pdf']);

    
end

