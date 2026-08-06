% Hunter Analysis - old plots and trials. 
%% Check for strikes in the post bug appearance, and zero-pad

% It's doing zero padding but not in the right times. 
% i'm missing some data, in a few hundred ms but i dont konw why. 

for i = 1:length(oeStartTrial)
   for j = 1:length(oeStrikes)
       % check for each strike is it's time is inside the win: meaning its
       % larger then the begining but smaller than the end
        if oeStartTrial(i)<oeStrikes(j) && (oeStartTrial(i)+win)>oeStrikes(j)
            fprintf('during the win of Trial num %d there was a strike, running zero padding.\n',i)
            %find the time of strike, cut before and padd with zeros
            startToStrike = (oeStrikes(j)-oeStartTrial(i))*fs/1000; %match the times diff
            % change to zero:
            postStart(1,i,round(startToStrike:end)) = 0;
        end
   end
      
end

%% plot traces
% plot the traces:

% plot the averave trace of each trail
% avarge:
meanPostStart = mean(postStart,1);
meanPreStart = mean(preStart,1);

%ploting
figure
subplot(2,1,1);
plotShifted(times, squeeze(meanPreStart).','verticalShift',0.05);
title('Before Bug Apperance')
subplot(2,1,2);
plotShifted(times, squeeze(meanPostStart(1,1:24,:)).','verticalShift',0.05);
title('After Bug Apperance')
xlabel ('Time(ms)')
sgtitle (sprintf("Traces of 25 trial, before/after bug apperance. Ch %d",defCh))
%% extract only relevant parts to see the effect:
% preStartSlim = preStart(:,[1:2,6:10, 12:15],1:(end-5000)); %exclude: 3, 4, 17, ~11
% postStartSlim = postStart(:,[5:14],1:(end-5000));
% 
% %figure;plotShifted(squeeze(meanPreStart)');title('pre');
% %figure;plotShifted(squeeze(meanPostStart)','verticalShift',0.02);title('post')
% 
% figure
% subplot(2,1,1);
% plotShifted(times(1:end-5000), squeeze(preStartSlim).','verticalShift',0.05);
% title('Before Bug Apperance')
% subplot(2,1,2);
% plotShifted(times(1:end-5000), squeeze(postStartSlim).','verticalShift',0.05);
% title('After Bug Apperance')
% xlabel('Time(ms)')
%% get and plot the max amplitude for each trial ch mean - before and after. 

maxPreStart = max(meanPreStart,[],3);
maxPostStart = max(meanPostStart,[],3); 
%meanMaxPre = mean(maxPreStart);
%meanMaxPost = mean(maxPostStart);
%stdMaxPre = std(maxPreStart);
%stdMaxPost = std(maxPostStart);

% swarmchart:
figure

swarmchart(ones(1,length(maxPreStart)),maxPreStart, "filled");
hold on
swarmchart(2*ones(1,length(maxPostStart)),maxPostStart, "filled");
legend({'before','After'})
title('Max Values for each trial')
ylabel('voltage')
hold off

% maybe add the mean and std..
 