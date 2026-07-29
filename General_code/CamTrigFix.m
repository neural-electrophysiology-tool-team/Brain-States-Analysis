% % cleaning up the camTrigs

SA.setCurrentRecording('');
vidTrig = SA.currentDataObj.getCamerasTrigger(7);

thresh= 15;

wrongT = find(diff(vidTrig)<thresh);
removeInd = wrongT(repmat([true,false],1,length(wrongT)/2));
newVidTrig = vidTrig;
newVidTrig(removeInd) = [];

isStart = [true,]

