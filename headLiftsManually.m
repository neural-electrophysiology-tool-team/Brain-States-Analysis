%% is there enough data to understand head lifts? 

% look at the angeles changesin the videos:
% check 5 recs:
% 1. PV161, night18 (probably not the right properties of the headstage)
recName = 'Animal=PV161,recNames=Night18';
SA.setCurrentRecording(recName);
LM = SA.getLizardMovements;
figure; plot(LM.t_static_ms, LM.angles(1,:)); 
figure; plot(LM.t_static_ms, LM.angles(2,:)); 
figure; plot(LM.t_static_ms, LM.angles(3,:)); 
diff1 = 0.9765-0.976;
diff2 = 2.1942-2.195;
diff3 = 2.1995-2.198; %biggest diff. tis is the critical channel.
% small bump in the middle of the sleep - invisible in the recording. (back
% camera). the head lift in the end is seen in the graph. at the same
% timings.
 %% 2. PV149, night12
recName = 'Animal=PV149,recNames=Night12';
SA.setCurrentRecording(recName);
LM = SA.getLizardMovements;
figure; plot(LM.t_static_ms, LM.angles(1,:)); 
figure; plot(LM.t_static_ms, LM.angles(2,:)); 
figure; plot(LM.t_static_ms, LM.angles(3,:)); 
% diff1 = ;
% diff2 = 2.1942-2.195;
% diff3 = 2.1995-2.198; %biggest diff between wake and sleep times. 
% .
% small bump in the middle of the sleep - invisible in the recording. (back
% camera). the head lift in the end is seen in the graph. at the same
% timings.
%% 3. PV159, night10 (probably not the right properties of the headstage)
recName = 'Animal=PV159,recNames=Night10';
SA.setCurrentRecording(recName);
LM = SA.getLizardMovements;
figure; plot(LM.t_static_ms, LM.angles(1,:)); 
figure; plot(LM.t_static_ms, LM.angles(2,:)); 
figure; plot(LM.t_static_ms, LM.angles(3,:)); 
% diff1 =  extreamlt small
diff2 = 2.1945-2.1955;% same and small
diff3 = 2.1995-2.1985;
% i don't understand which of the axis is the "right" one. 

%% 4. PV161, night13 - run the LM anallysis again. right parameters for headstage 
recName = 'Animal=PV161,recNames=Night13';
SA.setCurrentRecording(recName);
SA.getLizardMovements('overwrite',1)
LM = SA.getLizardMovements;
figure; plot(LM.t_static_ms, LM.angles(1,:)); 
figure; plot(LM.t_static_ms, LM.angles(2,:)); 
figure; plot(LM.t_static_ms, LM.angles(3,:)); % this is the right angle!
%  =  extreamlt small

%% 5. PV157, night16 -seems like the right headstage. 
recName = 'Animal=PV157,recNames=Night16';
SA.setCurrentRecording(recName);
% SA.getLizardMovements('overwrite',1)
LM = SA.getLizardMovements;
figure; plot(LM.t_static_ms, LM.angles(1,:)); 
figure; plot(LM.t_static_ms, LM.angles(2,:)); 
figure; plot(LM.t_static_ms, LM.angles(3,:)); % this is the right angle

%% 6. PV159, night34 -seems like the right headstage. 
recName = 'Animal=PV159,recNames=Night34';
SA.setCurrentRecording(recName);
SA.getLizardMovements('overwrite',1)
LM = SA.getLizardMovements;
figure; plot(LM.t_static_ms, LM.angles(1,:)); 
figure; plot(LM.t_static_ms, LM.angles(2,:)); 
figure; plot(LM.t_static_ms, LM.angles(3,:)); % this is the right angle

%% bottom line: 
