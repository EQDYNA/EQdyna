close all;
% selected stations to save acc seismograms; 
stCoorList = [0,15; 0,10; 0,5; 0,1; 0,0.1;
    5,15; 5,10; 5,5; 5,1; 5,0.1]; 
for i = 1: size(stCoorList,1)
    stNameList{i} = strcat('st.x',num2str(stCoorList(i,1)),'.y',num2str(stCoorList(i,1)));
end

dx = 500;
dt = 0.5*dx/6000; 
nIntervalEQdyna = 10; % EQdyna writes out results every nIntervalEQdyna steps.
dt = dt * nIntervalEQdyna;

figID = 1;
for k = 1: length(stNameList)
    stName = strcat(stNameList{k},'.accx.txt');
    acc = load(stName);
    gmStAccAnalysis(dt,acc,figID, stName); 
    figID = figID + 1;
    
    stName = strcat(stNameList{k},'.accy.txt');
    acc = load(stName);
    gmStAccAnalysis(dt,acc,figID, stName); 
    figID = figID + 1;
end