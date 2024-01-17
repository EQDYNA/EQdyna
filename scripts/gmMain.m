close all;
% main program to postproess GM data.
% It calls gmGetRSA.m and gmPlotRSA.m.
% Dunyu Liu <dliu@ig.utexas.edu>, 20240109.
addpath('D:\GitHub\EQdyna\gm_utils_steve');

processRawGM = 1; % 1: yes, 2: no;
np = 10;
dx = 500; % in m
dt = 0.5*dx/6000;
nIntervalEQdyna = 10; % EQdyna writes out results every nIntervalEQdyna steps.
dt       = dt * nIntervalEQdyna;
gamma    = 0.05; % critical damping fraction
T        = [0.5, 1, 3]; 
range    = [-40, 40, -20, 20];
xLimit   = 25; % x coor of one fault end in km.
% to calculate Rupture distance from a station.
RBinSize = 1;  % km
Rmin     = 1;  % km
Rmax     = 20; % km

% selected stations to save acc seismograms; 
stCoorList = [0,15; 0,10; 0,5; 0,1; 0,0.1;
    5,15; 5,10; 5,5; 5,1; 5,0.1]; 

if processRawGM == 1
    gmGetRSA(dx, dt, gamma, np, T, range, stCoorList);
end
gmPlotRSA(dx, np, T, range, xLimit, RBinSize, Rmin, Rmax);