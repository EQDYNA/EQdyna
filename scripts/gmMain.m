close all;
% main program to postproess GM data.
% It calls gmGetRSA.m and gmPlotRSA.m.
% Dunyu Liu <dliu@ig.utexas.edu>, 20240109.
addpath('D:\GitHub\EQdyna\gm_utils_steve');

processRawGM = 1; % 1: yes, 2: no;
np = 4;
dx = 500; % in m
dt = 0.5*dx/6000;
nIntervalEQdyna = 10; % EQdyna writes out results every nIntervalEQdyna steps.
dt = dt * nIntervalEQdyna;
gamma = 0.05; % critical damping fraction
T  = [0.5, 1, 3, 5]; 
range = [-40, 40, -20, 20];
xLimit = 20; % x coor of one fault end in km.
% to calculate Rupture distance from a station.

if processRawGM == 1
    gmGetRSA(dx, dt, gamma, np, T, range);
end
gmPlotRSA(dx, np, T, range, xLimit);