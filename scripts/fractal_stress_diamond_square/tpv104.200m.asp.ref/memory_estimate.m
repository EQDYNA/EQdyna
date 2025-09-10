clear all; close all;
% The script is created on 08/26/2021 to estimate the memory used by EQdyna
% Dunyu Liu (dliu@ig.utexas.edu)

x = 70e3; 
y = 50e3;
z = 35e3;
dx = 25;
numel = x*y*z/dx/dx/dx

numnp = numel;
maxm = numel;
neq = numnp*3;
mem_per_node = 192;

mem = 4*(maxm+11*numnp);
mem = mem + 4* (8+48+1+28+8*4+6*8)*numel;
mem = mem + 4*242*200*200*1;
mem = mem + 4*numel + 40*maxm;
mem = mem + 32*neq + 16*3*numnp;

mem = mem/1024/1024

n = mem/mem_per_node