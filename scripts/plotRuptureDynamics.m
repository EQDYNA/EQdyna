%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The script plots on fault rupture dynamics from EQdyna.
% Input: frt.txt* from EQdyna.
% Last update: 20230511.
% Author: Dunyu Liu (dliu@ig.utexas.edu).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%
% adjustable parameters
path  = './';
x0    = 25e3; % half length of the fault along strike, m. 
z0    = 15e3; % vertical dimension of the fault, m
dx    = 100;  % grid size, m
np    = 1000; % max # of CPUs used.
font  = 12;
line  = 0:0.5:15; % contour info for rupture time contours.
write_file = 0; % 1/0 write/not cplot.txt for SCEC TPV benchmark comparison.
%%
na = 2*x0/dx+1; % # of on-fault grids along strike.
ma = z0/dx+1;   % # of on-fault grids along dip.

x=-x0:dx:x0;
z=-z0:dx:0;
[xx,zz]=meshgrid(x,z);

moment = 0;
for me=0:np-1
    fname=strcat('frt.txt',num2str(me));
    if exist(fname, 'file')
        a=load(fname);
        [n,m] = size(a);
        for i = 1:n
            ii   = (a(i,1)--x0)/dx+1;      % on-fault node id along strike  
            jj   = (a(i,3)--z0)/dx+1;      % on-fault node id along dip
            rupt(jj*na+ii,1) = a(i,1);
            rupt(jj*na+ii,2) = -a(i,3);
            rupt(jj*na+ii,3) = a(i,4);
            t(jj,ii)         = a(i,4);     % rupture time, s  
            slip(jj,ii)      = a(i,5);     % final slip, m
            slipr(jj,ii)     = a(i,10);    % peak slip rate, m/s
            fslipr(jj,ii)    = a(i,11);    % final slip rate, m/s
            tstk(jj,ii)      = a(i,13)/1e6;% final shear stress in MPa
            tnrm(jj,ii)      = a(i,12)/1e6;% final effective normal stress in MPa 
            
            %shear_mod = calc_shear_mod(abs(zz(jj,ii)), TPV);
            %moment = moment + slip(jj,ii)*shear_mod*dx^2;
        end
        delete a;
    end
end
%moment = moment * 1e7; % in dyne`cm
%mag = 2/3*log10(moment) - 10.7

% plot rupture time contour.
peak_slip = max(max(slip));
h = figure(1);
set(h,'Position', [30 30 1000 700]);

% creating 6 panel figure.

% panel 1, slip + rupture time contour
subplot(3,2,1);
%contourf(xx/1e3,zz/1e3,slip,'LineStyle','None');
%colorbar;
hold on; 
contour(xx/1e3,zz/1e3,t,line,'k-'); 
axis equal;
caxis([0 peak_slip]);
title('Rupture time contour (s) & Slip (m)');
ylabel('Dip (km)');
set(gca, 'FontSize',font,'FontWeight', 'Bold');

peak_slipr = max(max(slipr));
subplot(3,2,2);
contourf(xx/1e3,zz/1e3,slipr,'LineStyle','None'); 
colorbar; 
axis equal;
caxis([0 peak_slipr]);
title('Slip Rate (m/s)');
set(gca, 'FontSize',font,'FontWeight', 'Bold');

subplot(3,2,3);
contourf(xx/1e3,zz/1e3,fslipr,'LineStyle','None'); 
colorbar; 
axis equal;

title('Final Slip Rate (m/s)');
xlabel('Strike (km)');
ylabel('Dip (km)');
set(gca, 'FontSize',font,'FontWeight', 'Bold');

subplot(3,2,4);
contourf(xx/1e3,zz/1e3,tstk,'LineStyle','None'); 
colorbar;
axis equal;
caxis([min(min(tstk)) max(max(tstk))]);
title('Shear Stress Change (MPa)');
xlabel('Strike (km)');
ylabel('Dip (km)');
set(gca, 'FontSize',font,'FontWeight', 'Bold');

subplot(3,2,5)
contourf(xx/1e3,zz/1e3,tnrm,'LineStyle','None'); 
colorbar; 
axis equal;
caxis([min(min(tnrm)) max(max(tnrm))]);
title('Normal Stress Change (MPa)');
xlabel('Strike (km)');
set(gca, 'FontSize',font,'FontWeight', 'Bold');
set(gcf, 'color', 'white');

% write cplot.txt in the format for SCEC TPV benchmark.
if write_file == 1
    fileID = fopen('cplot.txt', 'w');
    fprintf(fileID,'j k t\n#\n');
    for i=1:na*ma
        fprintf(fileID,'%12.8f %12.8f %12.8f\n', rupt(i,1),rupt(i,2),rupt(i,3));
    end
    fclose(fileID);
end