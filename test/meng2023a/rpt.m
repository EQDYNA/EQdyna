clear all;close all;
% The script is used to plot rupture time contour from EQdyna 3D's output
% frt.txt*.
% The script was restructed on 07/04/2021.
% Author: Dunyu Liu (dliu@ig.utexas.edu).
% Last updated on 07/04/2021. 
% Change log:
%   * Add parameter TPV. Currently, parameter sets for TPV=104/105 are given.

load('PurpleColormaps.mat');
mymap = rgb;
font = 12;

TPV = 28026;
write_file = 0;

[path, x0, z0, np, dx, dt] = model_path_initial(TPV);

na = 2*x0/dx + 1;
ma = z0/dx+1;

x=-x0:dx:x0;
z=-z0:dx:0;
[xx,zz]=meshgrid(x,z);
line = 0:0.5:15;


moment = 0;
for me=0:np-1
    fname=strcat(path,'frt.txt',num2str(me));
    %me
    if exist(fname, 'file')
        a=load(fname);
        [n,m]=size(a);
        for i=1:n
            ii=(a(i,1)--x0)/dx+1;
            jj=(a(i,3)--z0)/dx+1;
            rupt(jj*na+ii,1)=a(i,1);
            rupt(jj*na+ii,2)=-a(i,3);
            rupt(jj*na+ii,3)=a(i,4);
            t(jj,ii)=a(i,4); 
            slip(jj,ii) = a(i,5);
            slipr(jj,ii) = a(i,10);
            tstk(jj,ii) = a(i,11)/1e6;
            tnrm(jj,ii) = a(i,12)/1e6;
            tstk0(jj,ii) = a(i,13)/1e6;
            tdip0(jj,ii) = a(i,14)/1e6;
            tnrm0(jj,ii) = a(i,15)/1e6;
            tstk(jj,ii) = tstk(jj,ii) - tstk0(jj,ii);
            tnrm(jj,ii) = tnrm(jj,ii) - tnrm0(jj,ii);
            rat(jj,ii) = -tstk0(jj,ii)/tnrm0(jj,ii);
            shear_mod = calc_shear_mod(abs(zz(jj,ii)), TPV);
            moment = moment + slip(jj,ii)*shear_mod*dx^2;
             %if tstk(jj,ii)>-15
             %    tstk(jj,ii) = -15;
             %end
%             if tnrm(jj,ii)<-200
%                 tnrm(jj,ii)=-200;
%             end
            %if a(i,4)>500
            %    t(jj,ii) = -1;
            %end
        end
        delete a;
    end
end
moment = moment * 1e7; % in dyne`cm
mag = 2/3*log10(moment) - 10.7

% plot rupture time contour.
peak_slip = max(max(slip));
h = figure(1);
set(h,'Position', [30 30 1000 600]);

subplot(2,2,1);
contourf(xx/1e3,zz/1e3,slip,'LineStyle','None');
colorbar; colormap(mymap);
hold on; 
contour(xx/1e3,zz/1e3,t,line,'w-'); 
axis equal;
caxis([0 peak_slip]);
title('Rupture time contour (s) & Slip (m)');
ylabel('Dip (km)');
set(gca, 'FontSize',font,'FontWeight', 'Bold');

peak_slipr = max(max(slipr));
subplot(2,2,2);
contourf(xx/1e3,zz/1e3,slipr,'LineStyle','None'); 
colorbar; colormap(mymap);
axis equal;
caxis([0 peak_slipr]);
title('Slip Rate (m/s)');
set(gca, 'FontSize',font,'FontWeight', 'Bold');

subplot(2,2,3);
contourf(xx/1e3,zz/1e3,tstk,'LineStyle','None'); 
colorbar; %colormap(parula)
axis equal;
caxis([min(min(tstk)) max(max(tstk))]);
title('Shear Stress Change (MPa)');
xlabel('Strike (km)');
ylabel('Dip (km)');
set(gca, 'FontSize',font,'FontWeight', 'Bold');

subplot(2,2,4)
contourf(xx/1e3,zz/1e3,tnrm,'LineStyle','None'); 
colorbar; %colormap(parula)
axis equal;
caxis([min(min(tnrm)) max(max(tnrm))]);
title('Normal Stress Change (MPa)');
xlabel('Strike (km)');
set(gca, 'FontSize',font,'FontWeight', 'Bold');
set(gcf, 'color', 'white');

h = figure(6);
set(h,'Position', [30 30 1000 500]);
contourf(xx/1e3,zz/1e3,tstk0); colorbar;
axis equal;
title('Normal stress (MPa)');
xlabel('Strike (km)');
ylabel('Dip (km)');
set(gca, 'FontSize',12,'FontWeight', 'Bold');
set(gcf, 'color', 'white');
%caxis([-200 0]);

% Whether to write cplot.txt.
if write_file == 1
    
    fileID = fopen('cplot.txt', 'w');
    fprintf(fileID,'j k t\n#\n');
    for i=1:na*ma
        fprintf(fileID,'%12.8f %12.8f %12.8f\n', rupt(i,1),rupt(i,2),rupt(i,3));
    end
    fclose(fileID);
end