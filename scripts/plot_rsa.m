clear all; 
%close all;
% This script is used to plot RSA distribution from output of EQdyna.
% Created on 07/15/2021. 
% Author: Dunyu Liu (dliu@ig.utexas.edu).
% The script will call functions ResSpecRotD50.m and
% ResSpecTimeDomVectorized_3.m written by Steven Day.

TPV = 28026;
write_file = 0;

[path, x0, z0,np, dx, dt] = model_path_initial(TPV);

T = [0.5, 1, 3];
gamma = 0.05;
steps = 10;
dt = dt*steps;
NintMax = 5;
g = 9.8;

xs0 = -40; % km
xs1 = 40; % km
ys0 = -20; % km
ys1 = 20; % km

dx = dx/1000; % km

itag = 0;
for me=0:np-1
    fname=strcat(path,'/res_rsa/rsa',num2str(me),'.mat');
    if exist(fname, 'file')
        me
        itag = itag + 1;
        a = load(fname);
        tmp1 = a.coor;
        tmp2 = a.res_sa;
        if itag == 1
            res_sa = tmp2;
            coor = tmp1;
        else
            res_sa = [res_sa; tmp2;];
            coor = [coor; tmp1;];
        end
        clear a;
    end
end


[coor_new, id] =  unique(coor, 'rows');
for i=1:3
    tmp = res_sa(:,i);
    tmp1 = tmp(id);
    res_sa_new(:,i) = tmp1;
end

xx = xs0:dx:xs1; yy = ys0:dx:ys1;
xx = xx*1e3; yy = yy*1e3;
[x2,y2] = meshgrid(xx,yy);

h1=figure(1);
set(h1, 'Position', [10 10 500 800]);
for i = 1:3
    F = scatteredInterpolant(coor_new(:,1),coor_new(:,2), res_sa_new(:,i),'linear');
    sa_map = F(x2,y2);
    subplot(3,1,i)
    contourf(x2/1e3,y2/1e3,sa_map); colorbar;
    axis equal;
    if i == 1
        title('SA (g) at 0.5 s');
    elseif i == 2
        title('SA (g) at 1.0 s');ylabel('Fault-normal (km)');
    elseif i == 3
        title('SA (g) at 3.0 s');xlabel('Strike (km)');
    end
    set(gca, 'Fontsize', 12, 'Fontweight', 'bold');
end
set(gcf, 'color', 'white');



%%
[dmean, dmin, dmax, Rrup, period] = empirical_kevin;
k = [12,14,17]; % 0.5, 1, and 3 seconds.


h2=figure(2);
set(h2, 'Position', [10 10 600 900]);
for i = 1:3
    [stat, R] = sa_Rrupt(coor,res_sa(:,i));

    subplot(3,1,i)
    plot(Rrup,dmean(:,k(i)),'k-', 'LineWidth', 1); hold on;
    plot(Rrup,dmax(:,k(i)),'k:', 'LineWidth', 1); hold on;
    plot(Rrup,dmin(:,k(i)),'k:', 'LineWidth', 1); hold on;
    plot(R,stat(:,1),'o-','Linewidth',2);
    ylim([0.01 1]);
    if i == 1
        title('Spectral Acceleration(GMRotD50,g) at 0.5 s');
    elseif i == 2
        title('Spectral Acceleration(GMRotD50,g) at 1.0 s');
    elseif i == 3
        title('Spectral Acceleration(GMRotD50,g) at 3.0 s');
    end
    set(gca, 'YScale', 'log', 'XScale', 'log', 'color', 'white');
    set(gcf, 'color', 'white');
    set(gca, 'Fontsize', 12, 'Fontweight', 'bold');
end