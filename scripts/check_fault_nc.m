clear all; close all;
fname = 'fault.00041.nc'
a = ncread(fname ,'slip_rate');
a = a + 1e-50;
ss = ncread(fname ,'shear_strike');
sd = ncread(fname ,'shear_dip');
ns = ncread(fname ,'effective_normal');

dx = 250/1000;
x = -40:dx:40;
z = -40:dx:0;
[xx,zz] = meshgrid(x,z);

figure(1)
subplot(2,2,1)
contourf(xx,zz,log10(a'), 'LineColor','none'); colorbar;
title('log10(slip rate)');

subplot(2,2,2)
contourf(xx,zz,ss'/1e6,'LineColor','none'); colorbar;
title('Shear strike stress (MPa)');

subplot(2,2,3)
contourf(xx,zz,sd'/1e6, 'LineColor','none'); colorbar;
title('Shear dip stress (MPa)');

subplot(2,2,4)
contourf(xx,zz,ns'/1e6, 'LineColor','none'); colorbar;
title('Normal stress (MPa)');

[n,m] = size(a);
for i = 1:n
    for j = 1:m
        if a(i,j)>1e-6
            i,j
        end
    end
end
% peak = ncread('roughness.nc','peak');
% pypx = ncread('roughness.nc','pypx');
% pypz = ncread('roughness.nc','pypz');
% 
% figure(3)
% subplot(3,1,1)
% contourf(peak');
% subplot(3,1,2)
% contourf(pypx');
% subplot(3,1,3)
% contourf(pypz');