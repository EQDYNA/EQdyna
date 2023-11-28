clear all; close all;
% This script is used to plot RSA distribution from output of EQdyna.
% Created on 07/15/2021. 
% Author: Dunyu Liu (dliu@ig.utexas.edu).
% The script will call functions ResSpecRotD50.m and
% ResSpecTimeDomVectorized_3.m written by Steven Day.
me1 = 0;
me2 = 200;

TPV = 1043;
write_file = 0;

%[path,x0,z0,np,dx] = model_path_initial(TPV);
path = './'

T = [0.5, 1, 3];
gamma = 0.05;
steps = 10;
dt = 0.004*steps;
NintMax = 5;
g = 9.8;

x0 = -40; % km
x1 = 40; % km
y0 = -40; % km
y1 = 41; % km

dx = 100/1000; % km

for me=me1:me2
    ntag = 0;
    fname=strcat(path,'surface_coor.txt',num2str(me));
    if exist(fname, 'file')
        me
        loaddata = 1;
        a=load(fname);
        [n,m]=size(a);
        xtmp=a(1:n,1);
        ytmp=a(1:n,2);
        for i = 1:n
            if abs(xtmp(i)/1e3) < x1 && abs(ytmp(i)/1e3) < y1
                loaddata = 1;
            end
        end
        if loaddata == 1
            xcoor(1:n)=a(1:n,1);
            ycoor(1:n)=a(1:n,2);
            fname1=strcat(path,'gm',num2str(me));
            %b=load(fname1);
            fileID = fopen(fname1);
            C = fread(fileID, 'double');
            n1 = size(C,1);
            nt = n1/n/3;
            
            for i = 1:n % loop over stations
                %if abs(xcoor(initial_tag+i-1)/1e3) < 1 && (abs(ycoor(initial_tag+i-1)/1e3) < 1 ||abs(ycoor(initial_tag+i-1)/1e3-40) < 1)
                 if abs(xcoor(i)/1e3) < x1 && abs(ycoor(i)/1e3) < y1
                    ntag = ntag + 1;
                    coor(ntag,1) = xcoor(i);
                    coor(ntag,2) = ycoor(i);
                    for j = 1:nt % loop over time steps

                        timeseries_x(i,j) = C((j-1)*n*3 + (i-1)*3 + 1);
                        timeseries_y(i,j) = C((j-1)*n*3 + (i-1)*3 + 2); 
                    end
                    acc_x = vel_to_acc(timeseries_x(i,:)',dt);
                    acc_y = vel_to_acc(timeseries_y(i,:)',dt);
                    
                    sa_rotd50 = ResSpecRotD50(acc_x',acc_y',dt,T,gamma,NintMax);
                    res_sa(ntag,:) = sa_rotd50; 
                end
            end
            res_sa = res_sa/g;
            fname2=strcat(path,'rsa',num2str(me),'.mat');
            save(fname2,'res_sa','coor');
            delete b C xcoor ycoor;
        end
        delete a

    end
end

% xx = x0:dx:x1; yy = y0:dx:y1;
% xx = xx*1e3; yy = yy*1e3;
% [x2,y2] = meshgrid(xx,yy);
% F = scatteredInterpolant(coor(:,1),coor(:,2), res_sa(:,1),'linear','nearest');
% sa_map = F(x2,y2);
% %sa_map = griddata(coor(:,1),coor(:,2),res_sa(:,1),x2,y2);
% h=figure(1);
% contourf(x2,y2,sa_map); colorbar;
% %set(h,'LineColor','none');
% 
% [x2,y2] = meshgrid(xx,yy);
% F = scatteredInterpolant(coor(:,1),coor(:,2), res_sa(:,2),'linear','nearest');
% sa_map = F(x2,y2);
% %sa_map = griddata(coor(:,1),coor(:,2),res_sa(:,1),x2,y2);
% h=figure(2);
% contourf(x2,y2,sa_map); colorbar;
% %set(h,'LineColor','none');