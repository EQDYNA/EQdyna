clear all;close all;
np=100;
dx=100;
na=36e3/dx+1;
ma=18e3/dx+1;
for me=0:np-1
    fname=strcat('frt.txt',num2str(me));
    me
    if exist(fname, 'file')
        a=load(fname);
        [n,m]=size(a);
        for i=1:n
            ii=(a(i,1)--18e3)/dx+1;
            jj=(a(i,3)--18e3)/dx+1;
            rupt((ma+1-jj-1)*na+ii,1)=a(i,1);
            rupt((ma+1-jj-1)*na+ii,2)=-a(i,3);
            rupt((ma+1-jj-1)*na+ii,3)=a(i,4);
            t(ma+1-jj,ii)=a(i,4); 
            if a(i,4)>500
                t(ma+1-jj,ii) = -1;
            end
        end
        delete a;
    end
end

fileID = fopen('cplot.txt', 'w');
fprintf(fileID,'j k t\n#\n');
for i=1:na*ma
	fprintf(fileID,'%12.8f %12.8f %12.8f\n', rupt(i,1),rupt(i,2),rupt(i,3));
end
fclose(fileID);
x=-18e3:dx:18e3;
z=0:dx:18e3;
[xx,zz]=meshgrid(x,z);
line = 0:0.5:20;
figure(1)
contour(xx,zz,t,line);