function [acc] = vel_to_acc(vel,dt)
%VEL_TO_ACC Summary of this function goes here
%   Convert velocity time series to acceleration time series.
%   vel and acc are 1D columns.

n = size(vel,1);
acc(1,1) = 0;
for i = 2:n
    acc(i,1) = (vel(i,1) - vel(i-1,1))/dt; 
end

end

