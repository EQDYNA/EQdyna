function [stat,R] = sa_Rrupt(coor, sa)
%SA_RRUPT is the function to calculate the sa as a function of Rrupt, the
%closest distance to the fault from 1-30 km with a bin size of 1 km.
% Input: coor, coordiantes of the stations;
%   sa, psa at stations;
% Output:
%   stat, with a dimension of nr by 2, that stores the mean and std of sa
%   for differnet Rrupt.
%   R, with a dimension of nr, stores the Rrupt.
% First created on 07/22/2021.
% Author: Dunyu Liu, (dliu@ig.utexas.edu).

n = size(coor,1);
coor = coor/1e3;

R = 1:1:20;
nr = length(R);
sa_count = zeros(nr,1);
sa_rec = sa_count;
for i = 1:n
    if abs(coor(i,1))>25
        r(i) = (coor(i,1)^2 + coor(i,2)^2)^0.5;
    else
        r(i) = abs(coor(i,2));
    end
    
    intr = floor(r(i));
    if intr<nr+1 && intr>0
        sa_count(intr,1) = sa_count(intr,1) + 1;
        sa_rec(intr,sa_count(intr,1)) = sa(i);
    end
end

for i = 1:nr
    tmp = sa_rec(i,1:sa_count(i,1));
    stat(i,1) = mean(tmp);
    stat(i,2) = std(tmp);
    clear tmp;
end

end

