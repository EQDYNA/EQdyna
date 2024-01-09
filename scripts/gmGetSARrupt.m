function [stat, R] = gmGetSARrupt(coor, sa, xLimit)
    % gmGetSARrupt.m calculates SA as a function of Rupture distance Rrupt, the
    % closest distance from a station to the fault. Rrupt ranges from 1 to 20 km 
    % with a bin size of 1 km.
    % Input: 
    %   -coor, coordiantes of the stations;
    %   -sa, PSA at stations;
    %   -xLimit, x coordinate of one fault end; 
    % Output:
    %   -stat, with a dimension of nr by 2, that stores the mean and std of
    %       PSA for differnet Rrupt.
    %   -R, with a dimension of nr, stores the Rrupt.
    % Dunyu Liu <dliu@ig.utexas.edu>, 07/22/2021.
    
    nStation = size(coor,1);
    coor     = coor/1e3; % convert to km
    R        = 1:1:19; % in km
    nr       = length(R);
    sa_count = zeros(nr,1);
    sa_rec   = sa_count;
    
    for i = 1: nStation
        if abs(coor(i,1))> xLimit % x coordinate > xLimit km
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
    
    for i = 1: nr
        tmp = sa_rec(i,1:sa_count(i,1));
        stat(i,1) = mean(tmp);
        stat(i,2) = std(tmp);
        clear tmp;
    end
end

