function [RSAstats, R, numSA, r, sa_rec] = gmGetSARrupt(coor, sa, xLimit, RBinSize, Rmin, Rmax)
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
    R0       = Rmin:RBinSize:floor(Rmax/RBinSize)*RBinSize; % in km
    nr       = length(R0);
    % bin center R
    R        = Rmin+RBinSize/2 : RBinSize: (floor(Rmax/RBinSize)-0.5)*RBinSize;
    numSA    = zeros(nr-1); 
    sa_rec   = zeros(nr-1,1000);
    logSA    = sa_rec;
    RSAstats = zeros(nr-1,20);
    % 1:mean of logSA; 2:meanSqLog; 3:std of logSA; 4: min; 5:25 percentile;   
    %MeanSqLogSA = LogSA;

    for i = 1: nStation
        if coor(i,1)> xLimit % x coordinate > xLimit km
            r(i) = ((coor(i,1)-xLimit)^2 + coor(i,2)^2)^0.5;
        elseif coor(i,1)<-xLimit
            r(i) = ((coor(i,1)--xLimit)^2 + coor(i,2)^2)^0.5;
        else
            r(i) = abs(coor(i,2));
        end
        
        for iR = 1: nr-1
            if r(i)>=R0(iR) && r(i)<R0(iR+1) && abs(sa(i))>0
                numSA(iR) = numSA(iR) + 1;
                sa_rec(iR,numSA(iR)) = sa(i);
            end
        end
    end
    
    % calculate mean and std of SA in log unit.
    for iR = 1: nr-1
        sumLogSA = 0;
        sumSqLogSA = 0;
        for iS = 1: numSA(iR)
            logSA(iR,iS) = log(sa_rec(iR,iS));
            sumLogSA = sumLogSA + logSA(iR,iS);
            sumSqLogSA = sumSqLogSA + (log(sa_rec(iR,iS)))^2;
        end
        meanLogSA = sumLogSA/numSA(iR);
        meanSqLogSA = sumSqLogSA/(numSA(iR)-1);
        stdLogSA = ( meanSqLogSA - meanLogSA^2 * numSA(iR)/(numSA(iR)-1) )^0.5;
        RSAstats(iR,1) = exp(meanLogSA); % mean
        RSAstats(iR,2) = stdLogSA;
        RSAstats(iR,3) = exp(min(logSA(iR,1:numSA(iR))));
        RSAstats(iR,4) = exp(prctile(logSA(iR,1:numSA(iR)), 25));
        RSAstats(iR,5) = exp(prctile(logSA(iR,1:numSA(iR)), 75));
        RSAstats(iR,6) = exp(max(logSA(iR,1:numSA(iR))));
    end
end

