function gmCheck(minVs, dx, n, dt)
% gmCheck calculates the maximum freqeuency properly simulated.
%
minWaveLength = dx*n;
minPeriod     = minWaveLength/minVs
maxFreq       = 1./minPeriod

minPeriod2    = 2*dt;
nyquistFreq   = 1./minPeriod2

end

