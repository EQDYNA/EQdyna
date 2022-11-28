function [shear_mod] = calc_shear_mod(depth, TPV)
%CALC_SHEAR_MOD Summary of this function goes here
%   Created on 09/05/2021 by Dunyu Liu (dliu@ig.utexas.edu)
% The script is to calculate shear modulus given a velocity structure.

if depth<0 
    a = 'ERROR: Depth should be positive downward.'
end

if TPV >28000
    
depth = depth/1e3; % in km.
if depth < 0.03
    vs = 2.206*depth^0.272;
elseif depth<0.19
    vs = 3.542*depth^0.407;
elseif depth<4.0
    vs = 2.505*depth^0.199;
elseif depth<8.0
    vs = 2.927*depth^0.086;
else 
    vs = 2.927*8^0.086;
end
vp = max(1.4+1.14*vs, 1.68*vs);
rou = 2.4405 + 0.10271*vs;

vs = vs*1e3;
vp = vp*1e3;
rou = rou*1e3;
shear_mod = vs^2*rou;

else
    vs = 3464;
    rou = 2700;
    shear_mod = vs^2*rou;
end
end
