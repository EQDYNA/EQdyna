%
% Response Spectrum Function
%
%    S. Day  29dec2011
%
% Uses time-domain ODE solver
%
% Assumes equally-spaced time points
%
% Minimum period must not be less than the time step.
%
     function Sa=ResSpecTimeDomVectorized(acc,dtt,T,gamma,NintMax)
     
%   Function returns Sa(T), vector of same length as T
%
% ...........................................................
% .          FUNCTION ARGUMENTS                             .
% .                                                         .
% .     acc  = ground acceleration time series (row vector) .
% .     dtt  = time increment                               .
% .      T   = vector of periods                            .
% .  gamma   = critical damping fraction                    .
% .  NintMax = maximum interpolation factor to be permitted .
% ...........................................................
%
% Check  minimum period and set interpolation factor Nint
%
NintMaxLimit=100;
if (NintMax>NintMaxLimit);
    ['Nintmax cannot exceed' num2str(NintMaxLimit)]
    ['To override, reset "NintMaxLimit" in response spectrum function'];
end
if (min(T)<4*dtt/NintMax);
    [ 'MINIMUM PERIOD MUST BE GREATER THAN' num2str(4*dtt/NintMax) ]
    return
end
Nint=ceil(4*dtt/min(T))
%
%
% Eigenvalues
%
l1=-gamma+i*sqrt(1-gamma^2);   
l2=-gamma-i*sqrt(1-gamma^2); 
%
% Right (C for column) and Left (R for row) Eigenvectors
%
C=[l1,l2;1,1];   
R=(1/2i/sqrt(1-gamma^2))*[1,-l2;-1,l1];
%
dt=dtt/Nint;
om0=2*pi./T;
E11=exp(l1*dt*om0);
E22=exp(l2*dt*om0);
Eint11=(1/l1./om0).*(1-E11);
Eint22=(1/l2./om0).*(1-E22);
%
%  Initialize
%
  sample=ceil(0.5*max(T)/dtt); 
  Uold=zeros(2,length(om0));
  umax=zeros(1,length(om0));
  %aEx=0.5*([acc(2:length(acc)) zeros(1,sample)]+...
      %[acc(1:length(acc)-1) zeros(1,sample)]);
      aEx=[acc(1:length(acc)) zeros(1,sample)];
      for j=1:length(acc)+sample-1
      for k=1:Nint
          aExInterp(Nint*(j-1)+k)=aEx(j)*(Nint-k+1)/Nint+aEx(j+1)*(k-1)/Nint;
      end
      end
      aEx=aExInterp;
      Ntot=length(aExInterp)
% 
% Form propagator C*E*R and Integrate the ODE
%
for j=1:Ntot
  Unew(1,:)=R(1,1)*(aEx(j)*Eint11+Uold(1,:).*E11)+R(1,2)*Uold(2,:).*E11;
  Unew(2,:)=R(2,1)*(aEx(j)*Eint22+Uold(1,:).*E22)+R(2,2)*Uold(2,:).*E22;
  Unew=C*Unew;
  Uold=Unew;
  umax=max(abs([Unew(2,:);umax]));
%
% Save Pseudo Acceleration
%
end
Sa=om0.*umax;