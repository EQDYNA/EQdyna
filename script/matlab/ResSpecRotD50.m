%
% Response Spectrum Function for RotD50
%
%    S.M. Day, 30Jul2012 (based on ResSpecTimeDomVectorized_3,
%                         written by  Day  29dec2011)
%
% Uses time-domain ODE solver
%
% Assumes equally-spaced time points
%
% Minimum period must not be less than the time step.
%
     function Sa=ResSpecRotD50(accx,accy,dtt,T,gamma,NintMax)
%
% ...........................................................
% .          FUNCTION ARGUMENTS                             .
% .                                                         .
% .     accx  = x comp ground acc time series (row vector)  .  
% .     accy  = y comp ground acc time series (row vector)  .
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
Nint=ceil(4*dtt/min(T));
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
%  
for component=1:2
%
%  Initialize
%
  sample=ceil(0.5*max(T)/dtt); 
  if (component==1)
      acc=accx;
  else
      acc=accy;
  end
  Unew=zeros(2,length(om0));
  Uold=zeros(2,length(om0));
  
      aEx=[acc(1:length(acc)) zeros(1,sample)];
      for j=1:length(acc)+sample-1
      for k=1:Nint
          aExInterp(Nint*(j-1)+k)=aEx(j)*(Nint-k+1)/Nint+aEx(j+1)*(k-1)/Nint;
      end
      end
      aEx=aExInterp;
      Ntot=length(aExInterp);
   if(component==1)
     ux=zeros(Ntot,length(T));
   else
     uy=zeros(Ntot,length(T));
   end
% 
% Form propagator C*E*R and Integrate the ODE
%
  for j=1:Ntot
    Unew(1,:)=R(1,1)*(aEx(j)*Eint11+Uold(1,:).*E11)+R(1,2)*Uold(2,:).*E11;
    Unew(2,:)=R(2,1)*(aEx(j)*Eint22+Uold(1,:).*E22)+R(2,2)*Uold(2,:).*E22;
    Unew=C*Unew;
    Uold=Unew;
    if(component==1)
      ux(j,:)=Unew(2,:);
    else
      uy(j,:)=Unew(2,:);
    end
  end
  
end

%dtheta=pi/180;
%cth=cos(dtheta);
%sth=sin(dtheta);

GeoMax=zeros(90,length(T));
for j=1:90
 cth=cos((j-1)*pi/180);
 sth=sin((j-1)*pi/180);
osc1=cth*ux+sth*uy;
osc2=-sth*ux+cth*uy;
GeoMax(j,:)=sqrt(max(abs(osc2)).*max(abs(osc1)));
end
Sa=om0.*median(GeoMax);