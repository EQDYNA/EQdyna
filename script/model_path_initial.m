function [path,x0,z0,np,dx, dt] = model_path_initial(TPV)
%MODEL_PATH_INITIAL Summary of this function goes here
%   Detailed explanation goes here
dt = 0;
if TPV == 104
    np=40;
    dx=100;
    x0 = 18e3; z0 = 18e3;
    path = '../result/TPV104/';
elseif TPV == 105
    np=40;
    dx=100;
    x0 = 22e3; z0 = 22e3;
    path = '../result/TPV1053D/';na=2*x0/dx+1;
elseif TPV == 1041
    np=40;
    dx=80;
    x0 = 20.44e3; z0 = 20.4e3;
    path = '../result/TPV104_rough_fault/';na=2*x0/dx+1;    
elseif TPV == 1042
    np=144;
    dx=80;
    x0 = 20.44e3; z0 = 20.4e3;
    dt = 0.005;
    path = '../result/TPV104_rough_fault_kcut2belowNyquist/';na=2*x0/dx+1;   
elseif TPV == 1043
    np=1152;
    dx=80;
    x0 = 20.44e3; z0 = 20.4e3;
    dt = 0.005;
    path = '../result/P28_rough_kcut2_large/';na=2*x0/dx+1;       
elseif TPV == 2801
    np=2688;
    dx=50;
    x0 = 20e3; z0 = 20e3;
    path = '../result/TPV2801_50m/';na=2*x0/dx+1;     
elseif TPV == 2802
    np=3840;
    dx=25;
    x0 = 20e3; z0 = 20e3;
    path = '../result/TPV2801_25m_Grace/';na=2*x0/dx+1;      
elseif TPV == 28001
    np=960;
    dx=50;
    x0 = 25e3; z0 = 25e3;
    path = '../result/TPV2800_50m_a_Grace/';na=2*x0/dx+1; 
elseif TPV == 28021
    np=960;
    dx=50;
    dt = 0.004;
    x0 = 25e3; z0 = 25e3;
    path = '../result/TPV2802_50m_a1_Grace/';na=2*x0/dx+1;  
elseif TPV == 28022
    np=960;
    dx=50;
    dt = 0.004;
    x0 = 25e3; z0 = 25e3;
    path = '../result/TPV2802_50m_a2_Grace/';na=2*x0/dx+1;      
elseif TPV == 28023
    np=960;
    dx=50;
    dt = 0.004;
    x0 = 25e3; z0 = 25e3;
    path = '../result/TPV2802_50m_a3_Grace/';na=2*x0/dx+1;   
elseif TPV == 28024
    np=960;
    dx=50;
    dt = 0.004;
    x0 = 25e3; z0 = 25e3;
    path = '../result/TPV2802_50m_a4_Grace/';na=2*x0/dx+1;     
elseif TPV == 28025
    np=960;
    dx=50;
    dt = 0.004;
    x0 = 25e3; z0 = 25e3;
    path = '../result/TPV2802_50m_a5_Grace/';na=2*x0/dx+1;       
elseif TPV == 28026
    np=960;
    dx=50;
    dt = 0.004;
    x0 = 25e3; z0 = 25e3;
    path = '../result/TPV2802_50m_a6_Grace/';na=2*x0/dx+1;           
elseif TPV == 28031
    np=56;
    dx=100;
    x0 = 25e3; z0 = 25e3;
    path = '../result/TPV2802_100m_Frontera/';na=2*x0/dx+1;      
end    

