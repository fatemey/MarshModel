function [a, b] = BoxModel_parameters(t, y)
% plot_BoxModel_vars: plots different parameters
%
% Input
%           t : vector of time data in yr from BoxModel.
%           y : matrix of data from BoxModel. 1st vertical vector: tidal flat width (m),
%           2nd vector: tidal flat depth (m), 3rd vector: marsh depth (m),
%           and 4th vector: C_r (g/m3).
%           x_par : char vector of the x axis parameter of interest. It can
%                       be 't', 'd', 'h' or 'fetch'
%           y_par : char vector of the y axis parameter of interest. It can
%                       be 'tau', 'bed_e', 'wave', 'mar_e' (for margin erosion rate),
%                       'margin' (for margin accretion or erosion rate), 'ocean' (for ocean inout),
%                       'ratio' (for ocean to river input ratio). if margin
%                       is positive then erosion has happened, otherwise,
%                       accretion has happened.
%
% Output
%           a : vector of x-axis data of interset for plotting purposes
%           b : vector of y-axis data of interset for plotting purposes
%
% Last Update: 8/2/2017
%
%--------------------------------------------------------------------------------------------------

%-------------- Sediment input constants
C_o = 35 *10^-3;    % ocean concertation (kg/m3)
C_f = 15 *10^-3;   % river concentration (kg/m3)
Q_f = 20;          % river water discharge (m3/s)

%-------------- Erosion constants
k_0 = 1 *10^-3; % roughness (m)
v_w = 6;        % reference wind speed (m/s)

% -------------- Accretion constants
k_a = 2;        % margin accretion coefficient

%-------------- Vegetation properties
B_max = 1;      % maximum biomass density (kg/m2)
k_B = 2*10^-3 /365/24/60/60;    % vegetation characteristics (m3/s/kg)

%-------------- Basin properties
b_fm = 10 *10^3; % total basin width (both sides of the channel) (m)
L_E = 10 *10^3; % basin length (m)

%-------------- Tide Characteristics
T_T = 12 *60*60;   % tidal period (s) (= 12 hours)
H = 1.4 /2;          % tidal amplitude (range/2) (m)

%-------------- Sediment properties
rho_s = 1000;   % sediment bulk density (kg/m3)
omega_s = 0.5 *10^-3;   % settling velocity (m/s)

%-------------- Model constants
gamma = 9800;   % water specific weight (N/m3 or kg/m2/s2)

%-------------- Model assumptions
Q_f = Q_f/2;      % consider half of the discharge only for one side of the tidal platform (the same will be automatically considered below for Q_T)
% b_fm = b_fm/2;  % consider half of the basin only for one side of the tidal platform

%-------------- Model parameters
b_f = y(:,1);
d_f = y(:,2);
d_m = y(:,3);
C_r = y(:,4);
b_m = b_fm-b_f; % marsh width
Q_T = (d_f.*b_f+d_m.*b_m)*L_E/T_T-Q_f;
M_ocean2river = Q_T*C_o/Q_f/C_f; % Ocean Sediment Input (kg/s)
V_ocean2river = Q_T/Q_f; % Ocean Sediment Input (kg/s)
chi = b_f*2; % fetch (m)
h = (d_f+max(0,d_f-2*H))/2; % characteristic depth (m)
Mar_accr = k_a*omega_s*C_r/rho_s *365*24*60*60;
Bed_accr = C_r.*d_f/T_T/rho_s *365*24*60*60*1000;
z = H-d_m;       % elevation of marsh platform
r = -0.5*z/H+1;     % reproduction rate
m = 0.5*z/H;        % mortality rate
B = B_max*(1-m./r);  % steady state solution for biomass (kg/m2)
SOM = k_B*B *365*24*60*60*1000;  % organic matter production rate
SOM(z<0) = 0;

tau = zeros(size(h));
W = zeros(size(h)); % Wave Power Density (kg.m/s3)
for i=1:length(h)
    [ H_w, T_w] = WaveProps ( h(i), v_w, chi(i));   % compute significant height and peak period
    [ tau(i), k_w ] = ShearStress ( h(i), k_0, H_w, T_w);
    c_g = pi/k_w/T_w*(1+2*k_w*h(i)/sinh(2*k_w*h(i))); % wave group velocity (general form)
    W(i) = gamma*c_g*H_w^2/16; % wave power density (kg.m/s3)
end

%-------------- Plot the results
figure
clf
n = 3; m = 7;

subplot(n,m,1)
plot(t,b_f,'linewidth',2)
xlabel('Year')
ylabel('Tidal Flat Width (m)')
title('Solutions')

subplot(n,m,m+1)
plot(t,d_f,'linewidth',2)
xlabel('Year')
ylabel('Tidal Flat Depth (m)')

subplot(n,m,2*m+1)
plot(t,d_m,'linewidth',2)
xlabel('Year')
ylabel('Marsh Depth (m)')

subplot(n,m,2+2*m)
plot(t,tau,'linewidth',2)
xlabel('Year')
ylabel('\tau (PA)')

subplot(n,m,2+m)
plot(h*2,tau,'linewidth',2)
xlabel('Reference Depth (x2, m)')
ylabel('\tau (PA)')

subplot(n,m,2)
plot(chi/2,tau,'linewidth',2)
xlabel('Fetch (/2, m)')
ylabel('\tau (PA)')
title('Shear Stress')

subplot(n,m,3+2*m)
plot(t,W,'linewidth',2)
xlabel('Year')
ylabel('W (kg.m/s^3)')

subplot(n,m,3+m)
plot(h*2,W,'linewidth',2)
xlabel('Reference Depth (x2, m)')
ylabel('W (kg.m/s^3)')

subplot(n,m,3)
plot(chi/2,W,'linewidth',2)
xlabel('Fetch (/2, m)')
ylabel('W (kg.m/s^3)')
title('Wave Power Density')

subplot(n,m,4+2*m)
plot(t,Mar_accr,'linewidth',2)
xlabel('Year')
ylabel('Accretion Rate (m/yr)')

subplot(n,m,4+m)
plot(h*2,Mar_accr,'linewidth',2)
xlabel('Reference Depth (x2, m)')
ylabel('Accretion Rate (m/yr)')

subplot(n,m,4)
plot(chi/2,Mar_accr,'linewidth',2)
xlabel('Fetch (/2, m)')
ylabel('Accretion Rate (m/yr)')
title('Margin Accretion Rate')

subplot(n,m,5+2*m)
plot(t,M_ocean2river,'linewidth',2)
xlabel('Year')
ylabel('Qo.Co/Qf.Cf')

subplot(n,m,5+m)
plot(h*2,M_ocean2river,'linewidth',2)
xlabel('Reference Depth (x2, m)')
ylabel('Qo.Co/Qf.Cf')

subplot(n,m,5)
plot(chi/2,M_ocean2river,'linewidth',2)
xlabel('Fetch (/2, m)')
ylabel('Qo.Co/Qf.Cf')
title('Sediment Input Ratio')

subplot(n,m,6+2*m)
plot(t,V_ocean2river,'linewidth',2)
xlabel('Year')
ylabel('Qo/Qf')

subplot(n,m,6+m)
plot(h*2,V_ocean2river,'linewidth',2)
xlabel('Reference Depth (x2, m)')
ylabel('Qo/Qf')

subplot(n,m,6)
plot(chi/2,V_ocean2river,'linewidth',2)
xlabel('Fetch (/2, m)')
ylabel('Qo/Qf')
title('Tidal Prism Ratio')

subplot(n,m,7+2*m)
plot(t,SOM,'linewidth',2)
xlabel('Year')
ylabel('SOM Production (mm/yr)')

subplot(n,m,7+m)
plot(d_m,SOM,'linewidth',2)
xlabel('Marsh Depth (m)')
ylabel('SOM Production (mm/yr)')
title('SOM Production Rate')

subplot(n,m,7)
plot(t,Bed_accr,'linewidth',2)
xlabel('Year')
ylabel('Bed Accreation (mm/yr)')
title('TF Bed Accreation Rate')

% h_fig=gcf;
% set(h_fig,'PaperOrientation','portrait')
% set(h_fig,'PaperPosition', [0 0 14 6]) % [... ... max_width=7.5 max_height=9]
% tit='Co_95_parameters_ssshal';
% print(tit,'-dtiff','-r400')
% % movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\15 - Model Parameters relationships')
% movefile([tit,'.tif'],'/Users/Callisto/Files/Work/Marsh Model/Results/20 - Unstable Equlibrium Results through Optimization of Steady State')
% close all