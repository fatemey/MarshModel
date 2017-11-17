function [a, b] = plot_BoxModel_pars(t,y)
% plot_BoxModel_pars: plots different parameters
%
% Input
%           t : vector of time data in yr from BoxModel.
%           y : matrix of data from BoxModel. 1st vertical vector: tidal flat width (m),
%           2nd vector: tidal flat depth (m), 3rd vector: marsh depth (m),
%           and 4th vector: C_r (g/m3).
%           c : concentration 
%
% Output
%           a : vector of x-axis data of interset for plotting purposes
%           b : vector of y-axis data of interset for plotting purposes
%
% Last Update: 11/17/2017
%
%--------------------------------------------------------------------------------------------------

%-------------- Sediment input constants
C_o = 20 *10^-3;    % ocean concertation (kg/m3)
C_f = 15 *10^-3;   % river concentration (kg/m3)
Q_f = 20;          % river water discharge (m3/s)

%-------------- Erosion constants
k_0 = 1 *10^-3; % roughness (m)
v_w = 6;        % reference wind speed (m/s)

% -------------- Accretion constants
k_a = 2;        % margin accretion coefficient

%-------------- Basin properties
b_fm = 5 *10^3; % total basin width (both sides of the channel) (m)
L_E = 5 *10^3; % basin length (m)

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
Vol = (d_f.*b_f+d_m.*b_m)*L_E; % availble volume in the system to be filled with water
Q_T = max(Vol/T_T-Q_f,0);
M_ocean2river = Q_T*C_o/Q_f/C_f; % Ocean Sediment Input (kg/s)
V_ocean2river = Q_T/Q_f; % Ocean Sediment Input (kg/s)
chi = b_f*2; % fetch (m)
h = (d_f+max(0,d_f-2*H))/2; % characteristic depth (m)
Mar_accr = k_a*omega_s*C_r/rho_s *365*24*60*60;
Bed_accr = C_r.*d_f/T_T/rho_s *365*24*60*60*1000;

tau = zeros(size(h));
W = zeros(size(h)); % Wave Power Density (kg.m/s3)
for i=1:length(h)
    [ H_w, T_w] = WaveProps ( h(i), v_w, chi(i));   % compute significant height and peak period
    [ tau(i), k_w ] = ShearStress ( h(i), k_0, H_w, T_w);
    c_g = pi/k_w/T_w*(1+2*k_w*h(i)/sinh(2*k_w*h(i))); % wave group velocity (general form)
    W(i) = gamma*c_g*H_w^2/16; % wave power density (kg.m/s3)
end

%-------------- Plot the results
% figure
clf
n = 3; m = 4;

subplot(n,m,1)
plot(t,b_f,'linewidth',2)
xlabel('Year')
ylabel('Tidal Flat Width (m)')
title('TF Width')

subplot(n,m,2)
plot(t,d_f,'linewidth',2)
xlabel('Year')
ylabel('Tidal Flat Depth (m)')
title('TF Depth')

subplot(n,m,3)
plot(t,d_m,'linewidth',2)
xlabel('Year')
ylabel('Marsh Depth (m)')
title('M Depth')

subplot(n,m,4)
plot(t,C_r*1000,'linewidth',2)
xlabel('Year')
ylabel('Concentration (mg/l)')
title('Concentration')

subplot(n,m,5)
plot(t,tau,'linewidth',2)
xlabel('Year')
ylabel('\tau (PA)')
title('Shear Stress')

subplot(n,m,6)
plot(t,W,'linewidth',2)
xlabel('Year')
ylabel('W (kg.m/s^3)')
title('Wave Power Density')

subplot(n,m,7)
plot(t,M_ocean2river,'linewidth',2)
xlabel('Year')
ylabel('Qo.Co/Qr.Cr')
title('Sediment Input Ratio')

subplot(n,m,8)
plot(t,V_ocean2river,'linewidth',2)
xlabel('Year')
ylabel('Qo/Qr')
title('Tidal Prism Ratio')

subplot(n,m,9)
plot(t,Bed_accr,'linewidth',2)
xlabel('Year')
ylabel('Accreation Rate (mm/yr)')
title('TF Bed Accreation Rate')

subplot(n,m,10)
plot(t,Mar_accr,'linewidth',2)
xlabel('Year')
ylabel('Accretion Rate (m/yr)')
title('Margin Accretion Rate')

% h_fig=gcf;
% set(h_fig,'PaperOrientation','portrait')
% set(h_fig,'PaperPosition', [0 0 14 6]) % [... ... max_width=7.5 max_height=9]
% tit='Co_95_parameters_ssshal';
% print(tit,'-dtiff','-r400')
% % movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\15 - Model Parameters relationships')
% movefile([tit,'.tif'],'/Users/Callisto/Files/Work/Marsh Model/Results/20 - Unstable Equlibrium Results through Optimization of Steady State')
% close all