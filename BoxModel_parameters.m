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
C_o = 20 *10^-3;    % ocean concertation (kg/m3)
C_f = 15 *10^-3;    % river concentration (kg/m3)
Q_f = 20;           % river water discharge (m3/s)

%-------------- Erosion constants
k_0 = 1 *10^-3; % roughness (m)
v_w = 6;        % reference wind speed (m/s)

%-------------- Basin properties
b_fm = 10 *10^3; % total basin width (both sides of the channel) (m)
L_E = 10 *10^3;  % basin length (m)

%-------------- Tide Characteristics
T_T = 12 *60*60;   % tidal period (s) (= 12 hours)
H = 1.4 /2;        % tidal amplitude (range/2) (m)

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
h =  (d_f+max(0,d_f-2*H))/2; % characteristic depth (m)

tau = zeros(size(h));
W = zeros(size(h)); % Wave Power Density (kg.m/s3)
for i=1:length(h)
    [ H_w, T_w] = WaveProps ( h(i), v_w, chi(i));   % compute significant height and peak period
    [ tau(i), k_w ] = ShearStress ( h(i), k_0, H_w, T_w);
    c_g = pi/k_w/T_w*(1+2*k_w*h(i)/sinh(2*k_w*h(i))); % wave group velocity (general form)
    W(i) = gamma*c_g*H_w^2/16; % wave power density (kg.m/s3)
end

%-------------- Plot the results
figure(2)
clf
subplot(4,4,1)
plot(t,b_f,'linewidth',2)
xlabel('Year')
ylabel('Tidal Flat Width (m)')

subplot(4,4,2)
plot(t,d_f,'linewidth',2)
xlabel('Year')
ylabel('Tidal Flat Depth (m)')

subplot(4,4,3)
plot(t,d_m,'linewidth',2)
xlabel('Year')
ylabel('Marsh Depth (m)')

subplot(4,4,4)
plot(t,C_r,'linewidth',2)
xlabel('Year')
ylabel('Concentration (kg/m3)')

subplot(4,4,5)
plot(t,tau,'linewidth',2)
xlabel('Year')
ylabel('Shear Stress (PA)')

subplot(4,4,9)
plot(h,tau,'linewidth',2)
xlabel('Reference Depth (m)')
ylabel('Shear Stress (PA)')

subplot(4,4,13)
plot(chi,tau,'linewidth',2)
xlabel('Fetch (m)')
ylabel('Shear Stress (PA)')

subplot(4,4,6)
plot(t,W,'linewidth',2)
xlabel('Year')
ylabel('Wave Power Density (kg.m/s^3)')

subplot(4,4,10)
plot(h,W,'linewidth',2)
xlabel('Reference Depth (m)')
ylabel('Wave Power Density (kg.m/s^3)')

subplot(4,4,14)
plot(chi,W,'linewidth',2)
xlabel('Fetch (m)')
ylabel('Wave Power Density (kg.m/s^3)')

subplot(4,4,7)
plot(t,M_ocean2river,'linewidth',2)
xlabel('Year')
ylabel('Ocean to River Sediment Input')

subplot(4,4,11)
plot(h,M_ocean2river,'linewidth',2)
xlabel('Reference Depth (m)')
ylabel('Ocean to River Sediment Input')

subplot(4,4,15)
plot(chi,M_ocean2river,'linewidth',2)
xlabel('Fetch (m)')
ylabel('Ocean to River Sediment Input')

subplot(4,4,8)
plot(t,V_ocean2river,'linewidth',2)
xlabel('Year')
ylabel('Ocean to River Tidal Prism')

subplot(4,4,12)
plot(h,V_ocean2river,'linewidth',2)
xlabel('Reference Depth (m)')
ylabel('Ocean to River Tidal Prism')

subplot(4,4,16)
plot(chi,V_ocean2river,'linewidth',2)
xlabel('Fetch (m)')
ylabel('Ocean to River Tidal Prism')

% h_fig=gcf;
% set(h_fig,'PaperOrientation','portrait')
% set(h_fig,'PaperPosition', [0 0 7.5 7.5]) % [... ... max_width=7.5 max_height=9]
% tit='bf0';
% print(tit,'-dtiff','-r400')
% movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\15 - Model Parameters relationships')
% close all