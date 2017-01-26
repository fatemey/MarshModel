function MarshTidalFlat_BoxModel_IsotropicWind
% MarshTidalFlat_BoxModel_SA : Modeling 0d marsh and tidal flat time
% evolution using Matlab ode15s function for sensitivity analysis
%
% Purpose: Determining marsh and tidal flat depths and widths changes with
%          rising sea level. This model is solved using 4 equations and 4
%          unknowns.
%
%--------------------------------------------------------------------------------------------------
format compact
format longG
clear
clf
tic
kk=0;
C_o_ = [10*10^-3 20*10^-3 100*10^-3]; leg = {'10 g/m^3','20 g/m^3','100 g/m^3'}; tit = 'Ocean Concentration';
% C_f_ = [10*10^-3 20*10^-3 100*10^-3]; leg = {'0 g/m^3','15 g/m^3','100 g/m^3'}; tit = 'River Concentration';
% Q_f_ = [0 20 1000]; leg = {'0 m^3/s','20 m^3/s','1000 m^3/s'}; tit = 'River Discharge';
% tau_c_ = [.1 .3 .4 ]; leg = {'0.1 PA','0.3 PA','0.4 PA'}; tit = 'Critical Shear Stress';
% E_0_ = [10^-5 10^-4 10^-3]; leg = {'10^{-5} kg/m^2/s','10^{-4} kg/m^2/s','10^{-3} kg/m^2/s'}; tit = 'Bed Erosion Coefficient';
% k_e_ = [0.1/365/24/60/60 0.16/365/24/60/60 0.2/365/24/60/60]; leg = {'0.1 m^2/yr/W','0.16 m^2/yr/W','0.2 m^2/yr/W'}; tit = 'Margin Erodibility Coefficient';
% v_w_ = [0 6 10]; leg = {'0 m/s','6 m/s','10 m/s'}; tit = 'Wind Velocity';
% k_a_ = [1 2 4]; leg = {'1','2','4'}; tit ='Margin Accretion Coefficient';
% k_B_ = [1*10^-3/365/24/60/60 2*10^-3/365/24/60/60 4*10^-3/365/24/60/60]; leg = {'1 m3/yr/g','2 m3/yr/g','4 m3/yr/g'}; tit = 'Biomass Coefficient';
% R_ = [0 2*10^-3/365/24/60/60 10*10^-3/365/24/60/60]; leg = {'0 mm/yr','2 mm/yr','10 mm/yr'}; tit = 'Rate of Sea Level Rise';
% b_fm_ = [10^3 5*10^3 10*10^3];  leg = {'0.5 km','2.5 km','5 km'}; tit = 'Basin Width';
% L_E_ = [.15*10^3 15*10^3 150*10^3]; leg = {'0.15 km','15 km','150 km'}; tit = 'Estuary Length';
% b_f_0 = [1*10^3  5*10^3  9*10^3]/4; leg = {'0.25 km (0.1 L)','1.25 km (0.5 L)','2.25 km (0.9 L)'}; tit = 'Initial Tidal Flat Width';
% H_ = [1 1.4 2.4]/2; leg = {'1 m','1.4 m','2.4 m'}; tit = 'Tidal Range';

for i = 1 : length(leg)
    
%-------------- Set the time span
tyr = 300;  % solve for time tyr (years)
ts = tyr *365*24*60*60; % tyr in (s)
dt = 24*60*60; % time step in (s)
tspan = 0:dt:ts;

%-------------- Sediment input constants
% C_o = 20 *10^-3;    % ocean concertation (kg/m3)
C_o = C_o_(i);

C_f = 15 *10^-3;    % river concentration (kg/m3)
% C_f = C_f_(i);

Q_f = 20;         % river water discharge (m3/s)
% Q_f = Q_f_(i);

%-------------- Erosion constants
k_0 = 1 *10^-3; % roughness (m)
% k_0 = k_0_(i);

tau_c = 0.3;  % critical shear stress (Pa)
% tau_c = tau_c_(i);

E_0 = 10^-4;    % bed erosion coefficient (kg/m2/s)
% E_0 = E_0_(i);

k_e =  .16 /365/24/60/60;  % margin erodibility coefficient (m2/s/W)
% k_e = k_e_(i);

v_w = 6;        % reference wind speed (m/s)
% v_w = v_w_(i);
% v_w = v_w*(2*rand(length(tspan),1)-1);

% -------------- Accretion constants
k_a = 2;        % margin accretion coefficient
% k_a = k_a_(i);

%-------------- Vegetation properties
B_max = 1;      % maximum biomass density (kg/m2)

k_B = 2*10^-3 /365/24/60/60;    % vegetation characteristics (m3/s/kg)
% k_B = k_B_(i);

%-------------- Basin properties
b_fm = 5 *10^3; % total basin width (both sides of the channel) (m)
% b_fm = b_fm_(i);

L_E = 15 *10^3; % basin length (m)
% L_E = L_E_(i);

R = 2 *10^-3/365/24/60/60;   % sea level rise (m/s)
% R = R_(i);

%-------------- Tide Characteristics
H = 1.4 /2;          % tidal amplitude (range/2) (m)
% H = H_(i);          % tidal amplitude (range/2) (m)

T_T = 12 *60*60;   % tidal period (s) (= 12 hours)

%-------------- Sediment properties
rho_s = 1000;   % sediment bulk density (kg/m3)
omega_s = 0.5 *10^-3;   % settling velocity (m/s)

%-------------- Model constants
gamma = 9800;   % water specific weight (N/m3)
g = 9.81;       % gravitational acceleration (m/s2)

%-------------- Model assumptions
Q_f = Q_f/2;    % consider half of the discharge only for one side of the tidal platform (the same will be automatically considered below for Q_T)
b_fm = b_fm/2;  % consider half of the basin only for one side of the tidal platform

%-------------- Initial conditions, y0=[ b_f, d_f, d_m,u(=C_r*(b_f*d_f+b_m*d_m))]
y0(1) = b_fm/2;   % tidal flat width (m)
% y0(1) = b_f_0(i);

y0(2) = 1.0;         % tidal flat depth (m)
y0(3) = 0.4;         % marsh depth (m)
y0(4) = 0;

%-------------- Solve the system of differential equations
[t,y,te,ye,ie] = ode15s(@ode4marshtidalflat,tspan,y0); % or use ode15s/ode45/..23s
y(:,4) = y(:,4)./(y(:,1).*y(:,2)+y(:,3).*(b_fm-y(:,1))); % convert y(:,4) to C_r from the formula used before: y4=u (=C_r*(b_f*d_f+b_m*d_m)
ydata(:,:,i) = y; % save the solution for each SA value (3 values for the factor of interest)

end

%-------------- Plot Results
plot_results_BoxModel_SA(t,ydata,leg,tit)
% print(tit,'-dtiff','-r400')
% movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\8 - No wind condition')
% close all

timespent_min = toc/60

%======================= Nested Function =========================
    function dy = ode4marshtidalflat (t,y) %  y1=b_f, y2=d_f, y3=d_m, y4=u (=C_r*(b_f*d_f+b_m*d_m)

        v_w = v_w*(2*rand-1);
%         kk=kk+1
        %-------------- Setting width boundary limits
        local = y(1); % imposing a constraint for lower and upper limits of y(1)
        if (local < 0)
            local = 0 ;
        end
        if (local > b_fm)
            local = b_fm ;
        end
        
        %-------------- Model assumptions
        b_m = b_fm-local; % marsh width
        chi = 2*local;    % fetch
            
        %---------------------------------- Define the equations -----------------------------------
        
        %----------------------------- Marsh width changes equation ------------------------------
        
        %-------------- Compute margin erosion (m/s)
        if  chi<=0 || v_w<=0 % condition for no bed and margin erosion in case of a filled mudflat or no wind
            
            tau = 0;
            W = 0;
            
        else
            
            h = (y(2)+max(0,y(2)-2*H))/2;   % reference water depth
            [ H_w, T_w ] = WaveProps ( h, v_w, chi, g );   % compute significant height and peak period
            [ tau, k_w ] = ShearStress ( h, k_0, H_w, T_w, g );     % compute bed shear stress
            % c_g = sqrt(g*h);        % wave group velocity (shallow water)
            c_g = pi/k_w/T_w*(1+2*k_w*h/sinh(2*k_w*h)); % wave group velocity (general form)
            W = gamma*c_g*H_w^2/16; % wave power density (kg.m/s3)
            
        end
        
        B_e = k_e*W;    % margin erosion
        
        %-------------- Compute margin accretion (m/s)
        C_r = y(4)/(local*y(2)+b_m*y(3));
        B_a = k_a*omega_s*C_r/rho_s;    % margin accretion
        
        %-------------- Describe the equation for b_f (m)
        dy(1,1) = B_e - B_a;
        
        if (y(1) < 0 && dy(1,1) < 0) % imposing a constraint for lower and upper limits of y(1)
            dy(1,1) = 0 ;
        end
        if (y(1) > b_fm && dy(1,1)>0)
            dy(1,1) = 0 ;
        end

        %----------------------------- Tidal flat depth changes equation ---------------------------
        
        %-------------- Compute effective time when the tidal flat is submerged
        if y(2) < 2*H
            t_e = 1/2-1/pi*asin((H-y(2))/H);
        else
            t_e = 1;
        end
        
        %-------------- Compute the rate of bed erosion (m/s)
        TF_erosion = max(0,t_e*E_0/rho_s*(tau-tau_c)/tau_c);
        
        %-------------- Compute the rate of sediment accretion (m/s)
        TF_accretion = t_e*C_r*omega_s/rho_s;
        
        %-------------- Describe the equation for d_f (m)
        dy(2,1) = TF_erosion - TF_accretion + R;
        
        %----------------------------- Marsh depth changes equation ------------------------------
        
        %-------------- Compute the rate of sediment accretion (m/s)
        M_accretion = C_r*y(3)/T_T/rho_s;
        
        %-------------- Compute organice matter production (m/s)
        z = H-y(3);            % elevation of marsh platform
        r = -0.5*z/H+1;     % reproduction rate
        m = 0.5*z/H;         % mortality rate
        B = B_max*(1-m/r);  % steady state solution of B (kg/m2)
        O = k_B*B;             % organic matter production rate
        
        %-------------- Describe the equation for d_m (m)
        dy(3,1) = - M_accretion - O + R;
        
        %-------------------------------- Mass conservation equation ------------------------------
        
        %-------------- Compute tidal flat bed erosion (kg/s)
        bed_erosion = max(0,t_e*E_0*(tau-tau_c)/tau_c*local*L_E);
        
        %-------------- Compute marsh/tidal flat margin erosion/accretion (kg/s)
        margin = (B_e - B_a)*(y(2)-y(3))*L_E*rho_s;
        
        %-------------- Compute external sediment input (kg/s)
        Q_T = (y(2)*local+y(3)*b_m)*L_E/T_T-Q_f;
        ocean_in = Q_T*C_o;
        river_in = Q_f*C_f;
        
        %-------------- Compute deposition on tidal flat and marsh bed (kg/s)
        TF_deposition = t_e*C_r*local*omega_s*L_E;
        M_deposition = C_r*b_m*y(3)*L_E/T_T;
        
        %-------------- Compute export sediment to the ocean (kg/s)
        export = C_r*local*min(y(2),2*H)*L_E/T_T;
        
        %-------------- Describe the equation for C_r (kg/m3)
        var = bed_erosion + margin + ocean_in + river_in - TF_deposition - M_deposition - export;   % (kg/s)
        dy(4,1) =  var / L_E;
                 
    end

end