function [y, t] = BoxModel_SplittedResults (t1, t2, y0,v_w)
% BoxModel_SplittedResults : Modeling 0d marsh and tidal flat time
% evolution using Matlab ode15s function for sensitivity analysis. This
% function uses splitted time periods to be used in BoxModel_CombinedResults
%
% Purpose: Determining marsh and tidal flat depths and widths changes with
%          rising sea level. This model is solved using 4 equations and 4
%          unknowns.
%
%--------------------------------------------------------------------------------------------------

%-------------- Set the time span
t1s = t1 *365*24*60*60; % t1 in (s)
t2s = t2 *365*24*60*60; % t2 in (s)
dt = 12*60*60; % time step in (s)
tspan = t1s:dt:t2s;

%-------------- Sediment input constants
C_o = 20 *10^-3;    % ocean concertation (kg/m3)
C_f = 15 *10^-3;    % river concentration (kg/m3)
Q_f = 20;         % river water discharge (m3/s)

%-------------- Erosion constants
k_0 = 1 *10^-3; % roughness (m)
tau_c = 0.3;  % critical shear stress (Pa)
E_0 = 10^-4;    % bed erosion coefficient (kg/m2/s)
k_e =  0.16 /365/24/60/60;  % margin erodibility coefficient (m2/s/W)
% v_w = 6;        % reference wind speed (m/s)

% -------------- Accretion constants
k_a = 2;        % margin accretion coefficient

%-------------- Vegetation properties
B_max = 1;      % maximum biomass density (kg/m2)
k_B = 2*10^-3 /365/24/60/60;    % vegetation characteristics (m3/s/kg)

%-------------- Basin properties
b_fm = 5 *10^3; % total basin width (both sides of the channel) (m)
L_E = 15 *10^3; % basin length (m)
R = 2 *10^-3/365/24/60/60;   % sea level rise (m/s)

%-------------- Tide Characteristics
T_T = 12 *60*60;   % tidal period (s) (= 12 hours)
H = 1.4 /2;          % tidal amplitude (range/2) (m)

%-------------- Sediment properties
rho_s = 1000;   % sediment bulk density (kg/m3)
omega_s = 0.5 *10^-3;   % settling velocity (m/s)

%-------------- Model constants
gamma = 9800;   % water specific weight (N/m3)
g = 9.81;       % gravitational acceleration (m/s2)

%-------------- Model assumptions
Q_f = Q_f/2;    % consider half of the discharge only for one side of the tidal platform (the same will be automatically considered below for Q_T)
b_fm = b_fm/2;  % consider half of the basin only for one side of the tidal platform

%-------------- Solve the system of differential equations
[t, y] = ode15s(@ode4marshtidalflat,tspan,y0); % or use ode15s/ode45/..23s
y(:,4) = y(:,4)./(y(:,1).*y(:,2)+y(:,3).*(b_fm-y(:,1))); % convert y(:,4) to C_r from the formula used before: y4=u (=C_r*(b_f*d_f+b_m*d_m)



%======================= Nested Function =========================
    function dy = ode4marshtidalflat (t,y) %  y1=b_f, y2=d_f, y3=d_m, y4=u (=C_r*(b_f*d_f+b_m*d_m, why solving u instead of C_r? u is the variable on the left hand side of mass conservation equation.)
        
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
        h = (y(2)+max(0,y(2)-2*H))/2;   % reference water depth
        
        if  chi<=0 || v_w==0 % condition for no bed and margin erosion in case of a filled mudflat or no wind
            
            tau = 0;
            W = 0;
            
        else
            
            [ H_w, T_w ] = WaveProps ( h, v_w, chi, g );   % compute significant height and peak period
            [ tau, k_w ] = ShearStress ( h, k_0, H_w, T_w, g );     % compute bed shear stress
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
        TF_accretion = min(t_e*C_r*omega_s/rho_s ,t_e*C_r*h/T_T/rho_s);
        
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
        TF_deposition = min(t_e*C_r*local*omega_s*L_E, t_e*C_r*local*h*L_E/T_T);
        M_deposition = C_r*b_m*y(3)*L_E/T_T;
        
        %-------------- Compute export sediment to the ocean (kg/s)
        export = C_r*local*min(y(2),2*H)*L_E/T_T;
        
        %-------------- Describe the equation for C_r (kg/m3)
        var = bed_erosion + margin + ocean_in + river_in - TF_deposition - M_deposition - export;   % (kg/s)
        dy(4,1) =  var / L_E;
        
    end

end