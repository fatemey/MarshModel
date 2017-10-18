function BoxModel_SS_2var()
% BoxModel_SS_2var : plots the model equations (3 of them) based on 2
% variables (marsh and tidal flat depths)
%
%--------------------------------------------------------------------------
format compact
format longG
clf

b_f = 1;
%-------------- Sediment input constants
C_o = 60 *10^-3;    % ocean concertation (kg/m3)
C_f = 15 *10^-3;    % river concentration (kg/m3)
Q_f = 20;         % river water discharge (m3/s)

%-------------- Erosion constants
k_0 = 1 *10^-3; % roughness (m)
tau_c = 0.3;  % critical shear stress (Pa)
E_0 = 10^-4;    % bed erosion coefficient (kg/m2/s)
k_e = 0.16 /365/24/60/60;  % margin erodibility coefficient (m2/s/W)
v_w = 6;        % reference wind speed (m/s)

% -------------- Accretion constants
k_a = 2;        % margin accretion coefficient

%-------------- Vegetation properties
B_max = 1;      % maximum biomass density (kg/m2)
k_B = 2*10^-3 /365/24/60/60;    % vegetation characteristics (m3/s/kg)

%-------------- Basin properties
b_fm = 10 *10^3; % total basin width (both sides of the channel) (m)
L_E = 10 *10^3; % basin length (m)
R = 2 *10^-3/365/24/60/60;   % sea level rise (m/s)
b_r = 0; % river width (m)

%-------------- Tide Characteristics
T_T = 12 *60*60;   % tidal period (s) (= 12 hours)
H = 1.4/2;          % tidal amplitude (range/2) (m)

%-------------- Sediment properties
rho_s = 1000;   % sediment bulk density (kg/m3)
omega_s = 0.5 *10^-3;   % settling velocity (m/s)

%-------------- Model constants
gamma = 9800;   % water specific weight (N/m3 or kg/m2/s2)

%-------------- Model assumptions
Q_f = Q_f/2;    % consider half of the discharge only for one side of the tidal platform (the same will be automatically considered below for Q_T)


[X,Y] = meshgrid(0.71:.01:1.4,.01:.01:.7);
[Z1,Z2,Z3]=Equations(X,Y);

subplot(2,2,1)
surf(X,Y,Z1/1000);
xlabel('Tidal Flat Depth (m)');ylabel('Marsh Depth (m)');zlabel('Tidal Flat Width Rate (m/yr)')
set(findobj('type','axes'),'fontsize',15)

subplot(2,2,2)
surf(X,Y,Z2/1000);
xlabel('Tidal Flat Depth (m)');ylabel('Marsh Depth (m)');zlabel('Tidal Flat Depth Rate (m/yr)')
set(findobj('type','axes'),'fontsize',15)

subplot(2,2,3)
surf(X,Y,Z3/1000);
xlabel('Tidal Flat Depth (m)');ylabel('Marsh Depth (m)');zlabel('Basin Elevation Rate (m/yr)')
set(findobj('type','axes'),'fontsize',15)

    function [F1,F2,F3]=Equations(d_f,d_m)
        %-------------- Imposing a condition for tidal flat conversion to marsh in case of presence of new vegetation when tidal flat is above MSL
        flag_f2m = 0;     % showing that tidal flat is below MSL
        if d_f <= H
            flag_f2m = 1; % showing that tidal flat is above MSL
        end
        
        %-------------- Model assumptions
        b_m = b_fm-b_f; % marsh width
        if flag_f2m == 0
            chi = 2*b_f+b_r;    % fetch
        elseif flag_f2m == 1
            chi = b_r;
        end
        
        %-------------- Compute organice matter production (m/s)
        z = H-d_m;       % elevation of marsh platform
        if z >= 0           % condition for presence of vegetation when marsh is above MSL
            r = -0.5*z/H+1;     % reproduction rate
            m = 0.5*z/H;         % mortality rate
            B = B_max*(1-m./r);  % steady state solution for biomass (kg/m2)
            O = k_B*B;            % organic matter production rate
        else
            O = 0;
        end
        
        %-------------- Compute the concentration (kg/m^3)
        C_r = (R-O)*T_T*rho_s./d_m;
        
        %---------------------------------- Define the equations ----------------------------------
        
        %-------------------------------- Tidal flat width equation --------------------------------
        
        %-------------- Compute margin erosion (m/s)
        h = (d_f+max(0,d_f-2*H))/2;     % reference water depth
        
        if  flag_f2m==1 || chi<=b_r || v_w==0 % condition for no bed and margin erosion in case of a filled mudflat or no wind
            tau = 0; % bed shear stress
            W = 0;   % wave power density
        else
            [ H_w, T_w ] = WaveProps_2 ( h, v_w, chi );   % compute significant height and peak period
            [ tau, k_w ] = ShearStress_2 ( h, k_0, H_w, T_w );     % compute bed shear stress
            c_g = pi./k_w./T_w.*(1+2*k_w.*h./sinh(2*k_w.*h)); % wave group velocity (general form)
            W = gamma*c_g*H_w.^2/16; % wave power density (kg.m/s3)
        end
        
        B_e = k_e*W;    % margin erosion
        
        %-------------- Compute margin accretion (m/s)
        if flag_f2m == 0
            B_a = k_a*omega_s*C_r/rho_s;  % margin accretion
        else
            B_a = 0;
        end
        
        %-------------- Describe the equation for b_f (m)
        F1 = abs(B_e - B_a) ; % (m/s)
        F1 = F1 *1000*60*60*24*365; % (mm/s)
        
        %--------------------------------- Tidal flat depth equation -------------------------------
        if flag_f2m == 0
            
            %-------------- Compute submerged time when the tidal flat is covered with water
            if d_f < 2*H
                t_s = 1/2-1/pi*asin((H-d_f)/H);
            else
                t_s = 1;
            end
            
            %-------------- Compute the rate of bed erosion (m/s)
            TF_erosion = max(0,t_s*E_0/rho_s*(tau-tau_c)/tau_c);
            
            %-------------- Compute the rate of sediment accretion (m/s)
            TF_accretion = min(t_s*C_r*omega_s/rho_s ,C_r*d_f/T_T/rho_s);
            
            %-------------- Compute the rate of organic matter production in tidal flat (m/s)
            SOM = 0;
            
        elseif flag_f2m == 1
            
            %-------------- Compute the rate of bed erosion (m/s)
            TF_erosion = 0;
            
            %-------------- Compute the rate of sediment accretion (m/s)
            TF_accretion = C_r*d_f/T_T/rho_s;
            
            %-------------- Compute the rate of organic matter production in the new marsh (m/s)
            z_new = H-d_f;       % elevation of marsh platform
            r_new = -0.5*z_new/H+1;     % reproduction rate
            m_new = 0.5*z_new/H;         % mortality rate
            B_new = B_max*(1-m_new/r_new);  % steady state solution for biomass (kg/m2)
            SOM = k_B*B_new;        % organic matter production rate
            
        end
        
        %-------------- Describe the equation for d_f (m)
        F2 = TF_erosion - TF_accretion - SOM + R; % (m/s)
        F2 = F2 *1000*60*60*24*365; % (mm/s)
        
        %-------------------------------- Mass conservation equation ------------------------------
        
        %-------------- Compute tidal flat bed erosion (kg/s)
        if flag_f2m == 0
            bed_erosion = max(0,t_s*E_0*(tau-tau_c)/tau_c*b_f*L_E);
        else
            bed_erosion = 0;
        end
        
        %-------------- Compute external sediment input (kg/s)
        Q_T = (d_f*b_f+d_m*b_m)*L_E/T_T-Q_f;
        ocean_in = Q_T*C_o;
        river_in = Q_f*C_f;
        
        %-------------- Compute deposition on tidal flat (kg/s)
        if flag_f2m == 0
            TF_deposition = min(t_s*C_r*b_f*omega_s*L_E, C_r*b_f*d_f*L_E/T_T);
        else
            TF_deposition = C_r*b_f*d_f*L_E/T_T;
        end
        
        %-------------- Compute deposition on marsh (kg/s)
        M_deposition = C_r*b_m*d_m*L_E/T_T;
        
        %-------------- Compute export sediment to the ocean (kg/s)
        if flag_f2m == 0
            export = C_r*b_f*min(d_f,2*H)*L_E/T_T;
        else
            export = 0;
        end
        
        %-------------- Describe the equation for C_r (kg/m3)
        F3 = bed_erosion + ocean_in + river_in - TF_deposition - M_deposition - export;   % (kg/s)
        F3 = F3/rho_s/L_E/b_fm;   % (m/s)
        F3 = F3 *1000*60*60*24*365; %(mm/s)
        
    end

    function [ H_w, T_w ] = WaveProps_2 ( h, v_w, chi)
        % WaveProps_2 : Computing significant wave height and peak wave period
        %
        % Input
        %       h   : reference water depth
        %       v_w : wind velocity
        %       chi : fetch
        %
        % Output
        %       H_w : significant wave height
        %       T_w : peak wave period
        %
        %--------------------------------------------------------------------------
        g = 9.81;       % gravitational acceleration (m/s2)
        
        %-------------- Parameters
        A1 = .493*(g*h/v_w^2).^.75;
        B1 = 3.13*10^-3*(g*chi/v_w^2).^.57;
        A2 = .331*(g*h/v_w^2).^1.01;
        B2 = 5.215*10^-4*(g*chi/v_w^2).^.73;
        
        %-------------- Outpouts
        H_w = .2413*(tanh(A1).*tanh(B1./tanh(A1))).^.87*v_w^2/g;   % significant wave height
        T_w = 7.518*(tanh(A2).*tanh(B2./tanh(A2))).^.37*v_w/g;     % peak wave period
        
    end

    function [ tau, k_w ] = ShearStress_2 ( h, k_0, H_w, T_w )
        % WaveProps : Computing bed shear stress by waves
        %
        % Input
        %       h   : reference water depth
        %       k_0 : roughness
        %       H_w : significant wave height
        %       T_w : peak wave period
        %
        % Output
        %       tau : shear stress
        %
        %--------------------------------------------------------------------------
        
        %-------------- Constants
        rho_w = 1000;   % water density (kg/m3)
        g = 9.81;       % gravitational acceleratxion (m/s2)
        
        %-------------- Parameters
        k_w = 2*pi./T_w./sqrt(g*h);               % wave number (for k_w*h<0.1pi)
        f_w = 0.4*(H_w/k_0./sinh(k_w*h)).^-.75;   % friction factor
        u_w = pi*H_w./T_w./sinh(k_w*h);           % horizontal orbital velocity
        
        %-------------- Outpouts
        tau = rho_w.*f_w.*u_w.^2/2;    % wave shear stress
        
    end


end

