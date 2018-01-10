function BoxModel_SS_psdueo()
% BoxModel_SS_psdueo: Models 0d marsh and tidal flat time
% evolution by reducing 3 equations to 2, at equilibrium conditions.
%
% Last Update: 11/17/2017
%
%--------------------------------------------------------------------------------------------------
format compact
format longG
clear

%-------------- Set up shared variables with main functions

%-------------- Sediment input constants
C_o = 5 *10^-3;    % ocean concertation (kg/m3)
C_f = 90 *10^-3;    % river concentration (kg/m3)
Q_f = 50;         % river water discharge (m3/s)

%-------------- Erosion constants
k_0 = 1 *10^-3; % roughness (m)
tau_c = 0.3;  % critical shear stress (Pa)
E_0 = 10^-4;    % bed erosion coefficient (kg/m2/s)
k_e = 0.16 /365/24/60/60;  % margin erodibility coefficient (m2/s/W)
v_w = 2;        % reference wind speed (m/s)

% -------------- Accretion constants
k_a = 2;        % margin accretion coefficient

%-------------- Vegetation properties
B_max = 1;      % maximum biomass density (kg/m2)
k_B = 2*10^-3 /365/24/60/60;    % vegetation characteristics (m3/s/kg)

%-------------- Basin properties
b_fm = 1 *10^3; % total basin width (both sides of the channel) (m)
L_E = 4 *10^3; % basin length (m)
R = 2 *10^-3/365/24/60/60;   % sea level rise (m/s)
b_r = 0; % river width (m)

%-------------- Tide Characteristics
T_T = 12 *60*60;   % tidal period (s) (= 12 hours)
a = 2;            % tidal range (m)
H = a/2;          % tidal amplitude (m)

%-------------- Sediment properties
rho_s = 1000;   % sediment bulk density (kg/m3)
omega_s = 0.5 *10^-3;   % settling velocity (m/s)

%-------------- Model constants
gamma = 9800;   % water specific weight (N/m3 or kg/m2/s2)

%-------------- Model assumptions
% Q_f = Q_f/2;    % consider half of the discharge only for one side of the tidal platform (the same will be automatically considered below for Q_T)

%-------------- Initial conditions, x0 = [d_f,d_m]
x0(1) = H+H/2;        % tidal flat depth (m)
x0(2) = H-H/2;         % marsh depth (m)

%-------------- Solve the system
k = 1;
for j = 1 : 10 : b_fm
    
    j
    
    fun = @(x) Fun_BoxModel_SS_depth([j,x]);
    options = optimoptions('fsolve','Display','off');%,'FunctionTolerance',1e-18, 'MaxFunctionEvaluations', 10000,'MaxIterations',10000);%,'Algorithm','trust-region-dogleg','FunctionTolerance',1e-18));
    [x,fval] = fsolve(fun,x0,options);
    y = Fun_BoxModel_SS_width([j,x]);
    
    Sol(k,1:2+length(x)+length(fval)) = [j,x,y,fval];
    k = k+1;
    
%     if y > 0
%         disp('Tidal flat width is:')
%         disp(j)
%         break
%     end
    
end

%-------------- Plot Results
clf
% subplot(2,2,1)

scatter(Sol(:,1),Sol(:,4),'.')
xlabel('TF width (m)')
ylabel('Rate of TF Width (mm/yr)')
box on

% subplot(2,2,2)
%
% scatter(Sol(:,1),Sol(:,5),'.')
% xlabel('TF width (m)')
% ylabel('Rate of TF Depth (mm/yr)')
% box on
%
% subplot(2,2,3)
%
% scatter(Sol(:,1),Sol(:,6),'.')
% xlabel('TF width (m)')
% ylabel('Rate of M Depth (mm/yr)')
% box on
%
% subplot(3,2,4)
%
% scatter(Sol(:,1),Sol(:,2),'.')
% xlabel('TF width (m)')
% ylabel('TF Depth (m)')
% box on
%
% subplot(3,2,5)
%
% scatter(Sol(:,1),Sol(:,3),'.')
% xlabel('TF width (m)')
% ylabel('M Depth (m)')
% box on

% set(findobj('type','axes'),'fontsize',15)
% h_fig=gcf;
% set(h_fig,'PaperOrientation','portrait')
% set(h_fig,'PaperPosition', [0 0 7.5 6]) % [... ... max_width=7.5 max_height=9]
% tit='2';
% print(tit,'-dtiff','-r400')

%======================= Nested Function =========================
    function Feq = Fun_BoxModel_SS_depth (x)
        % function of depths equations at equilibrium
        
        b_f = x(1);
        d_f = x(2);
        d_m = x(3);
        
        %-------------- Imposing a condition for tidal flat conversion to marsh in case of presence of new vegetation when tidal flat is above MSL
        flag_f2m = 0;     % showing that tidal flat is below MSL
        %         if d_f <= H
        %             flag_f2m = 1; % showing that tidal flat is above MSL
        %         end
        
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
            B = B_max*(1-m/r);  % steady state solution for biomass (kg/m2)
            O = k_B*B;            % organic matter production rate
        else
            O = 0;
        end
        
        %-------------- Compute the concentration (kg/m^3)
        C_r = (R-O)*T_T*rho_s/d_m;
        
        %---------------------------------- Define the equations -----------------------------------
        
        
        %-------------------------------- Tidal flat width equation --------------------------------
        
        %-------------- Compute margin erosion (m/s)
        h = (d_f+max(0,d_f-2*H))/2;     % reference water depth
        
        if  flag_f2m==1 || chi<=b_r || v_w==0 % condition for no bed and margin erosion in case of a filled mudflat or no wind
            tau = 0; % bed shear stress
        else
            [ H_w, T_w ] = WaveProps ( h, v_w, chi );   % compute significant height and peak period
            [ tau, k_w ] = ShearStress ( h, k_0, H_w, T_w );     % compute bed shear stress
        end
        
        
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
        
        %-------------- Describe the equation for d_f (m/s)
        Feq(1) = TF_erosion - TF_accretion - SOM + R;
        
        
        %-------------------------------- Mass conservation equation ------------------------------
        
        %-------------- Compute tidal flat bed erosion (kg/s)
        if flag_f2m == 0
            bed_erosion = max(0,t_s*E_0*(tau-tau_c)/tau_c*b_f*L_E);
        else
            bed_erosion = 0;
        end
        
        %-------------- Compute external sediment input (kg/s)
        Vol = (d_f*b_f+d_m*b_m)*L_E; % availble volume in the system to be filled with water
        Q_T = max(Vol/T_T-Q_f,0);
        Q_f (Q_T==0) = Vol/T_T;
        
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
        
        %-------------- Describe the equation for C_r (kg/s)
        Feq(2) = bed_erosion + ocean_in + river_in - TF_deposition - M_deposition - export;   % (kg/s)
        Feq(2) = Feq(2)/rho_s/L_E/b_fm;   % (m/s)
        
        Feq = Feq *1000*60*60*24*365; %(mm/yr)
        
    end

    function F = Fun_BoxModel_SS_width (x)
        % function of width equations at equilibrium
        
        b_f = x(1);
        d_f = x(2);
        d_m = x(3);
        
        %-------------- Imposing a condition for tidal flat conversion to marsh in case of presence of new vegetation when tidal flat is above MSL
        flag_f2m = 0;     % showing that tidal flat is below MSL
        %         if d_f <= H
        %             flag_f2m = 1; % showing that tidal flat is above MSL
        %         end
        
        %-------------- Model assumptions
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
            B = B_max*(1-m/r);  % steady state solution for biomass (kg/m2)
            O = k_B*B;            % organic matter production rate
        else
            O = 0;
        end
        
        %-------------- Compute the concentration (kg/m^3)
        C_r = (R-O)*T_T*rho_s/d_m;
        
        %---------------------------------- Define the equations -----------------------------------
        
        
        %-------------------------------- Tidal flat width equation --------------------------------
        
        %-------------- Compute margin erosion (m/s)
        h = (d_f+max(0,d_f-2*H))/2;     % reference water depth
        
        if  flag_f2m==1 || chi<=b_r || v_w==0 % condition for no bed and margin erosion in case of a filled mudflat or no wind
            W = 0;   % wave power density
        else
            [ H_w, T_w ] = WaveProps ( h, v_w, chi );   % compute significant height and peak period
            [ tau, k_w ] = ShearStress ( h, k_0, H_w, T_w );     % compute bed shear stress
            c_g = pi/k_w/T_w*(1+2*k_w*h/sinh(2*k_w*h)); % wave group velocity (general form)
            W = gamma*c_g*H_w^2/16; % wave power density (kg.m/s3)
        end
        
        B_e = k_e*W;    % margin erosion
        
        %-------------- Compute margin accretion (m/s)
        if flag_f2m == 0
            B_a = k_a*omega_s*C_r/rho_s;  % margin accretion
        else
            B_a = 0;
        end
        
        %-------------- Describe the equation for b_f (m/s)
        F = B_e - B_a;
        F = F *1000*60*60*24*365; % (mm/yr)
        
    end

end