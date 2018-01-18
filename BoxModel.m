function BoxModel
% BoxModel: Models 0d marsh and tidal flat time
% evolution using Matlab ode15s function.
%
% Output
%           t : vector of time data in yr
%           y : matrix of data. 1st vertical vector: tidal flat width (m),
%           2nd vector: tidal flat depth (m), 3rd vector: marsh depth (m),
%           and 4th vector: C_r (g/m3).
%
% Purpose: Determining marsh and tidal flat depths and widths changes with
%               rising sea level. This model is solved using 4 equations and 4
%               unknowns.
%
% Last Update: 1/17/2018
%
%--------------------------------------------------------------------------------------------------
format compact
format longG
% clear
% clf

%-------------- Set the time span
tyr = 1000;  % solve for time tyr (years)
ts = tyr *365*24*60*60; % tyr in (s)
dt = 24*1000*60*60; % time step in (s)
tspan = 0:dt:ts;

%-------------- Sediment input constants
C_o = 5 *10^-3;    % ocean concertation (kg/m3)
C_f = 35 *10^-3;    % river concentration (kg/m3)
Q_f = 5;         % river water discharge (m3/s)

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

%**************************************************************************
% x=[ 0.01,  0.14,     100,   5000,   5000,    1.4,   6.0,     12,      6];
% C_o=x(1);C_f=x(2);Q_f=x(3);L_E=x(4);b_fm=x(5);a=x(6);R=x(7)*10^-3/365/24/60/60;T_T=x(8)*60*60;v_w=x(9);

C_o = 10*10^-3; % (0:5:100)*10^-3 , 10*10^-3
C_f = (62)*10^-3; % (0:5:100)*10^-3 , [5 50]*10^-3
Q_f = [ 5]; % (0:5:100), [5 50]
L_E = 5*1000;
b_fm = 5*1000;
a = 1.4; % [.5, 1:20] , 1.4
R = [6] *10^-3/365/24/60/60; % [0,2:4:18]
T_T = 12*60*60;
v_w = 6;

%**************************************************************************
H = a/2;          % tidal amplitude (m)

%-------------- Sediment properties
rho_s = 1000;   % sediment bulk density (kg/m3)
omega_s = 0.5 *10^-3;   % settling velocity (m/s)

%-------------- Model constants
gamma = 9800;   % water specific weight (N/m3 or kg/m2/s2)

%-------------- Model assumptions
% Q_f = Q_f/2;    % consider half of the discharge only for one side of the tidal platform (the same will be automatically considered below for Q_T)
% b_fm = b_fm/2;  % consider half of the basin only for one side of the tidal platform

%-------------- Initial conditions, y0=[ b_f, d_f, d_m,u (=C_r*(b_f*d_f+b_m*d_m))]
y0(1) =[147];%y1;%b_fm/2;      % tidal flat width (m)
y0(2) = H+H/2;        % tidal flat depth (m)
y0(3) = H-H/2;         % marsh depth (m)
y0(4) =C_o*(y0(1)*y0(2)+(b_fm-y0(1))*y0(3)); % u

%-------------- Solve the system of differential equations
[t, y] = ode15s(@ode4marshtidalflat,tspan,y0); % or use ode15s/ode23s/ode23tb
t = t /365/24/60/60; % convert time unit from s to yr for plotting purposes
y(:,4) = y(:,4)./(y(:,1).*y(:,2)+y(:,3).*(b_fm-y(:,1))); % convert y(:,4) to C_r from the equation used before: y4=u (=C_r*(b_f*d_f+b_m*d_m)

%-------------- Removing data corresponding to platform conversion and reaching to basin boundary limits
% ind = find(y(:,3)>H); % marsh conversion to tidal flat
% if ~isempty(ind) && length(ind)>1
%     y(ind(2):end,:)=[]; % retain only one value after conversion to remember it is a new tidal flat now
%     t(ind(2):end,:)=[];
% end

ind = find(y(:,2)<=H); % tidal flat conversion to marsh
if ~isempty(ind) && length(ind)>1
    y(ind(2):end,:)=[]; % retain only one value after conversion to remember it is a new marsh now
    t(ind(2):end,:)=[];
end

ind = find(y(:,1)>=b_fm); % tidal flat filling the basin
if ~isempty(ind) && length(ind)>1
    y(ind(2):end,:)=[]; % retain only one value after conversion to remember it is a new marsh now
    t(ind(2):end,:)=[];
end

ind = find(y(:,1)<=0); % marsh filling the basin
if ~isempty(ind) && length(ind)>1
    y(ind(2):end,:)=[]; % retain only one value after conversion to remember it is a new marsh now
    t(ind(2):end,:)=[];
end

%-------------- Plot Results
% plot_BoxModel(t,y)
plot_BoxModel_pars(t,y,C_o,C_f,Q_f,L_E,b_fm,a,T_T,v_w)

%======================= Nested Function =========================
    function dy = ode4marshtidalflat (t,y) 
        %  y1=b_f, y2=d_f, y3=d_m, y4=u (=C_r*(b_f*d_f+b_m*d_m, why solving u instead of C_r? u is the variable on the left hand side of mass conservation equation.)
        % solves the ODE system of equations
        
        %-------------- Setting width boundary limits
        local_bf = y(1); % imposing a constraint for lower and upper limits of y(1)
        if local_bf < 0
            local_bf = 0 ;
        end
        if local_bf > b_fm
            local_bf = b_fm ;
        end
        
        %-------------- Setting depth boundary limits
        local_df = y(2); % imposing a constraint for lower limit of y(2)
        if local_df < 0
            local_df = 0 ;
        end
        
        local_dm = y(3); % imposing a constraint for lower limit of y(3)
        if local_dm < 0
            local_dm = 0 ;
        end
        
        %-------------- Imposing a condition for tidal flat conversion to marsh in case of presence of new vegetation when tidal flat is above MSL
        flag_f2m = 0;     % showing that tidal flat is below MSL
        if local_df <= H
            flag_f2m = 1; % showing that tidal flat is above MSL
        end
        
        %-------------- Model assumptions
        b_m = b_fm-local_bf; % marsh width
        if flag_f2m == 0
            chi = 2*local_bf+b_r;    % fetch
        elseif flag_f2m == 1
            chi = b_r;
        end
        
        %---------------------------------- Define the equations -----------------------------------
        
        
        %---------------------------- Tidal flat width changes equation ----------------------------
        
        %-------------- Compute margin erosion (m/s)
        h = (local_df+max(0,local_df-2*H))/2;     % reference water depth
        
        if  flag_f2m==1 || chi<=b_r || v_w==0 % condition for no bed and margin erosion in case of a filled mudflat or no wind
            tau = 0; % bed shear stress
            W = 0;   % wave power density
        else
            [ H_w, T_w ] = WaveProps ( h, v_w, chi);   % compute significant height and peak period
            [ tau, k_w ] = ShearStress ( h, k_0, H_w, T_w);     % compute bed shear stress
            % c_g = sqrt(g*h);        % wave group velocity (shallow water)
            c_g = pi/k_w/T_w*(1+2*k_w*h/sinh(2*k_w*h)); % wave group velocity (general form)
            W = gamma*c_g*H_w^2/16; % wave power density (kg.m/s3)
        end
        
        B_e = k_e*W;    % margin erosion
        
        %-------------- Compute margin accretion (m/s)
        C_r = y(4)/(local_bf*local_df+b_m*local_dm); % concentration (kg/m3)
        
        if flag_f2m == 0
            B_a = k_a*omega_s*C_r/rho_s;  % margin accretion
        else
            B_a = 0;
        end
        
        %-------------- Describe the equation for b_f (m/s)
        dy(1,1) = B_e - B_a;
        
        if (y(1) < 0 && dy(1,1) < 0) % imposing a constraint for boundary limits of b_f (or y(1))
            dy(1,1) = 0 ;
        end
        if (y(1) > b_fm && dy(1,1) > 0)
            dy(1,1) = 0 ;
        end
        
        
        %----------------------------- Tidal flat depth changes equation ---------------------------
        if flag_f2m == 0
            
            %-------------- Compute submerged time when the tidal flat is covered with water
            if local_df < 2*H
                t_s = 1/2-1/pi*asin((H-local_df)/H);
            else
                t_s = 1;
            end
            
            %-------------- Compute the rate of bed erosion (m/s)
            TF_erosion = max(0,t_s*E_0/rho_s*(tau-tau_c)/tau_c);
            
            %-------------- Compute the rate of sediment accretion (m/s)
            TF_accretion = min(t_s*C_r*omega_s/rho_s ,C_r*local_df/T_T/rho_s);
            
            %-------------- Compute the rate of organic matter production in tidal flat (m/s)
            SOM = 0;
            
        elseif flag_f2m == 1
            
            %-------------- Compute the rate of bed erosion (m/s)
            TF_erosion = 0;
            
            %-------------- Compute the rate of sediment accretion (m/s)
            TF_accretion = C_r*local_df/T_T/rho_s;
            
            %-------------- Compute the rate of organic matter production in the new marsh (m/s)
            z_new = H-local_df;       % elevation of marsh platform
            r_new = -0.5*z_new/H+1;     % reproduction rate
            m_new = 0.5*z_new/H;         % mortality rate
            B_new = B_max*(1-m_new/r_new);  % steady state solution for biomass (kg/m2)
            SOM = k_B*B_new;        % organic matter production rate
            
        end
        
        %-------------- Describe the equation for d_f (m/s)
        dy(2,1) = TF_erosion - TF_accretion - SOM + R;
        
        
        %----------------------------- Marsh depth changes equation ------------------------------
        
        %-------------- Compute the rate of sediment accretion (m/s)
        M_accretion = C_r*local_dm/T_T/rho_s;
        
        %-------------- Compute organice matter production (m/s)
        z = H-local_dm;       % elevation of marsh platform
        if z >= 0           % condition for presence of vegetation when marsh is above MSL
            r = -0.5*z/H+1;     % reproduction rate
            m = 0.5*z/H;         % mortality rate
            B = B_max*(1-m/r);  % steady state solution for biomass (kg/m2)
            O = k_B*B;            % organic matter production rate
        else
            O = 0;
        end
        
        %-------------- Describe the equation for d_m (m/s)
        dy(3,1) = - M_accretion - O + R;
        
        %-------------------------------- Mass conservation equation ------------------------------
        
        %-------------- Compute tidal flat bed erosion (kg/s)
        if flag_f2m == 0
            bed_erosion = max(0,t_s*E_0*(tau-tau_c)/tau_c*local_bf*L_E);
        else
            bed_erosion = 0;
        end
        
        %-------------- Compute marsh/tidal flat margin erosion/accretion (kg/s)
        if flag_f2m == 0
            margin = (B_e - B_a)*(local_df-local_dm)*L_E*rho_s;
        else
            margin = 0;
        end
        
        %-------------- Compute external sediment input (kg/s)
        Vol = (local_df*local_bf+local_dm*b_m)*L_E; % availble volume in the system to be filled with water
        Q_T = max(Vol/T_T-Q_f,0);
        Q_f (Q_T==0) = Vol/T_T;
        
        ocean_in = Q_T*C_o;
        river_in = Q_f*C_f;
        
        %-------------- Compute deposition on tidal flat (kg/s)
        if flag_f2m == 0
            TF_deposition = min(t_s*C_r*local_bf*omega_s*L_E, C_r*local_bf*local_df*L_E/T_T);
        else
            TF_deposition = C_r*local_bf*local_df*L_E/T_T;
        end
        
        %-------------- Compute deposition on marsh (kg/s)
        M_deposition = C_r*b_m*local_dm*L_E/T_T;
        
        %-------------- Compute export sediment to the ocean (kg/s)
        if flag_f2m == 0
            export = C_r*local_bf*min(local_df,2*H)*L_E/T_T;
        else
            export = 0;
        end
        
        %-------------- Describe the equation for C_r (kg/s)
        var = bed_erosion + margin + ocean_in + river_in - TF_deposition - M_deposition - export;   % (kg/s)
        dy(4,1) =  var / L_E;
        
    end

end
