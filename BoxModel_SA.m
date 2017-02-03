function BoxModel_SA (parin)
% BoxModel_SA: Modeling 0d marsh and tidal flat time
% evolution using Matlab ode15s function for sensitivity analysis
%
% Purpose: Determining marsh and tidal flat depths and widths changes with
%          rising sea level. This model is solved using 4 equations and 4
%          unknowns.
%
%--------------------------------------------------------------------------------------------------
format compact
format longG

%-------------- Initial model adjustment for different runs
if nargin == 0 || strcmp(parin,'all') || unique(parin == 0)
    par_temp = 1:14;
    nj = length(par_temp);
else
    par = parin;
    nj = 1;
end

for j = 1 : nj
    
    tic
    
    if nj>1
        par = par_temp(j);
    end
    
    switch par
        case {'C_o' , 1}
            C_o_ = [10*10^-3 20*10^-3 100*10^-3]; leg = {'10 g/m^3','20 g/m^3','100 g/m^3'}; tit = 'Ocean Concentration';
        case {'C_f' , 2}
            C_f_ = [10*10^-3 20*10^-3 100*10^-3]; leg = {'0 g/m^3','15 g/m^3','100 g/m^3'}; tit = 'River Concentration';
        case {'Q_f' , 3}
            Q_f_ = [0 20 1000]; leg = {'0 m^3/s','20 m^3/s','1000 m^3/s'}; tit = 'River Discharge';
        case {'tau_c' , 4}
            tau_c_ = [.1 .3 .4 ]; leg = {'0.1 PA','0.3 PA','0.4 PA'}; tit = 'Critical Shear Stress';
        case {'E_0' , 5}
            E_0_ = [10^-5 10^-4 10^-3]; leg = {'10^{-5} kg/m^2/s','10^{-4} kg/m^2/s','10^{-3} kg/m^2/s'}; tit = 'Bed Erosion Coefficient';
        case {'k_e' , 6}
            k_e_ = [0.1/365/24/60/60 0.16/365/24/60/60 0.2/365/24/60/60]; leg = {'0.1 m^2/yr/W','0.16 m^2/yr/W','0.2 m^2/yr/W'}; tit = 'Margin Erodibility Coefficient';
        case {'v_w' , 7}
            v_w_ = [0 6 10]; leg = {'0 m/s','6 m/s','10 m/s'}; tit = 'Wind Velocity';
        case {'k_a' , 8}
            k_a_ = [1 2 4]; leg = {'1','2','4'}; tit ='Margin Accretion Coefficient';
        case {'k_B' , 9}
            k_B_ = [1*10^-3/365/24/60/60 2*10^-3/365/24/60/60 4*10^-3/365/24/60/60]; leg = {'1 m3/yr/g','2 m3/yr/g','4 m3/yr/g'}; tit = 'Biomass Coefficient';
        case {'R' , 10}
            R_ = [0 2*10^-3/365/24/60/60 10*10^-3/365/24/60/60]; leg = {'0 mm/yr','2 mm/yr','10 mm/yr'}; tit = 'Rate of Sea Level Rise';
        case {'b_fm' , 11}
            b_fm_ = [10^3 5*10^3 10*10^3];  leg = {'0.5 km','2.5 km','5 km'}; tit = 'Basin Width';
        case {'L_E' , 12}
            L_E_ = [.15*10^3 15*10^3 150*10^3]; leg = {'0.15 km','15 km','150 km'}; tit = 'Estuary Length';
        case {'b_f_0' , 13}
            b_f_0 = [1*10^3  5*10^3  9*10^3]/4; leg = {'0.25 km (0.1 L)','1.25 km (0.5 L)','2.25 km (0.9 L)'}; tit = 'Initial Tidal Flat Width';
        case {'H' , 14}
            H_ = [1 1.4 2.4]/2; leg = {'1 m','1.4 m','2.4 m'}; tit = 'Tidal Range';
    end
    
    %-------------- Start the loop for each run
    for i = 1 : length(leg)
        
        %-------------- Set the time span
        tyr = 1000;  % solve for time tyr (years)
        ts = tyr *365*24*60*60; % tyr in (s)
        dt = 12*60*60; % time step in (s)
        tspan = 0:dt:ts;
        
        %-------------- Sediment input constants
        C_o = 20 *10^-3;    % ocean concertation (kg/m3)
        C_f = 15 *10^-3;    % river concentration (kg/m3)
        Q_f = 20;         % river water discharge (m3/s)
        
        %-------------- Erosion constants
        k_0 = 1 *10^-3; % roughness (m)
        tau_c = 0.3;  % critical shear stress (Pa)
        E_0 = 10^-4;    % bed erosion coefficient (kg/m2/s)
        k_e =  0.16 /365/24/60/60;  % margin erodibility coefficient (m2/s/W)
        v_w = 6;        % reference wind speed (m/s)
        
        % -------------- Accretion constants
        k_a = 2;        % margin accretion coefficient
        
        %-------------- Vegetation properties
        B_max = 1;      % maximum biomass density (kg/m2)
        k_B = 2*10^-3 /365/24/60/60;    % vegetation characteristics (m3/s/kg)
        
        %-------------- Basin properties
        b_fm = 5 *10^3; % total basin width (both sides of the channel) (m)
        L_E = 15 *10^3; % basin length (m)
        R = 2 *10^-3/365/24/60/60;   % sea level rise (m/s)
        b_r = 50; % river width (m)
        
        %-------------- Tide Characteristics
        T_T = 12 *60*60;   % tidal period (s) (= 12 hours)
        H = 1.4 /2;          % tidal amplitude (range/2) (m)
        
        %-------------- Sediment properties
        rho_s = 1000;   % sediment bulk density (kg/m3)
        omega_s = 0.5 *10^-3;   % settling velocity (m/s)
        
        %-------------- Model constants
        gamma = 9800;   % water specific weight (N/m3)
        g = 9.81;       % gravitational acceleration (m/s2)
        
        %-------------- Initial conditions, y0=[ b_f, d_f, d_m,u(=C_r*(b_f*d_f+b_m*d_m))]
        y0(1) = b_fm/4;   % tidal flat width (m)
        y0(2) = 1.0;         % tidal flat depth (m)
        y0(3) = 0.4;         % marsh depth (m)
        y0(4) = 0*10^-3*(y0(1)*(y0(2)+y0(3)));
        
        %-------------- Parameter assignment
        
        switch par
            case {'C_o' , 1}
                C_o = C_o_(i);
            case {'C_f' , 2}
                C_f = C_f_(i);
            case {'Q_f' , 3}
                Q_f = Q_f_(i);
            case {'tau_c' , 4}
                tau_c = tau_c_(i);
            case {'E_0' , 5}
                E_0 = E_0_(i);
            case {'k_e' , 6}
                k_e = k_e_(i);
            case {'v_w' , 7}
                v_w = v_w_(i);
            case {'k_a' , 8}
                k_a = k_a_(i);
            case {'k_B' , 9}
                k_B = k_B_(i);
            case {'R' , 10}
                R = R_(i);
            case {'b_fm' , 11}
                b_fm = b_fm_(i);
            case {'L_E' , 12}
                L_E = L_E_(i);
            case {'b_f_0' , 13}
                y0(1) = b_f_0(i);
            case {'H' , 14}
                H = H_(i);
        end
        
        %-------------- Model assumptions
        Q_f = Q_f/2;    % consider half of the discharge only for one side of the tidal platform (the same will be automatically considered below for Q_T)
        b_fm = b_fm/2;  % consider half of the basin only for one side of the tidal platform
        
        %-------------- Solve the system of differential equations
        [t, y] = ode15s(@ode4marshtidalflat,tspan,y0); % or use ode15s/ode45/..23s
        y(:,4) = y(:,4)./(y(:,1).*y(:,2)+y(:,3).*(b_fm-y(:,1))); % convert y(:,4) to C_r from the formula used before: y4=u (=C_r*(b_f*d_f+b_m*d_m)
        ydata(:,:,i) = y; % save the solution for each SA value (3 values for the factor of interest)
        
    end
    
    %-------------- Plot Results
%     figure
    plot_BoxModel_SA(t,ydata,leg,tit)
        print(tit,'-dtiff','-r400')
        movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\12 - TF conversion to M corrected')
        close all
    
    timespent_min = toc/60
    
end

%======================= Nested Function =========================
    function dy = ode4marshtidalflat (t,y) %  y1=b_f, y2=d_f, y3=d_m, y4=u (=C_r*(b_f*d_f+b_m*d_m, why solving u instead of C_r? u is the variable on the left hand side of mass conservation equation.)
        % solves the ODE system of equations
        
        %-------------- Setting width boundary limits
        local = y(1); % imposing a constraint for lower and upper limits of y(1)
        if (local < 0)
            local = 0 ;
        end
        if (local > b_fm)
            local = b_fm ;
        end
        
        %-------------- Imposing a condition for when tidal flat turns into marsh
        flag = 0; % showing that tidal flat is below MSL
        if y(2) <= H
            flag = 1; % showing that tidal flat is above MSL
        end
        
        %-------------- Model assumptions
        b_m = b_fm-local; % marsh width
        if flag == 0
            chi = 2*local+b_r;    % fetch
        elseif flag == 1
            chi = b_r;
        end
        
        %---------------------------------- Define the equations -----------------------------------
        
        
        %---------------------------- Tidal flat width changes equation ----------------------------
        
        %-------------- Compute margin erosion (m/s)
        h = (y(2)+max(0,y(2)-2*H))/2;     % reference water depth
        
        if  flag==1 || chi<=b_r || v_w==0 % condition for no bed and margin erosion in case of a filled mudflat or no wind
            tau = 0; % bed shear stress
            W = 0;   % wave power density
        else
            [ H_w, T_w ] = WaveProps ( h, v_w, chi, g );   % compute significant height and peak period
            [ tau, k_w ] = ShearStress ( h, k_0, H_w, T_w, g );     % compute bed shear stress
            % c_g = sqrt(g*h);        % wave group velocity (shallow water)
            c_g = pi/k_w/T_w*(1+2*k_w*h/sinh(2*k_w*h)); % wave group velocity (general form)
            W = gamma*c_g*H_w^2/16; % wave power density (kg.m/s3)
        end
        
        B_e = k_e*W;    % margin erosion
        
        %-------------- Compute margin accretion (m/s)
        C_r = y(4)/(local*y(2)+b_m*y(3)); % concentration (kg/m3)
        
        if flag == 0
            B_a = k_a*omega_s*C_r/rho_s;  % margin accretion
        else
            B_a = 0;
        end
        
        %-------------- Describe the equation for b_f (m)
        dy(1,1) = B_e - B_a;
        
        if (y(1) < 0 && dy(1,1) < 0) % imposing a constraint for boundary limits of b_f (or y(1))
            dy(1,1) = 0 ;
        end
        if (y(1) > b_fm && dy(1,1) > 0)
            dy(1,1) = 0 ;
        end
        
        
        %----------------------------- Tidal flat depth changes equation ---------------------------
        if flag == 0
            
            %-------------- Compute submerged time when the tidal flat is covered with water
            if y(2) < 2*H
                t_s = 1/2-1/pi*asin((H-y(2))/H);
            else
                t_s = 1;
            end
            
            %-------------- Compute the rate of bed erosion (m/s)
            TF_erosion = max(0,t_s*E_0/rho_s*(tau-tau_c)/tau_c);
            
            %-------------- Compute the rate of sediment accretion (m/s)
            %         TF_accretion = t_e*C_r*omega_s/rho_s;
            TF_accretion = min(t_s*C_r*omega_s/rho_s ,C_r*y(2)/T_T/rho_s);
            %         if C_o == 10*10^-3 % condition to consider alpha based on D'Alpaos et al 2011 approach
            %             alpha = 4 *10^-3/365/24/60/60; % (m/s)
            %         elseif C_o == 20*10^-3
            %             alpha = 8 *10^-3/365/24/60/60;
            %         elseif C_o == 100*10^-3
            %             alpha = 38 *10^-3/365/24/60/60;
            %         end
            %         alpha = 0.38 * C_r;
            %         TF_accretion = alpha*dt*(y(2)/H); % based on D'Alpaos et al 2011 approach during one time step
            
            %-------------- Compute the rate of organice matter production in tidal flat (m/s)
            SOM = 0;
            
        elseif flag == 1
            
            %-------------- Compute the rate of bed erosion (m/s)
            TF_erosion = 0;
            
            %-------------- Compute the rate of sediment accretion (m/s)
            TF_accretion = C_r*y(2)/T_T/rho_s;
            
            %-------------- Compute the rate of organice matter production in the new marsh (m/s)
            z_new = H-y(2);            % elevation of marsh platform
            r_new = -0.5*z_new/H+1;     % reproduction rate
            m_new = 0.5*z_new/H;         % mortality rate
            B_new = B_max*(1-m_new/r_new);  % steady state solution for biomass (kg/m2)
            SOM = k_B*B_new;        % organic matter production rate
            
        end
        
        %-------------- Describe the equation for d_f (m)
        dy(2,1) = TF_erosion - TF_accretion - SOM + R;
        
        
        %----------------------------- Marsh depth changes equation ------------------------------
        
        %-------------- Compute the rate of sediment accretion (m/s)
        M_accretion = C_r*y(3)/T_T/rho_s;
        
        %-------------- Compute organice matter production (m/s)
        z = H-y(3);            % elevation of marsh platform
        r = -0.5*z/H+1;     % reproduction rate
        m = 0.5*z/H;         % mortality rate
        B = B_max*(1-m/r);  % steady state solution for biomass (kg/m2)
        O = k_B*B;             % organic matter production rate
        
        %-------------- Describe the equation for d_m (m)
        dy(3,1) = - M_accretion - O + R;
        
        %-------------------------------- Mass conservation equation ------------------------------
        
        %-------------- Compute tidal flat bed erosion (kg/s)
        if flag == 0
            bed_erosion = max(0,t_s*E_0*(tau-tau_c)/tau_c*local*L_E);
        else
            bed_erosion = 0;
        end
        
        %-------------- Compute marsh/tidal flat margin erosion/accretion (kg/s)
        if flag == 0
            margin = (B_e - B_a)*(y(2)-y(3))*L_E*rho_s;
        else
            margin = 0;
        end
        
        %-------------- Compute external sediment input (kg/s)
        Q_T = (y(2)*local+y(3)*b_m)*L_E/T_T-Q_f;
        ocean_in = Q_T*C_o;
        river_in = Q_f*C_f;
        
        %-------------- Compute deposition on tidal flat (kg/s)
        if flag == 0
            %         TF_deposition = t_e*C_r*local*omega_s*L_E;
            TF_deposition = min(t_s*C_r*local*omega_s*L_E, C_r*local*y(2)*L_E/T_T);
            %         TF_deposition = alpha*dt*(y(2)/H)*local*L_E*rho_s;
        else
            TF_deposition = C_r*local*y(2)*L_E/T_T;
        end
        
        %-------------- Compute deposition on marsh (kg/s)
        M_deposition = C_r*b_m*y(3)*L_E/T_T;
        
        %-------------- Compute export sediment to the ocean (kg/s)
        if flag == 0
            export = C_r*local*min(y(2),2*H)*L_E/T_T;
        else
            export = 0;
        end
        
        %-------------- Describe the equation for C_r (kg/m3)
        var = bed_erosion + margin + ocean_in + river_in - TF_deposition - M_deposition - export;   % (kg/s)
        dy(4,1) =  var / L_E;
        
    end

end
