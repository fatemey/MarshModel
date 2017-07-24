function tidalflatmarsh_convertion_linearized
% tidalflatmarsh_convertion_linearized: checks for tidal flat or marsh conversion to marsh or tidal flat based on different set of parametrs.
% here we are intersted in values from critical fetch settings based on BoxModel.
% The saved results can be used in plot_initialwidth function.
%
% Last Update: 5/23/2017
%
%--------------------------------------------------------------------------------------------------
format compact
format longG
clear
% clf

par_temp = 1 : 8;

for k = 2 : 2
    
    %-------------- Set the time span
    tyr = 500;  % solve for time tyr (years)
    ts = tyr *365*24*60*60; % tyr in (s)
    dt = 12*60*60; % time step in (s)
    tspan = 0:dt:ts;
    
    %-------------- Sediment input constants
    C_o = 20 *10^-3;    % ocean concertation (kg/m3)
    C_f = 600 *10^-3;    % river concentration (kg/m3)
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
    R = 3 *10^-3/365/24/60/60;   % sea level rise (m/s)
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
    %     b_fm = b_fm/2;  % consider half of the basin only for one side of the tidal platform
    
    par = par_temp(k);
    switch par
        case 1
            dat = callfun1; d = 10;
        case 2
            dat = callfun2; d = 5;
        case 3
            dat = callfun3; d = 5;
        case 4
            dat = callfun4; d = 1;
        case 5
            dat = callfun5; d = 1;
        case 6
            dat = callfun6; d = 10;
        case 7
            dat = callfun7; d = 5;
        case 8
            dat = callfun8; d = 5;
    end
    
    %-------------- Start the loop for each run
    for j = 1 : size(dat,1)
        
        j
        switch par
            case 1
                C_o = dat(j,1);
            case 2
                C_f = dat(j,1);
            case 3
                Q_f = dat(j,1)/2;
            case 4
                L_E = dat(j,1);
            case 5
                %  b_fm = dat(j,1)/2;
                b_fm = dat(j,1);
            case 6
                v_w = dat(j,1);
            case 7
                R = dat(j,1);
            case 8
                H = dat(j,1);
        end
        
        %-------------- Case 1: Tidal flat emergence
        %-------------- Initial conditions, y0=[ b_f, d_f, d_m,u (=C_r*(b_f*d_f+b_m*d_m))]
        y0(1) = dat(j,2)-d;      % tidal flat width (m)
        y0(2) = H+0.3;        % tidal flat depth (m)
        y0(3) = H-0.3;         % marsh depth (m)
        y0(4) =C_o*(y0(1)*y0(2)+(b_fm-y0(1))*y0(3)); % u
        
        %-------------- Solve the system of differential equations
        [t, y] = ode15s(@ode4marshtidalflat_linearized,tspan,y0); % or use ode15s/ode23s/ode23tb
        t = t /365/24/60/60; % convert time unit from s to yr for plotting purposes
        y(:,4) = y(:,4)./(y(:,1).*y(:,2)+y(:,3).*(b_fm-y(:,1))); % convert y(:,4) to C_r from the formula used before: y4=u (=C_r*(b_f*d_f+b_m*d_m)
        
        ind = find(y(:,3)>H); % remove data related to marsh conversion to tidal flat
        if ~isempty(ind) && length(ind)>1
            y(ind(2):end,:)=[]; % retain only one value afetr conversion to remember it is a new tidal flat now
            t(ind(2):end,:)=[];
        end
        
        ind = find(y(:,2)<=H); % remove data related to tidal flat to marsh conversion
        if ~isempty(ind) && length(ind)>1
            y(ind(2):end,:)=[]; % retain only one value afetr conversion to remember it is a new marsh now
            t(ind(2):end,:)=[];
        end
        
        width = y(:,1); % tidal falt width solution
        n = length(width);
        s = sign(width(n)-width(floor(n/2))); % -1 corresponds to TF contraction (or fully marsh) and +1 accounts for expansion (or fully tidal flat)
        
        %         emrg(j,1:4) = y(n,:);
        %         emrg(j,5) = L_E;
        %         emrg(j,6) = b_fm;
        %         emrg(j,7) = C_o;
        %         emrg(j,8) = C_f;
        
        if y(n,2) <= H % check tidal flat conversion to marsh
            if s == 1
                f(j,1) = 1; % corresponds to an expanding TF which converts to a marsh
            else
                f(j,1) = 2; % corresponds to a contracting TF which converts to a marsh
            end
        else
            f(j,1) = 0;
        end
        
        %-------------- Case 2: marsh drowning
        %-------------- Initial conditions, y0=[ b_f, d_f, d_m,u (=C_r*(b_f*d_f+b_m*d_m))]
        y0(1) = dat(j,2);      % tidal flat width (m)
        
        %-------------- Solve the system of differential equations
        [t, y] = ode15s(@ode4marshtidalflat_linearized,tspan,y0); % or use ode15s/ode23s/ode23tb
        t = t /365/24/60/60; % convert time unit from s to yr for plotting purposes
        y(:,4) = y(:,4)./(y(:,1).*y(:,2)+y(:,3).*(b_fm-y(:,1))); % convert y(:,4) to C_r from the formula used before: y4=u (=C_r*(b_f*d_f+b_m*d_m)
        
        ind = find(y(:,3)>H); % remove data related to marsh conversion to tidal flat
        if ~isempty(ind) && length(ind)>1
            y(ind(2):end,:)=[]; % retain only one value afetr conversion to remember it is a new tidal flat now
            t(ind(2):end,:)=[];
        end
        
        ind = find(y(:,2)<=H); % remove data related to tidal flat to marsh conversion
        if ~isempty(ind) && length(ind)>1
            y(ind(2):end,:)=[]; % retain only one value afetr conversion to remember it is a new marsh now
            t(ind(2):end,:)=[];
        end
        
        width = y(:,1); % tidal falt width solution
        n = length(width);
        s = sign(width(n)-width(floor(n/2))); % -1 corresponds to TF contraction (or fully marsh) and +1 accounts for expansion (or fully tidal flat)
        
        if y(n,3) > H % check marsh conversion to tidal flat
            if width(n) <= 0
                f(j,1) = 3; % corresponds to an expanding marsh which drwons
            elseif width(n) >= b_fm
                f(j,1) = 4; % corresponds to a contracting marsh which drwons
            elseif s == -1
                f(j,1) = 3; % corresponds to an expanding marsh which drwons
            elseif s == 1
                f(j,1) = 4; % corresponds to a contracting marsh which drwons
            end
            
        end
        
    end
    
    switch par
        case 1
            save('co_1_conversion.mat','f')
        case 2
            save('cf_l_conversion.mat','f')
        case 3
            save('qf_1_conversion.mat','f')
        case 4
            save('le_1_conversion.mat','f')
        case 5
            save('bfm_1_conversion.mat','f')
        case 6
            save('vw_1_conversion.mat','f')
        case 7
            save('R_1_conversion.mat','f')
        case 8
            save('H_1_conversion.mat','f')
    end
    
    clear f
    
end

%======================= Nested Functions =========================
    function dy = ode4marshtidalflat_linearized (t,y)
        %  y1=b_f, y2=d_f, y3=d_m, y4=u (=C_r*(b_f*d_f+b_m*d_m, why solving u instead of C_r? u is the variable on the left hand side of mass conservation equation.)
        % solves the ODE system of equations
        
        %-------------- Setting width boundary limits
        local = y(1); % imposing a constraint for lower and upper limits of y(1)
        if (local < 0)
            local = 0 ;
        end
        if (local > b_fm)
            local = b_fm ;
        end
        
        %-------------- Imposing a condition for tidal flat conversion to marsh in case of presence of new vegetation when tidal flat is above MSL
        flag_f2m = 0;     % showing that tidal flat is below MSL
        if y(2) <= H
            flag_f2m = 1; % showing that tidal flat is above MSL
        end
        
        %-------------- Model assumptions
        b_m = b_fm-local; % marsh width
        if flag_f2m == 0
            chi = 2*local+b_r;    % fetch
        elseif flag_f2m == 1
            chi = b_r;
        end
        
        %---------------------------------- Define the equations -----------------------------------
        
        
        %---------------------------- Tidal flat width changes equation ----------------------------
        
        %-------------- Compute margin erosion (m/s)
        h = (y(2)+max(0,y(2)-2*H))/2;     % reference water depth
        
        if  flag_f2m==1 || chi<=b_r || v_w==0 % condition for no bed and margin erosion in case of a filled mudflat or no wind
            tau = 0; % bed shear stress
            W = 0;   % wave power density
        else
            tau = 0.322485 + 0.0000323145 * local - 0.100764 * y(2);
            W = -13.1366 + 0.00286433 * local + 20.1048 * y(2);  
        end
        
        B_e = k_e*W;    % margin erosion
        
        %-------------- Compute margin accretion (m/s)
        C_r = y(4)/(local*y(2)+b_m*y(3)); % concentration (kg/m3)
        
        if flag_f2m == 0
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
        if flag_f2m == 0
            
            %-------------- Compute submerged time when the tidal flat is covered with water
            if y(2) < 2*H
%                 t_s = 1/2-1/pi*asin((H-y(2))/H); 
                t_s = 0.0374284 + 0.619389 * y(2); 
            else
                t_s = 1;
            end
            
            %-------------- Compute the rate of bed erosion (m/s)
            TF_erosion = max(0,t_s*E_0/rho_s*(tau-tau_c)/tau_c);
            
            %-------------- Compute the rate of sediment accretion (m/s)
            %         TF_accretion = t_e*C_r*omega_s/rho_s;
            TF_accretion = min(t_s*C_r*omega_s/rho_s ,C_r*y(2)/T_T/rho_s);
            
            %-------------- Compute the rate of organic matter production in tidal flat (m/s)
            SOM = 0;
            
        elseif flag_f2m == 1
            
            %-------------- Compute the rate of bed erosion (m/s)
            TF_erosion = 0;
            
            %-------------- Compute the rate of sediment accretion (m/s)
            TF_accretion = C_r*y(2)/T_T/rho_s;
            
            %-------------- Compute the rate of organic matter production in the new marsh (m/s)
            z_new = H-y(2);       % elevation of marsh platform
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
        z = H-y(3);       % elevation of marsh platform
        if z >= 0           % condition for no presence of vegetation when marsh is below MSL
            r = -0.5*z/H+1;     % reproduction rate
            m = 0.5*z/H;         % mortality rate
            B = B_max*(1-m/r);  % steady state solution for biomass (kg/m2)
            O = k_B*B;            % organic matter production rate
        else
            O = 0;
        end
        
        %-------------- Describe the equation for d_m (m)
        dy(3,1) = - M_accretion - O + R;
        
        %-------------------------------- Mass conservation equation ------------------------------
        
        %-------------- Compute tidal flat bed erosion (kg/s)
        if flag_f2m == 0
            bed_erosion = max(0,t_s*E_0*(tau-tau_c)/tau_c*local*L_E);
        else
            bed_erosion = 0;
        end
        
        %-------------- Compute marsh/tidal flat margin erosion/accretion (kg/s)
        if flag_f2m == 0
            margin = (B_e - B_a)*(y(2)-y(3))*L_E*rho_s;
        else
            margin = 0;
        end
        
        %-------------- Compute external sediment input (kg/s)
        Q_T = (y(2)*local+y(3)*b_m)*L_E/T_T-Q_f;
        ocean_in = Q_T*C_o;
        river_in = Q_f*C_f;
        
        %-------------- Compute deposition on tidal flat (kg/s)
        if flag_f2m == 0
            %         TF_deposition = t_e*C_r*local*omega_s*L_E;
            TF_deposition = min(t_s*C_r*local*omega_s*L_E, C_r*local*y(2)*L_E/T_T);
            %         TF_deposition = alpha*dt*(y(2)/H)*local*L_E*rho_s;
        else
            TF_deposition = C_r*local*y(2)*L_E/T_T;
        end
        
        %-------------- Compute deposition on marsh (kg/s)
        M_deposition = C_r*b_m*y(3)*L_E/T_T;
        
        %-------------- Compute export sediment to the ocean (kg/s)
        if flag_f2m == 0
            export = C_r*local*min(y(2),2*H)*L_E/T_T;
        else
            export = 0;
        end
        
        %-------------- Describe the equation for C_r (kg/m3)
        var = bed_erosion + margin + ocean_in + river_in - TF_deposition - M_deposition - export;   % (kg/s)
        dy(4,1) =  var / L_E;
        
    end

    function dat = callfun1
        % loading data suitable for static workspace
        load co_1 dat
    end

    function dat = callfun2
        load cf_l dat
    end

    function dat = callfun3
        load qf_1 dat
    end

    function dat = callfun4
        load le_1 dat
    end

    function dat = callfun5
        load bfm_1 dat
    end

    function dat = callfun6
        load vw_1 dat
    end

    function dat = callfun7
        load R_1 dat
    end

    function dat = callfun8
        load H_1 dat
    end

end
