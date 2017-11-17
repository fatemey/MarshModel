function Sol = BoxModel_SS()
% BoxModel_SS: Models 0d marsh and tidal flat time
% evolution at equilibrium conditions.
%
%
% Last Update: 11/17/2017
%
%--------------------------------------------------------------------------------------------------
format compact
format longG
clear

%-------------- Set up shared variables with main functions
load co_data_1000yr; bf_0=dat(:,2);
C_o_V = [5:5:100] *10^-3;

Sol = zeros(length(C_o_V),4);
for i = 1 : length(C_o_V)
    
    i
    %-------------- Sediment input constants
    C_o = C_o_V(i) ;    % ocean concertation (kg/m3)
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
    b_fm = 10*10^3; % total basin width (both sides of the channel) (m)
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
    
    %-------------- Initial conditions, x0=[b_f, d_f, d_m]
    x0(1) = bf_0(i)+200;%b_fm/2;      % tidal flat width (m)
    x0(2) = H+0.3;        % tidal flat depth (m)
    x0(3) = H-0.3;         % marsh depth (m)
    
    %-------------- Solve the system
    objfun = @Fun_BoxModel_SS;
    
    %     D = diag([1/b_fm,1/2/H,1/b_fm]); % used for scaling the problem ->width in thousands of m but depth in 10ths of m??!!
    %     objfun = @(x) Fun_BoxModel_SS(D*x);
    
    %     options = optimoptions('fminunc','Algorithm','quasi-newton','Display','final');
    %     [x,fval,exitflag,output,grad,hessian] = fminunc(objfun,x0',options);
    
    options = optimoptions('fsolve','Algorithm','trust-region','Display','final','FunctionTolerance',1e-18, 'MaxFunctionEvaluations', 10000,'MaxIterations',10000);
    [x,fval] = fsolve(objfun,x0,options);
    
    Sol(i,1:length(x)+length(fval)) = [x,fval];
    
end

% dat = [C_o_V', Sol];
% save('co_data_SS_opt_2.mat','dat')

%-------------- Plot Results
load co_data_1000yr; co=dat(:,1)*1000; bf_TM=dat(:,2);
bf=Sol(:,1); df=Sol(:,2); dm=Sol(:,3); feval=Sol(:,end);

clf
subplot(2,2,1)

hold on
scatter(co,bf_TM/1000,80,'k<','filled')
scatter(co,bf/1000,'ro','filled')
xlabel('Ocean Concentration (mg/l)')
ylabel('Tidal Flat Width (km)')
box on

subplot(2,2,2)

hold on
scatter(co,df*100,80,'b<','filled')
scatter(co,dm*100,'go','filled')
xlabel('Ocean Concentration (mg/l)')
ylabel('Depth (cm)')
legend('Tidal Flat','Marsh','location','northwest')
box on

subplot(2,2,3)

hold on
scatter(co,feval,'ko')
xlabel('Ocean Concentration (mg/l)')
ylabel('Tidal Flat Width Rate (mm/yr)')
box on

% set(findobj('type','axes'),'fontsize',15)
% h_fig=gcf;
% set(h_fig,'PaperOrientation','portrait')
% set(h_fig,'PaperPosition', [0 0 7.5 6]) % [... ... max_width=7.5 max_height=9]
% tit='2';
% print(tit,'-dtiff','-r400')

%======================= Nested Function =========================
    function F = Fun_BoxModel_SS (x)
        % solves the system of equations at equilibrium
        
        %-------------- Setting width boundary limits
        local_bf =  x(1); % imposing a constraint for lower and upper limits of x(1)
        if local_bf < 0
            local_bf = 0;
        end
        if local_bf > b_fm
            local_bf = b_fm;
        end
        
        %-------------- Setting depth boundary limits
        local_df =  x(2); % imposing a constraint for lower limit of x(2)
        if local_df < 0
            local_df = 0;
        end
        
        local_dm =  x(3); % imposing a constraint for lower limit of x(3)
        if local_dm < 0
            local_dm = 0;
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
        
        %-------------- Compute the concentration (kg/m^3)
        C_r = (R-O)*T_T*rho_s/local_dm;
        
        %---------------------------------- Define the equations -----------------------------------
        
        
        %-------------------------------- Tidal flat width equation --------------------------------
        
        %-------------- Compute margin erosion (m/s)
        h = (local_df+max(0,local_df-2*H))/2;     % reference water depth
        
        if  flag_f2m==1 || chi<=b_r || v_w==0 % condition for no bed and margin erosion in case of a filled mudflat or no wind
            tau = 0; % bed shear stress
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
        F(1) = B_e - B_a;
        
        
        %--------------------------------- Tidal flat depth equation -------------------------------
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
        F(2) = TF_erosion - TF_accretion - SOM + R;
        
        
        %-------------------------------- Mass conservation equation ------------------------------
        
        %-------------- Compute tidal flat bed erosion (kg/s)
        if flag_f2m == 0
            bed_erosion = max(0,t_s*E_0*(tau-tau_c)/tau_c*local_bf*L_E);
        else
            bed_erosion = 0;
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
        F(3) = bed_erosion + ocean_in + river_in - TF_deposition - M_deposition - export;   % (kg/s)
        F(3) = F(3)/rho_s/L_E/b_fm;   % (m/s)
        
        %-------------- Describe the objective function (m/s)
        F = F *1000*60*60*24*365; %(mm/yr)
        %         F = norm(F)^2;
        
    end

end