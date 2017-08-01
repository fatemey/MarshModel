function BoxModel_SS_opt()
% BoxModel_SS: Models 0d marsh and tidal flat time
% evolution using optimization and in equilibrium conditions.
% This approacche works better than BoxModel_SS
%
% Last Update: 7/5/2017
%
%--------------------------------------------------------------------------------------------------
format compact
format longG
clear
% clf

%-------------- Set up shared variables with OUTFUN
history.x = [];
history.fval = [];
% searchdir = [];

%-------------- Set up shared variables with main functions
load co_data; bf_0=dat(:,2); %bf_0(19:20)=[1500,1500];
C_o_V = 5:5:100;

Sol = zeros(length(C_o_V),4);
for i = 1 : length(C_o_V)
    i
    %-------------- Sediment input constants
    C_o = C_o_V(i) *10^-3;    % ocean concertation (kg/m3)
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
    b_fm = 5 *10^3; % total basin width (both sides of the channel) (m)
    L_E = 15 *10^3; % basin length (m)
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
    x0(1) = bf_0(i);%1570;%b_fm/2;      % tidal flat width (m)
    x0(2) = H+0.3;        % tidal flat depth (m)
    x0(3) = H-0.3;         % marsh depth (m)
    
    %-------------- Solve the system
    %     lb = [0,-Inf,-Inf];
    lb = [0,0,0];
    ub = [b_fm,Inf,Inf];
    objfun = @Fun_BoxModel_SS;
    confun = @Fun_BoxModel_SS_con;
    options = optimoptions('fmincon','Algorithm','interior-point','Display','final','StepTolerance',1e-200,'ConstraintTolerance',1e-15,'MaxFunctionEvaluations',100000,'MaxIterations',100000,'OptimalityTolerance',1e-20,'OutputFcn',@outfun); %'Algorithm','active-set','sqp'
%     options = optimoptions('fmincon','Algorithm','active-set','Display','final','StepTolerance',1e-300,'ConstraintTolerance',1e-100,'MaxFunctionEvaluations',1000000,'MaxIterations',1000000,'OptimalityTolerance',1e-30,'OutputFcn',@outfun); %'Algorithm','active-set','sqp'
    [x,fval,exitflag,output,lambda,grad,hessian] = fmincon(objfun,x0,[],[],[],[],lb,ub,confun,options);
    
    Sol(i,1:length(x)+length(fval)) = [x,fval];
%     [G,Geq] = Fun_BoxModel_SS_con(x);
%     Sol(i,length(x)+length(fval)+1:length(x)+length(fval)+2) = Geq;

end

dat = [C_o_V'*.001, Sol];
save('co_data_SS_opt.mat','dat')

%-------------- Plot Results
% figure(1)
% clf
% hold on
% yyaxis left
% bf=history.x(:,1);
% scatter(1:length(bf),bf)
% ylabel('Tidal Flat Width (m)')
% yyaxis right
% df=history.x(:,2);
% dm=history.x(:,3);
% scatter(1:length(bf),df)
% scatter(1:length(bf),dm,'filled')
% xlabel('Iteration')
% ylabel('Depth (m)')
% box on
% 
% figure(2)
% clf
% scatter(1:length(bf),history.fval,'k','o')
% xlabel('Iteration')
% ylabel('Tidal Flat Width Rate (mm/yr)')
% box on

figure(1)
clf
hold on
df=Sol(:,2);
dm=Sol(:,3);
scatter(C_o_V,df,60,'b<','filled')
scatter(C_o_V,dm,'go','filled')
xlabel('Ocean Concentration (mg/l)')
ylabel('Depth (m)')
legend('Tidal Flat','Marsh','location','northwest')
box on

figure(2)
clf
hold on
co=C_o_V; feval=Sol(:,end);
scatter(co,feval,'ko')
xlabel('Ocean Concentration (mg/l)')
ylabel('Tidal Flat Width Rate (mm/yr)')
box on

figure(3)
clf
hold on
load co_data; co1=dat(:,1)*1000; bf1=dat(:,2);
co2=C_o_V; bf2=Sol(:,1);
scatter(co1(1:7),bf1(1:7),55,'k','>','filled')
scatter(co1(8:end),bf1(8:end),55,'k','o','filled')
scatter(co2,bf2,'ro','filled')
xlabel('Ocean Concentration (mg/l)')
ylabel('Tidal Flat Width (m)')
box on

% figure(4)
% clf
% hold on
% g1=Sol(:,end-1);
% g2=Sol(:,end);
% scatter(C_o_V,g1,'b<','filled')
% scatter(C_o_V,g2,'go','filled')
% xlabel('Ocean Concentration (mg/l)')
% ylabel('Constraintents (mm/yr)')
% legend('C1','C2','location','northwest')
% box on

% set(findobj('type','axes'),'fontsize',15)
% h_fig=gcf;
% set(h_fig,'PaperOrientation','portrait')
% set(h_fig,'PaperPosition', [0 0 7.5 6]) % [... ... max_width=7.5 max_height=9]
% tit='3-initialguess';
% print(tit,'-dtiff','-r400')

%======================= Nested Function =========================
    function F = Fun_BoxModel_SS (x)
        % solves the system of equations at equilibrium
        
        b_f = x(1);
        d_f = x(2);
        d_m = x(3);
        
        %-------------- Imposing a condition for tidal flat conversion to marsh in case of presence of new vegetation when tidal flat is above MSL
        flag_f2m = 0;     % showing that tidal flat is below MSL
        if d_f <= H
            flag_f2m = 1; % showing that tidal flat is above MSL
        end
        
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
        
        %-------------- Describe the equation for b_f (m)
        F = abs(B_e - B_a) *1000*60*60*24*365;
        
    end

    function [G,Geq] = Fun_BoxModel_SS_con (x)
        % solves the system of equations at equilibrium
        
        b_f = x(1);
        d_f = x(2);
        d_m = x(3);
        
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
        
        %-------------- Describe the equation for d_f (m)
        Geq(1) = TF_erosion - TF_accretion - SOM + R;
        
        
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
        Geq(2) = bed_erosion + ocean_in + river_in - TF_deposition - M_deposition - export;   % (kg/s)
        
        Geq = Geq *1000*60*60*24*365;
        G = [];
        
    end

    function stop = outfun(x,optimValues,state)
        stop = false;
        
        switch state
            case 'init'
%                 hold on
            case 'iter'
                % Concatenate current point and objective function
                % value with history. x must be a row vector.
                history.fval = [history.fval; optimValues.fval];
                history.x = [history.x; x];
                % Concatenate current search direction with
                % searchdir.
                %                 searchdir = [searchdir;...
                %                     optimValues.searchdirection'];
                %                 plot(x(1),x(2),'o');
                %                 % Label points with iteration number and add title.
                %                 % Add .15 to x(1) to separate label from plotted 'o'
                %                 text(x(1)+.15,x(2),...
                %                     num2str(optimValues.iteration));
                %                 title('Sequence of Points Computed by fmincon');
            case 'done'
%                 hold off
            otherwise
        end
        
    end

end