function x = CriticalInitialWidth_LinearizedEquations
% Function CriticalInitialWidth_LinearizedEquations looks for a critical fetch value based on
% different values of a variable of interest. Functtion fetch_threshold is
% based on the function of BoxModel.
% after each run for the parameter of interest, save the x in the excle
% file, or save it as a mat file to plot later using the function plot_fetch.
%
% Last Update: 6/8/2017
%
%--------------------------------------------------------------------------------------------------
format compact
format longG
clear
% clf

par_temp = 1 : 8;

for k = 2 : 2
    
    k
    %-------------- Set the time span
    tyr = 500;  % solve for time tyr (years)
    ts = tyr *365*24*60*60; % tyr in (s)
    dt = 12*60*60; % time step in (s)
    tspan = 0:dt:ts;
    
    %-------------- Sediment input constants
    C_o = 20 *10^-3;    % ocean concertation (kg/m3)
    C_f = 15 *10^-3;    % river concentration (kg/m3)
    Q_f = 20;         % river water discharge (m3/s)
        
    %-------------- Vegetation properties
    B_max = 1;      % maximum biomass density (kg/m2)
    k_B = 2*10^-3 /365/24/60/60;    % vegetation characteristics (m3/s/kg)
    
    %-------------- Basin properties
    b_fm = 10 *10^3; % total basin width (both sides of the channel) (m)
    L_E = 10 *10^3; % basin length (m)
    R = 3 *10^-3/365/24/60/60;   % sea level rise (m/s)
    
    %-------------- Tide Characteristics
    T_T = 12 *60*60;   % tidal period (s) (= 12 hours)
    H = 1.4 /2;          % tidal amplitude (range/2) (m)
    
    %-------------- Sediment properties
    rho_s = 1000;   % sediment bulk density (kg/m3)
        
    %-------------- Model assumptions
    Q_f = Q_f/2;    % consider half of the discharge only for one side of the tidal platform (the same will be automatically considered below for Q_T)
    %  b_fm = b_fm/2;  % consider half of the basin only for one side of the tidal platform
    
    %-------------- Start the loop for each run
    par = par_temp(k);
    switch par
        case 1
            par_v = 5 *10^-3 : 5 *10^-3 : 100 *10^-3; % for C_o
            TF_width = 10 : 10 : b_fm-10;
        case 2
            par_v = 100 *10^-3: 100 *10^-3 : 300 *10^-3; % for C_f
%             par_v = 600*10^-3;
            TF_width = 100 : 100 : b_fm-100;
        case 3
            par_v = 0 : 100  : 1000; % for Q_f
            TF_width = 5 : 5 : b_fm-5;
        case 4
            par_v = 1 *10^3 : 2 *10^3 : 20 *10^3; % for L_E
            TF_width = 1 : 1 : b_fm-1;
        case 5
            par_v = 1 *10^3 : 2 *10^3 : 20 *10^3; % for b_fm
        case 6
            par_v = 0 : 2 : 20; % for v_w
            TF_width = 10 : 10 : b_fm-10;
        case 7
            par_v = 0 : 2 *10^-3/365/24/60/60 : 30 *10^-3/365/24/60/60; % for R
            TF_width = 5 : 5 : b_fm-5;
        case 8
            par_v = [1 : 1 : 10]/2; % for H
            TF_width = 5 : 5 : b_fm-5;
    end
    
    for j = 1 : length(par_v)
        
        j
        switch par
            case 1
                C_o = par_v(j);
            case 2
                C_f = par_v(j);
            case 3
                Q_f = par_v(j)/2;
            case 4
                L_E = par_v(j);
            case 5
                % b_fm = par_v(j)/2;
                b_fm = par_v(j);
                TF_width = 1 : 1 : b_fm-1;
            case 6
                v_w = par_v(j);
            case 7
                R = par_v(j);
            case 8
                H = par_v(j);
        end
        
        clear y t width_diff
        depthequil = zeros(length(par_v),2);
        
        for i = 1 : length(TF_width)
            
            i
            %-------------- Initial conditions, y0=[ b_f, d_f, d_m,u(=C_r*(b_f*d_f+b_m*d_m))]
            y0(1) = TF_width(i);
            y0(2) = H+0.3;         % tidal flat depth (m)
            y0(3) = H-0.3;         % marsh depth (m)
            y0(4) =C_o*(y0(1)*y0(2)+(b_fm-y0(1))*y0(3)); % u
            
            %-------------- Solve the system of differential equations
            [t, y] = ode15s(@ode4marshtidalflat_Linearized,tspan,y0); % or use ode15s/ode23s/ode23tb
            t = t /365/24/60/60; % convert time unit from s to yr for plotting purposes
            y(:,4) = y(:,4)./(y(:,1).*y(:,2)+y(:,3).*(b_fm-y(:,1))); % convert y(:,4) to C_r from the formula used before: y4=u (=C_r*(b_f*d_f+b_m*d_m)
            
            %-------------- Data removal
            ind = find(y(:,3)>H); % remove data related to marsh conversion to tidal flat
            if ~isempty(ind) && length(ind)>1
                y(ind(2):end,:)=[]; % retain only one value afetr conversion to remember it is a new tidal flat now
                t(ind(2):end,:)=[];
            end
            
            ind = find(y(:,2)<=H); % remove data related to tidal flat conversion to marsh
            if ~isempty(ind) && length(ind)>1
                y(ind(2):end,:)=[]; % retain only one value afetr conversion to remember it is a new marsh now
                t(ind(2):end,:)=[];
            end
            
            %-------------- Check tidal flat contraction/expansion
            width = y(:,1); % tidal falt width solution
            n = length(width);
            
            width_diff(i) = sign(width(n)-width(floor(n/2))); % -1 corresponds to TF contraction (or fully marsh) and +1 accounts for expansion (or fully tidal flat)
            
            %-------------- Check emergence, drowning and reaching to boundary limits
            f = 0; % flag for emergance and drowning
            if y(end,2) <= H % check whether if tidal flat has turned into a marsh
                width_diff(i) = -1;
                f = 1;
            end
            
            if width(n) >= b_fm % check whether if tidal flat has reached to its boundary limits
                width_diff(i) = 1;
            elseif width(n) <= 0
                width_diff(i) = -1;
            end
            
            if y(end,3) > H % check whether if marsh has drowned
                width_diff(i) = 1;
                f = 1;
            end
                        
            %-------------- Record the critical initial width
            n_width_diff = length(unique(width_diff));
            tidalflatd = y(floor(4*n/5):n,2);
            marshd = y(floor(4*n/5):n,3);
            
            if n_width_diff == 2
                x(j,1) = width(1);
                
                if f == 0
                    if max(diff(tidalflatd)) < 0.001
                        depthequil(j,1) = 1;
                    end
                    if max(diff(marshd)) < 0.001
                        depthequil(j,2) = 1;
                    end
                end
                
                break
                
            elseif i == length(TF_width) && unique(width_diff) == 1
                x(j,1) = TF_width(1);
                                
                if f == 0
                    if max(diff(tidalflatd)) < 0.001
                        depthequil(j,1) = 1;
                    end
                    if max(diff(marshd)) < 0.001
                        depthequil(j,2) = 1;
                    end
                end
                
                break
            elseif i == length(TF_width) && unique(width_diff) == -1
                x(j,1) = TF_width(end);
                                
                if f == 0
                    if max(diff(tidalflatd)) < 0.001
                        depthequil(j,1) = 1;
                    end
                    if max(diff(marshd)) < 0.001
                        depthequil(j,2) = 1;
                    end
                end
                
                break
            end
            
        end
        
    end
    
    %-------------- Save the results
    x_new = x;
    %     if length(par_v)>length(x)
    %         x_new(length(x)+1:length(par_v)) = -1;
    %     end
    
    dat = [par_v', x_new];
    
    save('cf_l.mat','dat','depthequil')
    
%     switch par
%         case 1
%             save('co_1.mat','dat','depthequil')
%         case 2
%             save('cf_1.mat','dat','depthequil')
%         case 3
%             save('qf_1.mat','dat','depthequil')
%         case 4
%             save('le_1.mat','dat','depthequil')
%         case 5
%             save('bfm_1.mat','dat','depthequil')
%         case 6
%             save('vw_1.mat','dat','depthequil')
%         case 7
%             save('R_1.mat','dat','depthequil')
%         case 8
%             save('H_1.mat','dat','depthequil')
%     end
    
    clear x dat x_new 
    
end

%======================= Nested Function =========================
    function dy = ode4marshtidalflat_Linearized (t,y) %  y1=b_f, y2=d_f, y3=d_m, y4=u (=C_r*(b_f*d_f+b_m*d_m, why solving u instead of C_r? u is the variable on the left hand side of mass conservation equation.)
        % solves the ODE system of equations
        
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
        
        %---------------------------------- Define the equations -----------------------------------
        C_r = y(4)/(local*y(2)+b_m*y(3)); % concentration (kg/m3)
        
        %---------------------------- Tidal flat width changes equation ----------------------------
        dy(1,1) = -7.*10^-8 + 3.*10^-12 *local + 1.4*10^-7 *y(2) + 2.*10^-9 *y(3);
        
        if (y(1) < 0 && dy(1,1) < 0) % imposing a constraint for boundary limits of b_f (or y(1))
            dy(1,1) = 0 ;
        end
        if (y(1) > b_fm && dy(1,1) > 0)
            dy(1,1) = 0 ;
        end
        
        %----------------------------- Tidal flat depth changes equation ---------------------------
        dy(2,1) = 8.*10^-9 + 1.*10^-12 *local - 6.*10^-9 *y(2) + 5.4*10^-10 *y(3);
        
        
        %----------------------------- Marsh depth changes equation ------------------------------
        M_accretion = C_r*y(3)/T_T/rho_s;
        
        %-------------- Compute organice matter production (m/s)
        z = H-y(3);         % elevation of marsh platform
        if z >= 0           % condition for presence of vegetation when marsh is above MSL
            r = -0.5*z/H+1; % reproduction rate
            m = 0.5*z/H;    % mortality rate
            B = B_max*(1-m/r);  % steady state solution for biomass (kg/m2)
            O = k_B*B;      % organic matter production rate
        else
            O = 0;
        end
        
        %-------------- Describe the equation for d_m (m)
        dy(3,1) = - M_accretion - O + R;
        
        %-------------------------------- Mass conservation equation ------------------------------
        var = 1.*10^-9 *local^2* L_E + b_fm *(7.*10^-8 + 0.00002 *C_o + 4.6*10^-8 *y(3))*y(3)* L_E + ...
            local* (8.*10^-6 + y(2)* (-6.*10^-6 + 0.00002 *C_o + 4.6*10^-8 *y(3)) + ...
            6.*10^-7 *y(3) - 0.00002 *C_o *y(3) - 4.6*10^-8 *y(3)^2)* L_E + C_f *Q_f - C_o *Q_f;
        dy(4,1) = var / L_E;
        
    end

end
