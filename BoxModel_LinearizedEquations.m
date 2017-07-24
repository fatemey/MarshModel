% function [t, y] = BoxModel_LinearizedEquations
function BoxModel_LinearizedEquations
% BoxModel_LinearizedEquations: Models linearized 0d marsh and tidal flat time
% evolution based on final linearized equations using Matlab ode23tb function.
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
% Last Update: 6/8/2017
%
%--------------------------------------------------------------------------------------------------
format compact
format longG
clear
% clf

%-------------- Set the time span
tyr = 1000;  % solve for time tyr (years)
ts = tyr *365*24*60*60; % tyr in (s)
dt = 12*60*60; % time step in (s)
tspan = 0:dt:ts;

%-------------- Sediment input constants
C_o = 20 *10^-3;    % ocean concertation (kg/m3)
C_f = 100 *10^-3;    % river concentration (kg/m3)
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
H = 1.4/2;          % tidal amplitude (range/2) (m)

%-------------- Sediment properties
rho_s = 1000;   % sediment bulk density (kg/m3)

%-------------- Model assumptions
Q_f = Q_f/2;    % consider half of the discharge only for one side of the tidal platform (the same will be automatically considered below for Q_T)
% b_fm = b_fm/2;  % consider half of the basin only for one side of the tidal platform

%-------------- Initial conditions, y0=[ b_f, d_f, d_m,u (=C_r*(b_f*d_f+b_m*d_m))]
y0(1) = 1000;%b_fm/2;      % tidal flat width (m)
y0(2) = H+0.3;        % tidal flat depth (m)
y0(3) = H-0.3;         % marsh depth (m)
y0(4) =C_o*(y0(1)*y0(2)+(b_fm-y0(1))*y0(3)); % u

%-------------- Solve the system of differential equations
[t, y] = ode15s(@ode4marshtidalflat_linearized,tspan,y0); % or use ode15s/ode23s/ode23tb
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

% fetch2depth = y(:,1)./y(:,2); % fetch to tidal flat depth ratio

%-------------- Plot Results
% figure(2)
clf
plot_BoxModel(t,y)
% tit = 'R_4-bf0_745';
% print(tit,'-dtiff','-r400')
% movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\15 - Model Parameters relationships')
% close all

% figure(2)
% clf
% scatter(t,fetch2depth,'.','k')
% xlabel('Year')
% ylabel('Fetch to Depth Ratio')
% box on
% set(findobj('type','axes'),'fontsize',15)

%======================= Nested Function =========================
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
