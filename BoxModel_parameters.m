function [a, b] = BoxModel_parameters(t, y, x_par, y_par)
% plot_BoxModel_vars: plots different parameters
%
% Input
%           t : vector of time data in yr
%           y : matrix of data. 1st vertical vector: tidal flat width (m),
%           2nd vector: tidal flat depth (m), 3rd vector: marsh depth (m),
%           and 4th vector: C_r (g/m3).
%           x_par : char vector of the x axis parameter of interest. It can
%                       be 't', 'd', 'h' or 'fetch'
%           y_par : char vector of the y axis parameter of interest. It can
%                       be 'tau', 'bed_e', 'mar_e' (for margin erosion rate),
%                       'margin' (for margin accretion or erosion rate), 'ocean' (for ocean inout),
%                       'ratio' (for ocean to river input ratio). if margin
%                       is positive it is the result of erosion, otherwise,
%                       it is due to accretion.
%
%--------------------------------------------------------------------------------------------------

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
gamma = 9800;   % water specific weight (N/m3)

%-------------- Model assumptions
Q_f = Q_f/2;    % consider half of the discharge only for one side of the tidal platform (the same will be automatically considered below for Q_T)
b_fm = b_fm/2;  % consider half of the basin only for one side of the tidal platform

switch x_par
    
    case 't'
        a = t;
        x_lab = 'Time (yr)';
        
    case 'd' % tidal flat depth
        a = y(:,2);
        x_lab = 'Tidal Flat Depth (m)';
        
        %-------------- Compute Fetch (m)
    case 'fetch'
        a = y(:,1)*2/1000;
        x_lab = 'Fetch (km)';
        
        %-------------- Compute Characteristic Depth (m)
    case 'h'
        a = (y(:,2)+max(0,y(:,2)-2*H))/2;
        x_lab = 'Reference Depth (m)';
        
end

switch y_par
    
    %-------------- Compute Shear Stress (PA)
    case 'tau'
        h = (y(:,2)+max(0,y(:,2)-2*H))/2;
        chi = y(:,1)*2;
        tau = zeros(size(h));
        for i=1:length(h)
            [ H_w, T_w] = WaveProps ( h(i), v_w, chi(i));   % compute significant height and peak period
            [ tau(i), k_w ] = ShearStress ( h(i), k_0, H_w, T_w);
        end
        %         load tau
        b = tau;
        y_lab = 'Shear Stress (PA)';
        
        %-------------- Compute Tidal Flat Bed Erosion (kg/s)
    case 'bed_e'
        h = (y(:,2)+max(0,y(:,2)-2*H))/2;
        chi = y(:,1)*2;
        tau = zeros(size(h));
        for i=1:length(h)
            [ H_w, T_w] = WaveProps ( h(i), v_w, chi(i));   % compute significant height and peak period
            [ tau(i), k_w ] = ShearStress ( h(i), k_0, H_w, T_w);
        end
        %         load tau
        t_s = ones(size(y,1),1);
        t_s(y(:,2) < 2*H) = 1/2-1/pi*asin((H-y(:,2))/H);
        
        b = zeros(size(y,1),1);
        b(y(:,2) > H) = max(0,t_s*E_0.*(tau-tau_c)./tau_c.*y(:,1)*L_E);
        y_lab = 'Bed Erosion Rate (kg/s)';
        
        %-------------- Compute Margin Erosion (kg/s)
    case 'mar_e'
        b_m = b_fm-y(:,1); % marsh width
        h = (y(:,2)+max(0,y(:,2)-2*H))/2;
        chi = y(:,1)*2;
        W = zeros(size(h));
        for i=1:length(h)
            [ H_w, T_w] = WaveProps ( h(i), v_w, chi(i));   % compute significant height and peak period
            [ tau, k_w ] = ShearStress ( h(i), k_0, H_w, T_w);
            c_g = pi/k_w/T_w*(1+2*k_w*h(i)/sinh(2*k_w*h(i))); % wave group velocity (general form)
            W(i) = gamma*c_g*H_w^2/16; % wave power density (kg.m/s3)
        end
        
        b = zeros(size(y,1),1);
        b(y(:,2) > H) = k_e*W;
        y_lab = 'Margin Erosion Rate (kg/s)';
        
        %-------------- Compute Margin Erosion/Accretion (kg/s)
    case 'margin'
        b_m = b_fm-y(:,1); % marsh width
        h = (y(:,2)+max(0,y(:,2)-2*H))/2;
        chi = y(:,1)*2;
        W = zeros(size(h));
        for i=1:length(h)
            [ H_w, T_w] = WaveProps ( h(i), v_w, chi(i));   % compute significant height and peak period
            [ tau, k_w ] = ShearStress ( h(i), k_0, H_w, T_w);
            c_g = pi/k_w/T_w*(1+2*k_w*h(i)/sinh(2*k_w*h(i))); % wave group velocity (general form)
            W(i) = gamma*c_g*H_w^2/16; % wave power density (kg.m/s3)
        end
        %         load W
        %         if  flag_f2m==1 || chi<=b_r || v_w==0 % condition for no bed and margin erosion in case of a filled mudflat or no wind
        %             tau = 0; % bed shear stress
        %             W = 0;   % wave power density
        %         else
        B_e = k_e*W;    % margin erosion
        C_r = y(:,4)./(y(:,1).*y(:,2)+b_m.*y(:,3)); % concentration (kg/m3)
        B_a = zeros(size(y,1),1);
        B_a(y(:,2) > H) = k_a*omega_s*C_r/rho_s;  % margin accretion
        
        b = zeros(size(y,1),1);
        b(y(:,2) > H) = (B_e - B_a).*(y(:,2)-y(:,3))*L_E*rho_s;
        y_lab = 'Margin Erosion/Accretion Rate (kg/s)';
        
        %-------------- Compute Ocean Sediment Input (kg/s)
    case 'ocean'
        b_m = b_fm-y(:,1); % marsh width
        Q_T = (y(:,2).*y(:,1)+y(:,3).*b_m)*L_E/T_T-Q_f;
        
        b = Q_T*C_o;
        y_lab = 'Ocean Sediment Input (kg/s)';
        
        %-------------- Compute Ocean Input to River Input(kg/s)
    case 'ratio'
        b_m = b_fm-y(:,1); % marsh width
        Q_T = (y(:,2).*y(:,1)+y(:,3).*b_m)*L_E/T_T-Q_f;
        ocean_in = Q_T*C_o;
        river_in = Q_f*C_f;
        
        b = ocean_in/river_in;
        y_lab = 'Ocean to River Input Ratio (kg/s)';
end

%-------------- Make a Plot
% figure
% plot_BoxModel(t,y)
figure
scatter(a,b,'k','.')
xlabel(x_lab)
ylabel(y_lab)
set(findobj('type','axes'),'fontsize',15)
box on
h_fig=gcf;
set(h_fig,'PaperOrientation','portrait')
set(h_fig,'PaperPosition', [0 0 7.5 6]) % [... ... max_width=7.5 max_height=9]
tit=[x_par, '-', y_par];
print(tit,'-dtiff','-r400')
movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\15 - Model Parameters relationships')
close all