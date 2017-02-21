function [ tau, k_w ] = ShearStress ( h, k_0, H_w, T_w )
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
g = 9.81;       % gravitational acceleration (m/s2)

%-------------- Parameters
% k_w_0 = 2*pi/T_w/sqrt(g*h);
k_w_0 = 1;
opts = optimoptions('fsolve','Display','off');
% opts = optimoptions('fsolve','Display','off','Algorithm','levenberg-marquardt');
% k_w = fsolve(@ (k_w) (2*pi/T_w)^2-g*k_w*tanh(k_w*h),k_w_0,opts); % wave number (general form)
k_w = 2*pi/T_w/sqrt(g*h);               % wave number (for k_w*h<0.1pi)
f_w = 0.4*(H_w/k_0/sinh(k_w*h))^-.75;   % friction factor
u_w = pi*H_w/T_w/sinh(k_w*h);           % horizontal orbital velocity

%-------------- Outpouts
tau = rho_w*f_w*u_w^2/2;    % wave shear stress

end

