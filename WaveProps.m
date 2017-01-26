function [ H_w, T_w ] = WaveProps ( h, v_w, chi, g )
% WaveProps : Computing significant wave height and peak wave period
%
% Input 
%       h   : reference water depth
%       v_w : wind velocity
%       chi : fetch
%
% Output 
%       H_w : significant wave height
%       T_w : peak wave period
%
%--------------------------------------------------------------------------

%-------------- Parameters
A1 = .493*(g*h/v_w^2)^.75;
B1 = 3.13*10^-3*(g*chi/v_w^2)^.57;
A2 = .331*(g*h/v_w^2)^1.01;
B2 = 5.215*10^-4*(g*chi/v_w^2)^.73;

%-------------- Outpouts
H_w = .2413*(tanh(A1)*tanh(B1/tanh(A1)))^.87*v_w^2/g;   % significant wave height
T_w = 7.518*(tanh(A2)*tanh(B2/tanh(A2)))^.37*v_w/g;     % peak wave period

end

