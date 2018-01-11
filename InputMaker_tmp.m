function [input,Co,k1,m] = InputMaker_tmp
%
% Last Update: 1/5/2018
%--------------------------------------------------------------------------------------------------
format compact
format longG

x=[0.015,  0.010,     90,   4000,   1000,    3.0,   10.0,     12,      2];
Co=[x(1)];Cf=x(2);Qf=x(3);LE=x(4);bfm=x(5);a=x(6);R=x(7)*10^-3/365/24/60/60;T=x(8)*60*60;vw=x(9);


k1 = length(Co);
k2 = length(Cf);
k3 = length(Qf);
k4 = length(LE);
k5 = length(bfm);
k6 = length(a);
k7 = length(R);
k8 = length(T);
k9 = length(vw);

m = k2*k3*k4*k5*k6*k7*k8*k9;

input = [Cf Qf LE bfm a R T vw 1];

end