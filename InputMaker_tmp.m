function [inmat,Co] = InputMaker_tmp
%
% Last Update: 1/5/2018
%--------------------------------------------------------------------------------------------------
format compact
format longG

% x=[ 0.01,  0.14,     100,   5000,   5000,    1.4,   6.0,     12,      6];
% Co=[x(1)];Cf=x(2);Qf=x(3);LE=x(4);bfm=x(5);a=x(6);R=x(7)*10^-3/365/24/60/60;T=x(8)*60*60;vw=x(9);

Co = 10  *10^-3; % 10
Cf = .05;%(0:5:100) *10^-3; % 5,100 (0:10:200) [5 50]
Qf = [(0:5:100)];%(0:5:100); % 10,100 (10:10:200) [5,50]
LE = 5  *1000;
bfm = 5 *1000;
a = 1.4;%(1:20); % 1.4
R = [6] *10^-3/365/24/60/60;
T = 12 *60*60;
vw = 6;

k1 = length(Co);
k2 = length(Cf);
k3 = length(Qf);
k4 = length(LE);
k5 = length(bfm);
k6 = length(a);
k7 = length(R);
k8 = length(T);
k9 = length(vw);

% m = k2*k3*k4*k5*k6*k7*k8*k9;
% 
% inmat = [Cf Qf LE bfm a R T vw 1];


m = k2*k3*k4*k5*k6*k7*k8*k9;
inmat = zeros(m,9);
n = 1;

for i9 =1 : k9
    for i2 = 1 : k2
        for i3 = 1 : k3
            for i4 = 1 : k4
                for i5 = 1 : k5
                    for i6 = 1 : k6
                        for i7 = 1 : k7
                            for i8 = 1 : k8
                                inmat(n,:) = [Cf(i2),Qf(i3),LE(i4),bfm(i5),a(i6),R(i7),T(i8),vw(i9),n];
                                n = n+1;
                            end
                        end
                    end
                end
            end
        end
    end
end

end