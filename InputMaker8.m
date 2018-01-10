function [input,Co,k1,m] = InputMaker8
%
% Last Update: 1/5/2018
%--------------------------------------------------------------------------------------------------
format compact
format longG

Co = ([5,15,25]) *10^-3;
Cf = ([10,50,90]) *10^-3;
Qf = [10,50,90];
LE = ([1,4,9]) * 1000;
bfm = ([1,4,9]) *1000;
a = ([1,2,3]);
R = ([6]) *10^-3/365/24/60/60;
T = 12 *60*60;
vw = 10;%[2,6,10];

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
input = zeros(m,8);
n = 1;

for i9 =1 : k9
    for i2 = 1 : k2
        for i3 = 1 : k3
            for i4 = 1 : k4
                for i5 = 1 : k5
                    for i6 = 1 : k6
                        for i7 = 1 : k7
                            for i8 = 1 : k8
                                input(n,:) = [Cf(i2),Qf(i3),LE(i4),bfm(i5),a(i6),R(i7),T(i8),vw(i9)];
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