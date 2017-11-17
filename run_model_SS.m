function run_model_SS

% Co = (10:10:40) *10^-3;
% Cf = (10:10:40) *10^-3;
% LE = (4:8) * 1000;
% Qf = 10:10:40;
% R = (2:6) *10^-3/365/24/60/60;
% T = [12,24] *60*60;
% H = [1.4,2.8] /2;
% bfm =  (5) *1000;
% vw = (2:4:10);

Co = (20) *10^-3;
Cf = (15) *10^-3;
LE = (5) *1000;
Qf = 20;
R = (1:.5:10) *10^-3/365/24/60/60;
T = ([12]) *60*60;
H = (1.4) /2;
bfm =  (5) *1000;
vw = (6);

k1 = length(Co);
k2 = length(Cf);
k3 = length(LE);
k4 = length(Qf);
k5 = length(R);
k6 = length(H);
k7 = length(T);
k8 = length(bfm);
k9 = length(vw);

dat = zeros(k1*k2*k3*k4*k5*k6*k7*k8*k9,11);
n = 1;
for i1 =1 : k1
    for i2 = 1 : k2
        for i3 = 1 : k3
            for i4 = 1 : k4
                for i5 = 1 : k5
                    for i6 = 1 : k6
                        for i7 = 1 : k7
                            for i8 = 1 : k8
                                for i9 = 1 : k9
                                    Sol = BoxModel_SS_2eq(Co(i1),Cf(i2),LE(i3),Qf(i4),R(i5),H(i6),T(i7),bfm(i8),vw(i9)); % Sol = [x,y,fval_x,fval_y]
                                    dat(n,1) = min(Sol(1),bfm(i8)); % b_f
                                    dat(n,2:3) = Sol(2:3); % d_f and d_m
                                    dat(n,4:12) = [Co(i1),Cf(i2),LE(i3),Qf(i4),R(i5),H(i6),T(i7),bfm(i8),vw(i9)];
                                    dat(n,13:15) = Sol(4:6); % fval_bf,fval_df, fval_d in: mm/yr
                                    n = n + 1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

save plotdatar dat

end


