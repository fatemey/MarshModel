function run_model_SS

Co = (10:10:40) *10^-3;
Cf = (10:10:40) *10^-3;
LE = (.6:.2:1.4) * 1000;
Qf = 10:10:40;
R_ = (2:6) *10^-3/365/24/60/60;

k1 = length(Co);
k2 = length(Cf);
k3 = length(LE);
k4 = length(Qf);
k5 = length(R_);
dat = zeros(k1*k2*k3*k4*k5,8);
n = 1;
for i1 =1 : k1
    for i2 = 1 : k2
        for i3 = 1 : k3
            for i4 = 1 : k4
                for i5 = 1 : k5
                    dat(n,1:3) = BoxModel_SS_2eq_org(Co(i1),Cf(i2),LE(i3),Qf(i4),R_(i5));
                    dat(n,4:8) = [Co(i1),Cf(i2),LE(i3),Qf(i4),R_(i5)];
                    n = n + 1;
                end
            end
        end
    end
end

save SS_results_small dat

end


