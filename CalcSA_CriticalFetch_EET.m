function CalcSA_CriticalFetch_EET
% CalcSA_CriticalFetch_EET calculates sensitivity using EET method
%
% Last Update: 11/27/2017
%--------------------------------------------------------------------------------------------------
format compact
format longG

par_temp =1:8;

L = 4;

Co = linspace(1,100,L)*10^-3;
Cf = linspace(1,100,L)*10^-3;
Qf = linspace(1,100,L);
LE = linspace(1,10,L)*10^3;
bfm = linspace(1,10,L)*10^3;
a =  linspace(1,3,L);
R = linspace(1,20,L)/10^3/365/24/60/60;
vw = linspace(0,10,L);

for k = 1 : length(par_temp)
    
    par = par_temp(k);
    
    switch par
        case 1
            load co_data_SA_EET
            delta = Co(3)-Co(1);
            range = [Co(1), Co(L)];
        case 2
            load cf_data_SA_EET
            delta = Cf(3)-Cf(1);
            range = [Cf(1), Cf(L)];
       case 3
            load qf_data_SA_EET
            delta = Qf(3)-Qf(1);
            range = [Qf(1), Qf(L)];
        case 4
            load le_data_SA_EET
            delta = LE(3)-LE(1);
               range = [LE(1), LE(L)];
     case 5
            load bfm_data_SA_EET
            delta = bfm(3)-bfm(1);
            range = [bfm(1), bfm(L)];
        case 6
            load a_data_SA_EET
            delta = a(3)-a(1);
              range = [a(1), a(L)];
      case 7
            load R_data_SA_EET
            delta = R(2)-R(1);
              range = [R(1), R(L)];
      case 8
            load vw_data_SA_EET
            delta = vw(3)-vw(1);
            range = [vw(1), vw(L)];
    end
    
    ee = zeros(2,2);
    for j = 1 : 2
        ee(j,:)  = (dat(j+2,1)-dat(j,1))/delta*range;
    end
    
    stat(k,1:2) = mean(ee); % mu
    stat(k,3:4) = sum(abs(ee))/2; % mu*
    stat(k,5:6) = sqrt(mean((ee-mean(ee)).^2)); % sigma^2
    
end

save SA_results_EET stat

end
