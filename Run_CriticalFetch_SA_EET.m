function Run_CriticalFetch_SA_EET
% Runs CriticalFetch for sensitivity analysis using EET method
%
% Last Update: 11/27/2017
%--------------------------------------------------------------------------------------------------
format compact
format longG

par_temp =1:8;

L = 4;

Co = linspace(1*10^-3,100 *10^-3,L);
Cf = linspace(1*10^-3,100 *10^-3,L);
Qf = linspace(1,100 ,L);
LE = linspace(1*10^3,10 *10^3,L);
bfm = linspace(1*10^3,10 *10^3,L);
a =  linspace(1,3 ,L);
R = linspace(1/10^3/365/24/60/60,20/10^3/365/24/60/60,L);
vw = linspace(0,10 ,L);

for k = 1 : length(par_temp)
    
    k
    
    Co_bar =ones(1,L)*50*10^-3;
    Cf_bar = ones(1,L)*50 *10^-3;
    Qf_bar =ones(1,L)*50 ;
    LE_bar = ones(1,L)*5 *10^3;
    bfm_bar = ones(1,L)*5 *10^3;
    a_bar =  ones(1,L)*2;
    R_bar = ones(1,L)*10/10^3/365/24/60/60;
    vw_bar = ones(1,L)*5;
    
    par = par_temp(k);
    
    switch par
        case 1
            Co_bar = Co;
        case 2
            Cf_bar = Cf;
        case 3
            Qf_bar = Qf;
        case 4
            LE_bar = LE;
        case 5
            bfm_bar = bfm;
        case 6
            a_bar = a;
        case 7
            R_bar = R;
        case 8
            vw_bar = vw;
    end
    
    dat = zeros(L,11);
    
    for j = 1 : L
        
        j
        
        SolSS = CriticaIFetchSS(Co_bar(j),Cf_bar(j),Qf_bar(j),LE_bar(j),bfm_bar(j),a_bar(j),R_bar(j),vw_bar(j),100); % Sol = [x,y,fval_x,fval_y]
        bf0 = min(SolSS(1),bfm_bar(j)-1); % b_f ; condition for when b_f hits the upper boundary
        if SolSS(3)>a_bar(j) || max(abs(SolSS(4:6)))>10^-3 % condition for when b_f hits the lower boundary OR SS solution is not relible
            bf0 = 1; % b_f ; condition for when b_f hits the lower boundary
        end
        
        SolTM = CriticaIFetchTM(Co_bar(j),Cf_bar(j),Qf_bar(j),LE_bar(j),bfm_bar(j),a_bar(j),R_bar(j),vw_bar(j),bf0,10);
        
        dat(j,1) = max(SolSS(1),SolTM(1)); % b_f
        dat(j,2:3) = SolTM(2:3); % d_f and d_m
        dat(j,4:11) = [Co_bar(j),Cf_bar(j),Qf_bar(j),LE_bar(j),bfm_bar(j),a_bar(j),R_bar(j) *10^3*365*24*60*60,vw_bar(j)];
        
    end
    
    switch par
        case 1
            save('co_data_SA_EET.mat','dat')
        case 2
            save('cf_data_SA_EET.mat','dat')
        case 3
            save('qf_data_SA_EET.mat','dat')
        case 4
            save('le_data_SA_EET.mat','dat')
        case 5
            save('bfm_data_SA_EET.mat','dat')
        case 6
            save('a_data_SA_EET.mat','dat')
        case 7
            save('R_data_SA_EET.mat','dat')
        case 8
            save('vw_data_SA_EET.mat','dat')
    end
    
end

end