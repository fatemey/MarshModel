function Run(ind1,ind2)
% Runs CriticalFetchTM and CriticalFetchSS for a method using both steady
% state and time marching approaches.
%
% Last Update: 1/16/2018
%--------------------------------------------------------------------------------------------------
format compact
format longG

width_inc = 10; % tidal flat width increment for time marching approach
bf0SS_small = 100;
bf0SS_large = 1000;
input = Input4graphs;
inmat = input(ind1:ind2,:);
m = ind2-ind1+1;

parpool()

fprintf('%6s, %7s, %6s, %6s, %7s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %4s, %4s\n','n','bf','df','dm','bfss','dfss','dmss','Co','Cf','Qf','LE','bfm','a','R','T','VW','Yr','MM','DD','Hr','Min','Sec');

parfor i = 1 : m
    
    outmat = zeros(1,6);
    Co=inmat(i,1);Cf=inmat(i,2);Qf=inmat(i,3);LE=inmat(i,4);bfm=inmat(i,5);a=inmat(i,6);R=inmat(i,7);T=inmat(i,8);vw=inmat(i,9);
    
    bf0SS = bf0SS_small;
    SolSS = CriticalFetchSS(Co,Cf,Qf,LE,bfm,a,R,T,vw,bf0SS); % Sol = [x,y,fval_x,fval_y]
    err_SS = max(abs(SolSS(4:6)));
    
    if err_SS < 10^-3 &&  SolSS(3) > a/2 % try a larger initial condition for width if needed
        bf0SS = bf0SS_large;
        SolSS_temp = CriticalFetchSS(Co,Cf,Qf,LE,bfm,a,R,T,vw,bf0SS); % Sol = [x,y,fval_x,fval_y]
        err_SS_temp = max(abs(SolSS_temp(4:6)));
        if err_SS_temp < 10^-3 &&  SolSS_temp(3) <= a/2
            SolSS = SolSS_temp;
            err_SS = err_SS_temp;
        end
    end
    
    if err_SS < 10^-3 && SolSS(3) <= a/2 && SolSS(2) > a/2 % set the initial depth conditions
        df0 = SolSS(2);
        dm0 = SolSS(3);
    else
        df0 = 3*a/4;
        dm0 = a/4;
    end
    
    if err_SS < 10^-3 % if SS solution is relible; setting the output
        outmat(4) = min(SolSS(1), bfm); % b_f
        outmat(5:6) = SolSS(2:3); % d_f and d_m
    else
        outmat(4:6) = [-1,-1,-1]; % showing that SS was not relible
    end
    
    if err_SS >= 10^-3 || SolSS(3) > a/2 % if the solution is not relible or marsh is drowned
        bf0TM = 1;
        SolTM = CriticalFetchTM(Co,Cf,Qf,LE,bfm,a,R,T,vw,bf0TM,df0,dm0,width_inc);
        outmat(1:3) = SolTM(1:3); % b_f, d_f and d_m
        
    elseif err_SS < 10^-3 && SolSS(1) >= bfm % if the solution is relible and b_f hits the upper boundary
        outmat(1) = bfm; % b_f
        outmat(2:3) = SolSS(2:3); % d_f and d_m
        
    else % if the solution is relible
        bf0TM = max(SolSS(1) - 5*width_inc, 1);
        SolTM = CriticalFetchTM(Co,Cf,Qf,LE,bfm,a,R,T,vw,bf0TM,df0,dm0,width_inc);
        outmat(1:3) = SolTM(1:3); % b_f, d_f and d_m
        
    end
    
    if outmat(1) == 0
        outmat(4) = 0;
    end
    
    input_data = [Co,Cf,Qf,LE,bfm,a,R*10^3*365*24*60*60,T/60/60,vw];
    t = datetime;
    time = datevec(t);
    
    fprintf('%6d, %7.2f, %6.2f, %6.2f, %7.2f, %6.2f, %6.2f, %6.3f, %6.3f, %6d, %6d, %6d, %6.1f, %6.1f, %6d, %6d, %6d, %6d, %6d, %6d, %4d, %4.0f\n',[inmat(i,10),outmat,input_data,time]);
    
end

end