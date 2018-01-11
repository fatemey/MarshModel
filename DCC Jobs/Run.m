function Run(ind1,ind2)
% Runs CriticalFetchTM and CriticalFetchSS for a method using both steady
% state and time marching approaches.
%
% Last Update: 1/10/2018
%--------------------------------------------------------------------------------------------------
format compact
format longG

width_inc = 10; % tidal flat width increment for time marching approach
bf0SS = 100;
[inmat,Co] = Input;
input = inmat(ind1:ind2,:);
k = length(Co);
m = ind2-ind1+1;

parpool()

fprintf('%6s, %6s, %7s, %6s, %6s, %7s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %4s, %4s\n','n','iCo','bf','df','dm','bfss','dfss','dmss','Co','Cf','Qf','LE','bfm','a','R','T','VW','Yr','MM','DD','Hr','Min','Sec');

parfor i = 1 : m
    
    flag = 0;
    dat_temp = zeros(k,6);
    
    for i1 = 1 : k
        
        SolSS = CriticalFetchSS(Co(i1),input(i,1),input(i,2),input(i,3),input(i,4),input(i,5),input(i,6),input(i,7),input(i,8),bf0SS); % Sol = [x,y,fval_x,fval_y]
        err_SS = max(abs(SolSS(4:6)));
        
        if err_SS < 10^-3 % if SS solution is relible
            dat_temp(i1,4) = min(SolSS(1), input(i,4)); % b_f
            dat_temp(i1,5:6) = SolSS(2:3); % d_f and d_m
        else
            dat_temp(i1,4:6) = [-1,-1,-1]; % showing that SS was not relible
        end
        
        if flag == 0 % if this is the first run for Co loop
            
            if err_SS >= 10^-3 || SolSS(3) > input(i,5)/2 % if the solution is not relible or marsh is drowned
                bf0TM = 1;
                SolTM = CriticalFetchTM(Co(i1),input(i,1),input(i,2),input(i,3),input(i,4),input(i,5),input(i,6),input(i,7),input(i,8),bf0TM,width_inc);
                dat_temp(i1,1:3) = SolTM(1:3); % b_f, d_f and d_m
                
            elseif err_SS < 10^-3 && SolSS(1) >= input(i,4) % if the solution is relible and b_f hits the upper boundary
                dat_temp(i1,1) = input(i,4); % b_f
                dat_temp(i1,2:3) = SolSS(2:3); % d_f and d_m
                
            else % if the solution is relible
                bf0TM = max(SolSS(1) - 5*width_inc, 1);
                SolTM = CriticalFetchTM(Co(i1),input(i,1),input(i,2),input(i,3),input(i,4),input(i,5),input(i,6),input(i,7),input(i,8),bf0TM,width_inc);
                dat_temp(i1,1:3) = SolTM(1:3); % b_f, d_f and d_m
                
            end
            
        else % if this is not the first run for Co loop
            
            if input(i,4) - dat_temp(i1-1,1) <= width_inc % if the width was close enough to the upper limit in the prevoius iteration
                dat_temp(i1,1:3) = dat_temp(i1-1,1:3);
                
            elseif err_SS < 10^-3 && SolSS(3) <= input(i,5)/2 % if the solution is relible and marsh is not drowned
                bf0TM = max(max(dat_temp(i1-1,1) - 1, 1), min(SolSS(1), input(i,4) - 1));
                SolTM = CriticalFetchTM(Co(i1),input(i,1),input(i,2),input(i,3),input(i,4),input(i,5),input(i,6),input(i,7),input(i,8),bf0TM,width_inc);
                dat_temp(i1,1:3) = SolTM(1:3); % b_f, d_f and d_m
                
            else
                bf0TM = max(dat_temp(i1-1,1) - 1, 1);
                SolTM = CriticalFetchTM(Co(i1),input(i,1),input(i,2),input(i,3),input(i,4),input(i,5),input(i,6),input(i,7),input(i,8),bf0TM,width_inc);
                dat_temp(i1,1:3) = SolTM(1:3); % b_f, d_f and d_m
                
            end
            
        end

        flag = 1;
        input_data = [Co(i1),input(i,1),input(i,2),input(i,3),input(i,4),input(i,5),input(i,6)*10^3*365*24*60*60,input(i,7)/60/60,input(i,8)];
        t = datetime;
        time = datevec(t);
        
        fprintf('%6d, %6d, %7.2f, %6.2f, %6.2f, %7.2f, %6.2f, %6.2f, %6.3f, %6.3f, %6d, %6d, %6d, %6.1f, %6.1f, %6d, %6d, %6d, %6d, %6d, %6d, %4d, %4.0f\n',[input(i,9),i1,dat_temp(i1,:),input_data,time]);
        
    end
    
end

end