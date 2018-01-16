function Run_CriticalFetch
% Runs CriticalFetchTM and CriticalFetchSS for a method using both steady
% state and time marching approaches.
%
% Last Update: 1/9/2018
%--------------------------------------------------------------------------------------------------
format compact
format longG

width_inc = 10; % tidal flat width increment for time marching approach
bf0SS_small = 100;
bf0SS_large = 1000;
[inmat,Co] = InputMaker_tmp;
% load data_cf_2 
% inmat(inmat(:,6)>2e-10,:)=[];
% Co=.01;
k = length(Co);
m = size(inmat,1);
output = zeros(m,15*k);

% parpool()

fprintf('%6s, %6s, %7s, %6s, %6s, %7s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %4s, %4s\n','n','iCo','bf','df','dm','bfss','dfss','dmss','Co','Cf','Qf','LE','bfm','a','R','T','VW','Yr','MM','DD','Hr','Min','Sec');

for i = 1 : m
    
    flag = 0;
    dat_temp = zeros(k,6);
    
    for i1 = 1 : k
        
        bf0SS = bf0SS_small;
        SolSS1 = CriticalFetchSS(Co(i1),inmat(i,1),inmat(i,2),inmat(i,3),inmat(i,4),inmat(i,5),inmat(i,6),inmat(i,7),inmat(i,8),bf0SS); % Sol = [x,y,fval_x,fval_y]
        err_SS1 = max(abs(SolSS1(4:6)));
        
        if err_SS1 < 10^-3 &&  SolSS1(3) > inmat(i,5)/2 % try a larger initial condition for width
            bf0SS = bf0SS_large;
            SolSS2 = CriticalFetchSS(Co(i1),inmat(i,1),inmat(i,2),inmat(i,3),inmat(i,4),inmat(i,5),inmat(i,6),inmat(i,7),inmat(i,8),bf0SS); % Sol = [x,y,fval_x,fval_y]
            err_SS2 = max(abs(SolSS2(4:6)));
            if err_SS2 < 10^-3 &&  SolSS2(3) <= inmat(i,5)/2
                SolSS = SolSS2;
                err_SS = err_SS2;
            else
                SolSS = SolSS1;
                err_SS = err_SS1;
            end
        else
            SolSS = SolSS1;
            err_SS = err_SS1;
        end
        
        if err_SS < 10^-3 && SolSS(3) <= inmat(i,5)/2 && SolSS(2) > inmat(i,5)/2
            df0 = SolSS(2);
            dm0 = SolSS(3);
        else
            df0 = inmat(i,5)/4*3;
            dm0 = inmat(i,5)/4;
        end
        
        if err_SS < 10^-3 % if SS solution is relible
            dat_temp(i1,4) = min(SolSS(1), inmat(i,4)); % b_f
            dat_temp(i1,5:6) = SolSS(2:3); % d_f and d_m
        else
            dat_temp(i1,4:6) = [-1,-1,-1]; % showing that SS was not relible
        end
        
        if flag == 0 % if this is the first run for Co loop
            
            if err_SS >= 10^-3 || SolSS(3) > inmat(i,5)/2 % if the solution is not relible or marsh is drowned
                bf0TM = 1;
                SolTM = CriticalFetchTM_t(Co(i1),inmat(i,1),inmat(i,2),inmat(i,3),inmat(i,4),inmat(i,5),inmat(i,6),inmat(i,7),inmat(i,8),bf0TM,df0,dm0,width_inc);
                dat_temp(i1,1:3) = SolTM(1:3); % b_f, d_f and d_m
                
            elseif err_SS < 10^-3 && SolSS(1) >= inmat(i,4) % if the solution is relible and b_f hits the upper boundary
                dat_temp(i1,1) = inmat(i,4); % b_f
                dat_temp(i1,2:3) = SolSS(2:3); % d_f and d_m
                
            else % if the solution is relible
                bf0TM = max(SolSS(1) - 5*width_inc, 1);
                SolTM = CriticalFetchTM_t(Co(i1),inmat(i,1),inmat(i,2),inmat(i,3),inmat(i,4),inmat(i,5),inmat(i,6),inmat(i,7),inmat(i,8),bf0TM,df0,dm0,width_inc);
                dat_temp(i1,1:3) = SolTM(1:3); % b_f, d_f and d_m
                
            end
            
        else % if this is not the first run for Co loop
            
            if inmat(i,4) - dat_temp(i1-1,1) <= width_inc % if the width was close enough to the upper limit in the prevoius iteration
                dat_temp(i1,1:3) = dat_temp(i1-1,1:3);
                
            elseif err_SS < 10^-3 && SolSS(3) <= inmat(i,5)/2 % if the solution is relible and marsh is not drowned
                bf0TM = max(max(dat_temp(i1-1,1) - 1, 1), min(SolSS(1), inmat(i,4) - 1));
                SolTM = CriticalFetchTM_t(Co(i1),inmat(i,1),inmat(i,2),inmat(i,3),inmat(i,4),inmat(i,5),inmat(i,6),inmat(i,7),inmat(i,8),bf0TM,df0,dm0,width_inc);
                dat_temp(i1,1:3) = SolTM(1:3); % b_f, d_f and d_m
                
            else
                bf0TM = max(dat_temp(i1-1,1) - 1, 1);
                SolTM = CriticalFetchTM_t(Co(i1),inmat(i,1),inmat(i,2),inmat(i,3),inmat(i,4),inmat(i,5),inmat(i,6),inmat(i,7),inmat(i,8),bf0TM,df0,dm0,width_inc);
                dat_temp(i1,1:3) = SolTM(1:3); % b_f, d_f and d_m
                
            end
            
        end
        
        if dat_temp(i1,1) == 0
            dat_temp(i1,4) = 0;
        end
        
        flag = 1;
        input_data = [Co(i1),inmat(i,1),inmat(i,2),inmat(i,3),inmat(i,4),inmat(i,5),inmat(i,6)*10^3*365*24*60*60,inmat(i,7)/60/60,inmat(i,8)];
        t = datetime;
        time = datevec(t);
        
        % fileID = fopen('C:\Users\fy23\Dropbox\Res_PC_2.txt','a');
        % fileID = fopen('/Users/Callisto/Dropbox/Res.txt','a');
        fprintf('%6d, %6d, %7.2f, %6.2f, %6.2f, %7.2f, %6.2f, %6.2f, %6.3f, %6.3f, %6d, %6d, %6d, %6.1f, %6.1f, %6d, %6d, %6d, %6d, %6d, %6d, %4d, %4.0f\n',[inmat(i,9),i1,dat_temp(i1,:),input_data,time]);
        % fcolse(fileID)
        
    end
    
        output(i,:) = reshape([dat_temp,input_data]',1,[]);
    
end

dat = reshape(output',15,[])';

% save SA_data_nau dat
figure
scatter(dat(:,9),dat(:,1))
hold on
scatter(dat(:,9),dat(:,4),'+')

end