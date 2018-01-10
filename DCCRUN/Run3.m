function Run3
% Runs CriticalFetchTM and CriticalFetchSS for a method using both steady
% state and time marching approaches.
%
% Last Update: 1/9/2018
%--------------------------------------------------------------------------------------------------
format compact
format longG

width_inc = 10; % tidal flat width increment for time marching approach
bf0SS = 100;
[input,Co,k1,m] = InputMaker3;
% output = zeros(m,15*k1);

parpool()

fprintf('%6s, %6s, %7s, %6s, %6s, %7s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %4s, %4s\n','i','iCo','bf','df','dm','bfss','dfss','dmss','Co','Cf','Qf','LE','bfm','a','R','T','VW','Yr','MM','DD','Hr','Min','Sec');

parfor i = 1 : m
    
    flag = 0;
    dat_temp = zeros(k1,15);
    
    for i1 = 1 : k1
        
        SolSS = CriticalFetchSS(Co(i1),input(i,1),input(i,2),input(i,3),input(i,4),input(i,5),input(i,6),input(i,7),input(i,8),bf0SS); % Sol = [x,y,fval_x,fval_y]
        err_SS = max(abs(SolSS(4:6)));
        
        if err_SS < 10^-3 % if SS solution is relible
            
            dat_temp(i1,4:6) = SolSS(1:3);
            
            if SolSS(1) >= input(i,4) % if b_f hits the upper boundary
                dat_temp(i1,1) = input(i,4); %b_f
                dat_temp(i1,2:3) = SolSS(2:3);% d_f and d_m
                dat_temp(i1,4) = input(i,4); % b_f
                
            else
                if flag == 0 % if this is the first run for Co loop
                    bf0TM = max(SolSS(1)-5*width_inc,1);
                    SolTM = CriticalFetchTM(Co(i1),input(i,1),input(i,2),input(i,3),input(i,4),input(i,5),input(i,6),input(i,7),input(i,8),bf0TM,width_inc);
                    dat_temp(i1,1:3) = SolTM(1:3); % b_f, d_f and d_m
                    
                else % if this is not the first run for Co loop
                    if input(i,4)-dat_temp(i1-1,1) <= width_inc
                        dat_temp(i1,1:6) = dat_temp(i1-1,1:6);
                    else
                        bf0TM = max(dat_temp(i1-1,1)-1,1);
                        SolTM = CriticalFetchTM(Co(i1),input(i,1),input(i,2),input(i,3),input(i,4),input(i,5),input(i,6),input(i,7),input(i,8),bf0TM,width_inc);
                        dat_temp(i1,1:3) =SolTM(1:3); % b_f, d_f and d_m
                    end
                    
                end
                
            end
            
        else % if SS solution is not relible
            
            dat_temp(i1,4:6) = [-1,-1,-1]; % showing that SS was not relible
            
            if flag == 0 % if this is the first run for Co loop
                bf0TM = 1;
                SolTM = CriticalFetchTM(Co(i1),input(i,1),input(i,2),input(i,3),input(i,4),input(i,5),input(i,6),input(i,7),input(i,8),bf0TM,width_inc);
                dat_temp(i1,1:3) =SolTM(1:3); % b_f, d_f and d_m
                
            else % if this is not the first run for Co loop
                bf0TM = max(dat_temp(i1-1,1)-1,1);
                SolTM = CriticalFetchTM(Co(i1),input(i,1),input(i,2),input(i,3),input(i,4),input(i,5),input(i,6),input(i,7),input(i,8),bf0TM,width_inc);
                dat_temp(i1,1:3) =SolTM(1:3); % b_f, d_f and d_m
                
            end
            
        end
        
        dat_temp(i1,7:15) = [Co(i1),input(i,1),input(i,2),input(i,3),input(i,4),input(i,5),input(i,6)*10^3*365*24*60*60,input(i,7)/60/60,input(i,8)];
        
        flag = 1;
        t = datetime;
        time = datevec(t);
        
        % fileID = fopen('C:\Users\fy23\Dropbox\Res_PC_2.txt','a');
        % fileID = fopen('/Users/Callisto/Dropbox/Res.txt','a');
        fprintf('%6d, %6d, %7.2f, %6.2f, %6.2f, %7.2f, %6.2f, %6.2f, %6.3f, %6.3f, %6d, %6d, %6d, %6.1f, %6.1f, %6d, %6d, %6d, %6d, %6d, %6d, %4d, %4.0f\n',[i,i1,dat_temp(i1,:),time]);
        % fcolse(fileID)
        
    end
    
    %     output(i,:) = reshape(dat_temp',1,[]);
    
end

% dat = reshape(output',15,[])';

% save SA_data_nau dat

end