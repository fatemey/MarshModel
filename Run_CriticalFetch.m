function Run_CriticalFetch
% Runs BoxModel_SS_2eq.
%
% Last Update: 12/28/2017
%--------------------------------------------------------------------------------------------------
format compact
format longG

width_inc = 10; % tidal flat width increment for time marching approach

Co = ([5,10,15]) *10^-3;
Cf = ([10,30,50]) *10^-3;
Qf = [10,50,90];
LE = ([1,5,9]) * 1000;
bfm = ([5]) *1000;
a = (1.4);
R = ([2:4:10]) *10^-3/365/24/60/60;
T = 12 *60*60;
vw = [6];

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

output = zeros(m,15*k1);


%%parpool(2)

parfor i = 1 : m
    
    flag = 0;
    dat_temp = zeros(k1,15);
    
    for i1 = 1 : k1
        
        SolSS = CriticalFetchSS(Co(i1),input(i,1),input(i,2),input(i,3),input(i,4),input(i,5),input(i,6),input(i,7),input(i,8),100); % Sol = [x,y,fval_x,fval_y]
        err_SS = max(abs(SolSS(4:6)));
        
        if flag == 1
            
            bf0 = max(dat_temp(i1-1,1)-width_inc,1);
            SolTM = CriticalFetchTM(Co(i1),input(i,1),input(i,2),input(i,3),input(i,4),input(i,5),input(i,6),input(i,7),input(i,8),bf0,width_inc);
            dat_temp(i1,1:3) =SolTM(1:3); % b_f, d_f and d_m
            if err_SS < 10^-3 % if SS solution is relible
                dat_temp(i1,4:6) = SolSS(1:3);
            else
                dat_temp(i1,4:6) = [-1,-1,-1]; % showing that SS was not relible
            end
            
        else
            
            if err_SS < 10^-3 % if SS solution is relible
                if  SolSS(3)>a(i6)/2 % if marsh is drowned
                    dat_temp(i1,1) = 0; % b_f
                    dat_temp(i1,2:3) = SolSS(2:3);% d_f and d_m
                elseif SolSS(1)>=bfm(i5) % if b_f hits the upper boundary
                    dat_temp(i1,1) = bfm(i5); %b_f
                    dat_temp(i1,2:3) = SolSS(2:3);% d_f and d_m
                else
                    bf0 = max(SolSS(1)-5*width_inc,1);
                    SolTM = CriticalFetchTM(Co(i1),input(i,1),input(i,2),input(i,3),input(i,4),input(i,5),input(i,6),input(i,7),input(i,8),bf0,width_inc);
                    dat_temp(i1,1:3) =SolTM(1:3); % b_f, d_f and d_m
                end
                dat_temp(i1,4:6) = SolSS(1:3);
            else % if SS solution is not relible
                bf0 = 1;
                SolTM = CriticalFetchTM(Co(i1),input(i,1),input(i,2),input(i,3),input(i,4),input(i,5),input(i,6),input(i,7),input(i,8),bf0,width_inc);
                dat_temp(i1,1:3) =SolTM(1:3); % b_f, d_f and d_m
                dat_temp(i1,4:6) = [-1,-1,-1]; % showing that SS was not relible
            end
            
        end
        
        dat_temp(i1,7:15) = [Co(i1),input(i,1),input(i,2),input(i,3),input(i,4),input(i,5),input(i,6)*10^3*365*24*60*60,input(i,7)/60/60,input(i,8)];
        
        t = datetime;
        time = datevec(t);
        
        % fileID = fopen('C:\Users\fy23\Dropbox\Res_PC.txt','a');
        fileID = fopen('/Users/Callisto/Dropbox/Res.txt','a');
        % fileID = 1;
        % fprintf(fileID,'%6s %6s %7s %6s %6s %7s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s \n','i','i1','b_f','d_f','d_m','bfs','dfs','dms','Co','Cf','Qf','LE','bfm','a','R','T','VW','time');
        fprintf(fileID,'%6d %6d %7.2f %6.2f %6.2f %7.2f %6.2f %6.2f %6.3f %6.3f %6d %6d %6d %6.1f %6.1f %6d %6d %6d%2d%2d%2d%2d%2.0f\n',[i,i1,dat_temp(i1,:),time]);
%         fcolse(fileID)

        flag = 1;
        
    end
    
    output(i,:) = reshape(dat_temp',1,[]);
    
end

dat = reshape(output,15,[])';
save SA_data dat

end