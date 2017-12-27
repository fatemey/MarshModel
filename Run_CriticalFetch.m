function Run_CriticalFetch
% Runs BoxModel_SS_2eq.
%
% Last Update: 12/22/2017
%--------------------------------------------------------------------------------------------------
format compact
format longG

% fileID = fopen('C:\Users\fy23\Dropbox\Res_PC.txt','a');
% fileID = fopen('/Users/Callisto/Dropbox/Res.txt','a');
fileID = 1;
fprintf(fileID,'%6s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s \n','n','b_f','d_f','d_m','Co','Cf','Qf','LE','bfm','a','R','T','VW','time');

width_inc = 10; % tidal flat width increment for time marching approach

Co = ([5,10,15,20]) *10^-3;
Cf = ([5:15:80]) *10^-3;
Qf = [5:15:80];
LE = ([1:2:9]) * 1000;
bfm = ([1:2:9]) *1000;
a = (1:1:3);
R = ([2:4:14]) *10^-3/365/24/60/60;
T = 12 *60*60;
vw = [2:4:10];

k1 = length(Co);
k2 = length(Cf);
k3 = length(Qf);
k4 = length(LE);
k5 = length(bfm);
k6 = length(a);
k7 = length(R);
k8 = length(T);
k9 = length(vw);

m = k1*k2*k3*k4*k5*k6*k7*k8*k9;
dat = zeros(m,12);
n = 1;
for i9 =1 : k9
    for i2 = 1 : k2
        for i3 = 1 : k3
            for i4 = 1 : k4
                for i5 = 1 : k5
                    for i6 = 1 : k6
                        for i7 = 1 : k7
                            for i8 = 1 : k8
                                
                                flag = 0;
                                
                                for i1 = 1 : k1
                                    
                                    if flag == 1
                                        bf0 = max(dat(n-1,1)-width_inc,1);
                                        SolTM = CriticalFetchTM(Co(i1),Cf(i2),Qf(i3),LE(i4),bfm(i5),a(i6),R(i7),T(i8),vw(i9),bf0,width_inc);
                                        dat(n,1:3) =SolTM(1:3); % b_f, d_f and d_m
                                    else
                                        SolSS = CriticalFetchSS(Co(i1),Cf(i2),Qf(i3),LE(i4),bfm(i5),a(i6),R(i7),T(i8),vw(i9),100); % Sol = [x,y,fval_x,fval_y]
                                        if max(abs(SolSS(4:6)))<10^-3 % if SS solution is relible
                                            if  SolSS(3)>a(i6)/2 % if marsh is drowned
                                                dat(n,1) = 0; % b_f
                                                dat(n,2:3) = SolSS(2:3);% d_f and d_m
                                            elseif SolSS(1)>=bfm(i5) % if b_f hits the upper boundary
                                                dat(n,1) = bfm(i5); %b_f
                                                dat(n,2:3) = SolSS(2:3);% d_f and d_m
                                            else
                                                bf0 = max(SolSS(1)-5*width_inc,1);
                                                SolTM = CriticalFetchTM(Co(i1),Cf(i2),Qf(i3),LE(i4),bfm(i5),a(i6),R(i7),T(i8),vw(i9),bf0,width_inc);
                                                dat(n,1:3) =SolTM(1:3); % b_f, d_f and d_m
                                            end
                                        else % if SS solution is not relible
                                            bf0 = 1;
                                            SolTM = CriticalFetchTM(Co(i1),Cf(i2),Qf(i3),LE(i4),bfm(i5),a(i6),R(i7),T(i8),vw(i9),bf0,width_inc);
                                            dat(n,1:3) =SolTM(1:3); % b_f, d_f and d_m
                                        end
                                    end
                                    
                                    dat(n,4:12) = [Co(i1),Cf(i2),Qf(i3),LE(i4),bfm(i5),a(i6),R(i7)*10^3*365*24*60*60,T(i8)/60/60,vw(i9)];
                                    
                                    t = datetime;
                                    time = datevec(t);
                                    fprintf(fileID,'%6d %6.2f %6.2f %6.2f %6.3f %6.3f %6d %6d %6d %6.1f %6.1f %6d %6d %6d%2d%2d%2d%2d%2.1f\n',[n,dat(n,:),time]);
                                    n = n + 1;
                                    flag = 1;
                                    
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

save SA_data dat

end