function Run_CriticalFetch_SSTM
% Runs BoxModel_SS_2eq.
%
% Last Update: 12/5/2017
%--------------------------------------------------------------------------------------------------
format compact
format longG

fileID = fopen('C:\Users\fy23\Dropbox\Res_PC.txt','a');
fprintf(fileID,'%5s %7s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s %6s \n','n/162','b_f','d_f','d_m','Co','Cf','Qf','LE','bfm','a','R','T','VW','time');

width_inc = 10; % tidal flat width increment for time marching approach

Co = (20) *10^-3;
Cf = (15) *10^-3;
Qf = [20];
LE = ([5]) * 1000;
bfm = (5) *1000;
a = (1.4);
R = ([2]) *10^-3/365/24/60/60;
T = 12 *60*60;
vw = 0:.5:15;

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
for i7 =1 : k7
    flag = 0;
    for i2 = 1 : k2
        for i3 = 1 : k3
            for i1 = 1 : k1
                for i5 = 1 : k5
                    for i6 = 1 : k6
                        for i4 = 1 : k4
                            for i8 = 1 : k8
                                for i9 = 1 : k9
                                    
                                    t1=clock;
                                    
                                    SolSS = CriticaIFetchSS(Co(i1),Cf(i2),Qf(i3),LE(i4),bfm(i5),a(i6),R(i7),T(i8),vw(i9),100); % Sol = [x,y,fval_x,fval_y]
                                    if ( SolSS(1)>=bfm(i5) || SolSS(3)>a(i6)/2 ) && max(abs(SolSS(4:6)))<10^-3 % condition for when b_f hits the upper boundary / marsh drowns and SS solution is relible
                                        dat(n,1) = min(SolSS(1),bfm(i5)); %b_f
                                        dat(n,2:3) = SolSS(2:3);% d_f and d_m
                                    else
                                        if width_inc>0
                                            if flag == 1
                                                bf0 =  max(dat(n-1,1)-2*width_inc,1);
                                            elseif max(abs(SolSS(4:6)))>10^-3 % condition for when SS solution is not relible
                                                bf0 = 1;
                                            else
                                                bf0 = min(max(SolSS(1)-20*width_inc,1),bfm(i5)-2*width_inc); % b_f ; condition for when b_f hits the upper boundary
                                            end
                                        else
                                            if flag == 1
                                                bf0 =  max(dat(n-1,1)-2*width_inc,1);
                                            elseif max(abs(SolSS(4:6)))>10^-3 % condition for when SS solution is not relible
                                                bf0 = bfm(i5)-1;
                                            else
                                                bf0 = max(SolSS(1),bfm(i5)-1); % b_f ; condition for when b_f hits the upper boundary
                                            end
                                        end
                                        SolTM = CriticaIFetchTM(Co(i1),Cf(i2),Qf(i3),LE(i4),bfm(i5),a(i6),R(i7),T(i8),vw(i9),bf0,width_inc);
                                        dat(n,1:3) =SolTM(1:3); % b_f, d_f and d_m
                                    end
                                    
                                    dat(n,4:12) = [Co(i1),Cf(i2),Qf(i3),LE(i4),bfm(i5),a(i6),R(i7)*10^3*365*24*60*60,T(i8)/60/60,vw(i9)];
                                    
                                    t =etime(clock,t1)/60;
                                    fprintf(fileID,'%5d %7.2f %6.2f %6.2f %6.2f %6.2f %6d %6d %6d %6.1f %6.1f %6d %6d %6.1f\n',[n,dat(n,:),t]);
                                    n = n + 1;
                                    flag = 0;
                                    
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

save vw_data_2Rs_2222 dat

end