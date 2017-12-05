function Run_CriticalFetch_SA
% Runs CriticalFetch for sensitivity analysis
%
% Last Update: 11/25/2017
%--------------------------------------------------------------------------------------------------
format compact
format longG

fileID = fopen('C:\Users\fy23\Dropbox\Res_PC.txt','w');
fprintf(fileID,'%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s \n','b_f','d_f','d_m','Co','Cf','Qf','LE','bfm','a','R','VW');

Co = [5,50,100] *10^-3;
Cf =  [5,50,100] *10^-3;
Qf = [5,50,100];
LE = [1,5,10] *10^3;
bfm = [1,5,10] *10^3;
a = [1,2,3];
R = [1,7,15] /10^3/365/24/60/60;
vw = [4,6,8];

k1 = length(Co);
k2 = length(Cf);
k3 = length(Qf);
k4 = length(LE);
k5 = length(bfm);
k6 = length(a);
k7 = length(R);
k8 = length(vw);

dat = zeros(k1*k2*k3*k4*k5*k6*k7*k8,11);
n = 1;
for i1 =2 : k1
    for i2 = 1 : k2
        for i3 = 1 : k3
            for i4 = 1 : k4
                for i5 = 1 : k5
                    for i6 = 1 : k6
                        for i7 = 1 : k7
                            for i8 = 1 : k8
                                
                                SolSS = CriticaIFetchSS(Co(i1),Cf(i2),Qf(i3),LE(i4),bfm(i5),a(i6),R(i7),vw(i8),100); % Sol = [x,y,fval_x,fval_y]
                                bf0 = min(SolSS(1),bfm(i8)-1); % b_f ; condition for when b_f hits the upper boundary
                                if SolSS(3)>a(i6) || max(abs(SolSS(4:6)))>10^-3 % condition for when b_f hits the lower boundary OR SS solution is not relible
                                    bf0 = 1; % b_f ; condition for when b_f hits the lower boundary
                                end
                                
                                SolTM = CriticaIFetchTM(Co(i1),Cf(i2),Qf(i3),LE(i4),bfm(i5),a(i6),R(i7),vw(i8),bf0);
                                
                                dat(n,1) = max(SolSS(1),SolTM(1)); % b_f
                                dat(n,2:3) = SolTM(2:3); % d_f and d_m
                                dat(n,4:11) = [Co(i1),Cf(i2),Qf(i3),LE(i4),bfm(i5),a(i6),R(i7) *10^3*365*24*60*60,vw(i8)];
                                
                                fprintf(fileID,'%12f %12f %12f %12f %12f %12d %12d %12d %12d %12d %12d \n',dat(n,:));
                                n = n + 1
                                
                            end
                        end
                    end
                end
            end
        end
    end
end

save data_CriticalFetch_SA dat
fclose(fileID);

end


