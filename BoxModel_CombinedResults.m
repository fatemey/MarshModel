function BoxModel_CombinedResults
% BoxModel_CombinedResults : This function splitts the time
% period and uses the function BoxModel_SplittedResults for
% each intervals and finally combines the results to overcome the problem
% of the code arose when the period was very long.
%
%--------------------------------------------------------------------------------------------------
format compact
format longG
clear
clf

tic

tyr = 1000;  % solve for time tyr (years)
t = 0:100:tyr;
b_fm = 5 *10^3; % total basin width (both sides of the channel) (m)

%-------------- Initial conditions, y0=[ b_f, d_f, d_m,u(=C_r*(b_f*d_f+b_m*d_m))]
y0(1) = b_fm/4;   % tidal flat width (m)
y0(2) = 1.0;         % tidal flat depth (m)
y0(3) = 0.4;         % marsh depth (m)
y0(4) = 0*10^-3*(y0(1)*(y0(2)+y0(3)));

ydata = [];
timedata = [];

for i = 1:length(t)-1
    
    i
    if rem(i,2) == 0
        v_w = 6;
    else
        v_w = 0;
    end
    
    [y, time] = MarshTidalFlat_BoxModel_SplittedResults (t(i), t(i+1), y0,v_w);
    y0 = y(end, :);
    y0(4) = y0(4)*(y0(1)*y0(2)+((b_fm/2-y0(1))*y0(3)));
    y(end,:) = [];
    time(end) = [];
    ydata = [ydata; y];
    timedata = [timedata; time];
    
end

timespent_min = toc/60

%-------------- Plot Results
plot_test(timedata,ydata)
tit = 'alt100yr';
print(tit,'-dtiff','-r400')
movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\11 - Different wind conditions')
close all

end