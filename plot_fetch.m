function plot_fetch
% Function plot_fetch plots fetch threshold values (computed from functin fetch_threshold) vs diffrent
% parameters. The size of the circles are the results from the function
% tidalflat_convertion.
%
%--------------------------------------------------------------------------------------------------
format compact
format longG
clear
% clf

%------------------------------------------------
subplot(3,3,1)
hold on
load co; a = dat(:,1)*1000; b = dat(:,2);
load co_conversion;
col = [f, f, f];
for i = 1 : length(a)
    scatter(a(i),b(i),10,'MarkerEdgeColor','k', 'MarkerFaceColor',col(i,:))
end
xlabel('Ocean Concentration (mg/l)')
ylabel('Critical Fetch (m)') 
box on

%------------------------------------------------
subplot(3,3,2)
hold on
load cf; a = dat(:,1)*1000; b = dat(:,2);
load cf_conversion
col = [f, f, f];
for i = 1 : length(a)
    scatter(a(i),b(i),10,'MarkerEdgeColor','k', 'MarkerFaceColor',col(i,:))
end
xlabel('River Concentration (mg/l)')
ylabel('Critical Fetch (m)') 
box on

%------------------------------------------------
subplot(3,3,3)
hold on
load qf; a = dat(:,1); b = dat(:,2);
load qf_conversion
col = [f, f, f];
for i = 1 : length(a)
    scatter(a(i),b(i),10,'MarkerEdgeColor','k', 'MarkerFaceColor',col(i,:))
end
xlabel('River Discharge (m^3/s)')
ylabel('Critical Fetch (m)') 
box on

%------------------------------------------------
subplot(3,3,4)
hold on
load le; a = dat(:,1)/1000; b = dat(:,2);
load le_conversion
col = [f, f, f];
for i = 1 : length(a)
    scatter(a(i),b(i),10,'MarkerEdgeColor','k', 'MarkerFaceColor',col(i,:))
end
xlabel('Estuary Length (km)')
ylabel('Critical Fetch (m)') 
box on

%------------------------------------------------
subplot(3,3,5)
hold on
load bfm; a = dat(:,1)/1000; b = dat(:,2);
load bfm_conversion
col = [f, f, f];
for i = 1 : length(a)
    scatter(a(i),b(i),10,'MarkerEdgeColor','k', 'MarkerFaceColor',col(i,:))
end
xlabel('Basin Width (km)')
ylabel('Critical Fetch (m)') 
box on

%------------------------------------------------
subplot(3,3,6)
hold on
load vw; a = dat(:,1); b = dat(:,2);
load vw_conversion
col = [f, f, f];
for i = 1 : length(a)
    scatter(a(i),b(i),10,'MarkerEdgeColor','k', 'MarkerFaceColor',col(i,:))
end
xlabel('Wind Speed (m/s)')
ylabel('Critical Fetch (m)') 
box on

%------------------------------------------------
subplot(3,3,7)
hold on
load R; a = dat(:,1)/10^-3*365*24*60*60; b = dat(:,2);
load R_conversion
col = [f, f, f];
for i = 1 : length(a)
    scatter(a(i),b(i),10,'MarkerEdgeColor','k', 'MarkerFaceColor',col(i,:))
end
xlabel('Rate of Sea Level Rise (mm/yr)')
ylabel('Critical Fetch (m)') 
box on
xlim([0 50])

%------------------------------------------------
subplot(3,3,8)
hold on
load H; a = dat(:,1)*2; b = dat(:,2);
load H_conversion
col = [f, f, f];
for i = 1 : length(a)
    scatter(a(i),b(i),10,'MarkerEdgeColor','k', 'MarkerFaceColor',col(i,:))
end
xlabel('Tidal Range (m)')
ylabel('Critical Fetch (m)') 
box on

%------------------------------------------------
set(findobj('type','axes'),'fontsize',10)
h_fig=gcf;
set(h_fig,'PaperOrientation','portrait')
set(h_fig,'PaperPosition', [0 0 7.5 6]) % [... ... max_width=7.5 max_height=9]
tit='threshold-1';
print(tit,'-dtiff','-r400')
movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\14 - Fetch threshold')
close all