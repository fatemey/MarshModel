function plot_fetch
% Function plot_fetch plots fetch threshold values (computed from functin fetch_threshold) vs diffrent
% parameters. 
%--------------------------------------------------------------------------------------------------

%------------------------------------------------
subplot(3,3,1)
load co; a = dat(:,1)*1000; b = dat(:,2);
scatter(a,b,'k','.')
xlabel('Ocean Concentration (mg/l)')
ylabel('Critical Fetch (m)') 
box on

%------------------------------------------------
subplot(3,3,2)
load cf; a = dat(:,1)*1000; b = dat(:,2);
scatter(a,b,'k','.')
xlabel('River Concentration (mg/l)')
ylabel('Critical Fetch (m)') 
box on

%------------------------------------------------
subplot(3,3,3)
load qf; a = dat(:,1); b = dat(:,2);
scatter(a,b,'k','.')
xlabel('River Discharge (m^3/s)')
ylabel('Critical Fetch (m)') 
box on

%------------------------------------------------
subplot(3,3,4)
load le; a = dat(:,1)/1000; b = dat(:,2);
scatter(a,b,'k','.')
xlabel('Estuary Length (km)')
ylabel('Critical Fetch (m)') 
box on

%------------------------------------------------
subplot(3,3,5)
load bfm; a = dat(:,1)/1000; b = dat(:,2);
scatter(a,b,'k','.')
xlabel('Basin Width (km)')
ylabel('Critical Fetch (m)') 
box on

%------------------------------------------------
subplot(3,3,6)
load vw; a = dat(:,1); b = dat(:,2);
scatter(a,b,'k','.')
xlabel('Wind Speed (m/s)')
ylabel('Critical Fetch (m)') 
box on

%------------------------------------------------
subplot(3,3,7)
load R; a = dat(:,1)/10^-3*365*24*60*60; b = dat(:,2);
scatter(a,b,'k','.')
xlabel('Rate of Sea Level Rise (mm/yr)')
ylabel('Critical Fetch (m)') 
box on
xlim([0 50])

%------------------------------------------------
subplot(3,3,8)
load H; a = dat(:,1)*2; b = dat(:,2);
scatter(a,b,'k','.')
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