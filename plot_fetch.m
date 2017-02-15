function plot_fetch
% Function plot_fetch plots fetch threshold values (computed from functin fetch_threshold) vs diffrent
% parameters. 
%--------------------------------------------------------------------------------------------------

%------------------------------------------------
subplot(3,3,1)
load x_co; a = dat(:,1); b = dat(:,2);
scatter(a,b,'k','.')
xlabel('Ocean Concentration (mg/l)')
ylabel('Threshold Width (m)')
box on

%------------------------------------------------
subplot(3,3,2)
load x_cf; a = dat(:,1); b = dat(:,2);
scatter(a,b,'k','.')
xlabel('River Concentration (mg/l)')
ylabel('Threshold Width (m)')
box on

%------------------------------------------------
subplot(3,3,3)
load x_qf; a = dat(:,1); b = dat(:,2);
scatter(a,b,'k','.')
xlabel('River Discharge (m^3/s)')
ylabel('Threshold Width (m)')
box on

%------------------------------------------------
subplot(3,3,4)
load x_le; a = dat(:,1); b = dat(:,2);
scatter(a,b,'k','.')
xlabel('Estuary Length (km)')
ylabel('Threshold Width (m)')
box on

%------------------------------------------------
subplot(3,3,5)
load x_bfm; a = dat(:,1); b = dat(:,2);
scatter(a,b,'k','.')
xlabel('Basin Width (km)')
ylabel('Threshold Width (m)')
box on

%------------------------------------------------
subplot(3,3,6)
load x_vw; a = dat(:,1); b = dat(:,2);
scatter(a,b,'k','.')
xlabel('Wind Speed (m/s)')
ylabel('Threshold Width (m)')
box on

%------------------------------------------------
subplot(3,3,7)
load x_r; a = dat(:,1); b = dat(:,2);
scatter(a,b,'k','.')
xlabel('Rate of Sea Level Rise (mm/yr)')
ylabel('Threshold Width (m)')
box on

%------------------------------------------------
subplot(3,3,8)
load x_2h; a = dat(:,1); b = dat(:,2);
scatter(a,b,'k','.')
xlabel('Tidal Range (m)')
ylabel('Threshold Width (m)')
box on

%------------------------------------------------
set(findobj('type','axes'),'fontsize',10)
h_fig=gcf;
set(h_fig,'PaperOrientation','portrait')
set(h_fig,'PaperPosition', [0 0 7.5 6]) % [... ... max_width=7.5 max_height=9]
tit='threshold';
print(tit,'-dtiff','-r400')
movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\14 - Fetch threshold')
close all