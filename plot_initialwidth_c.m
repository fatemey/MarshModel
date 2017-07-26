function plot_initialwidth_c
% Function plot_initialwidth plots critical initial width values (computed from functin CriticalInitialWidth) vs different
% parameters.
%
% Last Update: 7/24/2017
%
%--------------------------------------------------------------------------------------------------
format compact
format longG
clear
% clf

%------------------------------------------------
col = linspecer(7);
load co_data; a = dat(:,1)*1000; b = dat(:,2); dtf = dat(:,3); dm = dat(:,4); eqdtf = dat(:,5); eqdm = dat(:,6); conv = dat(:,7);
load co_data_SS_opt; b2 = dat(:,2); dtf2 = dat(:,3); dm2 = dat(:,4);
subplot(1,2,1)
hold on
scatter(a,b2,[],col(3,:),'^','filled')
gscatter(a,b,conv,'rk',[],[],'off')
xlabel('Ocean Concentration (mg/l)')
ylabel('Tidal Flat Width (m)') 
box on

subplot(1,2,2)
hold on
gscatter(a,dtf,eqdtf,col,[],[],'off')
gscatter(a,dm,eqdm,col,[],[],'off')
scatter(a,dtf2,'kx')
scatter(a,dm2,'k+')
plot([0,105],[.7,.7],'k')
xlabel('Ocean Concentration (mg/l)')
ylabel('Depth (m)') 
box on

%------------------------------------------------
set(findobj('type','axes'),'fontsize',10)
h_fig=gcf;
set(h_fig,'PaperOrientation','portrait')
set(h_fig,'PaperPosition', [0 0 7 3.5]) % [... ... max_width=7.5 max_height=9]
tit='equilibrium solutions';
print(tit,'-dtiff','-r400')
movefile([tit,'.tif'],'/Users/Callisto/Files/Work/Marsh Model/Results/19 - Unstable Equlibrium Results through Numerical Approach')
close all