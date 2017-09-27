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
clf

%------------------------------------------------
load co_data_1000yr; a = dat(:,1)*1000; b = dat(:,2); dtf = dat(:,3); dm = dat(:,4); eqdtf = dat(:,5); eqdm = dat(:,6); conv = dat(:,7);
load co_data_SS_opt; b2 = dat(:,2); dtf2 = dat(:,3); dm2 = dat(:,4);
subplot(1,2,1)
hold on
scatter(a,b2,80,[0, 153, 0]/255,'x','LineWidth',1.1)
gscatter(a,b,conv,'kk',[],[],'off')
xlabel('Ocean Concentration (mg/l)')
ylabel('Tidal Flat Width (m)') 
legend('Steady State (SS)','Time Marching (TM)','location','northwest')
box on

subplot(1,2,2)
hold on
% gscatter(a,dtf,eqdtf,'rk',[],[],'off')
% gscatter(a,dm,eqdm,'rk',[],[],'off')
scatter(a,dtf,20,[77, 148, 255]/255,'o','filled')
scatter(a,dtf2,90,[77, 148, 255]/255,'x','LineWidth',1.1)
scatter(a,dm,20,[.8 .1 0],'o','filled')
scatter(a,dm2,90,[.8 .1 0],'+','LineWidth',1.1)
plot([0,105],[.7,.7],'--k', 'LineWidth',.3)
text(90,.73,'MSL')
legend('Tidal Flat - SS','Tidal Flat - TM','Marsh - SS','Marsh - TM','location','best')
xlim([0,105])
xlabel('Ocean Concentration (mg/l)')
ylabel('Depth (m)') 
box on

%------------------------------------------------
set(findobj('type','axes'),'fontsize',10)
h_fig=gcf;
set(h_fig,'PaperOrientation','portrait')
set(h_fig,'PaperPosition', [0 0 12 5]) % [... ... max_width=7.5 max_height=9]
tit='equilibrium solutions-1kyr';
print(tit,'-dtiff','-r400')
movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\20 - Unstable Equlibrium Results through Optimization of Steady State')
close all