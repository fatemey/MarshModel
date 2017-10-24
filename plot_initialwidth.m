function plot_initialwidth
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
% load co_data_SS_opt; b2 = dat(:,2); dtf2 = dat(:,3); dm2 = dat(:,4);
load Sol_Co_SS; b2 = Sol(:,1); dtf2 = Sol(:,2); dm2 = Sol(:,3);
subplot(1,2,1) 
hold on
scatter(a,b2,80,[77, 148, 255]/255,'x','LineWidth',1.1)
gscatter(a,b,conv,[77, 148, 255]/255,[],[],'off')
xlabel('Ocean Concentration (mg/l)')
ylabel('Tidal Flat Width (m)') 
legend('Steady State (SS)','Time Marching (TM)','location','northwest')
box on

subplot(1,2,2)
hold on
% gscatter(a,dtf,eqdtf,'rk',[],[],'off')
% gscatter(a,dm,eqdm,'rk',[],[],'off')
p1=scatter(a(8:end),dtf(8:end),50,[0, 153, 0]/255,'o','filled');
scatter(a(1:7),dm(1:7),50,[0, 153, 0]/255,'o','filled')
p2=scatter(a,dtf2,90,[0, 153, 0]/255,'x','LineWidth',1.1);
p3=scatter(a,dm,20,[.8 .1 0],'o','filled');
p4=scatter(a,dm2,90,[.8 .1 0],'x','LineWidth',1.1);
plot([0,105],[.7,.7],'--k', 'LineWidth',.3)
text(90,.73,'MSL')
legend([p1,p2,p3,p4],'Tidal Flat - TM','Tidal Flat - SS','Marsh - TM','Marsh - SS','location','northwest')
xlim([0,105])
xlabel('Ocean Concentration (mg/l)')
ylabel('Depth (m)') 
box on

%------------------------------------------------
set(findobj('type','axes'),'fontsize',10)
h_fig=gcf;
set(h_fig,'PaperOrientation','portrait')
set(h_fig,'PaperPosition', [0 0 12 5]) % [... ... max_width=7.5 max_height=9]
tit='equilibrium solutions';
print(tit,'-dtiff','-r400')
movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\23 - Steady State Approach')
close all