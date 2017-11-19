function plot_fetch_TM
% Function plot_fetch_TM plots critical initial width values (computed from functin CriticalInitialWidth) vs difzrent
% parameters. The size of the circles are the results from the function
% tidalflat_convertion.
%
% Last Update: 3/28/2017
%
%--------------------------------------------------------------------------------------------------
format compact
format longG
clear
% clf

%------------------------------------------------
subplot(3,3,1)
hold on
load co_data; a = dat(:,1)*1000; b = dat(:,2)/1000;

scatter(a,b,20,'k','o','filled');
% scatter(a(f==0),b(f==0),20,'k','o','filled')
% scatter(a(f==1),b(f==1),20,'k','>','filled') 
% scatter(a(f==2),b(f==2),20,'k','<','filled')
% scatter(a(f==3),b(f==3),20,'k','s','filled')
% scatter(a(f==4),b(f==4),20,'k','s');

% ylim([0 5000])
xlabel('Ocean Concentration (mg/l)')
ylabel('Fetch Length (km)') 
box on
NumTicks = 5;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))

%------------------------------------------------
subplot(3,3,2)
hold on
load cf_data; a = dat(:,1)*1000; b = dat(:,2)/1000;

scatter(a,b,20,'k','o','filled');
% scatter(a(f==0),b(f==0),20,'k','o','filled')
% scatter(a(f==1),b(f==1),20,'k','>','filled') 
% scatter(a(f==2),b(f==2),20,'k','<','filled')
% scatter(a(f==3),b(f==3),20,'k','s','filled')
% scatter(a(f==4),b(f==4),20,'k','s');

% ylim([0 5000])
xlabel('River Concentration (mg/l)')
ylabel('Fetch Length (km)') 
box on
NumTicks = 5;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))

%------------------------------------------------
subplot(3,3,3)
hold on
load qf_data; a = dat(:,1); b = dat(:,2)/1000;

scatter(a,b,20,'k','o','filled');
% scatter(a(f==0),b(f==0),20,'k','o','filled')
% scatter(a(f==1),b(f==1),20,'k','>','filled') 
% scatter(a(f==2),b(f==2),20,'k','<','filled')
% scatter(a(f==3),b(f==3),20,'k','s','filled')
% scatter(a(f==4),b(f==4),20,'k','s');

% ylim([0 5000])
xlabel('River Discharge (m^3/s)')
ylabel('Fetch Length (km)') 
box on
NumTicks = 5;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))

%------------------------------------------------
subplot(3,3,4)
hold on
load le_data; a = dat(:,1)/1000; b = dat(:,2)/1000;

scatter(a,b,20,'k','o','filled');
% scatter(a(f==0),b(f==0),20,'k','o','filled')
% scatter(a(f==1),b(f==1),20,'k','>','filled') 
% scatter(a(f==2),b(f==2),20,'k','<','filled')
% scatter(a(f==3),b(f==3),20,'k','s','filled')
% scatter(a(f==4),b(f==4),20,'k','s');

% ylim([0 5000])
xlabel('Estuary Length (km)')
ylabel('Fetch Length (km)') 
box on
NumTicks = 5;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))

%------------------------------------------------
subplot(3,3,5)
hold on
load bfm_data; a = dat(:,1)/1000; b = dat(:,2)/1000;

scatter(a,b,20,'k','o','filled');
% scatter(a(f==0),b(f==0),20,'k','o','filled')
% scatter(a(f==1),b(f==1),20,'k','>','filled') 
% scatter(a(f==2),b(f==2),20,'k','<','filled')
% scatter(a(f==3),b(f==3),20,'k','s','filled')
% scatter(a(f==4),b(f==4),20,'k','s');

% ylim([0 5000])
xlabel('Basin Width (km)')
ylabel('Fetch Length (km)') 
box on
NumTicks = 5;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))

%------------------------------------------------
subplot(3,3,9)
hold on
load vw_data; a = dat(:,1); b = dat(:,2)/1000;

scatter(a,b,20,'k','o','filled');
% scatter(a(f==0),b(f==0),20,'k','o','filled')
% scatter(a(f==1),b(f==1),20,'k','>','filled') 
% scatter(a(f==2),b(f==2),20,'k','<','filled')
% scatter(a(f==3),b(f==3),20,'k','s','filled')
% scatter(a(f==4),b(f==4),20,'k','s');

% ylim([0 5000])
xlabel('Wind Speed (m/s)')
ylabel('Fetch Length (km)') 
box on
NumTicks = 5;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))

%------------------------------------------------
subplot(3,3,6)
hold on
load R_data; a = dat(:,1)/10^-3*365*24*60*60; b = dat(:,2)/1000;

scatter(a,b,20,'k','o','filled');
% scatter(a(f==0),b(f==0),20,'k','o','filled')
% scatter(a(f==1),b(f==1),20,'k','>','filled') 
% scatter(a(f==2),b(f==2),20,'k','<','filled')
% scatter(a(f==3),b(f==3),20,'k','s','filled')
% scatter(a(f==4),b(f==4),20,'k','s');

% ylim([0 5000])
xlabel('Rate of Sea Level Rise (mm/yr)')
ylabel('Fetch Length (km)') 
box on
NumTicks = 5;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))
L = get(gca,'YLim');
set(gca,'YTick',linspace(L(1),L(2),NumTicks))
% xlim([0 30])

%------------------------------------------------
subplot(3,3,7)
% plot([ 0 10],[5000,5000],'r--')
hold on
load H_data; a = dat(:,1)*2; b = dat(:,2)/1000;

scatter(a,b,20,'k','o','filled');
% l1 = scatter(a(f==0),b(f==0),20,'k','o','filled');
% l2 = scatter(a(f==1),b(f==1),20,'k','>','filled') ;
% l3 = scatter(a(f==2),b(f==2),20,'k','<','filled');
% l4 = scatter(a(f==3),b(f==3),20,'k','s','filled');
% l5 = scatter(a(f==4),b(f==4),20,'k','s');

% ylim([0 5000])
xlabel('Tidal Range (m)')
ylabel('Fetch Length (km)') 
box on
NumTicks = 6;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))

%------------------------------------------------
subplot(3,3,8)
hold on
load T_data; a = dat(:,1)/60/60; b = dat(:,2)/1000;

scatter(a,b,20,'k','o','filled');
% l1 = scatter(a(f==0),b(f==0),20,'k','o','filled');
% l2 = scatter(a(f==1),b(f==1),20,'k','>','filled') ;
% l3 = scatter(a(f==2),b(f==2),20,'k','<','filled');
% l4 = scatter(a(f==3),b(f==3),20,'k','s','filled');
% l5 = scatter(a(f==4),b(f==4),20,'k','s');

% ylim([0 5000])
xlabel('Tidal Period (hr)')
ylabel('Fetch Length (km)') 
box on
NumTicks = 6;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))

%------------------------------------------------
% % hL = legend([l1,l2,l3,l4],{'Normal Expansion/Contraction','TF Expansion & Emergence','TF Contraction & Emergence','Marsh Drowning'});
% hL = legend([l1,l2,l3,l4,l5],{'Normal Expansion/Contraction','TF Expansion & Emergence','TF Contraction & Emergence','Marsh Expansion & Drowning','Marsh Contraction & Drowning'});
% newPosition = [0.7 0.15 0.2 0.2];
% newUnits = 'normalized';
% set(hL,'Position', newPosition,'Units', newUnits);
% legend boxoff 

%------------------------------------------------
set(findobj('type','axes'),'fontsize',10)
h_fig=gcf;
set(h_fig,'PaperOrientation','portrait')
set(h_fig,'PaperPosition', [0 0 10 9]) % [... ... max_width=7.5 max_height=9]
tit='fetchvs9v_TM';
print(tit,'-dtiff','-r400')
movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\24 - Final Results')
close all