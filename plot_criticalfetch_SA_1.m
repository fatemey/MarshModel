function plot_criticalfetch_SA_1
% plots parameters for sensitivity analysis purposes
%--------------------------------------------------------------------------------------------------
format compact
format longG
clear
clf

load bf_SA
y = dat;
load bfm_SA
y(y>10000&dat==10000) = 10000; 
y(y>5000&dat==5000) = 5000; 
y(y>1000&dat==1000) = 1000;
y = y/1000;

%------------------------------------------------
subplot(3,3,1)
hold on
load co_SA
x = dat*1000;
scatter(x,y,[],'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1)

% xlim([5 95])
% ylim([0 2.5])
xlabel('Ocean Concentration (mg/l)')
ylabel('Fetch Length (km)')
box on

%------------------------------------------------
subplot(3,3,2)
hold on
load cf_SA
x = dat*1000;
scatter(x,y,[],'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1)

% ylim([0 1.5])
% xlim([0 250])
xlabel('River Concentration (mg/l)')
ylabel('Fetch Length (km)')
box on

%------------------------------------------------
subplot(3,3,3)
hold on
load qf_SA
x = dat;
scatter(x,y,[],'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1)

% xlim([0 250])
% ylim([.3 .8])
xlabel('River Discharge (m^3/s)')
ylabel('Fetch Length (km)')
box on

%------------------------------------------------
subplot(3,3,4)
hold on
load le_SA
x = dat/1000;
scatter(x,y,[],'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1)

% xlim([1 5])
% ylim([.3 .8])
xlabel('Estuary Length (km)')
ylabel('Fetch Length (km)')
box on

%------------------------------------------------
subplot(3,3,5)
hold on
load bfm_SA
x = dat/1000;
scatter(x,y,[],'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1)

% xlim([1 5])
% ylim([.3 .8])
xlabel('Basin Width (km)')
ylabel('Fetch Length (km)')
box on

%------------------------------------------------
subplot(3,3,6)
hold on
load R_SA
x = dat;
scatter(x,y,[],'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1)

% xlim([1 20])
% ylim([0 .8])
xlabel('Rate of Sea Level Rise (mm/yr)')
ylabel('Fetch Length (km)')
box on

%------------------------------------------------
subplot(3,3,7)
hold on
load a_SA
x = dat;
scatter(x,y,[],'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1)

% xlim([1 10])
% ylim([0 5])
xlabel('Tidal Range (m)')
ylabel('Fetch Length (km)')
box on

%------------------------------------------------
subplot(3,3,8)
hold on
load vw_SA
x = dat;
scatter(x,y,[],'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1)

% xlim([0 15])
% ylim([0 5])
xlabel('Wind Velocity (m/s)')
ylabel('Fetch Length (km)')
box on

%------------------------------------------------
set(findobj('type','axes'),'fontsize',10)
h_fig=gcf;
set(h_fig,'PaperOrientation','portrait')
set(h_fig,'PaperPosition', [0 0 10 9]) % [... ... max_width=7.5 max_height=9]
tit='SA_1';
print(tit,'-dtiff','-r400')
movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\25 - Sensitivity Results')
close all

end