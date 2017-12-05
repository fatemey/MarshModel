function plot_criticalfetch_SA_2
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

ran = rand(size(y))*10;

c = y;
[c_unique,I,J] = unique(c);
col_map = jet(length(c_unique));
col = col_map(J,:);

%------------------------------------------------
subplot(2,2,1)
hold on
load cf_SA
x1 = dat*1000+rand(size(y))*10;
load qf_SA
x2 = dat+rand(size(y))*10;
scatter(x1,x2,1,col,'filled')

xlabel('River Concentration (mg/l)')
ylabel('River Discharge (m^3/s)')
box on

%------------------------------------------------
subplot(2,2,2)
hold on
load le_SA
x1 = dat/1000+rand(size(y))*1;
load qf_SA
x2 = dat+rand(size(y))*10;
scatter(x1,x2,1,col,'filled')

xlabel('Estuary Length (km)')
ylabel('River Discharge (m^3/s)')
box on

%------------------------------------------------
subplot(2,2,3)
hold on
load le_SA
x1 = dat/1000+rand(size(y))*1;
load bfm_SA
x2 = dat/1000+rand(size(y))*1;
scatter(x1,x2,1,col,'filled')

xlabel('Estuary Length (km)')
ylabel('Basin Width (km)')
box on

%------------------------------------------------
subplot(2,2,4)
hold on
load vw_SA
x1 = dat+rand(size(y))*.5;
load bfm_SA
x2 = dat/1000+rand(size(y))*1;
scatter(x1,x2,1,col,'filled')

xlabel('Wind Velocity (m/s)')
ylabel('Basin Width (km)')
box on


%------------------------------------------------
set(findobj('type','axes'),'fontsize',10)
h_fig=gcf;
set(h_fig,'PaperOrientation','portrait')
set(h_fig,'PaperPosition', [0 0 6 6]) % [... ... max_width=7.5 max_height=9]
tit='SA_2';
print(tit,'-dtiff','-r400')
movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\25 - Sensitivity Results')
close all

end


%%
% make colorbar
% load bf_SA
% y = dat;
% load bfm_SA
% y(y>10000&dat==10000) = 10000; 
% y(y>5000&dat==5000) = 5000; 
% y(y>1000&dat==1000) = 1000;
% y = y/1000;
% v =linspace(0,max(y),5);
% subplot(2,2,3)
% colormap('jet')
% colorbar;
% c1=colorbar('Ticks',[0,1/4,1/2,3/4,1],'TickLabels',v)
% c1.Label.String = 'Critical Fetch (km)';
% set(findobj('type','axes'),'fontsize',10)
% h_fig=gcf;
% set(h_fig,'PaperOrientation','portrait')
% set(h_fig,'PaperPosition', [0 0 6 6]) % [... ... max_width=7.5 max_height=9]
% tit='SA_2_cbar';
% print(tit,'-dtiff','-r400')
% movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\25 - Sensitivity Results')
% close all