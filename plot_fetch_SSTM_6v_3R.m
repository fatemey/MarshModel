function plot_fetch_SSTM_6v_3R
% plots fetch vs other variables. dta from ss and tm approaches
%--------------------------------------------------------------------------------------------------
format compact
format longG
clear
clf

%------------------------------------------------
subplot(2,3,1)
hold on
load co_data; a1 = dat(:,1)*1000; b1 = dat(:,2)/1000;
load co_data_2Rs; n = size(dat,1);
a2 = dat(1:n/2,4)*1000; b2 = dat(1:n/2,1)/1000;
a3 = dat(n/2+1:end,4)*1000; b3 = dat(n/2+1:end,1)/1000;

plot(a1,b1,a2,b2,a3,b3,'linewidth',3)

xlim([5 95])
ylim([0 2.5])
xlabel('Ocean Concentration (mg/l)')
ylabel('Critical Fetch Length (km)')
box on

%------------------------------------------------
subplot(2,3,2)
hold on
load cf_data; a1 = dat(:,1)*1000; b1 = dat(:,2)/1000;
load cf_data_2Rs; n = size(dat,1);
a2 = dat(1:n/2,5)*1000; b2 = dat(1:n/2,1)/1000;
a3 = dat(n/2+1:end,5)*1000; b3 = dat(n/2+1:end,1)/1000;
% a2 = dat(1:n/2,5)*1000; b2 = linspace(min(dat(1:n/2,1)),max(dat(1:n/2,1)),n/2)/1000;
% a3 = dat(n/2+1:end,5)*1000; b3 = linspace(min(dat(n/2+1:end,1)),max(dat(n/2+1:end,1)),n/2)/1000;

plot(a1,b1,a2,b2,a3,b3,'linewidth',3)

xlim([0 250])
ylim([.5 1])
xlabel('River Concentration (mg/l)')
ylabel('Critical Fetch Length (km)')
box on

%------------------------------------------------
subplot(2,3,3)
hold on
load qf_data; a1 = dat(:,1); b1 = dat(:,2)/1000;
load qf_data_2Rs; n = size(dat,1);
a2 = dat(1:n/2,6); b2 = dat(1:n/2,1)/1000;
a3 = dat(n/2+1:end,6); b3 = dat(n/2+1:end,1)/1000;

plot(a1,b1,a2,b2,a3,b3,'linewidth',3)

xlim([0 250])
ylim([.5 1])
xlabel('River Discharge (m^3/s)')
ylabel('Critical Fetch Length (km)')
box on
% NumTicks = 5;
% L = get(gca,'XLim');
% set(gca,'XTick',linspace(L(1),L(2),NumTicks))

%------------------------------------------------
subplot(2,3,4)
hold on
load le_data; a1 = dat(:,1)/1000; b1 = dat(:,2)/1000;
load le_data_2Rs; n = size(dat,1);
a2 = dat(1:n/2,7)/1000; b2 = dat(1:n/2,1)/1000;
a3 = dat(n/2+1:end,7)/1000; b3 = dat(n/2+1:end,1)/1000;

plot(a1,b1,a2,b2,a3,b3,'linewidth',3)

xlim([1 5])
ylim([.4 .8])
xlabel('Estuary Length (km)')
ylabel('Critical Fetch Length (km)')
box on

%------------------------------------------------
subplot(2,3,5)
hold on
load bfm_data; a1 = dat(:,1)/1000; b1 = dat(:,2)/1000;
load bfm_data_2Rs; n = size(dat,1);
a2 = dat(1:n/2,8)/1000; b2 = dat(1:n/2,1)/1000;
a3 = dat(n/2+1:end,8)/1000; b3 = dat(n/2+1:end,1)/1000;

plot(a1,b1,a2,b2,a3,b3,'linewidth',3)

xlim([1 5])
ylim([.4 .8])
xlabel('Basin Width (km)')
ylabel('Critical Fetch Length (km)')
box on

%------------------------------------------------
subplot(2,3,6)
hold on
% load vw_data; a1 = dat(:,1); b1 = dat(:,2)/1000;
load vw_data_2Rs_2222; a1 = dat(:,12); b1 = dat(:,1)/1000;
load vw_data_2Rs; n = size(dat,1);
a2 = dat(1:n/2,12); b2 = dat(1:n/2,1)/1000;
a3 = dat(n/2+1:end,12)+1;a3(1)=0; b3 = dat(n/2+1:end,1)/1000;

plot(a1,b1,a2,b2,a3,b3,'linewidth',3)

xlim([0 15])
ylim([0 5])
xlabel('Wind Velocity (m/s)')
ylabel('Critical Fetch Length (km)')
box on

%------------------------------------------------
set(findobj('type','axes'),'fontsize',10)
h_fig=gcf;
set(h_fig,'PaperOrientation','portrait')
set(h_fig,'PaperPosition', [0 0 10 6]) % [... ... max_width=7.5 max_height=9]
tit='fetchvs6v_SSTM_3R_1111';
print(tit,'-dtiff','-r400')
% movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\24 - Final Results')
movefile([tit,'.tif'],'/Users/Callisto/Dropbox/Matlab')
close all

end