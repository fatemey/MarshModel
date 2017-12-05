function plot_fetch_SSTM
% plots fetch vs other variables. dta from ss and tm approaches
%--------------------------------------------------------------------------------------------------
format compact
format longG
clear
clf

%------------------------------------------------
subplot(3,3,1)
hold on
load plotdataco; as = dat(:,4)*1000; bs = dat(:,1)/1000;
load co_data; at = dat(:,1)*1000; bt = dat(:,2)/1000;
hs=area(as,bs);
hs.FaceColor=[.8 .8 .8];
ht=patch('XData',[5;at;95],'YData',[0;bt;0]);
hatchfill(ht,'single', 75, 4);

xlim([5 95])
ylim([0 2.5])
xlabel('Ocean Concentration (mg/l)')
ylabel('Fetch Length (km)')
box on

%------------------------------------------------
subplot(3,3,2)
hold on
load cf_data; at = dat(:,1)*1000; bt = dat(:,2)/1000;c = bt<.766;
load plotdatacf; as = dat(:,5)*1000; bs = dat(:,1)/1000;
hs=area(as,bs);
hs.FaceColor=[.8 .8 .8];
ht=patch('XData',[0;at;250],'YData',[0;bt;0]);
hatchfill(ht,'single', 75, 4);

ylim([0 1.5])
xlim([0 250])
xlabel('River Concentration (mg/l)')
ylabel('Fetch Length (km)')
box on

%------------------------------------------------
subplot(3,3,3)
hold on
load qf_data; at = dat(:,1); bt = dat(:,2)/1000;
load plotdataqf; as = dat(:,7); bs = dat(:,1)/1000;
hs=area(as,bs);
hs.FaceColor=[.8 .8 .8];
ht=patch('XData',[0;at;250],'YData',[0;bt;0]);
hatchfill(ht,'single', 75, 4);

xlim([0 250])
ylim([.3 .8])
xlabel('River Discharge (m^3/s)')
ylabel('Fetch Length (km)')
box on
% NumTicks = 5;
% L = get(gca,'XLim');
% set(gca,'XTick',linspace(L(1),L(2),NumTicks))

%------------------------------------------------
subplot(3,3,4)
hold on
load le_data; at = dat(:,1)/1000; bt = dat(:,2)/1000;
load plotdatal; as = dat(:,6)/1000; bs = dat(:,1)/1000;
hs=area(as,bs);
hs.FaceColor=[.8 .8 .8];
ht=patch('XData',[1;at;5],'YData',[0;bt;0]);
hatchfill(ht,'single', 75, 4);

xlim([1 5])
ylim([.3 .8])
xlabel('Estuary Length (km)')
ylabel('Fetch Length (km)')
box on
% NumTicks = 5;
% L = get(gca,'XLim');
% set(gca,'XTick',linspace(L(1),L(2),NumTicks))

%------------------------------------------------
subplot(3,3,5)
hold on
load bfm_data; at = dat(:,1)/1000; bt = dat(:,2)/1000;
load plotdatab; as = dat(:,11)/1000; bs = dat(:,1)/1000;
hs=area(as,bs);
hs.FaceColor=[.8 .8 .8];
ht=patch('XData',[1;at;5],'YData',[0;bt;0]);
hatchfill(ht,'single', 75, 4);

xlim([1 5])
ylim([.3 .8])
xlabel('Basin Width (km)')
ylabel('Fetch Length (km)')
box on
% NumTicks = 5;
% L = get(gca,'XLim');
% set(gca,'XTick',linspace(L(1),L(2),NumTicks))

%------------------------------------------------
subplot(3,3,6)
hold on
load R_data; at = dat(:,1)/10^-3*365*24*60*60; bt = dat(:,2)/1000;
load plotdatar; as = dat(:,8)/10^-3*365*24*60*60; bs = dat(:,1)/1000;
hs=area(as,bs);
hs.FaceColor=[.8 .8 .8];
ht=patch('XData',[1;at;20],'YData',[0;bt;0]);
hatchfill(ht,'single', 75, 4);

xlim([1 20])
ylim([0 .8])
xlabel('Rate of Sea Level Rise (mm/yr)')
ylabel('Fetch Length (km)')
box on

%------------------------------------------------
subplot(3,3,7)
hold on
load H_data; at = dat(:,1)*2; bt = dat(:,2)/1000;
load plotdatah; as = dat(:,9)*2; bs = dat(:,1)/1000;
hs=area(as,bs);
hs.FaceColor=[.8 .8 .8];
ht=patch('XData',[1;at;10],'YData',[0;bt;0]);
hatchfill(ht,'single', 75, 4);

xlim([1 10])
ylim([0 5])
xlabel('Tidal Range (m)')
ylabel('Fetch Length (km)')
box on
% NumTicks = 5;
% L = get(gca,'XLim');
% set(gca,'XTick',linspace(L(1),L(2),NumTicks))

%------------------------------------------------
subplot(3,3,8)
hold on
load T_data; at = dat(:,1)/60/60; bt = dat(:,2)/1000;
load plotdatat; as = dat(:,10)/60/60; bs = dat(:,1)/1000;
hs=area(as,bs);
hs.FaceColor=[.8 .8 .8];
ht=patch('XData',[12;at;24],'YData',[0;bt;0]);
hatchfill(ht,'single', 75, 4);

xlim([12 24])
ylim([.3 .8])
xlabel('Tidal Period (hr)')
ylabel('Fetch Length (km)')
box on
NumTicks = 2;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))

%------------------------------------------------
subplot(3,3,9)
hold on
load vw_data; at = dat(:,1); bt = dat(:,2)/1000;
load plotdatav; as = dat(:,12); bs = dat(:,1)/1000;
hs=area(as,bs);
hs.FaceColor=[.8 .8 .8];
ht=patch('XData',[0;at;15],'YData',[0;bt;0]);
hatchfill(ht,'single', 75, 4);

xlim([0 15])
ylim([0 5])
xlabel('Wind Velocity (m/s)')
ylabel('Fetch Length (km)')
box on

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
tit='fetchvs9v_SSTM';
print(tit,'-dtiff','-r400')
movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\24 - Final Results')
close all

end