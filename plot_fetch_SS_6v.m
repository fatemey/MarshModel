function plot_fetch_SS_6v
% plots fetch vs other variables
%--------------------------------------------------------------------------------------------------
format compact
format longG
clear
clf

%------------------------------------------------
subplot(2,3,1)
hold on
load plotdataco; a = dat(:,4)*1000; b = dat(:,1)/1000;
h=area(a,b);
h.FaceColor=[.8 .8 .8];
% plot(a(c),b(c),'r')
text(23,2,'Tidal Flat Dominance','fontsize',8)
text(40,.5,'Marsh Dominance','fontsize',8)
ht=text(58,1.3,'Critical Fetch Length','fontsize',7,'FontAngle', 'italic');
set(ht,'Rotation',41);
% NumTicks = 6;
% L = get(gca,'XLim');
% set(gca,'XTick',linspace(L(1),L(2),NumTicks))

xlim([5 95])
xlabel('Ocean Concentration (mg/l)')
ylabel('Fetch Length (km)')
box on

%------------------------------------------------
subplot(2,3,2)
hold on
load plotdatacf; a = dat(:,5)*1000; b = dat(:,1)/1000;
h=area(a,b);
h.FaceColor=[.8 .8 .8];
% plot(a(c),b(c),'r')
text(70,2,'Tidal Flat Dominance','fontsize',8)
text(70,.35,'Marsh Dominance','fontsize',8)
ht=text(120,.78,'Critical Fetch Length','fontsize',7,'FontAngle', 'italic');
set(ht,'Rotation',10);

xlim([0 250])
ylim([0 2.5])
xlabel('River Concentration (mg/l)')
ylabel('Fetch Length (km)')
box on
% NumTicks = 5;
% L = get(gca,'XLim');
% set(gca,'XTick',linspace(L(1),L(2),NumTicks))

%------------------------------------------------
subplot(2,3,3)
hold on
load plotdataqf; a = dat(:,7); b = dat(:,1)/1000;
h=area(a,b);
h.FaceColor=[.8 .8 .8];
% plot(a(c),b(c),'r')
text(70,.52,'Tidal Flat Dominance','fontsize',8)
text(70,.26,'Marsh Dominance','fontsize',8)
ht=text(30,.44,'Critical Fetch Length','fontsize',7,'FontAngle', 'italic');
set(ht,'Rotation',-16.8);

xlim([0 250])
ylim([.2 .6])
xlabel('River Discharge (m^3/s)')
ylabel('Fetch Length (km)')
box on
% NumTicks = 5;
% L = get(gca,'XLim');
% set(gca,'XTick',linspace(L(1),L(2),NumTicks))

%------------------------------------------------
subplot(2,3,4)
hold on
load plotdatal; a = dat(:,6)/1000; b = dat(:,1)/1000;
h=area(a,b);
h.FaceColor=[.8 .8 .8];
% plot(a(c),b(c),'r')
text(2,.475,'Tidal Flat Dominance','fontsize',8)
text(2,.3,'Marsh Dominance','fontsize',8)
ht=text(3,.433,'Critical Fetch Length','fontsize',7,'FontAngle', 'italic');
set(ht,'Rotation',3);

% xlim([1 19])
ylim([.25 .5])
xlabel('Estuary Length (km)')
ylabel('Fetch Length (km)')
box on
% NumTicks = 5;
% L = get(gca,'XLim');
% set(gca,'XTick',linspace(L(1),L(2),NumTicks))

%------------------------------------------------
subplot(2,3,5)
hold on
load plotdatab; a = dat(:,11)/1000; b = dat(:,1)/1000;
h=area(a,b);
h.FaceColor=[.8 .8 .8];
% plot(a(c),b(c),'r')
text(2,.475,'Tidal Flat Dominance','fontsize',8)
text(2,.3,'Marsh Dominance','fontsize',8)
ht=text(3,.421,'Critical Fetch Length','fontsize',7,'FontAngle', 'italic');
set(ht,'Rotation',9);

% xlim([1 19])
ylim([.25 .5])
xlabel('Basin Width (km)')
ylabel('Fetch Length (km)')
box on
% NumTicks = 5;
% L = get(gca,'XLim');
% set(gca,'XTick',linspace(L(1),L(2),NumTicks))

%------------------------------------------------
subplot(2,3,6)
hold on
load plotdatav; a = dat(:,12); b = dat(:,1)/1000;
h=area(a,b);
h.FaceColor=[.8 .8 .8];
text(6,4,'Tidal Flat Dominance','fontsize',8)
text(.9,.9,{'Marsh','Dominance'},'fontsize',8)
ht=text(8,.35,'Critical Fetch Length','fontsize',7,'FontAngle', 'italic');
set(ht,'Rotation',-1.5);

% ylim([.25 .5])
xlim([0 15])
xlabel('Wind Velocity (m/s)')
ylabel('Fetch Length (km)')
box on

%------------------------------------------------
set(findobj('type','axes'),'fontsize',10)
h_fig=gcf;
set(h_fig,'PaperOrientation','portrait')
set(h_fig,'PaperPosition', [0 0 10 6]) % [... ... max_width=7.5 max_height=9]
tit='fetchvs6v_SS';
print(tit,'-dtiff','-r400')
% movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\24 - Final Results')
movefile([tit,'.tif'],'/Users/Callisto/Dropbox/Matlab')
close all

end