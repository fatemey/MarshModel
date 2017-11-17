% plots fetch vs other variables
%--------------------------------------------------------------------------------------------------
format compact
format longG
clear
clf

%------------------------------------------------
subplot(3,3,1)
hold on
load co; a = dat(:,1)*1000; b = dat(:,2)/1000; c = a<=35;
% load plotdatac; a = dat(:,4)*1000; b = dat(:,1)/1000;
h=area(a,b);
h.FaceColor=[.8 .8 .8];
plot(a(c),b(c),'r')
text(23,2,'Tidal Flat Dominance','fontsize',8)
text(40,.5,'Marsh Dominance','fontsize',8)
ht=text(38,.91,'Critical Fetch Length','fontsize',7,'FontAngle', 'italic');
set(ht,'Rotation',41);
% NumTicks = 6;
% L = get(gca,'XLim');
% set(gca,'XTick',linspace(L(1),L(2),NumTicks))

xlim([5 100])
xlabel('Ocean Concentration (mg/l)')
ylabel('Fetch Length (km)')
box on

%------------------------------------------------
subplot(3,3,2)
hold on
load cf; a = dat(:,1)*1000; b = dat(:,2)/1000;c = b<.766;
% load plotdatacf; a = dat(:,5)*1000; b = dat(:,1)/1000;
h=area(a,b);
h.FaceColor=[.8 .8 .8];
plot(a(c),b(c),'r')
text(100,.94,'Tidal Flat Dominance','fontsize',8)
text(450,.735,'Marsh Dominance','fontsize',8)
ht=text(510,.786,'Critical Fetch Length','fontsize',7,'FontAngle', 'italic');
set(ht,'Rotation',59);

ylim([.7 1])
xlabel('River Concentration (mg/l)')
ylabel('Fetch Length (km)')
box on
NumTicks = 5;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))

%------------------------------------------------
subplot(3,3,3)
hold on
load plotdataq; a = dat(:,7); b = dat(:,1)/1000; c = dat(:,2)<=.7;
h=area(a,b);
h.FaceColor=[.8 .8 .8];
plot(a(c),b(c),'r')
text(45,.45,'Tidal Flat Dominance','fontsize',8)
text(45,.3,'Marsh Dominance','fontsize',8)
ht=text(70,.412,'Critical Fetch Length','fontsize',7,'FontAngle', 'italic');
set(ht,'Rotation',-16.5);

xlim([10 150])
ylim([.25 .5])
xlabel('River Discharge (m^3/s)')
ylabel('Fetch Length (km)')
box on
NumTicks = 5;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))

%------------------------------------------------
subplot(3,3,4)
hold on
load plotdatal; a = dat(:,6)/1000; b = dat(:,1)/1000;  c = dat(:,2)<=.7;
h=area(a,b);
h.FaceColor=[.8 .8 .8];
plot(a(c),b(c),'r')
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
subplot(3,3,5)
hold on
load plotdatab; a = dat(:,11)/1000; b = dat(:,1)/1000; c = dat(:,2)<=.7;
h=area(a,b);
h.FaceColor=[.8 .8 .8];
plot(a(c),b(c),'r')
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
subplot(3,3,6)
hold on
load plotdatar; a = dat(:,8)/10^-3*365*24*60*60; b = dat(:,1)/1000; c = dat(:,2)<=.7;
h=area(a,b);
h.FaceColor=[.8 .8 .8];
plot(a(c),b(c),'r')
text(3,.475,'Tidal Flat Dominance','fontsize',8)
text(3,.3,'Marsh Dominance','fontsize',8)
ht=text(1.2,.44,'Critical Fetch Length','fontsize',7,'FontAngle', 'italic');
set(ht,'Rotation',0);

ylim([.25 .5])
xlim([1 10])
xlabel('Rate of Sea Level Rise (mm/yr)')
ylabel('Fetch Length (km)')
box on

%------------------------------------------------
subplot(3,3,7)
hold on
load plotdatah; a = dat(:,9)*2; b = dat(:,1)/1000; c = dat(:,2)<=.7;
h=area(a,b);
plot(a(c),b(c),'r')
h.FaceColor=[.8 .8 .8];
text(3,.475,'Tidal Flat Dominance','fontsize',8)
text(3,.3,'Marsh Dominance','fontsize',8)
ht=text(5,.448,'Critical Fetch Length','fontsize',7,'FontAngle', 'italic');
set(ht,'Rotation',2);

xlim([1 10])
ylim([.25 .5])
xlabel('Tidal Range (m)')
ylabel('Fetch Length (km)')
box on
% NumTicks = 5;
% L = get(gca,'XLim');
% set(gca,'XTick',linspace(L(1),L(2),NumTicks))

%------------------------------------------------
subplot(3,3,8)
hold on
load plotdatat; a = dat(:,10)/60/60; b = dat(:,1)/1000; c = dat(:,2)<=.7;
h=area(a,b);
h.FaceColor=[.8 .8 .8];
plot(a(c),b(c),'r')
text(15,.475,'Tidal Flat Dominance','fontsize',8)
text(15,.3,'Marsh Dominance','fontsize',8)
ht=text(17,.433,'Critical Fetch Length','fontsize',7,'FontAngle', 'italic');
set(ht,'Rotation',-4);

xlim([12 24])
ylim([.25 .5])
xlabel('Tidal Period (hr)')
ylabel('Fetch Length (km)')
box on
NumTicks = 2;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))

%------------------------------------------------
subplot(3,3,9)
hold on
load plotdatav; a = dat(:,12); b = dat(:,1)/1000; c = dat(:,2)<=.7;
h=area(a,b);
h.FaceColor=[.8 .8 .8];
plot(a(c),b(c),'r')
text(6,4,'Tidal Flat Dominance','fontsize',8)
text(1.2,.9,{'Marsh','Dominance'},'fontsize',8)
ht=text(8,.35,'Critical Fetch Length','fontsize',7,'FontAngle', 'italic');
set(ht,'Rotation',-1.5);

% ylim([.25 .5])
xlim([1 15])
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
tit='fetchvs9v_2';
print(tit,'-dtiff','-r400')
movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\24 - Final Results')
close all