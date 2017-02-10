subplot(3,3,1)
a=[20 40 60 80 100]; b = [220 410 610 840 1110]; 
scatter(a,b,'k','o','filled')
xlabel('Ocean Concentration (g/L)')
ylabel('Threshold Width (m)')
box on
NumTicks = 6;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))

subplot(3,3,2)
a=0:20:100; b=215:5:240;
scatter(a,b,'k','o','filled')
xlabel('River Concentration (g/L)')
ylabel('Threshold Width (m)')
box on
NumTicks = 6;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))

subplot(3,3,3)
a=0:20:100; b =[220 220 220 220 215 215];
scatter(a,b,'k','o','filled')
xlabel('River Discharge (m^3/s)')
ylabel('Threshold Width (m)')
box on
NumTicks = 6;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))

subplot(3,3,4)
% a=5:5:25; b=ones(length(a),1)*220;
a=[.5,1,5:5:20]; b=[185;200;220;220;220;220];
scatter(a,b,'k','o','filled')
xlabel('Estuary Length (km)')
ylabel('Threshold Width (m)')
box on
NumTicks = 5;
% xlim([5 25])
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))

subplot(3,3,5)
a=[1,5:5:20]; b=[165 220 240 245 250];
scatter(a,b,'k','o','filled')
xlabel('Basin Width (km)')
ylabel('Threshold Width (m)')
box on
NumTicks = 5;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))

subplot(3,3,6)
a=[0 2 4 6 8]; b=[2500,2160,470,200,120]; 
scatter(a,b,'k','o','filled')
xlabel('Wind Speed (m/s)')
ylabel('Threshold Width (m)')
box on
NumTicks = 6;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))

subplot(3,3,7)
a=2:2:10; b=ones(length(a),1)*220;
scatter(a,b,'k','o','filled')
xlabel('Rate of Sea Level Rise (mm/yr)')
ylabel('Threshold Width (m)')
box on
NumTicks = 6;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))

subplot(3,3,8)
% a=[1 1.4 1.8]; b=[220 220 220];
a=[1 1.4 1.8]; b=[225 220 210];
scatter(a,b,'k','o','filled')
xlabel('Tidal Range (m)')
ylabel('Threshold Width (m)')
box on
NumTicks = 6;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))
xlim([1 2])

subplot(3,3,9)
a=[.1 .3 .5]; b=[0 220 220];
scatter(a,b,'k','o','filled')
xlabel('Critical Shear Stress (PA)')
ylabel('Threshold Width (m)')
box on
NumTicks = 6;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))

set(findobj('type','axes'),'fontsize',10)

h_fig=gcf;
set(h_fig,'PaperOrientation','portrait')
set(h_fig,'PaperPosition', [0 0 7.5 6]) % [... ... max_width=7.5 max_height=9]

tit='threshold';
print(tit,'-dtiff','-r400')
movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\13 - Animation results')
close all