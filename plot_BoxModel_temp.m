function plot_BoxModel_temp(t,y)
% Plots results vs time
% The plot includes 4 graphs (C_r, b_f, d_f, d_m) for 1 set of parameters.
%
%--------------------------------------------------------------------------
C_r = y(:,4); b_f = y(:,1)/1000; d_f = y(:,2); d_m = y(:,3);

left=.15; midh=.15; right=.05;
up=.1; midv=.05; bottom=left;
width=(1-left-midh-right)/2;
height=(1-up-midv-bottom)/2;
left1=left; left2=left+width+midh;
bottom1=midv+height+bottom; bottom2=bottom;

% tyr = ceil(t(end)/365/24/60/60);
% time=linspace(0,t(end),6);
% timec=num2cell(linspace(0,tyr,6));

axis off

axes('Position',[left1 bottom1 width height]);
plot(t,C_r*1000,'linewidth',2)
% scatter(t,C_r*1000,'o')
ylabel('SSC (g/m^3)')
% set(gca,'XTick',time)
set(gca,'XTickLabel',[])
% xlim([0 t(end)])

axes('Position',[left2 bottom1 width height]);
plot(t,b_f,'linewidth',2)
% scatter(t,b_f,'.')
ylabel('Tidal Flat Width (km)')
% set(gca,'XTick',time)
set(gca,'XTickLabel',[])
% xlim([0 t(end)])

axes('Position',[left1 bottom2 width height]);
plot(t,d_m,'linewidth',2)
% scatter(t,d_m,'.')
xlabel('Year')
ylabel('Marsh Depth (m)')
% set(gca,'XTick',time)
% set(gca,'XTickLabel',timec)
% xlim([0 t(end)])

axes('Position',[left2 bottom2 width height]);
plot(t,d_f,'linewidth',2)
% scatter(t,d_f,'.')
xlabel('Year')
ylabel('Tidal Flat Depth (m)')
% set(gca,'XTick',time)
% set(gca,'XTickLabel',timec)
% xlim([0 t(end)])

box on
set(findobj('type','axes'),'fontsize',15)

h_fig=gcf;
set(h_fig,'PaperOrientation','portrait')
set(h_fig,'PaperPosition', [0 0 7.5 6]) % [... ... max_width=7.5 max_height=9]