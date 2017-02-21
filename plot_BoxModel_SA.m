function plot_BoxModel_SA(t,y,leg,tit)
% Plots results vs time
% The plot includes 4 graphs (C_r, b_f, d_f, d_m) for diffrenet sets of parameters.
%
%--------------------------------------------------------------------------
C_r = squeeze(y(:,4,:)); b_f = squeeze(y(:,1,:))/1000; d_f = squeeze(y(:,2,:)); d_m = squeeze(y(:,3,:));

left=.15; midh=.15; right=.05;
up=.1; midv=.05; bottom=left;
width=(1-left-midh-right)/2;
height=(1-up-midv-bottom)/2;
left1=left; left2=left+width+midh;
bottom1=midv+height+bottom; bottom2=bottom;

title(tit)
axis off

axes('Position',[left1 bottom1 width height]);
plot(t,C_r*1000,'linewidth',2)
ylabel('SSC (mg/l)')
set(gca,'XTickLabel',[])
lgnd = legend(leg,'location','best','fontsize',10);
set(lgnd,'color','none')

axes('Position',[left2 bottom1 width height]);
plot(t,b_f,'linewidth',2)
ylabel('Tidal Flat Width (km)')
set(gca,'XTickLabel',[])

axes('Position',[left1 bottom2 width height]);
plot(t,d_m,'linewidth',2)
xlabel('Year')
ylabel('Marsh Depth (m)')

axes('Position',[left2 bottom2 width height]);
plot(t,d_f,'linewidth',2)
xlabel('Year')
ylabel('Tidal Flat Depth (m)')

box on
set(findobj('type','axes'),'fontsize',15)

h_fig=gcf;
set(h_fig,'PaperOrientation','portrait')
set(h_fig,'PaperPosition', [0 0 7.5 6]) % [... ... max_width=7.5 max_height=9]