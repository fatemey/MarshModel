function plot_results_BoxModel_SA(t,y,leg,tit)
% Plot results vs time
%
%--------------------------------------------------------------------------
C_r = squeeze(y(:,4,:)); b_f = squeeze(y(:,1,:))/1000; d_f = squeeze(y(:,2,:)); d_m = squeeze(y(:,3,:));

left=.15; midh=.15; right=.05;
up=.1; midv=.05; bottom=left;
width=(1-left-midh-right)/2;
height=(1-up-midv-bottom)/2;
left1=left; left2=left+width+midh;
bottom1=midv+height+bottom; bottom2=bottom;

tyr = ceil(t(end)/365/24/60/60);
time=linspace(0,t(end),6);
timec=num2cell(linspace(0,tyr,6));

title(tit)
axis off

axes('Position',[left1 bottom1 width height]);
plot(t,C_r*1000,'linewidth',2)
% scatter(t,C_r(:,2)*1000,'.')
ylabel('SSC (g/m^3)')
set(gca,'XTick',time)
set(gca,'XTickLabel',[])
xlim([0 t(end)])
lgnd = legend(leg,'location','best','fontsize',10);
set(lgnd,'color','none')
% NumTicks = 4;
% L = get(gca,'YLim');
% set(gca,'YTick',linspace(L(1),L(2),NumTicks))

axes('Position',[left2 bottom1 width height]);
plot(t,b_f,'linewidth',2)
% scatter(t,b_f(:,2),'.')
ylabel('Tidal Flat Width (km)')
set(gca,'XTick',time)
set(gca,'XTickLabel',[])
xlim([0 t(end)])

axes('Position',[left1 bottom2 width height]);
plot(t,d_m,'linewidth',2)
% scatter(t,d_m(:,2),'.')
xlabel('Year')
ylabel('Marsh Depth (m)')
set(gca,'XTick',time)
set(gca,'XTickLabel',timec)
xlim([0 t(end)])

axes('Position',[left2 bottom2 width height]);
plot(t,d_f,'linewidth',2)
% scatter(t,d_f(:,2),'.')
xlabel('Year')
ylabel('Tidal Flat Depth (m)')
set(gca,'XTick',time)
set(gca,'XTickLabel',timec)
xlim([0 t(end)])

box on
set(findobj('type','axes'),'fontsize',15)

h_fig=gcf;
set(h_fig,'PaperOrientation','portrait')
set(h_fig,'PaperPosition', [0 0 7.5 6]) % [... ... max_width=7.5 max_height=9]