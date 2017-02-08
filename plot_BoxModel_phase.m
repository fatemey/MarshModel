function plot_BoxModel_phase
% plot_BoxModel_phase: makes an animation of marsh and tidal flat depth and
% width changes in time
%
%--------------------------------------------------------------------------------------------------
clear
clf

%-------------- Set data values
% load data_base; tit = 'base data';
% load data_highcfqf; tit = 'highcfqf';
% load data_nowind; tit = 'no wind';
load data_lowwind; tit = 'low wind';

tyr = 1 : ceil(t(end)/365/24/60/60);
ts = tyr*365*24*60*60;

for i = 1 : length(ts)
    tyr_ind(i)=find (t==ts(i));
end

bf = y(tyr_ind, 1);
df = y(tyr_ind, 2);
dm = y(tyr_ind, 3);

b_fm = 5/2 *10^3;   % total basin width (both sides of the channel) (m)
H = 1.4/2;          % tidal range (m)

%-------------- PLot feature initiations
c1=[160,82,45]/255;     % brown
c2=[60,179,113]/255;    % green
c3=[173,223,255]/255;   % blue
axis([0 1 0 1.2])
set(gca,'XTickLabel',[])
set(gca,'YTickLabel',[])
set(gca,'XTick',[])
set(gca,'YTick',[])
% box on
axis off
hold on

for i = 1 : 500 %length(tyr)
    
    %-------------- Making plots
    h3=area([0,1],[1,1]); % plots marsh surface
    h1=area([0,bf(i)/b_fm],[1-df(i)/(2*H),1-df(i)/(2*H)]); % plots tidal flat surface
    h2=area([bf(i)/b_fm,1],[1-dm(i)/(2*H),1-dm(i)/(2*H)]); % plots marsh surface
    
    plot([0,1],[1,1],'--','color','k','LineWidth',2) % plots MHWL
    plot([0,1],[1/2,1/2],'--','color','k','LineWidth',2) % plots MSWL
    plot([0,1],[0,0],'--','color','k','LineWidth',2) % plots MLWL
%     plot([1,1],[0,1],'color','k','LineWidth',.5) % plots a vertical line on the right side of the graph
    
    h1.FaceColor = c1;
    h2.FaceColor = c2;
    h3.FaceColor = c3;
    
    text(.01,.98,'MHWL','fontsize',20)
    text(.01,.48,'MSL','fontsize',20)
    text(.01,.02,'MLWL','fontsize',20)
    
    caption = sprintf('Year = %d',i);
%     title(caption,'fontsize',25);
    text(.45,1.1,caption,'fontsize',25)
    
    %-------------- Start animation
    drawnow
    F(i) = getframe;
    children = get(gca, 'children');
    delete(children);
    
end

%-------------- Save animation
v = VideoWriter([tit,'.avi']);
v.FrameRate = 20;  % Default 30
v.Quality = 100;    % Default 75
open(v);
writeVideo(v,F)
close(v);
movefile([tit,'.avi'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\13 - Animation results')
% close all

end
