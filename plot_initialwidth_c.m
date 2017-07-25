function plot_initialwidth_c
% Function plot_initialwidth plots critical initial width values (computed from functin CriticalInitialWidth) vs different
% parameters.
%
% Last Update: 7/24/2017
%
%--------------------------------------------------------------------------------------------------
format compact
format longG
clear
% clf

%------------------------------------------------
load co_data; a = dat(:,1)*1000; b = dat(:,2); dtf = dat(:,3); dm = dat(:,4); eqdtf = dat(:,5); eqdm = dat(:,6); conv = dat(:,7);

subplot(1,2,1)
gscatter(a,b,conv,'krcy',[],[],'off')
xlabel('Ocean Concentration (mg/l)')
ylabel('Initial Width (m)') 
box on

subplot(1,2,2)
hold on
gscatter(a,dtf,eqdtf,'rb',[],[],'off')
gscatter(a,dm,eqdm,'rg',[],[],'off')
xlabel('Ocean Concentration (mg/l)')
ylabel('Depth (m)') 
box on
% NumTicks = 5;
% L = get(gca,'XLim');
% set(gca,'XTick',linspace(L(1),L(2),NumTicks))


%------------------------------------------------
% hL = legend([l1,l2,l3,l4],{'Normal Expansion/Contraction','TF Expansion & Emergence','TF Contraction & Emergence','Marsh Drowning'});
% hL = legend([l1,l2,l3,l4,l5],{'Normal Expansion/Contraction','TF Expansion & Emergence','TF Contraction & Emergence','Marsh Expansion & Drowning','Marsh Contraction & Drowning'});
% newPosition = [0.7 0.15 0.2 0.2];
% newUnits = 'normalized';
% set(hL,'Position', newPosition,'Units', newUnits);
% legend boxoff 

%------------------------------------------------
% set(findobj('type','axes'),'fontsize',10)
% h_fig=gcf;
% set(h_fig,'PaperOrientation','portrait')
% set(h_fig,'PaperPosition', [0 0 7.5 4]) % [... ... max_width=7.5 max_height=9]
% tit='Initial Width_largeCfQf';
% print(tit,'-dtiff','-r400')
% movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\14 - Fetch threshold')
% close all