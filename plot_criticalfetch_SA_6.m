function plot_criticalfetch_SA_6
% plots parameters for sensitivity analysis purposes
%--------------------------------------------------------------------------------------------------
format compact
format longG
clear
clf
hold on

load SA_results_EET

x = mean(stat(:,3:4),2);
y = mean(stat(:,5:6),2);
errh = stat(:,4)-stat(:,3);
errv = stat(:,6)-stat(:,5);

c = {'C_o';'C_r';'Q_r';'l_E';'b_B';'H';'R';'V_w'};
errorbar(x,y,errv,'k.')
errorbar(x,y,errh,'horizontal','k.')
gscatter(x,y,c);

% legend off

xlabel('Mean of EEs')
ylabel('Standard Deviation of EEs')
box on


%------------------------------------------------
set(findobj('type','axes'),'fontsize',10)
h_fig=gcf;
set(h_fig,'PaperOrientation','portrait')
set(h_fig,'PaperPosition', [0 0 3.2 3.2]) % [... ... max_width=7.5 max_height=9]
tit='SA_6';
print(tit,'-dtiff','-r400')
movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\25 - Sensitivity Results')
close all

end