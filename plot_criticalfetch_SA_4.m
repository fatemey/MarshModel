function plot_criticalfetch_SA_4
% plots parameters for sensitivity analysis purposes
%--------------------------------------------------------------------------------------------------
format compact
format longG
clear
clf

load SA_results_EET

x1 = stat(:,1);
x2 = stat(:,2);
y = sqrt(stat(:,3));

% x1 = log(abs(stat(:,1))).*sign(stat(:,1));
% x2 = log(abs(stat(:,2)));
% y = sqrt(log(abs(stat(:,3))));

c = {'C_o';'C_r';'Q_r';'l_E';'b_B';'H';'R';'V_w'};

%------------------------------------------------
subplot(1,2,1)
gscatter(x1,y,c);
legend off

xlabel('Mean of EEs')
ylabel('Standard Deviation of EEs')
box on

%------------------------------------------------
subplot(1,2,2)
gscatter(x2,y,c)

xlabel('Mean of |EEs|')
ylabel('Standard Deviation of EEs')
box on

%------------------------------------------------
set(findobj('type','axes'),'fontsize',10)
h_fig=gcf;
set(h_fig,'PaperOrientation','portrait')
set(h_fig,'PaperPosition', [0 0 6.5 3]) % [... ... max_width=7.5 max_height=9]
tit='SA_4';
print(tit,'-dtiff','-r400')
movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\25 - Sensitivity Results')
close all

end