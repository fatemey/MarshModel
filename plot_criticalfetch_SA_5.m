function plot_criticalfetch_SA_5
% plots parameters for sensitivity analysis purposes
%--------------------------------------------------------------------------------------------------
format compact
format longG
clear
clf

load SA_results_EET

% x1 = stat(:,1);
% x2 = stat(:,2);

x1 = log(abs(stat(:,1))).*sign(stat(:,1));
x2 = log(abs(stat(:,2)));

%------------------------------------------------
subplot(1,2,1)
bar(x1)

par_name= {'C_o';'C_r';'Q_r';'l_E';'b_B';'H';'R';'V_w'};
set(gca,'xticklabel',par_name)
xlabel('Parameters')
ylabel('Sensitivity (\mu)')
box on

%------------------------------------------------
subplot(1,2,2)
bar(x2)

set(gca,'xticklabel',par_name)
xlabel('Parameters')
ylabel('Sensitivity (\mu*)')
box on

%------------------------------------------------
set(findobj('type','axes'),'fontsize',10)
h_fig=gcf;
set(h_fig,'PaperOrientation','portrait')
set(h_fig,'PaperPosition', [0 0 8 3]) % [... ... max_width=7.5 max_height=9]
tit='SA_5_loglog';
print(tit,'-dtiff','-r400')
movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\25 - Sensitivity Results')
close all

end