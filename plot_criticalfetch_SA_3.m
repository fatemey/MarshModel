function plot_criticalfetch_SA_3
% plots parameters for sensitivity analysis purposes
%--------------------------------------------------------------------------------------------------
format compact
format longG
clear
clf

% load data_CriticalFetch_SA
load bf_SA
y = dat;
load bfm_SA
y(y>10000&dat==10000) = 10000; 
y(y>5000&dat==5000) = 5000; 
y(y>1000&dat==1000) = 1000;
y = y/1000;

dat = csvread('SA_DATA_.csv');
data = [dat(:,4:11),y];

%------------------------------------------------
% labels = {'Sepal Length','Sepal Width','Petal Length','Petal Width'};
% parallelcoords(meas,'Group',species,'Labels',labels)

 parallelcoords(data)
 
% xlabel('Ocean Concentration (mg/l)')
% ylabel('Fetch Length (km)')
box on

%------------------------------------------------
set(findobj('type','axes'),'fontsize',10)
h_fig=gcf;
set(h_fig,'PaperOrientation','portrait')
set(h_fig,'PaperPosition', [0 0 10 9]) % [... ... max_width=7.5 max_height=9]
tit='SA_3';
print(tit,'-dtiff','-r400')
movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\25 - Sensitivity Results')
close all

end