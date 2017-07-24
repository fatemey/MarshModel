function plot_nondpars
% plots nondimensional variables based on dat_pars.
%
%--------------------------------------------------------------------------------------------
format compact
format longG
clear

%-------------- Set variables
load dat_pars
X_c = dat(:,1);
C_o = dat(:,2);
C_f = dat(:,3);
Q_f = dat(:,4);
b_fm = dat(:,5);
L_E = dat(:,6);
A = b_fm.*L_E;
R = 2 *10^-3/365/24/60/60;  
H = 1.4/2;
T = 12 *60*60; 
rho_s = 1000;   % sediment bulk density (kg/m3)

%-------------- Plot 1
clf
scatter(C_o./C_f,X_c./b_fm,'k','o','filled')
xlabel('Ocean Concentration to River Concentration')
ylabel('Critical Initial Width to Basin Width')
box on
set(findobj('type','axes'),'fontsize',15)
h_fig=gcf;
set(h_fig,'PaperOrientation','portrait')
set(h_fig,'PaperPosition', [0 0 7.5 6]) % [... ... max_width=7.5 max_height=9]
tit='plot 1';
print(tit,'-dtiff','-r400')
movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\17 - Nondimensional analysis for critical fetch')

%-------------- Plot 1-2
clf
scatter(C_o./C_f,X_c./L_E,'k','o','filled')
xlabel('Ocean Concentration to River Concentration')
ylabel('Critical Initial Width to Estuary Length')
box on
set(findobj('type','axes'),'fontsize',15)
h_fig=gcf;
set(h_fig,'PaperOrientation','portrait')
set(h_fig,'PaperPosition', [0 0 7.5 6]) % [... ... max_width=7.5 max_height=9]
tit='plot 1-2';
print(tit,'-dtiff','-r400')
movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\17 - Nondimensional analysis for critical fetch')

%-------------- Plot 2
clf
scatter(R.*A./Q_f,X_c./b_fm,'k','o','filled')
xlabel('Sea Level Rise Rate x Area to River Discharge')
ylabel('Critical Initial Width to Basin Width')
box on
set(findobj('type','axes'),'fontsize',15)
h_fig=gcf;
set(h_fig,'PaperOrientation','portrait')
set(h_fig,'PaperPosition', [0 0 7.5 6]) % [... ... max_width=7.5 max_height=9]
tit='plot 2';
print(tit,'-dtiff','-r400')
movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\17 - Nondimensional analysis for critical fetch')

%-------------- Plot 2-2
clf
scatter(R.*A./Q_f,X_c./L_E,'k','o','filled')
xlabel('Sea Level Rise Rate x Area to River Discharge')
ylabel('Critical Initial Width to Estuary Length')
box on
set(findobj('type','axes'),'fontsize',15)
h_fig=gcf;
set(h_fig,'PaperOrientation','portrait')
set(h_fig,'PaperPosition', [0 0 7.5 6]) % [... ... max_width=7.5 max_height=9]
tit='plot 2-2';
print(tit,'-dtiff','-r400')
movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\17 - Nondimensional analysis for critical fetch')

%-------------- Plot 3/4
clf
var = X_c./sqrt(A);
n_clusters = 4;
circle_siz = ones(size(var));
circle_col = zeros(size(var));
col = zeros(length(var),3);
siz_step = linspace(0,max(var),n_clusters);
for i = 1 : n_clusters-1
    circle_siz(var>siz_step(i)&var<=siz_step(i+1)) = i*100;
    circle_col(var>siz_step(i)&var<=siz_step(i+1)) = i;
end
col(circle_col==0,:) = repmat([255,165,0]/255,sum(circle_col==0),1); % gold
col(circle_col==1,:) = repmat([240,128,128]/255,sum(circle_col==1),1); % light coral
col(circle_col==2,:) = repmat([50,205,50]/255,sum(circle_col==2),1); % lime green
col(circle_col==3,:) = repmat([0,0,205]/255,sum(circle_col==3),1); % medium blue
col(circle_col==4,:) = repmat([255,0,255]/255,sum(circle_col==4),1); % magenta 
col(circle_col==5,:) = repmat([112,128,144]/255,sum(circle_col==5),1); % slate gray

scatter(C_o./C_f,R.*A./Q_f,circle_siz,'k','o')
title('\chi_{c}/{\surd}A')
xlabel('C_{o}/C_{f}')
ylabel('R.A/Q_{f}')
box on
set(findobj('type','axes'),'fontsize',15)
h_fig=gcf;
set(h_fig,'PaperOrientation','portrait')
set(h_fig,'PaperPosition', [0 0 7.5 6]) % [... ... max_width=7.5 max_height=9]
tit='plot 3';
% print(tit,'-dtiff','-r400')
% movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\17 - Nondimensional analysis for critical fetch')

for i =1 : length(var)
    scatter(C_o(i)/C_f(i),R*A(i)/Q_f(i),circle_siz(i),'MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:),'MarkerFaceAlpha',.4,'MarkerEdgeAlpha',.4)
    hold on
end
title('\chi_{c}/{\surd}A')
xlabel('C_{o}/C_{f}')
ylabel('R.A/Q_{f}')
box on
set(findobj('type','axes'),'fontsize',15)
h_fig=gcf;
set(h_fig,'PaperOrientation','portrait')
set(h_fig,'PaperPosition', [0 0 7.5 6]) % [... ... max_width=7.5 max_height=9]
tit='plot 4';
% print(tit,'-dtiff','-r400')
% movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\17 - Nondimensional analysis for critical fetch')

%-------------- Plot 5
clf
scatter(C_o*H.*L_E.*b_fm/T./C_f./Q_f, X_c./b_fm,'k','o','filled')
xlabel('Ocean to River Input')
ylabel('Critical Initial Width to Basin Width')
box on
set(findobj('type','axes'),'fontsize',15)
h_fig=gcf;
set(h_fig,'PaperOrientation','portrait')
set(h_fig,'PaperPosition', [0 0 7.5 6]) % [... ... max_width=7.5 max_height=9]
tit='plot 5';
print(tit,'-dtiff','-r400')
movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\17 - Nondimensional analysis for critical fetch')

%-------------- Plot 6
clf
scatter(C_o*H.*L_E.*b_fm/T./C_f./Q_f, X_c./L_E,'k','o','filled')
xlabel('Ocean to River Input')
ylabel('Critical Initial Width to Estuary Length')
box on
set(findobj('type','axes'),'fontsize',15)
h_fig=gcf;
set(h_fig,'PaperOrientation','portrait')
set(h_fig,'PaperPosition', [0 0 7.5 6]) % [... ... max_width=7.5 max_height=9]
tit='plot 6';
print(tit,'-dtiff','-r400')
movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\17 - Nondimensional analysis for critical fetch')

%-------------- Plot 7
clf
scatter(((C_o*H.*L_E.*b_fm/T)+C_f.*Q_f)./R./L_E./b_fm./C_o, X_c./L_E,'k','o','filled')
xlabel('River + Ocean Input to Input for Created Space due to SLR ')
ylabel('Critical Initial Width to Estuary Length')
box on
set(findobj('type','axes'),'fontsize',15)
h_fig=gcf;
set(h_fig,'PaperOrientation','portrait')
set(h_fig,'PaperPosition', [0 0 7.5 6]) % [... ... max_width=7.5 max_height=9]
tit='plot 7';
print(tit,'-dtiff','-r400')
movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\17 - Nondimensional analysis for critical fetch')

%-------------- Plot 8
clf
scatter(((C_o*H.*L_E.*b_fm/T)+C_f.*Q_f)./R./L_E./b_fm./C_o, X_c./b_fm,'k','o','filled')
xlabel('River + Ocean Input to Input for Created Space due to SLR ')
ylabel('Critical Initial Width to Basin Width')
box on
set(findobj('type','axes'),'fontsize',15)
h_fig=gcf;
set(h_fig,'PaperOrientation','portrait')
set(h_fig,'PaperPosition', [0 0 7.5 6]) % [... ... max_width=7.5 max_height=9]
tit='plot 8';
print(tit,'-dtiff','-r400')
movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\17 - Nondimensional analysis for critical fetch')

%-------------- Plot 9
clf
scatter((H.*L_E.*b_fm+T*Q_f)./R./L_E./b_fm./T, X_c./L_E,'k','o','filled')
xlabel('River + Ocean Water Input to Water Input for Created Space due to SLR ')
ylabel('Critical Initial Width to Estuary Length')
box on
set(findobj('type','axes'),'fontsize',15)
h_fig=gcf;
set(h_fig,'PaperOrientation','portrait')
set(h_fig,'PaperPosition', [0 0 7.5 6]) % [... ... max_width=7.5 max_height=9]
tit='plot 9';
print(tit,'-dtiff','-r400')
movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\17 - Nondimensional analysis for critical fetch')

%-------------- Plot 10
clf
scatter((H.*L_E.*b_fm+T*Q_f)./R./L_E./b_fm./T, X_c./b_fm,'k','o','filled')
xlabel('River + Ocean Water Input to Water Input for Created Space due to SLR ')
ylabel('Critical Initial Width to Basin Width')
box on
set(findobj('type','axes'),'fontsize',15)
h_fig=gcf;
set(h_fig,'PaperOrientation','portrait')
set(h_fig,'PaperPosition', [0 0 7.5 6]) % [... ... max_width=7.5 max_height=9]
tit='plot 10';
print(tit,'-dtiff','-r400')
movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\17 - Nondimensional analysis for critical fetch')

close all
end