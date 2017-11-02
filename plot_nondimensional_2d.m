%-------------- Plot nondimensional analysis results
% loading data
clear

C_o = 20 *10^-3;    % ocean concertation (kg/m3)
C_f = 15 *10^-3;    % river concentration (kg/m3)
Q_f = 20/2;         % river water discharge (m3/s)
b_fm = 5 *10^3; % total basin width (both sides of the channel) (m)
L_E = 5 *10^3; % basin length (m)
R = 2 *10^-3/365/24/60/60;   % sea level rise (m/s)
T_T = 12 *60*60;   % tidal period (s) ( = 12 hours)
H = 1.4/2;          % tidal amplitude (range/2) (m)
rho_s = 1000;   % sediment bulk density (kg/m3)
par_v = [1/3, 2/3, 1, 4/3];

%%
Co_ = C_o;
n = length(Co_);
k = length (par_v);

load Sol_bfm_m
bfm_ = b_fm * par_v;

for i = 1 : n
    dat(:,1,i) = squeeze(Sol(i,1,:))/sqrt(Q_f/R);
    dat(:,2,i) = L_E*bfm_*H/T_T/Q_f;
    dat(:,3,i) = squeeze(Sol(i,1,:))./bfm_';
    dat(:,4,i) = Co_(i)*L_E*bfm_*H/T_T/Q_f/C_f;
    dat(:,5,i) = L_E*bfm_*H/Q_f^1.5*R^-1.5;
    dat(:,6,i) = (Co_(i)*(L_E*bfm_*H/T_T-Q_f)+C_f*Q_f)/R./bfm_/L_E/rho_s;
    dat(:,7,i) = squeeze(Sol(i,1,:));
    dat(:,8,i) = squeeze(Sol(i,2,:));
    dat(:,9,i) = squeeze(Sol(i,3,:));
end
pi_bfm_m = dat;
clear dat

load Sol_le_m
LE_ = L_E * par_v;

for i = 1 : n
    dat(:,1,i) = squeeze(Sol(i,1,:))/sqrt(Q_f/R);
    dat(:,2,i)  = LE_*b_fm*H/T_T/Q_f;
    dat(:,3,i) = squeeze(Sol(i,1,:))/b_fm;
    dat(:,4,i)  = Co_(i)*LE_*b_fm*H/T_T/Q_f/C_f;
    dat(:,5,i) = LE_*b_fm*H/Q_f^1.5*R^-1.5;
    dat(:,6,i) = (Co_(i)*(LE_*b_fm*H/T_T-Q_f)+C_f*Q_f)/R/b_fm./LE_/rho_s;
    dat(:,7,i) = squeeze(Sol(i,1,:));
    dat(:,8,i) = squeeze(Sol(i,2,:));
    dat(:,9,i) = squeeze(Sol(i,3,:));
end
pi_le_m = dat;
clear dat

R_ = R;
n = length(R_);

load Sol_bfm_d
bfm_ = b_fm ./ par_v;

for i = 1 : n
    dat(:,1,i) = squeeze(Sol(i,1,:))./sqrt(bfm_'*L_E);
    dat(:,2,i) = C_f*Q_f/C_o/R_(i)/L_E./bfm_;
    dat(:,3,i) = squeeze(Sol(i,1,:));
    dat(:,4,i) = squeeze(Sol(i,2,:));
    dat(:,5,i) = squeeze(Sol(i,3,:));
end
pi_bfm_d = dat;
clear dat

load Sol_le_d
LE_ = L_E ./ par_v;


for i = 1 : n
    dat(:,1,i) = squeeze(Sol(i,1,:))./sqrt(b_fm*LE_');
    dat(:,2,i) = C_f*Q_f/C_o/R_(i)./LE_/b_fm;
    dat(:,3,i) = squeeze(Sol(i,1,:));
    dat(:,4,i) = squeeze(Sol(i,2,:));
    dat(:,5,i) = squeeze(Sol(i,3,:));
end
pi_le_d = dat;
clear dat

load Sol_rho
rho_ = rho_s ./ par_v;


for i = 1 : n
    dat(:,1,i) = squeeze(Sol(i,1,:))/b_fm;
    dat(:,2,i) = (C_o*(b_fm*L_E*H/T_T-Q_f)+C_f*Q_f)/R/b_fm/L_E./rho_;
    dat(:,3,i) = squeeze(Sol(i,1,:));
    dat(:,4,i) = squeeze(Sol(i,2,:));
    dat(:,5,i) = squeeze(Sol(i,3,:));
end
pi_rho = dat;
clear dat

load Sol_r
R_ = R ./ par_v;

for i = 1 : n
    dat(:,1,i) = squeeze(Sol(i,1,:))/b_fm;
    dat(:,2,i) = (C_o*(b_fm*L_E*H/T_T-Q_f)+C_f*Q_f)./R_/b_fm/L_E/rho_s;
    dat(:,3,i) = squeeze(Sol(i,1,:));
    dat(:,4,i) = squeeze(Sol(i,2,:));
    dat(:,5,i) = squeeze(Sol(i,3,:));
end
pi_r = dat;
clear dat

%%
% chi_2
clf

subplot(1,2,1)

hold on
for i = 1 : n
    plot(pi_bfm_m(:,2,i),pi_bfm_m(:,1,i),'Linewidth',2);
end
xlabel('\pi_1 = Tidal Prism Ratio (Q_O/Q_R) , b_B')
ylabel('\pi =  Nondimensional Width (\chi/{\surd} Q_R/R)')
title('\pi_2 = Concentration Ratio Constant')
box on

subplot(1,2,2)

hold on
for i = 1 : n
    plot(pi_le_m(:,2,i),pi_le_m(:,1,i),'Linewidth',2);
end
xlabel('\pi_1 = Tidal Prism Ratio (Q_O/Q_R) , l_E')
ylabel('\pi =  Nondimensional Width (\chi/{\surd} Q_R/R)')
title('\pi_2 = Concentration Ratio Constant')
box on

% set(findobj('type','axes'),'fontsize',10)
% h_fig = gcf;
% set(h_fig,'PaperOrientation','portrait')
% set(h_fig,'PaperPosition', [0 0 8 3]) % [... ... max_width = 7.5 max_height = 9]
% tit = 'phi-2';
% print(tit,'-dtiff','-r400')
% movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\22 - Nondimensional Analysis')
% close all

%%
% chi_3
clf

subplot(1,2,1)

hold on
for i = 1 : n
    plot(pi_bfm_d(:,2,i),pi_bfm_d(:,1,i),'Linewidth',2);
end
xlabel('\pi_1 = Mass Ratio (C_R Q_R / C_O Q_{SLR}) , b_B')
ylabel('\pi = Nondimensional Width (\chi/\surd b_B l_E)')
title('\pi_2 = Water Rise Rate Ratio Constant')
box on

subplot(1,2,2)

hold on
for i = 1 : n
    plot(pi_le_d(:,2,i),pi_le_d(:,1,i),'Linewidth',2);
end
xlabel('\pi_1 = Mass Ratio (C_R Q_R / C_O Q_{SLR}) , l_E')
ylabel('\pi = Nondimensional Width (\chi/\surd b_B l_E)')
title('\pi_2 = Water Rise Rate Ratio Constant')
box on

% set(findobj('type','axes'),'fontsize',10)
% h_fig = gcf;
% set(h_fig,'PaperOrientation','portrait')
% set(h_fig,'PaperPosition', [0 0 8 3]) % [... ... max_width = 7.5 max_height = 9]
% tit = 'phi-3';
% print(tit,'-dtiff','-r400')
% movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\22 - Nondimensional Analysis')
% close all

%%
% chi_5
clf

subplot(1,2,1)

hold on
for i = 1 : n
    plot(pi_bfm_m(:,5,i),pi_bfm_m(:,1,i),'Linewidth',2);
end
xlabel('\pi_1 =  V_O/Q_R^{1/2} R^{-3/2}) , b_B')
ylabel('\pi =  Nondimensional Width (\chi/{\surd} Q_R/R)')
title('\pi_2 = Concentration Ratio Constant')
box on

subplot(1,2,2)

hold on
for i = 1 : n
    plot(pi_le_m(:,5,i),pi_le_m(:,1,i),'Linewidth',2);
end
xlabel('\pi_1 =  V_O/Q_R^{1/2} R^{-3/2} , l_E')
ylabel('\pi =  Nondimensional Width (\chi/{\surd} Q_R/R)')
title('\pi_2 = Concentration Ratio Constant')
box on

% set(findobj('type','axes'),'fontsize',10)
% h_fig = gcf;
% set(h_fig,'PaperOrientation','portrait')
% set(h_fig,'PaperPosition', [0 0 8 3]) % [... ... max_width = 7.5 max_height = 9]
% tit = 'phi-5';
% print(tit,'-dtiff','-r400')
% movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\22 - Nondimensional Analysis')
% close all

%%
% chi_7
clf

subplot(1,2,1)

hold on
for i = 1 : n
    plot(pi_bfm_m(:,4,i),pi_bfm_m(:,3,i),'Linewidth',2);
end
xlabel('\pi_1 = Mass Ratio (C_O Q_O / C_R Q_R) , b_B')
ylabel('\pi =  Nondimensional Width (\chi/b_B)')
title('\pi_2 = Mass Ratio (2) Constant')
box on

subplot(1,2,2)

hold on
for i = 1 : n
    plot(pi_le_m(:,4,i),pi_le_m(:,3,i),'Linewidth',2);
end
xlabel('\pi_1 = Mass Ratio (C_O Q_O / C_R Q_R) , l_E')
ylabel('\pi =  Nondimensional Width (\chi/b_B)')
title('\pi_2 = Mass Ratio (2) Constant')
box on

% set(findobj('type','axes'),'fontsize',10)
% h_fig = gcf;
% set(h_fig,'PaperOrientation','portrait')
% set(h_fig,'PaperPosition', [0 0 8 3]) % [... ... max_width = 7.5 max_height = 9]
% tit = 'phi-7';
% print(tit,'-dtiff','-r400')
% movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\22 - Nondimensional Analysis')
% close all

%%
% chi_8
clf

subplot(1,2,1)

hold on
for i = 1 : n
    plot(pi_bfm_m(:,2,i),pi_bfm_m(:,3,i),'Linewidth',2);
end
xlabel('\pi_1 = Discharge Ratio (Q_O / Q_R) , b_B')
ylabel('\pi =  Nondimensional Width (\chi/b_B)')
title('\pi_2 = Mass Ratio Constant')
box on

subplot(1,2,2)

hold on
for i = 1 : n
    plot(pi_le_m(:,2,i),pi_le_m(:,3,i),'Linewidth',2);
end
xlabel('\pi_1 = Discharge Ratio (Q_O / Q_R) , l_E')
ylabel('\pi =  Nondimensional Width (\chi/b_B)')
title('\pi_2 = Mass Ratio Constant')
box on

% set(findobj('type','axes'),'fontsize',10)
% h_fig = gcf;
% set(h_fig,'PaperOrientation','portrait')
% set(h_fig,'PaperPosition', [0 0 8 3]) % [... ... max_width = 7.5 max_height = 9]
% tit = 'phi-8';
% print(tit,'-dtiff','-r400')
% movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\22 - Nondimensional Analysis')
% close all

%%
% chi_9
clf

subplot(1,2,1)

hold on
for i = 1 : n
    plot(pi_r(:,2,i),pi_r(:,1,i),'Linewidth',2);
end
xlabel('\pi_1 = Mass Ratio (M_O + M_R / M_{SLR}) , R')
ylabel('\pi =  Nondimensional Width (\chi/b_B)')
title('\pi_2 = Discharge Ratio Constant')
box on

subplot(1,2,2)

hold on
for i = 1 : n
    plot(pi_rho(:,2,i),pi_rho(:,1,i),'Linewidth',2);
end
xlabel('\pi_1 = Mass Ratio (M_O + M_R / M_{SLR}) , \rho')
ylabel('\pi =  Nondimensional Width (\chi/b_B)')
title('\pi_2 = Discharge Ratio Constant')
box on

set(findobj('type','axes'),'fontsize',10)
h_fig = gcf;
set(h_fig,'PaperOrientation','portrait')
set(h_fig,'PaperPosition', [0 0 8 3]) % [... ... max_width = 7.5 max_height = 9]
tit = 'phi-9';
print(tit,'-dtiff','-r400')
movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\22 - Nondimensional Analysis')
close all

%%
% width and depths vs concentration
% Co_ = (10:10:80) *10^-3;

clf

subplot(4,3,1)

hold on
for i = 1 : k
    plot(Co_*1000,squeeze(pi_bfm_m(i,6,:)),'Linewidth',2)
end
xlabel('Concentration (mg/l)')
ylabel('TF Width (m)')
title('Basin Width')
box on

subplot(4,3,2)

hold on
for i = 1 : k
    plot(Co_*1000,squeeze(pi_cobfm(i,7,:)),'Linewidth',2)
end
plot([0,100],[.7,.7],'k--')
text(90,.72,'MSL')
xlabel('Concentration (mg/l)')
ylabel('Tidal Flat Depth (m)')
title('Basin Width')
box on

subplot(4,3,3)

hold on
for i = 1 : k
    plot(Co_*1000,squeeze(pi_cobfm(i,8,:)),'Linewidth',2)
end
plot([0,100],[.7,.7],'k--')
text(90,.72,'MSL')
xlabel('Concentration (mg/l)')
ylabel('Marsh Depth (m)')
title('Basin Width')
box on

subplot(4,3,4)

hold on
for i = 1 : k
    plot(Co_*1000,squeeze(pi_le_m(i,6,:)),'Linewidth',2)
end
xlabel('Concentration (mg/l)')
ylabel('TF Width (m)')
title('Estuary Length')
box on

subplot(4,3,5)

hold on
for i = 1 : k
    plot(Co_*1000,squeeze(pi_le_m(i,7,:)),'Linewidth',2)
end
plot([0,100],[.7,.7],'k--')
text(90,.72,'MSL')
xlabel('Concentration (mg/l)')
ylabel('Tidal Flat Depth (m)')
title('Estuary Length')
box on

subplot(4,3,6)

hold on
for i = 1 : k
    plot(Co_*1000,squeeze(pi_le_m(i,8,:)),'Linewidth',2)
end
plot([0,100],[.7,.7],'k--')
text(90,.72,'MSL')
xlabel('Concentration (mg/l)')
ylabel('Marsh Depth (m)')
title('Estuary Length')
box on

subplot(4,3,7)

hold on
for i = 1 : k
    plot(Co_*1000,squeeze(pi_coh(i,6,:)),'Linewidth',2)
end
xlabel('Concentration (mg/l)')
ylabel('TF Width (m)')
title('Tidal Range')
box on

subplot(4,3,8)

hold on
for i = 1 : k
    plot(Co_*1000,squeeze(pi_coh(i,7,:)),'Linewidth',2)
end
plot([0,100],[.7,.7],'k--')
text(90,.72,'MSL')
xlabel('Concentration (mg/l)')
ylabel('Tidal Flat Depth (m)')
title('Tidal Range')
box on

subplot(4,3,9)

hold on
for i = 1 : k
    plot(Co_*1000,squeeze(pi_coh(i,8,:)),'Linewidth',2)
end
plot([0,100],[.7,.7],'k--')
text(90,.72,'MSL')
xlabel('Concentration (mg/l)')
ylabel('Marsh Depth (m)')
title('Tidal Range')
box on

subplot(4,3,10)

hold on
for i = 1 : k
    plot(Co_*1000,squeeze(pi_coqf(i,6,:)),'Linewidth',2)
end
xlabel('Concentration (mg/l)')
ylabel('TF Width (m)')
title('River Discharge')
box on

subplot(4,3,11)

hold on
for i = 1 : k
    plot(Co_*1000,squeeze(pi_coqf(i,7,:)),'Linewidth',2)
end
plot([0,100],[.7,.7],'k--')
text(90,.72,'MSL')
xlabel('Concentration (mg/l)')
ylabel('Tidal Flat Depth (m)')
title('River Discharge')
box on

subplot(4,3,12)

hold on
for i = 1 : k
    plot(Co_*1000,squeeze(pi_coqf(i,8,:)),'Linewidth',2)
end
plot([0,100],[.7,.7],'k--')
text(90,.72,'MSL')
xlabel('Concentration (mg/l)')
ylabel('Marsh Depth (m)')
title('River Discharge')
box on

set(findobj('type','axes'),'fontsize',10)
h_fig = gcf;
set(h_fig,'PaperOrientation','portrait')
set(h_fig,'PaperPosition', [0 0 7 7]) % [... ... max_width = 7.5 max_height = 9]
tit = 'Co vs pars';
print(tit,'-dtiff','-r400')
movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\22 - Nondimensional Analysis')
close all

%%
% width and depths vs concentration
R_ = (1:8);

clf

subplot(4,3,1)

hold on
for i = 1 : k
    plot(R_,squeeze(pi_bfm_d(i,3,:)),'Linewidth',2)
end
xlabel('Sea Level Rise Rate (mm/yr)')
ylabel('TF Width (m)')
title('Basin Width')
box on

subplot(4,3,2)

hold on
for i = 1 : k
    plot(R_,squeeze(pi_bfm_d(i,4,:)),'Linewidth',2)
end
plot([0,10],[.7,.7],'k--')
text(8,.72,'MSL')
xlabel('Sea Level Rise Rate (mm/yr)')
ylabel('Tidal Flat Depth (m)')
title('Basin Width')
box on

subplot(4,3,3)

hold on
for i = 1 : k
    plot(R_,squeeze(pi_bfm_d(i,5,:)),'Linewidth',2)
end
plot([0,10],[.7,.7],'k--')
text(8,.72,'MSL')
xlabel('Sea Level Rise Rate (mm/yr)')
ylabel('Marsh Depth (m)')
title('Basin Width')
box on

subplot(4,3,4)

hold on
for i = 1 : k
    plot(R_,squeeze(pi_le_d(i,3,:)),'Linewidth',2)
end
xlabel('Sea Level Rise Rate (mm/yr)')
ylabel('TF Width (m)')
title('Estuary Length')
box on

subplot(4,3,5)

hold on
for i = 1 : k
    plot(R_,squeeze(pi_le_d(i,4,:)),'Linewidth',2)
end
plot([0,10],[.7,.7],'k--')
text(8,.72,'MSL')
xlabel('Sea Level Rise Rate (mm/yr)')
ylabel('Tidal Flat Depth (m)')
title('Estuary Length')
box on

subplot(4,3,6)

hold on
for i = 1 : k
    plot(R_,squeeze(pi_le_d(i,5,:)),'Linewidth',2)
end
plot([0,10],[.7,.7],'k--')
text(8,.72,'MSL')
xlabel('Sea Level Rise Rate (mm/yr)')
ylabel('Marsh Depth (m)')
title('Estuary Length')
box on

subplot(4,3,7)

hold on
for i = 1 : k
    plot(R_,squeeze(pi_rco(i,3,:)),'Linewidth',2)
end
xlabel('Concentration (mg/l)')
ylabel('TF Width (m)')
title('Ocean Concentration')
box on

subplot(4,3,8)

hold on
for i = 1 : k
    plot(R_,squeeze(pi_rco(i,4,:)),'Linewidth',2)
end
plot([0,10],[.7,.7],'k--')
text(8,.72,'MSL')
xlabel('Sea Level Rise Rate (mm/yr)')
ylabel('Tidal Flat Depth (m)')
title('Ocean Concentration')
box on

subplot(4,3,9)

hold on
for i = 1 : k
    plot(R_,squeeze(pi_rco(i,5,:)),'Linewidth',2)
end
plot([0,10],[.7,.7],'k--')
text(8,.72,'MSL')
xlabel('Sea Level Rise Rate (mm/yr)')
ylabel('Marsh Depth (m)')
title('Ocean Concentration')
box on

subplot(4,3,10)

hold on
for i = 1 : k
    plot(R_,squeeze(pi_rcf(i,3,:)),'Linewidth',2)
end
xlabel('Sea Level Rise Rate (mm/yr)')
ylabel('TF Width (m)')
title('River Concentration')
box on

subplot(4,3,11)

hold on
for i = 1 : k
    plot(R_,squeeze(pi_rcf(i,4,:)),'Linewidth',2)
end
plot([0,10],[.7,.7],'k--')
text(8,.72,'MSL')
xlabel('Sea Level Rise Rate (mm/yr)')
ylabel('Tidal Flat Depth (m)')
title('River Concentration')
box on

subplot(4,3,12)

hold on
for i = 1 : k
    plot(R_,squeeze(pi_rcf(i,5,:)),'Linewidth',2)
end
plot([0,10],[.7,.7],'k--')
text(8,.72,'MSL')
xlabel('Sea Level Rise Rate (mm/yr)')
ylabel('Marsh Depth (m)')
title('River Concentration')
box on

set(findobj('type','axes'),'fontsize',10)
h_fig = gcf;
set(h_fig,'PaperOrientation','portrait')
set(h_fig,'PaperPosition', [0 0 7 7]) % [... ... max_width = 7.5 max_height = 9]
tit = 'R vs pars';
print(tit,'-dtiff','-r400')
movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\22 - Nondimensional Analysis')
close all