%-------------- Plot nondimensional analysis results
% loading data
clear

load Sol_Co_R
C_f = 15 *10^-3;    % river concentration (kg/m3)
Q_f = 20 /2;         % river water discharge (m3/s)
L_E = 10 *10^3; % basin length (m)
b_fm = 10  *10^3; % total basin width (both sides of the channel) (m)
T_T = 12 *60*60;   % tidal period (s) (= 12 hours)
H = 1.4/2;          % tidal amplitude (range/2) (m)
Co_ = (5:5:95) *10^-3;
R_ = (1:2:9) *10^-3/365/24/60/60;
n = length(Co_);
k = length (R_); 

%%
% Co & R using average tidal volume
clf

subplot(2,3,2)

hold on
for i = 1 : n
    plot(ones(k,1)*L_E*b_fm*H/T_T/Q_f,squeeze(Sol(i,1,:)).*sqrt(R_')/sqrt(Q_f),'Linewidth',2);
end
xlabel('\pi_1 = Tidal Prism Ratio (Q_O/Q_R)')
ylabel('\pi = Tidal Flat Width Ratio (\chi/{\surd} (Q_R/R))')
title('\pi_2 = Concentration Ratio (C_O/C_R)')
box on

subplot(2,3,1)

hold on
for i = 1 : n
    plot(round(R_*b_fm*L_E/Q_f,4),squeeze(Sol(i,1,:))/b_fm,'Linewidth',2);
end
xlabel('\pi_1 = Tidal Prism Ratio (Q_{SLR}/Q_R)')
ylabel('\pi = Tidal Flat Width Ratio (\chi/b_B)')
title('\pi_2 = Concentration Ratio (C_O/C_R)')
box on

subplot(2,3,3)

hold on
for i = 1 : n
    plot(ones(k,1)*L_E*b_fm*H/T_T/Q_f,squeeze(Sol(i,1,:))/b_fm,'Linewidth',2);
end
xlabel('\pi_1 = Tidal Prism Ratio (Q_O/Q_R)')
ylabel('\pi = Tidal Flat Width Ratio (\chi/b_B)')
title('\pi_2 = Mass Ratio (C_O V_T/C_R Q_R T)')
box on

subplot(2,3,5)

hold on
for i = 1 : k
    plot(Co_/C_f,Sol(:,1,i)*sqrt(R_(i))/sqrt(Q_f),'Linewidth',2);
end
xlabel('\pi_1 = Concentration Ratio (C_O/C_R)')
ylabel('\pi = Tidal Flat Width Ratio (\chi/{\surd} (Q_R/R))')
title('\pi_2 = Tidal Prism Ratio (Q_O/Q_R)')
box on

subplot(2,3,4)

hold on
for i = 1 : k
    plot(Co_/C_f,Sol(:,1,i)/b_fm,'Linewidth',2)
end
title('\pi_2 = Tidal Prism Ratio (Q_{SLR}/Q_R)')
ylabel('\pi = Tidal Flat Width Ratio (\chi/b_B)')
xlabel('\pi_1 = Concentration Ratio (C_O/C_R)')
box on

subplot(2,3,6)

hold on
for i = 1 : k
    plot(Co_*L_E*b_fm*H/T_T/C_f/Q_f,Sol(:,1,i)/b_fm,'Linewidth',2);
end
xlabel('\pi_1 = Mass Ratio (C_O V_T/C_R Q_R T)')
ylabel('\pi = Tidal Flat Width Ratio (\chi/b_B)')
title('\pi_2 = Tidal Prism Ratio (Q_O/Q_R)')
box on

set(findobj('type','axes'),'fontsize',10)
h_fig=gcf;
set(h_fig,'PaperOrientation','portrait')
set(h_fig,'PaperPosition', [0 0 9 7]) % [... ... max_width=7.5 max_height=9]
tit='R_averageV';
print(tit,'-dtiff','-r400')
movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\22 - Nondimensional Analysis')
close all

%%
% Co & bfm using accurate tidal volume
clf

subplot(2,3,2)

hold on
for i = 1 : n
    plot(L_E*((b_fm-squeeze(Sol(i,1,:))).*squeeze(Sol(i,3,:))+squeeze(Sol(i,1,:)).*squeeze(Sol(i,2,:)))/T_T/Q_f,squeeze(Sol(i,1,:)).*sqrt(R_')/sqrt(Q_f),'Linewidth',2);
end
xlabel('\pi_1 = Tidal Prism Ratio (Q_O/Q_R)')
ylabel('\pi = Tidal Flat Width Ratio (\chi/{\surd} (Q_R/R))')
title('\pi_2 = Concentration Ratio (C_O/C_R)')
box on

subplot(2,3,1)

hold on
for i = 1 : n
    plot(round(R_*b_fm*L_E/Q_f,4),squeeze(Sol(i,1,:))/b_fm,'Linewidth',2);
end
xlabel('\pi_1 = Tidal Prism Ratio (Q_{SLR}/Q_R)')
ylabel('\pi = Tidal Flat Width Ratio (\chi/b_B)')
title('\pi_2 = Concentration Ratio (C_O/C_R)')
box on

subplot(2,3,3)

hold on
for i = 1 : n
    plot(L_E*((b_fm-squeeze(Sol(i,1,:))).*squeeze(Sol(i,3,:))+squeeze(Sol(i,1,:)).*squeeze(Sol(i,2,:)))/T_T/Q_f,squeeze(Sol(i,1,:))/b_fm,'Linewidth',2);
end
xlabel('\pi_1 = Tidal Prism Ratio (Q_O/Q_R)')
ylabel('\pi = Tidal Flat Width Ratio (\chi/b_B)')
title('\pi_2 = Mass Ratio (C_O V_T/C_R Q_R T)')
box on

subplot(2,3,5)

hold on
for i = 1 : k
    plot(Co_/C_f,Sol(:,1,i)*sqrt(R_(i))/sqrt(Q_f),'Linewidth',2);
end
xlabel('\pi_1 = Concentration Ratio (C_O/C_R)')
ylabel('\pi = Tidal Flat Width Ratio (\chi/{\surd} (Q_R/R))')
title('\pi_2 = Tidal Prism Ratio (Q_O/Q_R)')
box on

subplot(2,3,4)

hold on
for i = 1 : k
    plot(Co_/C_f,Sol(:,1,i)/b_fm,'Linewidth',2)
end
title('\pi_2 = Tidal Prism Ratio (Q_{SLR}/Q_R)')
ylabel('\pi = Tidal Flat Width Ratio (\chi/b_B)')
xlabel('\pi_1 = Concentration Ratio (C_O/C_R)')
box on

subplot(2,3,6)

hold on
for i = 1 : k
    plot(Co_'.*(L_E*((b_fm-Sol(:,1,i)).*Sol(:,3,i)+Sol(:,1,i).*Sol(:,2,i)))/T_T/C_f/Q_f,Sol(:,1,i)/b_fm,'Linewidth',2);
end
xlabel('\pi_1 = Mass Ratio (C_O V_T/C_R Q_R T)')
ylabel('\pi = Tidal Flat Width Ratio (\chi/b_B)')
title('\pi_2 = Tidal Prism Ratio (Q_O/Q_R)')
box on

set(findobj('type','axes'),'fontsize',10)
h_fig=gcf;
set(h_fig,'PaperOrientation','portrait')
set(h_fig,'PaperPosition', [0 0 9 7]) % [... ... max_width=7.5 max_height=9]
tit='R_accurateV';
print(tit,'-dtiff','-r400')
movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\22 - Nondimensional Analysis')
close all

%%
% width and depths vs concentration
clf

subplot(2,2,1)

hold on
for i = 1 : k
    plot(Co_*1000,Sol(:,1,i),'Linewidth',2)
end
xlabel('Concentration (mg/l)')
ylabel('TF Width (m)')
box on

subplot(2,2,2)

hold on
for i = 1 : k
    plot(Co_*1000,Sol(:,2,i),'Linewidth',2)
end
plot([0,100],[.7,.7],'k--')
text(90,.72,'MSL')
xlabel('Concentration (mg/l)')
ylabel('Tidal Flat Depth (m)')
box on

subplot(2,2,3)

hold on
for i = 1 : k
    plot(Co_*1000,Sol(:,3,i),'Linewidth',2)
end
plot([0,100],[.7,.7],'k--')
text(90,.72,'MSL')
xlabel('Concentration (mg/l)')
ylabel('Marsh Depth (m)')
box on

set(findobj('type','axes'),'fontsize',10)
h_fig=gcf;
set(h_fig,'PaperOrientation','portrait')
set(h_fig,'PaperPosition', [0 0 7 7]) % [... ... max_width=7.5 max_height=9]
tit='R';
print(tit,'-dtiff','-r400')
movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\22 - Nondimensional Analysis')
close all

