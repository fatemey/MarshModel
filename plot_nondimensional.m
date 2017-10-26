%-------------- Plot nondimensional analysis results
% loading data
clear

load Sol_Co_bfm
C_f = 15 *10^-3;    % river concentration (kg/m3)
Q_f = 20 /2;         % river water discharge (m3/s)
L_E = 10 *10^3; % basin length (m)
R = 2 *10^-3/365/24/60/60;   % sea level rise (m/s)
T_T = 12 *60*60;   % tidal period (s) (= 12 hours)
H = 1.4/2;          % tidal amplitude (range/2) (m)
Co_ = (5:5:95) *10^-3;
bfm_ = (5:1:10) *10^3;
n = length(Co_);
k = length (bfm_); 

%%
% Co & bfm using average tidal volume
clf

subplot(2,4,2)

hold on
for i = 1 : n
    plot(L_E*bfm_*H/T_T/Q_f,squeeze(Sol(i,1,:))*sqrt(R)/sqrt(Q_f),'Linewidth',2);
end
xlabel('\pi_1 = Tidal Prism Ratio (Q_O/Q_R)')
ylabel('\pi = Tidal Flat Width Ratio (\chi/{\surd} Q_R/R)')
title('\pi_2 = Concentration Ratio (C_O/C_R)')
box on

subplot(2,4,1)

hold on
for i = 1 : n
    plot(round(R*bfm_*L_E/Q_f,4),squeeze(Sol(i,1,:))./bfm_','Linewidth',2);
end
xlabel('\pi_1 = Tidal Prism Ratio (Q_{SLR}/Q_R)')
ylabel('\pi = Tidal Flat Width Ratio (\chi/b_B)')
title('\pi_2 = Concentration Ratio (C_O/C_R)')
box on

subplot(2,4,3)

hold on
for i = 1 : n
    plot(L_E*bfm_*H/T_T/Q_f,squeeze(Sol(i,1,:))./bfm_','Linewidth',2);
end
xlabel('\pi_1 = Tidal Prism Ratio (Q_O/Q_R)')
ylabel('\pi = Tidal Flat Width Ratio (\chi/b_B)')
title('\pi_2 = Mass Ratio (C_O Q_O/C_R Q_R)')
box on

subplot(2,4,6)

hold on
for i = 1 : k
    plot(Co_/C_f,Sol(:,1,i)*sqrt(R)/sqrt(Q_f),'Linewidth',2);
end
xlabel('\pi_1 = Concentration Ratio (C_O/C_R)')
ylabel('\pi = Tidal Flat Width Ratio (\chi/{\surd} Q_R/R)')
title('\pi_2 = Tidal Prism Ratio (Q_O/Q_R)')
box on

subplot(2,4,5)

hold on
for i = 1 : k
    plot(Co_/C_f,Sol(:,1,i)/bfm_(i),'Linewidth',2)
end
title('\pi_2 = Tidal Prism Ratio (Q_{SLR}/Q_R)')
ylabel('\pi = Tidal Flat Width Ratio (\chi/b_B)')
xlabel('\pi_1 = Concentration Ratio (C_O/C_R)')
box on

subplot(2,4,7)

hold on
for i = 1 : k
    plot(Co_*L_E*bfm_(i)*H/T_T/C_f/Q_f,Sol(:,1,i)/bfm_(i),'Linewidth',2);
end
xlabel('\pi_1 = Mass Ratio (C_O Q_O/C_R Q_R)')
ylabel('\pi = Tidal Flat Width Ratio (\chi/b_B)')
title('\pi_2 = Tidal Prism Ratio (Q_O/Q_R)')
box on

subplot(2,4,4)

hold on
for i = 1 : n
    plot(C_f*Q_f./(L_E*bfm_*R*Co_(i)),squeeze(Sol(i,1,:))./sqrt(bfm_'*L_E),'Linewidth',2);
end
xlabel('\pi_1 = Mass Ratio (C_R Q_R / C_O Q_{SLR})')
ylabel('\pi = Tidal Flat Width Ratio (\chi/\surd b_B l_E)')
title('\pi_2 = Water Rise Rate Ratio (H/T/R)')
box on

subplot(2,4,8)

hold on
for i = 1 : k
    plot(C_f*Q_f./(L_E*bfm_(i)*R*Co_),Sol(:,1,i)/sqrt(bfm_(i)*L_E),'Linewidth',2);
end
% for i = 1 : k
%     plot(ones(n,1)*H/T_T/R,Sol(:,1,i)/sqrt(bfm_(i)'*L_E),'Linewidth',2);
% end
xlabel('\pi_1 = Mass Ratio (C_R Q_R / C_O Q_{SLR})')
ylabel('\pi = Tidal Flat Width Ratio (\chi/\surd b_B l_E)')
title('\pi_2 = Water Rise Rate Ratio (H/T/R)')
box on

% set(findobj('type','axes'),'fontsize',10)
% h_fig=gcf;
% set(h_fig,'PaperOrientation','portrait')
% set(h_fig,'PaperPosition', [0 0 14 8]) % [... ... max_width=7.5 max_height=9]
% tit='bfm_averageV';
% print(tit,'-dtiff','-r400')
% movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\22 - Nondimensional Analysis')
% close all

%%
% Co & bfm using accurate tidal volume
clf

subplot(2,4,2)

hold on
for i = 1 : n
    plot(L_E*((bfm_'-squeeze(Sol(i,1,:))).*squeeze(Sol(i,3,:))+squeeze(Sol(i,1,:)).*squeeze(Sol(i,2,:)))/T_T/Q_f,squeeze(Sol(i,1,:))*sqrt(R)/sqrt(Q_f),'Linewidth',2);
end
xlabel('\pi_1 = Tidal Prism Ratio (Q_O/Q_R)')
ylabel('\pi = Tidal Flat Width Ratio (\chi/{\surd} Q_R/R)')
title('\pi_2 = Concentration Ratio (C_O/C_R)')
box on

subplot(2,4,1)

hold on
for i = 1 : n
    plot(round(R*bfm_*L_E/Q_f,4),squeeze(Sol(i,1,:))./bfm_','Linewidth',2);
end
xlabel('\pi_1 = Tidal Prism Ratio (Q_{SLR}/Q_R)')
ylabel('\pi = Tidal Flat Width Ratio (\chi/b_B)')
title('\pi_2 = Concentration Ratio (C_O/C_R)')
box on

subplot(2,4,3)

hold on
for i = 1 : n
    plot(L_E*((bfm_'-squeeze(Sol(i,1,:))).*squeeze(Sol(i,3,:))+squeeze(Sol(i,1,:)).*squeeze(Sol(i,2,:)))/T_T/Q_f,squeeze(Sol(i,1,:))./bfm_','Linewidth',2);
end
xlabel('\pi_1 = Tidal Prism Ratio (Q_O/Q_R)')
ylabel('\pi = Tidal Flat Width Ratio (\chi/b_B)')
title('\pi_2 = Mass Ratio (C_O Q_O/C_R Q_R)')
box on

subplot(2,4,6)

hold on
for i = 1 : k
    plot(Co_/C_f,Sol(:,1,i)*sqrt(R)/sqrt(Q_f),'Linewidth',2);
end
xlabel('\pi_1 = Concentration Ratio (C_O/C_R)')
ylabel('\pi = Tidal Flat Width Ratio (\chi/{\surd} Q_R/R)')
title('\pi_2 = Tidal Prism Ratio (Q_O/Q_R)')
box on

subplot(2,4,5)

hold on
for i = 1 : k
    plot(Co_/C_f,Sol(:,1,i)/bfm_(i),'Linewidth',2)
end
title('\pi_2 = Tidal Prism Ratio (Q_{SLR}/Q_R)')
ylabel('\pi = Tidal Flat Width Ratio (\chi/b_B)')
xlabel('\pi_1 = Concentration Ratio (C_O/C_R)')
box on

subplot(2,4,7)

hold on
for i = 1 : k
    plot(Co_'.*(L_E*((bfm_(i)-Sol(:,1,i)).*Sol(:,3,i)+Sol(:,1,i).*Sol(:,2,i)))/T_T/C_f/Q_f,Sol(:,1,i)/bfm_(i),'Linewidth',2);
end
xlabel('\pi_1 = Mass Ratio (C_O Q_O/C_R Q_R)')
ylabel('\pi = Tidal Flat Width Ratio (\chi/b_B)')
title('\pi_2 = Tidal Prism Ratio (Q_O/Q_R)')
box on

subplot(2,4,4)

hold on
for i = 1 : n
    plot(C_f*Q_f./(L_E*bfm_*R*Co_(i)),squeeze(Sol(i,1,:))./sqrt(bfm_'*L_E),'Linewidth',2);
end
xlabel('\pi_1 = Mass Ratio (C_R Q_R / C_O Q_{SLR})')
ylabel('\pi = Tidal Flat Width Ratio (\chi/\surd b_B l_E)')
title('\pi_2 = Water Rise Rate Ratio (H/T/R)')
box on

subplot(2,4,8)

hold on
for i = 1 : k
    plot(ones(n,1)*H/T_T/R,Sol(:,1,i)/sqrt(bfm_(i)*L_E),'Linewidth',2);
end
xlabel('\pi_1 = Water Rise Rate Ratio (H/T/R)')
ylabel('\pi = Tidal Flat Width Ratio (\chi/\surd b_B l_E)')
title('\pi_2 = Mass Ratio (C_R Q_R / C_O Q_{SLR})')
box on

% set(findobj('type','axes'),'fontsize',10)
% h_fig=gcf;
% set(h_fig,'PaperOrientation','portrait')
% set(h_fig,'PaperPosition', [0 0 14 8]) % [... ... max_width=7.5 max_height=9]
% tit='bfm_accurateV';
% print(tit,'-dtiff','-r400')
% movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\22 - Nondimensional Analysis')
% close all

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

% set(findobj('type','axes'),'fontsize',10)
% h_fig=gcf;
% set(h_fig,'PaperOrientation','portrait')
% set(h_fig,'PaperPosition', [0 0 7 7]) % [... ... max_width=7.5 max_height=9]
% tit='bfm';
% print(tit,'-dtiff','-r400')
% movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\22 - Nondimensional Analysis')
% close all
