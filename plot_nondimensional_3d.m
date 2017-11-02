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
Co_ = (10:10:80) *10^-3;
n1 = length(Co_);
k = length (par_v);

load Sol_bfm_co
bfm_ = b_fm * par_v;

for i = 1 : n1
    dat(:,1,i) = squeeze(Sol(i,1,:))/sqrt(Q_f/R);
    dat(:,2,i) = L_E*bfm_*H/T_T/Q_f;
    dat(:,3,i) = squeeze(Sol(i,1,:))./bfm_';
    dat(:,4,i) = Co_(i)*L_E*bfm_*H/T_T/Q_f/C_f;
    dat(:,5,i) = L_E*bfm_*H/Q_f^1.5*R^-1.5;
    dat(:,6,i) = (Co_(i)*(L_E*bfm_*H/T_T-Q_f)+C_f*Q_f)/R./bfm_/L_E/rho_s;
    dat(:,7,i) = ones(size(par_v))*Co_(i)/C_f;
    dat(:,8,i) = ones(size(par_v))*Co_(i)*H/T_T/R/rho_s;
%     dat(:,7,i) = squeeze(Sol(i,1,:));
%     dat(:,8,i) = squeeze(Sol(i,2,:));
%     dat(:,9,i) = squeeze(Sol(i,3,:));
end
datt=[];
for i = 1 : n1
    datt=vertcat(datt,dat(:,:,i));
end

pi_bfm_co = datt;
clear dat datt

load Sol_le_co
LE_ = L_E * par_v;

for i = 1 : n1
    dat(:,1,i) = squeeze(Sol(i,1,:))/sqrt(Q_f/R);
    dat(:,2,i)  = LE_*b_fm*H/T_T/Q_f;
    dat(:,3,i) = squeeze(Sol(i,1,:))/b_fm;
    dat(:,4,i)  = Co_(i)*LE_*b_fm*H/T_T/Q_f/C_f;
    dat(:,5,i) = LE_*b_fm*H/Q_f^1.5*R^-1.5;
    dat(:,6,i) = (Co_(i)*(LE_*b_fm*H/T_T-Q_f)+C_f*Q_f)/R/b_fm./LE_/rho_s;
    dat(:,7,i) = ones(size(par_v))*Co_(i)/C_f;
    dat(:,8,i) = ones(size(par_v))*Co_(i)*H/T_T/R/rho_s;
    %     dat(:,7,i) = squeeze(Sol(i,1,:));
%     dat(:,8,i) = squeeze(Sol(i,2,:));
%     dat(:,9,i) = squeeze(Sol(i,3,:));
end
datt=[];
for i = 1 : n1
    datt=vertcat(datt,dat(:,:,i));
end

pi_le_co = datt;
clear dat datt

R_ = (1:8) *10^-3/365/24/60/60;
n2 = length(R_);

load Sol_bfm_r
bfm_ = b_fm ./ par_v;

for i = 1 : n2
    dat(:,1,i) = squeeze(Sol(i,1,:))./sqrt(bfm_'*L_E);
    dat(:,2,i) = C_f*Q_f/C_o/R_(i)/L_E./bfm_;
    dat(:,3,i) = ones(size(par_v))*H/T_T/R_(i);
%     dat(:,3,i) = squeeze(Sol(i,1,:));
%     dat(:,4,i) = squeeze(Sol(i,2,:));
%     dat(:,5,i) = squeeze(Sol(i,3,:));
end
datt=[];
for i = 1 : n1
    datt=vertcat(datt,dat(:,:,i));
end

pi_bfm_r = datt;
clear dat datt

load Sol_le_r
LE_ = L_E ./ par_v;

for i = 1 : n2
    dat(:,1,i) = squeeze(Sol(i,1,:))./sqrt(b_fm*LE_');
    dat(:,2,i) = C_f*Q_f/C_o/R_(i)./LE_/b_fm;
    dat(:,3,i) = ones(size(par_v))*H/T_T/R_(i);
%     dat(:,3,i) = squeeze(Sol(i,1,:));
%     dat(:,4,i) = squeeze(Sol(i,2,:));
%     dat(:,5,i) = squeeze(Sol(i,3,:));
end
datt=[];
for i = 1 : n1
    datt=vertcat(datt,dat(:,:,i));
end

pi_le_r = datt;
clear dat datt

%%
% chi_2
figure(1)

subplot(1,2,1)

xx = pi_bfm_co(:,7);
yy = pi_bfm_co(:,2);
zz = pi_bfm_co(:,1);
xqq=linspace(min(xx),max(xx));
yqq=linspace(min(yy),max(yy));
[xq,yq] = meshgrid(xqq,yqq);
vq = griddata(xx,yy,zz,xq,yq);
mesh(xq,yq,vq)

xlabel('\pi_1 = Tidal Prism Ratio')
zlabel('\pi =  Nondimensional Width')
ylabel('\pi_2 = Concentration Ratio')
% h = rotate3d; 
% set(h,'ActionPreCallback','set(gcf,''windowbuttonmotionfcn'',@align_axislabel)'); 
% set(h,'ActionPostCallback','set(gcf,''windowbuttonmotionfcn'','''')'); 
% align_axislabel([], gca) % optional, this will call align_axislabels once before manually rotating the plot 
box on

subplot(1,2,2)

xx = pi_le_co(:,7);
yy = pi_le_co(:,2);
zz = pi_le_co(:,1);
xqq=linspace(min(xx),max(xx));
yqq=linspace(min(yy),max(yy));
[xq,yq] = meshgrid(xqq,yqq);
vq = griddata(xx,yy,zz,xq,yq);
mesh(xq,yq,vq)

xlabel('\pi_1 = Tidal Prism Ratio')
zlabel('\pi =  Nondimensional Width')
ylabel('\pi_2 = Concentration Ratio')
% h = rotate3d;  
% set(h,'ActionPreCallback','set(gcf,''windowbuttonmotionfcn'',@align_axislabel)');  
% set(h,'ActionPostCallback','set(gcf,''windowbuttonmotionfcn'','''')'); 
% align_axislabel([], gca) % optional, this will call align_axislabels once before manually rotating the plot 
box on

% set(findobj('type','axes'),'fontsize',7)
% h_fig = gcf;
% set(h_fig,'PaperOrientation','portrait')
% set(h_fig,'PaperPosition', [0 0 8 4]) % [... ... max_width = 7.5 max_height = 9]
% tit = 'phi-3d-2';
% print(tit,'-dtiff','-r400')
% movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\22 - Nondimensional Analysis')
%close all 

%%
% chi_3
figure(2)

subplot(1,2,1)

xx = pi_bfm_r(:,3);
yy = pi_bfm_r(:,2);
zz = pi_bfm_r(:,1);
xqq=linspace(min(xx),max(xx));
yqq=linspace(min(yy),max(yy));
[xq,yq] = meshgrid(xqq,yqq);
vq = griddata(xx,yy,zz,xq,yq);
mesh(xq,yq,vq)

xlabel('\pi_1 = Mass Ratio')
zlabel('\pi = Nondimensional Width')
ylabel('\pi_2 = Water Rise Rate Ratio')
% h = rotate3d;  
% set(h,'ActionPreCallback','set(gcf,''windowbuttonmotionfcn'',@align_axislabel)');  
% set(h,'ActionPostCallback','set(gcf,''windowbuttonmotionfcn'','''')'); 
% align_axislabel([], gca) % optional, this will call align_axislabels once before manually rotating the plot 
box on

subplot(1,2,2)

xx = pi_le_r(:,3);
yy = pi_le_r(:,2);
zz = pi_le_r(:,1);
xqq=linspace(min(xx),max(xx));
yqq=linspace(min(yy),max(yy));
[xq,yq] = meshgrid(xqq,yqq);
vq = griddata(xx,yy,zz,xq,yq);
mesh(xq,yq,vq)

xlabel('\pi_1 = Mass Ratio')
zlabel('\pi = Nondimensional Width')
ylabel('\pi_2 = Water Rise Rate Ratio')
% h = rotate3d;  
% set(h,'ActionPreCallback','set(gcf,''windowbuttonmotionfcn'',@align_axislabel)');  
% set(h,'ActionPostCallback','set(gcf,''windowbuttonmotionfcn'','''')'); 
% align_axislabel([], gca) % optional, this will call align_axislabels once before manually rotating the plot 
box on

% set(findobj('type','axes'),'fontsize',7)
% h_fig = gcf;
% set(h_fig,'PaperOrientation','portrait')
% set(h_fig,'PaperPosition', [0 0 8 4]) % [... ... max_width = 7.5 max_height = 9]
% tit = 'phi-3d-3';
% print(tit,'-dtiff','-r400')
% movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\22 - Nondimensional Analysis')
%close all 

%%
% chi_5
figure(3)

subplot(1,2,1)

xx = pi_bfm_co(:,7);
yy = pi_bfm_co(:,5);
zz = pi_bfm_co(:,1);
xqq=linspace(min(xx),max(xx));
yqq=linspace(min(yy),max(yy));
[xq,yq] = meshgrid(xqq,yqq);
vq = griddata(xx,yy,zz,xq,yq);
mesh(xq,yq,vq)

xlabel('\pi_1 =  Volume Ratio?')
zlabel('\pi =  Nondimensional Width')
ylabel('\pi_2 = Concentration Ratio')
% h = rotate3d;  
% set(h,'ActionPreCallback','set(gcf,''windowbuttonmotionfcn'',@align_axislabel)');  
% set(h,'ActionPostCallback','set(gcf,''windowbuttonmotionfcn'','''')'); 
% align_axislabel([], gca) % optional, this will call align_axislabels once before manually rotating the plot 
box on

subplot(1,2,2)

xx = pi_le_co(:,7);
yy = pi_le_co(:,5);
zz = pi_le_co(:,1);
xqq=linspace(min(xx),max(xx));
yqq=linspace(min(yy),max(yy));
[xq,yq] = meshgrid(xqq,yqq);
vq = griddata(xx,yy,zz,xq,yq);
mesh(xq,yq,vq)

xlabel('\pi_1 =  Volume Ratio?')
zlabel('\pi =  Nondimensional Width')
ylabel('\pi_2 = Concentration Ratio')
% h = rotate3d;  
% set(h,'ActionPreCallback','set(gcf,''windowbuttonmotionfcn'',@align_axislabel)');  
% set(h,'ActionPostCallback','set(gcf,''windowbuttonmotionfcn'','''')'); 
% align_axislabel([], gca) % optional, this will call align_axislabels once before manually rotating the plot 
box on

% set(findobj('type','axes'),'fontsize',7)
% h_fig = gcf;
% set(h_fig,'PaperOrientation','portrait')
% set(h_fig,'PaperPosition', [0 0 8 4]) % [... ... max_width = 7.5 max_height = 9]
% tit = 'phi-3d-5';
% print(tit,'-dtiff','-r400')
% movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\22 - Nondimensional Analysis')
%close all 

%%
% chi_7
figure(4)

subplot(1,2,1)

xx = pi_bfm_co(:,8);
yy = pi_bfm_co(:,4);
zz = pi_bfm_co(:,3);
xqq=linspace(min(xx),max(xx));
yqq=linspace(min(yy),max(yy));
[xq,yq] = meshgrid(xqq,yqq);
vq = griddata(xx,yy,zz,xq,yq);
mesh(xq,yq,vq)

xlabel('\pi_1 = Mass Ratio')
zlabel('\pi =  Nondimensional Width')
ylabel('\pi_2 = Mass Ratio')
% h = rotate3d;  
% set(h,'ActionPreCallback','set(gcf,''windowbuttonmotionfcn'',@align_axislabel)');  
% set(h,'ActionPostCallback','set(gcf,''windowbuttonmotionfcn'','''')'); 
% align_axislabel([], gca) % optional, this will call align_axislabels once before manually rotating the plot 
box on

subplot(1,2,2)

xx = pi_le_co(:,8);
yy = pi_le_co(:,4);
zz = pi_le_co(:,3);
xqq=linspace(min(xx),max(xx));
yqq=linspace(min(yy),max(yy));
[xq,yq] = meshgrid(xqq,yqq);
vq = griddata(xx,yy,zz,xq,yq);
mesh(xq,yq,vq)

xlabel('\pi_1 = Mass Ratio')
zlabel('\pi =  Nondimensional Width')
ylabel('\pi_2 = Mass Ratio')
% h = rotate3d;  
% set(h,'ActionPreCallback','set(gcf,''windowbuttonmotionfcn'',@align_axislabel)');  
% set(h,'ActionPostCallback','set(gcf,''windowbuttonmotionfcn'','''')'); 
% align_axislabel([], gca) % optional, this will call align_axislabels once before manually rotating the plot 
box on

% set(findobj('type','axes'),'fontsize',7)
% h_fig = gcf;
% set(h_fig,'PaperOrientation','portrait')
% set(h_fig,'PaperPosition', [0 0 8 4]) % [... ... max_width = 7.5 max_height = 9]
% tit = 'phi-3d-7';
% print(tit,'-dtiff','-r400')
% movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\22 - Nondimensional Analysis')
%close all 

%%
% chi_8
figure(5)

subplot(1,2,1)

xx = pi_bfm_co(:,8);
yy = pi_bfm_co(:,2);
zz = pi_bfm_co(:,3);
xqq=linspace(min(xx),max(xx));
yqq=linspace(min(yy),max(yy));
[xq,yq] = meshgrid(xqq,yqq);
vq = griddata(xx,yy,zz,xq,yq);
mesh(xq,yq,vq)

xlabel('\pi_1 = Discharge Ratio')
zlabel('\pi =  Nondimensional Width')
ylabel('\pi_2 = Mass Ratio')
% h = rotate3d;  
% set(h,'ActionPreCallback','set(gcf,''windowbuttonmotionfcn'',@align_axislabel)');  
% set(h,'ActionPostCallback','set(gcf,''windowbuttonmotionfcn'','''')'); 
% align_axislabel([], gca) % optional, this will call align_axislabels once before manually rotating the plot 
box on

subplot(1,2,2)

xx = pi_le_co(:,8);
yy = pi_le_co(:,2);
zz = pi_le_co(:,3);
xqq=linspace(min(xx),max(xx));
yqq=linspace(min(yy),max(yy));
[xq,yq] = meshgrid(xqq,yqq);
vq = griddata(xx,yy,zz,xq,yq);
mesh(xq,yq,vq)

xlabel('\pi_1 = Discharge Ratio')
zlabel('\pi =  Nondimensional Width')
ylabel('\pi_2 = Mass Ratio')
% h = rotate3d;  
% set(h,'ActionPreCallback','set(gcf,''windowbuttonmotionfcn'',@align_axislabel)');  
% set(h,'ActionPostCallback','set(gcf,''windowbuttonmotionfcn'','''')'); 
% align_axislabel([], gca) % optional, this will call align_axislabels once before manually rotating the plot 
box on

% set(findobj('type','axes'),'fontsize',7)
% h_fig = gcf;
% set(h_fig,'PaperOrientation','portrait')
% set(h_fig,'PaperPosition', [0 0 8 4]) % [... ... max_width = 7.5 max_height = 9]
% tit = 'phi-3d-8';
% print(tit,'-dtiff','-r400')
% movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\22 - Nondimensional Analysis')
%close all 

%%
% chi_9
figure(6)

subplot(1,2,1)

xx = pi_bfm_co(:,2);
yy = pi_bfm_co(:,6);
zz = pi_bfm_co(:,3);
xqq=linspace(min(xx),max(xx));
yqq=linspace(min(yy),max(yy));
[xq,yq] = meshgrid(xqq,yqq);
vq = griddata(xx,yy,zz,xq,yq);
mesh(xq,yq,vq)

xlabel('\pi_1 = Mass Ratio')
zlabel('\pi =  Nondimensional Width')
ylabel('\pi_2 = Discharge Ratio')
% h = rotate3d;  
% set(h,'ActionPreCallback','set(gcf,''windowbuttonmotionfcn'',@align_axislabel)');  
% set(h,'ActionPostCallback','set(gcf,''windowbuttonmotionfcn'','''')'); 
% align_axislabel([], gca) % optional, this will call align_axislabels once before manually rotating the plot 
box on

subplot(1,2,2)

xx = pi_le_co(:,2);
yy = pi_le_co(:,6);
zz = pi_le_co(:,3);
xqq=linspace(min(xx),max(xx));
yqq=linspace(min(yy),max(yy));
[xq,yq] = meshgrid(xqq,yqq);
vq = griddata(xx,yy,zz,xq,yq);
mesh(xq,yq,vq)

xlabel('\pi_1 = Mass Ratio')
zlabel('\pi =  Nondimensional Width')
ylabel('\pi_2 = Discharge Ratio')
% h = rotate3d;  
% set(h,'ActionPreCallback','set(gcf,''windowbuttonmotionfcn'',@align_axislabel)');  
% set(h,'ActionPostCallback','set(gcf,''windowbuttonmotionfcn'','''')'); 
% align_axislabel([], gca) % optional, this will call align_axislabels once before manually rotating the plot 
box on

% set(findobj('type','axes'),'fontsize',7)
% h_fig = gcf;
% set(h_fig,'PaperOrientation','portrait')
% set(h_fig,'PaperPosition', [0 0 8 4]) % [... ... max_width = 7.5 max_height = 9]
% tit = 'phi-3d-9';
% print(tit,'-dtiff','-r400')
% movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\22 - Nondimensional Analysis')
%close all 

%%
% width and depths vs concentration
% Co_ = (10:10:80) *10^-3;

% figure
% 
% subplot(4,3,1)
% 
% hold on
% for i = 1 : k
%     plot(Co_*1000,squeeze(pi_bfm_co(i,6,:)),'Linewidth',2)
% end
% xlabel('Concentration (mg/l)')
% zlabel('TF Width (m)')
% ylabel('Basin Width')
% h = rotate3d;  set(h,'ActionPreCallback','set(gcf,''windowbuttonmotionfcn'',@align_axislabel)');  set(h,'ActionPostCallback','set(gcf,''windowbuttonmotionfcn'','''')'); box on
% 
% subplot(4,3,2)
% 
% hold on
% for i = 1 : k
%     plot(Co_*1000,squeeze(pi_cobfm(i,7,:)),'Linewidth',2)
% end
% plot([0,100],[.7,.7],'k--')
% text(90,.72,'MSL')
% xlabel('Concentration (mg/l)')
% zlabel('Tidal Flat Depth (m)')
% ylabel('Basin Width')
% h = rotate3d;  set(h,'ActionPreCallback','set(gcf,''windowbuttonmotionfcn'',@align_axislabel)');  set(h,'ActionPostCallback','set(gcf,''windowbuttonmotionfcn'','''')'); box on
% 
% subplot(4,3,3)
% 
% hold on
% for i = 1 : k
%     plot(Co_*1000,squeeze(pi_cobfm(i,8,:)),'Linewidth',2)
% end
% plot([0,100],[.7,.7],'k--')
% text(90,.72,'MSL')
% xlabel('Concentration (mg/l)')
% zlabel('Marsh Depth (m)')
% ylabel('Basin Width')
% h = rotate3d;  set(h,'ActionPreCallback','set(gcf,''windowbuttonmotionfcn'',@align_axislabel)');  set(h,'ActionPostCallback','set(gcf,''windowbuttonmotionfcn'','''')'); box on
% 
% subplot(4,3,4)
% 
% hold on
% for i = 1 : k
%     plot(Co_*1000,squeeze(pi_le_co(i,6,:)),'Linewidth',2)
% end
% xlabel('Concentration (mg/l)')
% zlabel('TF Width (m)')
% ylabel('Estuary Length')
% h = rotate3d;  set(h,'ActionPreCallback','set(gcf,''windowbuttonmotionfcn'',@align_axislabel)');  set(h,'ActionPostCallback','set(gcf,''windowbuttonmotionfcn'','''')'); box on
% 
% subplot(4,3,5)
% 
% hold on
% for i = 1 : k
%     plot(Co_*1000,squeeze(pi_le_co(i,7,:)),'Linewidth',2)
% end
% plot([0,100],[.7,.7],'k--')
% text(90,.72,'MSL')
% xlabel('Concentration (mg/l)')
% zlabel('Tidal Flat Depth (m)')
% ylabel('Estuary Length')
% h = rotate3d;  set(h,'ActionPreCallback','set(gcf,''windowbuttonmotionfcn'',@align_axislabel)');  set(h,'ActionPostCallback','set(gcf,''windowbuttonmotionfcn'','''')'); box on
% 
% subplot(4,3,6)
% 
% hold on
% for i = 1 : k
%     plot(Co_*1000,squeeze(pi_le_co(i,8,:)),'Linewidth',2)
% end
% plot([0,100],[.7,.7],'k--')
% text(90,.72,'MSL')
% xlabel('Concentration (mg/l)')
% zlabel('Marsh Depth (m)')
% ylabel('Estuary Length')
% h = rotate3d;  set(h,'ActionPreCallback','set(gcf,''windowbuttonmotionfcn'',@align_axislabel)');  set(h,'ActionPostCallback','set(gcf,''windowbuttonmotionfcn'','''')'); box on
% 
% subplot(4,3,7)
% 
% hold on
% for i = 1 : k
%     plot(Co_*1000,squeeze(pi_coh(i,6,:)),'Linewidth',2)
% end
% xlabel('Concentration (mg/l)')
% zlabel('TF Width (m)')
% ylabel('Tidal Range')
% h = rotate3d;  set(h,'ActionPreCallback','set(gcf,''windowbuttonmotionfcn'',@align_axislabel)');  set(h,'ActionPostCallback','set(gcf,''windowbuttonmotionfcn'','''')'); box on
% 
% subplot(4,3,8)
% 
% hold on
% for i = 1 : k
%     plot(Co_*1000,squeeze(pi_coh(i,7,:)),'Linewidth',2)
% end
% plot([0,100],[.7,.7],'k--')
% text(90,.72,'MSL')
% xlabel('Concentration (mg/l)')
% zlabel('Tidal Flat Depth (m)')
% ylabel('Tidal Range')
% h = rotate3d;  set(h,'ActionPreCallback','set(gcf,''windowbuttonmotionfcn'',@align_axislabel)');  set(h,'ActionPostCallback','set(gcf,''windowbuttonmotionfcn'','''')'); box on
% 
% subplot(4,3,9)
% 
% hold on
% for i = 1 : k
%     plot(Co_*1000,squeeze(pi_coh(i,8,:)),'Linewidth',2)
% end
% plot([0,100],[.7,.7],'k--')
% text(90,.72,'MSL')
% xlabel('Concentration (mg/l)')
% zlabel('Marsh Depth (m)')
% ylabel('Tidal Range')
% h = rotate3d;  set(h,'ActionPreCallback','set(gcf,''windowbuttonmotionfcn'',@align_axislabel)');  set(h,'ActionPostCallback','set(gcf,''windowbuttonmotionfcn'','''')'); box on
% 
% subplot(4,3,10)
% 
% hold on
% for i = 1 : k
%     plot(Co_*1000,squeeze(pi_coqf(i,6,:)),'Linewidth',2)
% end
% xlabel('Concentration (mg/l)')
% zlabel('TF Width (m)')
% ylabel('River Discharge')
% h = rotate3d;  set(h,'ActionPreCallback','set(gcf,''windowbuttonmotionfcn'',@align_axislabel)');  set(h,'ActionPostCallback','set(gcf,''windowbuttonmotionfcn'','''')'); box on
% 
% subplot(4,3,11)
% 
% hold on
% for i = 1 : k
%     plot(Co_*1000,squeeze(pi_coqf(i,7,:)),'Linewidth',2)
% end
% plot([0,100],[.7,.7],'k--')
% text(90,.72,'MSL')
% xlabel('Concentration (mg/l)')
% zlabel('Tidal Flat Depth (m)')
% ylabel('River Discharge')
% h = rotate3d;  set(h,'ActionPreCallback','set(gcf,''windowbuttonmotionfcn'',@align_axislabel)');  set(h,'ActionPostCallback','set(gcf,''windowbuttonmotionfcn'','''')'); box on
% 
% subplot(4,3,12)
% 
% hold on
% for i = 1 : k
%     plot(Co_*1000,squeeze(pi_coqf(i,8,:)),'Linewidth',2)
% end
% plot([0,100],[.7,.7],'k--')
% text(90,.72,'MSL')
% xlabel('Concentration (mg/l)')
% zlabel('Marsh Depth (m)')
% ylabel('River Discharge')
% h = rotate3d;  set(h,'ActionPreCallback','set(gcf,''windowbuttonmotionfcn'',@align_axislabel)');  set(h,'ActionPostCallback','set(gcf,''windowbuttonmotionfcn'','''')'); box on
% 
% set(findobj('type','axes'),'fontsize',7)
% h_fig = gcf;
% set(h_fig,'PaperOrientation','portrait')
% set(h_fig,'PaperPosition', [0 0 7 7]) % [... ... max_width = 7.5 max_height = 9]
% tit = 'Co vs pars';
% print(tit,'-dtiff','-r400')
% movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\22 - Nondimensional Analysis')
% %close all 
% 
% %%
% % width and depths vs concentration
% R_ = (1:8);
% 
% figure
% 
% subplot(4,3,1)
% 
% hold on
% for i = 1 : k
%     plot(R_,squeeze(pi_bfm_r(i,3,:)),'Linewidth',2)
% end
% xlabel('Sea Level Rise Rate (mm/yr)')
% zlabel('TF Width (m)')
% ylabel('Basin Width')
% h = rotate3d;  set(h,'ActionPreCallback','set(gcf,''windowbuttonmotionfcn'',@align_axislabel)');  set(h,'ActionPostCallback','set(gcf,''windowbuttonmotionfcn'','''')'); box on
% 
% subplot(4,3,2)
% 
% hold on
% for i = 1 : k
%     plot(R_,squeeze(pi_bfm_r(i,4,:)),'Linewidth',2)
% end
% plot([0,10],[.7,.7],'k--')
% text(8,.72,'MSL')
% xlabel('Sea Level Rise Rate (mm/yr)')
% zlabel('Tidal Flat Depth (m)')
% ylabel('Basin Width')
% h = rotate3d;  set(h,'ActionPreCallback','set(gcf,''windowbuttonmotionfcn'',@align_axislabel)');  set(h,'ActionPostCallback','set(gcf,''windowbuttonmotionfcn'','''')'); box on
% 
% subplot(4,3,3)
% 
% hold on
% for i = 1 : k
%     plot(R_,squeeze(pi_bfm_r(i,5,:)),'Linewidth',2)
% end
% plot([0,10],[.7,.7],'k--')
% text(8,.72,'MSL')
% xlabel('Sea Level Rise Rate (mm/yr)')
% zlabel('Marsh Depth (m)')
% ylabel('Basin Width')
% h = rotate3d;  set(h,'ActionPreCallback','set(gcf,''windowbuttonmotionfcn'',@align_axislabel)');  set(h,'ActionPostCallback','set(gcf,''windowbuttonmotionfcn'','''')'); box on
% 
% subplot(4,3,4)
% 
% hold on
% for i = 1 : k
%     plot(R_,squeeze(pi_le_r(i,3,:)),'Linewidth',2)
% end
% xlabel('Sea Level Rise Rate (mm/yr)')
% zlabel('TF Width (m)')
% ylabel('Estuary Length')
% h = rotate3d;  set(h,'ActionPreCallback','set(gcf,''windowbuttonmotionfcn'',@align_axislabel)');  set(h,'ActionPostCallback','set(gcf,''windowbuttonmotionfcn'','''')'); box on
% 
% subplot(4,3,5)
% 
% hold on
% for i = 1 : k
%     plot(R_,squeeze(pi_le_r(i,4,:)),'Linewidth',2)
% end
% plot([0,10],[.7,.7],'k--')
% text(8,.72,'MSL')
% xlabel('Sea Level Rise Rate (mm/yr)')
% zlabel('Tidal Flat Depth (m)')
% ylabel('Estuary Length')
% h = rotate3d;  set(h,'ActionPreCallback','set(gcf,''windowbuttonmotionfcn'',@align_axislabel)');  set(h,'ActionPostCallback','set(gcf,''windowbuttonmotionfcn'','''')'); box on
% 
% subplot(4,3,6)
% 
% hold on
% for i = 1 : k
%     plot(R_,squeeze(pi_le_r(i,5,:)),'Linewidth',2)
% end
% plot([0,10],[.7,.7],'k--')
% text(8,.72,'MSL')
% xlabel('Sea Level Rise Rate (mm/yr)')
% zlabel('Marsh Depth (m)')
% ylabel('Estuary Length')
% h = rotate3d;  set(h,'ActionPreCallback','set(gcf,''windowbuttonmotionfcn'',@align_axislabel)');  set(h,'ActionPostCallback','set(gcf,''windowbuttonmotionfcn'','''')'); box on
% 
% subplot(4,3,7)
% 
% hold on
% for i = 1 : k
%     plot(R_,squeeze(pi_rco(i,3,:)),'Linewidth',2)
% end
% xlabel('Concentration (mg/l)')
% zlabel('TF Width (m)')
% ylabel('Ocean Concentration')
% h = rotate3d;  set(h,'ActionPreCallback','set(gcf,''windowbuttonmotionfcn'',@align_axislabel)');  set(h,'ActionPostCallback','set(gcf,''windowbuttonmotionfcn'','''')'); box on
% 
% subplot(4,3,8)
% 
% hold on
% for i = 1 : k
%     plot(R_,squeeze(pi_rco(i,4,:)),'Linewidth',2)
% end
% plot([0,10],[.7,.7],'k--')
% text(8,.72,'MSL')
% xlabel('Sea Level Rise Rate (mm/yr)')
% zlabel('Tidal Flat Depth (m)')
% ylabel('Ocean Concentration')
% h = rotate3d;  set(h,'ActionPreCallback','set(gcf,''windowbuttonmotionfcn'',@align_axislabel)');  set(h,'ActionPostCallback','set(gcf,''windowbuttonmotionfcn'','''')'); box on
% 
% subplot(4,3,9)
% 
% hold on
% for i = 1 : k
%     plot(R_,squeeze(pi_rco(i,5,:)),'Linewidth',2)
% end
% plot([0,10],[.7,.7],'k--')
% text(8,.72,'MSL')
% xlabel('Sea Level Rise Rate (mm/yr)')
% zlabel('Marsh Depth (m)')
% ylabel('Ocean Concentration')
% h = rotate3d;  set(h,'ActionPreCallback','set(gcf,''windowbuttonmotionfcn'',@align_axislabel)');  set(h,'ActionPostCallback','set(gcf,''windowbuttonmotionfcn'','''')'); box on
% 
% subplot(4,3,10)
% 
% hold on
% for i = 1 : k
%     plot(R_,squeeze(pi_rcf(i,3,:)),'Linewidth',2)
% end
% xlabel('Sea Level Rise Rate (mm/yr)')
% zlabel('TF Width (m)')
% ylabel('River Concentration')
% h = rotate3d;  set(h,'ActionPreCallback','set(gcf,''windowbuttonmotionfcn'',@align_axislabel)');  set(h,'ActionPostCallback','set(gcf,''windowbuttonmotionfcn'','''')'); box on
% 
% subplot(4,3,11)
% 
% hold on
% for i = 1 : k
%     plot(R_,squeeze(pi_rcf(i,4,:)),'Linewidth',2)
% end
% plot([0,10],[.7,.7],'k--')
% text(8,.72,'MSL')
% xlabel('Sea Level Rise Rate (mm/yr)')
% zlabel('Tidal Flat Depth (m)')
% ylabel('River Concentration')
% h = rotate3d;  set(h,'ActionPreCallback','set(gcf,''windowbuttonmotionfcn'',@align_axislabel)');  set(h,'ActionPostCallback','set(gcf,''windowbuttonmotionfcn'','''')'); box on
% 
% subplot(4,3,12)
% 
% hold on
% for i = 1 : k
%     plot(R_,squeeze(pi_rcf(i,5,:)),'Linewidth',2)
% end
% plot([0,10],[.7,.7],'k--')
% text(8,.72,'MSL')
% xlabel('Sea Level Rise Rate (mm/yr)')
% zlabel('Marsh Depth (m)')
% ylabel('River Concentration')
% h = rotate3d;  set(h,'ActionPreCallback','set(gcf,''windowbuttonmotionfcn'',@align_axislabel)');  set(h,'ActionPostCallback','set(gcf,''windowbuttonmotionfcn'','''')'); box on
% 
% set(findobj('type','axes'),'fontsize',7)
% h_fig = gcf;
% set(h_fig,'PaperOrientation','portrait')
% set(h_fig,'PaperPosition', [0 0 7 7]) % [... ... max_width = 7.5 max_height = 9]
% tit = 'R vs pars';
% print(tit,'-dtiff','-r400')
% movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\22 - Nondimensional Analysis')
% %close all 