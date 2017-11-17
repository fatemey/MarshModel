function plot_phi

format longG
clear

%----------------------------- Load Data ----------------------------
load SS_results

bf = dat(:,1);
df = dat(:,2);
dm = dat(:,3);
Co = dat(:,4);
Cf = dat(:,5);
LE = dat(:,6);
Qf =dat(:,7);
R = dat(:,8);
H = dat(:,9);
T = dat(:,10);
bfm = 5 *10^3;
% bfm = dat(:,11);
Vw = dat(:,12);
rho_s = 1000;

%--------------------- Select Data of Interest ---------------------
rem = H ~= .7 |  T ~= T(1)  ;% |  Vw ~= 10 ;%  | Co ~= .02;% | R~= R(1);% | Co ~= .02 ;%| LE ~= 4000 | Qf ~=40 | Cf ~= .01 | Co ~= .02 ; UniqueR=[1;13;25;37;49]
% rem =zeros(size(Co));

bfv = bf(~rem);
dfv = df(~rem);
dmv = dm(~rem);
Cov = Co(~rem);
Cfv = Cf(~rem);
LEv = LE(~rem);
Qfv = Qf(~rem);
Rv = R(~rem);
Hv = H(~rem);
Tv = T(~rem);
Vwv = Vw(~rem);

%--------------------- Compute Pi Terms ---------------------
C_ratio = Cov./Cfv;
M_ratio = Cov.*(bfm*LEv.*Hv./Tv-Qfv)./Cfv./Qfv; % mass of ocean sediment to river 
M_ratio_tot = (Cov.*(bfm*LEv.*Hv./Tv-Qfv)+Cfv.*Qfv)./Rv/bfm./LEv/rho_s; % mass total to mass filled in the space created by SLR
M_ratio_totacc = (Cov.*(LEv.*((dmv.*bfv)+((bfm-bfv).*dmv))./Tv-Qfv)+Cfv.*Qfv)./Rv/bfm./LEv/rho_s; % accurate mass total to mass filled in the space created by SLR
Q_ratio = (bfm*LEv.*Hv./Tv-Qfv)./Qfv; % ocean discharge to river
Q_ratio_acc = (LEv.*((dmv.*bfv)+((bfm-bfv).*dmv))./Tv-Qfv)./Qfv;  % accurate ocean discharge to river
Wind_ratio = Vwv./Rv;
Rise_ratio = Hv./Tv./Rv;
fetch_nond = bfv/bfm;

%--------------------- Select Pi Terms ---------------------
x = C_ratio;
% x = M_ratio;
% x = C_ratio./Q_ratio;
x = M_ratio_tot;

y = Q_ratio;
% y = fetch_nond;

z = fetch_nond;

% c = Rise_ratio;
c = Rv;
[c_unique,I,J] = unique(c);
col_map = winter(length(c_unique));
col = col_map(J,:);

%---------------------------- Plot ----------------------------
% figure
clf

% scatter(x,y)
% scatter(x,y,[],col)
% scatter3(x,y,z)
scatter3(x,y,z,[],col)

%--- add surface plot
% hold on
% xqq=linspace(min(x),max(x));
% yqq=linspace(min(y),max(y));
% [xq,yq] = meshgrid(xqq,yqq);
% vq = griddata(x,y,z,xq,yq);
% colormap('jet')
% % n=length(xqq);
% % R = linspace(63,255,n);
% % G = linspace(77,80,n);
% % B = linspace(184,8,n);
% % cmap =  [R', G', B']/255;
% % colormap(cmap);
% mesh(xq,yq,vq)

% zlim([0,1])

%------------------------ Select Plot Labels ------------------------
% xlabel('Concentration Ratio')
xlabel('Mass Ratio')
% xlabel('Concentration Ratio to Discharge Ratio')

% ylabel('Basin Size Ratio')
ylabel('Discharge Ratio')
% ylabel('Nondimensional Fetch')

zlabel('Nondimensional Fetch')

% title('Color : Discharge Ratio')
% title('Color : Water Level Rise Rate Ratio')
% title('Color : Wind Velocity to SLR Rate Ratio')

box on

%---------------------------- Save Plot ----------------------------
% set(findobj('type','axes'),'fontsize',10)
% h_fig=gcf;
% set(h_fig,'PaperOrientation','portrait')
% set(h_fig,'PaperPosition', [0 0 6 5]) % [... ... max_width=7.5 max_height=9]
% tit='nond2';
% print(tit,'-dtiff','-r400')
% movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\24 - Final Results')
% close all