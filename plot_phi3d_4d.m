function plot_phi3d_4d

figure
clear

load SS_results
rho_s = 1000;

%***
% Co = (10:10:40) *10^-3;
% Cf = (10:10:40) *10^-3;
% LE = (4:8) * 1000;
% Qf = 10:10:40;
% R_ = (2:6) *10^-3/365/24/60/60;
% T = (12:13) *60*60;
% H_ = [1.4,2.8] /2;
%***

chi = dat(:,1);
% chi = dat(:,1).*dat(:,2)./dat(:,3);
Co = dat(:,4);
Cf = dat(:,5);
LE = dat(:,6);
Qf =dat(:,7);
R_ = dat(:,8);
H_ = dat(:,9);
T = dat(:,10);
bfm = 5 *10^3;

for i = 1: 3
    
    hold on
    rem = H_ ~= .7 | T ~= T(1)|  R_~= R_(4*(i-1)+1);
    %     rem = H_ ~= 1.4 |  T ~= T(1) |  R_~= R_(9) | LE ~= 5000 | Qf ~=40 | Cf ~= .02 | Co ~= .04 ;
    chiv = chi(~rem);
    Cov = Co(~rem);
    Cfv = Cf(~rem);
    LEv = LE(~rem);
    Qfv = Qf(~rem);
    R_v = R_(~rem);
    H_v = H_(~rem);
    Tv = T(~rem);
    %
    %     %     x = (Cov.*(bfm*LEv.*H_v./Tv-Qfv)+Cfv.*Qfv)./R_v/bfm./LEv/rho_s;
    %         y = bfm*LEv.*H_v./Tv./Qfv;
    %     %     z = chiv/bfm;
    %     x= (Cov.*(bfm*LEv.*H_v./Tv-Qfv)+Cfv.*Qfv)./R_v/bfm./LEv/rho_s;
    % %     y = (bfm*LEv.*H_v./Tv-Qfv)./Qfv;
    %     z = chiv/bfm;
    %     c = H_v./Tv./R_v;
    
    x = Cov.*(bfm*LEv.*H_v./Tv-Qfv)./Cfv./Qfv;
    z = chiv/bfm;
    y = H_v./Tv./R_v;
    
    scatter3(x,y,z,'.')
    xqq=linspace(min(x),max(x));
    yqq=linspace(min(y),max(y));
    [xq,yq] = meshgrid(xqq,yqq);
    vq = griddata(x,y,z,xq,yq);
    n=length(xqq);
    R = linspace(63,255,n);  % Red
    G = linspace(77,80,n);   % Green
    B = linspace(184,8,n);  % Blue
    cmap =  [R', G', B']/255;
    %     colormap(cmap);
    colormap('cool')
%     mesh(xq,yq,vq)
%     grid off
    
end

xlabel('\pi_1 = Mass Ratio')
ylabel('\pi_2 = Discharge Ratio')
zlabel('\pi =  Tidal Flat Width Ratio')

box on

% set(findobj('type','axes'),'fontsize',10)
% h_fig=gcf;
% set(h_fig,'PaperOrientation','portrait')
% set(h_fig,'PaperPosition', [0 0 6 5]) % [... ... max_width=7.5 max_height=9]
% tit='nond2';
% print(tit,'-dtiff','-r400')
% movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\24 - Final Results')
% close all

end