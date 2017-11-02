function plot_phi1

figure
clear

load SS_results_small

chi = dat(:,1);
Co = dat(:,4);
Cf = dat(:,5);
LE = dat(:,6);
Qf =dat(:,7);
R_ = dat(:,8);
bfm = 5 *10^3;  



x = Co./Cf;
y = LE/bfm;
z = chi/bfm;

c = Qf./(R_*bfm^2);
 [c_unique,I,J] = unique(c);
 col_map = parula(length(c_unique));
 col = col_map(J,:);
scatter3(x,y,z,[],col)
hold on
xqq=linspace(min(x),max(x));
yqq=linspace(min(y),max(y));
[xq,yq] = meshgrid(xqq,yqq);
vq = griddata(x,y,z,xq,yq);
colormap('jet')
mesh(xq,yq,vq)

xlabel('\pi_1 = Concentration Ratio')
zlabel('\pi =  Tidal Flat Width Ratio')
ylabel('\pi_2 = Basin Size Ratio')
title('Color : \pi_3 = Discharge Ratio')
% ylabel('\pi_3 = Discharge Ratio')
% title('Color : \pi_2 = Basin Size Ratio')

box on

end