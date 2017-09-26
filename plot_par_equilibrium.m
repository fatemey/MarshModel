clear
clf

load par_equil_data_5kyr

col = linspecer(length(co));

hold on
for i = 1 : length(co)
    h(i) = scatter(depth(i,:),width(i,:)/1000,20,col(i,:),'o','filled');
end

xlabel('Tidal Flat Depth (m)')
ylabel('Tidal Flat Width (km)')
xlim([0.6 1.6])
legend([h(20),h(1)],'Co=100 mg/l','Co=5 mg/l','Location','NorthWest')
box on
% 
% h_fig=gcf;
% set(h_fig,'PaperOrientation','portrait')
% set(h_fig,'PaperPosition', [0 0 3 3]) % [... ... max_width=7.5 max_height=9]
% tit = 'width_depth_equil_5kyr_2';
% print(tit,'-dtiff','-r400')
% movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\20 - Unstable Equlibrium Results through Optimization of Steady State')
% close all