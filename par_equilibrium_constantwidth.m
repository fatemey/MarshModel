function par_equilibrium_constantwidth
% par_equilibrium_constantwidth: plots width vs depth equlibrium data for a specific Co
%
% Last Update: 9/20/2017
%
%--------------------------------------------------------------------------------------------------
clear

load co_data_1000yr_TM_SS
co = dat(:,1);
y1 = dat(:,2);
y2 = dat(:,3);
y3 = dat(:,4);

clear dat
dat = [co,y1];
for i =7 : length(co)
    i
    
    [t, y] = BoxModel(co(i),y1(i),y2(i),y3(i));
    dat(i,3) = y(end,2);
    dat(i,4) = y(end,3);
    
    plot_BoxModel(t,y)
    h_fig=gcf;
    set(h_fig,'PaperOrientation','portrait')
    set(h_fig,'PaperPosition', [0 0 7 7]) % [... ... max_width=7.5 max_height=9]
    tit = ['2_Co = ' ,num2str(co(i)*1000), ' mg l-1'];
    print(tit,'-dtiff','-r400')
    movefile([tit,'.tif'],'C:\Users\fy23\Fateme\Projects\Marsh Model\Results\21 - Equlibrium depth results')
    close all
    
end

save par_equildata_ctewidth_5kyr_2 dat

end

