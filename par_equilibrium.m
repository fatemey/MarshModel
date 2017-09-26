function par_equilibrium
% par_equilibrium: plots width vs depth equlibrium data for a specific Co
%
% Last Update: 9/20/2017
%
%--------------------------------------------------------------------------------------------------
clear

load co_data_1000yr_TM_SS_2
co = dat(:,1);
y1 = dat(:,2);
y2 = dat(:,3);
y3 = dat(:,4);

d_inc = .7:.1:1.4; %set of incerements for depth (m)

for i = 1 : length(co)
    
    i
    [t, y] = BoxModel(co(i),y1(i),y2(i),y3(i));
    
    %     plot_BoxModel(t,y)
    
    w = y(:,1); % tidal flat width
    d = y(:,2); % tidal flat depth
    
    %     scatter(w,d)
    
    for j =1 : length(d_inc)
        
        x = d>d_inc(j);
        unique_x = x(diff([-1;x])~=0);
        
        if ~ismember(0,x)
            depth(i,j) = NaN;
            width(i,j) = NaN;
            continue
        end
        
        if ~ismember(1,x)
            depth(i,j) = NaN;
            width(i,j) = NaN;
            continue
        end
        
        ind = find((diff([-1;x])~=0)==1);
        if length(unique_x)==2
            depth(i,j) = d(ind(2));
            width(i,j) = w(ind(2));
        elseif length(unique_x)==3
            depth(i,j) = d(ind(3));
            width(i,j) = w(ind(3));
        else
            i
            j
        end
        
    end
    
end

save par_equil_data_5kyr_2 co depth width

end

