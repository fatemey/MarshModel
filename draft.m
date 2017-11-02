% codes that are no longer needed but may be needed in the future
%---------------------------------------------------------------------
% time variant solution of SOM
% B=[];
%         syms B(t)
%         eqn = diff(B,t) == r*B*(1-B/B_max)-m*B;
%         cond = B(0) == B_max*y0(3)/H;
%         Sol(t) = dsolve(eqn,cond);
%         B = Sol(tspan);
%         Sol = dsolve(eqn); B = Sol(2);
%         B(t) = (B_max*y(3) + B_max*y(3)*tanh(B_max*y(3)*(t/(2*B_max*H) + atanh(y(3)/H)/(B_max*y(3)))))/(H + y(3)); %or B(t) = (B_max + B_max*tanh(B_max*(atanh(z^2 + z - 1)/B_max + t/(2*B_max*z))))/(z + 1);
%---------------------------------------------------------------------
% print graphs
% casen = 22;
print(figure(1),['Cr-casen',num2str(casen)],'-dtiff','-r400')
movefile(['Cr-casen',num2str(casen),'.tif'],'C:\Users\FY23\Fateme\ex\Work\Model\Results\res2')
print(figure(2),['bf-casen',num2str(casen)],'-dtiff','-r400')
movefile(['bf-casen',num2str(casen),'.tif'],'C:\Users\FY23\Fateme\ex\Work\Model\Results\res2')
print(figure(3),['dm-casen',num2str(casen)],'-dtiff','-r400')
movefile(['dm-casen',num2str(casen),'.tif'],'C:\Users\FY23\Fateme\ex\Work\Model\Results\res2')
print(figure(4),['df-casen',num2str(casen)],'-dtiff','-r400')
movefile(['df-casen',num2str(casen),'.tif'],'C:\Users\FY23\Fateme\ex\Work\Model\Results\res2')
close all
%---------------------------------------------------------------------
% set the number of axis ticks
NumTicks = 6;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))
%also
set(gca,'YTick',[0, 5:10 :51]) %label specification
set(gca,'YTickLabel','0.45|1|2|3|4|5')
%---------------------------------------------------------------------
% plotting
col(:,1) = ones(size(dat,1),1)/2;
col(:,3) = dat(:,1)/max(dat(:,1));
scatter(dat(:,2),dat(:,3),50,col,'o','filled')

n_clusters = 10;
circle_siz = ones(size(dat,1),1);
siz_step = linspace(0,max(fetch_c),n_clusters-1);
for i = 1 : length(siz_step)-1
    circle_siz(fetch_c>siz_step(i)&fetch_c<=siz_step(i+1)) = 30*i;
end

%----------------------------------------------------------------------
% linear solve for independent variables need for nondimensionla analysis
%       [alpha
% A *   beta      = B
%         gamma]
A=[ 1 0 0 ;... % each row corresponds to 1 variable/main dimension
       -3 3 1 ;...
      0 0 -1 ];
B=[ 0 ;...
       3 ;...
      -1 ];
format rat
linsolve(A,B)