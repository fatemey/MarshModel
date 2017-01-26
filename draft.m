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
