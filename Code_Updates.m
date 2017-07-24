% Code updates: 
% -----------------------------------------------------------
% 2/28/2017:
% Consider total b_fm instead of only half.
% remove the plot data after TF conversion to marsh.
%
% 3/10/2017:
% remove the plot data after marsh conversion to TF.
%
% 3/13/2017: (disregard because the changes were deleted.)
% Adding:  && length(ind)>1, to consider a case when conversion happens for
% more than one value.

% 3/28/2017: 
% function naming was changed.

% 4/11/2017
% correct data removal condition:
% substitute ~isnan(ind)  with ~isempty(ind) && length(ind)>1

% 7/17/17
% delete extra unused code lines
%         if C_o == 10*10^-3 % condition to consider alpha based on D'Alpaos et al 2011 approach
%             alpha = 4 *10^-3/365/24/60/60; % (m/s)
%         elseif C_o == 20*10^-3
%             alpha = 8 *10^-3/365/24/60/60;
%         elseif C_o == 100*10^-3
%             alpha = 38 *10^-3/365/24/60/60;
%         end
%         alpha = 0.38 * C_r;
%         TF_accretion = alpha*dt*(y(2)/H); % based on D'Alpaos et al 2011 approach during one time step
%         TF_deposition = alpha*dt*(y(2)/H)*local*L_E*rho_s;
