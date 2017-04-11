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
% correct data removal:
%substitue ~isnan(ind)  with ~isempty(ind) && length(ind)>1
