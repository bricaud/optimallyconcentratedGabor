function tight_win = compute_tight_window( window, n_step)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% USE: 
%   tight_win = compute_tight_window( window, n_step)
% 
% FUNCTION:
%   Compute a window for analysis and synthesis (use a strict frame)
% 
% 
% INPUT:
%       window: input window from which
%       n_step: step size between 2 frames
%      
%
% OUTPUT:
%       tight_win: window with strict frame
%
%
%
% NOTA: From the secret thesis of Florent Jaillet
%  
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENESIS S.A. - www.genesis.fr - January 2011 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Length of the window (no kidding)
window_length = length(window); 
%Number of overlaps
n_loop = floor( (window_length-1) / n_step); 

% Computation
window2 = window.^2;
tmp = window2;

for n = 1:n_loop,
    L_current = window_length - n*n_step;
    tmp(n*n_step+1:window_length) = tmp(n*n_step+1:window_length) + window2(1:L_current);
    tmp(1:L_current) = tmp(1:L_current) + window2(n*n_step+1:window_length);    
end;

tight_win = window ./ sqrt(tmp);


