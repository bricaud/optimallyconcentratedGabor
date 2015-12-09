function [yk] = compute_yk(tprec, xact, xprec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% USE: 
%   [yk_usual] = compute_yk(tprec, xact, xprec)
%   
% 
% FUNCTION:
%   Compute the next calculation point for the gradient descent
% 
% 
% INPUT: 
%       tprec: current value of t-parameter
%       xact: Current calculation point for the gradient
%       xprec: Last calculation point for the gradient
%
% OUTPUT:
%       yk_usual: Next calculation point for the gradient
%       
%
% NOTA:
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENESIS S.A. - www.genesis.fr - January 2011 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calcule le prochain point de calcul du gradient, à partir des 2 points de
% calcul précédents
yk = xact + ((tprec-1)/(compute_tk(tprec))).*(xact-xprec);
yk_usual = xact;
end
