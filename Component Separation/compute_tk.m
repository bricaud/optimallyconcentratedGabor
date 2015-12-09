function [tk]= compute_tk(tprec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% USE: 
%    [tk]= compute_tk(tprec)
%   
% 
% FUNCTION:
%   Compute the next t-parameter for FISTA algorithm
% 
% 
% INPUT: 
%       tprec: current value of t-parameter
%
%
% OUTPUT:
%       tk: Next value of t-parameter
%       
%
% NOTA:
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENESIS S.A. - www.genesis.fr - January 2011 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tk=(1+sqrt(1+4*(tprec^2)))/2;
end