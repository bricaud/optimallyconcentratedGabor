function [ alpha ] = update_opt( alpha,upd_vect)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% USE: 
%   [ alpha ] = update_opt( alpha,upd_vect)
%   
% 
% FUNCTION:
%   Performs a gradient descent step
% 
% 
% INPUT: 
%       vect: vector
%       upd_vect: Update vector (actually the gradient)
%
% OUTPUT:
%       alpha: vector updated
%       
%
% NOTA:
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENESIS S.A. - www.genesis.fr - January 2011 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha=alpha-upd_vect;
