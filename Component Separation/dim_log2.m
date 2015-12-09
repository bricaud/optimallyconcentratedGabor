function [y Ncoeffnonnul] =dim_log2(alpha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% USE: 
%   [y Ncoeffnonnul] = dim_log2(alpha)
% 
% FUNCTION:
%   Compute the logarithmic dimension of a complex vector 
% 
% INPUT:
%       alpha: Vector from the spectrogram
%    
% 
% OUTPUT:
%       y: logarithmic dimension
%       Ncoeffnonnul: number of non zero coefficients
%         
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENESIS S.A. - www.genesis.fr - January 2011 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alphanonnul = find( alpha );
Ncoeffnonnul = length(alphanonnul(:));

N=size(alpha,1)*size(alpha,2);
alpha( alpha == 0) = eps;

alpha=log2(abs(alpha(:)).^2);

y=sum(alpha)/N;




