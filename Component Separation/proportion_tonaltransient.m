function [Iton,Itrans]=proportion_tonaltransient(dim1,dim2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% USE: 
%  [Iton,Itrans]=proportion_tonaltransient(dim1,dim2)
%   
% 
% FUNCTION:
%   Estimate the proportion of tonal and transient coefficients in the
%   signal
% 
% 
% INPUT:
%       dim1: logarithmic dimension of tonal coefficients
%       dim2: logarithmic dimension of transient coefficients
%
% OUTPUT:
%       Iton   : estimated proportion of tonal coefficients
%       Itrans : estimated proportion of transient coefficients
%       
%
% NOTA:
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENESIS S.A. - www.genesis.fr - January 2011 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Iton=2^dim1/(2^dim1+2^dim2);
Itrans=2^dim2/(2^dim1+2^dim2);