function C=compute_stepsize(signal,spec_ton,spec_trans)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% USE: 
%   C=compute_stepsize(signal,spec_ton,spec_trans)
% 
% FUNCTION:
%   Compute an appropriate stepsize for the fista component separation
%   optimization
% 
% 
% INPUT:
%       signal: Temporal signal
%       spec_ton: Spectrogram (tonal view)
%       spec_trans: Spectrogram (transient view)
%      
%
% OUTPUT:
%       C: Stepsize for the gradient descent
%       
%
% NOTA: Computed according to the Gerd Teschke paper (Multi-frame
% representations in linear inverse problems with mixed multi-constraints)
%  
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENESIS S.A. - www.genesis.fr - January 2011 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coeff_ton=abs(spec_ton(:)).^2;
coeff_trans=abs(spec_trans(:)).^2;
B1=sqrt(sum(coeff_ton))/norm(signal(:),2);
B2=sqrt(sum(coeff_trans))/norm(signal(:),2);
C=1/(sqrt(sqrt(B1+B2)));
