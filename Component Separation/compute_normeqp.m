function [norme] = compute_normeqp(vect,p1,q1,group1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% USE: 
%   [norme] = compute_normeqp(vect,p1,q1,group1) 
%   
% 
% FUNCTION:
%   Compute the mixed norm qp of a vector 
% 
% 
% INPUT: 
%       vect: vector
%       signal_ton: current temporal tonal component
%       p1,q1: regularization norm
%       group1: groups for the mixed norm
%
% OUTPUT:
%       obj1: Regression term
%       obj2: Regularization term on the tonal part (norm 1)
%       obj3: Regularization term on the transient part (mixed norm 2,1)
%       
%
% NOTA: Ponderations should probably be taken in account in the calculus. TODO 
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENESIS S.A. - www.genesis.fr - January 2011 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
norme=0;
for iter = 1:size(group1),
    norme = norme + (sum((abs(vect(group1{iter}))).^p1))^(q1/p1);
end
    norme=norme^(1/q1);
