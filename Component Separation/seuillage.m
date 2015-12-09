function [retour] = seuillage (x, seuil, p, q, group)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% USE: 
%   [retour] = seuillage (x, seuil, p, q, group)
%   
% 
% FUNCTION:
%   Performs the (soft) shrinkage step of the FISTA algorithm
% 
% 
% INPUT: 
%       x: spectrogram of a view
%       seuil: threshold of the shrinkage
%       p,q: regularization norm
%       group: groups for (possibly) the mixed norm
%
% OUTPUT:
%       retour: shrinked spectrogram
%       
%
% NOTA: Norm 1,2 must be (re)retested before its use (for example in
% E-LASSO). See Matthieu Kowalski thesis for more details
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENESIS S.A. - www.genesis.fr - January 2011 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Save the signs of the coefficients of x
retour = zeros(size(x));
signe=sign(x);

%Soft shrinkage of x for norm L1
seuil=repmat(seuil,size(x,1),1);
if((p==1) && (q==1))
    val_seuillee = abs(x)-seuil;
    retour(val_seuillee>0)=signe(val_seuillee>0).*val_seuillee(val_seuillee>0);
end

%Soft shrinkage of x for mixed norm 2,1
if((p==2) && (q==1))
    for iter=1:size(group)
        group_iter = group{iter};
        x_iter = x(group_iter);
        seuil_iter = seuil(group_iter);
        val_seuillee_iter = (1 - seuil_iter./norm(x_iter(1:end), 2) );
        val_seuillee_iter(val_seuillee_iter<0)=0;
        retour(group_iter)= x_iter.*val_seuillee_iter;
    
    end
end

%Soft shrinkage of x for mixed norm 1,2
if((p==1) && (q==2))
    for iter=1:size(group)
        s = ((seuil/(1+(numel(group{iter})*seuil)))*norm(x(group{iter}),1));
        val_seuillee(group{iter}) = abs(x(group{iter}))-s;
        retour(val_seuillee>0)=signe(val_seuillee>0).*val_seuillee(val_seuillee>0);

    end
end
end
