function [pond_ton pond_trans] = compute_ponderation(signal_spec_ton,signal_spec_trans,nbsplits,taille_fenetre_ton,taille_fenetre_trans,noverlap_ratio,fe,noisepond,stdvalue)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% USE: 
%   [pond_ton pond_trans] = compute_ponderation(signal_spec_ton,signal_spec_trans,nbsplits,taille_fenetre_ton,taille_fenetre_trans,noverlap_ratio,fe,noisepond,stdvalue),
%   
% 
% FUNCTION:
%   Compute two regularization vector for LASSO/GLASSO problem, the first
%   one for the tonal part, second for transient. For each vector, there is
%   one coefficient by frame in the associated spectrogram
% 
% 
% INPUT:
%       signal_spec_ton: tonal representation of the signal (accurate in
%       frequency)
%       signal_spec_trans: transient representation of the signal (accurate in
%       time)
%       nbsplits: Number of different splits for which we will compute a
%       regularization parameter
%       taille_fenetre_ton: Size of window used to compute tonal spectrogram 
%       taille_fenetre_trans: Size of window used to compute transient
%       spectrogram 
%       noverlap_ratio: Overlap ratio for tonal and transient STFT
%       fe: sample frequency
%       noisepond: noise ponderation
%       stdvalue: standard divisor for tonal/transient estimator ponderation
%
% OUTPUT:
%       pond_ton   : vector of ponderation for tonal spectrogram
%       pond_trans : vector of ponderation for transient spectrogram
%       
%
% NOTA:
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENESIS S.A. - www.genesis.fr - January 2011 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialisation of ponderation vectors
pond_ton=zeros(1,size(signal_spec_ton,2));
pond_trans=zeros(1,size(signal_spec_trans,2));


%% Computation the ponderation

%Make the splits
for i=0:nbsplits-1

%In the tonal and the transient views, and for each split, compute the first and the last
%sample of the split
    debut_ton=1+ceil(i*size(signal_spec_ton,2)/nbsplits);    
    fin_ton=ceil((i+1)*size(signal_spec_ton,2)/nbsplits);
    debut_trans=1+ceil(i*size(signal_spec_trans,2)/nbsplits);    
    fin_trans=ceil((i+1)*size(signal_spec_trans,2)/nbsplits);

    %Average ratio of transient
    transient_ratio=0;
    tonal_ratio=0;
        
    %Compute the tonal/transitory ratio for each frame in the split
    for frame=debut_ton:fin_ton
    

        %Compute the center sample of the frame
        center_sample_tonal = ((1-noverlap_ratio)*(frame-1))* taille_fenetre_ton + 0.5*taille_fenetre_ton;
        
        %Compute the corresponding sample in the transient view
        pos_center_sample_trans=floor(center_sample_tonal/(taille_fenetre_trans*(1-noverlap_ratio)));
        
        %Compute the ratio for the current tonal frame (we chose enough
        %transient frames to get as much coefficients for the two views)
        
        y1 = dim_log2(signal_spec_ton(1:end-1,frame));
        y2 = dim_log2(reshape(signal_spec_trans(1:end-1,pos_center_sample_trans-taille_fenetre_ton/(2*taille_fenetre_trans)+1:pos_center_sample_trans+taille_fenetre_ton/(2*taille_fenetre_trans)),taille_fenetre_ton/2,1));
        [Iton,Itrans]=proportion_tonaltransient(y1,y2);
          
        transient_ratio = transient_ratio + Itrans;
        tonal_ratio = tonal_ratio + Iton;
        
    end;
        
    %Averaging transient ration
    transient_ratio = transient_ratio/(fin_ton-debut_ton+1);
    tonal_ratio = tonal_ratio/(fin_ton-debut_ton+1);
        
% Compute the ponderation according to the tonal/transient estimation

%The choice of stdvalue seems to work rather well with 1.85 well
pond_trans(1,debut_trans:fin_trans)=noisepond*(tonal_ratio/(stdvalue));
pond_ton(1,debut_ton:fin_ton)=noisepond*(1-tonal_ratio/(stdvalue));

end
end