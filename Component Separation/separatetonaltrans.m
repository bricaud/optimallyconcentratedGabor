function [tonale transitoire bruit] = separatetonaltrans( filename )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% USE: 
%   [tonale transitoire bruit] = separatetonaltrans(filename)
% 
% FUNCTION:
%   Decomposition of a wav signal into tonal, transient and noisy
%   components (indeed transient extraction and denoising of the remaining
%   signal)
% 
% 
% INPUT:
%       filename: Complete filename of the input wav signal
%    
% 
% OUTPUT:
%       tonale   : tonal component of the signal
%       transitoire: transient component of the signal
%       bruit: noisy component of the signal
%
%
% NOTA: Decomposition is done using a sparse regression, with a LASSO (norm
% 1) regularization for the tonal part and a G-LASSO regularization on the
% transient part. Optimization uses a slightly modified FISTA (Fast
% iterative shrinkage-threshold algorithm) strategy.
%  
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENESIS S.A. - www.genesis.fr - January 2011 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialisations

%Loading the signal
[signal sample_rate nbits]=wavread(filename);


%Initialization of Time Frequency parameters (we use a STFT)
%For transients, we choose short windows (in order to get an high accuracy
%in time). To obtain high frequency accuracy, for the tonal part we chose
%longer windows. 

%Size of windows for STFT (transient part): around 10 ms
parametres.taille_fenetre_trans = 2^round(log2(sample_rate/100));

%Size of windows for STFT (tonal part): around 100 ms
parametres.taille_fenetre_ton = 2^round(log2(sample_rate/10));

%Overlap ratio for STFT
parametres.noverlap_ratio = 0.5;
parametres.fe = sample_rate;


%Stopping condition of the optimization
epsilon=0.05;

%Noise ponderation
noisepond=0.8;

%NB of different ponderations by sec
splitbysec = 2;

%Standard divisor
stdvalue = 1.85;


%% Performing the separation on the the FIRST channel of the signal
[tonale transitoire bruit]=fista(signal(1:end,1),1,1,2,1,noisepond,epsilon,parametres, splitbysec, stdvalue);



%% Write the results

%Each component is written in the same repertory than the original signal

wavwrite(tonale,sample_rate,nbits,[filename,'_tonale.wav']);
wavwrite(transitoire,sample_rate,nbits,[filename,'_transitoire.wav']);
wavwrite(bruit,sample_rate,nbits,[filename,'_bruit.wav']);