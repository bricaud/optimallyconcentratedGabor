
function [sig_rec, time] = ispecgram(STFTcoeffs, FS, wind_rec, nOverlap)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% USE: 
%   [sig_rec, time] = ispecgram(STFTcoeffs, FS, wind_rec, nOverlap)
%   
% 
% FUNCTION:
%   Performs the inverse transform of a STFT
% 
% 
% INPUT:
%       STFTcoeffs: coefficient matrix
%       FS: sample frequency
%       Window: Synthesis window
%       nOverlap: overlap of the windows
%
% OUTPUT:
%       sig_rec: temporal signal
%       time: time vector
%       
%
% NOTA: STFTScoeffs are coefficients computed by spectrogram matlab
% standard function
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENESIS S.A. - www.genesis.fr - January 2011 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Sizes of coefficient matrix and of the window used to split the signal
NTime = size(STFTcoeffs, 2);
Nwind = length(wind_rec);
Nshift = Nwind - nOverlap; % number of samples between each IFFT

SigLength = Nwind + Nshift * (NTime-1);
sig_rec = zeros(SigLength, 1);

% Compute the first time segment

sig_rec(1:Nwind) = doIFFT(STFTcoeffs(:,1), Nwind).* wind_rec;



% IFFT by window
for cnt = 2 : NTime-1,
   % Index for the beginning and the end of current window
   nBegin = (cnt-1) * Nshift + 1;
   nEnd = nBegin + Nwind - 1; 
   % Computation of IFFT for the current segment(spectrum)
   sigtemp = doIFFT( STFTcoeffs(:,cnt), Nwind);
   % rebuild the signal (do not forget the overlap) 
   sig_rec(nBegin : nEnd) = sig_rec(nBegin : nEnd) + sigtemp .* wind_rec;
 
end;

%Last inverse window 
nBegin = (NTime-1) * Nshift + 1;
nEnd = nBegin + Nwind - 1;
% Computation of IFFT for the last segment(spectrum)
sig_end = doIFFT( STFTcoeffs(:,NTime), Nwind);
% rebuild the signal (do not forget the overlap) 
sig_rec(nBegin : nEnd) = sig_rec(nBegin : nEnd) + sig_end .* wind_rec;
  

%Temporal sampling
Nsig_rec = length(sig_rec);
time = linspace(0, (Nsig_rec-1)/FS, Nsig_rec);
   
%-----------------------------------------------------------------------------

function sig = doIFFT(inFreq, Nwind)
% realise la transformer de fourier inverse du morceau de signal contenu
% dans infreq a la longeur esxitger par la fenétre
%ENTREE  :   infreq : contient la transformer de fourrier d'un segment (localisation en temps)du signal
%            Nwind  : longeur du segment 

NFreq = length(inFreq);
%reconstituion du spectre complet grace à la periodciter de celui çi
temp = [inFreq(1:(NFreq-1)) ; conj(flipud(inFreq(2:NFreq)))];
temp = ifft(temp);
sig = temp(1:Nwind);

return
