function [frequencies, amplitudes, n] = findMainHarmonics(spectro, mode, threshold)
%For a given Time-Frequency representation, finds the predominant harmonic
%   present in it.
% 
%   findMainHarmonics(spectro, mode)
%       inputs:
%           spectro: a Time Frequency representation of the signal
%           mode: 'main' -> finds only the main harmonic
%                 'all' -> finds all harmonics
%           threshold: threshold under which spectro values are considered
%              to be zero 
%           N: maximum number of peaks to detect
%
%       outputs:
%           frequencies: the list of found frequencies
%           amplitudes: the list or corresponding amplitudes
%           n: the maximum number of found frequencies.
%               should be 1 if 'main'
%               should be size(frequencies,1) otherwise

spectro = spectro(1:size(spectro,1)/2,:);

%% Mode = 'main'
if strcmp(mode,'main'),
    [amplitudes,frequencies] = max(spectro);
    n = 1;    
    spectro(spectro == 0) = NaN;
%% Mode = 'all'
else
    frequencies = [];
    n = length(frequencies);
    for i = 1:size(spectro,2),
        [ampls,freqs] = findpeaks(spectro(:,i),'minpeakheight',threshold,'npeaks',10);
        nfreq = length(freqs);
        if isempty(freqs),
            freqs = NaN(1,n);
            ampls = NaN(1,n);
        elseif nfreq > n,
            frequencies = [frequencies ; NaN(nfreq-n , size(frequencies,2))];
            amplitudes = [amplitudes ; NaN(nfreq-n , size(amplitudes,2))];
            n = nfreq;
        elseif nfreq < n,
            freqs = [freqs , NaN(1,n-nfreq)];
            ampls = [ampls , NaN(1,n-nfreq)];
        end
        frequencies = [frequencies , freqs'];
        amplitudes = [amplitudes , ampls'];
    end
end


