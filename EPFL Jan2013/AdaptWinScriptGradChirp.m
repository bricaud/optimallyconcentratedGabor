% ADAPTWINSCRIPTGRADCHIRP
% ==============
%
% Script for optimizing a window in terms of the concentration of
% a target function, estimated by restricting the gabor transform of a
% target (noisy) signal to a given time-frequency region
%
% It uses a gradient method
%
% add a constraint on the length of the window

%set(0,'DefaultAxesFontSize',18)


%% Optimization gradient parameters
epsilon = 1e-3;
alpha = 1;
p = 2.5;
set(0,'DefaultAxesFontSize',18)

%% generate 'chirp' signal
%y=real(pchirp(12000,-0.298));
NN=12000;
tau=1:NN;
alpha=1; % chirp slope
y=exp(1i*2*pi*alpha*tau.^2/(2*NN)).';
%sgram(y)
 FS = 44100;

% fade in et fade out
% Lbis = length(ybis);
L = length(y);
sizeHanning = round(FS/200)*2;
hanni =pgauss(sizeHanning,78);
hanni = hanni/max(hanni);
y = y.*[hanni(end - sizeHanning/2 +1:end)' ones((L-sizeHanning),1)' hanni(1:sizeHanning/2)'  ]';

% zeros around
y = [zeros(2*size(y,1),1); y; zeros(2*size(y,1),1)];
%ytest = [zeros(2*size(ytest,1),1); ytest; zeros(2*size(ytest,1),1)];

% param. init.
L = length(y);
M = L/4;
a = 100;
Lgamma = L;
g = psech(Lgamma); 

g = gabtight(g, a, Lgamma);

figure;
w = dgt(y, g, a, M); w = 20*log10(abs(w));
imagesc(w); axis xy;
m = max(max(w)); caxis([m-60 m]);
title('Gabor transform with initial window');

z = y;

%% Optimization gradient ascent
[gamma, crit, entrop] = WinOptimGrad_variance(z, g, a, M, L, p, epsilon, alpha);%


%% Visualization
axisX = (((1:length(y)/a)*a - a)/FS) ;
axisY = (0:FS/M:FS-1);
w2 = dgt(y, gamma, a, M);
w2 = 20 * log10(abs(w2));
figure; imagesc(axisX,axisY, w2); axis xy; m = max(max(w2)); caxis([m-60 m]);
title('Gabor Transform of the chirp with optimal window');
xlabel('Time (s)');
ylabel('Frequency (Hz)');

gg = fftshift(hanning(1024));
w2 = dgt(y, gg, a, M);%
w2 = 20 * log10(abs(w2));
figure; imagesc(axisX,axisY, w2); axis xy; m = max(max(w2)); caxis([m-60 m]);
title('Gabor Transform of the chirp with Gaussian window');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
w2 = dgt(gamma, fftshift(hanning(1024)), a, M);%, length(sig2(1:150000)));
w2 = fftshift(w2);
w2 = 20 * log10(abs(w2));

figure; imagesc(axisX,axisY, w2); axis xy; m = max(max(w2)); caxis([m-60 m]);
title('Ambiguity function of optimal window');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
