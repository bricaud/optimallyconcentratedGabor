% DEMO_ADAPTWIN
% ==============
%
% Script for optimizing a window in terms of the concentration of
% a target function, estimated by restricting the gabor transform of a
% target (chirped gaussian) signal to a given time-frequency region
%
% It uses a gradient method
%
% add a constraint on the length of the window

set(0,'DefaultAxesFontSize',18)

%% Target signal
% Example 1: a bended and 'chirped' gaussian
NN=12000;
s=0.1; % shear or slope of the chirp
e=18; % excentricity of the gaussian
      % (length/width ratio)
y = chirplet (NN, s, e); %generate the chirp
z=exp(1i*2*pi*NN/8*(1:NN)/NN); %for frequency shift
y=y.*z; %shift to the positive frequency
cc=10; %curvature of the chirp
y=y.*exp(1i*2*exp(cc*(1:NN)/NN));
FS = 44100; %sampling rate for visualization


%% initialization of the Gabor window
L = length(y);
%lattice parameters of the Gabor transform
a = 50; % time step
M = L/4; % number of frequency channels
% initial window
Lgamma = L;
g = psech(Lgamma); % Example of window
g = gabtight(g, a, Lgamma); % Example of window
figure;
w = dgt(y, g, a, M); w = 20*log10(abs(w));
imagesc(w); axis xy;
m = max(max(w)); caxis([m-60 m]);
title('Gabor transform of the signal with initial window');


%% Optimization gradient ascent
epsilon = 1e-3; % error tolerance
alpha = 0.001; % step of the gradient optimization
p = 2.5; % order of the l^p norm
% call the optimization for the chirped gaussian family
[gamma, crit] = WinOptimgauss(y, g, a, M, L, p, epsilon, alpha);%


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
