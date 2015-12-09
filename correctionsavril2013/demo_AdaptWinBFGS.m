% DEMO_ADAPTWIN_BFGS
% ==============
%
% Script for optimizing a window in terms of the concentration of
% a target function, estimated by restricting the gabor transform of a
% target (chirped gaussian) signal to a given time-frequency region
%
% It uses a BFGS method


set(0,'DefaultAxesFontSize',18)

%% Target signal
% Example 1: a bended and 'chirped' gaussian
% NN=12000;
% s=0.2; % shear or slope of the chirp
% e=48; % excentricity of the gaussian
%       % (length/width ratio)
% y = chirplet (NN, s, e); %generate the chirp
% z=exp(1i*1*pi*NN/20*(1:NN)/NN); %for frequency shift
% y=0.5 * y.*z + 0.5 * chirplet (NN, -1, e).*z; %shift to the positive frequency
% cc=10; %curvature of the chirp
% y=y.*exp(1i*1*exp(cc*(1:NN)/NN));
% FS = 44100; %sampling rate for visualization
% y = y.';

[yfull FS] = wavread('RPMRaiseSquare10Harm.wav');
y = yfull(20001:160000);
L = length(y);
Lpad = 147000;
y = [y ; zeros(Lpad - L, 1)];
%ynoise = randn(10000,1);
%y = y + ynoise/5;

%% initialization of the Gabor window
L = length(y);
%lattice parameters of the Gabor transform
a = 1225; % time step
M = 36750; % number of frequency channels
Lt = [0 1];
% initial window
Lgamma = L;
g = psech(Lgamma); % Example of window
g = gabtight(g, a, Lgamma); % Example of window
figure;
w = dgt(y, g, a, M); 
w = 20*log10(abs(w));
imagesc(w); axis xy;
m = max(max(w)); caxis([m-60 m]);
title('Gabor transform of the signal with initial window');


%% Optimization gradient ascent
epsilon = 1e-5; % bound on the norm of the gradient (if less than epsilon then stop)
p = 2.5; % order of the l^p norm
% call the optimization for the chirped gaussian family
[gamma, s, d, crit] = WinOptimBFGS(y, a, M, Lt, L, p, epsilon);%


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
w2 = dgt(gamma, gamma, a, M);%, length(sig2(1:150000)));
w2 = fftshift(w2);
w2 = 20 * log10(abs(w2));

figure; imagesc(axisX,axisY, w2); axis xy; m = max(max(w2)); caxis([m-60 m]);
title('Ambiguity function of optimal window');
xlabel('Time (s)');
ylabel('Frequency (Hz)');


%Lattice adaptation
[Lpad, a_ad, M_ad, Lt] = optimalsampling_new(L, M/a, d,s);
ypad = [y; zeros(Lpad-length(y), 1)];
axisX = (((1:length(ypad)/a_ad)*a_ad - a_ad)/FS) ;
axisY = (0:FS/M_ad:FS-1);
w2 = dgt(ypad, gamma, a_ad, M_ad, 'lt', Lt);
w2 = 20 * log10(abs(w2));
figure; imagesc(axisX,axisY, w2); axis xy; m = max(max(w2)); caxis([m-60 m]);
title('Gabor Transform of the chirp with optimal window and optimal lattice');
xlabel('Time (s)');
ylabel('Frequency (Hz)');

