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

[yfull FS] = wavread('RPMRaiseLin2Harm.wav');
M = 20000; % number of frequency channels
L = floor(length(yfull)/M)*M;
y = yfull(1:40000);
L = length(y);
ynoise = randn(size(y),1)/2;
%y = y + ynoise;
Lpad = 147000;
%y = [y ; zeros(Lpad - L, 1)];
%% initialization of the Gabor window
L = length(y);
%lattice parameters of the Gabor transform
a = 200; % time step
Lt = [0 1];
% initial window
Lgamma = L;
% g = psech(Lgamma); % Example of window
% g = gabtight(g, a, Lgamma); % Example of window
% figure;
% w = dgt(y, g, a, M); 
% w = 20*log10(abs(w));
% initial = w;
% w = initial;
% imagesc(w); axis xy;
% %m = max(max(w)); caxis([m-60 m]);
% title('Gabor transform of the signal with initial window');
% % [frequencies, n] = findMainHarmonics(w,'main',m-30);
% hold on;

% rpmsquare = (0:0.0005:80).^2;
% rpmlin = 0:1/16:10000;
% harmonic2 = harmonic(2, 0.2, 0.9);
% [~,SquareOneHarmonic] = harmonic2.Synthesize(rpmlin,FS);
% SquareOneHarmonic = SquareOneHarmonic(20001:1167:160000);
% 
% plot(frequencies','k.');
% distance = measureDistanceToTrueHarmonics(frequencies,SquareOneHarmonic);

%% Optimization gradient ascent
epsilon = 1e-5; % bound on the norm of the gradient (if less than epsilon then stop)
p = 2.5; % order of the l^p norm
% call the optimization for the chirped gaussian family
[gamma, s, d, crit] = WinOptimBFGS(y, a, M, [0 1], L, p, epsilon);%

%% Visualization
axisX = (((1:length(y)/a)*a - a)/FS) ;
axisY = (0:FS/M:FS-1);
w2 = dgt(y, gamma, a, M);
w2 = 20 * log10(abs(w2));
optimal = w2;
w2=optimal;
figure; imagesc(axisX,axisY, w2); axis xy; m = max(max(w2)); caxis([m-60 m]);
title('Gabor Transform of the chirp with optimal window');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
hold on;

% % find main harmonics and plot
% [frequencies, n] = findMainHarmonics(w2,'main',m-30);
% plot(axisX,frequencies'*FS/M,'k.');
% distance = measureDistanceToTrueHarmonics(frequencies,SquareOneHarmonic);
% 
% estimate peak width and plot
% w2filt=filter(ones(1,10)/10,1,w2);
% widthopt = computePeakWidth(w2filt);
% for i=1:size(widthopt,1),
%     for j=1:size(widthopt,2),
%         f=frequencies(i,j);
%         w=widthopt(i,j);
%         line([j j]*a/FS-a/FS,[f-w f+w]*FS/M,'Color','k','LineWidth',2);
%     end
% end


gg = fftshift(gausswin(16384));
w2 = dgt(y, gg, a, M);%
w2 = 20 * log10(abs(w2));
gaussian = w2;
w2 = gaussian;
figure; imagesc(axisX,axisY, w2); axis xy; m = max(max(w2)); caxis([m-60 m]);
title('Gabor Transform of the chirp with Gaussian window');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
hold on;

% [frequencies, n] = findMainHarmonics(w2,'main',m-30);
% plot(axisX,frequencies'*FS/M,'k.');
% distance = measureDistanceToTrueHarmonics(frequencies,SquareOneHarmonic);

% estimate peak width and plot
% w2filt=filter(ones(1,10)/10,1,w2);
% widthgauss = computePeakWidth(w2filt);
% for i=1:size(widthgauss,1),
%     for j=1:size(widthgauss,2),
%         f=frequencies(i,j);
%         w=widthgauss(i,j);
%         line([j j]*a/FS-a/FS,[f-w f+w]*FS/M,'Color','w','LineWidth',2);
%     end
% end


% w2 = dgt(gamma, gamma, a, M);%, length(sig2(1:150000)));
% w2 = fftshift(w2);
% w2 = 20 * log10(abs(w2));
% figure; imagesc(axisX,axisY, w2); axis xy; m = max(max(w2)); caxis([m-60 m]);
% title('Ambiguity function of optimal window');
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% 
% %Lattice adaptation
% [Lpad, a_ad, M_ad, Lt] = optimalsampling_new(L, M/a, d,s);
% ypad = [y; zeros(Lpad-length(y), 1)];
% axisX = (((1:length(ypad)/a_ad)*a_ad - a_ad)/FS) ;
% axisY = (0:FS/M_ad:FS-1);
% w2 = dgt(ypad, gamma, a_ad, M_ad, 'lt', Lt);
% w2 = 20 * log10(abs(w2));
% figure; imagesc(axisX,axisY, w2); axis xy; m = max(max(w2)); caxis([m-60 m]);
% title('Gabor Transform of the chirp with optimal window and optimal lattice');
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% 


%% Comparisions

% With a classical Hamming window STFT 
w3 = spectrogram(y, a, 0, M);
w3 = 20 * log10(abs(w3));
figure; imagesc(axisX,axisY, w3); axis xy; m = max(max(w3)); caxis([m-60 m]);
title('Spectrogram of the chirp with Hamming window');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
% [frequencies, n] = findMainHarmonics(w3,'main',m-30);
hold on;
% plot(axisX,frequencies'*FS/M,'k.');
% distance = measureDistanceToTrueHarmonics(frequencies,SquareOneHarmonic);

% estimate peak width and plot
% w3filt=filter(ones(1,10)/10,1,w3);
% widthspectr = computePeakWidth(w3filt);
% for i=1:size(widthgauss,1),
%     for j=1:size(widthgauss,2),
%         f=frequencies(i,j);
%         w=widthspectr(i,j);
%         line([j j]*a/FS-a/FS,[f-w f+w]*FS/M*2,'Color','k','LineWidth',2);
%     end
% end

clear w w2 w3 m g initial
% Gabor with the optimal classical gaussian window


% Wigner Ville


% Quadratic generalisation, Altes, Bertrand