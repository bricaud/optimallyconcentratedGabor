clear all;
fileToLoad = '7. Multitone_example.wav.mat';

%% Data from optimal chirped window and optimal gaussian
% Load
load(['../results/WindowAdaptation/' fileToLoad]);
L1 = L;
% Compute spectros
gamma = fftshift(chirplet(L, s, d)).';
optimal = dgt(y, gamma, a, M);
optimal = 20 * log10(abs(optimal).^2);

g=gausswin(8192);
gg = fftshift(g)/sqrt(sum(g.^2));
gaussian = dgt(y, gg, a, M);%
gaussian = 20 * log10(abs(gaussian).^2);

% Find Main Harmonics and their amplitudes
m = max(max(optimal));
[~, amplitudesOptWin, ~] = findMainHarmonics(optimal,'main',m-30);
m = max(max(gaussian));
[~, amplitudesGauss, ~] = findMainHarmonics(gaussian,'main',m-30);

% Plot beautiful comparison figure
axisX = (((1:length(amplitudesOptWin))*a - a)/FS) ;
figure(); hold on; 
plot(axisX, amplitudesOptWin);
plot(axisX, amplitudesGauss,'r:');


%% Data from optimal chirped window on optimal lattice
% Load
load(['./' fileToLoad]);
Lt = [1 5];
% Compute spectros
gamma = fftshift(chirplet(L, s, d)).';
w2 = dgt(y, gamma, a, M, 'lt', Lt);
w2 = 20 * log10(abs(w2).^2);

% Find Main Harmonics and their amplitudes
m = max(max(w2));
[~, amplitudesOptWinLat, ~] = findMainHarmonics(w2,'main',m-30);

% Plot beautiful comparison figure
axisX = (((1:length(amplitudesOptWinLat))*a - a)/FS) ;
if L>L1,
    plot(axisX(1:length(amplitudesOptWinLat)*L1/L), amplitudesOptWinLat(1:length(amplitudesOptWinLat)*L1/L),'m*-');
else 
    plot(axisX(1:length(amplitudesOptWinLat)), amplitudesOptWinLat(1:length(amplitudesOptWinLat)),'m*-');
end
legend('optimal window','gaussian','optimal window and lattice');
title(fileToLoad);
