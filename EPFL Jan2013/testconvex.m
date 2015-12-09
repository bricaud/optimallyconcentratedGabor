
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
f=y;

%% initialization of the Gabor window
L = length(y);
%lattice parameters of the Gabor transform
a = 50; % time step
M = L/4; % number of frequency channels
% initial window
Lgamma = L;
g = psech(Lgamma); % Example of window
g = gabtight(g, a, Lgamma); % Example of window
epsilon = 1e-3; % error tolerance
alpha = 1; % step of the gradient optimization
p = 2.5; % order of the l^p norm
%% Initialization
% ------------------------
NbIterMax = 3000;
gamma = g;


%% penalization parameter for the support length of the window
sizeWindow = 4000; % window support 
weightWindow = fftshift(1-pgauss(sizeWindow,100)/max((pgauss(sizeWindow,100))));
weight = zeros(L,1);
weight(1:sizeWindow/2) = weightWindow(1:sizeWindow/2);
weight(end - sizeWindow/2: end) = weightWindow(end - sizeWindow/2 :  end);



%% Set negative frequencies of the signal to zero
ff = fft(f);
ff(ceil(length(ff))/2+1:end) = 0;
ff(2:ceil(length(ff))/2-1) = 2*ff(2:ceil(length(ff))/2-1);
f = ifft(ff);

figure;
w2 = dgt(f, gamma, a, M); w2 = 20 * log10(abs(w2));
imagesc(w2); axis xy; m = max(max(w2)); caxis([m-60 m]);
title('Gabor transform of initial signal analytic"-ized" with initial window');


%% Time Frequency centering of the signal
f = tfcenter(f);



%% gradient ascent
t = -(L-1)/2 : 1 : (L-1)/2;
t = t/sqrt(L);

s=0;
d=10;
new_s=s;
new_d=d;
sii=0;
dii=0;
for si= -50:1:50
    sii=sii+1
    dii=0;
for di = 10:10,
   dii=dii+1
    gamma2=fftshift(chirplet(L,si,di));
    % Optimization criterion (always increase)
    gamma2=gamma2.';
 
    tmp1 = abs(dgt(gamma2.*weight, f, a, M));
  
    crit1(sii,dii) = sum(tmp1(:).^p);
 
end;
end




