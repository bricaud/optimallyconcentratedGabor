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

[yfull FS] = wavread('RPMRaiseExp1Harm.wav');
M = 36750; % number of frequency channels
L = floor(length(yfull)/M)*M;
%L = 160000;
%L = length(yfull);
y = yfull(20000:L);
L = length(y);
Lpad = 147000;
y = [y ; zeros(Lpad - L, 1)];

%% initialization of the Gabor window
L = length(y);
Linit = L;
%Desired redundancy of the dgt transform
R = 6;
%initial window
Lgamma = L;
d = 600;
s = 0;
a = 0;
M = 0;
gamma = chirplet(L, s, d); % Example of window
%figure;

%On sauvegarde l'�volution des param�tres
L_list = [];
a_list = [];
M_list = [];    
Lt_list = [];
d_list = d;
s_list = s;
fobj_list = [];
fprintf(['L initial = ' int2str(Lgamma) '\n']);
max_it = 50;
for it = 1 : 50
    %lattice parameters of the Gabor transform
    [new_Lpad new_a new_M new_Lt] = optimalsampling_new(L, R, abs(d), s, 1.05*Linit);
    new_Lt(1) = mod(new_Lt(1), new_Lt(2));
    sentence = ['lattice n. ' int2str(it) ' : Lt = [' int2str(new_Lt) '] - L = ' int2str(new_Lpad) ...
        ' - a = ' int2str(new_a) ' - M = ' int2str(new_M) '\n\t d = ' num2str(abs(d)) ' - s = ' num2str(s) '\n'];
    fprintf(sentence)
    if (new_Lpad == L && new_a == a && new_M == M && new_Lt(1) == Lt(1) && new_Lt(2) == Lt(2))
        break;
    end;
    
    Lpad = new_Lpad;
    a = new_a;
    M = new_M;
    Lt = new_Lt;
    
    %Zero Padding
    y = [y ; zeros(Lpad - L, 1)];
    L = Lpad;
    
    % gamma = chirplet(L, s, d);
    % w = dgt(y, gamma, a, M, 'lt', Lt);
    % w = 20*log10(abs(w));
    % imagesc(w); axis xy;
    % m = max(max(w)); caxis([m-60 m]);
    % title('Gabor transform of the signal with current window and current lattice');
    
    %% Optimization gradient ascent
    epsilon = 1e-5; % bound on the norm of the gradient (if less than epsilon then stop)
    p = 2.5; % order of the l^p norm
    % call the optimization for the chirped gaussian family
    [gamma, s, d, crit, fobj] = WinOptimBFGSOptWinFast(y, a, M, Lt, L, p, epsilon, d, s);
    
    
    L_list = [L_list  L];
    a_list = [a_list  a];
    M_list = [M_list M];
    Lt_list = [Lt_list Lt];
    d_list =[d_list d];
    s_list = [s_list s];
    fobj_list = [fobj_list fobj];
end;


% %% Visualization
%  axisX = (((1:length(y)/a)*a - a)/FS) ;
%  axisY = (0:FS/M:FS-1);
% % w2 = dgt(y, gamma, a, M);
% % w2 = 20 * log10(abs(w2));
% % figure; imagesc(axisX,axisY, w2); axis xy; m = max(max(w2)); caxis([m-60 m]);
% % title('Gabor Transform of the chirp with optimal window');
% % xlabel('Time (s)');
% % ylabel('Frequency (Hz)');
% 
% gg = fftshift(hanning(1024));
% w2 = dgt(y, gg, a, M);%
% w2 = 20 * log10(abs(w2));
% figure; imagesc(axisX,axisY, w2); axis xy; m = max(max(w2)); caxis([m-60 m]);
% title('Gabor Transform of the chirp with Gaussian window');
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% w2 = dgt(gamma, gamma, a, M);%, length(sig2(1:150000)));
% w2 = fftshift(w2);
% w2 = 20 * log10(abs(w2));
% 
% figure; imagesc(axisX,axisY, w2); axis xy; m = max(max(w2)); caxis([m-60 m]);
% title('Ambiguity function of optimal window');
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% 

%Lattice adaptation
% axisX = (((1:length(y)/a)*a - a)/FS) ;
% axisY = (0:FS/M:FS-1);
% w2 = dgt(y, gamma, a, M, 'lt', Lt);
% w2 = 20 * log10(abs(w2));
% figure();imagescMatrixOnNonSeparableLattice(w2,Lt,a*size(w2,2)/FS,FS);
% title('Gabor Transform of the chirp with optimal window and optimal lattice');
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');

