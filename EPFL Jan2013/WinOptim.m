function [gamma, crit, entrop] = WinOptim(f, g, a, M, L, p, epsilon, alpha)
% WinOptim: find window optimally concentrating the DGT of a target signal
% ---------
%   usage: [gamma,crit,entrop] = WinOptim(f, g, a, M, L, p, epsilon, alpha)
%
%   find the optimal window with a gradient method
%
%   input:
%   f: target signal
%   g: initial window
%   a,M: gabor frame lattice parameters, 'a' is the time step, M the number
%   of freq. channels
%   p: p of the L^p norm
%   epsilon: stopping criterion
%   alpha: step of the gradient method
%
%   output:
%   gamma: optimal window (generally complex valued)
%
%
% gradient ascent for 
%
%    norm_Lp( Vg(f) ) + ...

%% Initialization
% ------------------------
NbIterMax = 3000;
gamma = g;
crit = zeros(NbIterMax, 1);
crit1 = zeros(NbIterMax, 1);
crit2 = zeros(NbIterMax, 1);
entrop = zeros(NbIterMax, 1);

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

figure;
w2 = dgt(f, gamma, a, M);
w2 = 20 * log10(abs(w2));
imagesc(w2); axis xy;  m = max(max(w2)); caxis([m-60 m]);
title('Gabor transform of initial signal analytic and TF centered');

%% gradient ascent

for iter = 1:NbIterMax,
    figure(2);
    riplot(gamma); title(['Window at loop ' num2str(iter)]);

    % 1st term gradient calculation
    tmp1 = abs(dgt(gamma.*weight, f, a, M));
    
    if (p==2)
        mask = - log(tmp1+eps);
    else
        mask = tmp1.^(p-2);
    end;
    
    grad1 = (p/2) * gabmul(gamma.*weight, mask, f, a);
        
    % gradient
    grad = grad1;% + grad2;
    
    % Optimization criterion (always increase)
    crit1(iter) = sum(tmp1(:).^p);
    %crit2_t(iter) = sum(tmp2_t .* conj(gamma));
    %crit2_f(iter) = sum(fft(tmp2_f) .* conj(fft(gamma)));
    %crit2(iter) = crit2_t(iter) + crit2_f(iter);
    crit(iter) = crit1(iter);%- lambda * crit2(iter);
   
    
    fmax = 0;
    %alphaCurr = max(norm(grad), 1 / norm(grad));
    %alphamax = alphaCurr;
   
    % update window 
    gamma2 = gamma + alpha*grad .* weight;
    %   window normalization
    gamma2 = gamma2 / norm(gamma2);
    tmp1 = abs(dgt(gamma2.*weight, f, a, M));
    fobj = sum(tmp1(:).^p);
    
    if (fobj > fmax)
           fmax = fobj;
    end

    gamma2 = gamma + alpha*grad .* weight;
    %   window normalization
    gamma2 = gamma2 / norm(gamma2);

    tmp1 = abs(dgt(gamma2.*weight, f, a, M));
    crit1(iter) = sum(tmp1(:).^p);
    crit(iter) = crit1(iter);%- lambda * crit2(iter);
    
    %riplot(gamma2);
    fprintf('Iteration %d: , error term:%e, crit:%e, c1:%e, alpha:%e\n', iter, norm(gamma-gamma2), crit(iter), crit1(iter), alpha);
    
    % stop condition
    if norm(gamma-gamma2)<epsilon
        gamma = gamma2;
        break;
    end;
    
    gamma = gamma2;
    
    figure(3);
    w2 = dgt(fftshift(gamma), gamma, a, M); w2 = 20 * log10(abs(w2));
    imagesc(w2); axis xy; title(['Ambiguity function window at loop ' num2str(iter)]);
    m = max(max(w2)); caxis([m-60 m]);
%    figure(4);
%    riplot(grad);
end;

% time-frequency centering
%riplot(gamma);
gamma = tfcenter(gamma);
%riplot(gamma);

figure;
riplot(gamma);
title('Optimal window Re + Im');

% figure;
% w2 = dgt(f, gamma, a, M, L); w2 = 20 * log10(abs(w2));
% imagesc(w2); axis xy; m = max(max(w2)); caxis([m-60 m]);
% title('Gabor Transform of TFcentered signal with optimal window');

crit = [crit crit1 crit2];

