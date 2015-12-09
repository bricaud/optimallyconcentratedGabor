function [gamma, crit] = WinOptimgauss(f, g, a, M, L, p, epsilon, alpha)
% WinOptimgauss: find window optimally concentrating the DGT of a target signal
% ---------
%   usage: [gamma,crit,entrop] = WinOptim(f, g, a, M, L, p, epsilon, alpha)
%
%   find the optimal window with a gradient method in the subclass of
%   chirped Gaussians
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
%   crit: value of the objective function
% 
% gradient ascent for 
%
%    norm_Lp( Vg(f) ) + constraint on the time length of the window

%% Initialization
% ------------------------
NbIterMax = 3000;
gamma = g;
crit = zeros(NbIterMax, 1);
epsilon_d = 1e-2;
epsilon_s = 1e-4;

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


%% gradient ascent initialization
t = -(L)/2 : 1 : (L-1)/2;
t = t/sqrt(L);

s=0;
d=10;
new_s=s;
new_d=d;
%Reinitization of alpha
alpha=1;
tmp1 = abs(dgt(gamma.*weight, f, a, M));

%% gradient ascent loop
for iter = 1:NbIterMax,
    figure(2);
    riplot(gamma); title(['Window at loop ' num2str(iter)]);
    
    %Last values of d and s
    old_s = s;
    old_d = d;
    
    %Increase alpha
    alpha = alpha * 10;
    
    % 1st term gradient calculation
    fcurrent = sum(tmp1(:).^p);
    
    if (p==2)
        mask = - log(tmp1+eps);
    else
        mask = tmp1.^(p-2);
    end;
    
    grad1 = (p/2) * gabmul(gamma.*weight, mask, f, a).*weight;
    % derivatives wrt s and d
    gprims=1i*pi*t.^2.*chirplet(L,s,d);
    gprimd=(pi*t.^2/d^2-1/(4*d)).*chirplet(L,s,d);
    % gradient
    grads = 2*real(fftshift(conj(gprims))*grad1) ;
    gradd = 2*real(fftshift(conj(gprimd))*grad1);
    % gradient update
    new_s=s+alpha*grads;
    new_d=d+alpha*gradd;
    % new window
    gamma2=fftshift(chirplet(L,new_s,new_d));
    gamma2=gamma2.';
    
    tmp1 = abs(dgt(gamma2.*weight, f, a, M));
    fobj = sum(tmp1(:).^p);
    
    
    % optimal step (alpha) search
    while (fobj < fcurrent && alpha > 1e-10)
        alpha = alpha/10;  
        new_s=s+alpha*grads;
        new_d=d+alpha*gradd;
        gamma2=fftshift(chirplet(L,new_s,new_d));
        gamma2=gamma2.';
 
        tmp1 = abs(dgt(gamma2.*weight, f, a, M));
        fobj = sum(tmp1(:).^p);     
    end

    fcurrent = fobj;
    alpha_s = 1;
    s = new_s;
    d = new_d
    % optimal step for s
    while (fobj <= fcurrent && alpha_s > 1e-10)
        alpha_s = alpha_s/10;  
        new_s=s+alpha_s*sign(grads);
        gamma2=fftshift(chirplet(L,new_s,d));
        gamma2=gamma2.';
 
        tmp1 = abs(dgt(gamma2.*weight, f, a, M));
        fobj = sum(tmp1(:).^p);     
    end
    
    s = new_s
    alpha_d = 1;
    fcurrent = fobj;
    %optimal step for d
    while (fobj <= fcurrent && alpha_d > 1e-10)
        alpha_d = alpha_d/10;  
        new_d=d+alpha_d*sign(gradd);
        gamma2=fftshift(chirplet(L,new_s,new_d));
        gamma2=gamma2.';
    
        tmp1 = abs(dgt(gamma2.*weight, f, a, M));
        fobj = sum(tmp1(:).^p);     
    end
 
    
    crit(iter) = sum(tmp1(:).^p);
   
    
    %riplot(gamma2);
    fprintf('Iteration %d: , error term d:%e, error term s:%e, crit:%e, alpha:%e\n', iter, old_s - s, old_d - d, crit(iter), alpha);
    
    % stop condition
    if (abs(old_s-new_s)<epsilon_s) && (abs(old_d-new_d)<epsilon_d)
        gamma = gamma2;
        s = new_s;
        d = new_d;
        break;
    end;
    
    gamma = gamma2;
    s=new_s;
    d=new_d;
    
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


