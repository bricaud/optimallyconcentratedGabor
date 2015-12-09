function [gamma, s, d, crit] = WinOptimBFGSOptLat(f, a, M, Lt, L, p, epsilon, d_init, s_init)
% WinOptimBFGS: find gaussian window optimally concentrating the DGT of a target signal
% by using a BFGS optimization scheme
% ---------
%   usage: [gamma, s, d, crit] = WinOptimBFGS(f, a, M, L, p, epsilon)
%
%   find the optimal window with a BFGS method in the subclass of
%   chirped Gaussians
%
%   input:
%   f: target signal
%   a,M: gabor frame lattice parameters, 'a' is the time step, M the number
%   of freq. channels
%   p: p of the L^p norm
%   epsilon: stopping criterion
%
%   output:
%   gamma: optimal window (generally complex valued)
%   s: shear of the optimal window
%   d: spread of the optimal window (std deviation)
%   crit: value of the objective function
% 
% BFGS descent for 
%
%    -norm_Lp( Vg(f) ) + constraint on the time length of the window

%% Initialization
% ------------------------
NbIterMax = 3000;
crit = zeros(NbIterMax, 1);

%Initializarion of estimated hessian matrix
Bk = eye(2);

%% penalization parameter for the support length of the window
% Not used anymore, the penalization is done thanks to the limitation to
% the class of chirped gaussian
% But we keep the code, maybe we will use it one day

% sizeWindow = 4000; % window support 
% weightWindow = fftshift(1-pgauss(sizeWindow,100)/max((pgauss(sizeWindow,100))));
% weight = zeros(L,1);
% weight(1:sizeWindow/2) = weightWindow(1:sizeWindow/2);
% weight(end - sizeWindow/2: end) = weightWindow(end - sizeWindow/2 :  end);

%No penalization, weight vector is set to one
weight = ones(L,1);


%% Set negative frequencies of the signal to zero
ff = fft(f);
ff(ceil(length(ff))/2+1:end) = 0;
ff(2:ceil(length(ff))/2-1) = 2*ff(2:ceil(length(ff))/2-1);
f = ifft(ff);

%% Time Frequency centering of the signal
f = tfcenter(f);


%% time vector
t = -(L)/2 : 1 : (L-1)/2;
t = t/sqrt(L);
%Initialization fot BFGS
s=s_init;
d=d_init;

%Initialisation of the window, and computation of the initial
%value of the objective function
gamma=fftshift(chirplet(L, s, d)).';
shearLat = Lt(1)/Lt(2);
channelsNumber =  round(M/(1+shearLat));

chirpGaussianNeg =  (exp(1i * pi * -shearLat * (t.^2))).';
chirpGaussianPos =  (exp(1i * pi * shearLat * (t.^2))).';

%Chirp the signal
f_init = f;
f = f .* chirpGaussianNeg; 

%Chirp the window
gamma = gamma .* chirpGaussianNeg;

tmp1 = abs(dgt(gamma.*weight, f, a, M * (Lt(2)-Lt(1))/Lt(2)));
fobj = -sum(tmp1(:).^p);

%Computation of the first gradient


% derivatives wrt s and d
gprims=1i*pi*t.^2.*chirplet(L,s,d);
gprimd=sign(d)*(pi*t.^2/d^2-1/(4*abs(d))).*chirplet(L,s,d);
    
% gradient computation
mask = tmp1.^(p-2); 
c=dgt(gamma .* weight,f,a, L * M / (L + M *a * shearLat));
c=c.*mask;
fprocessed=idgt(c,f,a);
fprocessed = fprocessed .* chirpGaussianPos;
grad1 = (p/2) * fprocessed .*weight;
grads = -2*real(fftshift(conj(gprims))*grad1) ;
gradd = -2*real(fftshift(conj(gprimd))*grad1);
nextgrad = [grads; gradd];


%% BFGS loop
%Initialisation of the update step
alpha_k = 1;
for iter = 1:NbIterMax,
    figure(2);
    riplot(gamma); title(['Window at loop ' num2str(iter)]);
    
    %Last values of d and s
    old_s = s;
    old_d = d;
        
    %Compute the update direction from the gradient computation and from
    %the hessian estimation
    grad = nextgrad;
    direction = -Bk \ grad;
    direction = normalize(direction);
    direction_s = direction(1);
    direction_d = direction(2);
        
    
    % linear step (alpha) search
    %Lower bound on the step
    alpha_l = 0;
    %Upper bound on the step
    alpha_r = 1e10;
    beta1 = 0.3;
    beta2 = 0.6;
    %Wolfe conditions
    condition1 = false;
    condition2 = false;
    
    oldfobj = fobj; 
    %Loop in order to find a step which meets the wolfe conditions
    while (((~condition1 || ~condition2) && alpha_r - alpha_l > 1e-10 ))

        % parameters update
        new_d=d+alpha_k*direction_d;
        %d must not be nul
        if (new_d==0)
            new_d= eps;
        end;
        
        new_s=s+alpha_k*direction_s;

        % new window
        chirpGaussianNeg = (exp(1i * pi * -new_s * (t.^2))).';
        chirpGaussianPos = (exp(1i * pi * new_s * (t.^2))).';

        %Chirp the signal
        f = f_init;
        f = f .* chirpGaussianNeg; 

        %Chirp the window
        new_gamma=fftshift(chirplet(L,new_s,new_d)).* chirpGaussianNeg.';
        new_gamma=new_gamma.';

        %%Checke the wolfe conditions
        tmp1 = abs(dgt(new_gamma.*weight, f, a, L * M / (L + M *a * shearLat)));
        fobj = -sum(tmp1(:).^p);
        
        %Condition 1
        if (fobj <= oldfobj + alpha_k * beta1 * grad' * direction);
            condition1 = true;
        else
            condition1 = false;
        end;
        
        % Computation of the gradient at new_s and new_d
        % derivatives wrt s and d
        gprims=1i*pi*t.^2.*chirplet(L,new_s,new_d);
        gprimd=sign(new_d)*(pi*t.^2/new_d^2-1/(4*abs(new_d))).*chirplet(L,new_s,new_d);
    
        % gradient
        mask = tmp1.^(p-2); 
        c=dgt(new_gamma.*weight, f, a, L * M / (L + M *a * shearLat));
        c=c.*mask;
        fprocessed=idgt(c,f,a);
        fprocessed = fprocessed .* chirpGaussianPos; 
        grad1 = -(p/2) * fprocessed .*weight;
        grads = 2*real(fftshift(conj(gprims))*grad1) ;
        gradd = 2*real(fftshift(conj(gprimd))*grad1);    
        nextgrad = [grads; gradd];
        
        %Condition 2
        if (nextgrad' * direction >=  beta2 * grad' * direction);
            condition2 = true;
        else
            condition2 = false;
        end;     
        
        %If one of the condition is not met, then update the step
        if (~condition1)
           alpha_r = alpha_k;
           alpha_k = (alpha_l + alpha_r)/2;
        else if (~condition2)
                alpha_l = alpha_k;
                if alpha_r < 1e10
                    alpha_k = (alpha_l + alpha_r)/2; 
                else
                    alpha_k = 100 * alpha_k;
                end
             end
        end
    end
    
    if alpha_r - alpha_l <= 1e-10
        fprintf('Failing to meet Wolfe conditions');
    end;
    
    %Estimate the hessian for the next step
    sk = alpha_k * direction;
    yk = nextgrad - grad;
    Bk = Bk + (yk * yk.')/(yk.' * sk) - (Bk*sk*sk.'*Bk)/(sk.' * Bk * sk);
    
    %Compute the objective function
    fobj = -sum(tmp1(:).^p);
    crit(iter) = fobj;
    
    fprintf('Iteration %d: , error term d:%e, error term s:%e, crit:%e, alpha:%e\n', iter, old_s - new_s, old_d - new_d, crit(iter), alpha_k);
    
    % stop condition
    if (norm(grad)<epsilon || alpha_k < 1e-10)
        gamma = new_gamma;
        s = new_s;
        d = new_d;
        break;
    end;
    
    %%Update the window and the parameters
    %and analyze the signal with the new window
    s=new_s;
    d=new_d;
    %tmp1 = dgt(gamma.*weight, f, a, M);
    chirpGaussianNeg = (exp(1i * pi * -s * (t.^2))).';
    chirpGaussianPos = (exp(1i * pi * s * (t.^2))).';
    gamma = new_gamma.*chirpGaussianPos;
    
    %Chirp the signal
    f = f_init;
    f = f .* chirpGaussianNeg; 

    %Chirp the window
    gammabis=gamma .* chirpGaussianNeg;

    tmp1 = dgt(gammabis.*weight, f, a,  L * M / (L + M *a * shearLat));
    tmp1 = abs(tmp1); 

    figure(3);
    w2 = dgt(fftshift(gammabis), gammabis, a, M); w2 = 20 * log10(abs(w2));
    imagesc(w2); axis xy; title(['Ambiguity function window at loop ' num2str(iter)]);
    m = max(max(w2)); caxis([m-60 m]);

end;

% time-frequency centering
gamma = tfcenter(gamma);

figure;
riplot(gamma);
title('Optimal window Re + Im');

